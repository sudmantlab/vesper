"""Command module for variant refinement functionality."""

from pathlib import Path
import logging
import time
import os
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor

from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn, TimeRemainingColumn
from vesper.processors.vcf import VCFProcessor, VCFWriter
from vesper.processors.reads import ReadProcessor
from vesper.utils.config import RefineConfig

def process_refine_chunk(variants: list, read_proc: ReadProcessor, chunk_idx: int, logger: logging.Logger) -> None:
    """Process a chunk of variants through the refinement pipeline.
    
    Args:
        variants: List of variants to process
        read_proc: ReadProcessor instance
        chunk_idx: Index of this chunk for logging
        logger: Logger instance
    """
    logger.info(f"Processing refinement chunk {chunk_idx} ({len(variants)} variants)")
    
    for variant in variants:
        support_reads, nonsupport_reads = read_proc.get_read_groups(variant.variant)
        variant.support_reads = support_reads
        variant.nonsupport_reads = nonsupport_reads
        
        variant._calculate_grouped_metrics()
        variant._calculate_confidence()
    
    # Save registry after each chunk (if modified)
    read_proc.save_registry()
    logger.info(f"Completed refinement chunk {chunk_idx}")

def run_refine(args, logger):
    """Run the refinement pipeline."""
    config = RefineConfig.from_args(args)
    
    if config.test_mode is not None:
        logger.info(f"Running in test mode (limited to {config.test_mode} variants)")
        print(f"{datetime.now().strftime('%m/%d/%Y %I:%M:%S %p')} - Running in test mode (limited to {config.test_mode} variants)")

    # Create output directory if it doesn't exist
    if not config.output_dir.exists():
        config.output_dir.mkdir(parents=True)
        logger.info(f"Created output directory: {config.output_dir}")
    
    logger.info(f"Loading variants from {config.vcf_input}")
    variants = []
    with VCFProcessor(config.vcf_input, test_mode=config.test_mode) as vcf_proc:
        variants = list(vcf_proc.instantiate_variants())
    print(f"{datetime.now().strftime('%m/%d/%Y %I:%M:%S %p')} - Loaded {len(variants)} variants")
    logger.info(f"Loaded {len(variants)} variants")

    # TODO: set config options for memory
    n_workers = config.threads
    chunk_size = max(10, len(variants) // (n_workers * 16)) # smaller chunks for more granular progress updates
    chunks = [variants[i:i + chunk_size] for i in range(0, len(variants), chunk_size)]
    
    logger.info(f"Processing {len(chunks)} chunks with {n_workers} workers")
    start_time = time.time()
    with Progress(TextColumn("[bold blue]{task.description}"),
                 BarColumn(complete_style="green"),
                 TaskProgressColumn(),
                 TimeElapsedColumn(),
                 TimeRemainingColumn()) as progress, \
         ReadProcessor(config.bam_file, registry_dir=config.output_dir / 'read_registry',
                      auto_load_registry=config.auto_load_registry,
                      force_new_registry=config.force_new_registry) as read_proc, \
         ThreadPoolExecutor(max_workers=n_workers) as executor:
        task = progress.add_task("[cyan]Refining variants...", total=len(variants))  
        
        futures = [
            executor.submit(
                process_refine_chunk, 
                chunk, 
                read_proc,
                idx,
                logger
            )
            for idx, chunk in enumerate(chunks)
        ]
        
        completed_variants = 0
        for future in futures:
            try:
                future.result()
                completed_variants += chunk_size
                progress.update(task, advance=chunk_size)
                logger.info(f"Progress: {completed_variants}/{len(variants)} variants completed ({completed_variants/len(variants)*100:.1f}%)")
                read_proc.save_registry(force=True)
            except Exception as e:
                logger.error(f"Error processing chunk: {str(e)}")
                raise
        
        elapsed = time.time() - start_time
        logger.info(f"Completed refinement in {elapsed:.2f} seconds")

        confidence_scores = [v.confidence for v in variants if v.confidence is not None]
        if confidence_scores:
            mean_conf = sum(confidence_scores) / len(confidence_scores)
            min_conf = min(confidence_scores)
            max_conf = max(confidence_scores)
            median_conf = sorted(confidence_scores)[len(confidence_scores)//2]
            
            confident_variants = [v for v in variants if v.confidence and v.confidence > 0.6]
            pct_confident = len(confident_variants) / len(variants) * 100
            
            logger.info(f"Confidence score statistics:")
            logger.info(f"    Mean: {mean_conf:.3f}")
            logger.info(f"    Median: {median_conf:.3f}") 
            logger.info(f"    Min: {min_conf:.3f}")
            logger.info(f"    Max: {max_conf:.3f}")
            logger.info(f"    High-confidence variants: {len(confident_variants)} ({pct_confident:.1f}%)")
        else:
            logger.warning("WARNING: No confidence scores calculated!")
            print(f"{datetime.now().strftime('%m/%d/%Y %I:%M:%S %p')} - WARNING: No confidence scores calculated!")

    print(f"{datetime.now().strftime('%m/%d/%Y %I:%M:%S %p')} - Completed refinement in {elapsed:.2f} seconds")
    print(f"{datetime.now().strftime('%m/%d/%Y %I:%M:%S %p')} - Confidence score statistics:")
    print(f"    Mean: {mean_conf:.3f}")
    print(f"    Median: {median_conf:.3f}")
    print(f"    Min: {min_conf:.3f}")
    print(f"    Max: {max_conf:.3f}")
    print(f"    High-confidence variants: {len(confident_variants)} ({pct_confident:.1f}%)")

    output_vcf_path = config.output_dir / config.vcf_input.name.replace('.vcf.gz', '.refined.vcf.gz')
    logger.info(f"Writing refined variants to {output_vcf_path}")
    print(f"{datetime.now().strftime('%m/%d/%Y %I:%M:%S %p')} - Writing refined variants to {output_vcf_path}")
    
    start_write_time = time.time()
    with Progress(
        TextColumn("[bold blue]{task.description}"),
        BarColumn(complete_style="green"), 
        TaskProgressColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn()
    ) as progress, VCFWriter(output_vcf_path) as vcf_writer:
        
        task = progress.add_task("[cyan]Writing refined variants...", total=len(variants))
        vcf_writer.write_header(variants)
        
        for i, variant in enumerate(variants):
            vcf_writer.write_record(variant)
            progress.update(task, advance=1)
            
            if (i + 1) % 1000 == 0:
                logger.info(f"Wrote {i + 1}/{len(variants)} variants")
    
    write_elapsed = time.time() - start_write_time
    logger.info(f"Completed writing VCF in {write_elapsed:.2f} seconds")
    
    logger.info(f"Creating tabix index for {output_vcf_path}")
    VCFWriter.create_tabix_index(output_vcf_path)
    
    logger.info(f"Refinement pipeline completed successfully")
    print(f"{datetime.now().strftime('%m/%d/%Y %I:%M:%S %p')} - Refinement pipeline completed successfully") 
"""Command module for variant annotation functionality."""

from pathlib import Path
import logging
import time
from datetime import datetime
import os
from concurrent.futures import ThreadPoolExecutor

from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn, TimeRemainingColumn
from vesper.processors.vcf import VCFProcessor, VCFWriter
from vesper.processors.annotations import BEDProcessor
from vesper.utils.config import AnnotateConfig

def process_annotate_chunk(variants: list, bed_proc: BEDProcessor, chunk_idx: int, proximal_span: int, logger: logging.Logger) -> None:
    """Process a chunk of variants through the annotation pipeline.
    
    Args:
        variants: List of variants to process
        bed_proc: BEDProcessor instance
        chunk_idx: Index of this chunk for logging
        proximal_span: Distance (+/-) to search for proximal features
        logger: Logger instance
    """
    logger.info(f"Processing annotation chunk {chunk_idx} ({len(variants)} variants)")
    
    for variant in variants:
        bed_proc._annotate_variant(variant, proximal_span=proximal_span)
    
    logger.info(f"Completed annotation chunk {chunk_idx}")

def run_annotate(args, logger):
    """Run the annotation pipeline."""
    config = AnnotateConfig.from_args(args)
    timestamp = datetime.now().strftime("%m/%d/%Y %I:%M:%S %p")
    if config.test_mode is not None:
        logger.info(f"Running in test mode (limited to {config.test_mode} variants)")
        print(f"{timestamp} - Running in test mode (limited to {config.test_mode} variants)")

    if not config.output_dir.exists():
        config.output_dir.mkdir(parents=True)
        logger.info(f"Created output directory: {config.output_dir}")
    
    logger.info(f"Loading variants from {config.vcf_input}")
    variants = []
    with VCFProcessor(config.vcf_input, test_mode=config.test_mode) as vcf_proc:
        variants = list(vcf_proc.instantiate_variants())
    logger.info(f"Loaded {len(variants)} variants")
    
    # TODO: set config options for threads, memory, etc.
    n_workers = min(32, max(4, os.cpu_count() * 2))
    chunk_size = max(1, len(variants) // (n_workers * 16)) # smaller chunks for more granular progress updates
    chunks = [variants[i:i + chunk_size] for i in range(0, len(variants), chunk_size)]
    
    logger.info(f"Processing {len(chunks)} variant chunks with {n_workers} workers")
    start_time = time.time()

    with Progress(TextColumn("[bold blue]{task.description}"),
                 BarColumn(complete_style="green"),
                 TaskProgressColumn(),
                 TimeElapsedColumn(),
                 TimeRemainingColumn()) as progress, \
         BEDProcessor(config.bed_file) as bed_proc, \
         ThreadPoolExecutor(max_workers=n_workers) as executor:
        task = progress.add_task("[cyan]Annotating variants...", total=len(variants))           

        futures = [
            executor.submit(
                process_annotate_chunk, 
                chunk, 
                bed_proc,
                idx,
                config.proximal_span,
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
            except Exception as e:
                logger.error(f"Error processing chunk: {str(e)}")
                raise
    
    elapsed = time.time() - start_time
    logger.info(f"Completed annotation in {elapsed:.2f} seconds")
    logger.info(f"Annotated {len(variants)} variants")
    logger.info(f"Mean overlapping features: {sum(len(v.overlapping_features) for v in variants)/len(variants):.1f}") 
    print(f"{timestamp} - Completed annotation in {elapsed:.2f} seconds")
    print(f"{timestamp} - Annotated {len(variants)} variants")
    print(f"{timestamp} - Mean overlapping features: {sum(len(v.overlapping_features) for v in variants)/len(variants):.1f}") 
    
    output_vcf_path = config.output_dir / config.vcf_input.name.replace('.vcf.gz', '.annotated.vcf.gz')
    logger.info(f"Writing annotated variants to {output_vcf_path}")
    print(f"{timestamp} - Writing annotated variants to {output_vcf_path}")
    
    start_write_time = time.time()
    with Progress(
        TextColumn("[bold blue]{task.description}"),
        BarColumn(complete_style="green"), 
        TaskProgressColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn()
    ) as progress, VCFWriter(output_vcf_path) as vcf_writer:
        
        task = progress.add_task("[cyan]Writing annotated variants...", total=len(variants))
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
    
    logger.info(f"Annotation pipeline completed successfully")
    print(f"{timestamp} - Annotation pipeline completed successfully") 
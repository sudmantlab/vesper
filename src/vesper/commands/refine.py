"""Command module for variant refinement functionality."""

from pathlib import Path
import logging
import time
import os
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor

from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn, TimeRemainingColumn
from vesper.processors.vcf import VCFProcessor
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
        # Get read groups
        support_reads, nonsupport_reads = read_proc.get_read_groups(variant.variant)
        variant.support_reads = support_reads
        variant.nonsupport_reads = nonsupport_reads
        
        # Calculate metrics and confidence
        variant._calculate_grouped_metrics()
        variant._calculate_confidence()
    
    # Save registry after each chunk (if modified)
    read_proc.save_registry()
    logger.info(f"Completed refinement chunk {chunk_idx}")

def run_refine(args, logger):
    """Run the refinement pipeline."""
    config = RefineConfig.from_args(args)
    
    timestamp = datetime.now().strftime("%m/%d/%Y %I:%M:%S %p")
    if config.test_mode is not None:
        logger.info(f"Running in test mode (limited to {config.test_mode} variants)")
        print(f"{timestamp} - Running in test mode (limited to {config.test_mode} variants)")

    # Create output directory if it doesn't exist
    if not config.output_dir.exists():
        config.output_dir.mkdir(parents=True)
        logger.info(f"Created output directory: {config.output_dir}")
    
    logger.info(f"Loading variants from {config.vcf_input}")
    variants = []
    with VCFProcessor(config.vcf_input, test_mode=config.test_mode) as vcf_proc:
        variants = list(vcf_proc.instantiate_variants())
    print(f"{timestamp} - Loaded {len(variants)} variants")
    logger.info(f"Loaded {len(variants)} variants")

    # Determine optimal chunk size and number of workers
    n_workers = min(32, max(4, os.cpu_count() * 2))
    chunk_size = max(1, len(variants) // (n_workers * 16)) # smaller chunks for more granular progress updates
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

        confident_variants = [v for v in variants if v.confidence > 0.6]
        logger.info(f"Found {len(confident_variants)} high-confidence variants")
        logger.info(f"Mean confidence score: {sum(v.confidence for v in variants)/len(variants):.3f}") 

    print(f"{timestamp} - Completed refinement in {elapsed:.2f} seconds")
    print(f"{timestamp} - Found {len(confident_variants)} high-confidence variants")
    print(f"{timestamp} - Mean confidence score: {sum(v.confidence for v in variants)/len(variants):.3f}") 
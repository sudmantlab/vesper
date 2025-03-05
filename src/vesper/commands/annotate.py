"""Command module for variant annotation functionality."""

from pathlib import Path
import logging
import time
from datetime import datetime
import os
from concurrent.futures import ThreadPoolExecutor

from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn, TimeRemainingColumn
from vesper.processors.vcf import VCFProcessor
from vesper.processors.annotations import BEDProcessor
from vesper.utils.config import AnnotateConfig

def process_annotate_chunk(variants: list, bed_proc: BEDProcessor, chunk_idx: int, logger: logging.Logger) -> None:
    """Process a chunk of variants through the annotation pipeline.
    
    Args:
        variants: List of variants to process
        bed_proc: BEDProcessor instance
        chunk_idx: Index of this chunk for logging
        logger: Logger instance
    """
    logger.info(f"Processing annotation chunk {chunk_idx} ({len(variants)} variants)")
    
    for variant in variants:
        # Annotate with genomic features
        bed_proc._annotate_variant(variant, proximal_span=500)
    
    logger.info(f"Completed annotation chunk {chunk_idx}")

def run_annotate(args, logger):
    """Run the annotation pipeline."""
    config = AnnotateConfig.from_args(args)
    timestamp = datetime.now().strftime("%m/%d/%Y %I:%M:%S %p")
    if config.test_mode is not None:
        logger.info(f"Running in test mode (limited to {config.test_mode} variants)")
        print(f"{timestamp} - Running in test mode (limited to {config.test_mode} variants)")

    # Create output directory if it doesn't exist
    if not config.output_dir.exists():
        config.output_dir.mkdir(parents=True)
        logger.info(f"Created output directory: {config.output_dir}")
    
    # Load variants from VCF
    logger.info(f"Loading variants from {config.vcf_input}")
    variants = []
    with VCFProcessor(config.vcf_input, test_mode=config.test_mode) as vcf_proc:
        variants = list(vcf_proc.instantiate_variants())
    logger.info(f"Loaded {len(variants)} variants")
    
    # Determine optimal chunk size and number of workers
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
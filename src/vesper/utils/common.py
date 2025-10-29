"""Common pipeline utilities shared between commands."""

from __future__ import annotations

import time
from pathlib import Path
from typing import List, Tuple
import logging

from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn, TimeRemainingColumn
from vesper.processors.vcf import VCFProcessor, VCFWriter


def load_variants(vcf_path: Path, test_mode: int = None, logger: logging.Logger = None) -> List:
    """Load variants from VCF file."""
    if logger:
        logger.info(f"Loading variants from {vcf_path}")
    
    variants = []
    with VCFProcessor(vcf_path, test_mode=test_mode) as vcf_proc:
        variants = list(vcf_proc.instantiate_variants())
    
    if logger:
        logger.info(f"Loaded {len(variants)} variants")
    
    return variants


def calculate_chunks(variants: List, n_threads: int) -> Tuple[List[List], int]:
    """Calculate optimal chunks for parallel processing."""
    chunk_size = max(10, len(variants) // (n_threads * 16))
    chunks = [variants[i:i + chunk_size] for i in range(0, len(variants), chunk_size)]
    return chunks, chunk_size


def setup_output_directory(output_dir: Path, logger: logging.Logger = None) -> None:
    """Create output directory if it doesn't exist."""
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
        if logger:
            logger.info(f"Created output directory: {output_dir}")


def write_vcf_with_progress(output_path: Path, variants: List, logger: logging.Logger = None, interval: int = 100) -> None:
    """Write VCF with progress bar and create tabix index."""
    if logger:
        logger.info(f"Writing variants to {output_path}")
    
    start_write_time = time.time()
    with Progress(
        TextColumn("[bold blue]{task.description}"),
        BarColumn(complete_style="green"), 
        TaskProgressColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn()
    ) as progress, VCFWriter(output_path) as vcf_writer:
        
        task = progress.add_task("[cyan]Writing variants...", total=len(variants))
        vcf_writer.write_header(variants)
        
        for i, variant in enumerate(variants):
            vcf_writer.write_record(variant)
            progress.update(task, advance=1)
            
            if (i + 1) % interval == 0 and logger:
                logger.debug(f"Wrote {i + 1}/{len(variants)} variants")
    
    write_elapsed = time.time() - start_write_time
    if logger:
        logger.info(f"Completed writing VCF in {write_elapsed:.2f} seconds")
        logger.info(f"Creating tabix index for {output_path}")
    
    VCFWriter.create_tabix_index(output_path)

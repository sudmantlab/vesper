"""Command module for variant annotation functionality."""

from pathlib import Path
import logging
import time
from datetime import datetime
import os
from concurrent.futures import ThreadPoolExecutor

from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn, TimeRemainingColumn
from vesper.processors.vcf import VCFProcessor, VCFWriter
from vesper.processors.annotations import BEDProcessor, GFFProcessor
from vesper.processors.repeatmasker import RepeatMaskerProcessor
from vesper.utils.config import AnnotateConfig

def process_annotate_chunk(variants: list, annotation_procs: list, repeatmasker_proc: RepeatMaskerProcessor, chunk_idx: int, proximal_span: int, logger: logging.Logger) -> None:
    """Process a chunk of variants through the annotation pipeline with multiple annotation processors.
    
    Args:
        variants: List of variants to process
        annotation_procs: List of AnnotationProcessor instances (BED or GFF)
        repeatmasker_proc: RepeatMasker processor instance
        chunk_idx: Index of this chunk for logging
        proximal_span: Distance (+/-) to search for proximal features
        logger: Logger instance
    """
    logger.info(f"Processing annotation chunk {chunk_idx} ({len(variants)} variants)")
    
    # Apply overlapping/proximal annotations first
    for variant in variants:
        for proc in annotation_procs:
            proc._annotate_variant(variant, proximal_span=proximal_span)
            
    # Process insertion sequences with RepeatMasker
    repeatmasker_proc.process_batch(variants)
    
    logger.info(f"Completed annotation chunk {chunk_idx}")

def run_annotate(args, logger):
    """Run the annotation pipeline."""
    config = AnnotateConfig.from_args(args)
    timestamp = datetime.now().strftime("%m/%d/%Y %I:%M:%S %p")
    if config.test_mode is not None:
        logger.info(f"Running in test mode (limited to {config.test_mode} variants)")
        print(f"{timestamp} - Running in test mode (limited to {config.test_mode} variants)")
        
    total_annotation_files = len(config.bed_files) + len(config.gff_files)
    logger.info(f"Using {total_annotation_files} annotation file(s):")
    if config.bed_files:
        logger.info(f"BED file(s) ({len(config.bed_files)}):")
        for i, bed_file in enumerate(config.bed_files):
            bed_name = config.bed_names[i] if config.bed_names and i < len(config.bed_names) else bed_file.stem
            logger.info(f"  - {bed_file} (source: '{bed_name}')")
    if config.gff_files:
        logger.info(f"GFF file(s) ({len(config.gff_files)}):")
        for i, gff_file in enumerate(config.gff_files):
            gff_name = config.gff_names[i] if config.gff_names and i < len(config.gff_names) else gff_file.stem
            logger.info(f"  - {gff_file} (source: '{gff_name}')")

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
         ThreadPoolExecutor(max_workers=n_workers) as executor:

        annotation_procs = []
        
        # Use BED names from config (validation already ensures names are available)
        for i, bed_file in enumerate(config.bed_files):
            task = progress.add_task("[cyan]Loading BED files...", total=len(config.bed_files))
            bed_name = config.bed_names[i]
            progress.update(task, description=f"Loading BED file: {bed_file.name} ({bed_name})")
            logger.info(f"Loading BED file: {bed_file.name} as '{bed_name}'")
            annotation_procs.append(BEDProcessor(bed_file, source_name=bed_name))
            progress.update(task, advance=1)
        
        # Use GFF names from config (validation already ensures names are available)
        for i, gff_file in enumerate(config.gff_files):
            task = progress.add_task("[cyan]Loading GFF files...", total=len(config.gff_files))
            gff_name = config.gff_names[i]
            progress.update(task, description=f"Loading GFF file: {gff_file.name} ({gff_name})")
            logger.info(f"Loading GFF file: {gff_file.name} as '{gff_name}'")
            annotation_procs.append(GFFProcessor(gff_file, source_name=gff_name))
            progress.update(task, advance=1)

    with Progress(TextColumn("[bold blue]{task.description}"),
                 BarColumn(complete_style="green"),
                 TaskProgressColumn(),
                 TimeElapsedColumn(),
                 TimeRemainingColumn()) as progress, \
         ThreadPoolExecutor(max_workers=n_workers) as executor:
        
        # Open all processors using context managers
        for proc in annotation_procs:
            proc.__enter__()
        repeatmasker_proc = RepeatMaskerProcessor(config.output_dir)
        repeatmasker_proc.__enter__()
        
        try:
            task = progress.add_task("[cyan]Annotating variants...", total=len(variants))           

            futures = [
                executor.submit(
                    process_annotate_chunk, 
                    chunk, 
                    annotation_procs,
                    repeatmasker_proc,
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
        finally:
            for proc in annotation_procs:
                proc.__exit__(None, None, None)
            repeatmasker_proc.__exit__(None, None, None)
    
    elapsed = time.time() - start_time
    logger.info(f"Completed annotation in {elapsed:.2f} seconds")
    logger.info(f"Annotated {len(variants)} variants using {total_annotation_files} annotation file(s)")
    
    if config.bed_files:
        logger.info(f"BED file(s) ({len(config.bed_files)}):")
        for i, bed_file in enumerate(config.bed_files):
            bed_name = config.bed_names[i]
            logger.info(f"  - {bed_file} (source: '{bed_name}')")
            
    if config.gff_files:
        logger.info(f"GFF file(s) ({len(config.gff_files)}):")
        for i, gff_file in enumerate(config.gff_files):
            gff_name = config.gff_names[i]
            logger.info(f"  - {gff_file} (source: '{gff_name}')")
            
    logger.info(f"Mean overlapping features: {sum(len(v.overlapping_features) for v in variants)/len(variants):.1f}") 
    logger.info(f"Mean proximal features: {sum(len(v.proximal_features) for v in variants)/len(variants):.1f}") 
    
    print(f"{timestamp} - Completed annotation in {elapsed:.2f} seconds")
    print(f"{timestamp} - Annotated {len(variants)} variants using {total_annotation_files} annotation file(s)")
    
    if config.bed_files:
        print(f"{timestamp} - BED file(s) ({len(config.bed_files)}):")
        for i, bed_file in enumerate(config.bed_files):
            bed_name = config.bed_names[i]
            print(f"  - {bed_file} (source: '{bed_name}')")
            
    if config.gff_files:
        print(f"{timestamp} - GFF file(s) ({len(config.gff_files)}):")
        for i, gff_file in enumerate(config.gff_files):
            gff_name = config.gff_names[i]
            print(f"  - {gff_file} (source: '{gff_name}')")
            
    print(f"{timestamp} - Mean overlapping features: {sum(len(v.overlapping_features) for v in variants)/len(variants):.1f}") 
    print(f"{timestamp} - Mean proximal features: {sum(len(v.proximal_features) for v in variants)/len(variants):.1f}") 
    
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
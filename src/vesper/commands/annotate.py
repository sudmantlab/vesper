"""Command module for variant annotation functionality."""

from pathlib import Path
import logging
import time
from datetime import datetime
import os
from glob import glob
from concurrent.futures import ThreadPoolExecutor

from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn, TimeRemainingColumn
from vesper.processors.annotations import BEDProcessor, GFFProcessor, TSVProcessor
from vesper.processors.repeatmasker import RepeatMaskerProcessor
from vesper.utils.config import AnnotateConfig
from vesper.utils.common import load_variants, calculate_chunks, setup_output_directory, write_vcf_with_progress

def process_annotate_chunk(variants: list, annotation_procs: list, chunk_idx: int, config: AnnotateConfig, logger: logging.Logger) -> None:
    """Process a chunk of variants through the annotation pipeline.
    
    Args:
        variants: List of variants to process
        annotation_procs: List of AnnotationProcessor instances (BED/GFF/TSV) - may be empty
        chunk_idx: Index of this chunk for logging
        config: AnnotateConfig instance
        logger: Logger instance
    """
    logger.info(f"Processing annotation chunk {chunk_idx} ({len(variants)} variants)")
    
    # Apply external annotations if provided
    if annotation_procs:
        for variant in variants:
            for proc in annotation_procs:
                proc._annotate_variant(variant, proximal_span=config.proximal_span)
            
    # Always process insertion sequences with RepeatMasker
    with RepeatMaskerProcessor(config.output_dir) as repeatmasker_proc:
        repeatmasker_proc.batch_analysis(
            variants, 
            chunk_idx=chunk_idx,
            n=config.repeatmasker_n,
        )
    
    repeatmasker_proc.__exit__(None, None, None)
    
    logger.info(f"Completed annotation chunk {chunk_idx}")

def run_annotate(args, logger):
    """Run the annotation pipeline."""
    config = AnnotateConfig.from_args(args)

    if config.test_mode is not None:
        logger.info(f"Running in test mode (limited to {config.test_mode} variants)")
        
    total_annotation_files = len(config.bed_files) + len(config.gff_files) + len(config.tsv_files)
    if total_annotation_files > 0:
        logger.info(f"Using {total_annotation_files} external annotation file(s):")
    else:
        logger.info("No external annotation files provided, will run RepeatMasker only")
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
    if config.tsv_files:
        logger.info(f"TSV file(s) ({len(config.tsv_files)}):")
        for i, tsv_file in enumerate(config.tsv_files):
            tsv_name = config.tsv_names[i] if config.tsv_names and i < len(config.tsv_names) else tsv_file.stem
            logger.info(f"  - {tsv_file} (source: '{tsv_name}')")

    setup_output_directory(config.output_dir, logger)
    variants = load_variants(config.vcf_input, config.test_mode, logger)
    
    # TODO: set config options for memory, etc.
    n_threads = config.threads
    chunks, chunk_size = calculate_chunks(variants, n_threads)
    
    
    start_time = time.time()
    logger.info(f"Starting annotation pipeline")
    logger.info(f"Using {n_threads} threads for parallel processing")
    logger.info(f"Processing {len(chunks)} variant chunks")
    logger.debug(f"Spawning {n_threads} threads for parallel processing")

    with Progress(TextColumn("[bold blue]{task.description}"),
                 BarColumn(complete_style="green"),
                 TaskProgressColumn(),
                 TimeElapsedColumn(),
                 TimeRemainingColumn()) as progress, \
         ThreadPoolExecutor(max_workers=n_threads) as executor:

        annotation_procs = []
        
        # Use BED names from config (validation already ensures names are available)
        for i, bed_file in enumerate(config.bed_files):
            task = progress.add_task("[cyan]Loading BED files...", total=len(config.bed_files))
            bed_name = config.bed_names[i]
            progress.update(task, description=f"Loading BED file: {bed_file.name} ({bed_name})")
            logger.info(f"Loading BED file: {bed_file.name} as '{bed_name}'")
            annotation_procs.append(BEDProcessor(bed_file, source_name=bed_name, rebuild=config.rebuild))
            progress.update(task, advance=1)
        
        # Use GFF names from config (validation already ensures names are available)
        for i, gff_file in enumerate(config.gff_files):
            task = progress.add_task("[cyan]Loading GFF files...", total=len(config.gff_files))
            gff_name = config.gff_names[i]
            progress.update(task, description=f"Loading GFF file: {gff_file.name} ({gff_name})")
            logger.info(f"Loading GFF file: {gff_file.name} as '{gff_name}'")
            annotation_procs.append(GFFProcessor(gff_file, source_name=gff_name, rebuild=config.rebuild))
            progress.update(task, advance=1)
            
        # Use TSV names from config (validation already ensures names are available)
        for i, tsv_file in enumerate(config.tsv_files):
            task = progress.add_task("[cyan]Loading TSV files...", total=len(config.tsv_files))
            tsv_name = config.tsv_names[i]
            progress.update(task, description=f"Loading TSV file: {tsv_file.name} ({tsv_name})")
            logger.info(f"Loading TSV file: {tsv_file.name} as '{tsv_name}'")
            
            # Generic TSV processor with first row as header
            annotation_procs.append(TSVProcessor(tsv_file, source_name=tsv_name, rebuild=config.rebuild))
                
            progress.update(task, advance=1)

    with Progress(TextColumn("[bold blue]{task.description}"),
                 BarColumn(complete_style="green"),
                 TaskProgressColumn(),
                 TimeElapsedColumn(),
                 TimeRemainingColumn()) as progress, \
         ThreadPoolExecutor(max_workers=n_threads) as executor:
        
        # Open all annotation processors using context managers (if any)
        # repeatmasker processor is opened within each batch
        for proc in annotation_procs:
            proc.__enter__()
        
        try:
            task = progress.add_task("[cyan]Annotating variants...", total=len(variants))           

            futures = [
                executor.submit(
                    process_annotate_chunk, 
                    chunk, 
                    annotation_procs,
                    idx,
                    config,
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
                    logger.debug(f"Progress: {completed_variants}/{len(variants)} variants completed ({completed_variants/len(variants)*100:.1f}%)")
                except Exception as e:
                    logger.error(f"Error processing chunk: {str(e)}")
                    raise
        finally:
            logger.debug("All threads completed, closing thread pool")
            for proc in annotation_procs:
                proc.__exit__(None, None, None)
    
    elapsed = time.time() - start_time
    logger.info(f"Completed annotation in {elapsed:.2f} seconds")
    if total_annotation_files > 0:
        logger.info(f"Annotated {len(variants)} variants with RepeatMasker and {total_annotation_files} external annotation file(s)")
    else:
        logger.info(f"Annotated {len(variants)} variants with RepeatMasker")
    
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
            
    if config.tsv_files:
        logger.info(f"TSV file(s) ({len(config.tsv_files)}):")
        for i, tsv_file in enumerate(config.tsv_files):
            tsv_name = config.tsv_names[i]
            logger.info(f"  - {tsv_file} (source: '{tsv_name}')")
            
    if total_annotation_files > 0:
        logger.info(f"Mean overlapping features: {sum(len(v.overlapping_features) for v in variants)/len(variants):.1f}") 
        logger.info(f"Mean proximal features: {sum(len(v.proximal_features) for v in variants)/len(variants):.1f}") 
    
    temp_jsons = glob(str(config.output_dir / 'repeatmasker' / '*.chunk_*.json'))
    repeatmasker_json_path = config.output_dir / config.vcf_input.name.replace('.vcf.gz', '.annotated.repeatmasker_output.json')
    logger.info(f"Merging {len(temp_jsons)} temporary JSON files into {repeatmasker_json_path}")
    logger.info(f"Saving RepeatMasker results to {repeatmasker_json_path}")
    RepeatMaskerProcessor.merge_temp_jsons(config.output_dir / 'repeatmasker', repeatmasker_json_path)
    
    output_vcf_path = config.output_dir / config.vcf_input.name.replace('.vcf.gz', '.annotated.vcf.gz')
    write_vcf_with_progress(output_vcf_path, variants, logger)
    
    logger.info(f"Annotation pipeline completed successfully")
    print(f"Wrote {len(variants)} variants to {output_vcf_path}")
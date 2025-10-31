"""Command module for variant annotation functionality."""

from __future__ import annotations

import logging
import time
from glob import glob
from concurrent.futures import ThreadPoolExecutor

from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn, TimeRemainingColumn
from vesper.processors.annotations import GFFProcessor
from vesper.processors.repeatmasker import RepeatMaskerProcessor
from vesper.utils.config import AnnotateConfig
from vesper.utils.common import load_variants, calculate_chunks, setup_output_directory, write_vcf_with_progress

def process_annotate_chunk(variants: list, annotation_procs: list, chunk_idx: int, config: AnnotateConfig, logger: logging.Logger) -> None:
    """Process a chunk of variants through the annotation pipeline.
    
    Args:
        variants: List of variants to process
        annotation_procs: List of GFFProcessor instances - may be empty
        chunk_idx: Index of this chunk for logging
        config: AnnotateConfig instance
        logger: Logger instance
    """
    logger.info(f"Processing annotation chunk {chunk_idx} ({len(variants)} variants)")
    
    # Batch process GFF annotations (if any)
    if annotation_procs:
        logger.debug(f"Starting batch GFF annotation for {len(variants)} variants against {len(annotation_procs)} databases")
        for proc in annotation_procs:
            proc.batch_annotate_variants(variants, config.proximal_span)
    
    # Same variants batch processed by RepeatMasker
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

    total_annotation_files = len(config.gff_files) if config.gff_files else 0
    if total_annotation_files > 0:
        logger.info(f"Using {total_annotation_files} external GFF annotation file(s):")
        for i, gff_file in enumerate(config.gff_files):
            gff_name = config.gff_names[i] if config.gff_names and i < len(config.gff_names) else gff_file.stem
            logger.info(f"  - {gff_file} (source: '{gff_name}')")
    else:
        logger.info("No external annotation files provided")

    setup_output_directory(config.output_dir, logger)
    variants = load_variants(config.vcf_input, logger)
    
    # TODO: set config options for memory, etc.
    n_threads = config.threads
    chunks, chunk_size = calculate_chunks(variants, n_threads)
    
    
    start_time = time.time()
    logger.info("Starting annotation pipeline")
    logger.info(f"Using {n_threads} threads for parallel processing")
    logger.info(f"Processing {len(chunks)} variant chunks")

    with Progress(TextColumn("[bold blue]{task.description}"),
                 BarColumn(complete_style="green"),
                 TaskProgressColumn(),
                 TimeElapsedColumn(),
                 TimeRemainingColumn()) as progress, \
         ThreadPoolExecutor(max_workers=n_threads) as executor:

        annotation_procs = []
        
        if config.gff_files:
            task = progress.add_task("[cyan]Loading annotations...", total=len(config.gff_files))
            for i, gff_file in enumerate(config.gff_files):
                gff_name = config.gff_names[i]
                logger.info(f"Loading annotations: {gff_file.name} as '{gff_name}'")
                annotation_procs.append(GFFProcessor(gff_file, source_name=gff_name, rebuild=config.rebuild))
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
    logger.info(f"Annotated {len(variants)} variants in {elapsed:.2f} seconds")
    
    if config.gff_files:
        logger.info(f"Annotated with file(s) ({len(config.gff_files)}):")
        for i, gff_file in enumerate(config.gff_files):
            gff_name = config.gff_names[i]
            logger.info(f"  - {gff_file} (source: '{gff_name}')")
    
    temp_jsons = glob(str(config.output_dir / 'repeatmasker' / '*.chunk_*.json'))
    repeatmasker_json_path = config.output_dir / config.vcf_input.name.replace('.vcf.gz', '.annotated.repeatmasker_output.json')
    logger.info(f"Merging {len(temp_jsons)} temporary JSON files into {repeatmasker_json_path}")
    logger.info(f"Saving RepeatMasker results to {repeatmasker_json_path}")
    RepeatMaskerProcessor.merge_temp_jsons(config.output_dir / 'repeatmasker', repeatmasker_json_path)
    
    output_vcf_path = config.output_dir / config.vcf_input.name.replace('.vcf.gz', '.annotated.vcf.gz')
    write_vcf_with_progress(output_vcf_path, variants, logger)
    
    logger.info("Annotation pipeline completed successfully")
    print(f"Wrote {len(variants)} variants to {output_vcf_path}")

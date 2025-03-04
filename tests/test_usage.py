#!/usr/bin/env python3
from pathlib import Path
import logging
import argparse
import time
from concurrent.futures import ThreadPoolExecutor
import os
from vesper.processors.vcf import VCFProcessor
from vesper.processors.reads import ReadProcessor
from vesper.processors.annotations import BEDProcessor
from vesper.models.variants import VariantAnalysis
from vesper.utils.config import Config

def display_logo():
    """Display the vesper logo from the utils/logo.txt file."""
    logo_path = Path(__file__).parents[1] / "src/vesper/utils/vesper_logo.txt"
    try:
        with open(logo_path, 'r') as f:
            logo = f.read()
            print(logo)
    except FileNotFoundError:
        # That's weird, but we'll just print some text instead
        print("vesper - Haplotype-aware structural variant calling")

def setup_logging():
    """Configure basic logging."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger('vesper')

def process_variant_chunk(variants: list, read_proc: ReadProcessor, bed_proc: BEDProcessor, 
                         chunk_idx: int, logger: logging.Logger) -> None:
    """Process a chunk of variants through the full pipeline.
    
    Args:
        variants: List of variants to process
        read_proc: ReadProcessor instance
        bed_proc: BEDProcessor instance
        chunk_idx: Index of this chunk for logging
        logger: Logger instance
    """
    logger.info(f"Processing chunk {chunk_idx} ({len(variants)} variants)")
    
    for variant in variants:
        # Step 1: Get read groups
        support_reads, nonsupport_reads = read_proc.get_read_groups(variant.variant)
        variant.support_reads = support_reads
        variant.nonsupport_reads = nonsupport_reads
        
        # Step 2: Annotate with genomic features
        bed_proc._annotate_variant(variant, proximal_span=500)
        
        # Step 3: Calculate metrics and confidence
        variant._calculate_grouped_metrics()
        variant._calculate_confidence()
    
    # Save registry after each chunk (if modified)
    read_proc.save_registry()
    logger.info(f"Completed chunk {chunk_idx}")

def process_variant_chunk(variants: list, read_proc: ReadProcessor, bed_proc: BEDProcessor, 
                         chunk_idx: int, logger: logging.Logger) -> None:
    """Process a chunk of variants through the full pipeline.
    
    Args:
        variants: List of variants to process
        read_proc: ReadProcessor instance
        bed_proc: BEDProcessor instance
        chunk_idx: Index of this chunk for logging
        logger: Logger instance
    """
    logger.info(f"Processing chunk {chunk_idx} ({len(variants)} variants)")
    
    for variant in variants:
        # Step 1: Get read groups
        support_reads, nonsupport_reads = read_proc.get_read_groups(variant.variant)
        variant.support_reads = support_reads
        variant.nonsupport_reads = nonsupport_reads
        
        # Step 2: Annotate with genomic features
        bed_proc._annotate_variant(variant, proximal_span=500)
        
        # Step 3: Calculate metrics and confidence
        variant._calculate_grouped_metrics()
        variant._calculate_confidence()
    
    # Save registry after each chunk (if modified)
    read_proc.save_registry()
    logger.info(f"Completed chunk {chunk_idx}")


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Vesper variant analysis pipeline")
    parser.add_argument("--vcf", "-v", required=True, help="Input VCF file (required)")
    parser.add_argument("--bam", "-b", required=True, help="Input BAM file (required)")
    parser.add_argument("--output-dir", "-o", required=True, help="Output directory (required)")
    parser.add_argument("--auto-load-registry", choices=['True', 'False'], default='True',
                       help="Whether to automatically load existing registry if found (default: True)")
    parser.add_argument("--force-new-registry", action="store_true",
                       help="Force creation of new registry even if one exists")
    parser.add_argument("--min-support", type=int, default=1,
                       help="Minimum supporting reads")
    parser.add_argument("--max-af", type=float, default=0.1,
                       help="Maximum allele frequency")
    parser.add_argument("--logging", help="Log directory")
    parser.add_argument("--debug", action="store_true",
                       help="Enable debug logging")
    return parser.parse_args()


def main():
    """Main entry point demonstrating the full variant analysis pipeline."""
    display_logo()
    args = parse_args()
    logger = setup_logging()
    
    # Create config from arguments
    config = Config.from_args(args)

    # This is the only BED for now, just hard code it
    bed_path = Path("annotations/hg38/GRCH38_repeatmasker.bed")

    # Create output directory if it doesn't exist
    if not config.output_dir.exists():
        config.output_dir.mkdir(parents=True)
        logger.info(f"Created output directory: {config.output_dir}")
    
    # Step 1: Load variants from VCF
    logger.info(f"Loading variants from {config.vcf_input}")
    variants = []
    with VCFProcessor(config.vcf_input) as vcf_proc:
        variants = list(vcf_proc.instantiate_variants(test_mode=True))
    logger.info(f"Loaded {len(variants)} variants")
    
    variants = variants[:100]  # analyze more variants since we're parallel
    logger.info(f"Continuing with a test set of {len(variants)} variants")

    # Determine optimal chunk size and number of workers
    n_workers = min(32, max(4, os.cpu_count() * 2))
    chunk_size = max(1, len(variants) // (n_workers * 2))
    chunks = [variants[i:i + chunk_size] for i in range(0, len(variants), chunk_size)]
    
    logger.info(f"Processing {len(chunks)} chunks with {n_workers} workers")
    start_time = time.time()
    
    # Process variants in parallel through the full pipeline
    with ReadProcessor(config.bam_file, registry_dir=config.output_dir,
                      auto_load_registry=config.auto_load_registry,
                      force_new_registry=config.force_new_registry) as read_proc, \
         BEDProcessor(bed_path) as bed_proc, \
         ThreadPoolExecutor(max_workers=n_workers) as executor:
        
        # Submit all chunks for processing
        futures = [
            executor.submit(
                process_variant_chunk, 
                chunk, 
                read_proc, 
                bed_proc,
                idx,
                logger
            )
            for idx, chunk in enumerate(chunks)
        ]
        
        # Wait for all chunks to complete
        completed = 0
        for future in futures:
            try:
                future.result()  # This will raise any exceptions from the threads
                completed += 1
                # Force save every 25% of chunks
                if completed % max(1, len(chunks) // 4) == 0:
                    read_proc.save_registry(force=True)
                    logger.info(f"Progress: {completed}/{len(chunks)} chunks completed ({completed/len(chunks)*100:.1f}%)")
            except Exception as e:
                logger.error(f"Error processing chunk: {str(e)}")
                raise
    
    elapsed = time.time() - start_time
    logger.info(f"Completed all processing in {elapsed:.2f} seconds")

    # Placeholder for writing output VCF

    # Print summary statistics
    confident_variants = [v for v in variants if v.confidence > 0.6]
    logger.info(f"\nAnalysis complete:")
    logger.info(f"Found {len(confident_variants)} high-confidence variants")
    logger.info(f"Mean confidence score: {sum(v.confidence for v in variants)/len(variants):.3f}")
    
    # Log some detailed stats for a few variants
    logger.info("Displaying stats for first 5 variants:")
    for variant in variants[:5]:
        logger.info(f"\nVariant {variant.variant.ID}:")
        logger.info(f"  Support reads: {len(variant.support_reads)}")
        logger.info(f"  Non-support reads: {len(variant.nonsupport_reads)}")
        logger.info(f"  Mean support mapQ: {variant.support_reads.mean_mapq:.1f}")
        logger.info(f"  Support soft-clip stats: {variant.support_reads.soft_clip_stats}")
        logger.info(f"  Overlapping features: {len(variant.overlapping_features)}")
        logger.info(f"  Confidence score: {variant.confidence:.3f}")


if __name__ == "__main__":
    main()

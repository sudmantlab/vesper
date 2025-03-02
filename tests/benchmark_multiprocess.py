
"""Stand-alone benchmark script for multiprocessing implementation that properly handles entry points."""

from pathlib import Path
import time
import multiprocessing
import logging
import argparse

from vesper.processors.variants import *
from vesper.processors.annotations import BEDProcessor
from vesper.processors.annotations_multiprocess import MultiprocessBEDProcessor
from vesper.models.variants import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Benchmark multiprocessing implementation for variant annotation')
    parser.add_argument('--batch-size', type=int,
                       help='Batch size for processing variants (suggested: 100)')
    parser.add_argument('--n-workers', type=int, default=multiprocessing.cpu_count() - 1,
                       help='Number of worker processes to spawn (default: cpu count - 1)')
    args = parser.parse_args()

    if not args.batch_size:
        print("Error: --batch-size is required")
        parser.print_help()
        exit(1)

    logging.basicConfig(level=logging.INFO)
    vcf_path = Path("tests/files/hg38/894.duplomap.vcf.gz")
    repeatmasker_path = Path("annotations/hg38/GRCH38_repeatmasker.bed")

    print(f"Loading variants from {vcf_path}...")
    with VCFProcessor(vcf_path) as vcf_proc:
        variants = list(vcf_proc.instantiate_variants(test_mode=False))
    print(f"Loaded {len(variants)} variants")

    variants_original = variants
    variants_multiprocess = [VariantAnalysis(variant=v.variant) for v in variants] # copy

    # Test 1: Original implementation
    print("\n--- TEST 1: ORIGINAL IMPLEMENTATION ---")
    start_time = time.time()
    with BEDProcessor(repeatmasker_path) as ann_proc:
        for variant in variants_original:
            variant.add_annotations(ann_proc)
    end_time = time.time()
    print(f"Original: Annotated in {end_time - start_time:.2f} seconds")

    # Test 2: Multiprocessing implementation
    print("\n--- TEST 2:MULTIPROCESSING IMPLEMENTATION ---")
    start_time = time.time()
    with MultiprocessBEDProcessor(repeatmasker_path) as ann_proc:
        print(f"Using {args.n_workers} workers")
        ann_proc.process_variants_multiprocess(
            variants_multiprocess, 
            proximal_span=500,
            batch_size=args.batch_size,
            n_workers=args.n_workers
        )
    end_time = time.time()
    print(f"Multiprocessing: Annotated in {end_time - start_time:.2f} seconds")

    print("\n--- VERIFYING RESULTS ---")
    total_orig_annotations = sum(len(v.overlapping_features) + len(v.proximal_features) 
                                for v in variants_original)
    total_multiprocess_annotations = sum(len(v.overlapping_features) + len(v.proximal_features) 
                                        for v in variants_multiprocess)

    print(f"Original implementation: {total_orig_annotations} total annotations")
    print(f"Multiprocessing implementation: {total_multiprocess_annotations} total annotations")
    
    if total_orig_annotations == total_multiprocess_annotations:
        print("✓ Results match! Implementations produce identical annotations.")
    else:
        print("⚠ Results differ(?!)")
        print("\nComparison of first 5 variants:")
        for i, (v1, v2) in enumerate(zip(variants_original[:5], variants_multiprocess[:5])):
            if v1.variant.sv_type != SVType.INS:
                continue
            print(f"\nVariant {i+1}:")
            print(f"  Original: {len(v1.overlapping_features)} overlapping, {len(v1.proximal_features)} proximal")
            print(f"  Multiprocessing: {len(v2.overlapping_features)} overlapping, {len(v2.proximal_features)} proximal")

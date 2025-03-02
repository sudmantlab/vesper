
"""Stand-alone benchmark script for the thread-pool based annotation implementation."""

from pathlib import Path
import time
import logging
import argparse

from vesper.processors.variants import *
from vesper.processors.annotations import BEDProcessor
from vesper.processors.annotations_threadpool import ThreadSafeBEDProcessor
from vesper.models.variants import *

parser = argparse.ArgumentParser(description='Benchmark thread pool implementation for variant annotation')
parser.add_argument('--n-workers', type=int,
                   help='Number of worker threads to use (suggested: 16)')
args = parser.parse_args()

if not args.n_workers:
    print("Error: --n-workers is required")
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
variants_threadpool = [VariantAnalysis(variant=v.variant) for v in variants] # copy

# Test 1: Original implementation
print("\n--- TEST 1: ORIGINAL IMPLEMENTATION ---")
start_time = time.time()
with BEDProcessor(repeatmasker_path) as ann_proc:
    for variant in variants_original:
        variant.add_annotations(ann_proc)
end_time = time.time()
print(f"Original: Annotated in {end_time - start_time:.2f} seconds")

# Test 2: Thread Pool implementation
print("\n--- TEST 2: THREAD POOL IMPLEMENTATION ---")
print(f"Using {args.n_workers} workers")
start_time = time.time()
with ThreadSafeBEDProcessor(repeatmasker_path) as ann_proc:
    ann_proc.process_variants_parallel(
        variants_threadpool,
        proximal_span=500,
        n_workers=args.n_workers
    )
end_time = time.time()
print(f"Thread Pool: Annotated in {end_time - start_time:.2f} seconds")

print("\n--- VERIFYING RESULTS ---")
total_orig_annotations = sum(len(v.overlapping_features) + len(v.proximal_features) 
                            for v in variants_original)
total_thread_annotations = sum(len(v.overlapping_features) + len(v.proximal_features) 
                              for v in variants_threadpool)

print(f"Original implementation: {total_orig_annotations} total annotations")
print(f"Thread pool implementation: {total_thread_annotations} total annotations")

if total_orig_annotations == total_thread_annotations:
    print("✓ Results match! Implementations produce identical annotations.")
else:
    print("⚠ Results differ(?!)")
    print("\nComparison of first 5 variants:")
    for i, (v1, v2) in enumerate(zip(variants_original[:5], variants_threadpool[:5])):
        if v1.variant.sv_type != SVType.INS:
            continue
        print(f"\nVariant {i+1}:")
        print(f"  Original: {len(v1.overlapping_features)} overlapping, {len(v1.proximal_features)} proximal")
        print(f"  Thread Pool: {len(v2.overlapping_features)} overlapping, {len(v2.proximal_features)} proximal")

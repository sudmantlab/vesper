# Refinement
`vesper refine` uses aligned PacBio HiFi reads to score each structural variant with simple heuristics derived from mapping quality and soft-clipping metrics. This document outlines the data flow so you can debug or extend the pipeline confidently.

## Overview

1. Load variants from an input VCF (often the output of `vesper annotate`).
2. Instantiate a `ReadProcessor` with the supporting BAM.
3. Split variants into chunks and process them concurrently with a `ThreadPoolExecutor`.
4. For each variant:
   - Fetch reads within ±1 kb (default window) of the variant position.
   - Identify supporting reads (from the `RNAMES` INFO field) versus non-supporting reads.
   - Build aligned read groups in-memory, preferring primary alignments when multiple segments share a read name.
   - Compute aggregate statistics (mean MAPQ, soft-clipping percentages).
   - Derive a confidence score and flags.
5. Compute simple summary stats and write the refined VCF plus tabix index.

Implementation lives in `src/vesper/commands/refine.py`, with supporting models in `src/vesper/models/variants.py` and `src/vesper/processors/reads.py`.

## Input Requirements

- **VCF:** Prefer the annotated VCF (`*.annotated.vcf.gz`). Must include `DR`, `DV`, and ideally `RNAMES` fields to define supporting reads.
- **BAM:** Coordinate-sorted PacBio HiFi alignments aligned to the same reference as the VCF. A `.bai` index must exist; if missing, `refine` attempts to create one via `samtools index`.
- **Output directory:** Directory where the refined VCF and logs will be stored.

## Command Reference

```bash
vesper refine \
  --vcf output/SAMPLE/SAMPLE.annotated.vcf.gz \
  --bam alignments/SAMPLE.hifi.sorted.bam \
  --output-dir output/SAMPLE \
  [--min-support 1] \
  [--max-af 0.1] \
  [--threads 4] \
  [--logging output/logs] \
  [--debug] \
  [--console-output]
```

### Argument Details

- `--vcf/-v` *(required)* – Input VCF (bgzipped, indexed).
- `--bam/-b` *(required)* – BAM with supporting reads.
- `--output-dir/-o` *(required)* – Destination for outputs (refined VCF, logs).
- `--min-support` – Placeholder for future filters; currently recorded in config only.
- `--max-af` – Placeholder for allele frequency; not yet used in scoring.
- `--threads` – Number of worker threads (default 4).
- `--logging`, `--debug`, `--console-output` – Same behavior as in the annotate pipeline.

## Confidence Scoring

Defined in `VariantAnalysis._calculate_confidence`:

1. Initialize `quality_multiplier = 1.0`.
2. Compute mean MAPQ for support vs non-support read groups. If `support_mapq / max(1, nonsupport_mapq) < 0.8`, multiply the confidence by 0.8 and add the `MAPQ_DIFF` flag.
3. Compute soft-clipping percentages. If the absolute difference exceeds 25 percentage points, multiply by 0.75 and add the `SOFTCLIP_DIFF` flag.
4. If all supporting reads are supplementary/secondary, confidence drops to 0 with `NO_PRIMARY_SUPPORT`.
5. If no metrics could be derived (e.g., missing reads), confidence is set to 0 with flag `NO_METRICS`.
6. Confidence ≤ 0.3 triggers a `LOW_CONFIDENCE` filter tag.

Confidence (0.0–1.0) is written to the INFO field, alongside `CONFIDENCE_FLAGS`, read depth metrics (`RMAPQ`, `NSMAPQ`, `RSFTCLIP`, `NSFTCLIP`), and the existing annotation payload.

## Outputs

Inside `--output-dir` expect:

- `logs/<timestamp>.refine.log` – Run log with chunk progress and summary statistics.
- `<input>.refined.vcf.gz` – Refined VCF (bgzipped).
- `<input>.refined.vcf.gz.tbi` – Tabix index.

The console prints a brief summary (`Wrote N variants to ...`) when the run completes.

## Performance Notes

- `--threads` controls how many variant chunks execute in parallel. The BAM fetch uses half of the available CPU cores (`os.cpu_count() // 2`) for faster IO.

## Troubleshooting

- **Missing RNAMES:** The pipeline still executes but supporting read groups will be empty, leading to `NO_METRICS`. Consider enriching the VCF upstream.
- **Index creation failure:** Ensure you have write permissions where the BAM resides. The error message from `pysam.samtools.index` is logged before the run aborts.
- **Unexpected zero confidence:** Inspect `CONFIDENCE_FLAGS` and the numeric metrics in the refined VCF. They correspond directly to the heuristics above.

## Related Resources

- [Annotate pipeline reference](annotate_pipeline.md)
- [Output file reference](outputs.md)

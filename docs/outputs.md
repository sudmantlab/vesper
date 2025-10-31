`vesper` generates a fair number of outputs. This is a non-exhaustive list and description of the generated files.

## Per-Sample Output Directory

After running `vesper annotate` and `vesper refine`, each sample folder under `output/<sample>/` should contain:

- `logs/` – Timestamped log files for each command invocation.
- `repeatmasker/` – Temporary FASTA/JSON files created during annotation (safe to remove post-run).
- `<sample>.annotated.repeatmasker_output.json` – Merged RepeatMasker batches, one JSON array entry per variant with the same schema as `REPEATMASKER_RESULTS` in the VCF.
- `<sample>.annotated.vcf.gz` & `.tbi` – Annotated VCF (pre-refinement).
- `<sample>.annotated.refined.vcf.gz` & `.tbi` – Final VCF with confidence metrics.

## Final VCF Contents

Key INFO fields to be aware of (populated by `VariantAnalysis.to_vcf_record`):

- `CONFIDENCE` – 0–1 confidence score. Derived from read quality comparisons.
- `CONFIDENCE_FLAGS` – Semicolon-separated list (`NONE`, `MAPQ_DIFF`, `SOFTCLIP_DIFF`, `NO_PRIMARY_SUPPORT`, `NO_METRICS`).
- `RMAPQ` / `NSMAPQ` – Mean mapping quality for supporting / non-supporting read groups.
- `RSFTCLIP` / `NSFTCLIP` – Percentage of soft-clipped bases in supporting / non-supporting groups.
- `OVERLAPPING` – JSON array of `GenomicInterval` features overlapping the variant (keys: `source`, `chrom`, `start`, `end`, plus metadata such as `feature_type`, `distance` when proximal).
- `PROXIMAL` – Similar JSON array for nearby features within the `--proximal-span`.
- `REPEATMASKER_RESULTS` – JSON array of `RepeatMaskerResult` dictionaries (fields listed below).

### RepeatMaskerResult Schema

Each object includes:

| Field | Meaning |
|-------|---------|
| `repeat_name` | Name of repeat family (e.g., `L1P1`, `(TA)n`) |
| `repeat_class` | RepeatMasker classification (`LINE/L1`, `Simple_repeat`, etc.) |
| `sw_score` | Smith-Waterman alignment score |
| `divergence`, `deletion`, `insertion` | Percent divergence/deletion/insertion |
| `query_start`, `query_end`, `query_left` | Match coordinates relative to the insertion sequence |
| `strand` | `'+'` or `'C'` (complement) |
| `repeat_start`, `repeat_end`, `repeat_left` | Coordinates within the repeat consensus |
| `match_length`, `match_coverage` | Number of matched bases and fractional coverage |

## Summarized TSV

Produced with `vesper summarize`, which wraps a `bcftools query` call. Purely a convenience function.

```bash
vesper summarize \
  --input output/<sample>/<sample>.annotated.refined.vcf.gz \
  --output output/<sample>/<sample>.summary.tsv
```

The optional `--sample-names` arguments can override the inferred sample IDs. The command streams rows as `bcftools` emits them, so arbitrarily large cohorts are supported.

### Column Definitions

| Column | Description |
|--------|-------------|
| `SAMPLE` | Specimen identifier (folder name in `output/`) |
| `CHROM` | Contig/chromosome (`Variant.CHROM`) |
| `POS` | 1-based position of the variant |
| `ID` | VCF record ID (e.g., `901_Sniffles2.INS.6S0`) |
| `QUAL` | QUAL field copied from the refined VCF |
| `FILTER` | Filter status (`PASS`, `NO_PRIMARY_SUPPORT`, etc.) |
| `CONFIDENCE` | Same as INFO/CONFIDENCE |
| `CONFIDENCE_FLAGS` | Same as INFO/CONFIDENCE_FLAGS (`NONE` when empty) |
| `RSFTCLIP` | Supporting read soft-clipping percentage (`INFO/RSFTCLIP`) |
| `NSFTCLIP` | Non-supporting read soft-clipping percentage (`INFO/NSFTCLIP`) |
| `RMAPQ` | Supporting read mean MAPQ (`INFO/RMAPQ`) |
| `NSMAPQ` | Non-supporting read mean MAPQ (`INFO/NSMAPQ`) |
| `SVLEN` | Structural variant length (`INFO/SVLEN`) |
| `REPEATMASKER_RESULTS` | JSON array of RepeatMasker hits (see schema above). `"."` denotes no hits. |
| `OVERLAPPING` | JSON array of overlapping features. `"."` denotes none. |

The head of the TSV (first 100 rows) includes real data for sample `901`, showcasing typical RepeatMasker JSON payloads, `CONFIDENCE_FLAGS` such as `MAPQ_DIFF`, and cases where filters like `NO_PRIMARY_SUPPORT` appear.
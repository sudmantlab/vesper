# Annotation
`vesper annotate` adds insertion annotation via RepeatMasker and contextual overlap/proximity information from provided annotation files.

## Overview

1. Load variants from a bgzipped VCF (`pysam.VariantFile`).
2. Preload GFF/GTF annotations into thread-safe SQLite caches (`GFFProcessor`) for rapid overlap/proximity queries.
3. Split variants into chunks (`calculate_chunks`) and process them in a thread pool.
4. For each chunk:
   - Batch query overlaps/proximal features for each configured annotation source.
   - Extract insertion sequences and run RepeatMasker in a temporary batch directory.
   - Attach RepeatMasker hits and annotation metadata to the `VariantAnalysis` objects.
5. Write `<sample>.annotated.vcf.gz`, merge RepeatMasker JSON batches, and tabix-index the output.

Refer to `src/vesper/commands/annotate.py` for the orchestration logic.

## Input Requirements

- **VCF:** Must be bgzipped and indexed (`.tbi`/`.csi`). Insertions should carry the full sequence in `ALT`.
- **GFF/GTF (optional):** Any number of external feature sets. Each file is converted to an adjacent `*.sqlite` cache on first use. Use `--rebuild` to force regeneration after updating the source file.
- **RepeatMasker:** CLI must be available (the pipeline calls `RepeatMasker -engine rmblast -gff -species human`).
- **Output directory:** Created if missing. RepeatMasker intermediates and logs live under this folder.

## Command Reference

```bash
vesper annotate \
  --vcf input/sample.vcf.gz \
  --output-dir output/SAMPLE \
  [--files annotations/features1.gff annotations/features2.gff] \
  [--names features1 features2] \
  [--proximal-span 100] \
  [--repeatmasker-n 0] \
  [--threads 4] \
  [--logging output/logs] \
  [--debug] \
  [--console-output] \
  [--rebuild]
```

### Argument Details

- `--vcf/-v` *(required)* – Input structural variant VCF.
- `--output-dir/-o` *(required)* – Destination for all artifacts.
- `--files/-f` – List of GFF/GTF annotation files. Paths can be gzipped or plain text.
- `--names/-n` – Display names matching the order of `--files`. Required if `--files` supplied.
- `--proximal-span` – +/- window (bp) used to report proximal features (default 100).
- `--repeatmasker-n` – Limit the number of RepeatMasker hits per variant. `0` (default) keeps all matches.
- `--threads/-t` – Number of worker threads for chunk processing (default 4).
- `--logging` – Override log directory (defaults to `<output-dir>/logs`).
- `--debug` – Enables verbose logging (root logger set to `DEBUG`).
- `--console-output` – Echo logs to stdout in addition to file.
- `--rebuild` – Rebuilds all GFF SQLite caches before running.

## Annotation Flow

- **Chunking:** `calculate_chunks` derives a chunk size from variant count and threads (minimum 10 records per chunk).
- **GFFProcessor:** Each annotation source is opened as a context manager. Features are indexed into SQLite with range indexes. `batch_annotate_variants` performs `UNION ALL` queries to attach overlapping and proximal features, avoiding duplicates.
- **RepeatMasker:** Each chunk gets its own temp directory (`output/repeatmasker/batch_<uuid>`). Insertions with literal sequences (not `<INS>`, `N`, or `.`) are written to `insertions.fa`. After `RepeatMasker` completes, `.fa.out` is parsed into `RepeatMaskerResult` dataclasses. Temporary JSON summaries (`*.chunk_*.json`) accumulate hits; they are merged at the end of the run.
- **Variant records:** Annotated features are stored in the `VariantAnalysis` object. During VCF writing, these become JSON strings in the `OVERLAPPING`, `PROXIMAL`, and `REPEATMASKER_RESULTS` INFO keys.

## Outputs

Within `--output-dir` you should see:

- `logs/<timestamp>.annotate.log` – Full run log (start/finish timestamps, chunk progress, errors).
- `repeatmasker/` – Temporary FASTA and JSON files per chunk. Safe to delete after confirming the merge.
- `<input>.annotated.repeatmasker_output.json` – Consolidated RepeatMasker output for downstream auditing.
- `<input>.annotated.vcf.gz` – Annotated VCF (bgzipped).
- `<input>.annotated.vcf.gz.tbi` – Tabix index regenerated after writing.

## Performance & Tuning

- `--threads` controls both RepeatMasker batches and GFF queries. Set according to available CPU (RepeatMasker itself is single-threaded per chunk).
- Large GFFs benefit from caching: expect the initial run to spend time building `*.sqlite`. Subsequent runs reuse the cache unless `--rebuild` is set.
- If variants lack `ALT` insertion sequences, RepeatMasker has nothing to annotate. Monitor the log for `Wrote 0 insertion sequences` warnings.

## Troubleshooting

- **Missing RepeatMasker binary:** The command raises `RuntimeError("RepeatMasker executable not found")`; add RepeatMasker to your `PATH`.
- **GFF parsing failures:** The loader validates mandatory columns and coordinate fields. Errors report the line number in the GFF; fix upstream and rerun with `--rebuild`.
- **Annotation duplication:** Proximal features already counted as overlaps are suppressed by design. If you expect duplicates, check the metadata JSON in `OVERLAPPING` and `PROXIMAL`.
- **No variants written:** Verify the VCF path and that the file has a proper `.tbi`/`.csi`. `load_variants` logs the total loaded count.

## Related Resources

- [Refine pipeline details](refine_pipeline.md)
- [Output file reference](outputs.md)
- RepeatMasker documentation: <https://www.repeatmasker.org/>

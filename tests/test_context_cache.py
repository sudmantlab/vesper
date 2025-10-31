from __future__ import annotations

from pathlib import Path

from vesper.models.variants import Variant, SVType, VariantAnalysis
from vesper.models.contexts import VariantContextCacheWriter, VariantContextCacheReader
from vesper.models.repeatmasker import RepeatMaskerResult
from vesper.models.reads import ReadGroup


class FakeRead:
    def __init__(
        self,
        name: str,
        sequence: str,
        strand: str = "+",
        mapq: int = 60,
        edit_distance: int = 0,
        is_secondary: bool = False,
        is_supplementary: bool = False,
    ) -> None:
        self.name = name
        self.sequence = sequence
        self.strand = strand
        self.mapq = mapq
        self.edit_distance = edit_distance
        self.cigar = f"{len(sequence)}M"
        self.cigartuples = [(0, len(sequence))]
        self.cigar_stats = {"M": len(sequence)}
        self.length = len(sequence)
        self.aligned_length = len(sequence)
        self.soft_clip_left = 0
        self.soft_clip_right = 0
        self.is_secondary = is_secondary
        self.is_supplementary = is_supplementary


def make_variant(variant_id: str, repeat_class: str | None = None) -> VariantAnalysis:
    variant = Variant(
        chrom="chr1",
        position=100,
        ID=variant_id,
        ref="A",
        alt="AT",
        qual=50.0,
        filter="PASS",
        info={},
        format="GT",
        samples=[{"DR": 10, "DV": 5}],
        sv_type=SVType.INS,
        sv_length=1,
        DR=10,
        DV=5,
        rnames=["read1"],
    )
    analysis = VariantAnalysis(variant=variant)
    read = FakeRead(name="read1", sequence="ATCGATCGAT")
    analysis.support_reads = ReadGroup([read])  # type: ignore[arg-type]
    analysis.nonsupport_reads = ReadGroup([])  # type: ignore[arg-type]
    analysis.metrics = {"support_mapq": 60}
    analysis.confidence = 0.85
    if repeat_class:
        analysis.repeatmasker_results = [
            RepeatMaskerResult(
                repeat_name="AluYa5",
                repeat_class=repeat_class,
                sw_score=1000,
                divergence=2.0,
                deletion=0.0,
                insertion=0.0,
                query_start=1,
                query_end=10,
                query_left=0,
                strand="+",
                repeat_start=1,
                repeat_end=280,
                repeat_left=0,
                match_length=280,
                match_coverage=1.0,
            )
        ]
    return analysis


def test_cache_write_and_read(tmp_path: Path) -> None:
    writer = VariantContextCacheWriter(tmp_path)
    analysis = make_variant("var1", "SINE/Alu")
    writer.add_variant(sample_id="sample1", analysis=analysis)
    writer.close()

    manifest_parquet = tmp_path / "context_cache.parquet"
    assert manifest_parquet.exists()
    reads_dir = tmp_path / "read_cache"
    assert reads_dir.exists()
    assert any(p.suffixes == [".jsonl", ".gz"] for p in reads_dir.iterdir())

    reader = VariantContextCacheReader(tmp_path)
    contexts = list(reader.iter())
    assert len(contexts) == 1
    ctx = contexts[0]
    assert ctx.sample_id == "sample1"
    assert ctx.variant.ID == "var1"
    assert ctx.metadata["repeat_class"] == "SINE/Alu"
    assert len(ctx.support_reads) == 1
    support = ctx.support_reads[0]
    assert support.sequence == "ATCGATCGAT"
    assert support.cigartuples == [(0, 10)]


def test_cache_filters(tmp_path: Path) -> None:
    writer = VariantContextCacheWriter(tmp_path)
    analysis1 = make_variant("var1", "SINE/Alu")
    analysis2 = make_variant("var2", "LINE/L1")
    analysis2.confidence = 0.4
    writer.add_variant(sample_id="sample1", analysis=analysis1)
    writer.add_variant(sample_id="sample1", analysis=analysis2)
    writer.close()

    reader = VariantContextCacheReader(tmp_path)
    alu_contexts = list(reader.iter(motif="Alu"))
    assert {ctx.variant.ID for ctx in alu_contexts} == {"var1"}

    confident = list(reader.iter(min_confidence=0.8))
    assert {ctx.variant.ID for ctx in confident} == {"var1"}

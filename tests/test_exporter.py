from __future__ import annotations

from pathlib import Path

import pandas as pd

from vesper.analysis.context_extractor import ContextExtractor
from vesper.analysis.exporter import export_breakpoint_tsv, iter_feature_rows, default_columns
from vesper.analysis.features import BASE_FEATURE_SET, ALU_FEATURE_SET
from vesper.models.contexts import BreakpointReadContext, ExtractedSupportRead, VariantContext
from vesper.models.variants import SVType, Variant


def build_variant(variant_id: str = "var1", sv_type: SVType = SVType.INS, sv_length: int = 4) -> Variant:
    return Variant(
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
        sv_type=sv_type,
        sv_length=sv_length,
        DR=10,
        DV=5,
        rnames=["read1"],
    )


def build_support(read_id: str, sequence: str, cigartuples) -> ExtractedSupportRead:
    return ExtractedSupportRead(
        read_id=read_id,
        strand="+",
        mapq=60,
        cigartuples=list(cigartuples),
        cigar="".join(f"{length}{op}" for op, length in [("M", cigartuples[0][1])] if False) if False else None,
        sequence=sequence,
        edit_distance=0,
        soft_clip_left=0,
        soft_clip_right=0,
        cigar_stats=None,
        length=len(sequence),
        aligned_length=len(sequence),
        is_secondary=False,
        is_supplementary=False,
    )


def test_export_breakpoint_tsv_base(tmp_path: Path):
    variant = build_variant()
    support = ExtractedSupportRead(
        read_id="read1",
        strand="+",
        mapq=60,
        cigartuples=[(0, 5), (1, 4), (0, 5)],
        cigar="5M4I5M",
        sequence="AAAATTTTGGGGG",
        edit_distance=0,
        soft_clip_left=0,
        soft_clip_right=0,
        cigar_stats=None,
        length=13,
        aligned_length=10,
        is_secondary=False,
        is_supplementary=False,
    )
    variant_ctx = VariantContext(
        sample_id="sample1",
        variant=variant,
        support_reads=(support,),
        confidence=0.9,
        metrics={},
        metadata={"repeat_query_start": 1, "repeat_query_end": 3},
    )

    extractor = ContextExtractor(context_size=2, strict=False)
    feature_set = BASE_FEATURE_SET
    rows = list(iter_feature_rows([variant_ctx], extractor, feature_set))

    output_path = tmp_path / "base.tsv"
    columns = default_columns(feature_set)
    df = export_breakpoint_tsv(rows, columns, output_path)

    assert output_path.exists()
    loaded = pd.read_csv(output_path, sep="\t")
    assert len(loaded) == 1
    for column in columns:
        assert column in df.columns
    assert loaded.loc[0, "read_id"] == "read1"

    df_dropped = export_breakpoint_tsv(rows, columns, tmp_path / "base_drop.tsv", drop_columns=["left_seq", "right_seq"])
    assert "left_seq" not in df_dropped.columns
    assert "right_seq" not in df_dropped.columns


def test_export_breakpoint_tsv_with_feature_set(tmp_path: Path):
    variant = build_variant(sv_length=4)
    support = ExtractedSupportRead(
        read_id="read1",
        strand="+",
        mapq=60,
        cigartuples=[(0, 5), (1, 4), (0, 5)],
        cigar="5M4I5M",
        sequence="AAAATTTTGGGGG",
        edit_distance=0,
        soft_clip_left=0,
        soft_clip_right=0,
        cigar_stats=None,
        length=13,
        aligned_length=10,
        is_secondary=False,
        is_supplementary=False,
    )
    variant_ctx = VariantContext(
        sample_id="sample1",
        variant=variant,
        support_reads=(support,),
        confidence=0.9,
        metrics={},
        metadata={"repeat_strand": "+", "repeat_class": "SINE/Alu", "repeat_query_start": 1, "repeat_query_end": 3},
    )

    extractor = ContextExtractor(context_size=2, strict=False)
    feature_set = ALU_FEATURE_SET
    rows = list(iter_feature_rows([variant_ctx], extractor, feature_set))

    df = export_breakpoint_tsv(rows, default_columns(feature_set), tmp_path / "alu.tsv")
    assert "tsd_present" in df.columns
    assert "poly_a_present" in df.columns

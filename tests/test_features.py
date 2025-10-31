from __future__ import annotations

from vesper.analysis.features import (
    ALU_FEATURE_SET,
    BASE_FEATURE_SET,
    FeatureSet,
    PolyATailDetector,
    TsdDetector,
    available_feature_sets,
    get_feature_set,
)
from vesper.models.contexts import BreakpointReadContext, VariantContext
from vesper.models.contexts import ExtractedSupportRead
from vesper.models.variants import SVType, Variant


def build_variant(sv_type: SVType = SVType.INS, sv_length: int = 4, alt: str = "AT") -> Variant:
    return Variant(
        chrom="chr1",
        position=100,
        ID="var1",
        ref="A",
        alt=alt,
        qual=50.0,
        filter="PASS",
        info={"SVTYPE": sv_type.name},
        format="GT",
        samples=[{"DR": 10, "DV": 5}],
        sv_type=sv_type,
        sv_length=sv_length,
        DR=10,
        DV=5,
        rnames=["read1"],
    )


def dummy_support() -> ExtractedSupportRead:
    return ExtractedSupportRead(
        read_id="read1",
        strand="+",
        mapq=60,
        cigartuples=[(0, 10)],
        cigar="10M",
        sequence="A" * 10,
        edit_distance=0,
        soft_clip_left=0,
        soft_clip_right=0,
        cigar_stats={"M": 10},
        length=10,
        aligned_length=10,
        is_secondary=False,
        is_supplementary=False,
    )


def test_feature_set_columns_union():
    feature_set = FeatureSet(name="custom", detectors=(TsdDetector(), PolyATailDetector()))
    columns = feature_set.columns()
    assert "tsd_present" in columns
    assert "poly_a_present" in columns
    assert "tsd_edit_distance" in columns
    assert columns[0:5] == [
        "sample_id",
        "variant_id",
        "read_id",
        "chrom",
        "position",
    ]


def test_tsd_detector_identifies_duplication():
    variant = build_variant(sv_length=4)
    variant_ctx = VariantContext(
        sample_id="sample1",
        variant=variant,
        support_reads=(dummy_support(),),
        confidence=0.9,
        metrics={},
        metadata={"repeat_class": "SINE/Alu"},
    )
    context = BreakpointReadContext(
        sample_id="sample1",
        variant_id="var1",
        read_id="read1",
        chrom="chr1",
        position=100,
        sv_type=SVType.INS,
        sv_length=variant.sv_length,
        strand="+",
        support_rank=0,
        ins_start=7,
        ins_end=11,
        left_seq="CCGTTAC",
        insert_seq="GGGG",
        right_seq="TTACAA",
        mapq=60,
        edit_distance=0,
        is_primary=True,
        metadata={},
    )
    detector = TsdDetector(min_len=3, max_len=4)
    result = detector.analyze(context, variant_ctx)
    assert result["tsd_present"] is True
    assert result["tsd_sequence"] == "TTAC"
    assert result["tsd_length"] == 4
    assert result["tsd_left_offset"] == -4
    assert result["tsd_right_offset"] == 0
    assert result["tsd_edit_distance"] == 0


def test_tsd_detector_absent_when_homopolymer():
    variant = build_variant(sv_length=4)
    variant_ctx = VariantContext(
        sample_id="sample1",
        variant=variant,
        support_reads=(dummy_support(),),
        confidence=0.9,
        metrics={},
        metadata={},
    )
    context = BreakpointReadContext(
        sample_id="sample1",
        variant_id="var1",
        read_id="read1",
        chrom="chr1",
        position=100,
        sv_type=SVType.INS,
        sv_length=variant.sv_length,
        strand="+",
        support_rank=0,
        ins_start=5,
        ins_end=8,
        left_seq="AAAAAA",
        insert_seq="GGG",
        right_seq="AAAAAA",
        mapq=60,
        edit_distance=0,
        is_primary=True,
        metadata={},
    )
    detector = TsdDetector(min_len=3, max_len=4)
    result = detector.analyze(context, variant_ctx)
    assert result["tsd_present"] is False
    assert result["tsd_edit_distance"] == 0


def test_tsd_detector_allows_mismatch():
    variant = build_variant(sv_length=6)
    variant_ctx = VariantContext(
        sample_id="sample1",
        variant=variant,
        support_reads=(dummy_support(),),
        confidence=0.9,
        metrics={},
        metadata={"repeat_class": "SINE/Alu"},
    )
    context = BreakpointReadContext(
        sample_id="sample1",
        variant_id="var1",
        read_id="read1",
        chrom="chr1",
        position=100,
        sv_type=SVType.INS,
        sv_length=variant.sv_length,
        strand="+",
        support_rank=0,
        ins_start=7,
        ins_end=13,
        left_seq="GGTACCTA",
        insert_seq="GGGGGG",
        right_seq="TTATCTA",
        mapq=60,
        edit_distance=0,
        is_primary=True,
        metadata={},
    )
    detector = TsdDetector(min_len=4, max_len=7, max_edit_distance=2)
    result = detector.analyze(context, variant_ctx)
    assert result["tsd_present"] is True
    assert result["tsd_edit_distance"] <= 2


def test_poly_a_detector_tracks_orientation():
    variant = build_variant(sv_length=12)
    variant_ctx = VariantContext(
        sample_id="sample1",
        variant=variant,
        support_reads=(dummy_support(),),
        confidence=0.9,
        metrics={},
        metadata={"repeat_strand": "+"},
    )
    context = BreakpointReadContext(
        sample_id="sample1",
        variant_id="var1",
        read_id="read1",
        chrom="chr1",
        position=100,
        sv_type=SVType.INS,
        sv_length=variant.sv_length,
        strand="+",
        support_rank=0,
        ins_start=5,
        ins_end=17,
        left_seq="AC",
        insert_seq="AAAATAAAAAA",
        right_seq="TT",
        mapq=60,
        edit_distance=0,
        is_primary=True,
        metadata={},
    )
    detector = PolyATailDetector(min_total_length=6, max_impurity=0.2)
    result = detector.analyze(context, variant_ctx)
    assert result["poly_a_present"] is True
    assert result["poly_a_length"] >= 6
    assert result["poly_a_sequence"].count("A") >= result["poly_a_length"] * 0.8


def test_poly_a_detector_minus_strand_uses_t():
    variant = build_variant(sv_length=12)
    variant_ctx = VariantContext(
        sample_id="sample1",
        variant=variant,
        support_reads=(dummy_support(),),
        confidence=0.9,
        metrics={},
        metadata={"repeat_strand": "C"},
    )
    context = BreakpointReadContext(
        sample_id="sample1",
        variant_id="var1",
        read_id="read1",
        chrom="chr1",
        position=100,
        sv_type=SVType.INS,
        sv_length=variant.sv_length,
        strand="+",
        support_rank=0,
        ins_start=5,
        ins_end=15,
        left_seq="AC",
        insert_seq="TTTTTTATTTTT",
        right_seq="GG",
        mapq=60,
        edit_distance=0,
        is_primary=True,
        metadata={},
    )
    detector = PolyATailDetector(min_total_length=6, max_impurity=0.2)
    result = detector.analyze(context, variant_ctx)
    assert result["poly_a_present"] is True
    assert set(result["poly_a_sequence"]) <= {"T", "A"}


def test_feature_sets_registry():
    assert BASE_FEATURE_SET.name in available_feature_sets()
    alu = get_feature_set("alu")
    assert alu is ALU_FEATURE_SET

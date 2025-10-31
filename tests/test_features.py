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
from vesper.models.contexts import BreakpointReadContext, VariantContext, ExtractedSupportRead
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


def make_support_read(left: str, insert: str, right: str, strand: str = "+") -> ExtractedSupportRead:
    sequence = "".join([left, insert, right])
    cigartuples = []
    if left:
        cigartuples.append((0, len(left)))
    if insert:
        cigartuples.append((1, len(insert)))
    if right:
        cigartuples.append((0, len(right)))
    return ExtractedSupportRead(
        read_id="read1",
        strand=strand,
        mapq=60,
        cigartuples=cigartuples or None,
        sequence=sequence,
        edit_distance=0,
        soft_clip_left=0,
        soft_clip_right=0,
        is_secondary=False,
        is_supplementary=False,
    )


def make_breakpoint_context(
    variant: Variant,
    left: str,
    insert: str,
    right: str,
    *,
    strand: str = "+",
    ctx_metadata: dict | None = None,
    variant_metadata: dict | None = None,
) -> tuple[BreakpointReadContext, VariantContext]:
    support = make_support_read(left, insert, right, strand=strand)
    context_metadata = {"context_status": "ok", "context_size": 30}
    if ctx_metadata:
        context_metadata.update(ctx_metadata)
    variant_metadata = dict(variant_metadata or {})
    context = BreakpointReadContext(
        sample_id="sample1",
        variant_id=variant.ID,
        read_id=support.read_id,
        chrom=variant.chrom,
        position=variant.position,
        sv_type=variant.sv_type,
        sv_length=variant.sv_length,
        strand=strand,
        support_rank=0,
        ins_start=max(0, len(left) - 1),
        ins_end=len(left) + len(insert),
        left_seq=left or ".",
        insert_seq=insert or ".",
        right_seq=right or ".",
        mapq=support.mapq,
        edit_distance=0,
        is_primary=True,
        metadata=context_metadata,
    )
    variant_ctx = VariantContext(
        sample_id="sample1",
        variant=variant,
        support_reads=(support,),
        confidence=0.9,
        metrics={},
        metadata=variant_metadata,
    )
    return context, variant_ctx


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
    variant = build_variant(sv_length=4, alt="GGGG")
    context, variant_ctx = make_breakpoint_context(
        variant,
        left="CCGTTAC",
        insert="GGGG",
        right="TTACAA",
    )
    detector = TsdDetector(min_len=3, max_len=4)
    result = detector.analyze(context, variant_ctx)
    assert result["tsd_present"] is True
    assert result["tsd_sequence"] == "TTA"
    assert result["tsd_length"] == 3
    assert result["tsd_left_offset"] == -3
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


def test_tsd_detector_requires_exact_match():
    variant = build_variant(sv_length=4, alt="GGGG")
    context, variant_ctx = make_breakpoint_context(
        variant,
        left="CCGTTAC",
        insert="GGGG",
        right="TTGCAA",  # one mismatch relative to left window
    )
    detector = TsdDetector(min_len=3, max_len=4, max_edit_distance=2)
    result = detector.analyze(context, variant_ctx)
    assert result["tsd_present"] is False


def test_poly_a_detector_tracks_orientation():
    variant = build_variant(sv_length=12, alt="AAAAAAAAAAAA")
    context, variant_ctx = make_breakpoint_context(
        variant,
        left="GGGGG",
        insert=variant.alt,
        right="GGGGG",
        ctx_metadata={"context_size": 30},
        variant_metadata={"repeat_strand": "+"},
    )
    detector = PolyATailDetector(min_total_length=6, max_impurity=0.2)
    result = detector.analyze(context, variant_ctx)
    assert result["poly_a_present"] is True
    assert result["poly_a_length"] == len(result["poly_a_sequence"])
    assert result["poly_a_sequence"].count("A") == result["poly_a_length"]


def test_poly_a_detector_minus_strand_uses_t():
    variant = build_variant(sv_length=12, alt="TTTTTTTTTTTT")
    context, variant_ctx = make_breakpoint_context(
        variant,
        left="GGGGG",
        insert=variant.alt,
        right="GGGGG",
        ctx_metadata={"context_size": 30},
        variant_metadata={"repeat_strand": "C"},
    )
    detector = PolyATailDetector(min_total_length=6, max_impurity=0.2)
    result = detector.analyze(context, variant_ctx)
    assert result["poly_a_present"] is True
    assert set(result["poly_a_sequence"]) == {"T"}


def test_feature_sets_registry():
    assert BASE_FEATURE_SET.name in available_feature_sets()
    alu = get_feature_set("alu")
    assert alu is ALU_FEATURE_SET

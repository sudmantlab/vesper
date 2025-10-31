from __future__ import annotations

import pytest

from vesper.analysis.context_extractor import (
    ContextExtractor,
    InsertionNotFound,
    MissingSequenceError,
)
from vesper.models.contexts import ExtractedSupportRead, VariantContext
from vesper.models.variants import SVType, Variant


def build_variant(alt: str, sv_length: int = 3, sv_type: SVType = SVType.INS) -> Variant:
    return Variant(
        chrom="chr1",
        position=100,
        ID="var1",
        ref="A",
        alt=alt,
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


def test_iter_read_contexts_cigar_hit():
    support = ExtractedSupportRead(
        read_id="read1",
        strand="+",
        mapq=60,
        cigartuples=[(0, 5), (1, 3), (0, 5)],
        sequence="AAAAATTTGGGGG",
        edit_distance=2,
        soft_clip_left=0,
        soft_clip_right=0,
        is_secondary=False,
        is_supplementary=False,
    )
    ctx = VariantContext(
        sample_id="sample1",
        variant=build_variant("TTT", sv_length=3),
        support_reads=(support,),
        confidence=0.9,
        metrics={},
        metadata={"repeat_class": "SINE/Alu", "repeat_query_start": 1, "repeat_query_end": 3},
    )
    extractor = ContextExtractor(context_size=2)

    contexts = list(extractor.iter_read_contexts(ctx))
    assert len(contexts) == 1
    result = contexts[0]
    assert result.left_seq == "AAT"
    assert result.insert_seq == "TTT"
    assert result.right_seq == "GG"
    assert result.ins_start == 5
    assert result.ins_end == 8
    assert result.metadata["repeat_class"] == "SINE/Alu"
    assert result.is_primary is True


def test_iter_read_contexts_sequence_alignment():
    support = ExtractedSupportRead(
        read_id="read1",
        strand="+",
        mapq=60,
        cigartuples=[(0, 12)],
        sequence="AAAAGGGAAATTT",
        edit_distance=1,
        soft_clip_left=0,
        soft_clip_right=0,
        is_secondary=True,
        is_supplementary=False,
    )
    ctx = VariantContext(
        sample_id="sample1",
        variant=build_variant("GGG", sv_length=3),
        support_reads=(support,),
        confidence=0.9,
        metrics={},
        metadata={},
    )
    extractor = ContextExtractor(context_size=2)
    contexts = list(extractor.iter_read_contexts(ctx))
    assert len(contexts) == 1
    result = contexts[0]
    assert result.insert_seq == "."
    assert result.left_seq == "."
    assert result.right_seq == "."
    assert result.metadata["context_status"] == "missing_event"
    assert result.is_primary is False


def test_iter_read_contexts_missing_sequence():
    support = ExtractedSupportRead(
        read_id="read1",
        strand="+",
        mapq=60,
        cigartuples=[(0, 5), (1, 3), (0, 5)],
        sequence=None,
    )
    ctx = VariantContext(
        sample_id="sample1",
        variant=build_variant("TTT", sv_length=3),
        support_reads=(support,),
        confidence=0.9,
        metrics={},
        metadata={},
    )
    extractor = ContextExtractor(context_size=2)
    with pytest.raises(MissingSequenceError):
        list(extractor.iter_read_contexts(ctx))


def test_iter_read_contexts_missing_insertion_non_strict():
    support = ExtractedSupportRead(
        read_id="read1",
        strand="+",
        mapq=60,
        cigartuples=[(0, 10)],
        sequence="AAAAAAAAAA",
    )
    ctx = VariantContext(
        sample_id="sample1",
        variant=build_variant("TTT", sv_length=3),
        support_reads=(support,),
        confidence=0.9,
        metrics={},
        metadata={},
    )
    extractor = ContextExtractor(context_size=2, strict=False)
    contexts = list(extractor.iter_read_contexts(ctx))
    assert len(contexts) == 1
    placeholder = contexts[0]
    assert placeholder.insert_seq == "."
    assert placeholder.metadata["context_status"] == "missing_event"

    extractor_strict = ContextExtractor(context_size=2, strict=True)
    with pytest.raises(InsertionNotFound):
        list(extractor_strict.iter_read_contexts(ctx))


def test_iter_read_contexts_deletion():
    support = ExtractedSupportRead(
        read_id="read1",
        strand="+",
        mapq=60,
        cigartuples=[(0, 5), (2, 4), (0, 5)],
        sequence="AAAACCCCTTTT",
    )
    variant = build_variant("<DEL>", sv_length=4, sv_type=SVType.DEL)
    ctx = VariantContext(
        sample_id="sample1",
        variant=variant,
        support_reads=(support,),
        confidence=0.9,
        metrics={},
        metadata={},
    )
    extractor = ContextExtractor(context_size=2)
    contexts = list(extractor.iter_read_contexts(ctx))
    assert len(contexts) == 1
    result = contexts[0]
    assert result.insert_seq == "."
    assert result.left_seq == "AC"
    assert result.right_seq == "CC"
    assert result.metadata["context_status"] == "ok"
    assert result.metadata["event_length"] == 4

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from .variants import Variant, SVType
from .reads import AlignedRead


@dataclass
class ExtractedSupportRead:
    """Lightweight snapshot of an aligned read used for breakpoint analysis."""

    read_id: str
    strand: str
    mapq: int
    cigartuples: Optional[List[Tuple[int, int]]] = None
    cigar: Optional[str] = None
    sequence: Optional[str] = None
    edit_distance: Optional[int] = None
    soft_clip_left: Optional[int] = None
    soft_clip_right: Optional[int] = None
    cigar_stats: Optional[Dict[str, int]] = None
    length: Optional[int] = None
    aligned_length: Optional[int] = None
    is_secondary: bool = False
    is_supplementary: bool = False

    @classmethod
    def from_aligned_read(cls, read: AlignedRead, include_sequence: bool = True) -> "ExtractedSupportRead":
        return cls(
            read_id=read.name,
            strand=read.strand,
            mapq=read.mapq,
            cigartuples=list(read.cigartuples) if read.cigartuples is not None else None,
            cigar=read.cigar,
            sequence=read.sequence if include_sequence else None,
            edit_distance=read.edit_distance,
            soft_clip_left=read.soft_clip_left,
            soft_clip_right=read.soft_clip_right,
            cigar_stats=dict(read.cigar_stats) if read.cigar_stats is not None else None,
            length=read.length,
            aligned_length=read.aligned_length,
            is_secondary=read.is_secondary,
            is_supplementary=read.is_supplementary,
        )


@dataclass
class VariantContext:
    """Variant plus supporting read snapshot used for breakpoint analysis."""

    sample_id: str
    variant: Variant
    support_reads: Tuple[ExtractedSupportRead, ...]
    confidence: Optional[float] = None
    metrics: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_variant_analysis(
        cls,
        sample_id: str,
        analysis: "VariantAnalysis",
        include_sequence: bool = True,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> "VariantContext":
        support_reads = tuple(
            ExtractedSupportRead.from_aligned_read(read, include_sequence=include_sequence)
            for read in analysis.support_reads.reads
        )
        ctx_metadata: Dict[str, Any] = dict(metadata or {})
        if getattr(analysis, "repeatmasker_results", None):
            if analysis.repeatmasker_results:
                primary_repeat = analysis.repeatmasker_results[0]
                ctx_metadata.setdefault("repeat_class", primary_repeat.repeat_class)
                ctx_metadata.setdefault("repeat_strand", primary_repeat.strand)
                ctx_metadata.setdefault("repeat_query_start", primary_repeat.query_start)
                ctx_metadata.setdefault("repeat_query_end", primary_repeat.query_end)
        return cls(
            sample_id=sample_id,
            variant=analysis.variant,
            support_reads=support_reads,
            confidence=analysis.confidence,
            metrics=dict(analysis.metrics),
            metadata=ctx_metadata,
        )


@dataclass
class BreakpointReadContext:
    sample_id: str
    variant_id: str
    read_id: str
    chrom: str
    position: int
    sv_type: SVType
    sv_length: int
    strand: str
    support_rank: int
    ins_start: Optional[int]
    ins_end: Optional[int]
    left_seq: str
    insert_seq: str
    right_seq: str
    mapq: int
    edit_distance: Optional[int] = None
    is_primary: bool = True
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        base = {
            "sample_id": self.sample_id,
            "variant_id": self.variant_id,
            "read_id": self.read_id,
            "chrom": self.chrom,
            "position": self.position,
            "sv_type": self.sv_type.name,
            "sv_length": self.sv_length,
            "strand": self.strand,
            "support_rank": self.support_rank,
            "ins_start": self.ins_start,
            "ins_end": self.ins_end,
            "left_seq": self.left_seq,
            "insert_seq": self.insert_seq,
            "right_seq": self.right_seq,
            "mapq": self.mapq,
            "edit_distance": self.edit_distance,
            "is_primary": self.is_primary,
        }
        base.update(self.metadata)
        return base

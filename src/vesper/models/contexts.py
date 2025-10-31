from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple
from pathlib import Path
from datetime import datetime, timezone
import gzip
import json

import pandas as pd

from .variants import Variant, SVType
from .reads import ReadGroup, AlignedRead


@dataclass
class ExtractedSupportRead:
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

    def to_dict(self) -> Dict[str, Any]:
        return {
            "read_id": self.read_id,
            "strand": self.strand,
            "mapq": self.mapq,
            "cigartuples": self.cigartuples,
            "cigar": self.cigar,
            "sequence": self.sequence,
            "edit_distance": self.edit_distance,
            "soft_clip_left": self.soft_clip_left,
            "soft_clip_right": self.soft_clip_right,
            "cigar_stats": self.cigar_stats,
            "length": self.length,
            "aligned_length": self.aligned_length,
            "is_secondary": self.is_secondary,
            "is_supplementary": self.is_supplementary,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ExtractedSupportRead":
        return cls(
            read_id=data["read_id"],
            strand=data["strand"],
            mapq=data["mapq"],
            cigartuples=[tuple(pair) for pair in data.get("cigartuples", [])] if data.get("cigartuples") is not None else None,
            cigar=data.get("cigar"),
            sequence=data.get("sequence"),
            edit_distance=data.get("edit_distance"),
            soft_clip_left=data.get("soft_clip_left"),
            soft_clip_right=data.get("soft_clip_right"),
            cigar_stats=dict(data["cigar_stats"]) if data.get("cigar_stats") is not None else None,
            length=data.get("length"),
            aligned_length=data.get("aligned_length"),
            is_secondary=data.get("is_secondary", False),
            is_supplementary=data.get("is_supplementary", False),
        )


def _serialize_variant(variant: Variant) -> Dict[str, Any]:
    return {
        "chrom": variant.chrom,
        "position": variant.position,
        "ID": variant.ID,
        "ref": variant.ref,
        "alt": variant.alt,
        "qual": variant.qual,
        "filter": variant.filter,
        "info": variant.info,
        "format": variant.format,
        "samples": variant.samples,
        "sv_type": variant.sv_type.name,
        "sv_length": variant.sv_length,
        "DR": variant.DR,
        "DV": variant.DV,
        "rnames": variant.rnames,
    }


def _deserialize_variant(data: Mapping[str, Any]) -> Variant:
    return Variant(
        chrom=data["chrom"],
        position=data["position"],
        ID=data["ID"],
        ref=data["ref"],
        alt=data["alt"],
        qual=data["qual"],
        filter=data["filter"],
        info=data["info"],
        format=data["format"],
        samples=list(data["samples"]),
        sv_type=SVType[data["sv_type"]],
        sv_length=data["sv_length"],
        DR=data["DR"],
        DV=data["DV"],
        rnames=list(data["rnames"]),
    )


@dataclass
class VariantContext:
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

    def to_manifest_record(self, support_path: Path) -> Dict[str, Any]:
        primary_repeat = self.metadata.get("repeat_class")
        return {
            "sample_id": self.sample_id,
            "variant_id": self.variant.ID,
            "chrom": self.variant.chrom,
            "position": self.variant.position,
            "sv_type": self.variant.sv_type.name,
            "sv_length": self.variant.sv_length,
            "confidence": self.confidence,
            "repeat_class": primary_repeat,
            "support_path": support_path.as_posix(),
            "n_support": len(self.support_reads),
            "created_at": datetime.now(timezone.utc).isoformat(),
        }

    def to_payload(self) -> Dict[str, Any]:
        return {
            "variant": _serialize_variant(self.variant),
            "support_reads": [s.to_dict() for s in self.support_reads],
            "confidence": self.confidence,
            "metrics": self.metrics,
            "metadata": self.metadata,
        }

    @classmethod
    def from_payload(cls, sample_id: str, payload: Mapping[str, Any]) -> "VariantContext":
        variant = _deserialize_variant(payload["variant"])
        support_reads = tuple(ExtractedSupportRead.from_dict(d) for d in payload["support_reads"])
        return cls(
            sample_id=sample_id,
            variant=variant,
            support_reads=support_reads,
            confidence=payload.get("confidence"),
            metrics=dict(payload.get("metrics", {})),
            metadata=dict(payload.get("metadata", {})),
        )


class VariantContextCacheWriter:
    MANIFEST_NAME = "context_cache.parquet"
    READS_DIR = "read_cache"

    def __init__(self, cache_dir: Path, include_sequence: bool = True):
        self.cache_dir = Path(cache_dir)
        self.include_sequence = include_sequence
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.reads_dir = self.cache_dir / self.READS_DIR
        self.reads_dir.mkdir(exist_ok=True)
        self._records: List[Dict[str, Any]] = []

    def add_variant(self, sample_id: str, analysis: "VariantAnalysis", metadata: Optional[Dict[str, Any]] = None) -> VariantContext:
        context = VariantContext.from_variant_analysis(
            sample_id=sample_id,
            analysis=analysis,
            include_sequence=self.include_sequence,
            metadata=metadata,
        )
        support_path = self._write_support_reads(context)
        record = context.to_manifest_record(support_path=support_path)
        self._records.append(record)
        return context

    def _write_support_reads(self, context: VariantContext) -> Path:
        filename = f"{context.sample_id}__{context.variant.ID}.jsonl.gz"
        path = self.reads_dir / filename
        payload = context.to_payload()
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            json.dump(payload, handle)
        return path.relative_to(self.cache_dir)

    def close(self) -> None:
        if not self._records:
            return
        manifest_path = self.cache_dir / self.MANIFEST_NAME
        df = pd.DataFrame.from_records(self._records)
        df.to_parquet(manifest_path, index=False)


class VariantContextCacheReader:
    MANIFEST_NAME = VariantContextCacheWriter.MANIFEST_NAME
    READS_DIR = VariantContextCacheWriter.READS_DIR

    def __init__(self, cache_dir: Path):
        self.cache_dir = Path(cache_dir)
        manifest_path = self.cache_dir / self.MANIFEST_NAME
        if not manifest_path.exists():
            raise FileNotFoundError(f"Manifest not found in {self.cache_dir}")
        self._manifest = pd.read_parquet(manifest_path)

    def iter(
        self,
        sample: Optional[str] = None,
        motif: Optional[str] = None,
        min_confidence: Optional[float] = None,
    ) -> Iterable[VariantContext]:
        df = self._manifest
        if sample is not None:
            df = df[df["sample_id"] == sample]
        if motif is not None and "repeat_class" in df.columns:
            df = df[df["repeat_class"].fillna("").str.contains(motif, na=False)]
        if min_confidence is not None and "confidence" in df.columns:
            df = df[df["confidence"].fillna(0) >= min_confidence]
        for _, row in df.iterrows():
            yield self._load_context(row)

    def _load_context(self, row: pd.Series) -> VariantContext:
        support_path = self.cache_dir / row["support_path"]
        with gzip.open(support_path, "rt", encoding="utf-8") as handle:
            payload = json.load(handle)
        context = VariantContext.from_payload(sample_id=row["sample_id"], payload=payload)
        return context


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

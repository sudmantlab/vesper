from __future__ import annotations

from pydantic import BaseModel, Field
from typing import List, Dict, Optional, Tuple


class ReadMetadata(BaseModel):
    """Metadata about an AlignedRead stored in the read metadata registry."""
    chrom: str
    start: int
    end: int
    mapq: int
    strand: str
    is_supplementary: bool
    is_secondary: bool
    edit_distance: int
    soft_clip_left: int
    soft_clip_right: int
    cigar_stats: Dict[str, int]

    # Optional additional fields to serialize
    cigartuples: Optional[List[Tuple[int, int]]] = None
    cigar: Optional[str] = None
    edit_distance: Optional[int] = None
    methylation_status: Optional[str] = None
    
    class Config:
        frozen = True

class VariantReadGroups(BaseModel):
    """ReadGroups (composed of AlignedReads) associated with a variant."""
    support: List[str] = Field(default_factory=list)
    nonsupport: List[str] = Field(default_factory=list)
    
    class Config:
        frozen = True

class RegistryMetadata(BaseModel):
    """Metadata about the entire read registry."""
    total_variants: int
    total_reads: int
    bam_path: str
    window_size: int
    processed_time: Optional[str] = None  # ISO format timestamp, serialized as string not super important to validate

    class Config:
        frozen = True 

class RepeatMasker(BaseModel):
    """Result of RepeatMasker analysis."""
    repeat_name: str
    repeat_class: str
    score: int
    divergence: float
    deletion: float
    insertion: float
    query_start: int
    query_end: int
    query_left: int
    strand: str
    repeat_start: int
    repeat_end: int
    repeat_left: int
    match_length: int
    match_coverage: float

    class Config:
        frozen = True

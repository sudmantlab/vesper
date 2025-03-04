from dataclasses import dataclass, field
from typing import List, Optional, ClassVar, Pattern, Dict, Any
import re
import numpy as np
from enum import Enum, auto

import logging

from ..processors.annotations import GenomicInterval
from .reads import AlignedRead, ReadGroup

class SVType(Enum):
    """Structural variant type classifications for downstream logic."""
    DEL = auto()  # deletion
    INS = auto()  # insertion
    DUP = auto()  # duplication 
    INV = auto()  # inversion
    BND = auto()  # breakend

@dataclass(frozen=True)
class Variant:
    """Immutable representation of a structural variant call from a VCF record."""
    ID: str
    contig: str
    position: int
    sv_type: SVType
    sv_length: int
    DR: int # DR field
    DV: int # DV field, may not match rnames
    rnames: List[str] = field(default_factory=list)
    quality: int = 0
    
    READ_PATTERN: ClassVar[Pattern] = re.compile(r'([^,;]+)(?:[,;]|$)')

    @classmethod
    def from_pysam_record(cls, record) -> "Variant":
        """Construct variant from pysam VariantRecord object.
        
        Args:
            record: pysam.VariantRecord object
        
        Returns:
            Variant
        """
        try:
            try:
                sv_type = SVType[record.info['SVTYPE']]
            except KeyError:
                raise ValueError(f"Unknown SVTYPE: {record.info['SVTYPE']}")
                
            # Calculate SV length
            sv_length = abs(record.info.get('SVLEN', 0))
            if sv_length == 0 and record.stop:
                sv_length = abs(record.stop - record.start)
            
            # Get supporting read names if present
            rnames = []
            if 'RNAMES' in record.info:
                rnames = record.info['RNAMES'][0].split(',')
                
            return cls(
                ID=record.id,
                contig=record.chrom,
                position=record.pos,
                sv_type=sv_type,
                sv_length=sv_length,
                DR=record.samples[0]['DR'],
                DV=record.samples[0]['DV'],
                rnames=rnames,
                quality=record.qual if record.qual is not None else 0
                )
        except (IndexError, KeyError, ValueError) as e:
            raise ValueError(f"Invalid VCF record: {str(e)}")

    def __post_init__(self):
        if self.position < 0:
            raise ValueError("Position must be non-negative.")
        if self.sv_length < 0:
            raise ValueError("sv_length must be non-negative")
            
    def __repr__(self) -> str:
        return f"Variant(contig={self.contig}, pos={self.position}, type={self.sv_type.name}, len={self.sv_length}, DR={self.DR}, DV={self.DV})"

@dataclass
class VariantAnalysis:
    """Mutable wrapper object for variant analysis and annotation."""
    variant: Variant
    support_reads: ReadGroup = field(default_factory=lambda: ReadGroup([]))
    nonsupport_reads: ReadGroup = field(default_factory=lambda: ReadGroup([]))
    metrics: Dict[str, Any] = field(default_factory=dict)
    overlapping_features: List[GenomicInterval] = field(default_factory=list)
    proximal_features: List[GenomicInterval] = field(default_factory=list)
    confidence: float = 0.0
    logger: logging.Logger = field(default_factory=lambda: logging.getLogger(__name__))

    def __repr__(self) -> str:
        return (f"\n{self.variant},\nconfidence={self.confidence},\nsupport_reads={len(self.support_reads)},\n"
                f"nonsupport_reads={len(self.nonsupport_reads)},\n"
                f"overlapping_features={len(self.overlapping_features)},\n"
                f"proximal_features={len(self.proximal_features)}")

    # This method is no longer used - reads are directly assigned in ReadProcessor.process_variants
    # Keeping for backwards compatibility
    def _add_read_groups(self, processor) -> None:
        """Add supporting and non-supporting read groups.
        
        Args:
            processor: ReadProcessor instance
        """
        support_reads, nonsupport_reads = processor.get_read_groups(self.variant)
        self.support_reads = support_reads
        self.nonsupport_reads = nonsupport_reads
    
    def _calculate_grouped_metrics(self) -> None:
        """Calculate metrics based on read groups."""
        if not self.support_reads.reads or not self.nonsupport_reads.reads:
            self.logger.warning(f"Missing reads for variant {self.variant.ID}, skipping metrics calculation")
            return

        # Get basic read counts
        self.metrics['n_support'] = len(self.support_reads)
        self.metrics['n_nonsupport'] = len(self.nonsupport_reads)
        
        # Get read quality metrics
        self.metrics['support_mapq'] = self.support_reads.mean_mapq
        self.metrics['nonsupport_mapq'] = self.nonsupport_reads.mean_mapq
        
        # Get soft-clipping metrics
        self.metrics['support_softclip'] = self.support_reads.soft_clip_stats
        self.metrics['nonsupport_softclip'] = self.nonsupport_reads.soft_clip_stats

        # Comparison metrics
        self.metrics['comparison'] = self.support_reads.compare_stats(self.nonsupport_reads)
            
    def _calculate_confidence(self) -> None:
        """Basic confidence calculation based on mapq differential and softclipping.
        
        Should be called after _calculate_grouped_metrics() and adding annotations.
        """
        if not self.metrics:
            self.confidence = 0.0
            return
            
        # Calculate mapq differential between support and nonsupport reads
        mapq_diff = self.metrics['comparison']['mapq_mean'] - self.metrics['comparison']['mapq_mean_other']
        
        # Calculate softclip differential between support and nonsupport reads
        support_pct = self.metrics['comparison']['softclip_stats']['pct_softclipped']
        nonsupport_pct = self.metrics['comparison']['softclip_stats_other']['pct_softclipped']
        softclip_diff = support_pct - nonsupport_pct
        
        # scaling
        mapq_factor = min(max(mapq_diff / 60.0, 0), 1)  # max mapq of 60
        softclip_factor = min(max(-softclip_diff / 50.0, 0), 1)  # less softclipping is better
        
        # TODO: fix weighting
        self.confidence = (mapq_factor * 0.7) + (softclip_factor * 0.3)
    

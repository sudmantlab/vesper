from dataclasses import dataclass, field
from typing import List, Optional, ClassVar, Pattern, Dict, Any
import re
import numpy as np
from enum import Enum, auto
import pysam
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
    chrom: str
    position: int
    ID: str
    ref: str
    alt: str
    qual: Optional[float]
    filter: Optional[str]
    info: Dict[str, Any]
    format: str
    samples: List[Dict[str, Any]]
    sv_type: SVType
    sv_length: int
    DR: int # DR field
    DV: int # DV field, may not match rnames
    rnames: List[str] = field(default_factory=list)
    
    READ_PATTERN: ClassVar[Pattern] = re.compile(r'([^,;]+)(?:[,;]|$)')

    @classmethod
    def from_pysam_record(cls, record) -> "Variant":
        """Construct Variant object from pysam VariantRecord object.
        
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
                
            # lead with first 8 required fields
            return cls(
                chrom=record.chrom, # equiv #CHROM
                position=record.pos, # equiv POS
                ID=record.id, 
                ref=record.ref,
                alt=record.alts[0], # assume single ALT
                qual=record.qual if record.qual is not None else 0,
                filter=record.filter,
                info=record.info,
                format=record.format,
                samples=record.samples,
                sv_type=sv_type, # equiv INFO: SVTYPE
                sv_length=sv_length, # equiv INFO: SVLEN
                DR=record.samples[0]['DR'],
                DV=record.samples[0]['DV'],
                rnames=rnames,
                )
        except (IndexError, KeyError, ValueError) as e:
            raise ValueError(f"Invalid VCF record: {str(e)}")

    def __post_init__(self):
        if self.position < 0:
            raise ValueError("Position must be non-negative.")
        if self.sv_length < 0:
            raise ValueError("sv_length must be non-negative")
            
    def __repr__(self) -> str:
        return f"Variant(chrom={self.chrom}, pos={self.position}, type={self.sv_type.name}, len={self.sv_length}, DR={self.DR}, DV={self.DV})"

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
            
        # Calculate mapq diff between support and nonsupport reads as a ratio between support and nonsupport
        # Under the assumption that support reads may have worse alignment quality
        mapq_ratio = max(1, self.metrics['comparison']['mapq_mean'] / self.metrics['comparison']['mapq_mean_other'])
        
        # Calculate softclip diff between support and nonsupport reads as a ratio between support and nonsupport
        # Under the same assumption as above
        support_pct = self.metrics['comparison']['softclip_stats']['pct_softclipped']
        nonsupport_pct = self.metrics['comparison']['softclip_stats_other']['pct_softclipped']
        softclip_ratio = max(1, support_pct / nonsupport_pct)
        
        # TODO: smarter weighting
        mapq_weight, softclip_weight = 0.3, 0.7
        weighted_score = (mapq_ratio ** mapq_weight) + (softclip_ratio ** softclip_weight)
        self.confidence = weighted_score
        
    def to_vcf_record(self) -> pysam.VariantRecord:
        """Convert VariantAnalysis to a pysam.VariantRecord object.
        
        """
        """Convert VariantAnalysis to a VCF line string.
        
        Returns:
            String representation of the variant in VCF format
        """
        chrom = self.variant.chrom
        pos = str(self.variant.position)
        id_field = self.variant.ID
        ref = self.variant.ref
        alt = self.variant.alt
        
        # Format QUAL and FILTER fields
        qual = str(self.variant.qual) if self.variant.qual is not None else "."
        filter_field = self.variant.filter
        
        # Build INFO field
        info = []
        info.append(f"SVTYPE={self.variant.sv_type.name}")
        info.append(f"SVLEN={self.variant.sv_length}")
        
        if self.variant.rnames:
            rnames_str = ",".join(self.variant.rnames)
            info.append(f"RNAMES={rnames_str}")
            
        info.append(f"CONFIDENCE={self.confidence:.2f}")
        
        # Add read group metrics to INFO field if available
        if self.metrics and 'comparison' in self.metrics:
            info.append(f"RMAPQ={self.metrics['comparison']['mapq_mean']:.0f}")
            info.append(f"NSMAPQ={self.metrics['comparison']['mapq_mean_other']:.0f}")
            info.append(f"RSFTCLIP={self.metrics['comparison']['softclip_stats']['pct_softclipped']:.1f}")
            info.append(f"NSFTCLIP={self.metrics['comparison']['softclip_stats_other']['pct_softclipped']:.1f}")
        
        # Add annotations to INFO field if available
        if self.overlapping_features:
            overlapping = []
            for key, interval in self.overlapping_features.items():
                if isinstance(interval, GenomicInterval):
                    flattened = interval.to_dict()
                    overlapping.append(str(flattened))
            if overlapping:
                info.append(f"OVERLAPPING={','.join(overlapping)}")

        if self.proximal_features:
            proximal = []
            for key, interval in self.proximal_features.items():
                if isinstance(interval, GenomicInterval):
                    flattened = interval.to_dict()
                    proximal.append(str(flattened))
            if proximal:
                info.append(f"PROXIMAL={','.join(proximal)}")
                
        info_field = ";".join(info)
        
        # Format FORMAT and sample fields
        format_field = "DR:DV"
        sample_field = f"{self.variant.DR}:{self.variant.DV}"
        
        # Combine all fields into a VCF record
        vcf_record = "\t".join([chrom, pos, id_field, ref, alt, qual, filter_field, info_field, format_field, sample_field])
        
        return vcf_record
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
            
        Raises:
            ValueError: If record is invalid or missing required fields
        """
        try:
            # Handle SVTYPE Enum conversion
            try:
                sv_type = SVType[record.info['SVTYPE']]
            except KeyError:
                raise ValueError(f"Unknown SVTYPE: {record.info['SVTYPE']}")
                
            # Calculate SV length - TODO (refine?): check against insertion seq length, NM edit distance
            sv_length = abs(record.info.get('SVLEN', 0))
            if sv_length == 0 and record.stop:
                sv_length = abs(record.stop - record.start)
            
            # Get supporting read names if present – TODO: handle compatibility with VCFs where this is not present
            rnames = []
            if 'RNAMES' in record.info:
                rnames = str(record.info['RNAMES'][0]).split(',')
            
            # handle filter field conversion
            # TODO: check cases where None in pysam
            filter_list = []
            for f in record.filter:
                if hasattr(f, 'name'):
                    filter_list.append(str(f.name))
                else:
                    filter_list.append(str(f))
            filter_str = ";".join(filter_list) if filter_list else "PASS"
            
            info_dict = {}
            for key, value in record.info.items():
                if isinstance(value, tuple): 
                    info_dict[key] = list(value)
                else:
                    info_dict[key] = value
            
            # TODO: unlikely to be needed given single sample processing, but keep and review later
            sample_dicts = []
            for sample in record.samples.values():
                sample_dict = {}
                for key, value in sample.items():
                    if isinstance(value, tuple):
                        sample_dict[key] = list(value)
                    else:
                        sample_dict[key] = value
                sample_dicts.append(sample_dict)
            
            return cls(
                chrom=str(record.chrom),
                position=int(record.pos),
                ID=str(record.id if record.id is not None else "."),
                ref=str(record.ref),
                alt=str(record.alts[0]),  # assume single ALT
                qual=float(record.qual) if record.qual is not None else 0.0,
                filter=filter_str,
                info=info_dict,
                format=str(record.format),
                samples=sample_dicts,
                sv_type=sv_type,
                sv_length=int(sv_length),
                DR=int(record.samples[0]['DR']),
                DV=int(record.samples[0]['DV']),
                rnames=rnames,
            )
        except (IndexError, KeyError, ValueError, TypeError) as e:
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
    confidence: Optional[float] = None
    logger: logging.Logger = field(default_factory=lambda: logging.getLogger(__name__))

    def __repr__(self) -> str:
        return (f"\n{self.variant},\nconfidence={self.confidence},\nsupport_reads={len(self.support_reads)},\n"
                f"nonsupport_reads={len(self.nonsupport_reads)},\n"
                f"overlapping_features={len(self.overlapping_features)},\n"
                f"proximal_features={len(self.proximal_features)}")
    
    def _calculate_grouped_metrics(self) -> None:
        """Calculate metrics based on read groups."""
        if not self.support_reads.reads or not self.nonsupport_reads.reads:
            self.logger.warning(f"Missing reads for variant {self.variant.ID}, skipping metrics calculation")
            return

        self.metrics['n_support'] = len(self.support_reads)
        self.metrics['n_nonsupport'] = len(self.nonsupport_reads)
        
        self.metrics['support_mapq'] = self.support_reads.mean_mapq
        self.metrics['nonsupport_mapq'] = self.nonsupport_reads.mean_mapq
        
        self.metrics['support_softclip'] = self.support_reads.soft_clip_stats
        self.metrics['nonsupport_softclip'] = self.nonsupport_reads.soft_clip_stats

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
        weighted_score = max(0, mapq_ratio * (1/softclip_ratio))
        self.confidence = weighted_score # TODO: add filtering as next method to mark low confidence variants
        
    def to_vcf_record(self) -> str:
        """Convert VariantAnalysis to a VCF line string.
        
        Returns:
            String representation of the variant in VCF format
        """
        chrom = self.variant.chrom
        pos = str(self.variant.position)
        id_field = self.variant.ID
        ref = self.variant.ref
        alt = self.variant.alt
        
        qual = str(self.variant.qual) if self.variant.qual is not None else "."
        filter_field = self.variant.filter
        
        # Start building INFO field
        info = []
        info.append(f"SVTYPE={self.variant.sv_type.name}")
        info.append(f"SVLEN={self.variant.sv_length}")
        
        if self.variant.rnames:
            rnames_str = ",".join(self.variant.rnames)
            info.append(f"RNAMES={rnames_str}")
            
        # Only add confidence if it has been calculated
        if self.confidence is not None:
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
            for interval in self.overlapping_features:
                if isinstance(interval, GenomicInterval):
                    flattened = interval.to_dict()
                    overlapping.append(str(flattened))
            if overlapping:
                info.append(f"OVERLAPPING={','.join(overlapping)}")

        if self.proximal_features:
            proximal = []
            for interval in self.proximal_features:
                if isinstance(interval, GenomicInterval):
                    flattened = interval.to_dict()
                    proximal.append(str(flattened))
            if proximal:
                info.append(f"PROXIMAL={','.join(proximal)}")
                
        info_field = ";".join(info) 
    
        format_field = "DR:DV"
        sample_field = f"{self.variant.DR}:{self.variant.DV}"
        
        vcf_record = "\t".join([chrom, pos, id_field, ref, alt, qual, filter_field, info_field, format_field, sample_field])
        
        return vcf_record
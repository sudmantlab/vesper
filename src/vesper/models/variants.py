from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, ClassVar, Pattern, Dict, Any
import re
from enum import Enum, auto
import logging
import json

from .interval import GenomicInterval
from .repeatmasker import RepeatMaskerResult
from .reads import ReadGroup

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
                if key in ['OVERLAPPING', 'PROXIMAL', 'REPEATMASKER_RESULTS']:
                    info_dict[key] = ','.join(value) # reconstruct json string
                elif isinstance(value, tuple): 
                    info_dict[key] = list(value)
                else:
                    info_dict[key] = str(value)
            
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
    logger: logging.Logger = field(default_factory=lambda: logging.getLogger(__name__))
    support_reads: ReadGroup = field(default_factory=lambda: ReadGroup([]))
    nonsupport_reads: ReadGroup = field(default_factory=lambda: ReadGroup([]))
    metrics: Dict[str, Any] = field(default_factory=dict)

    overlapping_features: List[GenomicInterval] = field(default_factory=list)
    proximal_features: List[GenomicInterval] = field(default_factory=list)
    repeatmasker_results: List[RepeatMaskerResult] = field(default_factory=list)
    confidence: Optional[float] = None
    confidence_flags: List[str] = field(default_factory=list)
    filter: Optional[str] = None

    def __post_init__(self):
        if 'OVERLAPPING' in self.variant.info:
            self.overlapping_features = [GenomicInterval.from_json(record) for record in json.loads(self.variant.info['OVERLAPPING'])]
        if 'PROXIMAL' in self.variant.info:
            self.proximal_features = [GenomicInterval.from_json(record) for record in json.loads(self.variant.info['PROXIMAL'])]
        if 'REPEATMASKER_RESULTS' in self.variant.info:
            self.repeatmasker_results = [RepeatMaskerResult.from_json(record) for record in json.loads(self.variant.info['REPEATMASKER_RESULTS'])]

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
            
    def _calculate_confidence(self) -> None:
        """Calculate confidence based on read quality metrics only.
        
        Should be called after _calculate_grouped_metrics().
        """
        if not self.metrics:
            self.logger.warning(f"Missing metrics for variant {self.variant.ID}, skipping confidence calculation")
            self.filter = 'NO_METRICS'
            self.confidence = 0.0
            return
        
        # immediately fail if all support reads are secondary/supplementary
        if all(read.is_secondary or read.is_supplementary for read in self.support_reads.reads):
            self.confidence = 0.0
            self.confidence_flags.append('NO_PRIMARY_SUPPORT')
            self.filter = 'NO_PRIMARY_SUPPORT'
            return
            
        quality_multiplier = 1.0
        # Calculate mapq diff between support and nonsupport reads as a ratio between support and nonsupport
        # Under the assumption that support reads may have worse alignment quality
        mapq_ratio = self.metrics['support_mapq'] / max(1, self.metrics['nonsupport_mapq'])
        if mapq_ratio < 0.8:
            quality_multiplier = quality_multiplier * 0.8
            self.confidence_flags.append('MAPQ_DIFF')
        
        # Calculate softclip diff between support and nonsupport reads as absolute difference
        support_pct = self.metrics['support_softclip']['pct_softclipped']
        nonsupport_pct = self.metrics['nonsupport_softclip']['pct_softclipped']
        softclip_diff = abs(support_pct - nonsupport_pct)
        if softclip_diff > 25:
            quality_multiplier = quality_multiplier * 0.75
            self.confidence_flags.append('SOFTCLIP_DIFF')

        # Calculate final confidence score based on read quality only
        self.confidence = quality_multiplier
        if self.confidence <= 0.3 and self.filter is None:
            self.filter = "LOW_CONFIDENCE"
        
    def to_vcf_record(self) -> str:
        """Convert VariantAnalysis to a VCF line string.
        
        Returns:
            String representation of the variant in VCF format
        """
        chrom = self.variant.chrom
        pos = self.variant.position
        id_field = self.variant.ID
        ref = self.variant.ref
        alt = self.variant.alt
        
        qual = self.variant.qual if self.variant.qual is not None else "."
        filter_field = self.filter if self.filter else self.variant.filter
        
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
            qual = round(self.confidence * float(qual), 0) # TODO: this is more or less a placeholder to modify the qual field
            if self.confidence_flags:
                info.append(f"CONFIDENCE_FLAGS={';'.join(self.confidence_flags)}")
            else:
                info.append("CONFIDENCE_FLAGS=NONE")
        
        # Add read group metrics to INFO field if available
        if self.metrics:
            info.append(f"RMAPQ={self.metrics['support_mapq']:.0f}")
            info.append(f"NSMAPQ={self.metrics['nonsupport_mapq']:.0f}")
            info.append(f"RSFTCLIP={self.metrics['support_softclip']['pct_softclipped']:.1f}")
            info.append(f"NSFTCLIP={self.metrics['nonsupport_softclip']['pct_softclipped']:.1f}")
        
        # Add annotations to INFO field if available
        # TODO: Make this a general method for all annotations
        if self.overlapping_features:
            overlapping = []
            for interval in self.overlapping_features:
                if isinstance(interval, GenomicInterval):
                    overlapping.append(interval.to_dict())
            if overlapping:
                info.append(f"OVERLAPPING={json.dumps(overlapping)}")

        if self.proximal_features:
            proximal = []
            for interval in self.proximal_features:
                if isinstance(interval, GenomicInterval):
                    proximal.append(interval.to_dict())
            if proximal:
                info.append(f"PROXIMAL={json.dumps(proximal)}")
        
        if self.repeatmasker_results:
            repeatmasker = []
            for result in self.repeatmasker_results: 
                # TODO: ensure this is compatible with both best/multiple hits
                if isinstance(result, RepeatMaskerResult):
                    repeatmasker.append(result.to_dict())
            if repeatmasker:
                info.append(f"REPEATMASKER_RESULTS={json.dumps(repeatmasker)}")
                
        info_field = ";".join(info) 
    
        format_field = "DR:DV"
        sample_field = f"{self.variant.DR}:{self.variant.DV}"
        
        record_fields = [chrom, pos, id_field, ref, alt, qual, filter_field, info_field, format_field, sample_field]
        vcf_record = "\t".join([str(field) for field in record_fields])
        
        return vcf_record

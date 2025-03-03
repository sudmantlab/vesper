from dataclasses import dataclass, field
from typing import List, Optional, ClassVar, Pattern, Dict
import re
import numpy as np
from enum import Enum, auto

from ..processors.annotations import GenomicInterval
from .reads import AlignedRead

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
    support_reads: List[AlignedRead] = field(default_factory=list)
    nonsupport_reads: List[AlignedRead] = field(default_factory=list)
    support_metrics: Optional[Dict] = None #todo
    nonsupport_metrics: Optional[Dict] = None #todo
    overlapping_features: List[GenomicInterval] = field(default_factory=list)
    proximal_features: List[GenomicInterval] = field(default_factory=list)
    confidence: float = 0.0

    def __repr__(self) -> str:
        return f"\n{self.variant},\nconfidence={self.confidence},\nsupport_reads={len(self.support_reads)},\nnonsupport_reads={len(self.nonsupport_reads)},\noverlapping_features={len(self.overlapping_features)},\nproximal_features={len(self.proximal_features)}"

    def add_locus_reads(self, processor: "ReadProcessor") -> None:
        """Add supporting and non-supporting reads to variant.
        
        Args:
            processor: ReadProcessor instance
        """
        support_reads, nonsupport_reads = processor.get_read_groups(self.variant)
        self.support_reads = support_reads
        self.nonsupport_reads = nonsupport_reads
    
    def calculate_confidence(self):
        """calculate confidence score based on metrics and annotations
        
        TODO: implement actual confidence calculation logic
        """
        self.confidence = 0.0
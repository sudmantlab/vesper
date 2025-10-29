from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple
from functools import cached_property

import pysam

@dataclass
class AlignedRead:
    """Represents an aligned read from a BAM file. Lazy-loaded to minimize memory overhead."""
    _read: pysam.AlignedSegment

    # eagerly loaded attributes
    name: str = field(init=False)
    chrom: str = field(init=False)
    start: int = field(init=False)
    end: int = field(init=False)
    is_supplementary: bool = field(init=False)
    is_secondary: bool = field(init=False)

    def __post_init__(self):
        # load minimal properties we always need
        self.name = self._read.query_name
        self.chrom = self._read.reference_name
        self.start = self._read.reference_start
        self.end = self._read.reference_end
        self.is_supplementary = self._read.is_supplementary
        self.is_secondary = self._read.is_secondary

    @classmethod
    def from_pysam_read(cls, read: pysam.AlignedSegment) -> 'AlignedRead':
        return cls(_read=read)

    @property
    def sequence(self) -> str:
        return self._read.query_sequence

    @property
    def length(self) -> int:
        return self._read.query_length

    @property
    def aligned_length(self) -> int:
        return self._read.reference_length

    @property
    def strand(self) -> str:
        return '-' if self._read.is_reverse else '+'

    @property
    def mapq(self) -> int:
        return self._read.mapping_quality

    @cached_property
    def cigartuples(self) -> List[Tuple[int, int]]:
        return self._read.cigartuples

    @property
    def cigar(self) -> str:
        return self._read.cigarstring

    @property
    def edit_distance(self) -> Optional[int]:
        return self._read.get_tag('NM') if self._read.has_tag('NM') else None

    @property
    def methylation(self) -> Optional[str]:
        return self._read.get_tag('MM') if self._read.has_tag('MM') else None

    @cached_property
    def soft_clip_left(self) -> int:
        return self.cigartuples[0][1] if self.cigartuples[0][0] == 4 else 0

    @cached_property
    def soft_clip_right(self) -> int:
        return self.cigartuples[-1][1] if self.cigartuples[-1][0] == 4 else 0

    @cached_property
    def cigar_stats(self):
        return self._calculate_cigar_stats(self.cigartuples)

    @staticmethod
    def _calculate_cigar_stats(cigartuples: List[Tuple[int, int]]) -> Dict[str, int]:
        op_map = {0:'M', 1:'I', 2:'D', 4:'S', 5:'H', 6:'P', 7:'=', 8:'X', 3:'N'}
        stats = {}
        for op, count in cigartuples:
            op_char = op_map[op]
            stats[op_char] = stats.get(op_char, 0) + count
        return stats

@dataclass
class ReadGroup:
    """A group of aligned reads with methods for calculating aggregate statistics.
    
    Attributes:
        reads: List of AlignedRead objects in this group
    """
    reads: List[AlignedRead]
    
    def __len__(self) -> int:
        """Return the number of reads in this group."""
        return len(self.reads)
    
    @cached_property
    def mean_mapq(self) -> float:
        """Calculate the mean mapping quality across all reads.
        
        Returns:
            Mean mapping quality value
        """
        return sum(r.mapq for r in self.reads) / len(self)
    
    @cached_property
    def soft_clip_stats(self) -> Dict[str, float]:
        """Calculate soft-clipping statistics across all reads.
        
        Returns:
            Dictionary containing:
                - pct_softclipped: Percentage of total bases that are soft-clipped
                - mean_left_clip: Average number of soft-clipped bases on left end
                - mean_right_clip: Average number of soft-clipped bases on right end
        """
        total_len = sum(r.length for r in self.reads)
        total_soft = sum(r.cigar_stats.get('S', 0) for r in self.reads)
        return {
            'pct_softclipped': (total_soft / total_len) * 100,
            'mean_left_clip': sum(r.soft_clip_left for r in self.reads) / len(self),
            'mean_right_clip': sum(r.soft_clip_right for r in self.reads) / len(self)
        }
    
    @cached_property
    def cigar_stats(self) -> Dict[str, int]:
        """Calculate aggregate CIGAR statistics across all reads.
        
        Returns:
            Dictionary mapping CIGAR operations to their total counts
        """
        total_cigar_stats = {}
        for r in self.reads:
            for op, count in r.cigar_stats.items():
                total_cigar_stats[op] = total_cigar_stats.get(op, 0) + count
        return total_cigar_stats

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple
from pathlib import Path

import pysam

@dataclass
class AlignedRead:
    read: pysam.AlignedSegment
    name: str = field(init=False)
    sequence: str = field(init=False)
    length: int = field(init=False)
    aligned_length: int = field(init=False)
    soft_clip_left: int = field(init=False)
    soft_clip_right: int = field(init=False)
    contig: str = field(init=False)
    start: int = field(init=False)
    end: int = field(init=False)
    strand: str = field(init=False)
    mapq: int = field(init=False)
    cigartuples: List[Tuple[int, int]] = field(init=False)
    cigar: str = field(init=False)
    edit_distance: int = field(init=False)
    is_supplementary: bool = field(init=False)
    is_secondary: bool = field(init=False)
    methylation: Optional[str] = field(init=False)
    cigar_stats: Dict[str, int] = field(init=False)

    def __post_init__(self):
        self.name = self.read.query_name
        self.sequence = self.read.query_sequence
        self.length = self.read.query_length
        self.aligned_length = self.read.reference_length
        self.contig = self.read.reference_name
        self.start = self.read.reference_start
        self.end = self.read.reference_end
        self.strand = '-' if self.read.is_reverse else '+'
        self.mapq = self.read.mapping_quality
        self.cigartuples = self.read.cigartuples
        self.cigar = self.read.cigarstring
        self.edit_distance = self.read.get_tag('NM') if self.read.has_tag('NM') else None
        self.is_supplementary = self.read.is_supplementary
        self.is_secondary = self.read.is_secondary
        self.methylation = self.read.get_tag('MM') if self.read.has_tag('MM') else None

        left_clip = right_clip = 0 # init clip values
        if self.cigartuples:  # code 4 = soft clipped
            if self.cigartuples[0][0] == 4:
                left_clip = self.cigartuples[0][1]
            if self.cigartuples[-1][0] == 4:
                right_clip = self.cigartuples[-1][1]
        self.soft_clip_left = left_clip
        self.soft_clip_right = right_clip

        # derived attributes
        self.cigar_stats = self._calculate_cigar_stats(self.cigartuples)

    @classmethod
    def from_pysam_read(cls, read: pysam.AlignedSegment) -> 'AlignedRead':
        return cls(read=read)

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
    reads: List[AlignedRead]
    
    def __len__(self) -> int:
        return len(self.reads)
    
    @property
    def mean_mapq(self) -> float:
        return sum(r.mapq for r in self.reads) / len(self)
    
    @property
    def soft_clip_stats(self) -> Dict[str, float]:
        total_len = sum(r.length for r in self.reads)
        total_soft = sum(r.cigar_stats.get('S', 0) for r in self.reads)
        return {
            'pct_softclipped': (total_soft / total_len) * 100,
            'mean_left_clip': sum(r.soft_clip_left for r in self.reads) / len(self),
            'mean_right_clip': sum(r.soft_clip_right for r in self.reads) / len(self)
        }
    
    def compare_stats(self, other: 'ReadGroup') -> Dict[str, float]:
        query, subject = self.soft_clip_stats, other.soft_clip_stats
        
        return {
            'mapq_differential': self.mean_mapq - other.mean_mapq,
            'softclip_pct_differential': query['pct_softclipped'] - subject['pct_softclipped'],
            'left_clip_differential': query['mean_left_clip'] - subject['mean_left_clip'], 
            'right_clip_differential': query['mean_right_clip'] - subject['mean_right_clip']
        }
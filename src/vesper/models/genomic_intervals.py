from dataclasses import dataclass
from typing import Dict, Any, Optional

@dataclass
class GenomicInterval:
    """A genomic interval with optional metadata."""
    chrom: str
    start: int
    end: int
    name: Optional[str] = None
    score: Optional[float] = None
    strand: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        d = {
            'chrom': self.chrom,
            'start': self.start,
            'end': self.end
        }
        if self.name is not None:
            d['name'] = self.name
        if self.score is not None:
            d['score'] = self.score
        if self.strand is not None:
            d['strand'] = self.strand
        if self.metadata is not None:
            d['metadata'] = self.metadata
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> 'GenomicInterval':
        """Create instance from dictionary representation."""
        return cls(
            chrom=d['chrom'],
            start=d['start'],
            end=d['end'],
            name=d.get('name'),
            score=d.get('score'),
            strand=d.get('strand'),
            metadata=d.get('metadata')
        ) 
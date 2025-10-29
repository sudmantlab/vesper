from __future__ import annotations

import json
from dataclasses import dataclass

@dataclass
class GenomicInterval:
    """Represents a genomic interval with optional metadata."""
    chrom: str
    start: int
    end: int
    source: str  # Required source identifier for the annotation
    metadata: dict = None

    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}
        if not self.source:
            raise ValueError("Source must be provided for GenomicInterval")

    def __repr__(self) -> str:
        return f"GenomicInterval(source={self.source}, chrom={self.chrom}, start={self.start}, end={self.end}, metadata={self.metadata})"

    @property
    def length(self) -> int:
        return self.end - self.start
    
    def to_dict(self) -> dict:
        """Convert interval to a flattened dictionary."""
        result = {
            'source': self.source,
            'chrom': self.chrom,
            'start': self.start,
            'end': self.end,
            'length': self.length
        }

        if self.metadata:
            result.update(self.metadata)      
        return result
    
    @staticmethod
    def from_json(json_record):
        """Convert json string/dict to GenomicInterval."""
        if isinstance(json_record, dict):
            loaded = json_record
        else: # string
            loaded = json.loads(json_record)
        metadata = {k: v for k, v in loaded.items() if k not in ['source', 'chrom', 'start', 'end', 'length']}
        return GenomicInterval(
            source=loaded['source'],
            chrom=loaded['chrom'],
            start=loaded['start'],
            end=loaded['end'],
            metadata=metadata
        )

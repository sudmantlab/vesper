from dataclasses import dataclass
import json
import logging

logger = logging.getLogger(__name__)

@dataclass
class RepeatMaskerResult:
    """Represents a RepeatMasker result for the given insertion sequence."""
    repeat_name: str
    repeat_class: str
    sw_score: int
    divergence: float
    deletion: float
    insertion: float
    query_start: int # ex. if 1, then the match started at the first base of the insertion sequence
    query_end: int
    query_left: int # remaining bases of the insertion sequence that are not part of the match
    strand: str
    repeat_start: int
    repeat_end: int
    repeat_left: int
    match_length: int
    match_coverage: float # motif coverage over read length: ex. 100% coverage if the whole insertion sequence is matched

    def to_dict(self) -> dict:
        """Convert annotation to a dictionary."""
        annotation = {
            'repeat_name': self.repeat_name,
            'repeat_class': self.repeat_class,
            'sw_score': self.sw_score,
            'divergence': self.divergence,
            'deletion': self.deletion,
            'insertion': self.insertion,
            'query_start': self.query_start,
            'query_end': self.query_end,
            'query_left': self.query_left,
            'strand': self.strand,
            'repeat_start': self.repeat_start,
            'repeat_end': self.repeat_end,
            'repeat_left': self.repeat_left,
            'match_length': self.match_length,
            'match_coverage': self.match_coverage,
        }
        return annotation
    
    @staticmethod
    def from_json(json_record):
        """Convert json string/dict to RepeatMaskerResult."""
        if isinstance(json_record, dict):
            loaded = json_record
        else: # string
            loaded = json.loads(json_record)
        return RepeatMaskerResult(**loaded) 
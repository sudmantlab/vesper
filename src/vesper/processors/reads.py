from dataclasses import dataclass
from typing import Iterator, Tuple, Optional
import pysam
from pathlib import Path
import logging

from ..models.variants import Variant
from ..models.reads import ReadGroup, AlignedRead


class ReadProcessor:
    """Extracts reads from a BAM file according to variant breakpoints."""

    def __init__(self, bam_path: Path):
        self.bam_path = bam_path
        self._bam: Optional[pysam.AlignmentFile] = None
        self.logger = logging.getLogger(__name__)
    
    def __enter__(self):
        """Open BAM file and check index."""
        # Check for index
        index_path = str(self.bam_path) + '.bai'
        if not Path(index_path).exists():
            self.logger.info(f"Index not found for {self.bam_path}, creating index...")
            try:
                pysam.samtools.index(str(self.bam_path))
                self.logger.info("Index created successfully")
            except Exception as e:
                self.logger.error(f"Failed to create BAM index: {e}")
                raise RuntimeError(f"Failed to create BAM index: {e}")
        
        self._bam = pysam.AlignmentFile(str(self.bam_path), 'rb')
        self.logger.debug(f"Opened BAM file: {self.bam_path}")
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._bam:
            self._bam.close()
            self.logger.debug(f"Closed BAM file: {self.bam_path}")
            self._bam = None
    
    def get_read_groups(self, variant: Variant) -> Tuple[ReadGroup, ReadGroup]:
        """Fetch supporting and non-supporting reads for the variant.
        
        Args:
            Variant
            
        Returns:
            Tuple of (supporting_reads, nonsupporting_reads) where each read group
            contains ALL alignments (primary, secondary, supplementary).
        """
        if not self._bam:
            raise RuntimeError("BAM file not opened. Use with-statement to open file.")
            
        support_reads = []
        nonsupport_reads = []
        
        for read in self._fetch_reads(variant):
            if read.name in variant.rnames:
                support_reads.append(read)
            else:
                nonsupport_reads.append(read)
                
        return ReadGroup(support_reads), ReadGroup(nonsupport_reads)
    
    def _fetch_reads(self, variant: Variant, window_size: int = 1000) -> Iterator[AlignedRead]:
        """Fetch all reads from region around variant (including non-primary alignments).
        
        Args:
            Variant
            window_size: Size of window around variant to fetch reads from
            
        Returns:
            Iterator of one or more AlignedReads.
        """
        if not self._bam:
            raise RuntimeError("BAM file not opened. Use with-statement to open file.")
            
        for read in self._bam.fetch(
            variant.contig,
            max(0, variant.position - window_size),
            variant.end + window_size
        ):
            yield AlignedRead.from_pysam_read(read) 
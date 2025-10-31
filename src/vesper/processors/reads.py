from __future__ import annotations

from typing import Tuple, Optional, Dict, List, TYPE_CHECKING

import pysam
from pathlib import Path
import logging
import os
import threading

from ..models.reads import ReadGroup, AlignedRead

if TYPE_CHECKING:
    from ..models.variants import Variant


class ReadProcessor:
    """Handles BAM access and produces read groups for variants."""

    def __init__(self, bam_path: Path, window_size: int = 1000):
        self.bam_path = bam_path
        self.window_size = window_size
        self.logger = logging.getLogger(__name__)
        self._bam: Optional[pysam.AlignmentFile] = None
        self._fetch_lock = threading.RLock()

    def __enter__(self):
        index_path = str(self.bam_path) + ".bai"
        if not Path(index_path).exists():
            self.logger.info(f"Index not found for {self.bam_path}, creating index...")
            try:
                pysam.samtools.index(str(self.bam_path))
                self.logger.info("Index created successfully")
            except Exception as exc:  # pragma: no cover
                self.logger.error(f"Failed to create BAM index: {exc}")
                raise RuntimeError(f"Failed to create BAM index: {exc}") from exc

        threads = max(1, os.cpu_count() // 2)
        self._bam = pysam.AlignmentFile(str(self.bam_path), "rb", threads=threads)
        self.logger.debug(f"Opened BAM file: {self.bam_path}")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._bam:
            self._bam.close()
            self.logger.debug(f"Closed BAM file: {self.bam_path}")
            self._bam = None

    def _fetch_reads(self, variant: Variant) -> Tuple[List[pysam.AlignedSegment], List[str]]:
        if not self._bam:
            raise RuntimeError("BAM file not opened")

        with self._fetch_lock:
            start = max(0, variant.position - self.window_size)
            end = min(self._bam.get_reference_length(variant.chrom), variant.position + self.window_size)
            all_reads = list(self._bam.fetch(variant.chrom, start, end))
        supporting_names = set(variant.rnames)
        return all_reads, supporting_names

    @staticmethod
    def _should_replace(current: AlignedRead, candidate: AlignedRead) -> bool:
        current_primary = not current.is_secondary and not current.is_supplementary
        candidate_primary = not candidate.is_secondary and not candidate.is_supplementary
        if candidate_primary and not current_primary:
            return True
        return False

    def get_read_groups(self, variant: Variant) -> Tuple[ReadGroup, ReadGroup]:
        all_reads, supporting_names = self._fetch_reads(variant)

        support_reads: Dict[str, AlignedRead] = {}
        nonsupport_reads: Dict[str, AlignedRead] = {}

        for read in all_reads:
            aligned_read = AlignedRead.from_pysam_read(read)
            target = support_reads if read.query_name in supporting_names else nonsupport_reads
            existing = target.get(aligned_read.name)
            if existing is None or self._should_replace(existing, aligned_read):
                target[aligned_read.name] = aligned_read

        return ReadGroup(list(support_reads.values())), ReadGroup(list(nonsupport_reads.values()))

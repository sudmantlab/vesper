"""Thread-based optimization for annotation processing.

This module provides a simplified, thread-safe annotation processor that uses
threading instead of multiprocessing to avoid the complexity of process spawning.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union, Dict, Tuple, Iterator
import logging
import sqlite3
import time
import threading
from concurrent.futures import ThreadPoolExecutor
from queue import Queue
import os

# Import after ensuring module is importable
try:
    from ..processors.annotations import AnnotationProcessor, BEDProcessor, GenomicInterval
except ImportError:
    # If running standalone, use relative import
    from .annotations import AnnotationProcessor, BEDProcessor, GenomicInterval


class ThreadSafeBEDProcessor(BEDProcessor):
    """Thread-safe implementation of BEDProcessor that uses thread-local sqlite connections."""
    
    def __init__(self, filepath: Union[str, Path]):
        super().__init__(filepath)
        self._local = threading.local()
        self._local.conn = None
        
    def _get_connection(self):
        """Get a thread-local connection to the database."""
        if not hasattr(self._local, 'conn') or self._local.conn is None:
            self._local.conn = sqlite3.connect(str(self.db_path))
            self._local.conn.row_factory = sqlite3.Row
            # Enable optimizations
            self._local.conn.execute("PRAGMA cache_size = -10000")  # 10MB cache
            self._local.conn.execute("PRAGMA temp_store = MEMORY")
        return self._local.conn
    
    def find_overlaps(self, chrom: str, pos: int) -> List[GenomicInterval]:
        """Thread-safe implementation of find_overlaps using thread-local connections."""
        conn = self._get_connection()
        
        query = '''
            SELECT * FROM intervals
            WHERE chrom = ?
            AND start <= ? AND end >= ?
        '''
        
        cursor = conn.execute(query, (chrom, pos, pos))
        
        return [GenomicInterval(
            chrom=row['chrom'],
            start=row['start'],
            end=row['end'],
            metadata={k: row[k] for k in row.keys() if k not in ('chrom', 'start', 'end')}
        ) for row in cursor.fetchall()]

    def find_proximal(self, chrom: str, pos: int, span: int = 100) -> List[GenomicInterval]:
        """Thread-safe implementation of find_proximal using thread-local connections."""
        conn = self._get_connection()
        
        start = max(1, pos - span)
        end = pos + span
        
        query = '''
            SELECT *, 
                   MIN(ABS(start - ?), ABS(end - ?)) as distance
            FROM intervals
            WHERE chrom = ?
            AND ((start BETWEEN ? AND ?) OR (end BETWEEN ? AND ?))
            ORDER BY distance
        '''
        
        cursor = conn.execute(query, (pos, pos, chrom, start, end, start, end))
        
        return [GenomicInterval(
            chrom=row['chrom'],
            start=row['start'],
            end=row['end'],
            metadata={k: row[k] for k in row.keys() if k not in ('chrom', 'start', 'end', 'distance')}
        ) for row in cursor.fetchall()]
    
    def process_single_variant(self, variant, proximal_span=500):
        """Process a single variant, adding annotations directly."""
        if variant.variant.sv_type.name != 'INS':
            return
            
        # Find directly overlapping features at the breakpoint
        overlaps = self.find_overlaps(
            variant.variant.contig,
            variant.variant.position
        )
        variant.overlapping_features.extend(overlaps)
        
        # Find proximal features at the breakpoint
        proximal = self.find_proximal(
            variant.variant.contig,
            variant.variant.position,
            proximal_span
        )
        
        # Add proximal features, excluding any that overlap
        seen = set()
        overlapping_keys = {(f.chrom, f.start, f.end) for f in variant.overlapping_features}
        
        for feature in proximal:
            feature_key = (feature.chrom, feature.start, feature.end)
            if feature_key not in seen and feature_key not in overlapping_keys:
                variant.proximal_features.append(feature)
                seen.add(feature_key)
    
    def process_variants_parallel(self, variants, proximal_span=500, n_workers=None):
        """Process variants in parallel using thread pool.
        
        This avoids multiprocessing complexity by using threads with thread-local
        connections, which is still effective for I/O bound operations.
        
        Args:
            variants: List of VariantAnalysis objects
            proximal_span: Distance for proximal feature search
            n_workers: Number of worker threads (defaults to 2x CPU count for I/O bound)
        """
        if n_workers is None:
            # For I/O bound, can use more threads than cores
            n_workers = min(32, max(4, os.cpu_count() * 2))
            
        # Filter to insertion variants first
        ins_variants = [v for v in variants if v.variant.sv_type.name == 'INS']
        
        self.logger.info(f"Processing {len(ins_variants)} variants with {n_workers} threads")
        start_time = time.time()
        
        # Process in batches using thread pool
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            # Submit all variants to the pool
            futures = [executor.submit(self.process_single_variant, variant, proximal_span) 
                      for variant in ins_variants]
            
            # Wait for completion and count processed
            completed = 0
            total = len(futures)
            for future in futures:
                # This will raise any exceptions from the thread
                future.result()
                completed += 1
                if completed % 100 == 0 or completed == total:
                    self.logger.info(f"Progress: {completed}/{total} variants processed")
        
        elapsed = time.time() - start_time
        self.logger.info(f"Completed processing {len(ins_variants)} variants in {elapsed:.2f} seconds")

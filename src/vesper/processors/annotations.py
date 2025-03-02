from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union
import logging
import sqlite3
import threading
import os
import time
from concurrent.futures import ThreadPoolExecutor


@dataclass
class GenomicInterval:
    """Represents a genomic interval with optional metadata."""
    chrom: str
    start: int
    end: int
    metadata: dict = None

    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}

    @property
    def length(self) -> int:
        return self.end - self.start

class AnnotationProcessor(ABC):
    """Base class for processing genomic annotation files."""
    
    def __init__(self, filepath: Union[str, Path]):
        self.logger = logging.getLogger(__name__)
        self.filepath = Path(filepath)
        self.db_path: Optional[Path] = None
        self.conn = None

    def __enter__(self):
        """Initialize database connection and build if needed."""
        if self.db_path is None:
            self.db_path = self.filepath.with_suffix(self.filepath.suffix + '.sqlite')
        
        if not self.db_path.exists():
            self.logger.info(f"Database not found, building new database for {self.filepath}")
            self._build_sqlite_db()
        
        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.row_factory = sqlite3.Row
        
        # check if database is properly initialized
        try:
            cursor = self.conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='intervals'")
            if cursor.fetchone() is None:
                self.logger.warning(f"Database exists but appears invalid, rebuilding: {self.db_path}")
                self.conn.close()
                self.db_path.unlink()
                self._build_sqlite_db()
                self.conn = sqlite3.connect(str(self.db_path))
                self.conn.row_factory = sqlite3.Row
        except sqlite3.Error as e:
            self.logger.error(f"Error verifying database: {e}")
            if self.conn:
                self.conn.close()
            self.db_path.unlink(missing_ok=True)
            raise RuntimeError(f"Failed to initialize database: {e}")
            
        return self
    
    def close(self):
        """Close the database connection."""
        if self.conn:
            self.conn.close()
            self.conn = None

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @abstractmethod
    def _build_sqlite_db(self) -> 'AnnotationProcessor':
        """Placeholder for format-specific DB building logic."""
        pass

    def find_overlaps(self, chrom: str, pos: int) -> List[GenomicInterval]:
        """Find all features overlapping the given position.
        
        Args:
            chrom: Chromosome/contig name
            pos: Position to query (1-based)
            
        Returns:
            List of GenomicInterval objects
        """
        if not self.conn:
            raise RuntimeError("Database connection not initialized. Use with-statement to open processor.")
            
        query = '''
            SELECT * FROM intervals
            WHERE chrom = ?
            AND start <= ? AND end >= ?
        '''
        
        cursor = self.conn.execute(query, (chrom, pos, pos))
        
        return [GenomicInterval(
            chrom=row['chrom'],
            start=row['start'],
            end=row['end'],
            metadata={k: row[k] for k in row.keys() if k not in ('chrom', 'start', 'end')}
        ) for row in cursor.fetchall()]

    def find_proximal(self, chrom: str, pos: int, span: int = 100) -> List[GenomicInterval]:
        """Find features within span distance of the position.
        
        Args:
            chrom: Chromosome name
            pos: Position (1-based)
            span: +/- interval to search
            
        Returns:
            List of GenomicInterval objects
        """
        if not self.conn:
            raise RuntimeError("Database connection not initialized. Use with-statement to open processor.")
            
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
        
        cursor = self.conn.execute(query, (pos, pos, chrom, start, end, start, end))
        
        return [GenomicInterval(
            chrom=row['chrom'],
            start=row['start'],
            end=row['end'],
            metadata={k: row[k] for k in row.keys() if k not in ('chrom', 'start', 'end', 'distance')}
        ) for row in cursor.fetchall()]

class BEDProcessor(AnnotationProcessor):
    """Thread-safe implementation of BED file processor that uses thread-local sqlite connections."""
    
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

    def _build_sqlite_db(self, rebuild: bool = False, batch_size: int = 100000) -> 'BEDProcessor':
        """Create a new SQLite database from a BED file if it doesn't already exist.
        
        Args:
            rebuild: If True, rebuild the database even if it exists
            batch_size: Number of records to insert in each batch
            
        Returns:
            Self for method chaining
        """
        if self.db_path is None:
            self.db_path = self.filepath.with_suffix(self.filepath.suffix + '.sqlite')
        
        # Check if database already exists and is valid
        if not rebuild and self.db_path.exists():
            try:
                # Test if we can connect and query
                conn = sqlite3.connect(str(self.db_path))
                cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='intervals'")
                if cursor.fetchone() is not None:
                    self.logger.info(f"Using existing database: {self.db_path}")
                    conn.close()
                    return self
                conn.close()
            except sqlite3.Error:
                self.logger.warning(f"Existing database {self.db_path} appears corrupt, rebuilding...")
                self.db_path.unlink(missing_ok=True)
        
        self.logger.info(f"Building new database from {self.filepath}")
        conn = sqlite3.connect(str(self.db_path))
        
        # Drop existing table if rebuilding
        conn.execute('DROP TABLE IF EXISTS intervals')
        
        conn.execute('''
            CREATE TABLE intervals (
                chrom TEXT,
                start INTEGER,
                end INTEGER,
                name TEXT,
                score REAL,
                strand TEXT
            )
        ''')
        conn.execute('CREATE INDEX idx_chrom_start ON intervals (chrom, start)')
        conn.execute('CREATE INDEX idx_chrom_end ON intervals (chrom, end)')

        self.logger.info(f"Indexing {self.filepath}...")
        batch = []

        with open(self.filepath, 'r') as f:
            for i, line in enumerate(f):
                if line.startswith('#') or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                    
                # Convert BED's 0-based to 1-based for consistency with VCF
                chrom = fields[0].strip(' ')
                start = int(fields[1]) + 1  # Convert to 1-based
                end = int(fields[2])
                name = fields[3] if len(fields) > 3 else None
                score = float(fields[4]) if len(fields) > 4 and fields[4] != '.' else None
                strand = fields[5] if len(fields) > 5 else None
                
                batch.append((chrom, start, end, name, score, strand))
                
                if len(batch) >= batch_size:
                    conn.executemany(
                        'INSERT INTO intervals VALUES (?, ?, ?, ?, ?, ?)',
                        batch
                    )
                    conn.commit()
                    batch = []
                    self.logger.info(f"Processed {i+1:,} lines...")
            
            # Insert remaining records
            if batch:
                conn.executemany(
                    'INSERT INTO intervals VALUES (?, ?, ?, ?, ?, ?)',
                    batch
                )
                conn.commit()
        
        self.logger.info(f"Completed indexing {self.filepath} -> {self.db_path}")
        conn.close()
        return self
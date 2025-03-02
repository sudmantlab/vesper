
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union, Dict, Tuple, Iterator
import logging
import sqlite3
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
import time
import os

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


# Worker function that must be at module level for multiprocessing
def _process_variant_batch(db_path, variants_data, proximal_span):
    """Process a batch of variants in a separate process.
    
    Args:
        db_path: Path to the SQLite database
        variants_data: List of (idx, contig, position, sv_type) tuples
        proximal_span: Distance for proximal feature search
        
    Returns:
        List of (idx, overlaps, proximals) tuples
    """
    # We need to set an exception handler for multiprocessing
    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("worker")
    
    try:
        # Connect to database in this process
        conn = sqlite3.connect(str(db_path))
        conn.row_factory = sqlite3.Row
        
        results = []
        
        for idx, contig, position, sv_type in variants_data:
            # Skip non-insertion variants if needed
            if sv_type != 'INS':
                continue
                
            # Find overlaps
            query_overlaps = '''
                SELECT * FROM intervals
                WHERE chrom = ?
                AND start <= ? AND end >= ?
            '''
            
            cursor = conn.execute(query_overlaps, (contig, position, position))
            overlaps = [_row_to_interval(row) for row in cursor.fetchall()]
            
            # Find proximals
            start = max(1, position - proximal_span)
            end = position + proximal_span
            
            query_proximal = '''
                SELECT *, 
                       MIN(ABS(start - ?), ABS(end - ?)) as distance
                FROM intervals
                WHERE chrom = ?
                AND ((start BETWEEN ? AND ?) OR (end BETWEEN ? AND ?))
                ORDER BY distance
            '''
            
            cursor = conn.execute(query_proximal, (position, position, contig, start, end, start, end))
            proximals = [_row_to_interval(row) for row in cursor.fetchall()]
            
            results.append((idx, overlaps, proximals))
        
        conn.close()
        return results
    except Exception as e:
        logger.error(f"Worker error: {str(e)}")
        raise


def _row_to_interval(row):
    """Convert a SQLite row to a GenomicInterval."""
    return GenomicInterval(
        chrom=row['chrom'],
        start=row['start'],
        end=row['end'],
        metadata={k: row[k] for k in row.keys() if k not in ('chrom', 'start', 'end', 'distance')}
    )


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
        
        # Check if database is properly initialized
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
    
    def process_variants_multiprocess(self, variants: List['VariantAnalysis'], 
                                    proximal_span: int = 500, 
                                    batch_size: int = 100,
                                    n_workers: int = None) -> None:
        """Process multiple variants using multiprocessing.
        
        Args:
            variants: List of VariantAnalysis objects to annotate
            proximal_span: Distance for proximal feature search
            batch_size: Number of variants to process per batch
            n_workers: Number of worker processes (defaults to CPU count)
        """
        if not self.db_path or not self.db_path.exists():
            raise RuntimeError("Database not initialized")
            
        if n_workers is None:
            n_workers = max(1, multiprocessing.cpu_count() - 1)
        
        # Must close this connection before forking processes
        self.close()
        
        # Create data for processing (only the necessary fields to minimize serialization)
        variant_data = []
        for idx, variant in enumerate(variants):
            if hasattr(variant, 'variant') and hasattr(variant.variant, 'sv_type'):
                variant_data.append((
                    idx,
                    variant.variant.contig,
                    variant.variant.position,
                    variant.variant.sv_type.name
                ))
        
        # Create batches
        batches = [variant_data[i:i+batch_size] for i in range(0, len(variant_data), batch_size)]
        
        # Process batches in parallel
        start_time = time.time()
        self.logger.info(f"Starting multiprocessing with {n_workers} workers for {len(variant_data)} variants")
        
        # Need to use the freeze_support and __name__ == '__main__' pattern for multiprocessing
        # The actual process pool is created in the main script
        
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = []
            for batch in batches:
                futures.append(
                    executor.submit(_process_variant_batch, self.db_path, batch, proximal_span)
                )
            
            # Get results as they complete
            completed = 0
            total_variants = 0
            
            for future in futures:
                try:
                    batch_results = future.result()
                    total_variants += len(batch_results)
                    
                    # Apply results to original variants
                    for idx, overlaps, proximals in batch_results:
                        variant = variants[idx]
                        
                        # Add overlapping features
                        variant.overlapping_features.extend(overlaps)
                        
                        # Add proximal features, excluding any that overlap
                        seen = set()
                        overlapping_keys = {(f.chrom, f.start, f.end) for f in variant.overlapping_features}
                        
                        for feature in proximals:
                            feature_key = (feature.chrom, feature.start, feature.end)
                            if feature_key not in seen and feature_key not in overlapping_keys:
                                variant.proximal_features.append(feature)
                                seen.add(feature_key)
                    
                    completed += 1
                    if completed % max(1, len(futures)//10) == 0:
                        elapsed = time.time() - start_time
                        self.logger.info(f"Processed {completed}/{len(futures)} batches, {total_variants} variants in {elapsed:.2f}s")
                
                except Exception as e:
                    self.logger.error(f"Error processing batch: {e}")
                    raise
        
        self.logger.info(f"Multiprocessing completed in {time.time() - start_time:.2f}s")
        
        # Reopen connection for this instance
        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.row_factory = sqlite3.Row


class MultiprocessBEDProcessor(AnnotationProcessor):
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
        conn.execute('CREATE INDEX idx_chrom_pos ON intervals (chrom, start, end)')

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


class GFF3Processor(AnnotationProcessor):
    def _build_sqlite_db(self) -> 'AnnotationProcessor':
        """Placeholder for format-specific DB building logic."""
        pass

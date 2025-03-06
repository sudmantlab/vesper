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
    
    def to_dict(self) -> dict:
        """Convert interval to a flattened dictionary."""
        result = {
            'chrom': self.chrom,
            'start': self.start,
            'end': self.end,
            'length': self.length
        }

        if self.metadata:
            result.update(self.metadata)      
        return result

class AnnotationProcessor(ABC):
    """Base class for processing genomic annotation files."""
    
    def __init__(self, filepath: Union[str, Path]):
        self.logger = logging.getLogger(__name__)
        self.filepath = Path(filepath)
        self.db_path: Optional[Path] = None
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

    def __enter__(self):
        """Initialize database connection and build if needed."""
        if self.db_path is None:
            self.db_path = self.filepath.with_suffix(self.filepath.suffix + '.sqlite')
        
        if not self.db_path.exists():
            self.logger.info(f"Database not found, building new database for {self.filepath}")
            self._build_sqlite_db()
        
        self.conn = self._get_connection()
        
        # check if database is properly initialized
        try:
            cursor = self.conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='intervals'")
            if cursor.fetchone() is None:
                self.logger.warning(f"Database exists but appears invalid, rebuilding: {self.db_path}")
                self.conn.close()
                self.db_path.unlink()
                self._build_sqlite_db(rebuild=True)
                self.conn = self._get_connection()
        except sqlite3.Error as e:
            self.logger.error(f"Error verifying database: {e}")
            if self.conn:
                self.conn.close()
            self.db_path.unlink(missing_ok=True)
            raise RuntimeError(f"Failed to initialize database: {e}")
            
        return self
    
    def close(self):
        """Close the database connection."""
        if hasattr(self._local, 'conn') and self._local.conn:
            self._local.conn.close()
            self._local.conn = None

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
        """Find features within span distance of the position.
        
        Args:
            chrom: Chromosome name
            pos: Position (1-based)
            span: +/- interval to search
            
        Returns:
            List of GenomicInterval objects
        """
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
    
    def _annotate_variant(self, variant, proximal_span) -> None:
        """Annotate a single variant."""
        overlaps = self.find_overlaps(
            variant.variant.chrom,
            variant.variant.position
        )
        variant.overlapping_features.extend(overlaps)
        
        proximal = self.find_proximal(
            variant.variant.chrom,
            variant.variant.position,
            proximal_span
        )
        
        seen = set() # avoid duplicates
        overlapping_keys = {(f.chrom, f.start, f.end) for f in variant.overlapping_features}
        
        for feature in proximal:
            feature_key = (feature.chrom, feature.start, feature.end)
            if feature_key not in seen and feature_key not in overlapping_keys:
                variant.proximal_features.append(feature)
                seen.add(feature_key)


class BEDProcessor(AnnotationProcessor):
    """Implementation of BED file processor."""
    
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

class GFFProcessor(AnnotationProcessor):
    """Implementation of GFF file processor."""
    def _build_sqlite_db(self, rebuild: bool = False, batch_size: int = 100000) -> 'GFFProcessor':
        """Create a new SQLite database from a GFF file if it doesn't already exist.
        
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
                if len(fields) < 9:  # GFF requires 9 fields
                    continue
             
                # GFFs are already 1-based
                chrom = fields[0].strip()
                source = fields[1].strip() 
                feature = fields[2].strip()
                start = int(fields[3])
                end = int(fields[4])
                score = float(fields[5]) if fields[5] != '.' else None
                strand = fields[6]
                frame = fields[7]

                # Parse attributes into key-value pairs
                metadata = {}
                if fields[8]:
                    attr_pairs = fields[8].strip().split(';')
                    for pair in attr_pairs:
                        if '=' in pair:
                            key, value = pair.split('=', 1)
                            metadata[key.strip()] = value.strip()
                
                batch.append((chrom, start, end, name, score, strand, frame, metadata))
                
                if len(batch) >= batch_size:
                    conn.executemany(
                        'INSERT INTO intervals VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
                        batch
                    )
                    conn.commit()
                    batch = []
                    self.logger.info(f"Processed {i+1:,} lines...")
            
            # Insert remaining records
            if batch:
                conn.executemany(
                    'INSERT INTO intervals VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
                    batch
                )
                conn.commit()
        
        self.logger.info(f"Completed indexing {self.filepath} -> {self.db_path}")
        conn.close()
        return self
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union, Dict, Any, Iterator
import logging
import sqlite3
import threading
import os
import time
import gzip
from concurrent.futures import ThreadPoolExecutor
import json

from ..models.intervals import GenomicInterval


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

class AnnotationProcessor(ABC):
    """Base class for processing genomic annotation files."""
    
    def __init__(self, filepath: Union[str, Path], source_name: str):
        self.logger = logging.getLogger(__name__)
        self.filepath = Path(filepath)
        self.source_name = source_name
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
        
        results = []
        for row in cursor.fetchall():
            metadata = json.loads(row['metadata']) if 'metadata' in row.keys() else {}
            results.append(GenomicInterval(
                source=self.source_name,
                chrom=row['chrom'],
                start=row['start'],
                end=row['end'],
                metadata=metadata
            ))
        return results

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
        
        results = []
        for row in cursor.fetchall():
            metadata = json.loads(row['metadata']) if 'metadata' in row.keys() else {}
            metadata['distance'] = row['distance']  # Add distance to metadata
            results.append(GenomicInterval(
                source=self.source_name,
                chrom=row['chrom'],
                start=row['start'],
                end=row['end'],
                metadata=metadata
            ))
        return results
    
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

    def iter_intervals(self) -> Iterator[GenomicInterval]:
        """Iterate over genomic intervals in the annotation file."""
        raise NotImplementedError

    def get_overlapping(self, chrom: str, start: int, end: int) -> List[GenomicInterval]:
        """Get intervals overlapping the given region."""
        return [
            interval for interval in self.iter_intervals()
            if (interval.chrom == chrom and
                interval.start <= end and
                interval.end >= start)
        ]

class BEDProcessor(AnnotationProcessor):
    """Implementation of BED file processor."""
    
    def _build_sqlite_db(self, rebuild: bool = False, batch_size: int = 100000) -> 'BEDProcessor':
        """Create a new SQLite database from a BED file if it doesn't already exist.
        
        Args:
            rebuild: If True, rebuild the database even if it exists
            batch_size: Number of records to insert in each batch
            
        Returns:
            Self for method chaining
            
        Raises:
            ValueError: If a BED field is empty when higher-numbered fields are present
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
        
        # Create table with only mandatory fields plus metadata JSON
        conn.execute('''
            CREATE TABLE intervals (
                chrom TEXT,
                start INTEGER,
                end INTEGER,
                metadata TEXT  -- JSON field for all optional columns
            )
        ''')
        conn.execute('CREATE INDEX idx_chrom_start ON intervals (chrom, start)')
        conn.execute('CREATE INDEX idx_chrom_end ON intervals (chrom, end)')

        self.logger.info(f"Indexing {self.filepath}...")
        batch = []

        if str(self.filepath).endswith(('.gz', '.gzip')):
            file_opener = gzip.open
            mode = 'rt'
        else:
            file_opener = open
            mode = 'r'

        with file_opener(self.filepath, mode) as f:
            for i, line in enumerate(f, start=1):  # 1-based line counting for error messages
                if line.startswith('#') or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 3: # BED requires 3 fields
                    raise ValueError(f"Line {i}: BED format requires at least 3 fields, found {len(fields)}")
                    
                for idx, field in enumerate(fields[:3]):
                    if not field.strip():
                        raise ValueError(f"Line {i}: Mandatory BED field {idx+1} is empty")
                    
                # Convert BED's 0-based to 1-based for consistency with VCF
                chrom = fields[0].strip(' ')
                start = int(fields[1]) + 1
                end = int(fields[2])
                
                optional = {3: ('name', str), 
                            4: ('score', int), 
                            5: ('strand', str), 
                            6: ('thickStart', int), 
                            7: ('thickEnd', int), 
                            8: ('itemRgb', str), 
                            9: ('blockCount', int), 
                            10: ('blockSizes', str), 
                            11: ('blockStarts', str)}
                metadata = {}
      
                for idx, value in enumerate(fields[3:], start=3):
                    if value == '': 
                        value = "." # placeholder/empty sub
                    elif idx <= 12:
                        metadata[optional[idx][0]] = optional[idx][1](value)
                    metadata[f'field{idx}'] = str(value) # cast generic name/value pairs as strings
                    
                metadata_json = json.dumps(metadata) if metadata else '{}'
                
                batch.append((chrom, start, end, metadata_json))
                
                if len(batch) >= batch_size:
                    conn.executemany(
                        'INSERT INTO intervals VALUES (?, ?, ?, ?)',
                        batch
                    )
                    conn.commit()
                    batch = []
                    self.logger.info(f"Processed {i:,} lines...")
            
            # Insert remaining records
            if batch:
                conn.executemany(
                    'INSERT INTO intervals VALUES (?, ?, ?, ?)',
                    batch
                )
                conn.commit()
        
        self.logger.info(f"Completed indexing {self.filepath} -> {self.db_path}")
        conn.close()
        return self

    def iter_intervals(self) -> Iterator[GenomicInterval]:
        """Iterate over intervals in BED format."""
        open_fn = gzip.open if str(self.filepath).endswith('.gz') else open
        with open_fn(self.filepath, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                    
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                
                interval = GenomicInterval(
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=fields[3] if len(fields) > 3 else None,
                    score=float(fields[4]) if len(fields) > 4 else None,
                    strand=fields[5] if len(fields) > 5 else None
                )
                
                yield interval

class GFFProcessor(AnnotationProcessor):
    """Implementation of GFF file processor."""
    def _build_sqlite_db(self, rebuild: bool = False, batch_size: int = 100000) -> 'GFFProcessor':
        """Create a new SQLite database from a GFF file if it doesn't already exist.
        
        Args:
            rebuild: If True, rebuild the database even if it exists
            batch_size: Number of records to insert in each batch
            
        Returns:
            Self for method chaining
            
        Raises:
            ValueError: If a GFF field is empty or malformed
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
                metadata TEXT  -- JSON field for all optional columns
            )
        ''')
        conn.execute('CREATE INDEX idx_chrom_start ON intervals (chrom, start)')
        conn.execute('CREATE INDEX idx_chrom_end ON intervals (chrom, end)')

        self.logger.info(f"Indexing {self.filepath}...")
        batch = []
        
        if str(self.filepath).endswith(('.gz', '.gzip')):
            file_opener = gzip.open
            mode = 'rt'
        else:
            file_opener = open
            mode = 'r'

        with file_opener(self.filepath, mode) as f:
            for i, line in enumerate(f, start=1):  # 1-based line counting for error messages
                if line.startswith('#') or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue
                    
                chrom = fields[0]
                start = int(fields[3]) - 1  # Convert to 0-based
                end = int(fields[4])
                strand = fields[6]
                
                # Parse attributes
                attrs = {}
                attr_str = fields[8]
                for pair in attr_str.strip(';').split(';'):
                    if not pair:
                        continue
                    try:
                        key, value = pair.strip().split('=', 1)
                    except ValueError:
                        try:
                            key, value = pair.strip().split(' ', 1)
                            value = value.strip('"')
                        except ValueError:
                            continue
                    attrs[key.strip()] = value.strip()
                
                interval = GenomicInterval(
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=attrs.get('ID') or attrs.get('gene_id'),
                    score=float(fields[5]) if fields[5] != '.' else None,
                    strand=strand if strand != '.' else None,
                    metadata={
                        'source': fields[1],
                        'type': fields[2],
                        'phase': fields[7],
                        'attributes': attrs
                    }
                )
                
                batch.append((chrom, start, end, json.dumps(interval.to_dict())))
                
                if len(batch) >= batch_size:
                    conn.executemany(
                        'INSERT INTO intervals VALUES (?, ?, ?, ?)',
                        batch
                    )
                    conn.commit()
                    batch = []
                    self.logger.info(f"Processed {i:,} lines...")
            
            # Insert remaining records
            if batch:
                conn.executemany(
                    'INSERT INTO intervals VALUES (?, ?, ?, ?)',
                    batch
                )
                conn.commit()
        
        self.logger.info(f"Completed indexing {self.filepath} -> {self.db_path}")
        conn.close()
        return self

    def iter_intervals(self) -> Iterator[GenomicInterval]:
        """Iterate over intervals in GFF/GTF format."""
        open_fn = gzip.open if str(self.filepath).endswith('.gz') else open
        with open_fn(self.filepath, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue
                    
                chrom = fields[0]
                start = int(fields[3]) - 1  # Convert to 0-based
                end = int(fields[4])
                strand = fields[6]
                
                # Parse attributes
                attrs = {}
                attr_str = fields[8]
                for pair in attr_str.strip(';').split(';'):
                    if not pair:
                        continue
                    try:
                        key, value = pair.strip().split('=', 1)
                    except ValueError:
                        try:
                            key, value = pair.strip().split(' ', 1)
                            value = value.strip('"')
                        except ValueError:
                            continue
                    attrs[key.strip()] = value.strip()
                
                interval = GenomicInterval(
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=attrs.get('ID') or attrs.get('gene_id'),
                    score=float(fields[5]) if fields[5] != '.' else None,
                    strand=strand if strand != '.' else None,
                    metadata={
                        'source': fields[1],
                        'type': fields[2],
                        'phase': fields[7],
                        'attributes': attrs
                    }
                )
                
                yield interval
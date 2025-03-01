from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union
import logging

import sqlite3


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

    # @abstractmethod
    # def _parse_file(self) -> Iterator[GenomicInterval]:
    #     """Placeholder for format-specific parsing logic."""
    #     pass

    # def load(self) -> None:
    #     """Load and index interval tree for annotations."""
    #     if self._loaded:
    #         return

    #     self.logger.info(f"Loading annotations from {self.filepath}")
    #     for interval in self._parse_file():
    #         if interval.chrom not in self._intervals:
    #             self._intervals[interval.chrom] = intervaltree.IntervalTree()
            
    #         self._intervals[interval.chrom].addi(interval.start, interval.end, interval)

    #     self._loaded = True
    #     self.logger.info(f"Loaded annotations from {len(self._intervals)} chromosomes")

    # def find_overlaps(self, chrom: str, start: int, end: int) -> list[GenomicInterval]:
    #     """Find all annotations overlapping the query interval."""
    #     if not self._loaded:
    #         self.load()

    #     if chrom not in self._intervals:
    #         return []

    #     if start != end:
    #         overlaps = self._intervals[chrom].overlap(start, end)
    #     else: # point query
    #         overlaps = self._intervals[chrom].at(start)

    #     return [interval.data for interval in overlaps]

    # def find_proximal(self, chrom: str, pos: int, span: int = 500) -> list[GenomicInterval]:
    #     """Find annotations within the +/- span of the query position."""
    #     if not self._loaded:
    #         self.load()

    #     if chrom not in self._intervals:
    #         return []

    #     # Find all intervals that overlap the span region around pos
    #     start = max(0, pos - span)
    #     end = pos + span
    #     overlaps = self._intervals[chrom].overlap(start, end)
        
    #     # Sort by distance to query position
    #     intervals = list(overlaps)
    #     intervals.sort(key=lambda x: min(abs(x.begin - pos), abs(x.end - pos)))
        
    #     return [interval.data for interval in intervals]


class BEDProcessor(AnnotationProcessor):
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

#     def _parse_file(self) -> Iterator[GenomicInterval]:
#         bed_columns = [('chrom', str), ('start', int), ('end', int), ('name', str), ('score', float), ('strand', str)]

#         df = pd.read_csv(
#             self.filepath,
#             sep='\t',
#             header=None,
#             comment='#',
#             na_values='.',
#             low_memory=False
#         )

#         if df.shape[1] < 3:
#             raise ValueError(f"BED file {self.filepath} must have at least 3 columns")
        
#         # assign names and dtypes according to existing columns
#         names, dtypes = zip(*bed_columns[:df.shape[1]])
#         df.columns = names
#         for col, dtype in zip(names, dtypes):
#             try:
#                 df[col] = df[col].astype(dtype)
#             except (ValueError, TypeError) as e:
#                 self.logger.error(f"Could not convert BED column '{col}' to type {dtype}: {str(e)}")
#                 raise ValueError(f"Failed to parse BED column '{col}' as {dtype}")
        
#         self.logger.info(f"Processing BED file with {df.shape[1]} columns")

#         for idx, row in df.iterrows():
#             try:
#                 if row.start + 1 == row.end: # avoid null interval while still converting to 1-based
#                     yield GenomicInterval(
#                         chrom=str(row.chrom).strip(' '),
#                         start=int(row.start) + 1,
#                         end=int(row.end) + 1,
#                         metadata={
#                             'name': row['name'] if 'name' in df.columns and pd.notna(row['name']) else None, # avoid using index name
#                             'score': float(row.score) if 'score' in df.columns and pd.notna(row.score) else None,
#                             'strand': str(row.strand) if 'strand' in df.columns and pd.notna(row.strand) else None
#                         }
#                     )
#                 else:
#                     # convert to 1-based for VCF compatibility
#                     yield GenomicInterval(
#                         chrom=str(row.chrom).strip(' '),
#                         start=int(row.start) + 1,
#                         end=int(row.end),
#                         metadata={
#                             'name': row['name'] if 'name' in df.columns and pd.notna(row['name']) else None, # avoid using index name
#                             'score': float(row.score) if 'score' in df.columns and pd.notna(row.score) else None,
#                             'strand': str(row.strand) if 'strand' in df.columns and pd.notna(row.strand) else None
#                         }
#                     )
#             except (ValueError, TypeError) as e:
#                 self.logger.error(f"Malformed data in row {idx + 1}: {str(e)}")
#                 raise ValueError(f"Failed to parse row {idx + 1} of BED file")


class GFF3Processor(AnnotationProcessor):
    def _build_sqlite_db(self) -> 'AnnotationProcessor':
        """Placeholder for format-specific DB building logic."""
        pass
    # def _parse_file(self) -> Iterator[GenomicInterval]:
    #     # create temp db for fast parsing
    #     db = gffutils.create_db(
    #         str(self.filepath),
    #         ':memory:',
    #         merge_strategy='create_unique',
    #         force=True
    #     )

    #     for feature in db.all_features():
    #         yield GenomicInterval(
    #             chrom=feature.seqid,
    #             start=feature.start - 1, # convert to 0-based
    #             end=feature.end,
    #             metadata={
    #                 'type': feature.featuretype,
    #                 'source': feature.source,
    #                 'phase': feature.frame,
    #                 'attributes': dict(feature.attributes)
    #             }
    #         )

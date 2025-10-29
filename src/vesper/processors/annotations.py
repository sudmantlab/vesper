from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Optional, Union, Dict, Tuple
import logging
import sqlite3
import threading
import time
import gzip
import json

from ..models.interval import GenomicInterval

class AnnotationProcessor(ABC):
    """Base class for processing genomic annotation files."""
    
    def __init__(self, filepath: Union[str, Path], source_name: str, rebuild: bool = False):
        self.logger = logging.getLogger(__name__)
        self.filepath = Path(filepath)
        self.source_name = source_name
        self.db_path: Optional[Path] = None
        self._local = threading.local()
        self._local.conn = None
        self.rebuild = rebuild
        
        if rebuild:
            self.logger.info(f"Rebuild flag set for {source_name} ({filepath.name})")

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
        
        if self.rebuild:
            if self.db_path.exists():
                self.logger.info(f"Rebuild flag set, rebuilding database for {self.filepath}")
                self.db_path.unlink()
            self._build_sqlite_db(rebuild=True)
        elif not self.db_path.exists():
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

    def _find_overlap(self, chrom: str, pos: int) -> List[GenomicInterval]:
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

    def _find_proximal(self, chrom: str, pos: int, span: int = 100) -> List[GenomicInterval]:
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
    
    def batch_find_overlaps(self, positions: List[Tuple[str, int]]) -> Dict[Tuple[str, int], List[GenomicInterval]]:
        """Find overlapping features for multiple positions in a single query."""
        if not positions:
            return {}
            
        conn = self._get_connection()
        
        # Build UNION ALL query for all positions
        union_parts = []
        params = []
        for chrom, pos in positions:
            union_parts.append("SELECT ?, ?, chrom, start, end, metadata FROM intervals WHERE chrom = ? AND start <= ? AND end >= ?")
            # params [0-1]: query variant at (chrom, pos) against union of dbs
            # params [2-4]: formulate overlap query on union of dbs where chrom = chrom, start <= pos, end >= pos
            params.extend([chrom, pos, chrom, pos, pos])
        
        query = " UNION ALL ".join(union_parts)
        cursor = conn.execute(query, params)
        
        results = {}
        for row in cursor.fetchall():
            query_chrom, query_pos = row[0], row[1]
            pos_key = (query_chrom, query_pos)
            
            if pos_key not in results:
                results[pos_key] = []
                
            metadata = json.loads(row[5]) if row[5] else {}

            # append interval: denote actual feature start/end points
            results[pos_key].append(GenomicInterval(
                source=self.source_name,
                chrom=row[2],
                start=row[3],
                end=row[4],
                metadata=metadata
            ))
        
        return results

    def batch_find_proximal(self, positions: List[Tuple[str, int]], span: int) -> Dict[Tuple[str, int], List[GenomicInterval]]:
        """Find proximal features for multiple positions in a single query."""
        if not positions:
            return {}
            
        conn = self._get_connection()
        
        # Build UNION ALL query for all positions
        union_parts = []
        params = []
        for chrom, pos in positions:
            proximal_start = max(1, pos - span)
            proximal_end = pos + span
            
            union_parts.append(
                "SELECT ?, ?, chrom, start, end, metadata, MIN(ABS(start - ?), ABS(end - ?)) as distance "
                "FROM intervals WHERE chrom = ? AND ((start BETWEEN ? AND ?) OR (end BETWEEN ? AND ?))"
            )
            params.extend([chrom, pos, pos, pos, chrom, proximal_start, proximal_end, proximal_start, proximal_end])
            # params [0-1]: getquery variant at (chrom, pos) against union of dbs
            # params [2-8]: formulate proximal query on union of dbs 
                # [2-3] (pos, pos): MIN(ABS(start - pos), ABS(end - pos)) sets proximity span
                # [4] (chrom): sets query chrom on union of dbs
                # [5-6]: (proximal_start, proximal_end): intersect between feature start and proximal span
                # [7-8]: (proximal_start, proximal_end): intersect between feature end and proximal span
        
        query = " UNION ALL ".join(union_parts) + " ORDER BY distance"
        cursor = conn.execute(query, params)
        
        results = {}
        for row in cursor.fetchall():
            query_chrom, query_pos = row[0], row[1]
            pos_key = (query_chrom, query_pos)
            
            if pos_key not in results:
                results[pos_key] = []
                
            metadata = json.loads(row[5]) if row[5] else {}
            metadata['distance'] = row[6]

            # append interval: denote actual feature start/end points
            results[pos_key].append(GenomicInterval(
                source=self.source_name,
                chrom=row[2],
                start=row[3],
                end=row[4],
                metadata=metadata
            ))
        
        return results

    def batch_annotate_variants(self, variants: List, proximal_span: int) -> None:
        """Annotate multiple variants in batch."""
        if not variants:
            return
            
        self.logger.info(f"Starting batch annotation for {self.source_name} with {len(variants)} variants")
        start_time = time.time()
        
        positions = [(v.variant.chrom, v.variant.position) for v in variants]
        
        overlap_results = self.batch_find_overlaps(positions)
        overlap_count = sum(len(features) for features in overlap_results.values())
        self.logger.debug(f"Batch overlap query returned {overlap_count} features")
        
        proximal_results = self.batch_find_proximal(positions, proximal_span)
        proximal_count = sum(len(features) for features in proximal_results.values())
        self.logger.debug(f"Batch proximal query returned {proximal_count} features within {proximal_span}bp")
        
        # Assign results back to variants
        for variant in variants:
            pos_key = (variant.variant.chrom, variant.variant.position)
            
            # Add overlaps
            overlaps = overlap_results.get(pos_key, [])
            variant.overlapping_features.extend(overlaps)
            
            # Add proximals (avoiding duplicates)
            proximals = proximal_results.get(pos_key, [])
            seen = set()
            overlapping_keys = {(f.chrom, f.start, f.end) for f in variant.overlapping_features}
            
            for feature in proximals:
                feature_key = (feature.chrom, feature.start, feature.end)
                if feature_key not in seen and feature_key not in overlapping_keys:
                    variant.proximal_features.append(feature)
                    seen.add(feature_key)
        
        elapsed = time.time() - start_time
        self.logger.info(f"Completed batch annotation for {len(variants)} variants in {elapsed:.3f}s")

    def annotate_variant(self, variant, proximal_span) -> None:
        """Annotate a single variant."""
        overlaps = self._find_overlap(
            variant.variant.chrom,
            variant.variant.position
        )
        variant.overlapping_features.extend(overlaps)
        
        proximal = self._find_proximal(
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
                if len(fields) < 9:  # GFF requires 9 fields
                    raise ValueError(f"Line {i}: GFF format requires 9 fields, found {len(fields)}")
             
                for idx, field in enumerate(fields[:3]):
                    if not field.strip():
                        raise ValueError(f"Line {i}: Mandatory GFF field {idx+1} is empty")
                
                chrom = fields[0].strip(' ')
                
                try:
                    start = int(fields[3])
                    end = int(fields[4])
                except ValueError:
                    raise ValueError(f"Line {i}: Invalid start/end coordinates: {fields[3]}, {fields[4]}")
                
                # write all non-mandatory fields into metadata json
                metadata = {
                    'source': fields[1].strip(),
                    'feature_type': fields[2].strip(),
                    'score': None if fields[5] == '.' else float(fields[5]),
                    'strand': fields[6],
                    'frame': fields[7]
                }

                attributes = fields[8].strip()
                if attributes:
                    parts = attributes.split(';')
                    for part in parts:
                        if not part.strip():
                            continue
                        if '=' in part:
                            key, value = part.split('=', 1)
                            metadata[key.strip()] = value.strip()
                
                for idx, value in enumerate(fields[9:], start=1):
                    if value == '':
                        value = "." # placeholder/empty sub
                    metadata[f'field{idx}'] = str(value) # cast generic name/value pairs as strings
                
                metadata_json = json.dumps(metadata)
                
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

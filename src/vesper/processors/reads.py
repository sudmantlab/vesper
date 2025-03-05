from dataclasses import dataclass, field
from typing import Iterator, Tuple, Optional, Dict, List, Set, Any
import pysam
from pathlib import Path
import logging
import os
import json
import threading
import datetime

from ..models.variants import Variant
from ..models.reads import ReadGroup, AlignedRead
from ..models.registry import ReadMetadata, VariantReadGroups, RegistryMetadata


class ReadProcessor:
    """Handles BAM file access and read group management with registry-based caching.
    
    The ReadProcessor maintains two registries:
    1. Read Registry: Maps read names to AlignedRead objects
    2. Variant Read Groups Registry: Maps variant IDs to lists of supporting/non-supporting read names
    
    These registries enable:
    - Efficient read group lookup for variants
    - Caching of read metadata and computations
    - Persistence of read data between runs
    """

    def __init__(self, bam_path: Path, window_size: int = 1000, registry_dir: Optional[Path] = None,
                 auto_load_registry: bool = True, force_new_registry: bool = False):
        """Initialize ReadProcessor.
        
        Args:
            bam_path: Path to BAM file
            window_size: Size of window around variant position to fetch reads
            registry_dir: Optional directory to save/load registry files
            auto_load_registry: Whether to automatically load existing registry if found
            force_new_registry: Force creation of new registry even if one exists
        """
        self.bam_path = bam_path
        self.window_size = window_size
        self.registry_dir = registry_dir
        self._bam: Optional[pysam.AlignmentFile] = None
        self.logger = logging.getLogger(__name__)
        
        # Registry data structures
        self._read_registry: Dict[str, AlignedRead] = {}
        self._variant_read_groups: Dict[str, Dict[str, List[str]]] = {}
        self._registry_lock = threading.RLock()
        self._registry_modified = False
        
        # Create registry directory if specified
        if registry_dir:
            os.makedirs(registry_dir, exist_ok=True)
            
            # Check for existing registry
            if not force_new_registry and auto_load_registry:
                metadata_path = registry_dir / "registry_metadata.json"
                if metadata_path.exists():
                    self.logger.info("Found existing registry, will load on BAM file open")
                    self._load_existing_registry = True
                else:
                    self._load_existing_registry = False
            else:
                self._load_existing_registry = False

    def __enter__(self):
        """Open BAM file, check/create index, and load registry if configured."""
        index_path = str(self.bam_path) + '.bai'
        if not Path(index_path).exists():
            self.logger.info(f"Index not found for {self.bam_path}, creating index...")
            try:
                pysam.samtools.index(str(self.bam_path))
                self.logger.info("Index created successfully")
            except Exception as e:
                self.logger.error(f"Failed to create BAM index: {e}")
                raise RuntimeError(f"Failed to create BAM index: {e}")
        
        threads = max(1, os.cpu_count() // 2)  # Use 1/2 of available cores for reading BAM
        self._bam = pysam.AlignmentFile(str(self.bam_path), 'rb', threads=threads)
        self.logger.debug(f"Opened BAM file: {self.bam_path}")
        
        # Load registry if configured
        if hasattr(self, '_load_existing_registry') and self._load_existing_registry:
            try:
                self.load_registry(self.registry_dir)
                self.logger.info("Successfully loaded existing registry")
            except Exception as e:
                self.logger.error(f"Failed to load registry: {e}")
                self.logger.info("Continuing with new registry")
                self._read_registry = {}
                self._variant_read_groups = {}
        
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close BAM file and save registry if specified."""
        if self._bam:
            self._bam.close()
            self.logger.debug(f"Closed BAM file: {self.bam_path}")
            self._bam = None
            
        # Save registry on exit if directory specified and modified
        if self.registry_dir and self._registry_modified:
            self._save_complete_registry()
            self._registry_modified = False
    
    def _fetch_reads(self, variant: Variant) -> Tuple[List[pysam.AlignedSegment], List[str]]:
        """Fetch reads overlapping a variant's region.
        
        Args:
            variant: Variant to fetch reads for
            
        Returns:
            Tuple of (all reads in region, names of supporting reads)
        """
        if not self._bam:
            raise RuntimeError("BAM file not opened")
            
        # Determine region to fetch
        start = max(0, variant.position - self.window_size)
        end = min(self._bam.get_reference_length(variant.contig), variant.position + self.window_size)
        
        # Get all reads in the region
        all_reads = list(self._bam.fetch(variant.contig, start, end))
        
        # Set of supporting read names from variant
        supporting_names = set(variant.rnames)
        
        return all_reads, supporting_names
    
    def _get_read(self, read: pysam.AlignedSegment) -> AlignedRead:
        """Get existing AlignedRead from registry or create a new one.
        
        Args:
            read: pysam AlignedSegment to wrap
            
        Returns:
            AlignedRead object (either cached or newly created)
        """
        read_name = read.query_name
        
        with self._registry_lock:
            if read_name in self._read_registry:
                return self._read_registry[read_name]
            else:
                aligned_read = AlignedRead.from_pysam_read(read)
                self._read_registry[read_name] = aligned_read
                return aligned_read
    
    def get_read_groups(self, variant: Variant) -> Tuple[ReadGroup, ReadGroup]:
        """Get supporting and non-supporting read groups for a variant.
        Creates and caches groups if not already in registry.
        
        Args:
            variant: Variant to get read groups for
            
        Returns:
            Tuple of (supporting ReadGroup, non-supporting ReadGroup)
        """
        with self._registry_lock:
            # Use cached read groups if available
            if variant.ID in self._variant_read_groups:
                groups = self._variant_read_groups[variant.ID]
                support_names = groups['support']
                nonsupport_names = groups['nonsupport']
                
                # Create ReadGroups from cached AlignedRead objects
                support_reads = [self._read_registry[name] for name in support_names]
                nonsupport_reads = [self._read_registry[name] for name in nonsupport_names]
                
                return ReadGroup(support_reads), ReadGroup(nonsupport_reads)
            
            # If not cached, fetch and process reads
            all_reads, supporting_names = self._fetch_reads(variant)
            
            support_reads = []
            nonsupport_reads = []
            
            for read in all_reads:
                aligned_read = self._get_read(read)
                
                if read.query_name in supporting_names:
                    support_reads.append(aligned_read)
                else:
                    nonsupport_reads.append(aligned_read)
            
            # Cache read groups and mark registry as modified
            self._variant_read_groups[variant.ID] = {
                'support': [r.name for r in support_reads],
                'nonsupport': [r.name for r in nonsupport_reads]
            }
            self._registry_modified = True
            
            return ReadGroup(support_reads), ReadGroup(nonsupport_reads)
    
    def save_registry(self, force: bool = False) -> None:
        """Save the registry to disk if it has been modified or if forced.
        
        Args:
            force: If True, save regardless of modification status
        """
        if not self.registry_dir:
            self.logger.warning("No registry directory specified, cannot save registry")
            return
            
        with self._registry_lock:
            if not force and not self._registry_modified:
                self.logger.debug("Registry not modified since last save, skipping")
                return
                
            self._save_complete_registry()
            self._registry_modified = False
    
    def _save_complete_registry(self) -> None:
        """Save the complete registry to disk.
        Uses pydantic models for validation and serialization.
        """
        with self._registry_lock:
            self.logger.info("Saving registry to disk...")
            
            # Save registry metadata
            registry_metadata = RegistryMetadata(
                total_variants=len(self._variant_read_groups),
                total_reads=len(self._read_registry),
                bam_path=str(self.bam_path),
                window_size=self.window_size,
                processed_time=datetime.datetime.now().isoformat()
            )
            
            registry_metadata_path = self.registry_dir / "registry_metadata.json"
            with open(registry_metadata_path, 'w') as f:
                json.dump(registry_metadata.model_dump(), f)
            
            # Save variant read groups
            variant_groups_dict = {}
            for variant_id, groups in self._variant_read_groups.items():
                variant_groups_dict[variant_id] = VariantReadGroups(**groups).model_dump()
                
            variant_groups_path = self.registry_dir / "variant_read_groups.json"
            with open(variant_groups_path, 'w') as f:
                json.dump(variant_groups_dict, f)
            
            # Save read metadata
            read_metadata = {}
            for name, read in self._read_registry.items():
                read_model = ReadMetadata(
                    contig=read.contig,
                    start=read.start,
                    end=read.end,
                    mapq=read.mapq,
                    strand=read.strand,
                    cigar=read.cigar,
                    is_supplementary=read.is_supplementary,
                    is_secondary=read.is_secondary,
                    edit_distance=read.edit_distance,
                    soft_clip_left=read.soft_clip_left,
                    soft_clip_right=read.soft_clip_right,
                    cigar_stats=read.cigar_stats
                )
                read_metadata[name] = read_model.model_dump()
                
            reads_path = self.registry_dir / "read_metadata.json"
            with open(reads_path, 'w') as f:
                json.dump(read_metadata, f)
                
            self.logger.info(f"Registry saved to {self.registry_dir}")
    
    def load_registry(self, registry_dir: Path) -> None:
        """Load a previously saved registry from disk.
        
        Args:
            registry_dir: Directory containing registry files
        """
        if not registry_dir.exists():
            raise FileNotFoundError(f"Registry directory does not exist: {registry_dir}")
            
        self.logger.info(f"Loading registry from {registry_dir}")
        
        # Load and validate metadata
        metadata_path = registry_dir / "registry_metadata.json"
        if not metadata_path.exists():
            raise FileNotFoundError(f"Registry metadata file not found: {metadata_path}")
        
        with open(metadata_path, 'r') as f:
            metadata_dict = json.load(f)
            metadata = RegistryMetadata(**metadata_dict)
            self.logger.info(f"Loading registry with {metadata.total_variants} variants and {metadata.total_reads} reads")
        
        # Load variant read groups
        variant_groups_path = registry_dir / "variant_read_groups.json"
        if not variant_groups_path.exists():
            raise FileNotFoundError(f"Variant read groups file not found: {variant_groups_path}")
            
        with open(variant_groups_path, 'r') as f:
            variant_groups_dict = json.load(f)
            
            self._variant_read_groups = {}
            for variant_id, groups_data in variant_groups_dict.items():
                validated_groups = VariantReadGroups(**groups_data)
                self._variant_read_groups[variant_id] = validated_groups.model_dump()
        
        # Load read metadata and reconstruct AlignedRead objects
        reads_path = registry_dir / "read_metadata.json"
        if not reads_path.exists():
            raise FileNotFoundError(f"Read metadata file not found: {reads_path}")
            
        with open(reads_path, 'r') as f:
            read_metadata = json.load(f)
            
            self._read_registry = {}
            for read_name, read_data in read_metadata.items():
                try:
                    validated_metadata = ReadMetadata(**read_data)
                    
                    if self._bam is None:
                        self.logger.warning("BAM file not open, can't fetch read")
                        continue
                        
                    # Find and wrap the read
                    for read in self._bam.fetch(validated_metadata.contig, 
                                              validated_metadata.start, 
                                              validated_metadata.end):
                        if read.query_name == read_name:
                            aligned_read = AlignedRead.from_pysam_read(read)
                            self._read_registry[read_name] = aligned_read
                            break
                            
                except Exception as e:
                    self.logger.warning(f"Skipping invalid read metadata for {read_name}: {e}")
        
        self.logger.info("Registry loaded successfully")
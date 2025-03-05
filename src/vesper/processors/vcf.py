from dataclasses import dataclass
from typing import Iterator, Tuple, Optional
import pysam
from pathlib import Path
import logging

from ..models.variants import Variant, SVType, VariantAnalysis
from ..models.reads import ReadGroup, AlignedRead


class VCFProcessor:
    """Handles reading and writing VCF files."""
    
    def __init__(self, vcf_path: Path, test_mode: Optional[int] = None):
        self.vcf_path = vcf_path
        self._vcf: Optional[pysam.VariantFile] = None
        self.logger = logging.getLogger(__name__)
        self.test_mode = test_mode
    
    def __enter__(self):
        """Open VCF file and check/create index."""
        # Check for both possible index types
        tbi_path, csi_path = str(self.vcf_path) + '.tbi', str(self.vcf_path) + '.csi'
        
        if not Path(tbi_path).exists() and not Path(csi_path).exists():
            self.logger.info(f"Index not found for {self.vcf_path}, creating index...")
            try:
                pysam.tabix_index(str(self.vcf_path), preset='vcf', force=True)
                self.logger.info("Index created successfully")
            except Exception as e:
                self.logger.error(f"Failed to create VCF index: {e}")
                raise RuntimeError(f"Failed to create VCF index: {e}")
        
        self._vcf = pysam.VariantFile(str(self.vcf_path), 'r')
        self.logger.debug(f"Opened VCF file: {self.vcf_path}")
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._vcf:
            self._vcf.close()
            self.logger.debug(f"Closed VCF file: {self.vcf_path}")
            self._vcf = None
    
    def instantiate_variants(self) -> Iterator[VariantAnalysis]:
        """Create VariantAnalysis objects from records in a VCF file.
        
        Args:
            test_mode: If an integer, only loads that many records from VCF.
                      If None, loads all records.
        
        Yields:
            VariantAnalysis objects for each record in the VCF
            
        Raises:
            RuntimeError: if VCF file not opened
            ValueError: if variant record is invalid
        """
        if not self._vcf:
            raise RuntimeError("VCF file not opened. Use with-statement to open file.")
            
        for i, record in enumerate(self._vcf):
            if self.test_mode is not None and i >= self.test_mode:
                break
                
            try:
                base_variant = Variant.from_pysam_record(record)
                analysis = VariantAnalysis(variant=base_variant)
                yield analysis
            except Exception as e:
                raise type(e)(f"Invalid variant record: {str(e)}") from e
    
    @staticmethod
    def create_vcf(output_path: Path, template_vcf: Optional[Path] = None) -> pysam.VariantFile:
        """Create a new VCF file, optionally copying header from template.
        
        Args:
            output_path: Path to create new VCF
            template_vcf: Optional path to VCF to copy header from
            TODO: Default use input VCF header
            
        Returns:
            pysam.VariantFile opened for writing
        """
        mode = 'w' if str(output_path).endswith('.vcf') else 'wz'  # Compressed if .gz
        
        if template_vcf:
            # Copy header from template
            with pysam.VariantFile(str(template_vcf)) as vcf:
                header = vcf.header
        else:
            # Create minimal header
            # TODO: This needs to be updated with dependent annotation information per variant
            header = pysam.VariantHeader()
            header.add_line('##fileformat=VCFv4.2')
            header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
            header.add_line('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">')
            header.add_line('##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names">')
            
            # Add required fields
            header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
            header.contigs.add('1')  # Minimal contig
            
        return pysam.VariantFile(str(output_path), mode, header=header)
from dataclasses import dataclass
from typing import Iterator, Tuple, Optional, List, TextIO
import pysam
from pathlib import Path
import logging
from datetime import datetime
from io import StringIO
import gzip
import tempfile
import os

from ..models.variants import Variant, SVType, VariantAnalysis
from ..models.reads import ReadGroup, AlignedRead


class VCFProcessor:
    """Handles reading VCF files and creating Variant + VariantAnalysis objects."""
    
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

class VCFWriter:
    """Handles writing VCF files from VariantAnalysis objects."""
    
    def __init__(self, output_path: Path, compress: bool = True):
        """Initialize a VCF writer.
        
        Args:
            output_path: Path to create new VCF
            compress: Whether to bgzip compress output VCF (default True)
        """
        self.output_path = output_path
        self.compress = compress
        self.file = None
        
        if compress:
            if not str(output_path).endswith('.gz'):
                raise ValueError("Output path must end with .gz for compressed VCF files!")
            else:
                self.output_path = output_path
            self.temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.vcf')
            self.temp_path = Path(self.temp_file.name)
        else:
            self.temp_file = None
            self.temp_path = None
    
    def __enter__(self):
        """Open the VCF file for writing."""
        if self.compress:
            # Open the temporary file for writing
            self.file = self.temp_file
        else:
            # Open the output file directly
            self.file = open(self.output_path, 'w')
        
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close the VCF file."""
        if self.file:
            self.file.close()
            
        if self.compress and self.temp_path:
            pysam.tabix_compress(str(self.temp_path), str(self.output_path), force=True) # force overwrite
            os.unlink(self.temp_path)
    
    def create_header(self, variants: List[VariantAnalysis], template_vcf: Optional[Path] = None) -> str:
        """Create a VCF header string, optionally copying from a template (i.e. input VCF).
        
        Args:
            variants: List of Variant objects to derive contig information from
            template_vcf: Optional path to VCF to copy header from
            
        Returns:
            String representation of the VCF header
        """
        buffer = StringIO()
        
        if template_vcf:
            # Read header lines from template VCF
            # Check if the file is gzipped
            is_gzipped = str(template_vcf).endswith('.gz')
            
            # Open with appropriate opener
            opener = gzip.open if is_gzipped else open
            with opener(template_vcf, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        self.write_line(buffer, line.strip())
                    else:
                        break
        else:
            # Create a new header
            contigs = set(v.variant.chrom for v in variants)

            self.write_header_line(buffer, "fileformat=VCFv4.2")
            self.write_header_line(buffer, "fileDate=" + datetime.now().strftime("%Y%m%d"))
            self.write_header_line(buffer, "source=vesper")
            self.write_header_line(buffer, "reference=file:///seq/references/placeholder.fasta")
            for contig in contigs:
                # TODO: Get real lengths of the contigs from the reference/BAM file!
                self.write_header_line(buffer, f"contig=<ID={contig},length=10000000>")
            self.write_header_line(buffer, "FILTER=<ID=PASS,Description=\"All filters passed\">")  
            self.write_header_line(buffer, "FILTER=<ID=REP,Description=\"High identity RepeatMasker overlap/proximal feature\">")
            self.write_header_line(buffer, "FILTER=<ID=NO_PRIMARY_SUPPORT,Description=\"No primary supporting reads\">")
            self.write_header_line(buffer, "FILTER=<ID=IN_HIGH_SEGDUP,Description=\"Variant is in a high identity (>98%) segmental duplication\">")
            self.write_header_line(buffer, "FILTER=<ID=REPEAT_OVERLAP_MATCH,Description=\"Insertion RepeatMasker annotation matches overlapped RepeatMasker annotation\">")
            self.write_header_line(buffer, "FILTER=<ID=LOW_CONFIDENCE_OTHER,Description=\"Variant is likely low confidence (see CONFIDENCE_FLAGS for details)\">")
            self.write_header_line(buffer, "INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">")
            self.write_header_line(buffer, "INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">")
            self.write_header_line(buffer, "INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Supporting read names\">")
            self.write_header_line(buffer, "INFO=<ID=RMAPQ,Number=1,Type=Integer,Description=\"Supporting read mapping qualities\">")
            self.write_header_line(buffer, "INFO=<ID=RSFTCLIP,Number=1,Type=Float,Description=\"Supporting read soft-clipping stats\">")
            self.write_header_line(buffer, "INFO=<ID=NSMAPQ,Number=1,Type=Float,Description=\"Nonsupporting read mapping qualities\">")
            self.write_header_line(buffer, "INFO=<ID=NSFTCLIP,Number=1,Type=Float,Description=\"Nonsupporting read soft-clipping stats\">")
            self.write_header_line(buffer, "INFO=<ID=CONFIDENCE,Number=1,Type=Float,Description=\"Confidence score for the variant\">")
            self.write_header_line(buffer, "INFO=<ID=OVERLAPPING,Number=.,Type=String,Description=\"Overlapping annotated features\">")
            self.write_header_line(buffer, "INFO=<ID=PROXIMAL,Number=.,Type=String,Description=\"Proximal annotated features\">")
            self.write_header_line(buffer, "INFO=<ID=REPEATMASKER_RESULTS,Number=.,Type=String,Description=\"Repetitive sequence motifs identified via RepeatMasker (insertion only)\">")
            self.write_header_line(buffer, "FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"Number of reference-supporting reads\">")
            self.write_header_line(buffer, "FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of variant-supporting reads\">")
            self.write_line(buffer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE")
        
        return buffer.getvalue().rstrip()
    
    def write_header(self, variants: List[VariantAnalysis], template_vcf: Optional[Path] = None) -> None:
        """Write the VCF header to the file.
        
        Args:
            variants: List of Variant objects to derive contig information from
            template_vcf: Optional path to VCF to copy header from
        """
        if not self.file:
            raise RuntimeError("VCF file not opened. Use with-statement to open file.")
            
        header = self.create_header(variants, template_vcf)
        for line in header.split('\n'):
            self.write_line(self.file, line)
    
    def write_line(self, file: TextIO, line: str) -> None:
        """Write any line to a VCF file.
        
        Args:
            file: Open file object for writing
            line: Line to write
        """
        file.write(line + "\n")
    
    def write_header_line(self, file: TextIO, line: str) -> None:
        """Convenience function to write a header line to a VCF file.
        
        Args:
            file: Open file object for writing
            line: Header line content (without the ## prefix)
        """
        if not line.startswith("#"):
            line = "##" + line
        self.write_line(file, line)

    def write_record(self, variant: VariantAnalysis) -> None:
        """Write a single VariantAnalysis object to the VCF file.
        
        Args:
            variant: VariantAnalysis object to write
        """
        if not self.file:
            raise RuntimeError("VCF file not opened. Use with-statement to open file.")
            
        vcf_record = variant.to_vcf_record()
        self.write_line(self.file, vcf_record)
    
    @staticmethod
    def create_tabix_index(vcf_path: Path) -> None:
        """Create a tabix index for a VCF file.
        
        Args:
            vcf_path: Path to the VCF file to index (must be bgzipped)
        """
        if not str(vcf_path).endswith('.gz'):
            raise ValueError("VCF file must be bgzipped to create a tabix index!")

        pysam.tabix_index(str(vcf_path), preset="vcf", force=True)
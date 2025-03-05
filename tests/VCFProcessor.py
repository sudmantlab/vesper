import pytest
import tempfile
import shutil
from pathlib import Path
import pysam
from vesper.processors.vcf import VCFProcessor
from vesper.models.variants import SVType

TEST_DATA_DIR = Path("tests/test_data")

@pytest.fixture(scope="session", autouse=True)
def setup_test_data_dir():
    """Create and cleanup the test data directory."""
    TEST_DATA_DIR.mkdir(parents=True, exist_ok=True)
    yield
    shutil.rmtree(TEST_DATA_DIR)

def create_bgzipped_vcf(content: str, filename: str) -> Path:
    """Helper to create a bgzipped VCF file with index."""
    # First create the raw VCF
    vcf_path = TEST_DATA_DIR / filename
    vcf_path.write_text(content)
    
    # BGZip it
    bgzip_path = TEST_DATA_DIR / f"{filename}.gz"
    pysam.tabix_compress(str(vcf_path), str(bgzip_path), force=True)
    
    # Index it
    pysam.tabix_index(str(bgzip_path), preset="vcf", force=True)
    
    # Remove the raw VCF
    vcf_path.unlink()
    
    return bgzip_path

@pytest.fixture
def minimal_vcf():
    """Create a minimal valid VCF file."""
    content = """##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chrZ>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t1000\tDEL_1000\tA\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=100\tGT:GQ:DR:DV\t0/1:60:10:5
chr1\t2000\tDUP_2000\tG\t<DUP>\t45\tPASS\tSVTYPE=DUP;SVLEN=200\tGT:GQ:DR:DV\t0/1:45:12:6
chr1\t5000\tINS_5000\tA\t<INS>\t45\tPASS\tSVTYPE=INS;SVLEN=100\tGT:GQ:DR:DV\t0/1:45:12:6
"""
    return create_bgzipped_vcf(content, "minimal.vcf")

@pytest.fixture
def broken_position_vcf():
    """Create a VCF with invalid position."""
    content = """##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chrZ>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t-5\tBROKEN_POS\tA\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=100\tGT:GQ:DR:DV\t0/1:60:10:5
"""
    return create_bgzipped_vcf(content, "broken_position.vcf")

@pytest.fixture
def broken_svtype_vcf():
    """Create a VCF with invalid SVTYPE."""
    content = """##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chrZ>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chrZ\t1000\tBROKEN_SVTYPE\tA\t<BROKEN>\t45\tPASS\tSVTYPE=BROKEN;SVLEN=100\tGT:GQ:DR:DV\t0/1:45:12:6
"""
    return create_bgzipped_vcf(content, "broken_svtype.vcf")

@pytest.fixture
def negative_svlen_vcf():
    """Create a VCF with negative SVLEN."""
    content = """##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chrZ>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chrZ\t1000\tNEG_SVLEN\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100\tGT:GQ:DR:DV\t0/1:45:12:6
"""
    return create_bgzipped_vcf(content, "negative_svlen.vcf")

@pytest.fixture
def empty_vcf():
    """Create an empty VCF file."""
    content = """##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chrZ>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""
    return create_bgzipped_vcf(content, "empty.vcf")

def test_read_minimal_vcf(minimal_vcf):
    """Test reading a minimal valid VCF file."""
    with VCFProcessor(minimal_vcf) as vcf_proc:
        variants = list(vcf_proc.instantiate_variants())
        
        assert len(variants) == 3
        
        # Check first variant
        v1 = variants[0].variant
        assert v1.ID == "DEL_1000"
        assert v1.contig == "chr1"
        assert v1.position == 1000
        assert v1.sv_type == SVType.DEL
        assert v1.sv_length == 100
        assert v1.DR == 10
        assert v1.DV == 5
        assert v1.quality == 60
        
        # Check second variant
        v2 = variants[1].variant
        assert v2.ID == "DUP_2000"
        assert v2.sv_type == SVType.DUP
        assert v2.sv_length == 200
        
        # Check third variant
        v3 = variants[2].variant
        assert v3.ID == "INS_5000"
        assert v3.sv_type == SVType.INS
        assert v3.sv_length == 100

def test_broken_position(broken_position_vcf):
    """Test that negative position raises OSError."""
    with pytest.raises(OSError):
        with VCFProcessor(broken_position_vcf) as vcf_proc:
            list(vcf_proc.instantiate_variants())

def test_broken_svtype(broken_svtype_vcf):
    """Test that invalid SVTYPE raises ValueError."""
    with pytest.raises(ValueError):
        with VCFProcessor(broken_svtype_vcf) as vcf_proc:
            list(vcf_proc.instantiate_variants())

def test_negative_svlen(negative_svlen_vcf):
    """Test that negative SVLEN is handled via absolute value."""
    with VCFProcessor(negative_svlen_vcf) as vcf_proc:
        variants = list(vcf_proc.instantiate_variants())
        assert len(variants) == 1
        assert variants[0].variant.sv_length == 100  # Should be absolute value

def test_empty_vcf(empty_vcf):
    """Test that empty VCF returns no variants."""
    with VCFProcessor(empty_vcf) as vcf_proc:
        variants = list(vcf_proc.instantiate_variants())
        assert len(variants) == 0 
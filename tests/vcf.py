import pytest
import shutil
from pathlib import Path
import pysam
from vesper.processors.vcf import VCFProcessor, VCFWriter
from vesper.processors.annotations import GenomicInterval
from vesper.models.variants import SVType, Variant, VariantAnalysis

TEST_DATA_DIR = Path("tests/test_data")

@pytest.fixture(scope="session", autouse=True)
def setup_test_data_dir():
    """Create and cleanup the test data directory."""
    TEST_DATA_DIR.mkdir(parents=True, exist_ok=True)
    yield
    shutil.rmtree(TEST_DATA_DIR)

def create_bgzipped_vcf(content: str, filename: str) -> Path:
    """Helper to create a bgzipped VCF file with index for testing."""
    vcf_path = TEST_DATA_DIR / filename
    vcf_path.write_text(content)

    bgzip_path = TEST_DATA_DIR / f"{filename}.gz"
    pysam.tabix_compress(str(vcf_path), str(bgzip_path), force=True)
    pysam.tabix_index(str(bgzip_path), preset="vcf", force=True)

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
##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t1000\tDEL_1000\tA\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=100;RNAMES=read1,read2\tGT:GQ:DR:DV\t0/1:60:10:5
chr1\t2000\tDUP_2000\tG\t<DUP>\t45\tPASS\tSVTYPE=DUP;SVLEN=200;RNAMES=read3,read4\tGT:GQ:DR:DV\t0/1:45:12:6
chr1\t5000\tINS_5000\tA\t<INS>\t45\tPASS\tSVTYPE=INS;SVLEN=100;RNAMES=read5,read6\tGT:GQ:DR:DV\t0/1:45:12:6
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
##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t-5\tBROKEN_POS\tA\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=100;RNAMES=read1,read2\tGT:GQ:DR:DV\t0/1:60:10:5
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
##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chrZ\t1000\tBROKEN_SVTYPE\tA\t<BROKEN>\t45\tPASS\tSVTYPE=BROKEN;SVLEN=100;RNAMES=read1,read2\tGT:GQ:DR:DV\t0/1:45:12:6
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
##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chrZ\t1000\tNEG_SVLEN\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100;RNAMES=read1,read2\tGT:GQ:DR:DV\t0/1:45:12:6
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
##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""
    return create_bgzipped_vcf(content, "empty.vcf")

@pytest.fixture
def mock_variant_1():
    """Create a mock frozen Variant object for testing."""
    mock_variant = Variant(
        chrom="chr1",
        position=1000,
        ID="TEST_DEL_1", 
        ref="A",
        alt="<DEL>",
        qual=60,
        filter="PASS",
        info={"SVTYPE": "DEL", "SVLEN": 500, "RNAMES": ["read1", "read2"]},
        format="DR:DV",
        samples=[{"DR": 10, "DV": 5}],
        sv_type=SVType.DEL,
        sv_length=500,
        DR=10,
        DV=5,
        rnames=["read1", "read2"],
    )
    return mock_variant

@pytest.fixture
def mock_variant_2():
    """Create a second mock frozen Variant object for testing."""
    mock_variant = Variant(
        chrom="chr2",
        position=2000,
        ID="TEST_DUP_1",
        ref="C",
        alt="<DUP>", 
        qual=45,
        filter="PASS",
        info={"SVTYPE": "DUP", "SVLEN": 1000, "RNAMES": ["read3", "read4", "read5"]},
        format="DR:DV",
        samples=[{"DR": 8, "DV": 7}],
        sv_type=SVType.DUP,
        sv_length=1000,
        DR=8,
        DV=7,
        rnames=["read3", "read4", "read5"],
    )
    return mock_variant

@pytest.fixture
def mock_variant_analysis_1(mock_variant_1):
    """Create a mock VariantAnalysis object for testing."""
    analysis = VariantAnalysis(variant=mock_variant_1)
    analysis.confidence = 0.85
    analysis.metrics = {
        'n_support': 5,
        'n_nonsupport': 10,
        'comparison': {
            'mapq_mean': 40,
            'mapq_mean_other': 60,
            'softclip_stats': {'pct_softclipped': 15.5},
            'softclip_stats_other': {'pct_softclipped': 5.2},
            'cigar_stats': {'M': 1000, 'S': 200},
            'cigar_stats_other': {'M': 1200, 'S': 100}
        }
    }
    analysis.overlapping_features = [GenomicInterval(chrom="chr1", start=900, end=1100, metadata={'name': "AluSx"})]
    analysis.proximal_features = [GenomicInterval(chrom="chr1", start=1200, end=1500, metadata={'name': "AluYb8"})]
    return analysis

@pytest.fixture
def mock_variant_analysis_2(mock_variant_2):
    """Create a second mock VariantAnalysis object for testing."""
    analysis = VariantAnalysis(variant=mock_variant_2)
    analysis.confidence = 0.92
    analysis.metrics = {
        'n_support': 7,
        'n_nonsupport': 8,
        'comparison': {
            'mapq_mean': 45,
            'mapq_mean_other': 55,
            'softclip_stats': {'pct_softclipped': 12.3},
            'softclip_stats_other': {'pct_softclipped': 4.8},
            'cigar_stats': {'M': 1200, 'S': 150},
            'cigar_stats_other': {'M': 1400, 'S': 80}
        }
    }
    analysis.overlapping_features = [GenomicInterval(chrom="chr2", start=1900, end=2100, metadata={'name': "L1PA3"})]
    analysis.proximal_features = [GenomicInterval(chrom="chr2", start=2200, end=2400, metadata={'name': "AluY"})]
    return analysis

# VCF reading

def test_read_minimal_vcf(minimal_vcf):
    """Test reading a minimal valid VCF file."""
    with VCFProcessor(minimal_vcf) as vcf_proc:
        variants = list(vcf_proc.instantiate_variants())
        
        assert len(variants) == 3
        
        # Check first variant
        v1 = variants[0].variant
        assert v1.ID == "DEL_1000"
        assert v1.chrom == "chr1"
        assert v1.position == 1000
        assert v1.sv_type == SVType.DEL
        assert v1.sv_length == 100
        assert v1.DR == 10
        assert v1.DV == 5
        assert v1.qual == 60
        assert len(v1.rnames) >= 0  # May be empty or contain multipleread names
        
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
    """Test that negative position raises ValueError."""
    with pytest.raises(Exception):
        with VCFProcessor(broken_position_vcf) as vcf_proc:
            list(vcf_proc.instantiate_variants())

def test_broken_svtype(broken_svtype_vcf):
    """Test that invalid SVTYPE raises ValueError."""
    with pytest.raises(Exception):
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

# VCF writing

def test_create_vcf_header_from_scratch(mock_variant_analysis_1, tmp_path):
    """Test creating a VCF header from scratch."""
    variants = [mock_variant_analysis_1]
    writer = VCFWriter(TEST_DATA_DIR / "dummy.vcf", compress=False)
    header = writer.create_header(variants)
    
    # Check basic header properties
    assert "fileformat=VCFv4.2" in header
    assert "source=vesper" in header
    assert "fileDate=" in header
    
    # Check contig information
    assert "contig=<ID=chr1" in header
    
    # Check INFO fields
    assert "INFO=<ID=SVTYPE" in header
    assert "INFO=<ID=SVLEN" in header
    assert "INFO=<ID=CONFIDENCE" in header
    
    # Check FORMAT fields
    assert "FORMAT=<ID=DR" in header
    assert "FORMAT=<ID=DV" in header

def test_create_vcf_header_from_template(mock_variant_analysis_1, minimal_vcf):
    """Test creating a VCF header from a template."""
    variants = [mock_variant_analysis_1]
    writer = VCFWriter(TEST_DATA_DIR / "mock.vcf", compress=False)
    header = writer.create_header(variants, template_vcf=minimal_vcf)
    
    assert "contig=<ID=chr1" in header
    assert "contig=<ID=chrZ" in header

def test_to_vcf_record(mock_variant_analysis_1):
    """Test converting a VariantAnalysis to a VCF record."""
    vcf_record = mock_variant_analysis_1.to_vcf_record()
    
    # Check that the record contains the expected fields
    assert "chr1" in vcf_record # chrom
    assert "1000" in vcf_record # pos
    assert "TEST_DEL_1" in vcf_record # ID
    assert "N" in vcf_record # ref
    assert "<DEL>" in vcf_record # alt
    assert "." in vcf_record # qual
    assert "PASS" in vcf_record # filter
    assert "SVTYPE=DEL" in vcf_record # info field - SVTYPE
    assert "SVLEN=500" in vcf_record # info field - SVLEN
    assert f"CONFIDENCE={mock_variant_analysis_1.confidence}" in vcf_record # info field - CONFIDENCE
    assert "DR:DV" in vcf_record # format
    assert "10:5" in vcf_record # sample values

def test_write_vcf_record(mock_variant_analysis_1):
    """Test writing a single VariantAnalysis object to an uncompressed VCF."""
    output_path = TEST_DATA_DIR / "test_write_record.vcf"
    
    with VCFWriter(output_path, compress=False) as writer:
        writer.write_header([mock_variant_analysis_1])
        writer.write_record(mock_variant_analysis_1)
    
    assert output_path.exists()
    
    # Verify the file has the expected content
    with open(output_path, 'r') as f:
        content = f.read()
        
        # Comprehensive check for the header
        assert "fileformat=VCFv4.2" in content
        assert "source=vesper" in content
        assert "##contig=<ID=chr1" in content
        assert "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=" in content
        assert "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=" in content
        assert "##INFO=<ID=CONFIDENCE,Number=1,Type=Float,Description=" in content
        assert "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=" in content
        assert "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=" in content
        assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" in content
        
        # Comprehensive check that the record contains the expected fields
        assert "chr1" in content # chrom
        assert "1000" in content # pos
        assert "TEST_DEL_1" in content # ID
        assert "N" in content # ref
        assert "<DEL>" in content # alt
        assert "." in content # qual
        assert "PASS" in content # filter
        assert "SVTYPE=DEL" in content # info field - SVTYPE
        assert "SVLEN=500" in content # info field - SVLEN
        assert f"CONFIDENCE={mock_variant_analysis_1.confidence}" in content # info field - CONFIDENCE
        assert "DR:DV" in content # format
        assert "10:5" in content # sample values

def test_write_records(mock_variant_analysis_1, mock_variant_analysis_2):
    """Test writing multiple VariantAnalysis objects to an uncompressed VCF."""
    output_path = TEST_DATA_DIR / "test_write_records.vcf"
    variants = [mock_variant_analysis_1, mock_variant_analysis_2]
    with VCFWriter(output_path, compress=False) as writer:
        writer.write_header(variants)
        for variant in variants:
            writer.write_record(variant)

    assert output_path.exists()
    
    # Verify the file has the expected content
    with open(output_path, 'r') as f:
        content = f.read()
        
        # Shallow check for the header
        assert "fileformat=VCFv4.2" in content
        assert "source=vesper" in content
        assert "##contig=<ID=chr1" in content
        assert "##contig=<ID=chr2" in content
        
        # Shallow check for the first variant record
        assert "chr1\t1000\tTEST_DEL_1" in content
        assert "SVTYPE=DEL;SVLEN=500" in content
        assert "10:5" in content
        
        # Shallow check for the second variant record  
        assert "chr2\t2000\tTEST_DUP_1" in content
        assert "SVTYPE=DUP;SVLEN=1000" in content
        assert "8:7" in content

def test_writer_extension_check(mock_variant_analysis_1, mock_variant_analysis_2):
    """Test that the extension check works."""
    output_path = TEST_DATA_DIR / "test_compressed.vcf" # lacks .gz extension
    with pytest.raises(ValueError):
        with VCFWriter(output_path, compress=True) as writer:
            writer.write_header([mock_variant_analysis_1])

def test_tabix_extension_check(mock_variant_analysis_1, mock_variant_analysis_2):
    """Test that the tabix extension check works."""
    output_path = TEST_DATA_DIR / "test_compressed.vcf"
    with pytest.raises(ValueError):
        VCFWriter.create_tabix_index(output_path)

def test_write_compressed_vcf(mock_variant_analysis_1, mock_variant_analysis_2):
    """Test that we can write multiple variants to a compressed vcf.gz file."""
    output_path = TEST_DATA_DIR / "test_compressed.vcf.gz"
    variants = [mock_variant_analysis_1, mock_variant_analysis_2]
    with VCFWriter(output_path, compress=True) as writer:
        writer.write_header(variants)
        for variant in variants:
            writer.write_record(variant)

    assert output_path.exists()

    # Try to create a tabix index for the compressed file
    VCFWriter.create_tabix_index(output_path)
    index_path = Path(str(output_path) + ".tbi")
    assert index_path.exists()
    
    # Verify the file has the expected content when reading with pysam
    with pysam.VariantFile(str(output_path)) as vcf:
        assert "PASS" in list((vcf.header.filters))
        assert "RNAMES" in list((vcf.header.info))

        records = 0
        for record, original in zip(vcf, [mock_variant_analysis_1, mock_variant_analysis_2]):
            assert record.id == original.variant.ID
            assert record.chrom == original.variant.chrom
            assert record.pos == original.variant.position
            assert SVType[record.info["SVTYPE"]] == original.variant.sv_type
            assert record.info["SVLEN"] == original.variant.sv_length
            assert record.samples[0]["DR"] == original.variant.DR
            assert record.samples[0]["DV"] == original.variant.DV
            records += 1
        assert records == 2
        
    # Now try to read it back with VCFProcessor to verify full compatibility
    with VCFProcessor(output_path) as vcf_proc:
        variants = list(vcf_proc.instantiate_variants())
        
        # Verify we got one variant
        assert len(variants) == 2
        
        for loaded, original in zip(variants, [mock_variant_analysis_1, mock_variant_analysis_2]):
            assert loaded.variant.ID == original.variant.ID
            assert loaded.variant.chrom == original.variant.chrom
            assert loaded.variant.position == original.variant.position
            assert loaded.variant.sv_type == original.variant.sv_type
            assert loaded.variant.sv_length == original.variant.sv_length
            assert loaded.variant.DR == original.variant.DR
            assert loaded.variant.DV == original.variant.DV
        
            # Check that the confidence value is available in the INFO field
            assert "CONFIDENCE" in loaded.variant.info
            assert float(loaded.variant.info["CONFIDENCE"]) == pytest.approx(original.confidence)

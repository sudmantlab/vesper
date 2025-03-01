import pytest
from pathlib import Path
from io import StringIO
from vesper.processors.variants import VCFProcessor, Variant
from vesper.processors.annotations import BEDProcessor, GenomicInterval
from vesper.models.variants import SVType, VariantAnalysis

@pytest.fixture(scope="session", autouse=True)
def test_files_dir():
    """Create and manage the test files directory."""
    path = Path("tests/files")
    path.mkdir(exist_ok=True)
    return path

### vcf tests

@pytest.fixture
def minimal_vcf(test_files_dir):
    """Create a basic VCF with simple valid content"""
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
    
    path = test_files_dir / "minimal.vcf"
    path.write_text(content)
    yield path
    path.unlink(missing_ok=True)

@pytest.fixture
def broken_position_vcf(test_files_dir):
    """Create a VCF with invalid position"""
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
    
    path = test_files_dir / "broken_position.vcf"
    path.write_text(content)
    yield path
    path.unlink(missing_ok=True)

@pytest.fixture
def broken_svtype_vcf(test_files_dir):
    """Create a VCF with invalid SVTYPE"""
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
    
    path = test_files_dir / "broken_svtype.vcf"
    path.write_text(content)
    yield path
    path.unlink(missing_ok=True)

@pytest.fixture
def negative_svlen_vcf(test_files_dir):
    """Create a VCF with negative SVLEN"""
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
    
    path = test_files_dir / "negative_svlen.vcf"
    path.write_text(content)
    yield path
    path.unlink(missing_ok=True)

@pytest.fixture
def empty_vcf(test_files_dir):
    """Create an empty VCF with just headers"""
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
    
    path = test_files_dir / "empty.vcf"
    path.write_text(content)
    yield path
    path.unlink(missing_ok=True)

def test_minimal_vcf_reading(minimal_vcf):
    """
    Validate that we can read a simple VCF without crashing out
    Note that VCFs are 1-based.
    """
    with VCFProcessor(minimal_vcf) as vcf:
        variants = list(vcf.read_variants())
        assert len(variants) == 3
        
        v = variants[0]
        assert v.contig == "chr1"
        assert v.position == 1000
        assert v.sv_type == SVType.DEL
        assert v.sv_length == 100
        
        v = variants[1]
        assert v.contig == "chr1"
        assert v.position == 2000
        assert v.sv_type == SVType.DUP
        assert v.sv_length == 200

        v = variants[2]
        assert v.contig == "chr1"
        assert v.position == 5000
        assert v.sv_type == SVType.INS
        assert v.sv_length == 100

def test_invalid_position_handling(broken_position_vcf):
    """Check failures from broken VCF upon pysam parse"""
    with pytest.raises(OSError):
        with VCFProcessor(broken_position_vcf) as vcf:
            variants = list(vcf.read_variants())

def test_invalid_svtype_handling(broken_svtype_vcf):
    """Check failures from broken VCF upon pysam parse"""
    with pytest.raises(ValueError):
        with VCFProcessor(broken_svtype_vcf) as vcf:
            variants = list(vcf.read_variants())

def test_negative_svlen_handling(negative_svlen_vcf):
    """Check appropriate handling for negative SVLEN"""
    with VCFProcessor(negative_svlen_vcf) as vcf:
        variants = list(vcf.read_variants())
        assert len(variants) == 1
        v = variants[0]
        assert v.sv_length == 100 # appropriate absolute value conversion
        assert v.sv_type == SVType.DEL

def test_empty_vcf(empty_vcf):
    """Test handling of empty VCF file"""
    with VCFProcessor(empty_vcf) as vcf:
        variants = list(vcf.read_variants())
        assert len(variants) == 0

### bed tests
@pytest.fixture
def minimal_bed(test_files_dir):
    """Create a small BED file for testing."""
    content = """
    chr1	900	1100	repeat1	1	+
    chr1	999	1000	SNP1	1	+
    chr1	1499	1500	SNP2	1	+
    chr1	2000	2500	repeat2	1	-
    """
    
    path = test_files_dir / "tiny.bed"
    path.write_text(content)
    yield path
    path.unlink(missing_ok=True)

def test_minimal_bed_reading(minimal_bed):
    """Validate that we can read a simple BED file without crashing out"""
    processor = BEDProcessor(minimal_bed)
    processor.load()
    total_intervals = sum(len(tree) for tree in processor._intervals.values())
    assert total_intervals == 4

def test_bed_overlap(minimal_bed):
    """
    Test BED self-loading and simple overlap behavior
    """
    processor = BEDProcessor(minimal_bed) # test self-loading
    overlaps = processor.find_overlaps("chr1", 900, 1100) # should self-load without .load() call
    assert len(overlaps) == 2

def test_point_bed_overlap(minimal_bed):
    """
    Test BED point overlap behavior
    Note internal conversion of BED to 1-based coordinates and handling of point records.
    """
    processor = BEDProcessor(minimal_bed)

    # correct SNP1 point query (overlaps repeat1 and SNP1)
    overlaps = processor.find_overlaps("chr1", 1000, 1000)
    assert len(overlaps) == 2

    # excludes SNP2 (exists at 1500) - malformed point query
    overlaps = processor.find_overlaps("chr1", 1499, 1500)
    assert len(overlaps) == 0

    # excludes SNP2 (exists at 1500) - off by one point query
    overlaps = processor.find_overlaps("chr1", 1499, 1499)
    assert len(overlaps) == 0

    # correct SNP2 only point query
    overlaps = processor.find_overlaps("chr1", 1500, 1500)
    assert len(overlaps) == 1
    
def test_no_bed_overlap(minimal_bed):
    """Test return on no overlaps"""
    processor = BEDProcessor(minimal_bed)
    overlaps = processor.find_overlaps("chr1", 3000, 3500)
    assert len(overlaps) == 0

def test_variant_annotation_overlap(minimal_vcf, minimal_bed):
    """Test exact overlap behavior between variants and annotations"""
    vcf_proc = VCFProcessor(minimal_vcf)
    bed_proc = BEDProcessor(minimal_bed)
    bed_proc.load()
    
    with vcf_proc as vcf:
        variants = list(vcf.read_variants())
        first_var, second_var, third_var = variants[0], variants[1], variants[2]
        
        # test 1000-1100 overlap
        overlaps = bed_proc.find_overlaps(
            first_var.contig,
            first_var.position,
            first_var.end
        )
        print(overlaps)
        assert len(overlaps) == 2

        # test 2000-2500 overlap
        overlaps = bed_proc.find_overlaps(
            second_var.contig,
            second_var.position,
            second_var.end
        )
        assert len(overlaps) == 1

        # test 5000-5100 overlap
        overlaps = bed_proc.find_overlaps(
            third_var.contig,
            third_var.position,
            third_var.end
        )
        assert len(overlaps) == 0

### VariantAnalysis tests

def test_variant_analysis_creation(minimal_vcf):
    """Test creation of VariantAnalysis objects linked to Variants"""
    
    with VCFProcessor(minimal_vcf) as vcf:
        variants = list(vcf.read_variants())
        
        # Test basic creation
        analysis = VariantAnalysis(variants[0])
        assert analysis.variant == variants[0]
        assert analysis.variant.sv_type == SVType.DEL
        assert analysis.variant.position == 1000
        assert len(analysis.overlapping_features) == 0
        assert len(analysis.proximal_features) == 0
        assert analysis.confidence == 0.0

def test_variant_analysis_annotation(minimal_vcf, minimal_bed):
    """Test adding annotations to VariantAnalysis objects"""
    
    with VCFProcessor(minimal_vcf) as vcf:
        variants = list(vcf.read_variants())
        first_var = variants[0]  # DEL at position 1000
        
        analysis = VariantAnalysis(first_var)
        bed_proc = BEDProcessor(minimal_bed)
        
        analysis.add_annotations(bed_proc)
        
        # Should overlap repeat1 and SNP1
        assert len(analysis.overlapping_features) == 2
        assert any(f.name == "repeat1" for f in analysis.overlapping_features)
        assert any(f.name == "SNP1" for f in analysis.overlapping_features)
        
        # Should have SNP2 as proximal (within default 500bp span)
        assert len(analysis.proximal_features) == 1
        assert analysis.proximal_features[0].name == "SNP2"

def test_variant_analysis_custom_span(minimal_vcf, minimal_bed):
    """Test VariantAnalysis with custom proximal span distances"""
    from vesper.models.variants import VariantAnalysis
    
    with VCFProcessor(minimal_vcf) as vcf:
        variants = list(vcf.read_variants())
        second_var = variants[1]  # DUP at position 2000
        
        analysis = VariantAnalysis(second_var)
        bed_proc = BEDProcessor(minimal_bed)
        
        analysis.add_annotations(bed_proc, proximal_span=1000)
        
        # Should have 1 overlapping feature (repeat2)
        assert len(analysis.overlapping_features) == 1
        assert analysis.overlapping_features[0].name == "repeat2"
        
        # Should have repeat1, SNP1, and SNP2 as proximal (within 1000bp non-overlapping span)
        proximal_names = {f.name for f in analysis.proximal_features}
        assert len(proximal_names) == 3
        assert proximal_names == {"repeat1", "SNP1", "SNP2"}

def test_variant_analysis_no_annotations(minimal_vcf, minimal_bed):
    """Test VariantAnalysis behavior with no overlapping or proximal features"""
    from vesper.models.variants import VariantAnalysis
    
    with VCFProcessor(minimal_vcf) as vcf:
        variants = list(vcf.read_variants())
        third_var = variants[2]  # INS at position 5000
        
        analysis = VariantAnalysis(third_var)
        bed_proc = BEDProcessor(minimal_bed)
        
        # Add annotations
        analysis.add_annotations(bed_proc)
        
        assert len(analysis.overlapping_features) == 0
        assert len(analysis.proximal_features) == 0
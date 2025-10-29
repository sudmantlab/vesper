import pytest
import tempfile
from pathlib import Path
import shutil
import pysam
import pysam.bcftools
import json
from vesper.processors.reads import ReadProcessor
from vesper.processors.vcf import VCFProcessor
from vesper.models.registry import RegistryMetadata

SOURCE_DATA_DIR = Path("tests/files/hg38")
TEST_DATA_DIR = Path("tests/test_data")

@pytest.fixture(scope="session", autouse=True)
def setup_test_data_dir():
    """Create and cleanup the test data directory."""
    TEST_DATA_DIR.mkdir(parents=True, exist_ok=True)
    yield
    shutil.rmtree(TEST_DATA_DIR)

@pytest.fixture(scope="session")
def test_region():
    """Define a small genomic region on chr21 for testing."""
    return {
        'chrom': 'chr21',
        'start': 10600000,
        'end': 10700000
    }

@pytest.fixture(scope="session")
def small_bam(test_region) -> Path:
    """Create a small BAM file from chr21."""
    source_bam = SOURCE_DATA_DIR / "894.duplomap.bam"
    if not source_bam.exists():
        raise FileNotFoundError(f"Source BAM file not found: {source_bam}")
    
    # Ensure test directory exists
    TEST_DATA_DIR.mkdir(parents=True, exist_ok=True)
    
    bam_path = TEST_DATA_DIR / "sample.bam"
    if not bam_path.exists():
        region = f"{test_region['chrom']}:{test_region['start']}-{test_region['end']}"
        pysam.view(
            str(source_bam),
            "-b",
            region,
            "-o", str(bam_path),
            catch_stdout=False
        )
        
        pysam.index(str(bam_path), catch_stdout=False)
        
        if not bam_path.exists():
            raise RuntimeError(f"Failed to create BAM file: {bam_path}")
        if not Path(str(bam_path) + '.bai').exists():
            raise RuntimeError(f"Failed to create BAM index: {bam_path}.bai")
    
    return bam_path

@pytest.fixture(scope="session")
def small_vcf(test_region) -> Path:
    """Create a small VCF file from chr21."""
    source_vcf = SOURCE_DATA_DIR / "894.duplomap.vcf.gz"
    if not source_vcf.exists():
        raise FileNotFoundError(f"Source VCF file not found: {source_vcf}")
    
    # Ensure test directory exists
    TEST_DATA_DIR.mkdir(parents=True, exist_ok=True)
    
    vcf_path = TEST_DATA_DIR / "sample.vcf.gz"
    if not vcf_path.exists():
        # Extract region directly to compressed VCF
        region = f"{test_region['chrom']}:{test_region['start']}-{test_region['end']}"
        try:
            pysam.bcftools.view(
                str(source_vcf),
                "-r", region,
                "-O", "z",  # Output compressed VCF
                "-o", str(vcf_path),
                "--write-index",
                catch_stdout=False
            )
            
            # Verify the file was created
            if not vcf_path.exists():
                raise RuntimeError(f"Failed to create VCF file: {vcf_path}")
            if not Path(str(vcf_path) + '.csi').exists():
                raise RuntimeError(f"Failed to create VCF index: {vcf_path}.csi")
                
        except Exception as e:
            # Clean up any partial files
            vcf_path.unlink(missing_ok=True)
            Path(str(vcf_path) + '.csi').unlink(missing_ok=True)
            raise RuntimeError(f"Failed to process VCF: {e}")
    
    return vcf_path

@pytest.fixture
def registry_dir():
    """Create a temporary directory for registry files."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield Path(temp_dir)

def test_registry_creation(small_bam: Path, small_vcf: Path, registry_dir: Path):
    """Test basic registry creation and saving."""
    with ReadProcessor(small_bam, registry_dir=registry_dir) as proc:
        # Process a few reads to populate registry
        with VCFProcessor(small_vcf) as vcf_proc:
            variant = next(vcf_proc.instantiate_variants())
            proc.get_read_groups(variant.variant)
        
        # Verify registry files are created on exit
        proc.save_registry(force=True)
        
        assert (registry_dir / "registry_metadata.json").exists()
        assert (registry_dir / "variant_read_groups.json").exists()
        assert (registry_dir / "read_metadata.json").exists()
        
        # Verify metadata content
        with open(registry_dir / "registry_metadata.json") as f:
            metadata = RegistryMetadata(**json.load(f))
            assert metadata.bam_path == str(small_bam)
            assert metadata.total_variants > 0
            assert metadata.total_reads > 0

def test_registry_loading(small_bam: Path, small_vcf: Path, registry_dir: Path):
    """Test loading an existing registry."""
    # First create a registry
    variant_id = None
    read_count = 0
    
    with ReadProcessor(small_bam, registry_dir=registry_dir) as proc:
        with VCFProcessor(small_vcf) as vcf_proc:
            variant = next(vcf_proc.instantiate_variants())
            variant_id = variant.variant.ID
            support, nonsupport = proc.get_read_groups(variant.variant)
            read_count = len(support.reads) + len(nonsupport.reads)
        proc.save_registry(force=True)
    
    # Now load the registry and verify contents
    with ReadProcessor(small_bam, registry_dir=registry_dir) as proc:
        # Verify variant read groups are loaded
        assert len(proc._variant_read_groups) > 0
        assert variant_id in proc._variant_read_groups
        
        # Verify reads are loaded
        assert len(proc._read_registry) == read_count
        
        # Verify we can still fetch the same variant
        with VCFProcessor(small_vcf) as vcf_proc:
            variant = next(vcf_proc.instantiate_variants())
            support2, nonsupport2 = proc.get_read_groups(variant.variant)
            
            # Should match original counts
            assert len(support2.reads) + len(nonsupport2.reads) == read_count

def test_force_new_registry(small_bam: Path, small_vcf: Path, registry_dir: Path):
    """Test forcing creation of new registry ignores existing one."""
    # Create initial registry
    with ReadProcessor(small_bam, registry_dir=registry_dir) as proc:
        with VCFProcessor(small_vcf) as vcf_proc:
            variant = next(vcf_proc.instantiate_variants())
            proc.get_read_groups(variant.variant)
        proc.save_registry(force=True)
    
    # Get initial metadata
    with open(registry_dir / "registry_metadata.json") as f:
        initial_metadata = json.load(f)
    
    # Force new registry
    with ReadProcessor(small_bam, registry_dir=registry_dir, force_new_registry=True) as proc:
        with VCFProcessor(small_vcf) as vcf_proc:
            variant = next(vcf_proc.instantiate_variants())
            proc.get_read_groups(variant.variant)
        proc.save_registry(force=True)
    
    # Verify metadata has changed
    with open(registry_dir / "registry_metadata.json") as f:
        new_metadata = json.load(f)
        assert new_metadata['processed_time'] != initial_metadata['processed_time']

def test_auto_load_registry_disabled(small_bam: Path, small_vcf: Path, registry_dir: Path):
    """Test that disabling auto-load creates new registry."""
    # Create initial registry
    with ReadProcessor(small_bam, registry_dir=registry_dir) as proc:
        with VCFProcessor(small_vcf) as vcf_proc:
            variant = next(vcf_proc.instantiate_variants())
            proc.get_read_groups(variant.variant)
        proc.save_registry(force=True)
    
    # Get initial metadata
    with open(registry_dir / "registry_metadata.json") as f:
        initial_metadata = json.load(f)
    
    # Create new processor with auto-load disabled
    with ReadProcessor(small_bam, registry_dir=registry_dir, auto_load_registry=False) as proc:
        with VCFProcessor(small_vcf) as vcf_proc:
            variant = next(vcf_proc.instantiate_variants())
            proc.get_read_groups(variant.variant)
        proc.save_registry(force=True)
    
    # Verify metadata has changed
    with open(registry_dir / "registry_metadata.json") as f:
        new_metadata = json.load(f)
        assert new_metadata['processed_time'] != initial_metadata['processed_time']

def test_registry_modification_tracking(small_bam: Path, small_vcf: Path, registry_dir: Path):
    """Test that registry modification is properly tracked."""
    with ReadProcessor(small_bam, registry_dir=registry_dir) as proc:
        # Initially not modified
        assert not proc._registry_modified
        
        # Process variant should modify registry
        with VCFProcessor(small_vcf) as vcf_proc:
            variant = next(vcf_proc.instantiate_variants())
            proc.get_read_groups(variant.variant)
        
        assert proc._registry_modified
        
        # Save should clear modified flag
        proc.save_registry()
        assert not proc._registry_modified
        
        # Processing same variant shouldn't modify (using cache)
        with VCFProcessor(small_vcf) as vcf_proc:
            variant = next(vcf_proc.instantiate_variants())
            proc.get_read_groups(variant.variant)
        
        assert not proc._registry_modified

def test_registry_recovery_from_invalid(small_bam: Path, small_vcf: Path, registry_dir: Path):
    """Test recovery when loading invalid registry."""
    # Create corrupted registry files
    (registry_dir / "registry_metadata.json").write_text("invalid json")
    (registry_dir / "variant_read_groups.json").write_text("{}")
    (registry_dir / "read_metadata.json").write_text("{}")
    
    # Should handle gracefully and create new registry
    with ReadProcessor(small_bam, registry_dir=registry_dir) as proc:
        with VCFProcessor(small_vcf) as vcf_proc:
            variant = next(vcf_proc.instantiate_variants())
            support, nonsupport = proc.get_read_groups(variant.variant)
            
        assert len(support.reads) + len(nonsupport.reads) > 0
        
        # Should be able to save new valid registry
        proc.save_registry(force=True)
        assert (registry_dir / "registry_metadata.json").exists()
        
        # Verify new metadata is valid
        with open(registry_dir / "registry_metadata.json") as f:
            metadata = RegistryMetadata(**json.load(f))
            assert metadata.total_variants > 0

import pytest
import sqlite3
from pathlib import Path
import shutil
from vesper.processors.annotations import BEDProcessor, GenomicInterval
from vesper.models.variants import Variant, VariantAnalysis, SVType

TEST_DATA_DIR = Path("tests/test_data")

@pytest.fixture(scope="session", autouse=True)
def setup_test_data_dir():
    """Create and cleanup the test data directory."""
    TEST_DATA_DIR.mkdir(parents=True, exist_ok=True)
    yield
    shutil.rmtree(TEST_DATA_DIR)

@pytest.fixture
def simple_bed():
    """Create a simple BED file with a few features."""
    content = """chr1\t1000\t2000\tfeature1\t100\t+
chr1\t5000\t6000\tfeature2\t200\t-
chr2\t1000\t3000\tfeature3\t150\t+
"""
    bed_path = TEST_DATA_DIR / "simple.bed"
    bed_path.write_text(content)
    yield bed_path
    # Cleanup any associated SQLite DB
    db_path = bed_path.with_suffix(bed_path.suffix + '.sqlite')
    if db_path.exists():
        db_path.unlink()

@pytest.fixture
def overlapping_features_bed():
    """Create a BED file with overlapping features."""
    content = """chr1\t1000\t3000\tfeature1\t100\t+
chr1\t2000\t4000\tfeature2\t200\t-
chr1\t2500\t3500\tfeature3\t150\t+
"""
    bed_path = TEST_DATA_DIR / "overlapping.bed"
    bed_path.write_text(content)
    yield bed_path
    # Cleanup any associated SQLite DB
    db_path = bed_path.with_suffix(bed_path.suffix + '.sqlite')
    if db_path.exists():
        db_path.unlink()

@pytest.fixture
def malformed_bed():
    """Create a malformed BED file to test error handling."""
    content = """chr1\t1000\tabc\tfeature1\t100\t+
chr1\t5000\t4000\tfeature2\t200\t-
"""
    bed_path = TEST_DATA_DIR / "malformed.bed"
    bed_path.write_text(content)
    yield bed_path
    # Cleanup any associated SQLite DB
    db_path = bed_path.with_suffix(bed_path.suffix + '.sqlite')
    if db_path.exists():
        db_path.unlink()

@pytest.fixture
def test_variants():
    """Create test variants for annotation."""
    return [
        VariantAnalysis(
            variant=Variant(
                ID="test1",
                chrom="chr1",
                position=2500,  # Should overlap with features
                sv_type=SVType.DEL,
                sv_length=100,
                DR=10,
                DV=5
            )
        ),
        VariantAnalysis(
            variant=Variant(
                ID="test2",
                chrom="chr1",
                position=4500,  # Should be proximal to feature2
                sv_type=SVType.INS,
                sv_length=50,
                DR=8,
                DV=4
            )
        ),
        VariantAnalysis(
            variant=Variant(
                ID="test3",
                chrom="chr2",
                position=5000,  # Should have no overlaps or proximal features
                sv_type=SVType.DUP,
                sv_length=200,
                DR=12,
                DV=6
            )
        )
    ]

def test_build_sqlite_db(simple_bed):
    """Test building SQLite DB from BED file."""
    with BEDProcessor(simple_bed) as bed_proc:
        # DB should be created automatically
        db_path = simple_bed.with_suffix(simple_bed.suffix + '.sqlite')
        assert db_path.exists()
        
        # Check DB contents
        conn = sqlite3.connect(str(db_path))
        cursor = conn.execute("SELECT COUNT(*) FROM intervals")
        assert cursor.fetchone()[0] == 3
        
        # Check schema
        cursor = conn.execute("PRAGMA table_info(intervals)")
        columns = {row[1] for row in cursor.fetchall()}
        expected_columns = {'chrom', 'start', 'end', 'name', 'score', 'strand'}
        assert columns == expected_columns
        
        conn.close()

def test_rebuild_db(simple_bed):
    """Test rebuilding DB when requested."""
    with BEDProcessor(simple_bed) as bed_proc:
        # First build
        db_path = simple_bed.with_suffix(simple_bed.suffix + '.sqlite')
        assert db_path.exists()
        
        # Modify DB to test rebuild
        conn = sqlite3.connect(str(db_path))
        conn.execute("DELETE FROM intervals")
        conn.commit()
        conn.close()
        
        # Rebuild
        bed_proc._build_sqlite_db(rebuild=True)
        
        # Check contents restored
        conn = sqlite3.connect(str(db_path))
        cursor = conn.execute("SELECT COUNT(*) FROM intervals")
        assert cursor.fetchone()[0] == 3
        conn.close()

def test_find_overlaps(overlapping_features_bed):
    """Test finding overlapping features."""
    with BEDProcessor(overlapping_features_bed) as bed_proc:
        # Query position with multiple overlaps
        overlaps = bed_proc.find_overlaps('chr1', 2750)
        assert len(overlaps) == 3
        
        # Check returned objects
        for feature in overlaps:
            assert isinstance(feature, GenomicInterval)
            assert feature.chrom == 'chr1'
            assert feature.start <= 2750 <= feature.end
            assert 'name' in feature.metadata

def test_find_proximal(overlapping_features_bed):
    """Test finding proximal features within a window."""
    with BEDProcessor(overlapping_features_bed) as bed_proc:
        # Query with 500bp window
        proximal = bed_proc.find_proximal('chr1', 4500, 500)
        assert len(proximal) == 1
        assert proximal[0].metadata['name'] == 'feature2'  # Only feature2 extends close enough

def test_malformed_bed_handling(malformed_bed):
    """Test handling of malformed BED files."""
    with pytest.raises(ValueError):
        with BEDProcessor(malformed_bed) as bed_proc:
            # Attempting to build DB should raise error
            bed_proc._build_sqlite_db(rebuild=True)

def test_query_nonexistent_region(simple_bed):
    """Test querying a region with no features."""
    with BEDProcessor(simple_bed) as bed_proc:
        overlaps = bed_proc.find_overlaps('chr3', 1000)
        assert len(overlaps) == 0
        
        proximal = bed_proc.find_proximal('chr3', 1000, 1000)
        assert len(proximal) == 0
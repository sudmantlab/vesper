"""Tests for RepeatMasker processor."""

import os
import tempfile
import subprocess
import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
from dataclasses import dataclass, field

from vesper.processors.repeatmasker import RepeatMaskerProcessor
from vesper.models.variants import Variant, VariantAnalysis


class MockVariant:
    """Mock Variant class for testing."""
    def __init__(self, chrom, pos, id=None, ref="A", alt="AALU", svlen=None):
        self.chrom = chrom
        self.pos = pos
        self.ID = id or f"{chrom}_{pos}"
        self.ref = ref
        self.alt = alt
        self.info = {}
        if svlen:
            self.info["SVLEN"] = svlen


class MockVariantAnalysis:
    """Mock VariantAnalysis class for testing."""
    def __init__(self, variant):
        self.variant = variant
        self.repeatmasker_results = None

@pytest.fixture
def mock_variants():
    """Create a list of mock variants for testing."""
    return [
        MockVariantAnalysis(MockVariant("chr1", 1000, "var1", "A", "AAAGGCTGGCCAGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG")),
        MockVariantAnalysis(MockVariant("chr2", 2000, "var2", "G", "GTTTTTTAAGACAGGGTTTTGCTGTGTTGCCCAGGCTGGTCTTGAACT")),
        MockVariantAnalysis(MockVariant("chr3", 3000, "var3", "C", "<INS>"))  # Disregard <INS> placeholder sequences
    ]

@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield Path(temp_dir)

@pytest.fixture
def sample_out_content():
    """Sample RepeatMasker .out file content."""
    return """  SW   perc perc perc  query                    position in query           matching       repeat              position in repeat
score   div. del. ins.  sequence                     begin     end    (left)    repeat                   begin  end      (left)   ID

 1482   0.0  0.0  0.0  var1                          1    282  (0)    +  AluYa5                      (0)     1    282    (0)      1
 1409   8.3  4.0  0.4  var1                        291    587  (0)    C  AluY                       (31)   251      1    (0)      2
  389  20.9  0.0  2.3  var2                          1    133  (0)    +  L1HS                    (6031)     1    133  (6031)      3
  215  14.9  2.9  5.3  var2                        152    288  (0)    +  L1PA2                   (4321)   134    270  (2143)      4
"""

@pytest.fixture
def sample_out_with_overlap_content():
    """Sample RepeatMasker .out file content with secondary annotations."""
    return """  SW   perc perc perc  query                    position in query           matching       repeat              position in repeat
score   div. del. ins.  sequence                     begin     end    (left)    repeat                   begin  end      (left)   ID

 1482   0.0  0.0  0.0  var1                          1    282  (0)    +  AluYa5                      (0)     1    282    (0)      1
 1409   8.3  4.0  0.4  var1                        291    587  (0)    C  AluY                       (31)   251      1    (0)      2
  389  20.9  0.0  2.3  var2                          1    133  (0)    +  L1HS                    (6031)     1    133  (6031)      3
  215  14.9  2.9  5.3  var2                        152    288  (0)    +  L1PA2                   (4321)   134    270  (2143)      4 
  225  19.0 10.3  0.0  var2                        200    300  (0)    C  L1PA3                     (10)    100      1    (0)      5 *
"""

@pytest.fixture
def sample_out_file(temp_dir, sample_out_content):
    """Create a sample RepeatMasker .out file for testing."""
    outfile = temp_dir / "insertions.fa.out"
    with open(outfile, "w") as f:
        f.write(sample_out_content)
    return outfile

@pytest.fixture
def sample_out_with_overlap_file(temp_dir, sample_out_with_overlap_content):
    """Create a sample RepeatMasker .out file with secondary annotations."""
    outfile = temp_dir / "insertions.fa.out"
    with open(outfile, "w") as f:
        f.write(sample_out_with_overlap_content)
    return outfile

def test_initialization(temp_dir):
    """Test RepeatMaskerProcessor initialization."""
    processor = RepeatMaskerProcessor(temp_dir)
    assert processor is not None
    assert processor.repeatmasker_dir == temp_dir / "repeatmasker"
    assert processor.repeatmasker_dir.exists()

def test_create_temp_dir(temp_dir):
    """Test creating a temporary directory."""
    processor = RepeatMaskerProcessor(temp_dir)
    temp_dir = processor._create_temp_dir()
    
    assert temp_dir.exists()
    assert temp_dir.parent == processor.repeatmasker_dir
    assert "batch_" in temp_dir.name

def test_write_insertion_fasta(temp_dir, mock_variants):
    """Test writing insertion sequences to FASTA."""
    processor = RepeatMaskerProcessor(temp_dir)
    temp_dir = processor._create_temp_dir()
    
    fasta_path = processor._write_insertion_fasta(mock_variants, temp_dir)
    
    assert fasta_path.exists()
    
    with open(fasta_path, "r") as f:
        content = f.read()
    
    assert ">var1" in content
    assert ">var2" in content
    assert ">var3" not in content  # No insertion for var3
    assert "AAAGGCTGGCCAGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG" in content
    assert "GTTTTTTAAGACAGGGTTTTGCTGTGTTGCCCAGGCTGGTCTTGAACT" in content

@patch("subprocess.run")
def test_run_repeatmasker(mock_run, temp_dir, mock_variants):
    """Test running RepeatMasker."""
    mock_process = MagicMock()
    mock_process.stderr = "RepeatMasker completed"
    mock_process.stdout = "RepeatMasker output"
    mock_run.return_value = mock_process
    
    processor = RepeatMaskerProcessor(temp_dir)
    temp_dir = processor._create_temp_dir()
    
    fasta_path = temp_dir / "insertions.fa"
    with open(fasta_path, "w") as f:
        f.write(">var1\nAAAAAAAAAAAAAAAAAA\n")
    
    for suffix in [".fa.out", ".fa.tbl", ".fa.cat.gz", ".fa.masked", ".fa.out.gff"]:
        with open(fasta_path.with_suffix(suffix), "w") as f:
            f.write("mock content")
    
    processor._run_repeatmasker(fasta_path, temp_dir)
    mock_run.assert_called_once()
    cmd_args = mock_run.call_args[0][0]
    assert "RepeatMasker" in cmd_args
    assert "-species" in cmd_args
    assert "human" in cmd_args

@patch("subprocess.run")
def test_repeatmasker_failure(mock_run, temp_dir):
    """Test handling RepeatMasker failures."""
    mock_run.side_effect = subprocess.CalledProcessError(
        returncode=1,
        cmd=["RepeatMasker"],
        stderr="Error running RepeatMasker",
        output=""
    )
    
    processor = RepeatMaskerProcessor(temp_dir)
    temp_dir = processor._create_temp_dir()
    
    fasta_path = temp_dir / "insertions.fa"
    with open(fasta_path, "w") as f:
        f.write(">var1\nAAAAAAAAAAAAAAAAAA\n")
    
    with pytest.raises(RuntimeError):
        processor._run_repeatmasker(fasta_path, temp_dir)

def test_parse_output(temp_dir, sample_out_file):
    """Test parsing RepeatMasker output with a fixture."""
    processor = RepeatMaskerProcessor(temp_dir)
    
    # Parse the sample output file
    results = processor._parse_output(temp_dir)
    
    # Check basic structure
    assert isinstance(results, dict)
    assert len(results) == 2  # 2 different variants
    
    # Check first variant
    assert "var1" in results
    var1_results = results["var1"]
    assert len(var1_results) == 2  # 2 annotations
    
    # Check specific values for first annotation
    first_annotation = var1_results[0]
    assert first_annotation["score"] == 1482
    assert first_annotation["divergence"] == 0.0
    assert first_annotation["repeat_name"] == "AluYa5"
    
    # Check second variant
    assert "var2" in results
    var2_results = results["var2"]
    assert len(var2_results) == 2  # 2 annotations
    
    # Check for L1HS annotation
    l1_annotation = var2_results[0]
    assert l1_annotation["repeat_name"] == "L1HS"
    assert l1_annotation["score"] == 389

def test_parse_output_with_overlap(temp_dir, sample_out_with_overlap_file):
    """Test parsing RepeatMasker output with secondary annotations."""
    processor = RepeatMaskerProcessor(temp_dir)
    results = processor._parse_output(temp_dir)
    
    assert isinstance(results, dict)
    assert len(results) == 2
    
    assert "var2" in results
    var2_results = results["var2"]
    assert len(var2_results) == 3
    
    assert var2_results[0]["is_secondary"] == False  # Non-overlapping
    assert var2_results[1]["is_secondary"] == False  # Overlaps, primary relative to third annotation
    assert var2_results[2]["is_secondary"] == True  # Overlaps, secondary relative to second annotation

def test_parse_and_simplify(temp_dir, sample_out_file):
    """Test parsing and simplifying RepeatMasker results with a fixture."""
    processor = RepeatMaskerProcessor(temp_dir)
    
    # check that all annotations exist
    results = processor._parse_output(temp_dir)
    assert len(results["var1"]) == 2
    assert len(results["var2"]) == 2

    # check that the top-scoring annotation is returned for each variant
    top_results = processor._parse_and_simplify(temp_dir)
    
    # Check basic structure
    assert isinstance(top_results, dict)
    assert len(top_results) == 2  # 2 variants
    
    # Check that we got the best match for each variant
    assert top_results["var1"]["repeat_name"] == "AluYa5"  # Higher score
    assert top_results["var2"]["repeat_name"] == "L1HS"    # Higher score
    
    # Check that we have the right fields
    assert "score" in top_results["var1"].keys()
    assert "divergence" in top_results["var1"].keys()
    assert "match_coverage" in top_results["var1"].keys()


def test_parse_output_real_data(temp_dir):
    """Test parsing RepeatMasker output with real data."""
    processor = RepeatMaskerProcessor(temp_dir)
    
    # mock from provided real data
    real_output = Path("tests/files/hg38/894.filt.alt.fa.out")
    temp_output = temp_dir / "insertions.fa.out"
    temp_output.write_text(real_output.read_text())
    results = processor._parse_output(temp_dir)
    assert isinstance(results, dict)
    
    # Check specific variant with two annotations
    variant_id = "894_Sniffles2.INS.C5S11"
    assert variant_id in results.keys()
    variant_results = results[variant_id]
    assert len(variant_results) == 2 
    
    first_hit = variant_results[0]
    assert first_hit["score"] == 2255
    assert first_hit["divergence"] == 10.0
    assert first_hit["deletion"] == 1.0
    assert first_hit["insertion"] == 1.0
    assert first_hit["query_start"] == 72
    assert first_hit["query_end"] == 381
    assert first_hit["query_left"] == 389
    assert first_hit["strand"] == "+"
    assert first_hit["repeat_name"] == "AluSg"
    assert first_hit["repeat_class"] == "SINE/Alu"
    assert first_hit["repeat_start"] == 1
    assert first_hit["repeat_end"] == 310
    assert first_hit["repeat_left"] == 0
    
    second_hit = variant_results[1]
    assert second_hit["score"] == 180
    assert second_hit["repeat_name"] == "(TGAG)n"
    assert second_hit["repeat_class"] == "Simple_repeat"

def test_overlapping_annotations_in_real_data(temp_dir):
    """Test that secondary annotations are correctly identified in real data."""
    processor = RepeatMaskerProcessor(temp_dir)
    
    # mock from provided real data
    real_output = Path("tests/files/hg38/894.filt.alt.fa.out")
    temp_output = temp_dir / "insertions.fa.out"
    temp_output.write_text(real_output.read_text())
    results = processor._parse_output(temp_dir)
    
    variant_id = "894_Sniffles2.INS.1CES4" # 19 total annotations in real data
    assert variant_id in results
    variant_results = results[variant_id]
    
    # get all secondary annotations
    secondary_annotations = [hit for hit in variant_results if hit["is_secondary"]]
    assert len(secondary_annotations) > 0
    print(secondary_annotations)
    
    # get specific secondary hit (flagged with *)
    secondary_hit = next(hit for hit in secondary_annotations
                          if hit["is_secondary"] and hit["repeat_name"] == "L1MEg")
    
    assert secondary_hit["score"] == 225
    assert secondary_hit["divergence"] == 19.0
    assert secondary_hit["deletion"] == 10.3
    assert secondary_hit["insertion"] == 0.0
    assert secondary_hit["query_start"] == 3446
    assert secondary_hit["query_end"] == 3503

@patch("vesper.processors.repeatmasker.RepeatMaskerProcessor._run_repeatmasker")
@patch("vesper.processors.repeatmasker.RepeatMaskerProcessor._parse_output")
@patch("vesper.processors.repeatmasker.RepeatMaskerProcessor._parse_and_simplify")
def test_batch_analysis(mock_simplify, mock_parse, mock_run, temp_dir, mock_variants):
    """Test batch analysis of variants with RepeatMasker."""
    # Mock the parsing results
    mock_parse.return_value = {
        "var1": [
            {"score": 1482, "divergence": 0.0, "repeat_name": "AluYa5", "match_length": 282, "is_secondary": False}
        ],
        "var2": [
            {"score": 389, "divergence": 20.9, "repeat_name": "L1HS", "match_length": 133, "is_secondary": False}
        ]
    }
    
    mock_simplify.return_value = {
        "var1": {"best_match": "AluYa5", "score": 1482, "all_hits": 1, "coverage": 0.95},
        "var2": {"best_match": "L1HS", "score": 389, "all_hits": 1, "coverage": 0.85}
    }
    
    processor = RepeatMaskerProcessor(temp_dir)
    processor.batch_analysis(mock_variants)
    
    # Check that variants were annotated
    assert mock_variants[0].repeatmasker_results is not None
    assert mock_variants[1].repeatmasker_results is not None
    assert mock_variants[2].repeatmasker_results is None  # No insertion
    
    # Check the structure of the annotated results
    assert "detailed" in mock_variants[0].repeatmasker_results
    assert "summary" in mock_variants[0].repeatmasker_results
    
def test_context_manager(temp_dir):
    """Test using the processor as a context manager."""
    with RepeatMaskerProcessor(temp_dir) as processor:
        assert processor is not None
        temp_dir = processor._create_temp_dir()
        assert temp_dir.exists()
    
    # Temp dir should be cleaned up after context exit
    assert not temp_dir.exists()
    
def test_no_run_on_empty_variants(temp_dir):
    """Test that batch_analysis does nothing with empty variant list."""
    processor = RepeatMaskerProcessor(temp_dir)
    
    # Call with empty list
    processor.batch_analysis([])
    
    # No temp dir should be created
    assert processor.current_temp_dir is None
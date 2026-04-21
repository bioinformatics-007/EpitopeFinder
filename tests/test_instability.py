import pytest
import os
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock
from io import StringIO

# Add the parent directory to sys.path to find instability.py
sys.path.append(str(Path(__file__).resolve().parent.parent / "modules"))

import instability  # Import the module

# Sample FASTA content with amino acid sequences
SAMPLE_FASTA = """>seq1 description1
MKFLV
>seq2 description2
AILWY
"""

@pytest.fixture
def temp_fasta_file(tmp_path):
    """Create a temporary FASTA file for testing."""
    fasta_path = tmp_path / "test.fasta"
    with open(fasta_path, "w") as f:
        f.write(SAMPLE_FASTA)
    return fasta_path

@pytest.fixture
def empty_fasta_file(tmp_path):
    """Create an empty FASTA file for testing."""
    fasta_path = tmp_path / "empty.fasta"
    with open(fasta_path, "w") as f:
        f.write("")
    return fasta_path

def test_print_status(capsys):
    """Test the print_status function for different status types."""
    instability.print_status("Test info message", "info")
    captured = capsys.readouterr()
    assert "Test info message" in captured.out
    assert "\033[94m" in captured.out  # Blue color for info

    instability.print_status("Test success message", "success")
    captured = capsys.readouterr()
    assert "Test success message" in captured.out
    assert "\033[92m" in captured.out  # Green color for success

    instability.print_status("Test error message", "error")
    captured = capsys.readouterr()
    assert "Test error message" in captured.out
    assert "\033[91m" in captured.out  # Red color for error

@patch("instability.ProteinAnalysis")
def test_protparm(mock_protein_analysis):
    """Test the protparm function with amino acid sequences."""
    # Mock ProteinAnalysis
    mock_instance = MagicMock()
    mock_instance.instability_index.return_value = 30.1234
    mock_protein_analysis.return_value = mock_instance

    # Test stable sequence
    ii, stab, stab_coff = instability.protparm("MKFLV")
    assert ii == 30.12  # Rounded to 2 decimal places
    assert stab == "stable"
    assert stab_coff == 1

    # Test unstable sequence
    mock_instance.instability_index.return_value = 50.5678
    ii, stab, stab_coff = instability.protparm("AILWY")
    assert ii == 50.57
    assert stab == "unstable"
    assert stab_coff == 0

@patch("instability.ProteinAnalysis")
def test_protparm_error(mock_protein_analysis):
    """Test protparm with an invalid amino acid sequence."""
    # No need to mock ValueError from ProteinAnalysis since is_valid_sequence raises it
    with pytest.raises(ValueError, match="Sequence contains invalid amino acid characters"):
        instability.protparm("XYZ")  # Invalid amino acids

@patch("instability.SeqIO.parse")
@patch("instability.protparm")
def test_run_instability_success(mock_protparm, mock_seqio_parse, temp_fasta_file, tmp_path):
    """Test run_instability with a valid FASTA file."""
    # Mock SeqIO.parse
    mock_record1 = MagicMock(id="seq1", description="description1", seq="MKFLV")
    mock_record2 = MagicMock(id="seq2", description="description2", seq="AILWY")
    mock_seqio_parse.return_value = [mock_record1, mock_record2]

    # Mock protparm
    mock_protparm.side_effect = [
        (30.12, "stable", 1),  # seq1
        (50.57, "unstable", 0)  # seq2
    ]

    output_dir = tmp_path / "output"
    status = instability.run_instability(
        temp_fasta_file,
        output_dir=str(output_dir),
        binary_output_file="binary.txt",  # Ignored by function
        value_output_file="value.txt"     # Ignored by function
    )

    assert status == 0
    # Check output files
    seq1_csv = output_dir / "seq1_instability.csv"
    seq2_csv = output_dir / "seq2_instability.csv"
    combined_csv = output_dir / "combined_instability_batch_1.csv"
    assert seq1_csv.exists()
    assert seq2_csv.exists()
    assert combined_csv.exists()

    # Verify contents of seq1_instability.csv
    with open(seq1_csv) as f:
        seq1_content = f.read()
        assert "seq1,description1,MKFLV,30.12,stable,1" in seq1_content

    # Verify contents of seq2_instability.csv
    with open(seq2_csv) as f:
        seq2_content = f.read()
        assert "seq2,description2,AILWY,50.57,unstable,0" in seq2_content

    # Verify contents of combined_instability_batch_1.csv
    with open(combined_csv) as f:
        combined_content = f.read()
        assert "seq1,description1,MKFLV,30.12,stable,1" in combined_content
        assert "seq2,description2,AILWY,50.57,unstable,0" in combined_content

@patch("instability.SeqIO.parse")
def test_run_instability_empty_fasta(mock_seqio_parse, empty_fasta_file, tmp_path, capsys):
    """Test run_instability with an empty FASTA file."""
    mock_seqio_parse.return_value = []
    output_dir = tmp_path / "output"
    status = instability.run_instability(
        empty_fasta_file,
        output_dir=str(output_dir),
        binary_output_file="binary.txt",
        value_output_file="value.txt"
    )
    assert status == 1
    captured = capsys.readouterr()
    assert "No valid sequences found in FASTA file" in captured.out

def test_run_instability_nonexistent_fasta(tmp_path, capsys):
    """Test run_instability with a nonexistent FASTA file."""
    nonexistent_file = tmp_path / "nonexistent.fasta"
    output_dir = tmp_path / "output"
    status = instability.run_instability(
        nonexistent_file,
        output_dir=str(output_dir),
        binary_output_file="binary.txt",
        value_output_file="value.txt"
    )
    assert status == 1
    captured = capsys.readouterr()
    assert "General error during instability analysis: Input FASTA file not found" in captured.out

@patch("instability.SeqIO.parse")
def test_run_instability_file_error(mock_seqio_parse, temp_fasta_file, tmp_path):
    """Test run_instability when file writing fails."""
    mock_seqio_parse.side_effect = IOError("File error")
    output_dir = tmp_path / "output"
    status = instability.run_instability(
        temp_fasta_file,
        output_dir=str(output_dir),
        binary_output_file="binary.txt",
        value_output_file="value.txt"
    )
    assert status == 1

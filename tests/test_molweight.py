import pytest
import pandas as pd
from pathlib import Path
from unittest.mock import patch
from io import StringIO
import logging
import tempfile
import shutil
import sys

# Add parent directory to sys.path to access modules/
sys.path.append(str(Path(__file__).resolve().parent.parent))

from modules.molweight import (
    print_status,
    is_valid_sequence,
    protparm,
    validate_paths,
    verify_output_files,
    run_molwt
)

# Set up logging capture
@pytest.fixture
def capture_logs():
    """Fixture to capture log output with consistent formatting."""
    logger = logging.getLogger('MolWt')
    logger.handlers = []  # Clear existing handlers
    log_capture = StringIO()
    handler = logging.StreamHandler(log_capture)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    yield log_capture
    handler.close()

# Fixture for temporary directory
@pytest.fixture
def temp_dir():
    """Create a temporary directory for tests, cleaned up afterward."""
    temp_dir = tempfile.mkdtemp()
    yield Path(temp_dir)
    shutil.rmtree(temp_dir, ignore_errors=True)

# Fixture for sample FASTA content
@pytest.fixture
def sample_fasta_content():
    """Provide valid FASTA content for testing."""
    return """>sp|P12345|TEST1 Test protein 1
ACDEFGHIKLMNPQRSTVWY
>sp|Q67890|TEST2 Test protein 2
LMNPQRSTVWYACDEFGHIK
"""

# Test print_status
def test_print_status(capture_logs, capsys):
    """Test print_status for colored console output and logging."""
    msg = "Test message"
    print_status(msg, "info")
    captured = capsys.readouterr()
    assert "\033[94mTest message\033[0m" in captured.out
    log_output = capture_logs.getvalue()
    assert "- MolWt - Test message" in log_output

# Test is_valid_sequence
def test_is_valid_sequence_valid():
    """Test is_valid_sequence with valid amino acid sequence."""
    assert is_valid_sequence("ACDEFGHIKLMNPQRSTVWY") is True

def test_is_valid_sequence_invalid():
    """Test is_valid_sequence with invalid characters."""
    assert is_valid_sequence("ACDEXYZ") is False

def test_is_valid_sequence_empty():
    """Test is_valid_sequence with empty sequence."""
    assert is_valid_sequence("") is True  # Empty sequence is technically valid per function

# Test protparm
def test_protparm_valid():
    """Test protparm with valid sequence."""
    seq = "ACDEFGHIKLMNPQRSTVWY"
    mol_weight = protparm(seq)
    assert isinstance(mol_weight, float)
    assert mol_weight > 0

def test_protparm_invalid(capture_logs):
    """Test protparm with invalid sequence."""
    seq = "ACDEXYZ"
    with pytest.raises(ValueError) as exc:
        protparm(seq)
    assert "Sequence contains invalid amino acid characters" in str(exc.value)
    log_output = capture_logs.getvalue()
    assert "Error computing molecular weight" in log_output

# Test validate_paths
def test_validate_paths_valid(sample_fasta_content, temp_dir):
    """Test validate_paths with valid input FASTA and output directory."""
    fasta_file = temp_dir / "test.fasta"
    with open(fasta_file, "w", newline='') as f:
        f.write(sample_fasta_content)
    output_dir = temp_dir / "output"
    input_fasta, out_dir, out_file = validate_paths(fasta_file, output_dir)
    assert input_fasta == fasta_file
    assert out_dir == output_dir
    assert out_file is None

def test_validate_paths_missing_fasta(temp_dir, capture_logs):
    """Test validate_paths with missing FASTA file."""
    fasta_file = temp_dir / "nonexistent.fasta"
    output_dir = temp_dir / "output"
    with pytest.raises(FileNotFoundError) as exc:
        validate_paths(fasta_file, output_dir)
    assert "Input FASTA file not found" in str(exc.value)

def test_validate_paths_output_file(sample_fasta_content, temp_dir):
    """Test validate_paths with specified output file."""
    fasta_file = temp_dir / "test.fasta"
    with open(fasta_file, "w", newline='') as f:
        f.write(sample_fasta_content)
    output_dir = temp_dir / "output"
    output_file = temp_dir / "output/result.csv"
    input_fasta, out_dir, out_file = validate_paths(fasta_file, output_dir, output_file)
    assert input_fasta == fasta_file
    assert out_dir == output_dir
    assert out_file == output_file

# Test verify_output_files
def test_verify_output_files(temp_dir, capsys):
    """Test verify_output_files with various file states."""
    file1 = temp_dir / "file1.csv"
    file2 = temp_dir / "file2.csv"
    file3 = temp_dir / "file3.csv"
    with open(file1, "w") as f:
        f.write("data")
    file2.touch()
    verify_output_files([file1, file2, file3])
    captured = capsys.readouterr()
    assert "Output file found and non-empty: " in captured.out
    assert "Output file is empty: " in captured.out
    assert "Output file not found: " in captured.out

# Test run_molwt
def test_run_molwt_success(sample_fasta_content, temp_dir, capture_logs):
    """Test run_molwt with valid FASTA file."""
    fasta_file = temp_dir / "test.fasta"
    with open(fasta_file, "w", newline='') as f:
        f.write(sample_fasta_content)
    output_dir = temp_dir / "output"
    result = run_molwt(fasta_file, output_dir, batch_idx=1)
    assert result == 0
    combined_file = output_dir / "combined_molwt_batch_1.csv"
    assert combined_file.exists()
    df = pd.read_csv(combined_file)
    assert len(df) == 2
    assert set(df["Sequence_ID"]) == {"P12345", "Q67890"}
    log_output = capture_logs.getvalue()
    assert "Validated paths" in log_output
    assert "Molecular weight processing completed" in log_output

def test_run_molwt_invalid_sequence(temp_dir, capture_logs, capsys):
    """Test run_molwt with FASTA containing invalid sequence."""
    fasta_file = temp_dir / "test.fasta"
    invalid_fasta = """>sp|P12345|TEST1 Test protein 1
ACDEXYZ
"""
    with open(fasta_file, "w", newline='') as f:
        f.write(invalid_fasta)
    output_dir = temp_dir / "output"
    result = run_molwt(fasta_file, output_dir, batch_idx=1)
    assert result == 1
    combined_file = output_dir / "combined_molwt_batch_1.csv"
    assert not combined_file.exists()
    captured = capsys.readouterr()
    assert "Invalid amino acid characters in P12345" in captured.out
    log_output = capture_logs.getvalue()
    assert "Invalid amino acid characters in P12345" in log_output
    assert "No valid sequences found in FASTA file" in log_output

def test_run_molwt_empty_fasta(temp_dir, capture_logs, capsys):
    """Test run_molwt with empty FASTA file."""
    fasta_file = temp_dir / "empty.fasta"
    fasta_file.touch()
    output_dir = temp_dir / "output"
    result = run_molwt(fasta_file, output_dir, batch_idx=1)
    assert result == 1
    combined_file = output_dir / "combined_molwt_batch_1.csv"
    assert not combined_file.exists()
    captured = capsys.readouterr()
    assert "No valid sequences found in FASTA file" in captured.out
    log_output = capture_logs.getvalue()
    assert "No valid sequences found in FASTA file" in log_output

# Test main
@patch("sys.argv", ["molweight.py", "-i", "test.fasta", "-d", "output"])
def test_main_valid_args(temp_dir, sample_fasta_content, capture_logs, capsys):
    """Test main function with valid arguments."""
    input_file = temp_dir / "test.fasta"
    output_dir = temp_dir / "output"
    with open(input_file, "w", newline='') as f:
        f.write(sample_fasta_content)
    # Patch sys.argv to point to the temporary input file
    patched_argv = ["molweight.py", "-i", str(input_file), "-d", str(output_dir)]
    with patch("sys.argv", patched_argv):
        from modules.molweight import main
        try:
            main()
        except SystemExit as exc:
            assert exc.code == 0
    log_output = capture_logs.getvalue()
    assert "Validated paths" in log_output
    # Verify output files
    combined_file = output_dir / "combined_molwt_batch_1.csv"
    assert combined_file.exists()
    df = pd.read_csv(combined_file)
    assert len(df) == 2
    assert set(df["Sequence_ID"]) == {"P12345", "Q67890"}

@patch("sys.argv", ["molweight.py", "-i", "nonexistent.fasta", "-d", "output"])
def test_main_invalid_input(temp_dir, capture_logs, capsys):
    """Test main function with invalid input file."""
    output_dir = temp_dir / "output"
    try:
        from modules.molweight import main
        main()
    except SystemExit as exc:
        assert exc.code == 1
    captured = capsys.readouterr()
    assert "Input FASTA file not found" in captured.out
    log_output = capture_logs.getvalue()
    assert "Input FASTA file not found" in log_output

import pytest
import os
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock
import logging
import tempfile
import subprocess
from Bio import SeqIO
import sys
import shutil

# Adjust sys.path to include the directory containing algpred.py
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'modules')))
from algpred import print_status, is_csv_empty, get_fasta_sequences, validate_fasta_file, verify_output_files, run_algpred

# Configure logging for tests
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("TestAlgpred")

@pytest.fixture
def temp_fasta_file():
    """Create a temporary FASTA file with UniProt-style IDs."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">sp|P12345|TEST1\nACDEFGHIK\n>sp|Q67890|TEST2\nLMNPQRSTV\n")
    yield f.name
    os.unlink(f.name)

@pytest.fixture
def temp_dir():
    """Create a temporary directory."""
    temp_dir = tempfile.mkdtemp()
    yield Path(temp_dir)
    shutil.rmtree(temp_dir, ignore_errors=True)

@pytest.fixture
def caplog(caplog):
    """Fixture to capture log messages."""
    caplog.set_level(logging.DEBUG)
    return caplog

def test_print_status(capsys, caplog):
    """Test print_status function with different status types."""
    print_status("Test info message", "info")
    captured = capsys.readouterr()
    assert "Test info message" in captured.out
    assert "INFO" in caplog.text
    assert "Test info message" in caplog.text

    print_status("Test success message", "success")
    captured = capsys.readouterr()
    assert "Test success message" in captured.out
    assert "INFO" in caplog.text
    assert "Test success message" in caplog.text

    print_status("Test warning message", "warning")
    captured = capsys.readouterr()
    assert "Test warning message" in captured.out
    assert "WARNING" in caplog.text
    assert "Test warning message" in caplog.text

    print_status("Test error message", "error")
    captured = capsys.readouterr()
    assert "Test error message" in captured.out
    assert "ERROR" in caplog.text
    assert "Test error message" in caplog.text

def test_is_csv_empty(tmp_path):
    """Test is_csv_empty function."""
    non_existent = tmp_path / "non_existent.csv"
    assert is_csv_empty(non_existent)

    empty_csv = tmp_path / "empty.csv"
    empty_csv.write_text("")
    assert is_csv_empty(empty_csv)

    headers_csv = tmp_path / "headers.csv"
    pd.DataFrame(columns=["A", "B"]).to_csv(headers_csv, index=False)
    assert is_csv_empty(headers_csv)

    data_csv = tmp_path / "data.csv"
    pd.DataFrame({"A": [1, 2], "B": [3, 4]}).to_csv(data_csv, index=False)
    assert not is_csv_empty(data_csv)

def test_get_fasta_sequences(temp_fasta_file):
    """Test get_fasta_sequences function."""
    sequences = get_fasta_sequences(temp_fasta_file)
    assert len(sequences) == 2
    assert sequences[0] == ("P12345", "ACDEFGHIK")
    assert sequences[1] == ("Q67890", "LMNPQRSTV")

    empty_fasta = temp_fasta_file + ".empty"
    with open(empty_fasta, 'w') as f:
        pass
    assert get_fasta_sequences(empty_fasta) == []

    invalid_fasta = temp_fasta_file + ".invalid"
    with open(invalid_fasta, 'w') as f:
        f.write("not a fasta file")
    assert get_fasta_sequences(invalid_fasta) == []

def test_validate_fasta_file(temp_fasta_file, tmp_path, capsys, caplog):
    """Test validate_fasta_file function."""
    # Test valid FASTA file
    is_valid, clean_fasta, temp_dir, error_msg = validate_fasta_file(temp_fasta_file)
    assert is_valid
    assert clean_fasta is not None
    assert os.path.isfile(clean_fasta)
    assert temp_dir is not None
    assert error_msg is None
    captured = capsys.readouterr()
    assert "FASTA file validated and cleaned" in captured.out
    assert "Temporary directory created" in caplog.text

    with open(clean_fasta, 'r') as f:
        content = f.read()
    assert ">P12345\nACDEFGHIK\n>Q67890\nLMNPQRSTV" in content

    # Test non-existent FASTA file
    non_existent = tmp_path / "non_existent.fasta"
    is_valid, clean_fasta, temp_dir, error_msg = validate_fasta_file(non_existent)
    assert not is_valid
    assert clean_fasta is None
    assert temp_dir is None
    assert "Input FASTA file not found" in error_msg

    # Test empty FASTA file
    empty_fasta = tmp_path / "empty.fasta"
    empty_fasta.write_text("")
    is_valid, clean_fasta, temp_dir, error_msg = validate_fasta_file(empty_fasta)
    assert not is_valid
    assert clean_fasta is None
    assert temp_dir is None
    assert "FASTA file is empty" in error_msg

    # Test invalid sequence (short sequence)
    short_fasta = tmp_path / "short.fasta"
    with open(short_fasta, 'w') as f:
        f.write(">sp|P12345|TEST1\nAC\n")
    is_valid, clean_fasta, temp_dir, error_msg = validate_fasta_file(short_fasta)
    assert not is_valid
    assert clean_fasta is None
    assert temp_dir is None
    assert "No valid sequences found" in error_msg
    assert "Sequence too short" in caplog.text

    # Test invalid characters
    invalid_fasta = tmp_path / "invalid.fasta"
    with open(invalid_fasta, 'w') as f:
        f.write(">sp|P12345|TEST1\nACDXFG\n")
    is_valid, clean_fasta, temp_dir, error_msg = validate_fasta_file(invalid_fasta)
    assert not is_valid
    assert clean_fasta is None
    assert temp_dir is None
    assert "No valid sequences found" in error_msg
    assert "Invalid characters" in caplog.text

def test_verify_output_files(tmp_path, capsys):
    """Test verify_output_files function."""
    file1 = tmp_path / "file1.csv"
    file1.write_text("data")
    file2 = tmp_path / "file2.csv"
    file2.write_text("")
    file3 = tmp_path / "file3.csv"

    verify_output_files([file1, file2, file3])
    captured = capsys.readouterr()
    assert "Output file found and non-empty: {}".format(file1) in captured.out
    assert "Output file is empty: {}".format(file2) in captured.out
    assert "Output file not found: {}".format(file3) in captured.out

@patch("algpred.verify_output_files")
@patch("algpred.get_fasta_sequences")
@patch("algpred.validate_fasta_file")
@patch("subprocess.run")
@patch("uuid.uuid4")
def test_run_algpred_success(mock_uuid, mock_subprocess_run, mock_validate_fasta_file, mock_get_fasta_sequences, mock_verify_output_files, temp_fasta_file, temp_dir, capsys, caplog):
    """Test run_algpred with successful AlgPred2 execution."""
    logger.debug("Starting test_run_algpred_success")

    # Mock validate_fasta_file
    mock_validate_fasta_file.return_value = (True, temp_fasta_file, temp_dir, None)
    # Mock get_fasta_sequences
    mock_get_fasta_sequences.return_value = [("P12345", "ACDEFGHIK"), ("Q67890", "LMNPQRSTV")]
    # Mock verify_output_files
    mock_verify_output_files.return_value = None

    # Mock UUID
    uuid_mock = MagicMock()
    uuid_mock.hex = "fixed-uuid-for-test"
    mock_uuid.return_value = uuid_mock
    temp_output_uuid = temp_dir / f"algpred_temp_fixed-uuid-for-test.csv"

    # Mock subprocess.run
    mock_subprocess_run.return_value = MagicMock(
        returncode=0,
        stdout="AlgPred2 completed successfully",
        stderr=""
    )

    # Create UUID output file
    pd.DataFrame({
        "Sequence": ["ACDEFGHIK"],
        "Score": [0.8],
        "Prediction": ["Allergen"]
    }).to_csv(temp_output_uuid, index=False)

    # Mock os.path.isfile and os.path.exists for input and output files
    def isfile_side_effect(path):
        return str(path) == str(temp_fasta_file) or str(path) == str(temp_output_uuid)

    def exists_side_effect(path):
        return str(path) == str(temp_fasta_file) or str(path) == str(temp_dir) or str(path) == str(temp_output_uuid) or str(path).endswith("P12345_algpred.csv") or str(path).endswith("Q67890_algpred.csv") or str(path).endswith("combined_algpred_batch_1.csv")

    with patch("os.path.isfile", side_effect=isfile_side_effect), \
         patch("os.path.exists", side_effect=exists_side_effect), \
         patch("algpred.is_csv_empty", return_value=False):
        exit_code = run_algpred(
            input_fasta=temp_fasta_file,
            output_dir=temp_dir,
            batch_idx=1,
            model="1"
        )

    logger.debug(f"Exit code: {exit_code}")
    captured = capsys.readouterr()
    logger.debug(f"Captured stdout: {captured.out}")

    assert exit_code == 0
    assert "AlgPred2 completed (Model 1, Threshold 0.3)" in captured.out
    assert "Combined output saved to" in captured.out
    assert "AlgPred2 processing completed" in caplog.text

@patch("algpred.validate_fasta_file")
@patch("subprocess.run")
def test_run_algpred_no_output(mock_subprocess_run, mock_validate_fasta_file, temp_fasta_file, temp_dir, capsys, caplog):
    """Test run_algpred when AlgPred2 produces no output."""
    mock_subprocess_run.return_value = MagicMock(
        returncode=0,
        stdout="AlgPred2 completed successfully",
        stderr=""
    )

    # Mock validate_fasta_file
    mock_validate_fasta_file.return_value = (True, temp_fasta_file, temp_dir, None)

    # Mock os.path.isfile and os.path.exists for input file
    def isfile_side_effect(path):
        return str(path) == temp_fasta_file

    def exists_side_effect(path):
        return str(path) == temp_fasta_file or str(path) == str(temp_dir) or str(path).endswith("P12345_algpred.csv") or str(path).endswith("Q67890_algpred.csv") or str(path).endswith("combined_algpred_batch_1.csv")

    with patch("os.path.isfile", side_effect=isfile_side_effect), \
         patch("os.path.exists", side_effect=exists_side_effect), \
         patch("algpred.is_csv_empty", return_value=True):
        exit_code = run_algpred(
            input_fasta=temp_fasta_file,
            output_dir=temp_dir,
            batch_idx=1,
            model="1"
        )

    assert exit_code == 0
    captured = capsys.readouterr()
    assert "No allergens detected for" in captured.out
    assert "Combined output saved to" in captured.out
    assert "No allergens detected" in caplog.text

@patch("algpred.verify_output_files")
@patch("algpred.get_fasta_sequences")
@patch("algpred.validate_fasta_file")
@patch("subprocess.run")
@patch("uuid.uuid4")
def test_run_algpred_failure(mock_uuid, mock_subprocess_run, mock_validate_fasta_file, mock_get_fasta_sequences, mock_verify_output_files, temp_fasta_file, temp_dir, capsys, caplog):
    """Test run_algpred with AlgPred2 failure and model fallback."""
    logger.debug("Starting test_run_algpred_failure")

    # Mock validate_fasta_file
    mock_validate_fasta_file.return_value = (True, temp_fasta_file, temp_dir, None)
    # Mock get_fasta_sequences
    mock_get_fasta_sequences.return_value = [("P12345", "ACDEFGHIK"), ("Q67890", "LMNPQRSTV")]
    # Mock verify_output_files
    mock_verify_output_files.return_value = None

    # Mock UUID
    uuid_mock = MagicMock()
    uuid_mock.hex = "fixed-uuid-for-test"
    mock_uuid.return_value = uuid_mock
    temp_output_uuid = temp_dir / f"algpred_temp_fixed-uuid-for-test.csv"

    # Mock subprocess.run with fallback logic
    def side_effect(*args, **kwargs):
        if " -m 2 " in args[0]:
            pd.DataFrame({
                "Sequence": ["ACDEFGHIK"],
                "Score": [0.8],
                "Prediction": ["Allergen"]
            }).to_csv(temp_output_uuid, index=False)
            return MagicMock(returncode=0, stdout="AlgPred2 completed successfully", stderr="")
        raise subprocess.CalledProcessError(1, "algpred2", stderr="Failed")

    mock_subprocess_run.side_effect = side_effect

    # Mock os.path.isfile and os.path.exists for input and output files
    def isfile_side_effect(path):
        return str(path) == str(temp_fasta_file) or str(path) == str(temp_output_uuid)

    def exists_side_effect(path):
        return str(path) == str(temp_fasta_file) or str(path) == str(temp_dir) or str(path) == str(temp_output_uuid) or str(path).endswith("P12345_algpred.csv") or str(path).endswith("Q67890_algpred.csv") or str(path).endswith("combined_algpred_batch_1.csv")

    with patch("os.path.isfile", side_effect=isfile_side_effect), \
         patch("os.path.exists", side_effect=exists_side_effect), \
         patch("algpred.is_csv_empty", return_value=False):
        exit_code = run_algpred(
            input_fasta=temp_fasta_file,
            output_dir=temp_dir,
            batch_idx=1,
            model="1"
        )

    logger.debug(f"Exit code: {exit_code}")
    captured = capsys.readouterr()
    logger.debug(f"Captured stdout: {captured.out}")

    assert exit_code == 0
    assert "All thresholds failed with Model 1. Retrying with Model 2" in captured.out
    assert "AlgPred2 completed (Model 2, Threshold 0.3)" in captured.out
    assert "Combined output saved to" in captured.out
    assert "Retrying with Model 2" in caplog.text

@patch("subprocess.run")
def test_run_algpred_invalid_fasta(mock_subprocess_run, temp_dir, capsys, caplog):
    """Test run_algpred with invalid FASTA file."""
    invalid_fasta = temp_dir / "invalid.fasta"
    invalid_fasta.write_text("not a fasta file")

    exit_code = run_algpred(
        input_fasta=invalid_fasta,
        output_dir=temp_dir,
        batch_idx=1,
        model="1"
    )

    assert exit_code == 1
    captured = capsys.readouterr()
    assert "Invalid FASTA file" in captured.out
    assert "FASTA file does not start with '>'" in caplog.text
    mock_subprocess_run.assert_not_called()

if __name__ == "__main__":
    pytest.main([__file__])

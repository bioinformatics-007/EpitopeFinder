import pytest
import os
import pandas as pd
from pathlib import Path
from unittest.mock import patch, mock_open, MagicMock
from io import StringIO
import logging
import tempfile
import shutil
import sys
import subprocess

# Add parent directory to sys.path to access modules/
sys.path.append(str(Path(__file__).resolve().parent.parent))

from modules.mhc_ii import find_python_interpreter, print_status, log_memory_usage, validate_fasta, get_valid_alleles, suggest_alternative_methods, run_mhc2_prediction, parse_results, run_mhc2

# Set up logging capture
@pytest.fixture
def capture_logs():
    """Fixture to capture log output with consistent formatting."""
    logger = logging.getLogger('VaxElan')
    logger.handlers = []  # Clear existing handlers
    log_capture = StringIO()
    handler = logging.StreamHandler(log_capture)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
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

# Test find_python_interpreter
@patch("shutil.which")
@patch("subprocess.run")
def test_find_python_interpreter_success(mock_subprocess, mock_which, capture_logs):
    """Test find_python_interpreter with a valid Python 3 interpreter."""
    mock_which.return_value = "/usr/bin/python3"
    mock_subprocess.return_value = MagicMock(stdout="Python 3.8.0\n", stderr="")
    result = find_python_interpreter()
    assert result == "/usr/bin/python3"
    log_output = capture_logs.getvalue()
    assert "Using Python 3 interpreter: /usr/bin/python3" in log_output

@patch("shutil.which")
def test_find_python_interpreter_none(mock_which, capture_logs):
    """Test find_python_interpreter when no interpreter is found."""
    mock_which.return_value = None
    result = find_python_interpreter()
    assert result is None
    log_output = capture_logs.getvalue()
    assert "No Python 3 interpreter found" in log_output

# Test print_status
def test_print_status(capture_logs, capsys):
    """Test print_status for colored console output and logging."""
    msg = "Test message"
    print_status(msg, "info")
    captured = capsys.readouterr()
    assert "\033[94mTest message\033[0m" in captured.out
    log_output = capture_logs.getvalue()
    assert "INFO - Test message" in log_output

# Test log_memory_usage
@patch("psutil.Process")
def test_log_memory_usage(mock_process):
    """Test log_memory_usage returns memory in MB."""
    mock_process.return_value.memory_info.return_value.rss = 104857600  # 100 MB
    mem_usage = log_memory_usage()
    assert mem_usage == 100.0

# Test validate_fasta
def test_validate_fasta_valid(sample_fasta_content, temp_dir, capture_logs):
    """Test validate_fasta with a valid FASTA file."""
    fasta_file = temp_dir / "test.fasta"
    with open(fasta_file, "w", newline='') as f:
        f.write(sample_fasta_content)
    assert validate_fasta(fasta_file) is True
    log_output = capture_logs.getvalue()
    assert log_output == ""  # No errors logged

def test_validate_fasta_empty(temp_dir, capture_logs):
    """Test validate_fasta with an empty FASTA file."""
    fasta_file = temp_dir / "empty.fasta"
    fasta_file.touch()
    assert validate_fasta(fasta_file) is False
    log_output = capture_logs.getvalue()
    assert "Input FASTA" in log_output
    assert "is empty" in log_output

def test_validate_fasta_invalid_chars(sample_fasta_content, temp_dir, capture_logs):
    """Test validate_fasta with invalid characters."""
    fasta_file = temp_dir / "test.fasta"
    invalid_fasta = sample_fasta_content.replace("ACDEFGHIKLMNPQRSTVWY", "ACDEXYZ")
    with open(fasta_file, "w", newline='') as f:
        f.write(invalid_fasta)
    assert validate_fasta(fasta_file) is False
    log_output = capture_logs.getvalue()
    assert "Invalid amino acid characters" in log_output

def test_validate_fasta_dna(sample_fasta_content, temp_dir, capture_logs):
    """Test validate_fasta with DNA-like sequences."""
    fasta_file = temp_dir / "test.fasta"
    dna_fasta = sample_fasta_content.replace("ACDEFGHIKLMNPQRSTVWY", "ACGTACGT")
    with open(fasta_file, "w", newline='') as f:
        f.write(dna_fasta)
    assert validate_fasta(fasta_file) is False
    log_output = capture_logs.getvalue()
    assert "contains DNA-like sequences" in log_output

# Test get_valid_alleles
@patch("subprocess.run")
@patch("modules.mhc_ii.find_python_interpreter")
def test_get_valid_alleles(mock_find_python, mock_subprocess, capture_logs):
    """Test get_valid_alleles with valid output."""
    mock_find_python.return_value = "/usr/bin/python3"
    mock_subprocess.return_value = MagicMock(
        stdout="HLA-DRB1*01:01\nHLA-DRB1*03:01\n----\n",
        stderr=""
    )
    alleles = get_valid_alleles("consensus3")
    assert alleles == ["HLA-DRB1*01:01", "HLA-DRB1*03:01"]
    log_output = capture_logs.getvalue()
    assert "Raw allele output for method consensus3" in log_output

@patch("subprocess.run")
@patch("modules.mhc_ii.find_python_interpreter")
def test_get_valid_alleles_empty(mock_find_python, mock_subprocess, capture_logs):
    """Test get_valid_alleles with empty output."""
    mock_find_python.return_value = "/usr/bin/python3"
    mock_subprocess.return_value = MagicMock(stdout="----\n", stderr="")
    alleles = get_valid_alleles("consensus3")
    assert alleles == [
        'HLA-DRB1*01:01', 'HLA-DRB1*01:02', 'HLA-DRB1*03:01', 'HLA-DRB1*04:01',
        'HLA-DRB1*04:04', 'HLA-DRB1*04:05', 'HLA-DRB1*07:01', 'HLA-DRB1*08:01',
        'HLA-DRB1*08:02', 'HLA-DRB1*09:01', 'HLA-DRB1*11:01', 'HLA-DRB1*11:04',
        'HLA-DRB1*12:01', 'HLA-DRB1*13:01', 'HLA-DRB1*13:02', 'HLA-DRB1*15:01',
        'HLA-DRB1*15:02', 'HLA-DRB3*01:01', 'HLA-DRB3*02:01', 'HLA-DRB4*01:01',
        'HLA-DRB5*01:01', 'HLA-DQA1*05:01/DQB1*02:01', 'HLA-DQA1*05:01/DQB1*03:01',
        'HLA-DQA1*03:01/DQB1*03:02', 'HLA-DQA1*04:01/DQB1*04:02',
        'HLA-DQA1*01:01/DQB1*05:01', 'HLA-DQA1*01:02/DQB1*06:02',
        'HLA-DPA1*02:01/DPB1*01:01', 'HLA-DPA1*01:03/DPB1*02:01',
        'HLA-DPA1*03:01/DPB1*04:02', 'HLA-DPA1*02:01/DPB1*05:01'
    ]
    log_output = capture_logs.getvalue()
    assert "No alleles parsed for method consensus3" in log_output

# Test suggest_alternative_methods
@patch("subprocess.run")
@patch("modules.mhc_ii.find_python_interpreter")
def test_suggest_alternative_methods(mock_find_python, mock_subprocess, capture_logs):
    """Test suggest_alternative_methods with valid methods."""
    mock_find_python.return_value = "/usr/bin/python3"
    mock_subprocess.side_effect = [
        MagicMock(stdout="consensus3\nnetmhciipan_el\n", stderr=""),  # methods
        MagicMock(stdout="HLA-DRB1*01:01\nHLA-DRB1*03:01\n", stderr=""),  # alleles for consensus3
        MagicMock(stdout="HLA-DQA1*05:01/DQB1*02:01\n", stderr="")  # alleles for netmhciipan_el
    ]
    suggestion = suggest_alternative_methods()
    assert "Try alternative methods:" in suggestion
    assert "- consensus3: Supports alleles like HLA-DRB1*01:01, HLA-DRB1*03:01" in suggestion
    assert "- netmhciipan_el: Supports alleles like HLA-DQA1*05:01/DQB1*02:01" in suggestion

# Test run_mhc2_prediction
@patch("subprocess.run")
@patch("modules.mhc_ii.find_python_interpreter")
def test_run_mhc2_prediction_success(mock_find_python, mock_subprocess, sample_fasta_content, temp_dir, capture_logs):
    """Test run_mhc2_prediction with successful execution."""
    fasta_file = temp_dir / "test.fasta"
    with open(fasta_file, "w", newline='') as f:
        f.write(sample_fasta_content)
    mock_find_python.return_value = "/usr/bin/python3"
    mock_subprocess.return_value = MagicMock(stdout="peptide,score\nP12345,5.0\n", stderr="")
    with patch("os.path.isfile", return_value=True):
        with patch("os.access", return_value=True):
            result = run_mhc2_prediction("consensus3", fasta_file, "HLA-DRB1*01:01", temp_dir)
    assert result == "peptide,score\nP12345,5.0\n"
    log_output = capture_logs.getvalue()
    assert "Running prediction for allele: HLA-DRB1*01:01, method: consensus3" in log_output

@patch("subprocess.run")
@patch("modules.mhc_ii.find_python_interpreter")
def test_run_mhc2_prediction_missing_tool(mock_find_python, mock_subprocess, sample_fasta_content, temp_dir, capture_logs):
    """Test run_mhc2_prediction with missing MHC-II tool."""
    fasta_file = temp_dir / "test.fasta"
    with open(fasta_file, "w", newline='') as f:
        f.write(sample_fasta_content)
    mock_find_python.return_value = "/usr/bin/python3"
    with patch("os.path.isfile", return_value=False):
        result = run_mhc2_prediction("consensus3", fasta_file, "HLA-DRB1*01:01", temp_dir)
    assert result is None
    log_output = capture_logs.getvalue()
    assert "MHC-II tool" in log_output
    assert "not found" in log_output

# Test parse_results
def test_parse_results_consensus3(capture_logs):
    """Test parse_results with consensus3 output."""
    response_text = "peptide,consensus_percentile_rank\nP12345,2.5\nQ67890,15.0\n"
    df = parse_results(response_text, "consensus3", "HLA-DRB1*01:01", score_threshold=10)
    assert not df.empty
    assert df["Method"].iloc[0] == "consensus3"
    assert df["Allele"].iloc[0] == "HLA-DRB1*01:01"
    assert df["score"].iloc[0] == 2.5
    assert len(df) == 1  # Only P12345 passes threshold
    log_output = capture_logs.getvalue()
    assert "Using score column: consensus_percentile_rank" in log_output

def test_parse_results_empty(capture_logs):
    """Test parse_results with empty response."""
    df = parse_results("", "consensus3", "HLA-DRB1*01:01", score_threshold=10)
    assert df.empty
    log_output = capture_logs.getvalue()
    assert "No results for HLA-DRB1*01:01" in log_output

# Test run_mhc2
@patch("modules.mhc_ii.get_valid_alleles")
@patch("modules.mhc_ii.run_mhc2_prediction")
@patch("modules.mhc_ii.validate_fasta")
def test_run_mhc2_success(mock_validate_fasta, mock_run_mhc2_prediction, mock_get_valid_alleles, sample_fasta_content, temp_dir, capture_logs):
    """Test run_mhc2 with successful predictions."""
    fasta_file = temp_dir / "test.fasta"
    with open(fasta_file, "w", newline='') as f:
        f.write(sample_fasta_content)
    output_file = temp_dir / "mhcii_out.csv"
    mock_validate_fasta.return_value = True
    mock_get_valid_alleles.return_value = ["HLA-DRB1*01:01"]
    mock_run_mhc2_prediction.return_value = "peptide,consensus_percentile_rank\nP12345,2.5\n"
    with patch("pandas.DataFrame.to_csv"):
        result = run_mhc2("cons", fasta_file, str(temp_dir), str(output_file), score_threshold=10)
    assert result == 0
    log_output = capture_logs.getvalue()
    assert "Starting run_mhc2 with method_code: cons" in log_output
    assert "Predictions complete" in log_output

@patch("modules.mhc_ii.validate_fasta")
def test_run_mhc2_invalid_fasta(mock_validate_fasta, temp_dir, capture_logs):
    """Test run_mhc2 with invalid FASTA file."""
    fasta_file = temp_dir / "test.fasta"
    output_file = temp_dir / "mhcii_out.csv"
    mock_validate_fasta.return_value = False
    result = run_mhc2("cons", fasta_file, str(temp_dir), str(output_file), score_threshold=10)
    assert result == 1
    log_output = capture_logs.getvalue()
    assert "Starting run_mhc2 with method_code: cons" in log_output

# Test main
@patch("sys.argv", ["mhc_ii.py", "-m", "cons", "-i", "test.fasta", "-o", "mhcii_out.csv", "-d", "output"])
@patch("modules.mhc_ii.run_mhc2")
def test_main_valid_args(mock_run_mhc2, temp_dir, capture_logs):
    """Test main function with valid arguments."""
    input_file = temp_dir / "test.fasta"
    input_file.touch()
    output_dir = temp_dir / "output"
    mock_run_mhc2.return_value = 0
    with pytest.raises(SystemExit) as exc:
        from modules.mhc_ii import main
        main()
    assert exc.value.code == 0
    log_output = capture_logs.getvalue()
    assert "Running mhc_ii.py with args" in log_output

@patch("sys.argv", ["mhc_ii.py", "-m", "invalid", "-i", "test.fasta", "-o", "mhcii_out.csv"])
def test_main_invalid_method(temp_dir, capture_logs):
    """Test main function with invalid method code."""
    input_file = temp_dir / "test.fasta"
    input_file.touch()
    with pytest.raises(SystemExit) as exc:
        from modules.mhc_ii import main
        main()
    assert exc.value.code == 1
    log_output = capture_logs.getvalue()
    assert "Invalid method code: invalid" in log_output

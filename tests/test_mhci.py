import pytest
import os
import sys
import subprocess
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock
from io import StringIO
import logging

# Debug sys.path
print("Initial sys.path:", sys.path)
sys.path.append(str(Path(__file__).resolve().parent.parent / "modules"))
print("Updated sys.path:", sys.path)

import mhc_i  # Import the module

# Sample FASTA content with amino acid sequences
SAMPLE_FASTA = """>pep1
MKFLV
>pep2
AILWY
"""

# Sample MHC-I tool output
SAMPLE_MHC_OUTPUT = """
Peptide,Allele,Score
MKFLV,HLA-A*01:01,1000
AILWY,HLA-A*01:01,6000
"""

@pytest.fixture
def temp_fasta_file(tmp_path):
    """Create a temporary FASTA file with amino acid sequences."""
    fasta_path = tmp_path / "peptides.fasta"
    with open(fasta_path, "w") as f:
        f.write(SAMPLE_FASTA)
    return fasta_path

@pytest.fixture
def setup_logging():
    """Set up logging for tests."""
    logger = logging.getLogger('VaxElan')
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    logger.addHandler(handler)
    yield
    logger.handlers = []

def test_print_status(capsys, setup_logging):
    """Test the print_status function for different status types."""
    mhc_i.print_status("Test info message", "info")
    captured = capsys.readouterr()
    assert "Test info message" in captured.out
    assert "\033[94m" in captured.out  # Blue color for info

    mhc_i.print_status("Test warning message", "warning")
    captured = capsys.readouterr()
    assert "Test warning message" in captured.out
    assert "\033[93m" in captured.out  # Yellow color for warning

    mhc_i.print_status("Test error message", "error")
    captured = capsys.readouterr()
    assert "Test error message" in captured.out
    assert "\033[91m" in captured.out  # Red color for error

@patch("psutil.Process")
def test_log_memory_usage(mock_process, setup_logging):
    """Test log_memory_usage function."""
    mock_mem_info = MagicMock(rss=100 * 1024 * 1024, vms=200 * 1024 * 1024)  # 100 MB RSS, 200 MB VMS
    mock_process.return_value.memory_info.return_value = mock_mem_info
    mem_usage = mhc_i.log_memory_usage()
    assert mem_usage == 100.0  # RSS in MB

@patch("psutil.virtual_memory")
@patch("psutil.Process")
@patch("builtins.open", new_callable=MagicMock)
@patch("os.makedirs")
@patch("subprocess.run")
@patch("os.path.isfile")
def test_run_mhc_predictor_success(mock_isfile, mock_subprocess, mock_makedirs, mock_open, mock_process, mock_virtual_memory, temp_fasta_file, tmp_path, setup_logging):
    """Test run_mhc_predictor with successful execution."""
    mock_isfile.return_value = True
    mock_subprocess.return_value = MagicMock(stdout=SAMPLE_MHC_OUTPUT, stderr="", returncode=0)
    mock_mem_info = MagicMock(rss=100 * 1024 * 1024, vms=200 * 1024 * 1024)
    mock_process.return_value.memory_info.return_value = mock_mem_info
    mock_virtual_memory.return_value.total = 8 * 1024 * 1024 * 1024  # 8 GB
    result = mhc_i.run_mhc_predictor(
        method="netmhccons",
        allele="HLA-A*01:01",
        length=9,
        input_file=str(temp_fasta_file),
        output_dir=str(tmp_path),
        timeout=3600
    )
    assert result == SAMPLE_MHC_OUTPUT
    mock_subprocess.assert_called_once()
    mock_makedirs.assert_called_once_with(os.path.join(str(tmp_path), "mhci_debug"), exist_ok=True)
    debug_file_path = os.path.join(str(tmp_path), "mhci_debug", "mhc1_HLA-A_01:01_9.txt")
    mock_open.assert_any_call(debug_file_path, "w")

@patch("psutil.virtual_memory")
@patch("psutil.Process")
@patch("builtins.open", new_callable=MagicMock)
@patch("os.makedirs")
@patch("subprocess.run")
@patch("os.path.isfile")
def test_run_mhc_predictor_timeout(mock_isfile, mock_subprocess, mock_makedirs, mock_open, mock_process, mock_virtual_memory, temp_fasta_file, tmp_path, setup_logging):
    """Test run_mhc_predictor with timeout."""
    mock_isfile.return_value = True
    mock_subprocess.side_effect = subprocess.TimeoutExpired(cmd="test", timeout=3600)
    mock_mem_info = MagicMock(rss=100 * 1024 * 1024, vms=200 * 1024 * 1024)
    mock_process.return_value.memory_info.return_value = mock_mem_info
    mock_virtual_memory.return_value.total = 8 * 1024 * 1024 * 1024  # 8 GB
    result = mhc_i.run_mhc_predictor(
        method="netmhccons",
        allele="HLA-A*01:01",
        length=9,
        input_file=str(temp_fasta_file),
        output_dir=str(tmp_path),
        timeout=3600
    )
    assert result is None
    mock_makedirs.assert_called_once_with(os.path.join(str(tmp_path), "mhci_debug"), exist_ok=True)
    debug_error_file_path = os.path.join(str(tmp_path), "mhci_debug", "mhc1_HLA-A_01:01_9.txt.error")
    mock_open.assert_any_call(debug_error_file_path, "w")

@patch("os.path.isfile")
def test_run_mhc_predictor_missing_tool(mock_isfile, temp_fasta_file, tmp_path, setup_logging):
    """Test run_mhc_predictor with missing MHC-I tool."""
    mock_isfile.side_effect = lambda x: x != str(mhc_i.MHC_I_TOOL)
    result = mhc_i.run_mhc_predictor(
        method="netmhccons",
        allele="HLA-A*01:01",
        length=9,
        input_file=str(temp_fasta_file),
        output_dir=str(tmp_path)
    )
    assert result is None

@patch("pandas.read_csv")
def test_parse_results_success(mock_read_csv, setup_logging):
    """Test parse_results with valid output."""
    mock_df = pd.DataFrame({
        "Peptide": ["MKFLV", "AILWY"],
        "Score": [1000, 6000]
    })
    mock_read_csv.return_value = mock_df
    df = mhc_i.parse_results(
        SAMPLE_MHC_OUTPUT,
        method="netmhccons",
        allele="HLA-A*01:01",
        length=9,
        score_threshold=5000
    )
    assert not df.empty
    assert list(df.columns[:3]) == ["Method", "Allele", "Peptide_Length"]
    assert df["Method"].iloc[0] == "netmhccons"
    assert df["Allele"].iloc[0] == "HLA-A*01:01"
    assert df["Peptide_Length"].iloc[0] == 9
    assert len(df) == 1  # Only MKFLV (score <= 5000)

def test_parse_results_empty(setup_logging):
    """Test parse_results with empty or invalid output."""
    df = mhc_i.parse_results("", method="netmhccons", allele="HLA-A*01:01", length=9)
    assert df.empty
    df = mhc_i.parse_results("ERROR: Invalid data", method="netmhccons", allele="HLA-A*01:01", length=9)
    assert df.empty

@patch("mhc_i.run_mhc_predictor")
@patch("mhc_i.parse_results")
@patch("pandas.DataFrame.to_csv")
@patch("os.makedirs")
@patch("psutil.virtual_memory")
@patch("os.cpu_count")
def test_run_mhc1_success(
    mock_cpu_count, mock_virtual_memory, mock_makedirs, mock_to_csv,
    mock_parse_results, mock_run_mhc_predictor, temp_fasta_file, tmp_path, setup_logging
):
    """Test run_mhc1 with successful predictions."""
    mock_cpu_count.return_value = 4
    mock_virtual_memory.return_value.total = 8 * 1024 * 1024 * 1024  # 8 GB
    mock_run_mhc_predictor.return_value = SAMPLE_MHC_OUTPUT
    mock_df = pd.DataFrame({
        "Method": ["netmhccons"],
        "Allele": ["HLA-A*01:01"],
        "Peptide_Length": [9],
        "Peptide": ["MKFLV"],
        "score": [1000]
    })
    mock_parse_results.return_value = mock_df
    output_file = tmp_path / "output.csv"
    status = mhc_i.run_mhc1(
        method_code="d",  # netmhccons
        input_file=str(temp_fasta_file),
        output_dir=str(tmp_path),
        output_file=str(output_file),
        score_threshold=5000
    )
    assert status == 0
    mock_run_mhc_predictor.assert_called()
    mock_to_csv.assert_called_once()

@patch("mhc_i.run_mhc_predictor")
@patch("os.makedirs")
@patch("psutil.virtual_memory")
@patch("os.cpu_count")
def test_run_mhc1_invalid_method(
    mock_cpu_count, mock_virtual_memory, mock_makedirs, mock_run_mhc_predictor, temp_fasta_file, tmp_path, setup_logging
):
    """Test run_mhc1 with invalid method code."""
    mock_cpu_count.return_value = 4
    mock_virtual_memory.return_value.total = 8 * 1024 * 1024 * 1024  # 8 GB
    status = mhc_i.run_mhc1(
        method_code="x",  # Invalid method
        input_file=str(temp_fasta_file),
        output_dir=str(tmp_path),
        output_file=str(tmp_path / "output.csv")
    )
    assert status == 1
    mock_run_mhc_predictor.assert_not_called()

@patch("mhc_i.run_mhc_predictor")
@patch("mhc_i.parse_results")
@patch("os.makedirs")
@patch("psutil.virtual_memory")
@patch("os.cpu_count")
def test_run_mhc1_all_fail(
    mock_cpu_count, mock_virtual_memory, mock_makedirs, mock_parse_results, mock_run_mhc_predictor, temp_fasta_file, tmp_path, setup_logging
):
    """Test run_mhc1 when all predictions fail."""
    mock_cpu_count.return_value = 4
    mock_virtual_memory.return_value.total = 8 * 1024 * 1024 * 1024  # 8 GB
    mock_run_mhc_predictor.return_value = None
    mock_parse_results.return_value = pd.DataFrame()
    status = mhc_i.run_mhc1(
        method_code="d",  # netmhccons
        input_file=str(temp_fasta_file),
        output_dir=str(tmp_path),
        output_file=str(tmp_path / "output.csv")
    )
    assert status == 2
    mock_run_mhc_predictor.assert_called()

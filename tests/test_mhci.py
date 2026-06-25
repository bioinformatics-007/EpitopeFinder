import pytest
import os
import sys
import subprocess
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock
from io import StringIO
import json

# Add the modules directory to sys.path
sys.path.append(str(Path(__file__).resolve().parent.parent / "modules"))
import mhc_i

SAMPLE_FASTA = """>pep1
MKFLV
>pep2
AILWY
"""

SAMPLE_MHC_OUTPUT = """peptide\tallele\tscore
MKFLV\tHLA-A*01:01\t1000
AILWY\tHLA-A*01:01\t6000
"""

@pytest.fixture
def temp_fasta_file(tmp_path):
    """Create a temporary FASTA file with amino acid sequences."""
    fasta_path = tmp_path / "peptides.fasta"
    with open(fasta_path, "w") as f:
        f.write(SAMPLE_FASTA)
    return str(fasta_path)

def test_create_iedb_json(temp_fasta_file):
    """Test create_iedb_json creates a valid JSON file."""
    method = "netmhcpan_ba"
    alleles = ["HLA-A*01:01", "HLA-B*07:02"]
    
    json_path = mhc_i.create_iedb_json(temp_fasta_file, method, alleles)
    assert os.path.exists(json_path)
    
    with open(json_path, 'r') as f:
        data = json.load(f)
        
    assert data["input_sequence_text_file_path"] == os.path.abspath(temp_fasta_file)
    assert data["alleles"] == "HLA-A*01:01,HLA-B*07:02"
    assert data["predictors"][0]["method"] == "netmhcpan_ba"
    
    os.remove(json_path)

@patch("subprocess.run")
def test_run_mhc1_success(mock_subprocess, temp_fasta_file, tmp_path):
    """Test run_mhc1 with successful execution."""
    mock_subprocess.return_value = MagicMock(returncode=0, stdout=SAMPLE_MHC_OUTPUT, stderr="")
    
    output_file = str(tmp_path / "mhci_out.csv")
    exit_code = mhc_i.run_mhc1(
        method_code="d",  # netmhcpan_ba
        input_file=temp_fasta_file,
        output_file=output_file,
        score_threshold=5000
    )
    
    assert exit_code == 0
    assert os.path.exists(output_file)
    
    # Check that it filtered out the score > 5000
    df = pd.read_csv(output_file)
    assert "Method" in df.columns
    assert "peptide" in df.columns
    assert len(df) == 1
    assert df["peptide"].iloc[0] == "MKFLV"
    assert df["Method"].iloc[0] == "netmhcpan_ba"

@patch("subprocess.run")
def test_run_mhc1_empty_output(mock_subprocess, temp_fasta_file, tmp_path):
    """Test run_mhc1 when tool returns no data."""
    mock_subprocess.return_value = MagicMock(returncode=0, stdout="  \n ", stderr="")
    
    output_file = str(tmp_path / "mhci_out.csv")
    exit_code = mhc_i.run_mhc1(
        method_code="d",
        input_file=temp_fasta_file,
        output_file=output_file,
        score_threshold=5000
    )
    
    assert exit_code == 0
    assert os.path.exists(output_file)
    
    df = pd.read_csv(output_file)
    assert "Method" in df.columns
    assert "peptide" in df.columns
    assert len(df) == 0

@patch("subprocess.run")
def test_run_mhc1_subprocess_error(mock_subprocess, temp_fasta_file, tmp_path):
    """Test run_mhc1 when tool fails."""
    mock_subprocess.return_value = MagicMock(returncode=1, stdout="", stderr="Tool crashed")
    
    output_file = str(tmp_path / "mhci_out.csv")
    exit_code = mhc_i.run_mhc1(
        method_code="d",
        input_file=temp_fasta_file,
        output_file=output_file,
        score_threshold=5000
    )
    
    assert exit_code == 1

def test_run_mhc1_invalid_method(temp_fasta_file, tmp_path):
    """Test run_mhc1 with invalid method code."""
    output_file = str(tmp_path / "mhci_out.csv")
    exit_code = mhc_i.run_mhc1(
        method_code="x",  # Invalid
        input_file=temp_fasta_file,
        output_file=output_file
    )
    
    assert exit_code == 1

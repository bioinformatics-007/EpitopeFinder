import sys
import os
import pytest
import tempfile
from unittest.mock import patch, MagicMock
from pathlib import Path

# Adjust import path as needed
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'modules')))
import bcell


@pytest.fixture
def fasta_file(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nMVLSPADKTNVKAA\n")
    return fasta

@pytest.fixture
def empty_fasta_file(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text(">seq1\n")
    return fasta

@pytest.fixture
def invalid_fasta_file(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    fasta.write_text("")
    return fasta

def test_read_fasta_valid(fasta_file):
    seq = bcell.read_fasta(str(fasta_file))
    assert seq == "MVLSPADKTNVKAA"

def test_read_fasta_empty_raises(empty_fasta_file):
    with pytest.raises(ValueError):
        bcell.read_fasta(str(empty_fasta_file))

def test_read_fasta_missing_file():
    with pytest.raises(FileNotFoundError):
        bcell.read_fasta("nonexistent.fasta")

@patch("bcell.subprocess.run")
def test_query_iedb_success(mock_run, fasta_file):
    mock_proc = MagicMock()
    mock_proc.stdout = "Predicted peptides\n1 5 MVLSP 5.0 0.9\n"
    mock_proc.stderr = ""
    mock_proc.returncode = 0
    mock_run.return_value = mock_proc

    # Temporarily override IEDB_TOOL_PATH to avoid file check
    original_path = bcell.IEDB_TOOL_PATH
    bcell.IEDB_TOOL_PATH = str(fasta_file)  # any existing file path
    try:
        output = bcell.query_iedb("MVLSPADKTNVKAA", method="Emini")
        assert "Predicted peptides" in output
    finally:
        bcell.IEDB_TOOL_PATH = original_path

@patch("bcell.subprocess.run")
def test_query_iedb_calledprocesserror(mock_run, fasta_file):
    mock_run.side_effect = bcell.subprocess.CalledProcessError(
        1, "cmd", output="", stderr="error"
    )

    # Override IEDB_TOOL_PATH to bypass file check
    original_path = bcell.IEDB_TOOL_PATH
    bcell.IEDB_TOOL_PATH = str(fasta_file)
    try:
        with pytest.raises(RuntimeError):
            bcell.query_iedb("MVLSPADKTNVKAA")
    finally:
        bcell.IEDB_TOOL_PATH = original_path

def test_parse_output_correctness():
    sample_output = """
Predicted peptides
1 5 MVLSP 5.0 0.9
Position  Residue  Start  End  Peptide  Score
1 M 1 5 MVLSP 0.9
"""
    peptides, epitopes = bcell.parse_output(sample_output)
    assert len(peptides) == 1
    assert peptides[0][2] == "MVLSP"
    assert len(epitopes) == 1
    assert epitopes[0][1] == "M"

@patch("bcell.query_iedb")
@patch("bcell.read_fasta")
def test_run_bcell_success(mock_read_fasta, mock_query_iedb, tmp_path):
    mock_read_fasta.return_value = "MVLSPADKTNVKAA"
    mock_query_iedb.return_value = (
        "Predicted peptides\n1 5 MVLSP 5.0 0.9\n"
        "Position  Residue  Start  End  Peptide  Score\n"
        "1 M 1 5 MVLSP 0.9\n"
    )

    output_dir = tmp_path / "output"
    ret = bcell.run_bcell(
        input_fasta="dummy.fasta",
        output_dir=str(output_dir),
        output_file="result.csv",
        method="Emini"
    )
    assert ret == 0
    output_file = output_dir / "result.csv"
    assert output_file.exists()
    content = output_file.read_text()
    assert "Peptide" in content
    assert "Epitope" in content

def test_main_cli_invokes_run_bcell(monkeypatch, tmp_path):
    test_fasta = tmp_path / "test.fasta"
    test_fasta.write_text(">seq1\nMVLSPADKTNVKAA\n")

    # Patch run_bcell to just check call args
    called = {}

    def fake_run_bcell(input_fasta, output_dir, output_file, method):
        called["input_fasta"] = input_fasta
        called["output_dir"] = output_dir
        called["output_file"] = output_file
        called["method"] = method
        return 0

    monkeypatch.setattr(bcell, "run_bcell", fake_run_bcell)

    sys.argv = ["prog", "-i", str(test_fasta), "-d", "outdir", "-o", "myout.csv", "-m", "Emini"]
    bcell.main()

    assert called["input_fasta"] == str(test_fasta)
    assert called["output_dir"] == "outdir"
    assert called["output_file"] == "myout.csv"
    assert called["method"] == "Emini"

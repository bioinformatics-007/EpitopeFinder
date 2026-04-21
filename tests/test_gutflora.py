import pytest
import os
import sys
import subprocess
from pathlib import Path
import tempfile
import pandas as pd
from unittest.mock import patch, mock_open, MagicMock
import xml.etree.ElementTree as ET
import logging
import re

# Add parent directory to sys.path to allow absolute imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules.gutflora import parse_fasta, write_temp_fasta, run_blastp, parse_blast_output, save_output, cleanup_files, run_gutflora, print_status, ROOT_DIR, GUT_DB

# Configure logging for tests
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("GutfloraTest")

@pytest.fixture
def temp_fasta_file(tmp_path):
    """Create a temporary FASTA file for testing."""
    fasta_content = """>prot1 Test Protein 1
ACDEFGHIKLMNPQRSTVWY
>prot2 Test Protein 2
ACDEFGHIKLMNPQRSTVWY*"""
    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(fasta_content)
    return fasta_path

@pytest.fixture
def empty_fasta_file(tmp_path):
    """Create an empty FASTA file for testing."""
    fasta_path = tmp_path / "empty.fasta"
    fasta_path.write_text("")
    return fasta_path

@pytest.fixture
def invalid_fasta_file(tmp_path):
    """Create a FASTA file with invalid sequence for testing."""
    fasta_content = """>prot1 Test Protein 1
XYZ123"""
    fasta_path = tmp_path / "invalid.fasta"
    fasta_path.write_text(fasta_content)
    return fasta_path

@pytest.fixture
def blast_xml_content():
    """Mock BLAST XML content."""
    return """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_id>subject1</Hit_id>
          <Hit_def>Subject Protein 1</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_identity>80</Hsp_identity>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_evalue>1e-5</Hsp_evalue>
              <Hsp_bit-score>200</Hsp_bit-score>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_id>subject2</Hit_id>
          <Hit_def>Subject Protein 2</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_identity>20</Hsp_identity>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_evalue>0.1</Hsp_evalue>
              <Hsp_bit-score>50</Hsp_bit-score>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""

def test_parse_fasta(temp_fasta_file, caplog):
    """Test parsing a valid FASTA file."""
    caplog.set_level(logging.DEBUG)
    sequences = parse_fasta(temp_fasta_file)
    assert len(sequences) == 2
    assert sequences[0] == ("prot1", "Test Protein 1", "ACDEFGHIKLMNPQRSTVWY")
    assert sequences[1] == ("prot2", "Test Protein 2", "ACDEFGHIKLMNPQRSTVWY*")
    assert "Parsed sequence: prot1" in caplog.text

def test_parse_fasta_empty(empty_fasta_file, caplog):
    """Test parsing an empty FASTA file."""
    caplog.set_level(logging.ERROR)
    with pytest.raises(ValueError, match="Empty FASTA file"):
        parse_fasta(empty_fasta_file)
    assert "FASTA file is empty" in caplog.text

def test_parse_fasta_invalid_sequence(invalid_fasta_file, caplog):
    """Test parsing a FASTA file with invalid sequence."""
    caplog.set_level(logging.ERROR)
    with pytest.raises(ValueError, match="Invalid sequence"):
        parse_fasta(invalid_fasta_file)
    assert "Invalid sequence for >prot1 Test Protein 1" in caplog.text

def test_parse_fasta_no_header(tmp_path, caplog):
    """Test parsing a FASTA file with sequence data before header."""
    fasta_path = tmp_path / "no_header.fasta"
    fasta_path.write_text("ACDEFGHIKLMNPQRSTVWY\n>prot1 Test Protein")
    caplog.set_level(logging.ERROR)
    with pytest.raises(ValueError, match="Sequence data found before header"):
        parse_fasta(fasta_path)
    assert "Sequence data found before header" in caplog.text

def test_parse_fasta_no_sequence(tmp_path, caplog):
    """Test parsing a FASTA file with no sequence data."""
    fasta_path = tmp_path / "no_sequence.fasta"
    fasta_path.write_text(">prot1 Test Protein")
    caplog.set_level(logging.ERROR)
    with pytest.raises(ValueError, match="No sequence for header"):
        parse_fasta(fasta_path)
    assert "No sequence found for header" in caplog.text

def test_write_temp_fasta(tmp_path, caplog):
    """Test writing a temporary FASTA file."""
    caplog.set_level(logging.INFO)
    temp_file = tmp_path / "temp.fasta"
    write_temp_fasta("prot1", "Test Protein", "ACDEFGHIKLMNPQRSTVWY", temp_file)
    assert temp_file.exists()
    content = temp_file.read_text()
    assert content == ">prot1 Test Protein\nACDEFGHIKLMNPQRSTVWY\n"
    assert f"Created temporary FASTA file: {temp_file}" in caplog.text

def test_write_temp_fasta_error(tmp_path, caplog):
    """Test error handling in write_temp_fasta."""
    caplog.set_level(logging.ERROR)
    temp_file = tmp_path / "unwritable" / "temp.fasta"
    with pytest.raises(Exception):
        write_temp_fasta("prot1", "Test Protein", "ACDEFGHIKLMNPQRSTVWY", temp_file)
    assert "Error writing temporary FASTA file" in caplog.text

@patch("subprocess.run")
def test_run_blastp_success(mock_run, tmp_path, caplog):
    """Test running BLASTP successfully."""
    caplog.set_level(logging.INFO)
    db_path = tmp_path / "db"
    (db_path.with_suffix(".phr")).touch()
    query_file = "query.fasta"
    out_xml = "blast.xml"
    mock_run.return_value = MagicMock(returncode=0, stdout="BLASTP completed", stderr="")
    status = run_blastp(db_path, query_file, out_xml, tmp_path)
    assert status == 0
    assert "Running BLASTP" in caplog.text
    assert "BLASTP STDOUT: BLASTP completed" in caplog.text

@patch("subprocess.run")
def test_run_blastp_db_not_found(mock_run, tmp_path, caplog):
    """Test BLASTP when database is not found."""
    caplog.set_level(logging.ERROR)
    db_path = tmp_path / "db"
    query_file = "query.fasta"
    out_xml = "blast.xml"
    with pytest.raises(FileNotFoundError, match="BLAST database not found"):
        run_blastp(db_path, query_file, out_xml, tmp_path)
    assert "BLAST database not found" in caplog.text

@patch("subprocess.run")
def test_run_blastp_failure(mock_run, tmp_path, caplog):
    """Test BLASTP failure."""
    caplog.set_level(logging.ERROR)
    db_path = tmp_path / "db"
    (db_path.with_suffix(".phr")).touch()
    query_file = "query.fasta"
    out_xml = "blast.xml"
    mock_run.side_effect = subprocess.CalledProcessError(1, ["blastp"], stderr="BLASTP error")
    status = run_blastp(db_path, query_file, out_xml, tmp_path)
    assert status == 1
    assert "BLASTP failed" in caplog.text
    assert "STDERR: BLASTP error" in caplog.text

@patch("xml.etree.ElementTree.parse")
@patch("os.path.exists")
def test_parse_blast_output(mock_exists, mock_parse, blast_xml_content, caplog):
    """Test parsing BLAST XML output."""
    caplog.set_level(logging.INFO)
    mock_exists.return_value = True
    mock_tree = MagicMock()
    mock_root = ET.fromstring(blast_xml_content)
    mock_tree.getroot.return_value = mock_root
    mock_parse.return_value = mock_tree
    xml_path = "blast.xml"
    matches = parse_blast_output(xml_path)
    assert len(matches) == 1  # Only one match with identity >= 30% and evalue <= 1e-2
    assert matches[0]["subject_id"] == "subject1"
    assert matches[0]["identity_percent"] == 80.0
    assert matches[0]["evalue"] == 1e-5
    assert matches[0]["bit_score"] == 200
    assert f"Parsed 1 matches from {xml_path}" in caplog.text

def test_parse_blast_output_not_found(tmp_path, caplog):
    """Test parsing non-existent BLAST output."""
    caplog.set_level(logging.ERROR)
    xml_path = tmp_path / "nonexistent.xml"
    matches = parse_blast_output(xml_path)
    assert matches == []
    assert f"BLAST output file not found: {xml_path}" in caplog.text

def test_save_output(tmp_path, caplog):
    """Test saving results to CSV."""
    caplog.set_level(logging.INFO)
    output_file = tmp_path / "output.csv"
    results = [
        ("prot1", "Test Protein 1", [
            {"subject_id": "subject1", "subject_def": "Subject 1", "identity_percent": 80.0, "evalue": 1e-5, "bit_score": 200}
        ]),
        ("prot2", "Test Protein 2", [])
    ]
    save_output(output_file, results)
    assert output_file.exists()
    df = pd.read_csv(output_file)
    assert len(df) == 2
    assert df.iloc[0]["Sequence_ID"] == "prot1"
    assert df.iloc[0]["Cross_Reactivity_Detected"] == True
    assert "subject1 (Subject 1) - Identity: 80.0% | E-value: 1e-05 | Bit-score: 200" in df.iloc[0]["Matches"]
    assert df.iloc[1]["Sequence_ID"] == "prot2"
    assert df.iloc[1]["Cross_Reactivity_Detected"] == False
    assert pd.isna(df.iloc[1]["Matches"])
    assert f"Results saved to {output_file}" in caplog.text

def test_cleanup_files(tmp_path, caplog):
    """Test cleaning up temporary files."""
    caplog.set_level(logging.INFO)
    temp_file = tmp_path / "temp.fasta"
    temp_file.write_text("test content")
    cleanup_files([str(temp_file)])
    assert not temp_file.exists()
    assert f"Cleaned up {temp_file}" in caplog.text

def test_cleanup_files_not_found(tmp_path, caplog):
    """Test cleaning up non-existent files."""
    caplog.set_level(logging.INFO)
    temp_file = tmp_path / "ntest.fasta"
    cleanup_files([str(temp_file)])
    # No log assertion, as cleanup_files may not log for non-existent files

@patch("modules.gutflora.parse_fasta")
@patch("modules.gutflora.write_temp_fasta")
@patch("modules.gutflora.run_blastp")
@patch("modules.gutflora.parse_blast_output")
@patch("modules.gutflora.save_output")
@patch("modules.gutflora.cleanup_files")
@patch("os.path.exists")
def test_run_gutflora_check(mock_exists, mock_cleanup, mock_save, mock_parse_blast, mock_run_blastp, mock_write_fasta, mock_parse_fasta, temp_fasta_file, tmp_path, caplog):
    """Test the main execution of run_gutflora function."""
    caplog.set_level(logging.DEBUG)
    mock_parse_fasta.return_value = [("prot1", "Test Protein", "ACDEFGHIKLMNPQRSTVWY")]
    mock_run_blastp.return_value = 0
    mock_parse_blast.return_value = [{"subject_id": "subject1", "subject_def": "Subject 1", "identity_percent": 80.0, "evalue": 1e-5, "bit_score": 200}]
    mock_write_fasta.return_value = None
    mock_exists.return_value = True
    output_file = tmp_path / "output.csv"
    gut_db = tmp_path / "gut_db"
    (gut_db.with_suffix(".phr")).touch()
    
    status = run_gutflora(temp_fasta_file, output_file, gut_db)
    assert status == 0
    assert "Parsing FASTA file" in caplog.text
    assert "Found 1 sequences" in caplog.text
    assert "Processing sequence 1: prot1" in caplog.text
    assert "Cross-reactivity detected" in caplog.text
    mock_write_fasta.assert_called()
    mock_run_blastp.assert_called()
    mock_parse_blast.assert_called()
    mock_save.assert_called()
    mock_cleanup.assert_called()

@patch("modules.gutflora.parse_fasta")
def test_run_gutflora_input_not_found(mock_parse_fasta, tmp_path, caplog):
    """Test run_gutflora with missing input file."""
    caplog.set_level(logging.ERROR)
    input_file = tmp_path / "nonexistent.fasta"
    output_file = tmp_path / "output.csv"
    status = run_gutflora(input_file, output_file)
    assert status == 1
    assert "Input FASTA file not found" in caplog.text

@patch("modules.gutflora.parse_fasta")
def test_run_gutflora_db_not_found(mock_parse_fasta, tmp_path, caplog):
    """Test run_gutflora with missing database."""
    caplog.set_level(logging.ERROR)
    input_file = tmp_path / "test.fasta"
    input_file.touch()
    output_file = tmp_path / "output.csv"
    gut_db = tmp_path / "gut_db"
    mock_parse_fasta.return_value = [("prot1", "Test Protein", "ACDEFGHIKLMNPQRSTVWY")]
    status = run_gutflora(input_file, output_file, gut_db)
    assert status == 1
    assert "Gut microbiome database not found" in caplog.text

def test_print_status(capfd):
    """Test print_status function for different status levels."""
    print_status("Test info message", "info")
    out, _ = capfd.readouterr()
    assert "Test info message" in out

    print_status("Test warning message", "warning")
    out, _ = capfd.readouterr()
    assert "Test warning message" in out

    print_status("Test error message", "error")
    out, _ = capfd.readouterr()
    assert "Test error message" in out

    print_status("Test success message", "success")
    out, _ = capfd.readouterr()
    assert "Test success message" in out

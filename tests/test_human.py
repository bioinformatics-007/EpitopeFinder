import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'modules')))
import pytest
import subprocess
import tempfile
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import logging
from human_1 import validate_fasta_file, validate_paths, check_blastp, parse_blast_xml, run_human, SWISSPROT_DB

# Configure logging for tests
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("TestHumanBlast")

@pytest.fixture
def temp_fasta_file():
    """Create a temporary FASTA file for testing."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">seq1\nATCG\n>seq2\nGCTA\n")
    yield f.name
    os.unlink(f.name)

@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield Path(temp_dir)

@pytest.fixture
def mock_blast_db(tmp_path):
    """Create mock BLAST database files."""
    db_path = tmp_path / "swissprot"
    for ext in [".phr", ".pin", ".psq"]:
        (db_path.with_suffix(ext)).touch()
    return db_path

def test_validate_fasta_file(temp_fasta_file):
    """Test validate_fasta_file function."""
    sequences = validate_fasta_file(temp_fasta_file)
    assert len(sequences) == 2
    assert sequences[0].id == "seq1"
    assert str(sequences[0].seq) == "ATCG"
    assert sequences[1].id == "seq2"
    assert str(sequences[1].seq) == "GCTA"

def test_validate_fasta_file_empty(tmp_path):
    """Test validate_fasta_file with an empty file."""
    empty_fasta = tmp_path / "empty.fasta"
    empty_fasta.write_text("")
    with pytest.raises(ValueError, match="FASTA file is empty or contains no valid sequences"):
        validate_fasta_file(empty_fasta)

def test_validate_fasta_file_invalid(tmp_path):
    """Test validate_fasta_file with an invalid FASTA file."""
    invalid_fasta = tmp_path / "invalid.fasta"
    invalid_fasta.write_text("not a fasta file")
    with pytest.raises(ValueError, match=r"Invalid FASTA file.*"):
        validate_fasta_file(invalid_fasta)

def test_validate_paths(temp_fasta_file, tmp_path, mock_blast_db):
    """Test validate_paths function."""
    input_fasta = Path(temp_fasta_file)
    output_file = tmp_path / "output.txt"
    db_path = mock_blast_db
    output_dir = tmp_path / "test_output"

    # Test with default output_dir
    input_fasta, output_file, db_path, output_dir = validate_paths(input_fasta, output_file, db_path, None)
    assert input_fasta == Path(temp_fasta_file).resolve()
    assert output_file == output_dir / "combined.txt"
    assert db_path == mock_blast_db.resolve()
    assert output_dir == tmp_path / "human_blast"

    # Test with specified output_dir
    input_fasta, output_file, db_path, output_dir = validate_paths(input_fasta, output_file, db_path, tmp_path)
    assert output_dir == tmp_path / "human_blast"

def test_validate_paths_missing_fasta(tmp_path, mock_blast_db):
    """Test validate_paths with missing FASTA file."""
    with pytest.raises(FileNotFoundError, match="Input FASTA file not found"):
        validate_paths(tmp_path / "nonexistent.fasta", tmp_path / "output.txt", mock_blast_db, None)

def test_validate_paths_missing_db_files(tmp_path):
    """Test validate_paths with missing BLAST database files."""
    input_fasta = tmp_path / "test.fasta"
    input_fasta.write_text(">seq1\nATCG\n")
    db_path = tmp_path / "swissprot"
    (db_path.with_suffix(".phr")).touch()  # Only one DB file
    with pytest.raises(FileNotFoundError, match=r"Database files missing: .pin, .psq"):
        validate_paths(input_fasta, tmp_path / "output.txt", db_path, None)

@patch("subprocess.run")
def test_check_blastp(mock_subprocess_run):
    """Test check_blastp function."""
    mock_subprocess_run.return_value = MagicMock(stdout="/usr/bin/blastp\n")
    with patch("os.access", return_value=True):
        blastp_path = check_blastp()
        assert blastp_path == "/usr/bin/blastp"
        mock_subprocess_run.assert_called_once_with(
            ["which", "blastp"], capture_output=True, text=True, check=True
        )

@patch("subprocess.run")
def test_check_blastp_not_found(mock_subprocess_run):
    """Test check_blastp when blastp is not found."""
    mock_subprocess_run.side_effect = subprocess.CalledProcessError(1, ["which", "blastp"])
    with pytest.raises(FileNotFoundError, match="blastp not found"):
        check_blastp()

@patch("subprocess.run")
def test_check_blastp_not_executable(mock_subprocess_run):
    """Test check_blastp when blastp is not executable."""
    mock_subprocess_run.return_value = MagicMock(stdout="/usr/bin/blastp\n")
    with patch("os.access", return_value=False):
        with pytest.raises(PermissionError, match="blastp executable is not executable"):
            check_blastp()

def test_parse_blast_xml(tmp_path):
    """Test parse_blast_xml function."""
    xml_content = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>seq1</BlastOutput_query-def>
  <BlastOutput_query-len>100</BlastOutput_query-len>
  <BlastOutput_db>swissprot</BlastOutput_db>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>1e-5</Parameters_expect>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>F</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_query-def>seq1</Iteration_query-def>
      <Iteration_hits>
        <Hit>
          <Hit_id>hit1</Hit_id>
          <Hit_hsps>
            <Hsp>
              <Hsp_identity>90</Hsp_identity>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_gaps>2</Hsp_gaps>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>100</Hsp_hit-to>
              <Hsp_evalue>1e-10</Hsp_evalue>
              <Hsp_score>200</Hsp_score>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""
    xml_file = tmp_path / "blast_results.xml"
    xml_file.write_text(xml_content)

    df = parse_blast_xml(xml_file, identity_thresh=80.0)
    assert len(df) == 1
    assert df.iloc[0]["qseqid"] == "seq1"
    assert df.iloc[0]["sseqid"] == "hit1"
    assert df.iloc[0]["identity"] == 90.0
    assert df.iloc[0]["length"] == 100
    assert df.iloc[0]["mismatch"] == 10
    assert df.iloc[0]["gapopen"] == 2
    assert df.iloc[0]["qstart"] == 1
    assert df.iloc[0]["qend"] == 100
    assert df.iloc[0]["sstart"] == 1
    assert df.iloc[0]["send"] == 100
    assert df.iloc[0]["evalue"] == 1e-10
    assert df.iloc[0]["bitscore"] == 200

def test_parse_blast_xml_below_threshold(tmp_path):
    """Test parse_blast_xml with hits below identity threshold."""
    xml_content = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>seq1</BlastOutput_query-def>
  <BlastOutput_query-len>100</BlastOutput_query-len>
  <BlastOutput_db>swissprot</BlastOutput_db>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>1e-5</Parameters_expect>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>F</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_query-def>seq1</Iteration_query-def>
      <Iteration_hits>
        <Hit>
          <Hit_id>hit1</Hit_id>
          <Hit_hsps>
            <Hsp>
              <Hsp_identity>20</Hsp_identity>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_gaps>2</Hsp_gaps>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>100</Hsp_hit-to>
              <Hsp_evalue>1e-5</Hsp_evalue>
              <Hsp_score>50</Hsp_score>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""
    xml_file = tmp_path / "blast_results.xml"
    xml_file.write_text(xml_content)

    df = parse_blast_xml(xml_file, identity_thresh=30.0)
    assert df.empty

@patch("subprocess.run")
@patch("human_1.check_blastp")
@patch("human_1.SeqIO.write")
@patch("human_1.parse_blast_xml")
@patch("pandas.DataFrame.to_csv")
@patch("pandas.DataFrame.to_excel")
def test_run_human(mock_to_excel, mock_to_csv, mock_parse_blast_xml, mock_seqio_write, mock_check_blastp, mock_subprocess_run, temp_fasta_file, mock_blast_db, tmp_path):
    """Test run_human function with mocked BLASTP and file operations."""
    # Mock check_blastp to return a path
    mock_check_blastp.return_value = "/usr/bin/blastp"

    # Mock subprocess.run to simulate successful BLASTP execution
    mock_subprocess_run.return_value = MagicMock(stdout="BLASTP success", stderr="")

    # Mock parse_blast_xml to return a sample DataFrame
    sample_df = pd.DataFrame({
        "qseqid": ["seq1"],
        "sseqid": ["hit1"],
        "identity": [90.0],
        "length": [100],
        "mismatch": [10],
        "gapopen": [2],
        "qstart": [1],
        "qend": [100],
        "sstart": [1],
        "send": [100],
        "evalue": [1e-10],
        "bitscore": [200]
    })
    mock_parse_blast_xml.return_value = sample_df

    # Mock os.access to allow file operations
    with patch("os.access", return_value=True):
        exit_code = run_human(
            fasta_file=temp_fasta_file,
            evalue_thresh=1e-5,
            identity_thresh=30.0,
            output_file=tmp_path / "output.txt",
            outfmt=5,
            output_dir=tmp_path
        )
    
    assert exit_code == 0
    mock_check_blastp.assert_called_once()
    mock_subprocess_run.assert_called_once()
    mock_seqio_write.assert_called_once()
    mock_parse_blast_xml.assert_called_once()
    mock_to_csv.assert_called()
    mock_to_excel.assert_called()

@patch("subprocess.run")
@patch("human_1.check_blastp")
def test_run_human_blastp_failure(mock_check_blastp, mock_subprocess_run, temp_fasta_file, mock_blast_db, tmp_path):
    """Test run_human when BLASTP fails."""
    mock_check_blastp.return_value = "/usr/bin/blastp"
    mock_subprocess_run.side_effect = subprocess.CalledProcessError(1, ["blastp"], stderr="BLASTP error")
    with patch("os.access", return_value=True):
        exit_code = run_human(
            fasta_file=temp_fasta_file,
            evalue_thresh=1e-5,
            identity_thresh=30.0,
            output_file=tmp_path / "output.txt",
            outfmt=5,
            output_dir=tmp_path
        )
    assert exit_code == 1
    mock_check_blastp.assert_called_once()
    mock_subprocess_run.assert_called_once()

if __name__ == "__main__":
    pytest.main([__file__])

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
import tempfile
import subprocess
from unittest.mock import patch, MagicMock
from pathlib import Path
from modules.psortb import (
    validate_fasta,
    read_fasta,
    run_psortb_tool,
    parse_psortb_output,
    run_psortb,
    print_status,
    find_psortb_executable,
    check_blast,
    setup_bioperl
)

# Define ROOT_DIR for tests
ROOT_DIR = Path(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Sample amino acid sequence for testing
SAMPLE_FASTA = """>seq1
MKLYNLKDHNEQVSFAKEFLQDKHTFDIDWGKLSVDEHLKSYQQLAEHLTKQKSS
>seq2
MRLKTIQFDHLEAIEKAFKTAEAWYEAGDVQNNVGLSGNAFKKALDDILTKVNYD
"""

# Invalid FASTA files for testing
INVALID_FASTA_EMPTY = ""
INVALID_FASTA_NO_HEADER = "MKLYNLKDHNEQVSFAKEFLQDKHTFDIDWGKLSVDEHLKSYQQLAEHLTKQKSS"
INVALID_FASTA_NO_SEQUENCE = ">seq1\n"

# Mock PSORTb output for testing
MOCK_PSORTB_OUTPUT = """SeqID\tFinal_Localization\tScore
seq1\tCytoplasmic\t0.95
seq2\tMembrane\t0.90
"""

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        yield Path(tmpdirname)

@pytest.fixture
def sample_fasta_file(temp_dir):
    """Create a temporary FASTA file with sample amino acid sequences."""
    fasta_path = temp_dir / "test.fasta"
    with open(fasta_path, "w") as f:
        f.write(SAMPLE_FASTA)
    return fasta_path

@pytest.fixture
def invalid_fasta_files(temp_dir):
    """Create temporary invalid FASTA files."""
    empty_fasta = temp_dir / "empty.fasta"
    no_header_fasta = temp_dir / "no_header.fasta"
    no_sequence_fasta = temp_dir / "no_sequence.fasta"

    with open(empty_fasta, "w") as f:
        f.write(INVALID_FASTA_EMPTY)
    with open(no_header_fasta, "w") as f:
        f.write(INVALID_FASTA_NO_HEADER)
    with open(no_sequence_fasta, "w") as f:
        f.write(INVALID_FASTA_NO_SEQUENCE)

    return {
        "empty": empty_fasta,
        "no_header": no_header_fasta,
        "no_sequence": no_sequence_fasta
    }

def test_validate_fasta_valid(sample_fasta_file):
    """Test validate_fasta with a valid FASTA file."""
    assert validate_fasta(sample_fasta_file) is True

def test_validate_fasta_invalid(invalid_fasta_files):
    """Test validate_fasta with various invalid FASTA files."""
    with pytest.raises(ValueError, match="FASTA file is empty"):
        validate_fasta(invalid_fasta_files["empty"])
    
    with pytest.raises(ValueError, match="FASTA file contains no valid headers"):
        validate_fasta(invalid_fasta_files["no_header"])
    
    with pytest.raises(ValueError, match="FASTA file contains no valid sequence data"):
        validate_fasta(invalid_fasta_files["no_sequence"])

def test_read_fasta(sample_fasta_file):
    """Test read_fasta with a valid FASTA file."""
    sequences = read_fasta(sample_fasta_file)
    assert len(sequences) == 2
    assert sequences[0][0] == "seq1"
    assert sequences[0][1] == "MKLYNLKDHNEQVSFAKEFLQDKHTFDIDWGKLSVDEHLKSYQQLAEHLTKQKSS"
    assert sequences[1][0] == "seq2"
    assert sequences[1][1] == "MRLKTIQFDHLEAIEKAFKTAEAWYEAGDVQNNVGLSGNAFKKALDDILTKVNYD"

def test_read_fasta_file_not_found(temp_dir):
    """Test read_fasta with a non-existent file."""
    non_existent_file = temp_dir / "non_existent.fasta"
    with pytest.raises(FileNotFoundError, match="Input FASTA file not found"):
        read_fasta(non_existent_file)

@patch("modules.psortb.shutil.which")
@patch("modules.psortb.subprocess.run")
@patch("modules.psortb.os.path.isfile")
@patch("modules.psortb.os.makedirs")
@patch("modules.psortb.check_blast")
@patch("modules.psortb.setup_bioperl")
@patch("modules.psortb.find_psortb_executable")
def test_run_psortb_tool(mock_find_psortb_executable, mock_setup_bioperl, mock_check_blast, mock_makedirs, mock_isfile, mock_subprocess_run, mock_shutil_which, sample_fasta_file, temp_dir):
    """Test run_psortb_tool with mocked dependencies."""
    # Mock dependencies
    mock_shutil_which.return_value = "/usr/bin/blastp"
    mock_check_blast.return_value = True
    mock_isfile.return_value = True
    mock_subprocess_run.return_value = MagicMock(stdout="PSORTb ran successfully", stderr="")
    mock_find_psortb_executable.return_value = "/usr/bin/psortb"
    
    # Patch the global PSORTB_EXECUTABLE if used
    with patch("modules.psortb.PSORTB_EXECUTABLE", "/usr/bin/psortb"):
        # Run the function
        result = run_psortb_tool(sample_fasta_file, temp_dir, organism_type="n")
    
        # Verify command
        expected_command = [
            "/usr/bin/psortb",
            "-i", str(sample_fasta_file),
            "-r", str(temp_dir),
            "-n",
            "-o", "long"
        ]
        mock_subprocess_run.assert_called_once_with(
            expected_command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        assert result == "PSORTb ran successfully"

def test_run_psortb_tool_invalid_organism(sample_fasta_file, temp_dir):
    """Test run_psortb_tool with invalid organism type."""
    with pytest.raises(ValueError, match="Invalid organism type: x"):
        run_psortb_tool(sample_fasta_file, temp_dir, organism_type="x")

def test_run_psortb_tool_file_not_found(temp_dir):
    """Test run_psortb_tool with non-existent FASTA file."""
    non_existent_file = temp_dir / "non_existent.fasta"
    with pytest.raises(FileNotFoundError, match="Input FASTA file not found"):
        run_psortb_tool(non_existent_file, temp_dir)

def test_parse_psortb_output(temp_dir):
    """Test parse_psortb_output with a mock PSORTb output file."""
    output_file = temp_dir / "mock_psortb_output.txt"
    with open(output_file, "w") as f:
        f.write(MOCK_PSORTB_OUTPUT)
    
    header, predictions = parse_psortb_output(output_file)
    
    assert header == ["SeqID", "Final_Localization", "Score"]
    assert len(predictions) == 2
    assert predictions[0] == {"SeqID": "seq1", "Final_Localization": "Cytoplasmic", "Score": "0.95"}
    assert predictions[1] == {"SeqID": "seq2", "Final_Localization": "Membrane", "Score": "0.90"}

def test_parse_psortb_output_file_not_found(temp_dir):
    """Test parse_psortb_output with a non-existent file."""
    non_existent_file = temp_dir / "non_existent.txt"
    with pytest.raises(FileNotFoundError, match="PSORTb output file not found"):
        parse_psortb_output(non_existent_file)

@patch("modules.psortb.run_psortb_tool")
@patch("modules.psortb.find_latest_output")
@patch("modules.psortb.read_fasta")
@patch("modules.psortb.parse_psortb_output")
def test_run_psortb(mock_parse_output, mock_read_fasta, mock_find_latest_output, mock_run_psortb_tool, sample_fasta_file, temp_dir):
    """Test run_psortb with mocked dependencies."""
    # Mock dependencies
    mock_read_fasta.return_value = [
        ("seq1", "MKLYNLKDHNEQVSFAKEFLQDKHTFDIDWGKLSVDEHLKSYQQLAEHLTKQKSS"),
        ("seq2", "MRLKTIQFDHLEAIEKAFKTAEAWYEAGDVQNNVGLSGNAFKKALDDILTKVNYD")
    ]
    mock_run_psortb_tool.return_value = "PSORTb ran successfully"
    mock_find_latest_output.return_value = temp_dir / "mock_psortb_output.txt"
    mock_parse_output.return_value = (
        ["SeqID", "Final_Localization", "Score"],
        [
            {"SeqID": "seq1", "Final_Localization": "Cytoplasmic", "Score": "0.95"},
            {"SeqID": "seq2", "Final_Localization": "Membrane", "Score": "0.90"}
        ]
    )
    
    # Run the function
    output_csv = temp_dir / "output.csv"
    result = run_psortb(sample_fasta_file, temp_dir, "output.csv", "n")
    
    # Verify results
    assert result == 0
    assert os.path.exists(output_csv)
    
    # Check CSV content
    with open(output_csv, "r") as f:
        lines = f.readlines()
        assert lines[0].strip() == "SeqID,Final_Localization,Score"
        assert lines[1].strip() == "seq1,Cytoplasmic,0.95"
        assert lines[2].strip() == "seq2,Membrane,0.90"

def test_find_psortb_executable_not_found():
    """Test find_psortb_executable when no executable is found."""
    with patch("modules.psortb.os.path.isfile", return_value=False):
        with patch("modules.psortb.shutil.which", return_value=None):
            result = find_psortb_executable()
            assert result is None

def test_check_blast_not_found():
    """Test check_blast when BLAST+ is not found."""
    with patch("modules.psortb.shutil.which", return_value=None):
        assert check_blast() is False

def test_setup_bioperl():
    """Test setup_bioperl sets PERL5LIB correctly."""
    with patch("modules.psortb.os.path.isdir", return_value=True):
        with patch("modules.psortb.os.environ", {"PERL5LIB": ""}):
            setup_bioperl()
            expected_path = str(ROOT_DIR / "tools/psortb/bioperl-live")
            assert expected_path in os.environ["PERL5LIB"].split(":")

if __name__ == "__main__":
    pytest.main()

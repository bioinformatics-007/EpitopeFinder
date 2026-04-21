
import pytest
import pandas as pd
import subprocess
from unittest.mock import patch, mock_open, MagicMock
from pathlib import Path
import tempfile
import os
import logging
import sys

# Add parent directory to sys.path to resolve imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules.iapred import parse_fasta, run_iapred, BASE_DIR, IAPRED_TOOL_PATH

# Configure logging for tests
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

# Fixture to create a temporary FASTA file
@pytest.fixture
def temp_fasta_file():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        content = ">Sequence1|UniProt1\nATGCTAGCT\n>Sequence2|UniProt2\nGCTAGCTAG\n"
        f.write(content)
        f.flush()
        logger.debug(f"Created temp FASTA file {f.name} with content:\n{content}")
        # Verify file is readable
        with open(f.name, 'r') as f_check:
            read_content = f_check.read()
            logger.debug(f"Read back content from {f.name}:\n{read_content}")
            assert read_content == content, f"File content mismatch in {f.name}"
        yield f.name
    logger.debug(f"Removing temp FASTA file {f.name}")
    os.unlink(f.name)

# Fixture to create a temporary directory
@pytest.fixture
def temp_output_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)

# Fixture to mock IAPRED_TOOL_PATH and models directory
@pytest.fixture
def mock_iapred_path():
    with patch("modules.iapred.IAPRED_TOOL_PATH", Path(tempfile.gettempdir()) / "tools" / "IApred") as mock_path:
        (mock_path / "models").mkdir(parents=True, exist_ok=True)
        (mock_path / "IApred.py").touch()
        yield mock_path

# Test FASTA file content
def test_fasta_file_content(temp_fasta_file):
    with open(temp_fasta_file, 'r') as f:
        content = f.read()
    expected = ">Sequence1|UniProt1\nATGCTAGCT\n>Sequence2|UniProt2\nGCTAGCTAG\n"
    assert content == expected, f"FASTA file content mismatch. Expected:\n{expected}\nGot:\n{content}"

# Test parse_fasta directly
def test_parse_fasta_direct(temp_fasta_file):
    try:
        sequences = list(parse_fasta(temp_fasta_file))
        logger.debug(f"Parsed sequences: {sequences}")
        assert len(sequences) == 2, f"Expected 2 sequences, got {len(sequences)}"
        assert sequences[0] == ("Sequence1|UniProt1", "UniProt1", "ATGCTAGCT")
        assert sequences[1] == ("Sequence2|UniProt2", "UniProt2", "GCTAGCTAG")
    except Exception as e:
        logger.error(f"parse_fasta failed: {str(e)}")
        raise

# Test parse_fasta function
def test_parse_fasta_valid(temp_fasta_file):
    sequences = list(parse_fasta(temp_fasta_file))
    assert len(sequences) == 2
    assert sequences[0] == ("Sequence1|UniProt1", "UniProt1", "ATGCTAGCT")
    assert sequences[1] == ("Sequence2|UniProt2", "UniProt2", "GCTAGCTAG")

def test_parse_fasta_empty():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write("")
        f.flush()
        with pytest.raises(ValueError, match="No valid sequences found in FASTA file"):
            list(parse_fasta(f.name))
    os.unlink(f.name)

def test_parse_fasta_invalid():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write("Invalid content")
        f.flush()
        with pytest.raises(ValueError, match="No valid sequences found in FASTA file"):
            list(parse_fasta(f.name))
    os.unlink(f.name)

# Test run_iapred function
@patch("modules.iapred.parse_fasta")
@patch("subprocess.run")
@patch("shutil.copytree")
@patch("os.remove")
@patch("shutil.rmtree")
@patch("pandas.DataFrame.to_csv")
def test_run_iapred_success(mock_to_csv, mock_rmtree, mock_remove, mock_copytree, mock_subprocess, mock_parse_fasta, temp_fasta_file, temp_output_dir, mock_iapred_path):
    # Mock parse_fasta to return expected sequences
    mock_parse_fasta.return_value = [
        ("Sequence1|UniProt1", "UniProt1", "ATGCTAGCT"),
        ("Sequence2|UniProt2", "UniProt2", "GCTAGCTAG")
    ]
    # Mock subprocess.run to simulate successful IAPred execution
    mock_subprocess.return_value = MagicMock(returncode=0, stdout="IAPred completed", stderr="")
    # Create a mock IAPred output CSV
    mock_csv_content = "Header,Start,End,Prediction\nSequence1|UniProt1,1,9,Positive\nSequence2|UniProt2,1,9,Negative\n"
    with patch("builtins.open", mock_open(read_data=mock_csv_content)):
        with patch("pandas.read_csv", return_value=pd.DataFrame({
            "Header": ["Sequence1|UniProt1", "UniProt2"],
            "Start": [1, 1],
            "End": [9, 9],
            "Prediction": ["Positive", "Negative"]
        })):
            with patch("os.path.isfile", return_value=True):  # Mock file existence
                logger.debug(f"Calling run_iapred with FASTA file: {temp_fasta_file}")
                result = run_iapred(temp_fasta_file, temp_output_dir, batch_idx=1)
                assert result == 0
                logger.debug(f"to_csv called with: {mock_to_csv.call_args_list}")
                # Check to_csv calls instead of exists()
                expected_calls = [
                    f"{temp_output_dir}/combined_iapred_batch_1.csv",
                    f"{temp_output_dir}/UniProt1_iapred.csv",
                    f"{temp_output_dir}/UniProt2_iapred.csv"
                ]
                to_csv_calls = [call[0][0] for call in mock_to_csv.call_args_list]
                for expected in expected_calls:
                    assert any(expected in call for call in to_csv_calls), f"Expected to_csv call for {expected}"

def test_run_iapred_missing_iapred_script(temp_fasta_file, temp_output_dir):
    # Mock IAPRED_TOOL_PATH to point to a non-existent IApred.py
    with patch("modules.iapred.IAPRED_TOOL_PATH", Path(tempfile.gettempdir()) / "nonexistent"):
        result = run_iapred(temp_fasta_file, temp_output_dir)
        assert result == 1

def test_run_iapred_missing_models_dir(temp_fasta_file, temp_output_dir, mock_iapred_path):
    # Remove models directory
    (mock_iapred_path / "models").rmdir()
    result = run_iapred(temp_fasta_file, temp_output_dir)
    assert result == 1

@patch("subprocess.run")
def test_run_iapred_subprocess_failure(mock_subprocess, temp_fasta_file, temp_output_dir, mock_iapred_path):
    # Mock subprocess.run to simulate failure
    mock_subprocess.side_effect = subprocess.CalledProcessError(1, cmd=["python3"], stderr="IAPred error")
    result = run_iapred(temp_fasta_file, temp_output_dir)
    assert result == 1

@patch("modules.iapred.parse_fasta")
@patch("subprocess.run")
def test_run_iapred_empty_output(mock_subprocess, mock_parse_fasta, temp_fasta_file, temp_output_dir, mock_iapred_path):
    # Mock parse_fasta to return no sequences
    mock_parse_fasta.return_value = []
    # Mock subprocess.run
    mock_subprocess.return_value = MagicMock(returncode=0, stdout="IAPred completed", stderr="")
    with patch("pandas.read_csv", return_value=pd.DataFrame({
        "Header": pd.Series([], dtype="object"),
        "Start": [],
        "End": [],
        "Prediction": [],
        "UniProt_ID": pd.Series([], dtype="object"),
        "Epitope_Sequence": pd.Series([], dtype="object")
    })):
        with patch("os.path.isfile", return_value=True):
            with patch("pandas.DataFrame.to_csv"):
                logger.debug(f"Calling run_iapred with FASTA file: {temp_fasta_file}")
                result = run_iapred(temp_fasta_file, temp_output_dir)
                # Expect non-zero return code for empty sequences
                assert result == 1

@patch("modules.iapred.parse_fasta")
@patch("subprocess.run")
@patch("shutil.copytree")
@patch("os.remove")
@patch("shutil.rmtree")
@patch("pandas.DataFrame.to_csv")
def test_run_iapred_specific_output_file(mock_to_csv, mock_rmtree, mock_remove, mock_copytree, mock_subprocess, mock_parse_fasta, temp_fasta_file, temp_output_dir, mock_iapred_path):
    # Mock parse_fasta to return expected sequences
    mock_parse_fasta.return_value = [
        ("Sequence1|UniProt1", "UniProt1", "ATGCTAGCT")
    ]
    # Mock subprocess.run to simulate successful IAPred execution
    mock_subprocess.return_value = MagicMock(returncode=0, stdout="IAPred completed", stderr="")
    # Create a mock IAPred output CSV
    mock_csv_content = "Header,Start,End,Prediction\nSequence1|UniProt1,1,9,Positive\n"
    output_file = temp_output_dir / "custom_output.csv"
    with patch("builtins.open", mock_open(read_data=mock_csv_content)):
        with patch("pandas.read_csv", return_value=pd.DataFrame({
            "Header": ["Sequence1|UniProt1"],
            "Start": [1],
            "End": [9],
            "Prediction": ["Positive"]
        })):
            with patch("os.path.isfile", return_value=True):  # Mock file existence
                logger.debug(f"Calling run_iapred with FASTA file: {temp_fasta_file}")
                result = run_iapred(temp_fasta_file, temp_output_dir, batch_idx=1, output_file=str(output_file), uniprot_id="UniProt1")
                assert result == 0
                logger.debug(f"to_csv called with: {mock_to_csv.call_args_list}")
                # Check to_csv call for output_file
                to_csv_calls = [call[0][0] for call in mock_to_csv.call_args_list]
                assert any(str(output_file) in call for call in to_csv_calls), f"Expected to_csv call for {output_file}"


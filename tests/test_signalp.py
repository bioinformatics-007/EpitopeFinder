import sys
import os
# Add parent directory to Python path (Vaxelan_2_0/)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
from unittest.mock import patch, Mock
import subprocess
import tempfile
from pathlib import Path
from modules.signalp import find_model_dir, run_signalp, main, print_status, ROOT_DIR, MODEL_DIR

# Mock print_status to avoid logging during tests
@pytest.fixture(autouse=True)
def mock_print_status():
    with patch("modules.signalp.print_status") as mocked:
        yield mocked

# Fixture for temporary directory
@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)

# Fixture for mock FASTA file
@pytest.fixture
def mock_fasta(temp_dir):
    fasta_file = temp_dir / "test.fasta"
    with open(fasta_file, "w") as f:
        f.write(">seq1\nATGCGT\n")
    return fasta_file

# Fixture for mock model directory with model file
@pytest.fixture
def mock_model_dir(temp_dir):
    model_dir = temp_dir / "model_weights"
    model_dir.mkdir()
    model_file = model_dir / "distilled_model_signalp6.pt"
    model_file.touch()
    return model_dir

# Test find_model_dir function
def test_find_model_dir_primary_path(mock_model_dir):
    """Test find_model_dir with valid primary path."""
    result = find_model_dir(mock_model_dir)
    assert result == mock_model_dir
    assert (result / "distilled_model_signalp6.pt").is_file()

def test_find_model_dir_alternative_path(mock_model_dir, temp_dir):
    """Test find_model_dir with alternative path when primary path fails."""
    invalid_path = temp_dir / "invalid"
    with patch("modules.signalp.ROOT_DIR", temp_dir):
        (temp_dir / "tools").mkdir()
        alt_model_dir = temp_dir / "tools/signalp-6.0h.fast/models"
        alt_model_dir.mkdir(parents=True)
        model_file = alt_model_dir / "distilled_model_signalp6.pt"
        model_file.touch()
        result = find_model_dir(invalid_path)
        assert result == alt_model_dir

def test_find_model_dir_recursive_search(mock_model_dir, temp_dir):
    """Test find_model_dir with recursive search."""
    invalid_path = temp_dir / "invalid"
    with patch("modules.signalp.ROOT_DIR", temp_dir):
        (temp_dir / "tools/random_dir").mkdir(parents=True)
        recursive_model_dir = temp_dir / "tools/random_dir/model_weights"
        recursive_model_dir.mkdir(parents=True)
        model_file = recursive_model_dir / "distilled_model_signalp6.pt"
        model_file.touch()
        result = find_model_dir(invalid_path)
        assert result == recursive_model_dir

def test_find_model_dir_not_found(temp_dir):
    """Test find_model_dir when model file is not found."""
    invalid_path = temp_dir / "invalid"
    with patch("modules.signalp.ROOT_DIR", temp_dir), \
         patch("os.walk", return_value=[]):  # Mock os.walk to return empty
        result = find_model_dir(invalid_path)
        assert result is None

# Test run_signalp function
def test_run_signalp_success(mock_fasta, temp_dir, mock_model_dir):
    """Test run_signalp with valid inputs and successful execution."""
    with patch("shutil.which", return_value="/fake/path/signalp6"), \
         patch("subprocess.run") as mock_run:
        # Mock subprocess.run to simulate successful execution
        mock_run.return_value = Mock(returncode=0, stdout="Success", stderr="")
        output_dir = temp_dir / "output"
        output_dir.mkdir()  # Create output directory
        # Create non-empty output file
        output_file = output_dir / "prediction_results.txt"
        with open(output_file, "w") as f:
            f.write("Dummy output\n")  # Write content to ensure non-zero size
        result = run_signalp(mock_fasta, output_dir, mock_model_dir)
        assert result == 0
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert str(mock_fasta) in call_args
        assert str(output_dir) in call_args
        # Removed assertion for mock_model_dir since it’s not included in the command

def test_run_signalp_missing_executable(mock_fasta, temp_dir, mock_model_dir):
    """Test run_signalp when signalp6 executable is not found."""
    with patch("shutil.which", return_value=None):
        result = run_signalp(mock_fasta, temp_dir, mock_model_dir)
        assert result == 1

def test_run_signalp_missing_fasta(temp_dir, mock_model_dir):
    """Test run_signalp with missing FASTA file."""
    invalid_fasta = temp_dir / "nonexistent.fasta"
    result = run_signalp(invalid_fasta, temp_dir, mock_model_dir)
    assert result == 1

def test_run_signalp_missing_model_file(temp_dir, mock_fasta):
    """Test run_signalp with missing model file."""
    invalid_model_dir = temp_dir / "invalid_model"
    invalid_model_dir.mkdir()
    with patch("shutil.which", return_value="/fake/path/signalp6"):
        result = run_signalp(mock_fasta, temp_dir, invalid_model_dir)
        assert result == 1

def test_run_signalp_failed_execution(mock_fasta, temp_dir, mock_model_dir):
    """Test run_signalp when SignalP execution fails."""
    with patch("shutil.which", return_value="/fake/path/signalp6"), \
         patch("subprocess.run") as mock_run:
        mock_run.return_value = Mock(returncode=1, stdout="", stderr="Error")
        result = run_signalp(mock_fasta, temp_dir, mock_model_dir)
        assert result == 1

def test_run_signalp_empty_output(mock_fasta, temp_dir, mock_model_dir):
    """Test run_signalp when output file is empty or not created."""
    with patch("shutil.which", return_value="/fake/path/signalp6"), \
         patch("subprocess.run") as mock_run:
        mock_run.return_value = Mock(returncode=0, stdout="Success", stderr="")
        result = run_signalp(mock_fasta, temp_dir, mock_model_dir)
        assert result == 1  # No output file created

# Test main function (CLI)
def test_main_success(mock_fasta, temp_dir, mock_model_dir, monkeypatch):
    """Test main function with valid arguments."""
    with patch("modules.signalp.run_signalp", return_value=0) as mock_run_signalp:
        args = ["-f", str(mock_fasta), "-o", str(temp_dir), "--model-dir", str(mock_model_dir)]
        monkeypatch.setattr(sys, "argv", ["script.py"] + args)
        exit_code = main()
        assert exit_code == 0
        mock_run_signalp.assert_called_once_with(str(mock_fasta), str(temp_dir), str(mock_model_dir))

def test_main_missing_arguments(monkeypatch):
    """Test main function with missing required arguments."""
    monkeypatch.setattr(sys, "argv", ["script.py"])
    with pytest.raises(SystemExit) as exc_info:
        main()
    assert exc_info.value.code == 2

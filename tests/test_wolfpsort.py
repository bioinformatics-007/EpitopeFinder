import pytest
import os
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock
import tempfile
import subprocess
import sys
from io import StringIO
import shutil

# Adjust sys.path to include the directory containing wolfpsort.py
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'modules')))
try:
    from wolfpsort import print_status, validate_file, validate_directory, run_wolf_psort, ROOT_DIR
except ModuleNotFoundError as e:
    raise ModuleNotFoundError(f"Could not find 'wolfpsort' module. Ensure 'wolfpsort.py' is in ~/Downloads/Vaxelan_2_0/modules. Error: {e}")

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
def capsys(capsys):
    """Fixture to capture stdout/stderr."""
    return capsys

@pytest.fixture
def mock_wolf_psort_dir(tmp_path):
    """Create a mock WoLF PSORT directory structure."""
    wolf_dir = tmp_path / "WoLFPSort-master"
    bin_dir = wolf_dir / "bin"
    bin_dir.mkdir(parents=True)
    # Create mock scripts
    (bin_dir / "runWolfPsortSummary").write_text("")
    (bin_dir / "runWolfPsortHtmlTables").write_text("")
    return wolf_dir

def test_print_status(capsys):
    """Test print_status function with different status types."""
    print_status("Test info message", "info")
    captured = capsys.readouterr()
    assert "Test info message" in captured.out
    assert "\033[94m" in captured.out  # Blue color for info

    print_status("Test success message", "success")
    captured = capsys.readouterr()
    assert "Test success message" in captured.out
    assert "\033[92m" in captured.out  # Green color for success

    print_status("Test warning message", "warning")
    captured = capsys.readouterr()
    assert "Test warning message" in captured.out
    assert "\033[93m" in captured.out  # Yellow color for warning

    print_status("Test error message", "error")
    captured = capsys.readouterr()
    assert "Test error message" in captured.out
    assert "\033[91m" in captured.out  # Red color for error

def test_validate_file(temp_fasta_file, tmp_path):
    """Test validate_file function."""
    # Valid file
    assert validate_file(temp_fasta_file, "Test FASTA") == Path(temp_fasta_file)

    # Non-existent file
    non_existent = tmp_path / "non_existent.fasta"
    with pytest.raises(FileNotFoundError, match="Test file not found or not readable"):
        validate_file(non_existent, "Test file")

    # Non-readable file
    non_readable = tmp_path / "non_readable.fasta"
    non_readable.write_text("test")
    os.chmod(non_readable, 0o000)  # Remove read permissions
    with pytest.raises(FileNotFoundError, match="Test file not found or not readable"):
        validate_file(non_readable, "Test file")
    os.chmod(non_readable, 0o600)  # Restore permissions for cleanup

def test_validate_directory(tmp_path):
    """Test validate_directory function."""
    # Existing directory
    existing_dir = tmp_path / "existing"
    existing_dir.mkdir()
    assert validate_directory(existing_dir, "Test dir") == existing_dir.resolve()

    # Non-existent directory (should be created)
    new_dir = tmp_path / "new_dir"
    assert validate_directory(new_dir, "Test dir") == new_dir.resolve()
    assert new_dir.exists()

    # Directory creation failure (e.g., no permissions)
    with patch("pathlib.Path.mkdir", side_effect=OSError("Permission denied")):
        with pytest.raises(OSError, match="Cannot create Test dir"):
            validate_directory(tmp_path / "no_perm", "Test dir")

@patch("subprocess.run")
@patch("os.chdir")
def test_run_wolf_psort_success(mock_chdir, mock_subprocess_run, temp_fasta_file, temp_dir, mock_wolf_psort_dir, capsys):
    """Test run_wolf_psort with successful execution."""
    mock_subprocess_run.return_value = MagicMock(
        returncode=0,
        stdout="P12345\tnucleus\t0.8\nQ67890\tcytoplasm\t0.6\n",
        stderr=""
    )

    output_file = temp_dir / "output.txt"
    exit_code = run_wolf_psort(
        input_fasta=temp_fasta_file,
        output_file=output_file,
        wolf_psort_path=mock_wolf_psort_dir,
        organism_type="fungi",
        html_output=False
    )

    assert exit_code == 0
    captured = capsys.readouterr()
    assert "Running WoLF PSORT analysis on" in captured.out
    assert "WoLF PSORT analysis completed" in captured.out
    assert "Created individual UniProt ID files and combined_results.txt" in captured.out
    assert (temp_dir / "combined_results.txt").exists()
    assert (temp_dir / "P12345.txt").exists()
    assert (temp_dir / "Q67890.txt").exists()
    with open(temp_dir / "combined_results.txt", "r") as f:
        content = f.read()
        assert "P12345\tnucleus\t0.8" in content
        assert "Q67890\tcytoplasm\t0.6" in content

@patch("subprocess.run")
@patch("os.chdir")
def test_run_wolf_psort_html_output(mock_chdir, mock_subprocess_run, temp_fasta_file, temp_dir, mock_wolf_psort_dir, capsys):
    """Test run_wolf_psort with HTML output mode."""
    mock_subprocess_run.return_value = MagicMock(
        returncode=0,
        stdout="P12345\tnucleus\t0.8\nQ67890\tcytoplasm\t0.6\n",
        stderr=""
    )

    output_file = temp_dir / "output.html"
    exit_code = run_wolf_psort(
        input_fasta=temp_fasta_file,
        output_file=output_file,
        wolf_psort_path=mock_wolf_psort_dir,
        organism_type="plant",
        html_output=True
    )

    assert exit_code == 0
    captured = capsys.readouterr()
    assert "Running WoLF PSORT analysis on" in captured.out
    assert "WoLF PSORT analysis completed" in captured.out
    assert "Created individual UniProt ID files and combined_results.txt" in captured.out
    assert (temp_dir / "combined_results.txt").exists()
    assert (temp_dir / "P12345.txt").exists()
    assert (temp_dir / "Q67890.txt").exists()

@patch("subprocess.run")
@patch("os.chdir")
def test_run_wolf_psort_invalid_organism(mock_chdir, mock_subprocess_run, temp_fasta_file, temp_dir, mock_wolf_psort_dir):
    """Test run_wolf_psort with invalid organism type."""
    with pytest.raises(ValueError, match="Organism type must be one of"):
        run_wolf_psort(
            input_fasta=temp_fasta_file,
            output_file=temp_dir / "output.txt",
            wolf_psort_path=mock_wolf_psort_dir,
            organism_type="invalid"
        )
    mock_subprocess_run.assert_not_called()

@patch("subprocess.run")
@patch("os.chdir")
def test_run_wolf_psort_missing_script(mock_chdir, mock_subprocess_run, temp_fasta_file, temp_dir):
    """Test run_wolf_psort with missing WoLF PSORT script."""
    invalid_wolf_dir = temp_dir / "invalid_wolf"
    bin_dir = invalid_wolf_dir / "bin"
    bin_dir.mkdir(parents=True)  # Create bin directory but no scripts
    with pytest.raises(FileNotFoundError, match="WoLF PSORT script"):
        run_wolf_psort(
            input_fasta=temp_fasta_file,
            output_file=temp_dir / "output.txt",
            wolf_psort_path=invalid_wolf_dir
        )
    mock_subprocess_run.assert_not_called()

@patch("subprocess.run")
@patch("os.chdir")
def test_run_wolf_psort_subprocess_failure(mock_chdir, mock_subprocess_run, temp_fasta_file, temp_dir, mock_wolf_psort_dir, capsys):
    """Test run_wolf_psort with subprocess failure."""
    mock_subprocess_run.side_effect = subprocess.CalledProcessError(
        1, ["./runWolfPsortSummary", "fungi"], stderr="Subprocess failed"
    )

    exit_code = run_wolf_psort(
        input_fasta=temp_fasta_file,
        output_file=temp_dir / "output.txt",
        wolf_psort_path=mock_wolf_psort_dir,
        organism_type="fungi"
    )

    assert exit_code == 1
    captured = capsys.readouterr()
    assert "Error during WoLF PSORT analysis" in captured.out

@patch("subprocess.run")
@patch("os.chdir")
def test_run_wolf_psort_invalid_fasta(mock_chdir, mock_subprocess_run, temp_dir, mock_wolf_psort_dir, capsys):
    """Test run_wolf_psort with invalid FASTA file."""
    invalid_fasta = temp_dir / "invalid.fasta"
    invalid_fasta.write_text("not a fasta file")

    mock_subprocess_run.side_effect = subprocess.CalledProcessError(
        1, ["./runWolfPsortSummary", "fungi"], stderr="Invalid FASTA format"
    )

    exit_code = run_wolf_psort(
        input_fasta=invalid_fasta,
        output_file=temp_dir / "output.txt",
        wolf_psort_path=mock_wolf_psort_dir
    )

    assert exit_code == 1
    captured = capsys.readouterr()
    assert "Error during WoLF PSORT analysis" in captured.out

if __name__ == "__main__":
    pytest.main([__file__])

import pytest
import os
import subprocess
from pathlib import Path
from unittest.mock import patch, mock_open, MagicMock

# Add the modules directory to sys.path using a relative path
import sys
sys.path.append(str(Path(__file__).resolve().parent.parent / 'modules'))
from netchop import run_netchop, run_all_methods, print_status, ROOT_DIR, NETCHOP_SCRIPT_PATH

# Sample amino acid FASTA content for testing
SAMPLE_FASTA = """>protein_1
ACDEFGHIKLMNPQRSTVWY
>protein_2
GAVLIPFYWSTCMNQDEHKR
"""


@pytest.fixture
def temp_dir(tmp_path):
    """Fixture to create a temporary directory."""
    return tmp_path

@pytest.fixture
def sample_fasta_file(temp_dir):
    """Fixture to create a sample amino acid FASTA file."""
    fasta_path = temp_dir / "input.fasta"
    with open(fasta_path, "w") as f:
        f.write(SAMPLE_FASTA)
    return str(fasta_path)

@pytest.fixture
def mock_subprocess():
    """Fixture to mock subprocess.run."""
    with patch("subprocess.run") as mock_run:
        yield mock_run

def test_print_status(capsys):
    """Test the print_status function."""
    print_status("Test message", "info")
    captured = capsys.readouterr()
    assert "Test message" in captured.out
    assert "\033[94m" in captured.out  # Check for blue color code


def test_run_netchop_success(sample_fasta_file, temp_dir, mock_subprocess):
    """Test run_netchop with successful execution on amino acid FASTA."""
    mock_stdout = """#amino_acid\tprediction_score
A\t0.95
C\t0.90
"""
    mock_subprocess.return_value = MagicMock(returncode=0, stdout=mock_stdout, stderr="")
    
    output_dir = str(temp_dir)
    with patch("netchop.plot_results"):  # Mock plot_results to avoid matplotlib rendering in tests
        result = run_netchop(sample_fasta_file, output_dir=output_dir, method="netchop")
    
    assert result == 0
    assert os.path.exists(os.path.join(output_dir, "netchop_out.csv"))
    
    with open(os.path.join(output_dir, "netchop_out.csv"), "r") as f:
        content = f.read()
        assert "amino_acid\tprediction_score" in content
        assert "A\t0.95" in content

def test_run_netchop_input_not_found(temp_dir):
    """Test run_netchop with non-existent input file."""
    with pytest.raises(FileNotFoundError):
        run_netchop(str(temp_dir / "nonexistent.fasta"), output_dir=str(temp_dir), method="netchop")

def test_run_netchop_script_not_found(sample_fasta_file, temp_dir):
    """Test run_netchop when NetChop script is not found."""
    with patch("os.path.isfile", side_effect=lambda x: False if x == NETCHOP_SCRIPT_PATH else True):
        with pytest.raises(FileNotFoundError):
            run_netchop(sample_fasta_file, output_dir=str(temp_dir), method="netchop")

def test_run_netchop_subprocess_error(sample_fasta_file, temp_dir, mock_subprocess):
    """Test run_netchop when subprocess fails."""
    mock_subprocess.side_effect = subprocess.CalledProcessError(1, cmd=["python"], stderr="Subprocess error")
    
    with pytest.raises(RuntimeError, match="netchop execution failed: Subprocess error"):
        run_netchop(sample_fasta_file, output_dir=str(temp_dir), method="netchop")

@patch("netchop.plot_results")
def test_run_netchop_with_threshold(mock_plot_results, sample_fasta_file, temp_dir, mock_subprocess):
    """Test run_netchop with a threshold value."""
    mock_subprocess.return_value = MagicMock(returncode=0, stdout="Peptide,Score\nACDEFGHIK,0.95", stderr="")
    
    run_netchop(sample_fasta_file, output_dir=str(temp_dir), method="netchop", threshold=0.5)
    
    assert mock_subprocess.called
    args = mock_subprocess.call_args[0][0]
    assert "--threshold" in args
    assert "0.5" in args

def test_run_netchop_custom_output_file(sample_fasta_file, temp_dir, mock_subprocess):
    """Test run_netchop with a custom output file name."""
    mock_stdout = """#amino_acid\tprediction_score
A\t0.95
C\t0.90
"""
    mock_subprocess.return_value = MagicMock(returncode=0, stdout=mock_stdout, stderr="")
    
    output_dir = str(temp_dir)
    custom_output = "custom_netchop"
    with patch("netchop.plot_results"):
        run_netchop(sample_fasta_file, output_dir=output_dir, method="netchop", output_file=custom_output)
    
    assert os.path.exists(os.path.join(output_dir, "custom_netchop_netchop.csv"))

@patch("netchop.plot_results")
def test_run_all_methods_success(mock_plot_results, sample_fasta_file, temp_dir, mock_subprocess):
    """Test run_all_methods with successful execution for all methods."""
    mock_subprocess.return_value = MagicMock(returncode=0, stdout="Peptide,Score\nACDEFGHIK,0.95", stderr="")
    
    results = run_all_methods(sample_fasta_file, output_dir=str(temp_dir))
    
    assert results == {"netchop": 0, "netctlpan": 0}
    assert os.path.exists(os.path.join(temp_dir, "netchop_out.csv"))
    assert os.path.exists(os.path.join(temp_dir, "netctlpan_out.csv"))

@patch("netchop.plot_results")
def test_run_all_methods_partial_failure(mock_plot_results, sample_fasta_file, temp_dir, mock_subprocess):
    """Test run_all_methods when one method fails."""
    mock_subprocess.side_effect = [
        MagicMock(returncode=0, stdout="Peptide,Score\nACDEFGHIK,0.95", stderr=""),  # netchop succeeds
        subprocess.CalledProcessError(1, cmd=["python"], stderr="netctlpan error")  # netctlpan fails
    ]
    
    results = run_all_methods(sample_fasta_file, output_dir=str(temp_dir))
    
    assert results == {"netchop": 0, "netctlpan": -1}
    assert os.path.exists(os.path.join(temp_dir, "netchop_out.csv"))
    assert not os.path.exists(os.path.join(temp_dir, "netctlpan_out.csv"))

import pytest
from pathlib import Path
import sys
from unittest.mock import patch, MagicMock

# Add modules directory to sys.path for import
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "modules"))

from clbtope import run_clbtope

@pytest.fixture
def sample_fasta(tmp_path):
    """
    Create a sample FASTA file in a temporary directory for testing.
    """
    fasta_file = tmp_path / "sample.fasta"
    fasta_file.write_text(
        ">Test_protein\n"
        "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYTPGAGDP\n"
    )
    return fasta_file

@pytest.fixture
def output_file(tmp_path):
    """
    Define output file path in a temporary directory.
    """
    return tmp_path / "clbtope_results.csv"

@pytest.fixture
def mock_dependencies(output_file):
    """
    Mock subprocess.run and simulate file creation.
    """
    with patch("clbtope.subprocess.run") as mock_run:
        # Mock subprocess behavior
        mock_process = MagicMock()
        mock_process.returncode = 0
        mock_process.stdout = "Mocked ClbTope run successful."
        mock_process.stderr = ""
        mock_run.return_value = mock_process

        # Simulate output file creation
        output_file.write_text("Mock output content\n")
        yield mock_run

def test_run_clbtope(sample_fasta, output_file, mock_dependencies):
    """
    Test the run_clbtope function with mocked subprocess calls.
    """
    # Run the function
    ret_code = run_clbtope(
        input_file=str(sample_fasta),
        output_file=str(output_file),
    )

    # Assertions
    assert ret_code == 0, "run_clbtope should return 0 on success"
    assert output_file.is_file(), "Output file should be created"
    assert output_file.stat().st_size > 0, "Output file should not be empty"

    # Verify that the subprocess was called correctly
    mock_dependencies.assert_called_once()
    command_args = mock_dependencies.call_args[0][0]
    assert "clbtope.py" in " ".join(command_args), "The ClbTope script should be called"


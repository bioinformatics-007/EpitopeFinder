import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
import os
import sys

# Adjust import path to find your module correctly
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "modules"))

from deeptmhmm import run_deeptmhmm

@pytest.fixture
def sample_fasta(tmp_path):
    fasta_file = tmp_path / "sample.fasta"
    fasta_file.write_text(">Test_protein\nMKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYTPGAGDP\n")
    return fasta_file

@pytest.fixture
def output_dir(tmp_path):
    return tmp_path / "output"

@pytest.fixture
def biolib_path(tmp_path):
    biolib_exe = tmp_path / "biolib"
    biolib_exe.write_text("#!/bin/bash\necho 'Fake biolib'\n")
    biolib_exe.chmod(0o755)
    return biolib_exe

@pytest.fixture
def setup_biolib_results(tmp_path):
    # Create biolib_results in tmp_path to match run_deeptmhmm's expectation
    results_dir = tmp_path / "biolib_results"
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "dummy_output.txt").write_text("dummy content")
    return results_dir

@pytest.fixture
def mock_subprocess_run():
    with patch("deeptmhmm.subprocess.run") as mock_run:
        mock_process = MagicMock()
        mock_process.returncode = 0
        mock_process.stdout = "Mocked DeepTMHMM run successful."
        mock_process.stderr = ""
        mock_run.return_value = mock_process
        yield mock_run

@pytest.fixture
def mock_shutil_move(tmp_path, output_dir):
    with patch("deeptmhmm.shutil.move") as mock_move:
        def fake_move(src, dst):
            # Simulate moving biolib_results to deeptmhmm_results
            dest_dir = output_dir / "deeptmhmm_results"
            dest_dir.mkdir(parents=True, exist_ok=True)
            (dest_dir / "dummy_output.txt").write_text("dummy content")
        mock_move.side_effect = fake_move
        yield mock_move

@pytest.fixture
def set_cwd(tmp_path):
    # Set current working directory to tmp_path for the test
    original_cwd = os.getcwd()
    os.chdir(tmp_path)
    yield
    os.chdir(original_cwd)

def test_run_deeptmhmm(
    tmp_path,  # Add tmp_path to the parameter list
    sample_fasta,
    output_dir,
    biolib_path,
    setup_biolib_results,
    mock_subprocess_run,
    mock_shutil_move,
    set_cwd  # Ensure cwd is tmp_path
):
    # Call the function to test
    ret_code = run_deeptmhmm(
        input_fasta=str(sample_fasta),
        output_dir=str(output_dir),
        biolib_path=str(biolib_path)
    )

    # Assertions
    assert ret_code == 0, "run_deeptmhmm should return 0 on success"
    
    results_dir = output_dir / "deeptmhmm_results"
    assert results_dir.is_dir(), f"{results_dir} directory should be created"
    assert (results_dir / "dummy_output.txt").is_file(), "Output file should be moved"

    mock_subprocess_run.assert_called_once()
    mock_shutil_move.assert_called_once_with(
        str(tmp_path / "biolib_results"),
        str(output_dir / "deeptmhmm_results")
    )

    # Check the subprocess command
    called_args = mock_subprocess_run.call_args[0][0]
    assert "DTU/DeepTMHMM" in called_args

import sys
import os

# Add the parent directory of 'tests/' to the sys.path so 'modules' can be found
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules.netctl import run_netctl, NETCTL_SCRIPT_PATH
import pytest
from unittest.mock import patch
import subprocess

TEST_INPUT_FASTA = "test_input.fasta"
TEST_OUTPUT_DIR = "test_output"
TEST_OUTPUT_FILE = "netctl_test_out.txt"

def test_run_netctl_error_in_subprocess():
    with patch("os.path.isfile", return_value=True), \
         patch("subprocess.run") as mock_run:
        mock_run.side_effect = subprocess.CalledProcessError(
            returncode=1,
            cmd="mock command",
            stderr="Error occurred"
        )
        with pytest.raises(RuntimeError) as excinfo:
            run_netctl(
                input_fasta=TEST_INPUT_FASTA,
                output_dir=TEST_OUTPUT_DIR,
                output_file=TEST_OUTPUT_FILE,
                netctl_script_path=NETCTL_SCRIPT_PATH,
                supertype="A2"
            )
        assert "NetCTL execution failed" in str(excinfo.value)


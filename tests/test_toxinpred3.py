import pytest
import subprocess
from unittest.mock import patch, mock_open, MagicMock
from modules.toxin_epitope import run_toxinpred3

# Paths used in the tests
VALID_MODEL_PATH = "/home/hp-lapi/Downloads/Vaxelan_2_0/tools/toxinpred3/toxinpred3.0_model.pkl"  # Target model path
TOXINPRED_SCRIPT_PATH = "/home/hp-lapi/Downloads/Vaxelan_2_0/tools/toxinpred3/toxinpred3.py"  # Script path
OUTPUT_FILE_PATH = "/home/hp-lapi/Downloads/Vaxelan_2_0/Results/output.csv"
INPUT_FILE_PATH = "/home/hp-lapi/Downloads/Vaxelan_2_0/modules/test_fasta.fasta"
ZIPPED_MODEL_PATH = "/home/hp-lapi/Downloads/Vaxelan_2_0/tools/toxinpred3/model/toxinpred3.0_model.pkl.zip"
UNZIPPED_MODEL_PATH = "/home/hp-lapi/Downloads/Vaxelan_2_0/tools/toxinpred3/model/toxinpred3.0_model.pkl"  # Source path after unzip

# Example valid FASTA content for mocking file reading
VALID_FASTA_CONTENT = (
    ">sequence1\n"
    "GLSQGNENV\n"
    ">sequence2\n"
    "HAYQGALAL\n"
)

@pytest.fixture
def mock_file_open():
    with patch("builtins.open", mock_open(read_data=VALID_FASTA_CONTENT)) as mock_file:
        yield mock_file

@pytest.fixture
def mock_subprocess_run():
    with patch("subprocess.run") as mock_run:
        yield mock_run

def test_run_toxinpred3_success(mock_file_open, mock_subprocess_run):
    with patch("os.path.isfile") as mock_isfile, \
         patch("os.path.exists") as mock_exists, \
         patch("os.makedirs") as mock_makedirs, \
         patch("os.symlink") as mock_symlink, \
         patch("zipfile.ZipFile") as mock_zipfile, \
         patch("shutil.copy") as mock_copy, \
         patch("os.path.getsize") as mock_getsize:

        # Mock file existence checks
        def isfile_side_effect(path):
            return path in (
                VALID_MODEL_PATH,  # Target path after copy
                TOXINPRED_SCRIPT_PATH,
                OUTPUT_FILE_PATH,  # Ensure output file is mocked as existing
                INPUT_FILE_PATH,
                ZIPPED_MODEL_PATH,
                UNZIPPED_MODEL_PATH
            )

        mock_isfile.side_effect = isfile_side_effect
        mock_exists.return_value = True
        mock_zipfile.return_value.__enter__.return_value.extractall.return_value = None
        mock_copy.return_value = None
        mock_getsize.return_value = 1024  # Mock non-empty file size for output file

        # Mock subprocess.run for successful run
        mock_subprocess_run.return_value = MagicMock(stdout="Success output", stderr="", returncode=0)

        # Run function
        run_toxinpred3(INPUT_FILE_PATH, OUTPUT_FILE_PATH)

        # Check subprocess was called correctly
        expected_cmd = [
            "python3",
            TOXINPRED_SCRIPT_PATH,
            "-i", INPUT_FILE_PATH,
            "-o", OUTPUT_FILE_PATH,
            "-t", "0.38",
            "-m", "2",
            "-d", "2"
        ]
        mock_subprocess_run.assert_called_once()
        args_passed = mock_subprocess_run.call_args[0][0]
        assert args_passed == expected_cmd, f"Expected command {expected_cmd}, got {args_passed}"

def test_run_toxinpred3_input_file_not_found(mock_file_open, mock_subprocess_run):
    with patch("os.path.isfile") as mock_isfile:
        # Only model and script files exist; input file does not
        def isfile_side_effect(path):
            return path in (VALID_MODEL_PATH, ZIPPED_MODEL_PATH, UNZIPPED_MODEL_PATH, TOXINPRED_SCRIPT_PATH)
        mock_isfile.side_effect = isfile_side_effect

        with pytest.raises(SystemExit):
            run_toxinpred3("/nonexistent/input.fasta", OUTPUT_FILE_PATH)

def test_run_toxinpred3_invalid_fasta_format(mock_subprocess_run):
    invalid_fasta_content = "Invalid content missing fasta header"
    with patch("builtins.open", mock_open(read_data=invalid_fasta_content)):
        with pytest.raises(SystemExit):
            run_toxinpred3(INPUT_FILE_PATH, OUTPUT_FILE_PATH)

def test_run_toxinpred3_subprocess_failure(mock_file_open, mock_subprocess_run):
    with patch("os.path.isfile") as mock_isfile, \
         patch("os.path.exists") as mock_exists, \
         patch("zipfile.ZipFile") as mock_zipfile, \
         patch("shutil.copy") as mock_copy:

        def isfile_side_effect(path):
            return path in (
                VALID_MODEL_PATH,
                TOXINPRED_SCRIPT_PATH,
                OUTPUT_FILE_PATH,
                INPUT_FILE_PATH,
                ZIPPED_MODEL_PATH,
                UNZIPPED_MODEL_PATH
            )
        mock_isfile.side_effect = isfile_side_effect
        mock_exists.return_value = True
        mock_zipfile.return_value.__enter__.return_value.extractall.return_value = None
        mock_copy.return_value = None

        # Simulate subprocess failure
        mock_subprocess_run.side_effect = subprocess.CalledProcessError(1, "command", output="", stderr="error")

        # Expect RuntimeError since the function raises it on subprocess failure
        with pytest.raises(RuntimeError):
            run_toxinpred3(INPUT_FILE_PATH, OUTPUT_FILE_PATH)

if __name__ == "__main__":
    pytest.main()

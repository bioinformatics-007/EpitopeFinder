import os
import subprocess
import tempfile
from pathlib import Path
import shutil
import warnings

# Suppress scipy RuntimeWarning if necessary
warnings.filterwarnings("ignore", category=RuntimeWarning, module="scipy._lib.messagestream")

# Define root directory
try:
    ROOT_DIR = Path(__file__).resolve().parent.parent  # Aims for Vaxelan_2_0/
except NameError:
    ROOT_DIR = Path.cwd()  # Fallback to current working directory


# Define default paths with fallbacks
CLBTOPE_SCRIPT = os.environ.get("CLBTOPE_SCRIPT", ROOT_DIR / "tools/clbtope/clbtope/clbtope.py")
BLAST_PATH = os.environ.get("BLAST_PATH", shutil.which("blastp") or "/usr/bin/blastp")
BLAST_DATABASE = os.environ.get("BLAST_DATABASE", ROOT_DIR / "tools/clbtope/clbtope/Database")

def find_clbtope_script(base_path: Path) -> Path:
    """Search for clbtope.py in possible locations, prioritizing the provided path."""
    base_path = base_path.resolve()

    provided_path = ROOT_DIR / "tools/clbtope/clbtope/clbtope.py"
    if provided_path.is_file():
        return provided_path

    possible_paths = [
        base_path / "tools/clbtope/clbtope/clbtope.py",
        ROOT_DIR / "clbtope/clbtope.py",
        ROOT_DIR / "tools/clbtope/clbtope.py",
        base_path / "clbtope.py",
        Path.cwd() / "tools/clbtope/clbtope/clbtope.py",
    ]

    for path in possible_paths:
        if path.is_file():
            return path

    return provided_path

def find_blast_database(base_path: Path) -> Path:
    """Search for BLAST database directory containing train.* files."""
    base_path = base_path.resolve()

    provided_path = ROOT_DIR / "tools/clbtope/clbtope/Database"
    if (provided_path / "train.pdb").is_file():
        return provided_path

    possible_paths = [
        base_path / "tools/clbtope/clbtope/Database",
        ROOT_DIR / "clbtope/Database",
        ROOT_DIR / "tools/clbtope/Database",
        Path.cwd() / "tools/clbtope/clbtope/Database",
        base_path / "Database",
    ]

    for path in possible_paths:
        if (path / "train.pdb").is_file():
            return path
    return provided_path

def run_clbtope(
    input_file: str,
    output_file: str,
    clbtope_script: Path = CLBTOPE_SCRIPT,
    blast_path: str = BLAST_PATH,
    blast_database: Path = BLAST_DATABASE
) -> int:
    """
    Run ClbTope to predict B-cell epitopes using BLASTP.

    Args:
        input_file: Path to the input FASTA file.
        output_file: Path to the output file (e.g., CSV).
        clbtope_script: Path to clbtope.py script.
        blast_path: Path to BLASTP executable.
        blast_database: Path to BLAST database.

    Returns:
        int: 0 if successful, non-zero (1) otherwise.
    """
    try:
        # Resolve paths
        clbtope_script = find_clbtope_script(Path(clbtope_script).parent.parent)
        blast_database = find_blast_database(Path(blast_database).parent.parent)

        # Validate input files
        for path, name in [
            (input_file, "Input FASTA file"),
            (clbtope_script, "ClbTope script"),
            (blast_path, "BLASTP executable")
        ]:
            if not Path(path).is_file():
                raise FileNotFoundError(f"{name} not found: {path}")

        # Verify BLAST database files
        required_db_extensions = [".pdb", ".phr", ".pin", ".pog", ".pos", ".pot", ".psd", ".psi", ".psq", ".ptf", ".pto"]
        blast_db_path = Path(blast_database) / "train"
        missing_files = [str(blast_db_path.with_suffix(ext)) for ext in required_db_extensions if not blast_db_path.with_suffix(ext).is_file()]
        if missing_files:
            raise FileNotFoundError(f"BLAST database files missing: {', '.join(missing_files)}")

        # Handle output path
        output_path = Path(output_file)
        if output_path.is_dir():
            output_file = str(output_path / "clbtope_results.csv")
        else:
            output_path.parent.mkdir(parents=True, exist_ok=True)

        # Create temporary directory for .env file
        with tempfile.TemporaryDirectory() as temp_dir:
            envfile_path = Path(temp_dir) / "clbtope_envfile.env"
            with open(envfile_path, "w") as envfile:
                envfile.write(f"BLAST:{blast_path}\n")
                envfile.write(f"BLAST_DATABASE:{blast_db_path}\n")

            # Construct and run command
            command = [
                "python3", str(clbtope_script),
                "-i", str(input_file),
                "-o", str(output_file),
                "-e", str(envfile_path)
            ]
            

            result = subprocess.run(command, text=True, capture_output=True)
            
            if result.returncode != 0:
               
                return 1

            # Validate output file
            output_path = Path(output_file)
            if output_path.is_file() and output_path.stat().st_size > 0:
                return 0
            else:
                return 1

    except subprocess.CalledProcessError as e:
        
        return 1
    except Exception as e:
       
        return 1

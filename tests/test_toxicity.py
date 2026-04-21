import sys
from pathlib import Path
import os
import pytest
import pandas as pd

# Add project root (parent of tests/) to sys.path so we can import modules.toxicity
sys.path.append(str(Path(__file__).resolve().parent.parent))

from modules.toxicity import (
    find_model_file,
    get_blastp_path,
    reformat_output,
    MERCI_Processor,
    Merci_after_processing,
    BLAST_processor,
    hybrid,
    run_toxinpred,
)

# Define paths for testing relative to project root
ROOT_DIR = Path(__file__).resolve().parent.parent
TEST_FASTA = ROOT_DIR / "tests" / "test_sequences.fasta"
TEST_OUTPUT_DIR = ROOT_DIR / "tests" / "output"
MODEL_FILE = ROOT_DIR / "tools" / "toxinpred2" / "RF_model"

@pytest.fixture(scope="module", autouse=True)
def setup_module():
    """Setup output directory before tests and clean it afterwards."""
    TEST_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Create a simple dummy fasta file if not exists
    if not TEST_FASTA.exists():
        with open(TEST_FASTA, "w") as f:
            f.write(">seq1\nMKTLLILALVGA\n>seq2\nGAVLILMVL\n")
    yield
    # Cleanup test outputs after tests
    for file in TEST_OUTPUT_DIR.glob("*"):
        if file.is_file():
            file.unlink()
    try:
        TEST_OUTPUT_DIR.rmdir()
    except OSError:
        pass  # Directory not empty or other error, ignore

def test_find_model_file():
    model_path = find_model_file(str(MODEL_FILE))
    assert model_path is not None
    assert model_path.exists()

def test_get_blastp_path():
    try:
        blastp_path = get_blastp_path()
        assert blastp_path is not None
        assert os.path.exists(blastp_path), f"BLASTP binary not found at: {blastp_path}"
    except FileNotFoundError as e:
        pytest.skip(f"Skipping test due to missing BLASTP binary: {str(e)}")

def test_reformat_output():
    output_file = TEST_OUTPUT_DIR / "mock_output.csv"
    df = pd.DataFrame({
        'Subject': ['seq1', 'seq2'],
        'ML Score': [0.7, 0.4],
        'Prediction': ['Toxin', 'Non-Toxin']
    })
    df.to_csv(output_file, index=False)

    fasta_file = TEST_FASTA
    status = reformat_output(str(output_file), str(fasta_file), str(TEST_OUTPUT_DIR), display=2)
    assert status == 0
    combined_output = TEST_OUTPUT_DIR / "toxinpred" / "combined.csv"
    assert combined_output.exists()

def test_MERCI_Processor():
    # Create a mock MERCI output text file with tab-separated format
    merci_file = TEST_OUTPUT_DIR / "mock_merci.txt"
    with open(merci_file, 'w') as f:
        f.write("seq1\tCOVERAGE: 100%\nseq2\tCOVERAGE: 0%\n")

    seq_names = ['>seq1', '>seq2']
    merci_processed = TEST_OUTPUT_DIR / "processed_merci.csv"
    
    # Call the MERCI_Processor function
    MERCI_Processor(str(merci_file), str(merci_processed), seq_names)
    
    # Check if the processed file exists and is not empty
    assert merci_processed.exists(), f"Processed MERCI file not created: {merci_processed}"
    assert merci_processed.stat().st_size > 0, "Processed MERCI file is empty"

    # Check DataFrame structure
    df1 = pd.read_csv(merci_processed)
    if df1.empty:
        pytest.fail("Processed DataFrame is empty.")
    
    # Ensure expected columns exist
    assert all(col in df1.columns for col in ['Name', 'Hits', 'Prediction']), \
        f"Expected 'Name', 'Hits', 'Prediction' columns, got {df1.columns}"
    
    # Convert Hits to numeric if they are not already
    df1['Hits'] = df1['Hits'].str.extract('(\d+)').astype(float)  # Extract numeric values
    assert pd.api.types.is_numeric_dtype(df1['Hits']), \
        f"Hits column should be numeric, got values: {df1['Hits'].tolist()}"

def test_Merci_after_processing():
    processed_merci = TEST_OUTPUT_DIR / "processed_merci.csv"
    df = pd.DataFrame({
        'Name': ['seq1', 'seq2'],
        'Hits': [1.0, 0.0],
        'Prediction': ['Toxin', 'Non-Toxin']
    })
    df.to_csv(processed_merci, index=False)

    final_merci = TEST_OUTPUT_DIR / "final_merci.csv"
    Merci_after_processing(str(processed_merci), str(final_merci))
    assert final_merci.exists()
    assert final_merci.stat().st_size > 0

def test_BLAST_processor():
    blast_result = TEST_OUTPUT_DIR / "mock_blast.txt"
    df = pd.DataFrame({
        0: ['seq1', 'seq1', 'seq2'],
        1: ['P_query1_1', 'N_query2_1', 'P_query3_1']
    })
    df.to_csv(blast_result, sep="\t", index=False, header=False)

    seq_names = ['>seq1', '>seq2']
    blast_processed = TEST_OUTPUT_DIR / "processed_blast.csv"
    BLAST_processor(str(blast_result), str(blast_processed), seq_names)
    assert blast_processed.exists()
    assert blast_processed.stat().st_size > 0

def test_hybrid():
    ml_output = TEST_OUTPUT_DIR / "ml_output.csv"
    pd.DataFrame({'ML Score': [0.7, 0.4]}).to_csv(ml_output, index=False)

    seq_names = ['>seq1', '>seq2']
    merci_output = TEST_OUTPUT_DIR / "final_merci.csv"
    pd.DataFrame({'Subject': ['seq1', 'seq2'], 'MERCI Score': [0.5, 0]}).to_csv(merci_output, index=False)

    blast_output = TEST_OUTPUT_DIR / "processed_blast.csv"
    pd.DataFrame({'Subject': ['seq1', 'seq2'], 'BLAST Score': [0.5, 0]}).to_csv(blast_output, index=False)

    final_output = TEST_OUTPUT_DIR / "final_output.csv"
    hybrid(str(ml_output), seq_names, str(merci_output), str(blast_output), threshold=0.6, final_output=str(final_output))
    assert final_output.exists()
    assert final_output.stat().st_size > 0

def test_run_toxinpred():
    status = run_toxinpred(str(TEST_FASTA), str(TEST_OUTPUT_DIR))
    assert status in [0, 1], f"Unexpected status code: {status}"


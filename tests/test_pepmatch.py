import sys
import os
import pytest

# Add project root (one level up) to sys.path so 'modules' can be imported
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules import pepmatch  # Now this import should work

def test_dummy():
    # A simple dummy test to check pytest is working
    assert True

# Example: test validate_fasta function with a minimal valid FASTA string
def test_validate_fasta(tmp_path):
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(">seq1\nACDEFGHIKLMNPQRSTVWY\n")

    assert pepmatch.validate_fasta(str(fasta_file)) is True

def test_validate_fasta_invalid(tmp_path):
    fasta_file = tmp_path / "invalid.fasta"
    fasta_file.write_text("No header line\nACDEFGHIKLMNPQRSTVWY\n")

    assert pepmatch.validate_fasta(str(fasta_file)) is False


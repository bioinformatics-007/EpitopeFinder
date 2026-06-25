import sys
import os
import pytest

# Add project root (one level up) to sys.path so 'modules' can be imported
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules import pepmatch  # Now this import should work

# Test tile_sequence function
def test_tile_sequence(tmp_path):
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(">seq1\nACDEFGHIK\n")
    tiled_fasta = tmp_path / "tiled.fasta"
    
    count = pepmatch.tile_sequence(str(fasta_file), str(tiled_fasta), window_size=5)
    assert count == 5  # ACDEF, CDEFG, DEFGH, EFGHI, FGHIK
    assert os.path.exists(tiled_fasta)

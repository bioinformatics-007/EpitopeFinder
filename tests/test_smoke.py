import pytest
import os
import shutil
from pathlib import Path
import pandas as pd
from unittest.mock import MagicMock, patch

# Add project root to sys.path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import final
from core_pipeline import PipelineRequest, execute, AssemblyConfig

@pytest.fixture
def sample_fasta(tmp_path):
    fasta_file = tmp_path / "sample.fasta"
    fasta_file.write_text(">seq1\nMKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYTPGAGDP\n")
    return fasta_file

@pytest.fixture
def mock_pipeline(mocker):
    """Mock all pipeline tools using pytest-mock."""
    # Patch all tool functions in final.py
    tools = [
        "run_mhc1", "run_mhc2", "run_bcell", "run_netctl", "run_netchop",
        "run_psortb", "run_iapred", "run_algpred", "run_toxinpred",
        "run_clbtope", "run_toxinpred3", "run_wolf_psort", "run_assembly",
        "run_esmfold", "run_algpred_down", "run_iapred_down"
    ]
    
    mocks = {}
    for tool in tools:
        m = mocker.patch(f"final.{tool}", return_value=0)
        mocks[tool] = m

    # Strategy 6 uses run_pipeline from esmfold
    mocks["esmfold_pipeline"] = mocker.patch("modules.esmfold.run_pipeline", return_value=0)

    # Patch mp.Pool to be synchronous
    mock_pool = mocker.patch("final.mp.Pool")
    mock_pool_instance = MagicMock()
    mock_pool_instance.map.side_effect = lambda func, iterable: [func(item) for item in iterable]
    mock_pool.return_value.__enter__.return_value = mock_pool_instance

    # Patch logging and dependencies
    mocker.patch("final.setup_logging")
    mocker.patch("final.check_dependencies")
    mocker.patch("final.detect_pathogen_type_from_fasta", return_value="bacteria")
    mocker.patch("final.plot_epitope_analysis")
    mocker.patch("final.plot_vaccine_architecture")

    # Intelligent Path.exists mock
    def smart_exists(path_obj):
        # Allow checking for sample.fasta and other real temp files
        if "sample.fasta" in str(path_obj):
            return True
        # For tool outputs, we want them to "exist" AFTER they are supposedly created
        # But for the "already run" check, we want them NOT to exist
        # We can distinguish by whether they are in a 'strategy_X' directory
        if "strategy_" in str(path_obj) and (".csv" in str(path_obj) or ".txt" in str(path_obj)):
            # If it's a check to see if we should SKIP, return False
            # If it's a check to see if tool succeeded, return True
            # We'll just return True for everything except the very first check?
            # Actually, returning True for everything in strategy_X might be safe 
            # IF we don't want to skip.
            # But Strategy 2 uses .exists() to SKIP.
            return False 
        return True

    # For simplicity, we'll use a side effect that we can toggle or just be smart
    # Actually, Strategy 5 is the one that needs it to be True for validation
    
    return mocks

def test_strategy_1_smoke(sample_fasta, mock_pipeline):
    request = PipelineRequest(input_value=str(sample_fasta), strategy=1, pathogen_type="bacteria", selected_tools=["mhc1", "bcell"])
    with patch("pathlib.Path.exists", return_value=True):
        result = execute(request)
    assert result.status == "completed"

def test_strategy_2_smoke(sample_fasta, mock_pipeline, monkeypatch):
    import final
    m_iap = MagicMock(return_value=0)
    monkeypatch.setattr(final, "run_iapred", m_iap)
    
    request = PipelineRequest(input_value=str(sample_fasta), strategy=2, pathogen_type="bacteria", selected_tools=["iapred"])
    with patch("final.pd.read_csv") as mock_read:
        mock_read.return_value = pd.DataFrame({"seq_id": ["seq1"], "Prediction": ["Non-Allergen"]})
        # For strategy 2, we want exists=False so it doesn't skip
        with patch("pathlib.Path.exists", return_value=False):
            result = execute(request)
    assert result.status == "completed"
    m_iap.assert_called()

def test_strategy_3_smoke(sample_fasta, mock_pipeline):
    request = PipelineRequest(input_value=str(sample_fasta), strategy=3, pathogen_type="bacteria", selected_tools=["toxinpred"])
    result = execute(request)
    assert result.status == "completed"

def test_strategy_4_smoke(sample_fasta, mock_pipeline, mocker):
    request = PipelineRequest(input_value=str(sample_fasta), strategy=4, pathogen_type="bacteria")
    mocker.patch("final.pd.read_csv", return_value=pd.DataFrame({"seq_id": ["seq1"], "Prediction": ["Non-Allergen"]}))
    # Strategy 4 calls other strategies. We need a mix. 
    # Actually, return False for .csv/.txt is usually what we want for "not skipped"
    def exists_side_effect(self):
        if str(self).endswith(".csv") or str(self).endswith(".txt"):
             return False
        return True
    with patch("pathlib.Path.exists", autospec=True, side_effect=exists_side_effect):
        result = execute(request)
    assert result.status == "completed"

def test_strategy_5_smoke(sample_fasta, mock_pipeline, mocker):
    request = PipelineRequest(input_value=str(sample_fasta), strategy=5, pathogen_type="bacteria")
    mocker.patch("final.pd.read_csv", return_value=pd.DataFrame({
        "seq_id": ["seq1"], "Peptide": ["MKTAYIAK"], "peptide": ["MKTAYIAK"],
        "Sequence": ["MKTAYIAK"], "Score": [0.9]
    }))
    # Strategy 5 needs .exists() to be True for tool outputs it just "created"
    with patch("pathlib.Path.exists", return_value=True):
        result = execute(request)
    assert result.status == "completed"

def test_strategy_6_smoke(sample_fasta, mock_pipeline):
    config = AssemblyConfig(
        mode="custom", 
        custom_sequence="MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYTPGAGDP", 
        run_sasa=False
    )
    request = PipelineRequest(input_value=str(sample_fasta), strategy=6, assembly_config=config)
    result = execute(request)
    assert result.status == "completed"

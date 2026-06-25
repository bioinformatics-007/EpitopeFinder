import pytest
import os
import shutil
import tempfile
import threading
from pathlib import Path
from unittest.mock import patch, MagicMock
import requests

# Add project root to sys.path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core_pipeline import PipelineRequest, execute, PipelineResult
import final

@pytest.fixture
def mock_all_tools(monkeypatch):
    """Bypass any tool execution logic so strategies return immediately."""
    tools = [
        "run_mhc1", "run_mhc2", "run_bcell", "run_netctl", "run_netchop",
        "run_psortb", "run_iapred", "run_algpred", "run_toxinpred",
        "run_clbtope", "run_toxinpred3", "run_wolf_psort", "run_assembly",
        "run_esmfold", "run_algpred_down", "run_iapred_down"
    ]
    for t in tools:
        monkeypatch.setattr(final, t, lambda *args, **kwargs: 0)
    monkeypatch.setattr(final, "check_dependencies", lambda: None)
    monkeypatch.setattr(final, "setup_logging", lambda *args: MagicMock())
    monkeypatch.setattr(final, "plot_epitope_analysis", lambda *args, **kwargs: None)
    monkeypatch.setattr(final, "plot_vaccine_architecture", lambda *args, **kwargs: None)

def test_invalid_fasta_handling(mock_all_tools, tmp_path):
    """Verify that an invalid FASTA file (e.g. junk content) fails validation gracefully."""
    invalid_file = tmp_path / "invalid.fasta"
    invalid_file.write_text("This is not a fasta file, it has no headers.")
    
    request = PipelineRequest(
        input_value=str(invalid_file),
        strategy=1,
        pathogen_type="bacteria"
    )
    result = execute(request)
    assert result.status == "failed"
    assert "Invalid FASTA file" in result.error

def test_empty_fasta_handling(mock_all_tools, tmp_path):
    """Verify that an empty FASTA input fails validation gracefully."""
    empty_file = tmp_path / "empty.fasta"
    empty_file.write_text("")
    
    request = PipelineRequest(
        input_value=str(empty_file),
        strategy=1,
        pathogen_type="bacteria"
    )
    result = execute(request)
    assert result.status == "failed"
    assert "Invalid FASTA file" in result.error or "No valid sequences" in result.error

def test_uniprot_fetch_failure(mock_all_tools, monkeypatch):
    """Verify that UniProt ID retrieval failure is handled gracefully without crashing."""
    def mock_get(*args, **kwargs):
        raise requests.RequestException("UniProt API is down")
        
    monkeypatch.setattr(requests.Session, "get", mock_get)
    
    request = PipelineRequest(
        input_value="INVALID_ID",
        strategy=1,
        pathogen_type="bacteria"
    )
    result = execute(request)
    assert result.status == "failed"
    assert "Failed to fetch UniProt ID" in result.error

def test_concurrency_no_collision(mock_all_tools, tmp_path, monkeypatch):
    """Verify that multiple concurrent executions of the pipeline do not interfere or collide."""
    # We want to run execute in parallel threads
    results = []
    
    fasta_content = ">seq1\nMKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYTPGAGDP\n"
    fasta_file = tmp_path / "input.fasta"
    fasta_file.write_text(fasta_content)
    
    # We want to mock Path.exists to return True so the pipeline succeeds
    monkeypatch.setattr(Path, "exists", lambda self: True)
    
    def run_worker():
        request = PipelineRequest(
            input_value=str(fasta_file),
            strategy=1,
            pathogen_type="bacteria"
        )
        res = execute(request)
        results.append(res)
        
    threads = []
    for _ in range(5):
        t = threading.Thread(target=run_worker)
        threads.append(t)
        t.start()
        
    for t in threads:
        t.join()
        
    # Check that all jobs finished
    assert len(results) == 5
    # Check that each job had a unique results directory to prevent file collisions
    res_dirs = [r.results_dir for r in results]
    assert len(set(res_dirs)) == 5

def test_redis_broker_disconnection_graceful_error(monkeypatch):
    """Verify that a Redis broker connection failure during job submission is handled gracefully returning 503."""
    from fastapi.testclient import TestClient
    from backend.main import app
    from backend.routes.jobs import run_pipeline

    def mock_delay(*args, **kwargs):
        # Celery raises ConnectionError when broker is down
        raise ConnectionError("Could not connect to Redis broker at localhost:6379")

    monkeypatch.setattr(run_pipeline, "delay", mock_delay)

    client = TestClient(app)
    
    response = client.post(
        "/api/jobs/submit",
        json={
            "input_value": "seq1.fasta",
            "strategy": 1,
            "pathogen_type": "bacteria"
        }
    )
    assert response.status_code == 503
    assert "broker is currently unreachable" in response.json()["detail"]

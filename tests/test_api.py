import pytest
import os
import json
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock
from fastapi.testclient import TestClient

# Add project root to sys.path so we can import backend
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from backend.main import app

client = TestClient(app)

def test_health():
    """Test health check endpoint."""
    response = client.get("/api/health")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["version"] == "2.0.0"
    assert "redis" in data["services"]
    assert data["services"]["redis"]["status"] == "ok"

def test_get_strategies():
    """Test strategies configuration endpoint."""
    response = client.get("/api/config/strategies")
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)
    assert len(data) == 6
    assert data[0]["name"] == "Epitope Prediction"

def test_get_mhci_methods():
    """Test MHC-I methods configuration endpoint."""
    response = client.get("/api/config/mhci-methods")
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)
    assert any(m["key"] == "f" for m in data)

def test_get_mhcii_methods():
    """Test MHC-II methods configuration endpoint."""
    response = client.get("/api/config/mhcii-methods")
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)
    assert any(m["key"] == "1" for m in data)

@patch("backend.routes.jobs.run_pipeline.delay")
def test_submit_job(mock_delay, tmp_path):
    """Test job submission endpoint."""
    # Mock JOBS_DIR and _write_meta to avoid real file I/O if possible, 
    # but let's just mock the task dispatch for now.
    payload = {
        "input_value": "P01234",
        "strategy": 1,
        "pathogen_type": "bacteria",
        "mhci_method": "f",
        "mhcii_method": "1",
        "selected_tools": ["mhc1", "bcell"]
    }
    
    # We need to mock _write_meta because it tries to write to the real jobs/ directory
    with patch("backend.routes.jobs._write_meta") as mock_write:
        response = client.post("/api/jobs/submit", json=payload)
        assert response.status_code == 200
        data = response.json()
        assert "job_id" in data
        assert data["status"] == "pending"
        mock_delay.assert_called_once()
        mock_write.assert_called_once()

def test_get_job_status_not_found():
    """Test status polling for non-existent job."""
    with patch("backend.routes.jobs._read_meta", return_value={}):
        response = client.get("/api/jobs/job_unknown/status")
        assert response.status_code == 404
        assert "not found" in response.json()["detail"]

def test_get_job_status_success():
    """Test status polling for an existing job."""
    mock_meta = {
        "job_id": "job_test_123",
        "status": "running",
        "progress_pct": 45.5,
        "current_tool": "MHC-I",
        "failed_tools": [],
        "error": ""
    }
    with patch("backend.routes.jobs._read_meta", return_value=mock_meta):
        response = client.get("/api/jobs/job_test_123/status")
        assert response.status_code == 200
        data = response.json()
        assert data["job_id"] == "job_test_123"
        assert data["status"] == "running"
        assert data["progress_pct"] == 45.5

def test_get_job_results_not_found():
    """Test results retrieval for non-existent job."""
    with patch("backend.routes.jobs._read_meta", return_value={}):
        response = client.get("/api/jobs/job_unknown/results")
        assert response.status_code == 404

def test_get_job_results_pending():
    """Test results retrieval for a job that is not yet completed."""
    mock_meta = {
        "job_id": "job_test_123",
        "status": "running",
        "failed_tools": []
    }
    with patch("backend.routes.jobs._read_meta", return_value=mock_meta):
        response = client.get("/api/jobs/job_test_123/results")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "running"
        assert "outputs" not in data or len(data["outputs"]) == 0

def test_get_job_results_completed(tmp_path):
    """Test results retrieval for a completed job."""
    # Create a dummy result file
    result_file = tmp_path / "summary.csv"
    result_file.write_text("id,score\nseq1,0.9")
    
    mock_meta = {
        "job_id": "job_test_123",
        "status": "completed",
        "results_dir": str(tmp_path),
        "outputs": {
            "summary.csv": str(result_file)
        },
        "failed_tools": []
    }
    with patch("backend.routes.jobs._read_meta", return_value=mock_meta):
        with patch("os.path.exists", return_value=True):
            with patch("os.path.getsize", return_value=100):
                response = client.get("/api/jobs/job_test_123/results")
                assert response.status_code == 200
                data = response.json()
                assert data["status"] == "completed"
                assert len(data["outputs"]) == 1
                assert data["outputs"][0]["relative_path"] == "summary.csv"
                assert "/results/summary.csv" in data["outputs"][0]["download_url"]

@patch("backend.routes.jobs.run_pipeline.delay")
def test_submit_job_with_file(mock_delay, tmp_path):
    """Test submitting a job with a FASTA file upload."""
    fasta_content = ">seq1\nATCG\n"
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)
    
    with patch("backend.routes.jobs._write_meta") as mock_write:
        with open(fasta_file, "rb") as f:
            response = client.post(
                "/api/jobs/submit-with-file",
                files={"file": ("test.fasta", f, "application/octet-stream")},
                data={
                    "strategy": 1,
                    "pathogen_type": "virus",
                    "mhci_method": "f",
                    "mhcii_method": "1",
                    "selected_tools": "mhc1,bcell"
                }
            )
        assert response.status_code == 200
        data = response.json()
        assert "job_id" in data
        assert data["status"] == "pending"
        mock_delay.assert_called_once()
        mock_write.assert_called_once()

def test_download_result_file_success(tmp_path):
    """Test downloading a specific result file."""
    result_file = tmp_path / "summary.csv"
    result_file.write_text("test data")
    
    mock_meta = {
        "job_id": "job_test_123",
        "outputs": {
            "summary.csv": str(result_file)
        }
    }
    with patch("backend.routes.jobs._read_meta", return_value=mock_meta):
        with patch("os.path.exists", side_effect=lambda p: str(p) == str(result_file)):
            response = client.get("/api/jobs/job_test_123/results/summary.csv")
            assert response.status_code == 200
            assert response.content == b"test data"
            assert response.headers["content-type"] == "application/octet-stream"

def test_download_result_file_not_found():
    """Test downloading a non-existent result file."""
    mock_meta = {
        "job_id": "job_test_123",
        "outputs": {}
    }
    with patch("backend.routes.jobs._read_meta", return_value=mock_meta):
        response = client.get("/api/jobs/job_test_123/results/missing.csv")
        assert response.status_code == 404

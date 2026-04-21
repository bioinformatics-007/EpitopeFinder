"""
VaxElan 2.0 — Job management endpoints.

POST /api/jobs/submit      → Submit a new pipeline job
GET  /api/jobs/{id}/status  → Poll job progress
GET  /api/jobs/{id}/results → Retrieve results when completed
GET  /api/jobs/{id}/results/{filename} → Download a specific output file
"""
import os
import uuid
from pathlib import Path

from fastapi import APIRouter, HTTPException, UploadFile, File, Form
from fastapi.responses import FileResponse

from backend.models import (
    JobSubmitRequest,
    JobSubmitResponse,
    JobStatusResponse,
    JobResultsResponse,
    OutputFile,
)
from backend.tasks import run_pipeline, _read_meta, _write_meta, JOBS_DIR

router = APIRouter(prefix="/api/jobs", tags=["Jobs"])

# Directory for uploaded FASTA files
UPLOAD_DIR = Path(__file__).resolve().parent.parent.parent / "uploads"
UPLOAD_DIR.mkdir(exist_ok=True)


# ── POST /submit ─────────────────────────────────────────────────

@router.post("/submit", response_model=JobSubmitResponse)
async def submit_job(request: JobSubmitRequest):
    """Submit a new VaxElan pipeline job for asynchronous execution."""
    job_id = f"job_{uuid.uuid4().hex[:12]}"

    _write_meta(job_id, {"job_id": job_id, "status": "pending", "progress_pct": 0.0, "current_tool": "Queued", "failed_tools": [], "error": ""})

    # Serialize for Celery (Pydantic → dict)
    request_dict = request.model_dump()

    # Dispatch to Celery worker
    run_pipeline.delay(job_id, request_dict)

    return JobSubmitResponse(job_id=job_id, status="pending")


@router.post("/submit-with-file", response_model=JobSubmitResponse)
async def submit_job_with_file(
    file: UploadFile = File(...),
    strategy: int = Form(...),
    pathogen_type: str = Form("bacteria"),
    mhci_method: str = Form("f"),
    mhcii_method: str = Form("nmel"),
    selected_tools: str = Form(""),   # comma-separated
):
    """Submit a job with an uploaded FASTA file."""
    job_id = f"job_{uuid.uuid4().hex[:12]}"

    _write_meta(job_id, {"job_id": job_id, "status": "pending", "progress_pct": 0.0, "current_tool": "Queued", "failed_tools": [], "error": ""})

    # Save uploaded file
    upload_path = UPLOAD_DIR / f"{job_id}_{file.filename}"
    with open(upload_path, "wb") as f:
        content = await file.read()
        f.write(content)

    tools = [t.strip() for t in selected_tools.split(",") if t.strip()] or None

    request_dict = {
        "input_value": str(upload_path),
        "strategy": strategy,
        "pathogen_type": pathogen_type,
        "mhci_method": mhci_method,
        "mhcii_method": mhcii_method,
        "selected_tools": tools,
        "pre_predicted_fastas": None,
        "assembly_config": None,
    }

    run_pipeline.delay(job_id, request_dict)
    return JobSubmitResponse(job_id=job_id, status="pending")


# ── GET /status ──────────────────────────────────────────────────

@router.get("/{job_id}/status", response_model=JobStatusResponse)
async def get_job_status(job_id: str):
    """Poll the status and progress of a running job."""
    meta = _read_meta(job_id)
    if not meta:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    return JobStatusResponse(
        job_id=meta["job_id"],
        status=meta.get("status", "unknown"),
        progress_pct=meta.get("progress_pct", 0.0),
        current_tool=meta.get("current_tool", ""),
        failed_tools=meta.get("failed_tools", []),
        error=meta.get("error", ""),
    )


# ── GET /results ─────────────────────────────────────────────────

@router.get("/{job_id}/results", response_model=JobResultsResponse)
async def get_job_results(job_id: str):
    """Retrieve the output files once a job is completed."""
    meta = _read_meta(job_id)
    if not meta:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if meta.get("status") != "completed":
        return JobResultsResponse(
            job_id=job_id,
            status=meta.get("status", "unknown"),
            failed_tools=meta.get("failed_tools", []),
        )

    outputs = []
    results_dir = meta.get("results_dir", "")
    for rel_path, abs_path in meta.get("outputs", {}).items():
        size = os.path.getsize(abs_path) if os.path.exists(abs_path) else 0
        outputs.append(OutputFile(
            relative_path=rel_path,
            download_url=f"/api/jobs/{job_id}/results/{rel_path}",
            size_bytes=size,
        ))

    return JobResultsResponse(
        job_id=job_id,
        status="completed",
        results_dir=results_dir,
        outputs=outputs,
        failed_tools=meta.get("failed_tools", []),
    )


# ── GET /results/{filename} ─────────────────────────────────────

@router.get("/{job_id}/results/{file_path:path}")
async def download_result_file(job_id: str, file_path: str):
    """Download a specific output file from a completed job."""
    meta = _read_meta(job_id)
    if not meta:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    abs_path = meta.get("outputs", {}).get(file_path)
    if not abs_path or not os.path.exists(abs_path):
        raise HTTPException(status_code=404, detail=f"File {file_path} not found")

    return FileResponse(
        path=abs_path,
        filename=Path(abs_path).name,
        media_type="application/octet-stream",
    )

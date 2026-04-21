"""
VaxElan 2.0 — Celery tasks wrapping core_pipeline.execute().

Each submitted job becomes a Celery task. Job metadata (status,
progress, outputs) is persisted as a JSON file in the results
directory so that the API can poll without hitting Celery.
"""
import json
import logging
import sys
from pathlib import Path
from datetime import datetime

from backend.celery_app import celery_app

# Ensure project root is importable
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

logger = logging.getLogger("vaxelan.tasks")

# ── Job metadata helpers ─────────────────────────────────────────

JOBS_DIR = _PROJECT_ROOT / "jobs"
JOBS_DIR.mkdir(exist_ok=True)


def _job_meta_path(job_id: str) -> Path:
    return JOBS_DIR / f"{job_id}.json"


def _read_meta(job_id: str) -> dict:
    path = _job_meta_path(job_id)
    if path.exists():
        return json.loads(path.read_text())
    return {}


def _write_meta(job_id: str, meta: dict):
    path = _job_meta_path(job_id)
    path.write_text(json.dumps(meta, indent=2, default=str))


# ── Celery Task ──────────────────────────────────────────────────

@celery_app.task(bind=True, name="vaxelan.run_pipeline")
def run_pipeline(self, job_id: str, request_dict: dict):
    """
    Execute a VaxElan pipeline strategy asynchronously.

    Parameters
    ----------
    job_id : str
        Unique job identifier.
    request_dict : dict
        Serialized PipelineRequest fields.
    """
    import sys
    from pathlib import Path
    _PROJECT_ROOT = Path(__file__).resolve().parent.parent
    if str(_PROJECT_ROOT) not in sys.path:
        sys.path.insert(0, str(_PROJECT_ROOT))
        
    from core_pipeline import PipelineRequest, AssemblyConfig, execute

    # Update status → running
    meta = {
        "job_id": job_id,
        "status": "running",
        "submitted_at": datetime.utcnow().isoformat(),
        "request": request_dict,
        "progress_pct": 0.0,
        "current_tool": "",
        "failed_tools": [],
        "results_dir": "",
        "outputs": {},
        "error": "",
    }
    _write_meta(job_id, meta)

    try:
        # Reconstruct AssemblyConfig if present
        assembly_cfg = None
        if request_dict.get("assembly_config"):
            assembly_cfg = AssemblyConfig(**request_dict["assembly_config"])

        request = PipelineRequest(
            input_value=request_dict["input_value"],
            strategy=request_dict["strategy"],
            pathogen_type=request_dict.get("pathogen_type", "bacteria"),
            mhci_method=request_dict.get("mhci_method", "f"),
            mhcii_method=request_dict.get("mhcii_method", "nmel"),
            selected_tools=request_dict.get("selected_tools"),
            pre_predicted_fastas=request_dict.get("pre_predicted_fastas"),
            assembly_config=assembly_cfg,
        )

        # ── Fix: Daemonic processes are not allowed to have children ──
        # Celery workers are daemonic by default, which prevents final.py 
        # from using mp.Pool(). We override this flag for the task process.
        import multiprocessing
        multiprocessing.current_process().daemon = False

        result = execute(request)

        # Update status → completed / failed
        meta["status"] = result.status
        meta["results_dir"] = result.results_dir
        meta["failed_tools"] = result.failed_tools
        meta["outputs"] = result.outputs
        meta["error"] = result.error
        meta["progress_pct"] = 100.0 if result.status == "completed" else meta["progress_pct"]
        meta["completed_at"] = datetime.utcnow().isoformat()
        _write_meta(job_id, meta)

        return {"job_id": job_id, "status": result.status}

    except BaseException as e:
        # Catch SystemExit too — pipeline modules call sys.exit() which raises
        # SystemExit(BaseException), bypassing a plain `except Exception` and
        # killing the Celery worker process with exitcode 1.
        err_msg = f"[{type(e).__name__}] {e}"
        logger.error(f"Task {job_id} failed: {err_msg}", exc_info=True)
        meta["status"] = "failed"
        meta["error"] = err_msg
        meta["completed_at"] = datetime.utcnow().isoformat()
        _write_meta(job_id, meta)
        return {"job_id": job_id, "status": "failed", "error": err_msg}

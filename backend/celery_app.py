"""
VaxElan 2.0 — Celery application configuration.

Uses Redis as both broker and result backend.
The Redis URL is read from the CELERY_BROKER_URL environment variable
(defaults to localhost for local development).
"""
import os
import sys
from pathlib import Path
from celery import Celery

# ── Ensure project root is in sys.path ───────────────────────────
_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

BROKER_URL = os.getenv("CELERY_BROKER_URL", "redis://localhost:6379/0")
RESULT_BACKEND = os.getenv("CELERY_RESULT_BACKEND", BROKER_URL)

celery_app = Celery(
    "vaxelan_worker",
    broker=BROKER_URL,
    backend=RESULT_BACKEND,
    include=["backend.tasks"],
)

celery_app.conf.update(
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    timezone="UTC",
    enable_utc=True,
    # Long-running bioinformatics tasks — generous timeouts
    task_soft_time_limit=7200,   # 2 hours soft
    task_time_limit=10800,       # 3 hours hard
    worker_prefetch_multiplier=1,
    task_acks_late=True,
)

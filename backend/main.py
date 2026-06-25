"""
EpitopeFinder — FastAPI Application Entry Point.

Run with:
    uvicorn backend.main:app --host 0.0.0.0 --port 8000 --reload
"""
from dotenv import load_dotenv
load_dotenv()

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from backend.routes import jobs, config

app = FastAPI(
    title="EpitopeFinder API",
    description="Programmatic access to the EpitopeFinder vaccine design pipeline.",
    version="2.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
)

# ── CORS ─────────────────────────────────────────────────────────
# Allow the Next.js frontend (local dev and production) to call the API.
import os

allowed_origins = [
    "http://localhost:3000",       # Next.js dev server
    "http://127.0.0.1:3000",
    "http://localhost:3001",       # Custom port for this session
    "http://127.0.0.1:3001",
    "https://*.vercel.app",        # Vercel preview deployments
]
env_origins = os.getenv("CORS_ALLOWED_ORIGINS")
if env_origins:
    allowed_origins = [o.strip() for o in env_origins.split(",") if o.strip()]

app.add_middleware(
    CORSMiddleware,
    allow_origins=allowed_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Routers ──────────────────────────────────────────────────────
app.include_router(jobs.router)
app.include_router(config.router)


# ── Production Logging & Error Diagnostics ────────────────────────
import logging
from fastapi import Request
from fastapi.responses import JSONResponse

# Create file handler for tracebacks
try:
    os.makedirs("Results", exist_ok=True)
    file_handler = logging.FileHandler("Results/backend.log")
    file_handler.setLevel(logging.ERROR)
    file_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    file_handler.setFormatter(file_formatter)
    logging.getLogger().addHandler(file_handler)
except Exception:
    pass

@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception):
    logging.error(f"Global exception caught on {request.url.path}: {exc}", exc_info=True)
    return JSONResponse(
        status_code=500,
        content={"detail": "An internal server error occurred."}
    )


# ── Health check ─────────────────────────────────────────────────
import redis
from backend.celery_app import BROKER_URL

@app.get("/api/health", tags=["Health"])
async def health():
    redis_status = "ok"
    redis_error = None
    try:
        # Simple ping check to verify connection
        r = redis.Redis.from_url(BROKER_URL, socket_timeout=2.0)
        r.ping()
    except Exception as e:
        redis_status = "disconnected"
        redis_error = str(e)

    status = "ok" if redis_status == "ok" else "degraded"
    content = {
        "status": status,
        "version": "2.0.0",
        "services": {
            "redis": {
                "status": redis_status,
                "error": redis_error
            }
        }
    }
    
    if status == "degraded":
        return JSONResponse(status_code=503, content=content)
    return content

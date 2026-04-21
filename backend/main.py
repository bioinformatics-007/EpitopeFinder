"""
VaxElan 2.0 — FastAPI Application Entry Point.

Run with:
    uvicorn backend.main:app --host 0.0.0.0 --port 8000 --reload
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from backend.routes import jobs, config

app = FastAPI(
    title="VaxElan 2.0 API",
    description="Programmatic access to the VaxElan vaccine design pipeline.",
    version="2.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
)

# ── CORS ─────────────────────────────────────────────────────────
# Allow the Next.js frontend (local dev and production) to call the API.
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",       # Next.js dev server
        "http://127.0.0.1:3000",
        "http://localhost:3001",       # Custom port for this session
        "http://127.0.0.1:3001",
        "https://*.vercel.app",        # Vercel preview deployments
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Routers ──────────────────────────────────────────────────────
app.include_router(jobs.router)
app.include_router(config.router)


# ── Health check ─────────────────────────────────────────────────
@app.get("/api/health", tags=["Health"])
async def health():
    return {"status": "ok", "version": "2.0.0"}

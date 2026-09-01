#!/bin/bash
set -e

# ── Conda environment activation ──────────────────────────────────────────────
CONDA_BASE=$(conda info --base 2>/dev/null || echo "/home/amity/miniconda3")
if [ -d "$CONDA_BASE" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate epitopefinder_web_env
fi

CONDA_PYTHON="$CONDA_BASE/envs/epitopefinder_web_env/bin/python"
CONDA_UVICORN="$CONDA_BASE/envs/epitopefinder_web_env/bin/uvicorn"
CONDA_CELERY="$CONDA_BASE/envs/epitopefinder_web_env/bin/celery"

# ── Override Redis URL for bare-metal (not Docker-compose) ────────────────────
export CELERY_BROKER_URL=redis://localhost:6379/0
export CELERY_RESULT_BACKEND=redis://localhost:6379/0
export CLBTOPE_DB=$(pwd)/tools/clbtope/clbtope/Database

# ── Start Redis (native first, Docker fallback) ───────────────────────────────
echo "Starting Redis..."
if ! redis-cli ping > /dev/null 2>&1; then
    echo "Redis not running. Starting via Docker..."
    docker start epitopepred_redis 2>/dev/null || \
        docker run -d --name epitopepred_redis -p 6379:6379 redis:7-alpine 2>/dev/null || \
        echo "WARNING: Could not start Redis via Docker. Ensure Redis is running manually."
else
    echo "Redis already running on localhost:6379 ✓"
fi

# ── Start Backend ─────────────────────────────────────────────────────────────
echo "Starting Backend (port 8000)..."
if [ -f backend.pid ]; then kill $(cat backend.pid) 2>/dev/null || true; rm -f backend.pid; fi
setsid nohup "$CONDA_UVICORN" backend.main:app --host 0.0.0.0 --port 8000 > backend.log 2>&1 &
echo $! > backend.pid
disown $! 2>/dev/null || true

# ── Start Celery Worker ───────────────────────────────────────────────────────
echo "Starting Celery Worker..."
if [ -f worker.pid ]; then kill $(cat worker.pid) 2>/dev/null || true; rm -f worker.pid; fi
setsid nohup "$CONDA_CELERY" -A backend.celery_app worker --loglevel=info > worker.log 2>&1 &
echo $! > worker.pid
disown $! 2>/dev/null || true

# ── Start Frontend ────────────────────────────────────────────────────────────
echo "Starting Frontend (port 3000)..."
if [ -f frontend.pid ]; then kill $(cat frontend.pid) 2>/dev/null || true; rm -f frontend.pid; fi
cd frontend
export NEXT_PUBLIC_API_URL=http://localhost:8000
setsid nohup npm run dev -- -H 0.0.0.0 -p 3000 > ../frontend.log 2>&1 &
echo $! > ../frontend.pid
disown $! 2>/dev/null || true
cd ..

echo ""
echo "All services started! Waiting for backend to become ready..."
sleep 5
curl -s http://localhost:8000/api/health | python3 -c "import sys,json; d=json.load(sys.stdin); print('Backend health:', d.get('status','unknown'))" 2>/dev/null || echo "Backend still starting up — check backend.log"
echo ""
echo "  Backend:  http://localhost:8000"
echo "  Frontend: http://localhost:3000"
echo "  API Docs: http://localhost:8000/docs"

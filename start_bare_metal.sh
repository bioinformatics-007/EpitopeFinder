#!/bin/bash
echo "Starting Redis..."
docker start epitopefinder_redis || docker run -d --name epitopefinder_redis -p 6379:6379 redis:7-alpine

# Auto-detect Conda base if possible, otherwise rely on current environment
CONDA_BASE=$(conda info --base 2>/dev/null)
if [ ! -z "$CONDA_BASE" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate epitopefinder_web_env
fi
export CLBTOPE_DB=$(pwd)/tools/clbtope/clbtope/Database

echo "Starting Backend..."
nohup uvicorn backend.main:app --host 0.0.0.0 --port 8000 > backend.log 2>&1 &
echo $! > backend.pid

echo "Starting Celery Worker..."
nohup celery -A backend.celery_app worker --loglevel=info > worker.log 2>&1 &
echo $! > worker.pid

echo "Starting Frontend..."
cd frontend
export NEXT_PUBLIC_API_URL=http://localhost:8000

nohup npm run dev -- -H 0.0.0.0 -p 3001 > ../frontend.log 2>&1 &
echo $! > ../frontend.pid

echo "All services started!"

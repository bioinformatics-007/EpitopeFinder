#!/bin/bash
echo "Starting Redis..."
docker start vaxelan_redis || docker run -d --name vaxelan_redis -p 6379:6379 redis:7-alpine

export PATH=/home/amity/miniconda3/envs/vaxelan_web_env/bin:$PATH
export CLBTOPE_DB=$(pwd)/tools/clbtope/clbtope/Database

echo "Starting Backend..."
nohup python3 -m uvicorn backend.main:app --host 0.0.0.0 --port 8000 > backend.log 2>&1 &
echo $! > backend.pid

echo "Starting Celery Worker..."
nohup celery -A backend.celery_app worker --loglevel=info > worker.log 2>&1 &
echo $! > worker.pid

echo "Starting Frontend..."
cd frontend
export NEXT_PUBLIC_API_URL=http://localhost:8000
nohup npm run dev -- -p 3001 > ../frontend.log 2>&1 &
echo $! > ../frontend.pid

echo "All services started!"

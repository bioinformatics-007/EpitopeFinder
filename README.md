EpitopeFinder

EpitopeFinder is a web-based bioinformatics platform for multi-epitope vaccine design against bacterial, viral, protozoal, and fungal pathogens.

It integrates 30+ bioinformatics tools into a unified pipeline with a modern web interface, enabling end-to-end vaccine design from epitope prediction to final construct assembly.

Features
6 vaccine design strategies (MHC, B-cell, physicochemical, full pipeline, assembly modes)
30+ integrated bioinformatics tools (NetMHCpan, BepiPred, SignalP, PSORTb, ESMFold, etc.)
Asynchronous task processing using Celery + Redis
Interactive web dashboard for results visualization and download
Full-stack architecture (Next.js + FastAPI)
 Architecture
Frontend (Next.js)
        │
        ▼
Backend API (FastAPI)
        │
        ▼
Celery Worker ── Redis Queue ── Bioinformatics Pipeline

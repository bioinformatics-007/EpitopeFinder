# VaxElan 2.0

**VaxElan** is a web-based bioinformatics platform for multi-epitope vaccine design against bacterial, viral, protozoal, and fungal pathogens.

It combines 30+ bioinformatics tools into a single automated pipeline with a modern web interface, supporting six distinct vaccine design strategies — from MHC epitope prediction and B-cell epitope mapping to full multi-epitope vaccine assembly.

---

## Architecture

```
         User Browser
              │
         Next.js UI  (port 3001 / 3000)
              │  REST API
         FastAPI  (port 8000)
              │  Celery task queue
         Redis  (port 6379)
              │
         Celery Worker  ─► core_pipeline.py ─► modules/
```

---

## Features

- **6 Vaccine Design Strategies**
  - Strategy 1: MHC-I / MHC-II / B-cell epitope prediction
  - Strategy 2: Physicochemical profiling (allergenicity, instability, solubility)
  - Strategy 3: Subcellular localisation & signal peptide analysis
  - Strategy 4: Full combined pipeline (Strategies 1–3)
  - Strategy 5: Assembly from pre-predicted epitope lists
  - Strategy 6: Multi-epitope vaccine assembly with adjuvants
- **30+ Integrated Tools** — NetMHCpan, NetChop, BepiPred, ESMFold, PSORTb, SignalP, WoLF PSORT, AlgPred, IAPred, ToxinPred, and more.
- **Asynchronous Job Processing** via Celery + Redis.
- **Interactive Results Dashboard** with downloadable outputs.

---

## Installation

### Prerequisites

| Requirement | Version |
|---|---|
| Python | 3.10 |
| Conda / Miniconda | Any |
| Node.js | v18+ |
| Docker | 24+ (optional, for Redis or full stack) |

### 1. Clone the Repository

```bash
git clone https://github.com/<your-username>/Vaxelan_2_0.git
cd Vaxelan_2_0
```

### 2. Set Up the Python Environment

```bash
conda env create -f vaxelan_web_env.yml
conda activate vaxelan_web_env
```

> **Full CLI environment** (includes all ML/Bioinformatics libraries):
> ```bash
> conda env create -f vaxelan_new_env.yml
> conda activate vaxelan_new
> pip install fastapi==0.111.0 uvicorn[standard]==0.29.0 celery[redis]==5.4.0 pydantic==2.7.2 python-multipart==0.0.9
> ```

### 3. Install Frontend Dependencies

```bash
cd frontend
npm install
cd ..
```

### 4. Download Required Large Files (Not in Repository)

The following files are too large for GitHub and must be downloaded separately:

| File / Directory | Description |
|---|---|
| `tools/` | 30+ bioinformatics tools (NetMHCpan, NetChop, PSORTb, etc.) |
| `model_sav/` | Pre-trained ML model files |
| `modules/*.db` | UniProt & SwissProt sequence databases |
| `modules/*.pkl` | Pre-computed k-mer index files |

See **[Tool Installation Guide](docs/tool_installation_guide.md)** for detailed download and setup instructions.

---

## Quick Start

### Option A: Docker Compose (Recommended)

```bash
docker compose up --build -d
```

- **Frontend:** http://localhost:3000
- **API Docs:** http://localhost:8000/docs

### Option B: Bare-Metal

See **[Web Execution Guide](docs/web_execution_guide.md)** for full instructions.

```bash
# Terminal 1: Redis
docker run -d --name vaxelan_redis -p 6379:6379 redis:7-alpine

# Terminal 2: Backend
conda activate vaxelan_web_env
export CLBTOPE_DB=$(pwd)/tools/clbtope/clbtope/Database
uvicorn backend.main:app --host 0.0.0.0 --port 8000

# Terminal 3: Celery Worker
conda activate vaxelan_web_env
export CLBTOPE_DB=$(pwd)/tools/clbtope/clbtope/Database
celery -A backend.celery_app worker --loglevel=info

# Terminal 4: Frontend
cd frontend
export NEXT_PUBLIC_API_URL=http://localhost:8000
npm run dev -- -p 3001
```

---

## Documentation

| Document | Description |
|---|---|
| [Web Execution Guide](docs/web_execution_guide.md) | Full Docker + bare-metal setup |
| [Shortest Run Guide](docs/shortest_run_guide.md) | Quickstart on a pre-configured system |
| [Tool Installation Guide](docs/tool_installation_guide.md) | How to install all 30+ bioinformatics tools |
| [CLI Execution Guide](docs/cli_execution_guide.md) | Running the pipeline from the terminal |

---

## License

This project is licensed under the [MIT License](LICENSE).

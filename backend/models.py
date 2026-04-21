"""
VaxElan 2.0 — Pydantic schemas for API requests and responses.
These mirror (and extend) the dataclasses in core_pipeline.py but
are proper Pydantic v2 models for FastAPI auto-validation + docs.
"""
from __future__ import annotations
from typing import Optional
from pydantic import BaseModel, Field


# ── Request Schemas ──────────────────────────────────────────────

class AssemblyConfigSchema(BaseModel):
    """Parameters for Strategy 6 vaccine assembly."""
    mode: str = Field("assemble", description="'assemble' from CSVs or 'custom' for pre-assembled sequence")
    # Assemble mode fields
    bcell_csv_path: str = ""
    ctl_csv_path: str = ""
    htl_csv_path: str = ""
    bcell_count: int = 0
    ctl_count: int = 0
    htl_count: int = 0
    assembly_order: str = "1"
    add_adjuvant: bool = False
    add_his_tag: bool = False
    # Custom mode fields
    custom_sequence: str = ""
    custom_fasta_path: str = ""
    # SASA
    run_sasa: bool = False
    sasa_csv_path: str = ""


class JobSubmitRequest(BaseModel):
    """POST /api/jobs/submit body."""
    input_value: str = Field(..., description="Path to FASTA file (server-side) or UniProt ID")
    strategy: int = Field(..., ge=1, le=6, description="Strategy number 1-6")
    pathogen_type: str = Field("bacteria", description="bacteria / virus / protozoa / fungi")
    mhci_method: str = Field("f", description="MHC-I method key (a-j)")
    mhcii_method: str = Field("nmel", description="MHC-II method key or code")
    selected_tools: Optional[list[str]] = Field(None, description="Tool names for strategies 1-3")
    pre_predicted_fastas: Optional[dict[str, str]] = Field(None, description="Strategy 5: {bcell, mhci, mhcii} FASTA paths")
    assembly_config: Optional[AssemblyConfigSchema] = Field(None, description="Strategy 6 assembly parameters")


# ── Response Schemas ─────────────────────────────────────────────

class JobSubmitResponse(BaseModel):
    job_id: str
    status: str = "pending"


class ToolStatus(BaseModel):
    name: str
    status: str  # pending / running / completed / failed


class JobStatusResponse(BaseModel):
    job_id: str
    status: str  # pending / running / completed / failed
    progress_pct: float = 0.0
    current_tool: str = ""
    failed_tools: list[str] = []
    error: str = ""


class OutputFile(BaseModel):
    relative_path: str
    download_url: str
    size_bytes: int = 0


class JobResultsResponse(BaseModel):
    job_id: str
    status: str
    results_dir: str = ""
    outputs: list[OutputFile] = []
    failed_tools: list[str] = []


# ── Config Schemas ───────────────────────────────────────────────

class MethodOption(BaseModel):
    key: str
    name: str


class StrategyInfo(BaseModel):
    number: int
    name: str
    description: str
    available_tools: list[str] = []

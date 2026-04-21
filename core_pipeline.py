#!/usr/bin/env python3
"""
core_pipeline.py — Programmatic API for VaxElan 2.0

This module provides a clean, non-interactive interface to the VaxElan
vaccine design pipeline.  It is consumed by:
  • The FastAPI backend  (backend/tasks.py)
  • Unit / integration tests
  • Any future scripted workflows

The original CLI entry-point (final.py → main()) is left untouched so
that ``python final.py`` keeps working exactly as before.
"""

import os
import sys
import shutil
import tempfile
import logging
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, field
from typing import Optional

from Bio import SeqIO

# ---------------------------------------------------------------------------
# Pydantic-style dataclasses (stdlib — no extra dep until FastAPI is wired)
# ---------------------------------------------------------------------------

@dataclass
class AssemblyConfig:
    """Configuration for Strategy 6 vaccine assembly."""
    mode: str = "assemble"                  # "assemble" or "custom"
    # Mode = assemble
    bcell_csv_path: str = ""
    ctl_csv_path: str = ""
    htl_csv_path: str = ""
    bcell_count: int = 0
    ctl_count: int = 0
    htl_count: int = 0
    assembly_order: str = "1"               # "1" - "6"
    add_adjuvant: bool = False
    add_his_tag: bool = False
    # Mode = custom
    custom_sequence: str = ""
    custom_fasta_path: str = ""
    # SASA
    run_sasa: bool = False
    sasa_csv_path: str = ""


@dataclass
class PipelineRequest:
    """All parameters needed to run a VaxElan strategy."""
    input_value: str                        # FASTA file path or UniProt ID
    strategy: int                           # 1-6
    pathogen_type: str = "bacteria"         # bacteria / virus / protozoa / fungi
    mhci_method: str = "f"                  # a-j  (default: netmhcpan_el)
    mhcii_method: str = "nmel"              # method code string
    selected_tools: Optional[list] = None   # Tool names for strategies 1-3
    # Strategy 5 specific
    pre_predicted_fastas: Optional[dict] = None  # {bcell: path, mhci: path, mhcii: path}
    # Strategy 6 specific
    assembly_config: Optional[AssemblyConfig] = None


@dataclass
class PipelineResult:
    """Structured result returned after pipeline execution."""
    job_id: str = ""
    status: str = "pending"                 # pending / running / completed / failed
    results_dir: str = ""
    failed_tools: list = field(default_factory=list)
    outputs: dict = field(default_factory=dict)
    error: str = ""


# ---------------------------------------------------------------------------
# Imports from the existing VaxElan codebase
# ---------------------------------------------------------------------------
# We add the project root to sys.path so that ``from modules.xxx`` works
# regardless of the working directory.
_PROJECT_ROOT = Path(__file__).resolve().parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

# These are imported lazily inside execute() to avoid circular imports and
# to tolerate missing optional native tools during unit-testing.

# ---------------------------------------------------------------------------
# Helper: resolve input to a local FASTA file
# ---------------------------------------------------------------------------

def _resolve_input(input_value: str, logger) -> tuple:
    """Return (input_file: Path, temp_dir: Path | None, is_uniprot: bool)."""
    from final import fetch_uniprot_sequence, validate_fasta_file

    temp_dir = Path(tempfile.mkdtemp(prefix="vaxelan_"))

    # Case 1: raw FASTA content pasted by user (starts with '>')
    if input_value.strip().startswith(">"):
        input_file = temp_dir / "input.fasta"
        input_file.write_text(input_value.strip())
        logger.info(f"[_resolve_input] Detected raw FASTA content, wrote to {input_file}")
        return input_file, temp_dir, False

    # Case 2: path to an existing file
    if os.path.isfile(input_value):
        input_file = Path(input_value)
        if not validate_fasta_file(input_file):
            raise ValueError(f"Invalid FASTA file: {input_value}")
        return input_file, temp_dir, False

    # Case 3: assume UniProt accession ID (short, no whitespace)
    is_uniprot = True
    fasta_content = fetch_uniprot_sequence(input_value.strip(), logger)
    input_file = temp_dir / f"{input_value.strip()}.fasta"
    with open(input_file, "w") as f:
        f.write(fasta_content)
    return input_file, temp_dir, is_uniprot


# ---------------------------------------------------------------------------
# Main programmatic entry-point
# ---------------------------------------------------------------------------

def execute(request: PipelineRequest) -> PipelineResult:
    """
    Run a VaxElan strategy **without any interactive prompts**.

    Parameters
    ----------
    request : PipelineRequest
        Fully-populated request object.

    Returns
    -------
    PipelineResult
        Structured result with paths to every output file.
    """
    # Lazy imports — keeps module import lightweight
    from final import (
        setup_logging,
        check_dependencies,
        split_fasta_into_batches,
        validate_pathogen_type,
        strategy_1,
        strategy_2,
        strategy_3,
        strategy_4,
        strategy_5,
        strategy_6,
        BATCH_SIZE,
    )

    # ----- Prepare result & directories -----
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_base = _PROJECT_ROOT / f"Results_{timestamp}"
    results_dir = results_base / f"strategy_{request.strategy}"
    results_dir.mkdir(parents=True, exist_ok=True)

    job_id = f"job_{timestamp}"
    result = PipelineResult(
        job_id=job_id,
        status="running",
        results_dir=str(results_dir),
    )

    logger = setup_logging(results_base)
    logger.info(f"[core_pipeline] Starting job {job_id} — strategy {request.strategy}")

    temp_dir = None
    try:
        check_dependencies()

        # ---- Strategy 6 has its own input handling ----
        if request.strategy == 6:
            failed = strategy_6(
                results_dir=results_dir,
                logger=logger,
                # Programmatic overrides (handled by refactored strategy_6)
                _mode=request.assembly_config.mode if request.assembly_config else "custom",
                _assembly_config=request.assembly_config,
            )
            result.failed_tools = failed or []
            result.status = "failed" if failed else "completed"
            return result

        # ---- Strategies 1-5: need input file + pathogen type ----
        input_file, temp_dir, is_uniprot = _resolve_input(request.input_value, logger)

        batch_files, batch_dirs, uniprot_mapping = split_fasta_into_batches(
            input_file, BATCH_SIZE, temp_dir
        )
        if not batch_files:
            raise ValueError("No valid sequences found in input.")

        confirmed_type = validate_pathogen_type(
            request.input_value if is_uniprot else str(input_file), request.pathogen_type, is_uniprot
        )

        failed: list = []

        if request.strategy == 1:
            failed = strategy_1(
                input_file, results_dir, confirmed_type, logger,
                batch_files, batch_dirs, uniprot_mapping,
                failed_tools=[],
                selected_tools=request.selected_tools or [
                    "mhc1", "mhc2", "netctl", "netchop", "bcell", "psortb"
                ],
                _mhci_method=request.mhci_method,
                _mhcii_method=request.mhcii_method,
            )

        elif request.strategy == 2:
            failed = strategy_2(
                input_file, results_dir, request.input_value, is_uniprot,
                logger, batch_files, batch_dirs, uniprot_mapping,
                failed_tools=[],
                selected_tools=request.selected_tools or [
                    "iapred", "algpred", "instability", "molwt", "wolfpsort"
                ],
            )

        elif request.strategy == 3:
            selected_tool = None
            if request.selected_tools and len(request.selected_tools) == 1:
                selected_tool = request.selected_tools[0]
            failed = strategy_3(
                input_file, confirmed_type, results_dir, logger,
                batch_files, batch_dirs, uniprot_mapping,
                selected_tool=selected_tool,
                failed_tools=[],
            )

        elif request.strategy == 4:
            failed = strategy_4(
                input_file, confirmed_type, results_dir, logger,
                batch_files, batch_dirs, uniprot_mapping,
                failed_tools=[],
            )

        elif request.strategy == 5:
            failed = strategy_5(
                input_file, confirmed_type, results_dir, logger,
                batch_files, batch_dirs, uniprot_mapping,
                failed_tools=[],
                _pre_predicted_fastas=request.pre_predicted_fastas,
                _mhci_method=request.mhci_method,
                _mhcii_method=request.mhcii_method,
            )

        result.failed_tools = failed or []
        result.status = "completed"

    except BaseException as e:
        # Catch sys.exit() (raises SystemExit which is a BaseException, not Exception)
        logger.error(f"[core_pipeline] Job {job_id} failed: {e}", exc_info=True)
        result.status = "failed"
        result.error = f"[{type(e).__name__}] {e}"

    finally:
        if temp_dir and temp_dir.exists():
            shutil.rmtree(temp_dir, ignore_errors=True)
            logger.info(f"Cleaned up temp dir: {temp_dir}")

    # ---- Collect output file paths ----
    if result.status == "completed" and Path(result.results_dir).exists():
        for fpath in Path(result.results_dir).rglob("*"):
            if fpath.is_file():
                rel = fpath.relative_to(result.results_dir)
                result.outputs[str(rel)] = str(fpath)

    logger.info(f"[core_pipeline] Job {job_id} finished — status={result.status}")
    return result

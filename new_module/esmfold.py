from __future__ import annotations

import requests
import os
import logging
import shutil
import subprocess
from pathlib import Path
import numpy as np
import platform
import sys
import time
import matplotlib.pyplot as plt
from typing import Union, Optional, Dict, Any


def run_esmfold(
    sequence: str,
    output_dir: Union[str, Path],
    logger: Optional[logging.Logger] = None
) -> Dict[str, Any]:
    """
    Predict 3D structure using ESMFold, generate pLDDT plot,
    and create PyMOL script with special coloring for linkers.
    """
    if logger is None:
        logger = logging.getLogger('VaxElan')

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pdb_path = output_dir / "vaccine_3D_model.pdb"
    pml_path = output_dir / "visualize_structure.pml"
    plot_path = output_dir / "plddt_per_residue.png"

    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"

    try:
        clean_seq = "".join(c for c in sequence.strip().upper() if c.isalpha() and c in 'ACDEFGHIKLMNPQRSTVWY')
        if not clean_seq:
            raise ValueError("Empty or invalid sequence after cleaning")

        print(f"\033[94mRequesting 3D structure from ESMAtlas... (length: {len(clean_seq)} aa)\033[0m")

        response = requests.post(url, data=clean_seq.encode('utf-8'), timeout=200)

        if response.status_code != 200:
            error_msg = f"ESM Atlas API error {response.status_code} - {response.text[:300]}"
            logger.error(error_msg)
            print(f"\033[91m{error_msg}\033[0m")
            return {"status": -1, "error": error_msg}

        pdb_content = response.text

        with open(pdb_path, "w", encoding="utf-8") as f:
            f.write(pdb_content)

        # ── Extract per-residue pLDDT ────────────────────────────────────────
        residue_plddts = []
        last_res = None

        for line in pdb_content.splitlines():
            if line.startswith("ATOM") and " CA " in line:
                try:
                    res_num = int(line[22:26].strip())
                    val_str = line[60:66].strip()
                    if val_str:
                        val = float(val_str)
                        plddt_val = val if 0 <= val <= 1 else (val / 100 if 0 <= val <= 100 else None)
                        if plddt_val is not None:
                            if last_res is None or res_num != last_res:
                                residue_plddts.append(plddt_val)
                                last_res = res_num
                except:
                    continue

        if not residue_plddts:
            plddt_values = [float(line[60:66].strip()) for line in pdb_content.splitlines()
                            if line.startswith("ATOM") and len(line) >= 66 and line[60:66].strip()]
            residue_plddts = [v/100 if v > 1 else v for v in plddt_values if 0 <= v <= 100]

        avg_plddt_raw = np.mean(residue_plddts) if residue_plddts else 0.0
        avg_plddt_percent = avg_plddt_raw * 100

        print(f"\033[92mStructure predicted. Average pLDDT: {avg_plddt_percent:.1f}% (0-100 scale)\033[0m")

        # ── Generate per-residue pLDDT plot ──────────────────────────────────
        if residue_plddts:
            plt.figure(figsize=(12, 6), dpi=120)
            residues = np.arange(1, len(residue_plddts) + 1)
            plddt_percent = np.array(residue_plddts) * 100

            plt.plot(residues, plddt_percent, color='royalblue', linewidth=2.5, label='pLDDT')
            plt.axhspan(90, 100, facecolor='green', alpha=0.12)
            plt.axhspan(70, 90,  facecolor='lime',   alpha=0.12)
            plt.axhspan(50, 70,  facecolor='yellow', alpha=0.12)
            plt.axhspan(0,  50,  facecolor='orange', alpha=0.12)

            plt.axhline(y=avg_plddt_percent, color='red', linestyle='--',
                        label=f'Average: {avg_plddt_percent:.1f}%')

            plt.title('Per-Residue Model Confidence (pLDDT) - ESMFold')
            plt.xlabel('Residue Position')
            plt.ylabel('pLDDT (%)')
            plt.ylim(0, 105)
            plt.grid(True, alpha=0.3)
            plt.legend(loc='lower right')
            plt.tight_layout()
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"\033[92mPer-residue pLDDT plot saved → {plot_path}\033[0m")

        create_pymol_script(str(pdb_path), str(pml_path), sequence.upper())
        launched = launch_pymol_from_pymol_env(str(pml_path), logger)

        return {
            "status": 0,
            "pdb_file": str(pdb_path),
            "avg_confidence": round(float(avg_plddt_percent), 2),
            "pml_file": str(pml_path),
            "plddt_plot": str(plot_path) if residue_plddts else None,
            "pymol_launched": launched
        }

    except Exception as e:
        logger.exception("Error in ESMFold")
        print(f"\033[91mError: {str(e)}\033[0m")
        return {"status": -1, "error": str(e)}


def create_pymol_script(pdb_file: str, pml_file: str, full_sequence: str) -> None:
    """
    Generate PyMOL script with:
    - pLDDT-based cartoon coloring
    - Distinct colors for specific linkers (AYY, GPGPG, KK)
    """
    script = f"""\
# VaxElan generated PyMOL script - {os.path.basename(pdb_file)}
# pLDDT confidence coloring + special linker highlighting

load {pdb_file}, vaccine

# Background and basic view
bg_color white
hide everything
show cartoon, vaccine
spectrum b, red_orange_yellow_green_blue, vaccine, minimum=0, maximum=100
set cartoon_fancy_helices, 1
set ray_opaque_background, on
set antialias, 3

# ── Special coloring for linkers ─────────────────────────────────────────
# Note: selections are approximate - may need manual adjustment for junctions

# AYY linker (magenta)
select linker_ayy, vaccine and (resn ALA+TYR+TYR or (resi {find_linker_positions(full_sequence, "AYY")}))

# GPGPG linker (cyan)
select linker_gpgpg, vaccine and (resn GLY+PRO+GLY+PRO+GLY or (resi {find_linker_positions(full_sequence, "GPGPG")}))

# KK linker (yellow)
select linker_kk, vaccine and (resn LYS+LYS or (resi {find_linker_positions(full_sequence, "KK")}))

# Apply colors and style to linkers
color magenta, linker_ayy
color cyan,    linker_gpgpg
color yellow,  linker_kk

# Show linkers as sticks for better visibility
show sticks, linker_ayy or linker_gpgpg or linker_kk
set stick_radius, 0.18, linker_ayy or linker_gpgpg or linker_kk

# Zoom
zoom vaccine, buffer=4

# Optional: hide selection names in viewer
hide labels
"""

    with open(pml_file, "w", encoding="utf-8") as f:
        f.write(script)

    print(f"\033[92mPyMOL script with linker coloring created → {pml_file}\033[0m")


def find_linker_positions(sequence: str, linker: str) -> str:
    """
    Find residue numbers (1-based) where linker sequence appears.
    Returns comma-separated list for PyMOL resi selection.
    """
    positions = []
    linker_len = len(linker)
    for i in range(len(sequence) - linker_len + 1):
        if sequence[i:i+linker_len] == linker:
            # Starting residue number (1-based)
            start_res = i + 1
            # Add all residues of this linker instance
            for j in range(linker_len):
                positions.append(str(start_res + j))
    return ",".join(positions) if positions else "none"


def launch_pymol_from_pymol_env(pml_file: str, logger: logging.Logger) -> bool:
    """Launch PyMOL using 'conda run' from the 'pymol_env' environment."""
    pymol_env_name = "pymol_env"  # ← Change if your env name is different!

    print(f"  → Attempting to launch PyMOL from conda env: {pymol_env_name}")

    cmd = [
        "conda", "run",
        "--no-capture-output",
        "-n", pymol_env_name,
        "pymol", "-qx", pml_file
    ]

    try:
        print("  → Running: conda run -n pymol_env pymol -qx ...")
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True
        )

        time.sleep(4.0)

        if proc.poll() is None:
            print("  → PyMOL launched successfully from 'pymol_env'!")
            return True
        else:
            stdout, stderr = proc.communicate(timeout=5)
            print(f"  → conda run failed (return code {proc.returncode}):")
            if stderr:
                print(f"     Error: {stderr.strip()}")
            return False

    except Exception as e:
        print(f"  → Failed to launch via conda run: {e}")

    # Fallback...
    print("  → Falling back to direct launch with DISPLAY/XAUTHORITY...")
    pymol_bin = shutil.which("pymol")
    if not pymol_bin:
        _show_manual_instructions(pml_file)
        return False

    cmd_fallback = [pymol_bin, "-qx", pml_file]
    env = os.environ.copy()

    display_candidates = [os.environ.get('DISPLAY', ':0'), ':0', ':1', ':0.0', ':10.0']
    xauth_candidates = [
        os.environ.get('XAUTHORITY'),
        os.path.expanduser("~/.Xauthority"),
        f"/run/user/{os.getuid()}/gdm/Xauthority",
        f"/tmp/.Xauthority"
    ]

    for disp in display_candidates:
        for xauth in [x for x in xauth_candidates if x and os.path.exists(x)]:
            env['DISPLAY'] = disp
            env['XAUTHORITY'] = xauth
            print(f"     Trying fallback: DISPLAY={disp} XAUTHORITY={xauth}")

            try:
                proc = subprocess.Popen(
                    cmd_fallback,
                    env=env,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    start_new_session=True
                )
                time.sleep(3.0)
                if proc.poll() is None:
                    print(f"  → Success with fallback DISPLAY={disp}")
                    return True
            except:
                continue

    _show_manual_instructions(pml_file)
    logger.warning("PyMOL GUI auto-launch failed after all attempts")
    return False


def _show_manual_instructions(pml_file: str) -> None:
    print("\n" + "═"*70)
    print(" AUTO-LAUNCH FAILED")
    print(" PyMOL is installed in a separate env ('pymol_env').")
    print(" Please open the structure manually:")
    print()
    print(f"   1. conda activate pymol_env")
    print(f"   2. pymol {pml_file}")
    print(f"      OR")
    print(f"      pymol {Path(pml_file).parent / 'vaccine_3D_model.pdb'}")
    print()
    print("═"*70 + "\n")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    dummy_seq = "MADEAPKPRMKKKKGVKRKRT" * 5 + "GPGPG" + "AYY" * 3 + "KK" * 2
    result = run_esmfold(dummy_seq, "./test_output")
    print(result)

#/usr/bin/env python3
"""
Automate dataset preparation and Nextflow runs for every species in TODOLIST.

Look at TODOLIST 
Create a analysis (species name) dir in analysis/ 
For each one. run prepare_dataset.py to download genomes and prepare samplesheet in to the dir
Then run nextflow main.nf with that samplesheet. Be sure to give them a work dir in their species dir so they dont clash
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple


ROOT_DIR = Path(__file__).resolve().parent.parent
TODO_PATH = ROOT_DIR / "TODOLIST"
ANALYSIS_DIR = ROOT_DIR / "analysis"
PREPARE_SCRIPT = ROOT_DIR / "prepare_scripts" / "prepare_dataset.py"
MAIN_NF_PATH = ROOT_DIR / "main.nf"


def slugify(value: str) -> str:
    """Return a filesystem-friendly slug."""
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return slug.strip("_") or "species"


def read_todolist() -> List[Tuple[str, str]]:
    """Parse species and taxid pairs from TODOLIST."""
    if not TODO_PATH.exists():
        raise FileNotFoundError(f"TODOLIST not found: {TODO_PATH}")
    entries: List[Tuple[str, str]] = []
    pattern = re.compile(r"^(?P<name>.+?)\s*\((?P<taxid>\d+)\)$")
    for line in TODO_PATH.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        match = pattern.match(stripped)
        if not match:
            raise ValueError(f"Cannot parse line in TODOLIST: {stripped}")
        entries.append((match.group("name").strip(), match.group("taxid")))
    return entries


def run_command(cmd: List[str]) -> None:
    """Execute a subprocess command in the repository root."""
    print(f"\n$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=ROOT_DIR)


def run_prepare_dataset(species: str, taxid: str) -> Path:
    """Download assemblies and build the samplesheet for the species."""
    cmd = [
        "pixi",
        "run",
        "python",
        str(PREPARE_SCRIPT),
        "--taxid",
        taxid,
        "--outdir",
        str(ANALYSIS_DIR),
        species,
    ]
    run_command(cmd)
    species_dir = ANALYSIS_DIR / slugify(species)
    samplesheet = species_dir / "samplesheet.csv"
    if not samplesheet.exists():
        raise FileNotFoundError(f"Samplesheet was not created: {samplesheet}")
    return species_dir


def run_nextflow(species_dir: Path) -> None:
    """Launch Nextflow using the species-specific samplesheet and work dir."""
    samplesheet = species_dir / "samplesheet.csv"
    work_dir = species_dir / "work"
    work_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "nextflow",
        "run",
        str(MAIN_NF_PATH),
        "--samplesheet",
        str(samplesheet),
        "-work-dir",
        str(work_dir),
    ]
    run_command(cmd)


def main() -> None:
    species_list = read_todolist()
    if not species_list:
        print("No species found in TODOLIST.")
        return
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    for species, taxid in species_list:
        print(f"\n=== Processing {species} (taxid {taxid}) ===")
        species_dir = run_prepare_dataset(species, taxid)
        run_nextflow(species_dir)


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        print(f"Command failed with exit code {exc.returncode}: {' '.join(exc.cmd)}", file=sys.stderr)
        sys.exit(exc.returncode)
    except Exception as exc:  # noqa: BLE001 - want to bubble up friendly errors
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
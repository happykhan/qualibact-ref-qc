#!/usr/bin/env python3
"""
Prepare Singularity images for every container declared in the Nextflow modules.

The script scans the modules directory, extracts `container 'image:tag'` entries,
and downloads corresponding Singularity (`.sif`) files into a user-specified
directory. It also configures APPTAINER cache/tmp locations so the downloads
remain self-contained for shared clusters (e.g. BMRC).
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, Set

DIRECT_CONTAINER_PATTERN = re.compile(r"container\s+['\"]([^'\"]+)['\"]")
CONSTANT_CONTAINER_PATTERN = re.compile(r"containerImage\s*=\s*['\"]([^'\"]+)['\"]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pull Singularity images for containers declared in modules/*.nf"
    )
    parser.add_argument(
        "--modules-dir",
        default="modules",
        help="Directory containing Nextflow modules (default: %(default)s).",
    )
    parser.add_argument(
        "--singularity-dir",
        default="singularity",
        help="Directory where .sif images will be stored (default: %(default)s).",
    )
    parser.add_argument(
        "--cache-dir",
        default="apptainer-cache",
        help="Directory to use for APPTAINER_CACHEDIR (default: %(default)s).",
    )
    parser.add_argument(
        "--tmp-dir",
        default="tmp",
        help="Directory to use for APPTAINER_TMPDIR (default: %(default)s).",
    )
    parser.add_argument(
        "--singularity-command",
        default="singularity",
        help="Singularity/Apptainer executable to use (default: %(default)s).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the commands without executing Singularity pulls.",
    )
    return parser.parse_args()


def collect_containers(modules_dir: Path) -> Dict[str, Set[Path]]:
    containers: Dict[str, Set[Path]] = {}
    for nf_file in sorted(modules_dir.glob("*.nf")):
        try:
            text = nf_file.read_text()
        except OSError as exc:
            print(f"[WARN] Unable to read {nf_file}: {exc}", file=sys.stderr)
            continue
        for pattern in (DIRECT_CONTAINER_PATTERN, CONSTANT_CONTAINER_PATTERN):
            for match in pattern.finditer(text):
                image = match.group(1).strip()
                if not image:
                    continue
                containers.setdefault(image, set()).add(nf_file)
        # If both patterns match the same literal, we only add once due to set.
    return containers


def sanitize_image_name(image: str) -> str:
    # Replace slashes/colons to create a reproducible filename.
    sanitized = image.replace("/", "_").replace(":", "_")
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", sanitized)
    return f"{sanitized}.sif"


def ensure_directories(*paths: Path) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def pull_image(
    image: str,
    sif_path: Path,
    singularity_cmd: str,
    env: Dict[str, str],
    dry_run: bool,
) -> None:
    cmd = [
        singularity_cmd,
        "pull",
        str(sif_path),
        f"docker://{image}",
    ]
    print(" ".join(cmd))
    if dry_run:
        return
    subprocess.run(cmd, check=True, env=env)


def main() -> None:
    args = parse_args()

    modules_dir = Path(args.modules_dir).resolve()
    singularity_dir = Path(args.singularity_dir).resolve()
    cache_dir = Path(args.cache_dir).resolve()
    tmp_dir = Path(args.tmp_dir).resolve()

    if not modules_dir.exists():
        sys.exit(f"Modules directory not found: {modules_dir}")

    ensure_directories(singularity_dir, cache_dir, tmp_dir)

    containers = collect_containers(modules_dir)
    if not containers:
        sys.exit(f"No container declarations found in {modules_dir}")

    env = os.environ.copy()
    env["APPTAINER_CACHEDIR"] = str(cache_dir)
    env["APPTAINER_TMPDIR"] = str(tmp_dir)

    print(f"export APPTAINER_CACHEDIR={cache_dir}")
    print(f"export APPTAINER_TMPDIR={tmp_dir}")

    for image, sources in sorted(containers.items()):
        sif_name = sanitize_image_name(image)
        sif_path = singularity_dir / sif_name
        if sif_path.exists():
            if args.dry_run:
                print(f"[INFO] Would remove existing {sif_path}")
            else:
                print(f"[INFO] Removing existing {sif_path}")
                try:
                    sif_path.unlink()
                except OSError as exc:
                    print(f"[ERROR] Failed to remove {sif_path}: {exc}", file=sys.stderr)
                    sys.exit(1)
        print(f"[INFO] Pulling {image} referenced in {', '.join(str(s) for s in sorted(sources))}")
        try:
            pull_image(image, sif_path, args.singularity_command, env, args.dry_run)
        except subprocess.CalledProcessError as exc:
            print(f"[ERROR] Singularity pull failed for {image}: {exc}", file=sys.stderr)
            sys.exit(1)

    print(f"[DONE] Images stored in {singularity_dir}")


if __name__ == "__main__":
    main()

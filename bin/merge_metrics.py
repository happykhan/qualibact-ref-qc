#!/usr/bin/env python3
"""
Merge assembly stats, base composition, and CheckM outputs; generate summaries/plots.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def strip_ext(name: str) -> str:
    for ext in (".fasta", ".fa", ".fna"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def load_assembly_stats(files: Iterable[Path]) -> Dict[str, Dict[str, str]]:
    records = {}
    for path in files:
        rows = read_tsv(path)
        if not rows:
            continue
        row = rows[0]
        sample_id = strip_ext(row.get("filename", Path(path).name))
        records[sample_id] = row
    return records


def load_base_comp(files: Iterable[Path]) -> Dict[str, Dict[str, str]]:
    records = {}
    for path in files:
        rows = read_tsv(path)
        if not rows:
            continue
        row = rows[0]
        sample_id = strip_ext(row.get("Filename", Path(path).name))
        records[sample_id] = row
    return records


def load_checkm(files: Iterable[Path]) -> Dict[str, Dict[str, str]]:
    records = {}
    for path in files:
        rows = read_tsv(path)
        if not rows:
            continue
        row = rows[0]
        sample_id = strip_ext(row.get("Name", Path(path).stem))
        records[sample_id] = row
    return records


def compute_gc(row: Dict[str, str]) -> float | None:
    try:
        a = int(row.get("A", 0))
        t = int(row.get("T", 0))
        g = int(row.get("G", 0))
        c = int(row.get("C", 0))
    except ValueError:
        return None
    total = a + t + g + c
    if total == 0:
        return None
    return (g + c) / total


def percentile(data: List[float], pct: float) -> float:
    if not data:
        return math.nan
    data = sorted(data)
    k = (len(data) - 1) * pct
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return data[int(k)]
    d0 = data[int(f)] * (c - k)
    d1 = data[int(c)] * (k - f)
    return d0 + d1


def summarize(values: List[float]) -> Dict[str, float]:
    if not values:
        return {k: math.nan for k in ["min", "max", "q1", "mean", "median", "q3", "iqr"]}
    q1 = percentile(values, 0.25)
    q3 = percentile(values, 0.75)
    return {
        "min": min(values),
        "max": max(values),
        "q1": q1,
        "mean": statistics.fmean(values),
        "median": statistics.median(values),
        "q3": q3,
        "iqr": q3 - q1,
    }


def histogram(data: List[float], bins: int = 20) -> List[Tuple[float, float, int]]:
    if not data:
        return []
    mn, mx = min(data), max(data)
    if mn == mx:
        return [(mn, mx, len(data))]
    width = (mx - mn) / bins
    edges = [mn + i * width for i in range(bins + 1)]
    counts = [0] * bins
    for v in data:
        idx = min(int((v - mn) / width), bins - 1)
        counts[idx] += 1
    return [(edges[i], edges[i + 1], counts[i]) for i in range(bins)]


def load_metadata_genome_sizes(path: Path | None) -> List[float]:
    if not path or not path.exists():
        return []
    entries = json.loads(path.read_text())
    sizes = []
    for e in entries:
        stats = e.get("assembly_stats") or {}
        size = stats.get("total_sequence_length")
        if size is None:
            continue
        try:
            sizes.append(float(size))
        except (TypeError, ValueError):
            continue
    return sizes


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assembly-stats", nargs="+", required=True, type=Path)
    ap.add_argument("--base-comp", nargs="+", required=True, type=Path)
    ap.add_argument("--checkm", nargs="+", required=True, type=Path)
    ap.add_argument("--metadata-json", type=Path)
    ap.add_argument("--species", default="species")
    ap.add_argument("--samplesheet", type=Path)
    ap.add_argument("--outdir", type=Path, default=Path("."))
    args = ap.parse_args()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    asm_stats = load_assembly_stats(args.assembly_stats)
    base_comp = load_base_comp(args.base_comp)
    checkm = load_checkm(args.checkm)

    all_ids = set(asm_stats) | set(base_comp) | set(checkm)

    merged_path = outdir / "merged_metrics.tsv"
    with merged_path.open("w", newline="") as handle:
        fieldnames = [
            "sample_id",
            "total_length",
            "gc_fraction",
            "genome_size_checkm",
            "gc_content_checkm",
            "cds_count",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for sid in sorted(all_ids):
            row = {"sample_id": sid}
            stats_row = asm_stats.get(sid, {})
            bc_row = base_comp.get(sid, {})
            ch_row = checkm.get(sid, {})

            row["total_length"] = stats_row.get("total_length")
            row["genome_size_checkm"] = ch_row.get("Genome_Size")
            row["gc_content_checkm"] = ch_row.get("GC_Content")
            row["cds_count"] = ch_row.get("Total_Coding_Sequences")

            gc_val = compute_gc(bc_row)
            row["gc_fraction"] = f"{gc_val:.5f}" if gc_val is not None else ""

            writer.writerow(row)

    # Numeric series
    def to_floats(values):
        out = []
        for v in values:
            try:
                out.append(float(v))
            except (TypeError, ValueError):
                continue
        return out

    genome_sizes_assembly = to_floats([asm_stats[sid].get("total_length") for sid in asm_stats])
    genome_sizes_checkm = to_floats([checkm[sid].get("Genome_Size") for sid in checkm])
    gc_vals = [compute_gc(base_comp[sid]) for sid in base_comp]
    gc_vals = [v for v in gc_vals if v is not None]
    cds_vals = to_floats([checkm[sid].get("Total_Coding_Sequences") for sid in checkm])
    ref_sizes = load_metadata_genome_sizes(args.metadata_json)

    # Summary CSV
    summary_rows = []
    summary_rows.append({"metric": "genome_size_assembly", **summarize(genome_sizes_assembly)})
    summary_rows.append({"metric": "genome_size_checkm", **summarize(genome_sizes_checkm)})
    summary_rows.append({"metric": "gc_fraction", **summarize(gc_vals)})
    summary_rows.append({"metric": "cds_count", **summarize(cds_vals)})
    if ref_sizes:
        summary_rows.append({"metric": "genome_size_reference", **summarize(ref_sizes)})

    with (outdir / "summary.csv").open("w", newline="") as handle:
        fieldnames = ["metric", "min", "max", "q1", "mean", "median", "q3", "iqr"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    # Species metrics (ranges)
    def rounded_range(values: List[float]) -> Tuple[float, float] | None:
        if not values:
            return None
        return (round(min(values), -3), round(max(values), -3))

    species_rows = []
    if genome_sizes_assembly:
        gr = rounded_range(genome_sizes_assembly)
        gcr = rounded_range(gc_vals) if gc_vals else None
        species_rows.append(
            {
                "source": "assemblies",
                "genome_size_min": gr[0] if gr else "",
                "genome_size_max": gr[1] if gr else "",
                "gc_min": gcr[0] if gcr else "",
                "gc_max": gcr[1] if gcr else "",
            }
        )
    if ref_sizes:
        gr = rounded_range(ref_sizes)
        species_rows.append(
            {
                "source": "reference",
                "genome_size_min": gr[0] if gr else "",
                "genome_size_max": gr[1] if gr else "",
                "gc_min": "",
                "gc_max": "",
            }
        )

    with (outdir / f"{args.species}_metrics.csv").open("w", newline="") as handle:
        fieldnames = ["source", "genome_size_min", "genome_size_max", "gc_min", "gc_max"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in species_rows:
            writer.writerow(row)

    # Histogram data
    hist_data = []
    bins = histogram(genome_sizes_assembly + ref_sizes, bins=20) if (genome_sizes_assembly or ref_sizes) else []
    ref_set = histogram(ref_sizes, bins=len(bins)) if ref_sizes and bins else []
    ref_counts = {(round(b[0], 6), round(b[1], 6)): b[2] for b in ref_set}
    with (outdir / "genome_size_histogram.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["bin_start", "bin_end", "assemblies", "reference"])
        for b in bins:
            key = (round(b[0], 6), round(b[1], 6))
            writer.writerow([b[0], b[1], b[2], ref_counts.get(key, 0)])

    # CDS vs genome size scatter data
    with (outdir / "cds_vs_genome_size.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", "cds_count", "genome_size"])
        for sid, row in checkm.items():
            try:
                cds = float(row.get("Total_Coding_Sequences"))
                gsize = float(row.get("Genome_Size"))
            except (TypeError, ValueError):
                continue
            writer.writerow([sid, cds, gsize])

    # Optional plots if matplotlib is available
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return

    if bins:
        fig, ax = plt.subplots()
        asm_counts = [b[2] for b in bins]
        bin_edges = [b[0] for b in bins] + [bins[-1][1]]
        ax.hist(
            [genome_sizes_assembly, ref_sizes],
            bins=bin_edges,
            label=["assemblies", "reference"],
            alpha=0.7,
        )
        ax.set_xlabel("Genome size")
        ax.set_ylabel("Count")
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "genome_size_histogram.png", dpi=200)
        plt.close(fig)

    if checkm:
        xs = []
        ys = []
        for row in checkm.values():
            try:
                ys.append(float(row.get("Total_Coding_Sequences")))
                xs.append(float(row.get("Genome_Size")))
            except (TypeError, ValueError):
                continue
        if xs and ys:
            fig, ax = plt.subplots()
            ax.scatter(xs, ys, alpha=0.6)
            ax.set_xlabel("Genome size")
            ax.set_ylabel("Total coding sequences")
            fig.tight_layout()
            fig.savefig(outdir / "cds_vs_genome_size.png", dpi=200)
            plt.close(fig)


if __name__ == "__main__":
    main()

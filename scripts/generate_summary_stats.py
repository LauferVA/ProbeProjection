#!/usr/bin/env python3
"""
generate_summary_stats.py — Produce Table 1 and Table 3 data for ProbeProjection.

Reads pipeline output (manifest.json, per-profile coverage and diversity JSONs)
and generates summary statistics for:
  - Table 1: Technology Summary (profiles, probes, coverage, disparity by tech)
  - Table 3: HPRC Cohort Composition (superpopulation sample sizes and power)

When pipeline data is unavailable, the --demo flag generates tables using known
catalog metadata with placeholder values for projection results.

Usage:
    python generate_summary_stats.py [OPTIONS]
    python generate_summary_stats.py --demo
    python generate_summary_stats.py --help

Options:
    --results-dir DIR   Path to pipeline results directory (default: ../results/)
    --output-dir DIR    Path for summary output (default: ../results/summaries/)
    --format FMT        Output format: json, tsv, or both (default: both)
    --demo              Generate tables with catalog metadata and placeholders

Dependencies: standard library.  Uses scipy.stats.norm for power calculations
if available; otherwise falls back to a manual normal approximation.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import sys
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Import power computation from sibling module if available, else inline
# ---------------------------------------------------------------------------
try:
    # If run from the scripts directory or the module is on sys.path
    _script_dir = pathlib.Path(__file__).resolve().parent
    if str(_script_dir) not in sys.path:
        sys.path.insert(0, str(_script_dir))
    from power_analysis import min_detectable_delta, SUPERPOPULATIONS, P0, TARGET_POWER
except ImportError:
    # Inline fallback: minimal power computation
    try:
        from scipy.stats import norm as _norm

        def _ppf(p: float) -> float:
            return float(_norm.ppf(p))
    except ImportError:
        def _ppf(p: float) -> float:
            if p <= 0.0 or p >= 1.0:
                raise ValueError(f"p must be in (0, 1), got {p}")
            if p < 0.5:
                return -_ppf(1.0 - p)
            t = math.sqrt(-2.0 * math.log(1.0 - p))
            c0, c1, c2 = 2.515517, 0.802853, 0.010328
            d1, d2, d3 = 1.432788, 0.189269, 0.001308
            return t - (c0 + c1 * t + c2 * t * t) / (
                1.0 + d1 * t + d2 * t * t + d3 * t * t * t
            )

    P0 = 0.95
    TARGET_POWER = 0.80
    SUPERPOPULATIONS = {
        "AFR": {"full_name": "African / African American", "samples": 53, "haplotypes": 106},
        "NFE": {"full_name": "Non-Finnish European",       "samples": 35, "haplotypes": 70},
        "EAS": {"full_name": "East Asian",                  "samples": 24, "haplotypes": 48},
        "SAS": {"full_name": "South Asian",                 "samples": 19, "haplotypes": 38},
        "AMR": {"full_name": "Admixed American",            "samples": 5,  "haplotypes": 10},
    }

    def min_detectable_delta(
        n_target: int, p0: float, alpha: float, power: float, n_ref: int | None = None
    ) -> float:
        z_alpha = _ppf(1.0 - alpha / 2.0)
        z_beta = _ppf(power)
        lo, hi = 0.0, p0
        for _ in range(200):
            mid = (lo + hi) / 2.0
            p1 = p0 - mid
            if p1 <= 0.0:
                hi = mid
                continue
            p_bar = (p0 + p1) / 2.0
            if n_ref is not None:
                se_h0 = math.sqrt(p_bar * (1.0 - p_bar) * (1.0 / n_ref + 1.0 / n_target))
            else:
                se_h0 = math.sqrt(p_bar * (1.0 - p_bar) / n_target)
            if n_ref is not None:
                se_h1 = math.sqrt(p0 * (1.0 - p0) / n_ref + p1 * (1.0 - p1) / n_target)
            else:
                se_h1 = math.sqrt(p1 * (1.0 - p1) / n_target)
            if se_h0 < 1e-15 or se_h1 < 1e-15:
                hi = mid
                continue
            ncp = mid / se_h1
            achieved = ncp - z_alpha * se_h0 / se_h1
            if achieved < z_beta:
                lo = mid
            else:
                hi = mid
        return min((lo + hi) / 2.0, 1.0)


# ---------------------------------------------------------------------------
# Known catalog metadata for demo mode
# ---------------------------------------------------------------------------
# From the spec and manifest: 16 technologies, 1628 total profiles, 625 GRCh38-eligible.
# Profile and probe counts come from the normalized profile catalog.
DEMO_TECHNOLOGY_CATALOG: List[Dict[str, Any]] = [
    {"technology": "FISH",           "profiles": 319, "total_probes": 1276,    "probe_size_min_bp": 50000,   "probe_size_max_bp": 1000000},
    {"technology": "NGS_Panel",      "profiles": 184, "total_probes": 92000,   "probe_size_min_bp": 120,     "probe_size_max_bp": 500},
    {"technology": "HLA",            "profiles": 38,  "total_probes": 380,     "probe_size_min_bp": 500,     "probe_size_max_bp": 50000},
    {"technology": "PGx",            "profiles": 23,  "total_probes": 690,     "probe_size_min_bp": 200,     "probe_size_max_bp": 100000},
    {"technology": "CMA",            "profiles": 15,  "total_probes": 6850000, "probe_size_min_bp": 25,      "probe_size_max_bp": 60},
    {"technology": "MLPA",           "profiles": 12,  "total_probes": 600,     "probe_size_min_bp": 50,      "probe_size_max_bp": 100},
    {"technology": "OGM",            "profiles": 8,   "total_probes": 500000,  "probe_size_min_bp": 6,       "probe_size_max_bp": 6},
    {"technology": "Methylation",    "profiles": 7,   "total_probes": 850000,  "probe_size_min_bp": 50,      "probe_size_max_bp": 50},
    {"technology": "Karyotype",      "profiles": 5,   "total_probes": 100,     "probe_size_min_bp": 5000000, "probe_size_max_bp": 10000000},
    {"technology": "SNP_Array",      "profiles": 4,   "total_probes": 750000,  "probe_size_min_bp": 25,      "probe_size_max_bp": 50},
    {"technology": "Sanger",         "profiles": 3,   "total_probes": 45,      "probe_size_min_bp": 200,     "probe_size_max_bp": 1000},
    {"technology": "ddPCR",          "profiles": 3,   "total_probes": 30,      "probe_size_min_bp": 60,      "probe_size_max_bp": 150},
    {"technology": "RT_PCR",         "profiles": 2,   "total_probes": 24,      "probe_size_min_bp": 60,      "probe_size_max_bp": 150},
    {"technology": "WGS",            "profiles": 1,   "total_probes": 0,       "probe_size_min_bp": None,    "probe_size_max_bp": None},
    {"technology": "WES",            "profiles": 1,   "total_probes": 22000,   "probe_size_min_bp": 120,     "probe_size_max_bp": 300},
    {"technology": "RNA_Seq",        "profiles": 1,   "total_probes": 0,       "probe_size_min_bp": None,    "probe_size_max_bp": None},
]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_manifest(results_dir: pathlib.Path) -> Optional[Dict]:
    """Load manifest.json from the results directory."""
    manifest_path = results_dir / "manifest.json"
    if not manifest_path.exists():
        print(f"Warning: manifest not found at {manifest_path}", file=sys.stderr)
        return None
    with open(manifest_path) as f:
        return json.load(f)


def discover_technologies(results_dir: pathlib.Path) -> List[str]:
    """List technology subdirectories under per_profile/."""
    per_profile_dir = results_dir / "per_profile"
    if not per_profile_dir.is_dir():
        return []
    return sorted(
        d.name
        for d in per_profile_dir.iterdir()
        if d.is_dir()
    )


def load_per_profile_data(
    results_dir: pathlib.Path, technology: str
) -> Tuple[List[Dict], List[Dict]]:
    """Load all coverage and diversity JSONs for a technology.

    Returns (coverage_list, diversity_list).
    """
    tech_dir = results_dir / "per_profile" / technology
    coverage_files: List[Dict] = []
    diversity_files: List[Dict] = []

    if not tech_dir.is_dir():
        return coverage_files, diversity_files

    for fpath in sorted(tech_dir.iterdir()):
        if not fpath.name.endswith(".json"):
            continue
        try:
            with open(fpath) as f:
                data = json.load(f)
        except (json.JSONDecodeError, OSError) as exc:
            print(f"Warning: could not read {fpath}: {exc}", file=sys.stderr)
            continue

        if fpath.name.endswith("_coverage.json"):
            coverage_files.append(data)
        elif fpath.name.endswith("_diversity.json"):
            diversity_files.append(data)

    return coverage_files, diversity_files


# ---------------------------------------------------------------------------
# Table 1: Technology Summary
# ---------------------------------------------------------------------------
def compute_table1_from_data(
    results_dir: pathlib.Path,
    manifest: Optional[Dict],
) -> List[Dict]:
    """Build Table 1 rows from actual pipeline output files."""
    technologies = discover_technologies(results_dir)
    if not technologies:
        return []

    rows: List[Dict] = []
    for tech in technologies:
        coverage_list, diversity_list = load_per_profile_data(results_dir, tech)
        if not coverage_list and not diversity_list:
            continue

        n_profiles = len(coverage_list)
        total_probes = sum(c.get("total_probes", 0) for c in coverage_list)
        probes_universal = sum(c.get("probes_universal", 0) for c in coverage_list)
        pct_universal = (probes_universal / total_probes * 100) if total_probes > 0 else None

        # Disparity from diversity files
        total_evaluated = sum(d.get("total_probes_evaluated", 0) for d in diversity_list)
        total_disparate = sum(d.get("probes_with_significant_disparity", 0) for d in diversity_list)
        pct_gap = (total_disparate / total_evaluated * 100) if total_evaluated > 0 else None

        # Most affected superpopulation: aggregate mean_coverage_by_population
        pop_coverage_sums: Dict[str, float] = defaultdict(float)
        pop_coverage_counts: Dict[str, int] = defaultdict(int)
        for d in diversity_list:
            mcbp = d.get("mean_coverage_by_population", {})
            for pop, cov in mcbp.items():
                pop_coverage_sums[pop] += cov
                pop_coverage_counts[pop] += 1

        most_affected = None
        if pop_coverage_sums:
            pop_means = {
                pop: pop_coverage_sums[pop] / pop_coverage_counts[pop]
                for pop in pop_coverage_sums
            }
            most_affected = min(pop_means, key=lambda p: pop_means[p])

        # Probe size range: not directly in coverage JSON; compute from
        # per-probe results if available.  For now, mark as not available
        # unless we can infer it.
        probe_size_min = None
        probe_size_max = None

        rows.append({
            "technology": tech.upper(),
            "n_profiles": n_profiles,
            "total_probes": total_probes,
            "probe_size_min_bp": probe_size_min,
            "probe_size_max_bp": probe_size_max,
            "probe_size_range": _format_size_range(probe_size_min, probe_size_max),
            "pct_universal_coverage": round(pct_universal, 2) if pct_universal is not None else None,
            "pct_population_specific_gap": round(pct_gap, 2) if pct_gap is not None else None,
            "most_affected_superpopulation": most_affected,
        })

    return rows


def compute_table1_demo() -> List[Dict]:
    """Build Table 1 rows in demo mode using catalog metadata and placeholders."""
    rows: List[Dict] = []
    for entry in DEMO_TECHNOLOGY_CATALOG:
        tech = entry["technology"]
        n_profiles = entry["profiles"]
        total_probes = entry["total_probes"]
        size_min = entry["probe_size_min_bp"]
        size_max = entry["probe_size_max_bp"]

        rows.append({
            "technology": tech,
            "n_profiles": n_profiles,
            "total_probes": total_probes,
            "probe_size_min_bp": size_min,
            "probe_size_max_bp": size_max,
            "probe_size_range": _format_size_range(size_min, size_max),
            "pct_universal_coverage": None,  # placeholder: requires pipeline
            "pct_population_specific_gap": None,  # placeholder: requires pipeline
            "most_affected_superpopulation": None,  # placeholder: requires pipeline
        })
    return rows


def _format_size_range(lo: Optional[int], hi: Optional[int]) -> str:
    """Human-readable probe size range."""
    if lo is None and hi is None:
        return "N/A"
    if lo == hi:
        return _format_bp(lo)
    return f"{_format_bp(lo)}-{_format_bp(hi)}"


def _format_bp(n: Optional[int]) -> str:
    """Format a base-pair count with appropriate unit."""
    if n is None:
        return "N/A"
    if n >= 1_000_000:
        return f"{n / 1_000_000:.1f}Mb"
    if n >= 1_000:
        return f"{n / 1_000:.0f}kb"
    return f"{n}bp"


# ---------------------------------------------------------------------------
# Table 3: HPRC Cohort Composition
# ---------------------------------------------------------------------------
def compute_table3(manifest: Optional[Dict]) -> List[Dict]:
    """Build Table 3 rows from manifest cohort metadata or SUPERPOPULATIONS fallback."""
    # Try to pull from manifest
    if manifest is not None:
        rm = manifest.get("resource_metadata", {})
        pop_samples = rm.get("registry_population_counts", {})
        pop_haplotypes = rm.get("registry_haplotype_counts", {})
    else:
        pop_samples = {}
        pop_haplotypes = {}

    rows: List[Dict] = []
    total_samples = 0
    total_haplotypes = 0

    for code, info in SUPERPOPULATIONS.items():
        n_samples = pop_samples.get(code, info["samples"])
        n_haplotypes = pop_haplotypes.get(code, info["haplotypes"])

        # Compute minimum detectable disparity at 80% power, uncorrected
        delta = min_detectable_delta(
            n_target=n_haplotypes,
            p0=P0,
            alpha=0.05,
            power=TARGET_POWER,
        )

        statistical_note = f"~{delta * 100:.0f}% min detectable disparity"
        if n_haplotypes <= 10:
            statistical_note = f"Underpowered; Fisher exact only; ~{delta * 100:.0f}% min detectable"
        elif n_haplotypes <= 40:
            statistical_note = f"Limited power; {statistical_note}"
        elif n_haplotypes <= 50:
            statistical_note = f"Moderate power; {statistical_note}"
        elif n_haplotypes <= 80:
            statistical_note = f"Well powered; {statistical_note}"
        else:
            statistical_note = f"Best powered; {statistical_note}"

        total_samples += n_samples
        total_haplotypes += n_haplotypes

        rows.append({
            "superpopulation_code": code,
            "superpopulation_name": info["full_name"],
            "samples": n_samples,
            "haplotypes": n_haplotypes,
            "min_detectable_disparity_uncorrected": round(delta, 4),
            "min_detectable_disparity_pct": round(delta * 100, 1),
            "statistical_note": statistical_note,
        })

    # Add total row
    rows.append({
        "superpopulation_code": "Total",
        "superpopulation_name": "",
        "samples": total_samples,
        "haplotypes": total_haplotypes,
        "min_detectable_disparity_uncorrected": None,
        "min_detectable_disparity_pct": None,
        "statistical_note": "",
    })

    return rows


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------
def write_json_file(data: Any, path: pathlib.Path) -> None:
    """Write data as formatted JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Wrote {path}")


def write_tsv_file(rows: List[Dict], path: pathlib.Path) -> None:
    """Write rows as a TSV file."""
    if not rows:
        print(f"Skipping {path} (no rows)", file=sys.stderr)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0].keys())
    with open(path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            values = []
            for c in columns:
                v = row.get(c)
                values.append("" if v is None else str(v))
            f.write("\t".join(values) + "\n")
    print(f"Wrote {path}")


# ---------------------------------------------------------------------------
# Console output
# ---------------------------------------------------------------------------
def print_table1(rows: List[Dict]) -> None:
    """Print Table 1 to stdout."""
    print()
    print("=" * 100)
    print("TABLE 1: Technology Summary")
    print("=" * 100)

    header = (
        f"{'Technology':<16} {'Profiles':>8} {'Probes':>12} "
        f"{'Size Range':<15} {'% Universal':>11} {'% Pop Gap':>9} "
        f"{'Most Affected':>14}"
    )
    print(header)
    print("-" * 100)

    for r in rows:
        pct_u = f"{r['pct_universal_coverage']:.1f}" if r["pct_universal_coverage"] is not None else "[pending]"
        pct_g = f"{r['pct_population_specific_gap']:.1f}" if r["pct_population_specific_gap"] is not None else "[pending]"
        most = r["most_affected_superpopulation"] or "[pending]"
        probes_str = f"{r['total_probes']:,}" if r["total_probes"] else "0"

        print(
            f"{r['technology']:<16} {r['n_profiles']:>8} {probes_str:>12} "
            f"{r['probe_size_range']:<15} {pct_u:>11} {pct_g:>9} "
            f"{most:>14}"
        )
    print("-" * 100)
    print()


def print_table3(rows: List[Dict]) -> None:
    """Print Table 3 to stdout."""
    print()
    print("=" * 90)
    print("TABLE 3: HPRC Cohort Composition")
    print("=" * 90)

    header = (
        f"{'Code':<6} {'Full Name':<28} {'Samples':>8} "
        f"{'Haplotypes':>10} {'Min Disp %':>10}  {'Note'}"
    )
    print(header)
    print("-" * 90)

    for r in rows:
        disp = f"{r['min_detectable_disparity_pct']:.1f}" if r["min_detectable_disparity_pct"] is not None else ""
        print(
            f"{r['superpopulation_code']:<6} {r['superpopulation_name']:<28} "
            f"{r['samples']:>8} {r['haplotypes']:>10} "
            f"{disp:>10}  {r['statistical_note']}"
        )
    print("-" * 90)
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate summary statistics for Table 1 (Technology Summary) "
            "and Table 3 (HPRC Cohort Composition) from ProbeProjection "
            "pipeline output."
        ),
    )
    parser.add_argument(
        "--results-dir",
        type=pathlib.Path,
        default=None,
        help=(
            "Path to pipeline results directory. "
            "Default: <script_dir>/../results/"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=None,
        help=(
            "Path for summary output files. "
            "Default: <script_dir>/../results/summaries/"
        ),
    )
    parser.add_argument(
        "--format",
        choices=["json", "tsv", "both"],
        default="both",
        help="Output format (default: both).",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help=(
            "Generate tables with known catalog metadata and placeholder "
            "values for projection results that require pipeline execution."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    script_dir = pathlib.Path(__file__).resolve().parent
    results_dir = args.results_dir or (script_dir.parent / "results")
    output_dir = args.output_dir or (script_dir.parent / "results" / "summaries")

    # Load manifest
    manifest = load_manifest(results_dir)

    # ---- Table 1 ----
    if args.demo:
        print("Running in --demo mode: using catalog metadata with placeholders.")
        table1_rows = compute_table1_demo()
    else:
        table1_rows = compute_table1_from_data(results_dir, manifest)
        if not table1_rows:
            print(
                "No per-profile data found. Use --demo for placeholder output, "
                "or ensure pipeline results exist under per_profile/.",
                file=sys.stderr,
            )

    # ---- Table 3 ----
    table3_rows = compute_table3(manifest)

    # ---- Print to stdout ----
    if table1_rows:
        print_table1(table1_rows)
    print_table3(table3_rows)

    # ---- Write files ----
    is_demo = args.demo

    table1_doc = {
        "table": "Table 1: Technology Summary",
        "source": "demo catalog metadata" if is_demo else "pipeline per-profile output",
        "note": (
            "Columns pct_universal_coverage, pct_population_specific_gap, and "
            "most_affected_superpopulation require pipeline execution and are "
            "null in demo mode."
        ) if is_demo else None,
        "rows": table1_rows,
    }
    table3_doc = {
        "table": "Table 3: HPRC Cohort Composition",
        "source": "manifest.json + power analysis",
        "reference_coverage_assumed": P0,
        "power_target": TARGET_POWER,
        "alpha_base": 0.05,
        "rows": table3_rows,
    }

    fmt = args.format

    if fmt in ("json", "both"):
        write_json_file(table1_doc, output_dir / "table1_technology_summary.json")
        write_json_file(table3_doc, output_dir / "table3_cohort_composition.json")

    if fmt in ("tsv", "both"):
        write_tsv_file(table1_rows, output_dir / "table1_technology_summary.tsv")
        write_tsv_file(table3_rows, output_dir / "table3_cohort_composition.tsv")

    print("\nDone.")


if __name__ == "__main__":
    main()

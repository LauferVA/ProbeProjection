#!/usr/bin/env python3
"""
power_analysis.py — Minimum detectable disparity for ProbeProjection cohort sizes.

Computes the minimum detectable absolute disparity (%) between a reference
population coverage rate and each HPRC superpopulation, using a two-proportion
z-test framework at 80% power across several multiple-testing scenarios.

Usage:
    python power_analysis.py [--output-dir DIR] [--format {json,tsv,both}]
    python power_analysis.py --help

Outputs:
    results/summaries/power_analysis.json
    results/summaries/power_analysis.tsv
    Formatted table to stdout.

Dependencies: standard library.  Uses scipy.stats.norm if available; falls back
to a rational approximation of the inverse normal CDF otherwise.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import sys
from typing import Dict, List, Tuple

# ---------------------------------------------------------------------------
# Normal distribution utilities — scipy preferred, manual fallback
# ---------------------------------------------------------------------------
try:
    from scipy.stats import norm as _norm

    def _ppf(p: float) -> float:
        """Inverse CDF of the standard normal distribution."""
        return float(_norm.ppf(p))

except ImportError:
    # Rational approximation (Abramowitz & Stegun 26.2.23, |error| < 4.5e-4).
    def _ppf(p: float) -> float:
        """Inverse CDF of the standard normal via rational approximation."""
        if p <= 0.0 or p >= 1.0:
            raise ValueError(f"p must be in (0, 1), got {p}")
        if p < 0.5:
            return -_ppf(1.0 - p)
        t = math.sqrt(-2.0 * math.log(1.0 - p))
        c0, c1, c2 = 2.515517, 0.802853, 0.010328
        d1, d2, d3 = 1.432788, 0.189269, 0.001308
        return t - (c0 + c1 * t + c2 * t * t) / (1.0 + d1 * t + d2 * t * t + d3 * t * t * t)


# ---------------------------------------------------------------------------
# HPRC cohort data
# ---------------------------------------------------------------------------
SUPERPOPULATIONS: Dict[str, Dict] = {
    "AFR": {"full_name": "African / African American", "samples": 53, "haplotypes": 106},
    "NFE": {"full_name": "Non-Finnish European",       "samples": 35, "haplotypes": 70},
    "EAS": {"full_name": "East Asian",                  "samples": 24, "haplotypes": 48},
    "SAS": {"full_name": "South Asian",                 "samples": 19, "haplotypes": 38},
    "AMR": {"full_name": "Admixed American",            "samples": 5,  "haplotypes": 10},
}

# Reference coverage rate assumed for the "no-disparity" baseline.
P0 = 0.95

# Multiple-testing scenarios: label -> number of tests.
CORRECTION_SCENARIOS: List[Tuple[str, int]] = [
    ("Uncorrected",            1),
    ("Bonferroni (4, FISH)",   4),
    ("Bonferroni (500, NGS)",  500),
    ("Bonferroni (6.85M, CMA)", 6_850_000),
]

# Target power.
TARGET_POWER = 0.80


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------
def min_detectable_delta(
    n_target: int,
    p0: float,
    alpha: float,
    power: float,
    n_ref: int | None = None,
) -> float:
    """Return the minimum detectable absolute disparity (p0 - p1) for a
    two-proportion z-test.

    Parameters
    ----------
    n_target : int
        Number of haplotypes in the target (smaller) group.
    p0 : float
        Coverage rate in the reference group (assumed known / large-sample).
    alpha : float
        Significance level (after any multiple-testing correction).
    power : float
        Desired statistical power (probability of detecting a true effect).
    n_ref : int or None
        Number of haplotypes in the reference group.  If None, the reference
        is treated as infinitely large (conservative for the target group).

    Returns
    -------
    float
        Minimum detectable absolute disparity as a proportion (0-1).
        Returns 1.0 if the computation yields a non-physical result.
    """
    z_alpha = _ppf(1.0 - alpha / 2.0)
    z_beta = _ppf(power)

    # Iterative solver: find p1 such that the non-centrality parameter
    # equals z_beta.  Binary search on delta = p0 - p1.
    lo, hi = 0.0, p0  # delta in [0, p0] so p1 in [0, p0]
    for _ in range(200):
        mid = (lo + hi) / 2.0
        p1 = p0 - mid

        if p1 <= 0.0:
            hi = mid
            continue

        p_bar = (p0 + p1) / 2.0  # pooled proportion (equal group weight approx)

        # Standard error under H0 (pooled)
        if n_ref is not None:
            se_h0 = math.sqrt(p_bar * (1.0 - p_bar) * (1.0 / n_ref + 1.0 / n_target))
        else:
            se_h0 = math.sqrt(p_bar * (1.0 - p_bar) / n_target)

        # Standard error under H1
        if n_ref is not None:
            se_h1 = math.sqrt(p0 * (1.0 - p0) / n_ref + p1 * (1.0 - p1) / n_target)
        else:
            se_h1 = math.sqrt(p1 * (1.0 - p1) / n_target)

        if se_h0 < 1e-15 or se_h1 < 1e-15:
            hi = mid
            continue

        # Achieved power for this delta
        ncp = mid / se_h1
        achieved = ncp - z_alpha * se_h0 / se_h1

        if achieved < z_beta:
            lo = mid
        else:
            hi = mid

    delta = (lo + hi) / 2.0
    return min(delta, 1.0)


def build_power_table(
    p0: float = P0,
    power: float = TARGET_POWER,
) -> List[Dict]:
    """Build the full power-analysis table.

    Returns a list of dicts, one per (superpopulation, correction_scenario).
    """
    rows: List[Dict] = []
    for code, info in SUPERPOPULATIONS.items():
        n_hap = info["haplotypes"]
        for scenario_label, n_tests in CORRECTION_SCENARIOS:
            alpha_corrected = 0.05 / n_tests
            delta = min_detectable_delta(
                n_target=n_hap,
                p0=p0,
                alpha=alpha_corrected,
                power=power,
            )
            rows.append({
                "superpopulation": code,
                "full_name": info["full_name"],
                "samples": info["samples"],
                "haplotypes": n_hap,
                "correction_scenario": scenario_label,
                "n_tests": n_tests,
                "alpha_corrected": round(alpha_corrected, 12),
                "reference_coverage": p0,
                "power": power,
                "min_detectable_disparity": round(delta, 4),
                "min_detectable_disparity_pct": round(delta * 100, 2),
            })
    return rows


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------
def print_table(rows: List[Dict]) -> None:
    """Print a formatted table to stdout."""
    # Group by superpopulation for readability.
    header = (
        f"{'Pop':>5}  {'Haplotypes':>10}  "
        f"{'Uncorrected':>12}  {'Bonf(4)':>10}  "
        f"{'Bonf(500)':>10}  {'Bonf(6.85M)':>12}"
    )
    separator = "-" * len(header)

    print()
    print("Minimum Detectable Absolute Disparity (%) at 80% Power")
    print(f"Reference coverage = {P0 * 100:.0f}%, two-proportion z-test")
    print()
    print(header)
    print(separator)

    pops_seen: set = set()
    for row in rows:
        code = row["superpopulation"]
        if code not in pops_seen:
            pops_seen.add(code)
            # Collect all scenarios for this pop
            pop_rows = [r for r in rows if r["superpopulation"] == code]
            vals = {r["correction_scenario"]: r["min_detectable_disparity_pct"] for r in pop_rows}
            n_hap = row["haplotypes"]
            note = " *" if n_hap <= 10 else ""
            print(
                f"{code:>5}  {n_hap:>10}  "
                f"{vals.get('Uncorrected', 'N/A'):>11.1f}%  "
                f"{vals.get('Bonferroni (4, FISH)', 'N/A'):>9.1f}%  "
                f"{vals.get('Bonferroni (500, NGS)', 'N/A'):>9.1f}%  "
                f"{vals.get('Bonferroni (6.85M, CMA)', 'N/A'):>11.1f}%{note}"
            )

    print(separator)
    print()
    print("* AMR (n=10 haplotypes) is underpowered for all but extreme effects.")
    print("  Results should use Fisher's exact test and be flagged LOW_POWER.")
    print()


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------
def write_json(rows: List[Dict], path: pathlib.Path) -> None:
    """Write rows as a JSON file."""
    output = {
        "description": "Power analysis: minimum detectable disparity by superpopulation",
        "methodology": "Two-proportion z-test, binary search for delta at target power",
        "reference_coverage": P0,
        "target_power": TARGET_POWER,
        "base_alpha": 0.05,
        "rows": rows,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Wrote {path}")


def write_tsv(rows: List[Dict], path: pathlib.Path) -> None:
    """Write rows as a TSV file."""
    columns = [
        "superpopulation", "full_name", "samples", "haplotypes",
        "correction_scenario", "n_tests", "alpha_corrected",
        "reference_coverage", "power",
        "min_detectable_disparity", "min_detectable_disparity_pct",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            f.write("\t".join(str(row[c]) for c in columns) + "\n")
    print(f"Wrote {path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute minimum detectable disparity for HPRC superpopulation "
            "sample sizes under several multiple-testing correction scenarios."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=None,
        help=(
            "Directory for output files. "
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
        "--reference-coverage",
        type=float,
        default=P0,
        help=f"Assumed reference coverage rate (default: {P0}).",
    )
    parser.add_argument(
        "--power",
        type=float,
        default=TARGET_POWER,
        help=f"Target statistical power (default: {TARGET_POWER}).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    output_dir = args.output_dir
    if output_dir is None:
        output_dir = pathlib.Path(__file__).resolve().parent.parent / "results" / "summaries"

    rows = build_power_table(p0=args.reference_coverage, power=args.power)

    print_table(rows)

    if args.format in ("json", "both"):
        write_json(rows, output_dir / "power_analysis.json")
    if args.format in ("tsv", "both"):
        write_tsv(rows, output_dir / "power_analysis.tsv")


if __name__ == "__main__":
    main()

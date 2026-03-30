#!/usr/bin/env python3
"""
generate_bed_files.py

Read ProbeProjection pipeline output (per-profile coverage and diversity JSONs)
and produce GA4GH BED v1.0 files with population-stratified coverage metrics.

Output BED schema (BED4+ with 9 additional columns):
    chrom  chromStart  chromEnd  name  score  product_id  technology
    coverage_AFR  coverage_EAS  coverage_NFE  coverage_SAS  coverage_AMR
    disparity_magnitude  significant

Usage:
    python generate_bed_files.py --results-dir ../results --output-dir ../distribution/bed
    python generate_bed_files.py --demo          # write 20 simulated probes
"""

import argparse
import json
import os
import pathlib
import random
import sys

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
POPULATIONS = ["AFR", "EAS", "NFE", "SAS", "AMR"]

BED_HEADER_LINES = [
    "#GA4GH BED v1.0",
    "#chrom\tchromStart\tchromEnd\tname\tscore\tproduct_id\ttechnology"
    "\tcoverage_AFR\tcoverage_EAS\tcoverage_NFE\tcoverage_SAS\tcoverage_AMR"
    "\tdisparity_magnitude\tsignificant",
]

# Mapping from chromosome name to sort key (chr1..chr22, chrX, chrY, chrM, then
# anything else alphabetically).
_CHR_ORDER = {f"chr{i}": i for i in range(1, 23)}
_CHR_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25})


def chrom_sort_key(chrom: str) -> tuple:
    """Return a tuple suitable for sorting chromosomes in karyotypic order."""
    return (_CHR_ORDER.get(chrom, 99), chrom)


# ---------------------------------------------------------------------------
# Data extraction helpers
# ---------------------------------------------------------------------------

def _extract_probe_records(coverage_path: pathlib.Path,
                           diversity_path: pathlib.Path,
                           technology: str) -> list[dict]:
    """
    Parse one profile's coverage + diversity JSONs and return a list of
    per-probe BED-row dicts.

    Each coverage JSON contains ``probe_results`` -- one entry per
    (probe, haplotype) pair.  We aggregate across haplotypes to compute
    per-population coverage rates, then join with the diversity JSON for
    disparity information.
    """
    with open(coverage_path) as fh:
        cov = json.load(fh)
    with open(diversity_path) as fh:
        div = json.load(fh)

    profile_id = cov["profile_id"]
    probe_results = cov.get("probe_results", [])

    if not probe_results:
        return []

    # -- Build per-probe, per-population haplotype counts ----------------
    # We need the ancestry registry to know which sample belongs to which
    # population.  The diversity JSON carries ``sample_size_per_population``
    # and ``population_haplotype_counts``, but the individual sample->pop
    # mapping is NOT in the JSON.  Therefore we compute an *overall*
    # coverage rate per probe (fraction of haplotypes that are mapped) and
    # use the diversity JSON's ``mean_coverage_by_population`` for the
    # per-population breakdown.

    # Unique probes in this profile
    probe_ids = sorted({r["probe_id"] for r in probe_results})

    # Per-probe aggregation
    probe_agg: dict[str, dict] = {}
    for pid in probe_ids:
        rows = [r for r in probe_results if r["probe_id"] == pid]
        total = len(rows)
        mapped = sum(1 for r in rows if r.get("is_mapped", False))

        # Derive chromosome and interval.  If any row has a non-null
        # mapped_interval we use that; otherwise fall back to the
        # chromosome field (start/end unknown -> use 0/0 placeholder).
        chrom = rows[0].get("chromosome", "chrUn")
        chrom_start = 0
        chrom_end = 0
        for r in rows:
            mi = r.get("mapped_interval")
            if mi is not None:
                # mapped_interval is expected as [start, end] or
                # {"start": ..., "end": ...}
                if isinstance(mi, list) and len(mi) >= 2:
                    chrom_start = int(mi[0])
                    chrom_end = int(mi[1])
                    break
                elif isinstance(mi, dict):
                    chrom_start = int(mi.get("start", 0))
                    chrom_end = int(mi.get("end", 0))
                    break

        overall_rate = mapped / total if total > 0 else 0.0

        probe_agg[pid] = {
            "chrom": chrom,
            "chromStart": chrom_start,
            "chromEnd": chrom_end,
            "name": pid,
            "overall_rate": overall_rate,
            "product_id": profile_id,
            "technology": technology.upper(),
        }

    # -- Join per-population coverage from diversity JSON -----------------
    mean_cov_by_pop = div.get("mean_coverage_by_population", {})
    disparity_mag = div.get("max_population_disparity", 0.0)

    # Build a set of probe IDs flagged as significantly disparate
    disparate_ids = set()
    for entry in div.get("disparate_probes", []):
        if isinstance(entry, dict):
            disparate_ids.add(entry.get("probe_id", ""))
        elif isinstance(entry, str):
            disparate_ids.add(entry)

    records = []
    for pid, agg in probe_agg.items():
        # For single-probe profiles the diversity JSON's mean_coverage is
        # the per-population rate for that probe.  For multi-probe profiles
        # it is an average, so this is an approximation until we have
        # per-probe per-population breakdowns.
        pop_coverage = {}
        for pop in POPULATIONS:
            pop_coverage[pop] = mean_cov_by_pop.get(pop, 0.0)

        score = min(int(round(agg["overall_rate"] * 1000)), 1000)

        records.append({
            "chrom": agg["chrom"],
            "chromStart": agg["chromStart"],
            "chromEnd": agg["chromEnd"],
            "name": agg["name"],
            "score": score,
            "product_id": agg["product_id"],
            "technology": agg["technology"],
            "coverage_AFR": round(pop_coverage.get("AFR", 0.0), 4),
            "coverage_EAS": round(pop_coverage.get("EAS", 0.0), 4),
            "coverage_NFE": round(pop_coverage.get("NFE", 0.0), 4),
            "coverage_SAS": round(pop_coverage.get("SAS", 0.0), 4),
            "coverage_AMR": round(pop_coverage.get("AMR", 0.0), 4),
            "disparity_magnitude": round(disparity_mag, 4),
            "significant": 1 if pid in disparate_ids else 0,
        })

    return records


def _sort_bed_records(records: list[dict]) -> list[dict]:
    """Sort records by chromosome (karyotypic) then chromStart."""
    return sorted(records, key=lambda r: (
        chrom_sort_key(r["chrom"]),
        r["chromStart"],
        r["chromEnd"],
        r["name"],
    ))


def _format_bed_line(rec: dict) -> str:
    """Format one record as a tab-delimited BED line."""
    fields = [
        rec["chrom"],
        str(rec["chromStart"]),
        str(rec["chromEnd"]),
        rec["name"],
        str(rec["score"]),
        rec["product_id"],
        rec["technology"],
        f'{rec["coverage_AFR"]:.4f}',
        f'{rec["coverage_EAS"]:.4f}',
        f'{rec["coverage_NFE"]:.4f}',
        f'{rec["coverage_SAS"]:.4f}',
        f'{rec["coverage_AMR"]:.4f}',
        f'{rec["disparity_magnitude"]:.4f}',
        str(rec["significant"]),
    ]
    return "\t".join(fields)


def _write_bed(path: pathlib.Path, records: list[dict]) -> None:
    """Write a sorted BED file with header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    sorted_recs = _sort_bed_records(records)
    with open(path, "w") as fh:
        for hdr in BED_HEADER_LINES:
            fh.write(hdr + "\n")
        for rec in sorted_recs:
            fh.write(_format_bed_line(rec) + "\n")


# ---------------------------------------------------------------------------
# Pipeline results processing
# ---------------------------------------------------------------------------

def process_results(results_dir: pathlib.Path,
                    output_dir: pathlib.Path,
                    technology_filter: list[str] | None) -> None:
    """
    Walk results/per_profile/{technology}/ and generate BED files.

    One BED per technology plus a combined all-technologies file.
    """
    per_profile_dir = results_dir / "per_profile"

    if not per_profile_dir.is_dir():
        print(f"[INFO] per_profile directory not found: {per_profile_dir}")
        print("       Run the ProbeProjection pipeline first to generate results,")
        print("       or use --demo to produce an example BED file.")
        return

    # Discover technologies (subdirectories of per_profile/)
    tech_dirs = sorted([
        d for d in per_profile_dir.iterdir()
        if d.is_dir()
    ])

    if not tech_dirs:
        print(f"[INFO] No technology subdirectories found under {per_profile_dir}")
        return

    all_records: list[dict] = []
    technologies_processed = 0

    for tech_dir in tech_dirs:
        tech_name = tech_dir.name.upper()

        if technology_filter and tech_name not in [t.upper() for t in technology_filter]:
            continue

        # Find coverage/diversity JSON pairs
        coverage_files = sorted(tech_dir.glob("*_coverage.json"))
        if not coverage_files:
            print(f"[INFO] No coverage files found for technology {tech_name}")
            continue

        tech_records: list[dict] = []
        for cov_path in coverage_files:
            # Derive diversity path from coverage path
            div_path = cov_path.parent / cov_path.name.replace(
                "_coverage.json", "_diversity.json"
            )
            if not div_path.exists():
                print(f"[WARN] Missing diversity file: {div_path}")
                continue

            try:
                recs = _extract_probe_records(cov_path, div_path, tech_name)
                tech_records.extend(recs)
            except (json.JSONDecodeError, KeyError) as exc:
                print(f"[WARN] Failed to parse {cov_path.name}: {exc}")
                continue

        if tech_records:
            bed_path = output_dir / f"{tech_name.lower()}.probe_coverage.bed"
            _write_bed(bed_path, tech_records)
            all_records.extend(tech_records)
            technologies_processed += 1
            print(f"[OK]   {bed_path.name}: {len(tech_records)} probe records")

    if all_records:
        combined_path = output_dir / "all_technologies.probe_coverage.bed"
        _write_bed(combined_path, all_records)
        print(f"[OK]   {combined_path.name}: {len(all_records)} probe records (combined)")

    # -- Summary statistics -----------------------------------------------
    print()
    print("=== Summary ===")
    print(f"Technologies processed : {technologies_processed}")
    print(f"Total probe records    : {len(all_records)}")
    if all_records:
        scores = [r["score"] for r in all_records]
        significant = sum(1 for r in all_records if r["significant"])
        print(f"Score range            : {min(scores)} - {max(scores)}")
        print(f"Significant disparities: {significant} / {len(all_records)}")
        # Per-population mean coverage
        for pop in POPULATIONS:
            key = f"coverage_{pop}"
            vals = [r[key] for r in all_records]
            mean_val = sum(vals) / len(vals) if vals else 0.0
            print(f"Mean coverage {pop:>3s}      : {mean_val:.4f}")
    print(f"Output directory       : {output_dir}")


# ---------------------------------------------------------------------------
# Demo mode
# ---------------------------------------------------------------------------

_DEMO_GENES = [
    ("chr1",  10000,  15000, "GENE_A"),
    ("chr1",  50000,  55000, "GENE_B"),
    ("chr2",  20000,  25000, "GENE_C"),
    ("chr2",  80000,  85000, "GENE_D"),
    ("chr3",  30000,  35000, "GENE_E"),
    ("chr5",  10000,  14000, "GENE_F"),
    ("chr7",  45000,  50000, "GENE_G"),
    ("chr7",  90000,  95000, "GENE_H"),
    ("chr8",  12000,  17000, "GENE_I"),
    ("chr9",  60000,  65000, "GENE_J"),
    ("chr10", 25000,  30000, "GENE_K"),
    ("chr11", 35000,  40000, "GENE_L"),
    ("chr12", 55000,  60000, "GENE_M"),
    ("chr13", 15000,  20000, "GENE_N"),
    ("chr15", 40000,  45000, "GENE_O"),
    ("chr17", 39000000, 39100000, "ERBB2"),
    ("chr19", 10000,  15000, "GENE_Q"),
    ("chr20", 50000,  55000, "GENE_R"),
    ("chr22", 20000,  25000, "GENE_S"),
    ("chrX",  30000,  35000, "GENE_T"),
]


def _generate_demo(output_dir: pathlib.Path) -> None:
    """Generate a small demo BED file with 20 simulated probes."""
    rng = random.Random(42)
    technologies = ["NGS", "FISH", "MICROARRAY", "RTPCR"]
    records = []

    for chrom, start, end, gene in _DEMO_GENES:
        tech = rng.choice(technologies)
        product_id = f"DEMO_{gene}_{tech}"
        probe_id = f"{product_id}_probe_{gene}"

        # Simulate population-specific coverage rates
        base_rate = rng.uniform(0.5, 1.0)
        pop_cov = {}
        for pop in POPULATIONS:
            noise = rng.gauss(0, 0.08)
            pop_cov[pop] = max(0.0, min(1.0, base_rate + noise))

        rates = list(pop_cov.values())
        disparity = max(rates) - min(rates)
        overall = sum(rates) / len(rates)
        significant = 1 if disparity > 0.15 else 0

        records.append({
            "chrom": chrom,
            "chromStart": start,
            "chromEnd": end,
            "name": probe_id,
            "score": min(int(round(overall * 1000)), 1000),
            "product_id": product_id,
            "technology": tech,
            "coverage_AFR": round(pop_cov["AFR"], 4),
            "coverage_EAS": round(pop_cov["EAS"], 4),
            "coverage_NFE": round(pop_cov["NFE"], 4),
            "coverage_SAS": round(pop_cov["SAS"], 4),
            "coverage_AMR": round(pop_cov["AMR"], 4),
            "disparity_magnitude": round(disparity, 4),
            "significant": significant,
        })

    demo_path = output_dir / "demo.probe_coverage.bed"
    _write_bed(demo_path, records)
    print(f"[DEMO] Wrote {len(records)} simulated probe records to {demo_path}")
    print()

    # Print first few lines as a preview
    with open(demo_path) as fh:
        lines = fh.readlines()
    print("Preview (first 8 lines):")
    for line in lines[:8]:
        print("  " + line.rstrip())


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate GA4GH BED v1.0 files with population-stratified coverage "
            "metrics from ProbeProjection pipeline output."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  %(prog)s --results-dir ../results --output-dir ../distribution/bed\n"
            "  %(prog)s --technologies FISH,NGS\n"
            "  %(prog)s --demo\n"
        ),
    )
    parser.add_argument(
        "--results-dir",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parent.parent / "results",
        help="Path to pipeline results directory (default: ../results/)",
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parent.parent / "distribution" / "bed",
        help="Path for BED output files (default: ../distribution/bed/)",
    )
    parser.add_argument(
        "--technologies",
        type=str,
        default=None,
        help="Comma-separated technology filter (default: all)",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Generate a small example BED file with 20 simulated probes",
    )
    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    tech_filter = None
    if args.technologies:
        tech_filter = [t.strip() for t in args.technologies.split(",")]

    if args.demo:
        _generate_demo(args.output_dir)
        return

    process_results(args.results_dir, args.output_dir, tech_filter)


if __name__ == "__main__":
    main()

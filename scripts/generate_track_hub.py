#!/usr/bin/env python3
"""
generate_track_hub.py

Generate UCSC Track Hub configuration files from ProbeProjection BED files.

Produces:
    hub.txt         -- hub descriptor
    genomes.txt     -- genome assembly reference
    hg38/trackDb.txt -- track database with per-technology subtracks
    description.html -- human-readable data description

Usage:
    python generate_track_hub.py --bed-dir ../distribution/bed
    python generate_track_hub.py --hub-url https://hgwdev.gi.ucsc.edu/~user/ProbeProjection
"""

import argparse
import pathlib
import sys

# ---------------------------------------------------------------------------
# Technology colour palette (R,G,B strings for UCSC trackDb)
# ---------------------------------------------------------------------------
TECHNOLOGY_COLOURS = {
    "FISH":       "70,130,180",   # steel blue
    "NGS":        "60,150,80",    # forest green
    "MICROARRAY": "180,100,50",   # burnt sienna
    "RTPCR":      "140,80,160",   # muted purple
    "SEQUENCING": "60,150,80",    # forest green (alias)
    "MLPA":       "200,140,60",   # gold
    "PCR":        "100,100,180",  # slate blue
    "IHC":        "160,80,80",    # muted red
    "ISH":        "70,130,180",   # steel blue (alias)
}

DEFAULT_COLOUR = "100,100,100"   # grey for unlisted technologies


# ---------------------------------------------------------------------------
# Hub file generators
# ---------------------------------------------------------------------------

def _write_hub_txt(output_dir: pathlib.Path,
                   hub_name: str,
                   hub_url: str) -> pathlib.Path:
    """Write hub.txt descriptor."""
    path = output_dir / "hub.txt"
    lines = [
        f"hub {hub_name}",
        "shortLabel Probe Projection Coverage",
        "longLabel Population-Stratified Clinical Assay Probe Coverage via HPRC Pangenome",
        "genomesFile genomes.txt",
        "email probe-projection@placeholder.org",
        "descriptionUrl description.html",
        "",
    ]
    path.write_text("\n".join(lines))
    return path


def _write_genomes_txt(output_dir: pathlib.Path) -> pathlib.Path:
    """Write genomes.txt pointing to GRCh38 (hg38)."""
    path = output_dir / "genomes.txt"
    lines = [
        "genome hg38",
        "trackDb hg38/trackDb.txt",
        "",
    ]
    path.write_text("\n".join(lines))
    return path


def _discover_bed_files(bed_dir: pathlib.Path) -> list[pathlib.Path]:
    """Return sorted list of *.probe_coverage.bed files in bed_dir."""
    if not bed_dir.is_dir():
        return []
    beds = sorted(bed_dir.glob("*.probe_coverage.bed"))
    # Exclude the combined file for per-technology subtracks; it will be
    # added as its own track.
    return beds


def _technology_from_bed(bed_path: pathlib.Path) -> str:
    """Derive technology name from BED filename."""
    # e.g.  fish.probe_coverage.bed  ->  FISH
    return bed_path.name.split(".")[0].upper()


def _count_bed_records(bed_path: pathlib.Path) -> int:
    """Count non-header, non-empty lines in a BED file."""
    count = 0
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            count += 1
    return count


def _write_track_db(output_dir: pathlib.Path,
                    bed_files: list[pathlib.Path],
                    hub_name: str,
                    hub_url: str) -> pathlib.Path:
    """
    Write hg38/trackDb.txt with a composite track and per-technology
    subtracks.
    """
    hg38_dir = output_dir / "hg38"
    hg38_dir.mkdir(parents=True, exist_ok=True)
    path = hg38_dir / "trackDb.txt"

    lines: list[str] = []

    # -- Composite track --------------------------------------------------
    lines.extend([
        "track probeProjection",
        "compositeTrack on",
        "shortLabel Probe Coverage",
        "longLabel ProbeProjection: Population-Stratified Clinical Assay Probe Coverage",
        "type bigBed 5 +",
        "visibility dense",
        "subGroup1 technology Technology \\",
    ])

    # Build subGroup1 enum from discovered technologies
    per_tech = [f for f in bed_files
                if _technology_from_bed(f) != "ALL_TECHNOLOGIES"]
    tech_names = [_technology_from_bed(f) for f in per_tech]
    all_present = any(_technology_from_bed(f) == "ALL_TECHNOLOGIES"
                      for f in bed_files)

    subgroup_entries = []
    for tn in tech_names:
        tag = tn.lower()
        subgroup_entries.append(f"    {tag}={tn}")
    if all_present:
        subgroup_entries.append("    all=All_Technologies")
    lines.append(" \\\n".join(subgroup_entries))
    lines.append("")

    # -- Per-technology subtracks -----------------------------------------
    for bed_path in bed_files:
        tech = _technology_from_bed(bed_path)
        tag = tech.lower()
        colour = TECHNOLOGY_COLOURS.get(tech, DEFAULT_COLOUR)
        n_records = _count_bed_records(bed_path)

        is_combined = (tech == "ALL_TECHNOLOGIES")
        short = "All Technologies" if is_combined else tech
        long_label = (
            f"Combined probe coverage across all technologies ({n_records} probes)"
            if is_combined else
            f"{tech} probe coverage ({n_records} probes)"
        )
        subgroup_val = "all" if is_combined else tag
        # The bigDataUrl assumes BED files have been converted to bigBed
        # and placed alongside the trackDb.
        bigbed_name = bed_path.stem + ".bb"

        lines.extend([
            f"    track {tag}ProbeProjection",
            f"    parent probeProjection on",
            f"    shortLabel {short}",
            f"    longLabel {long_label}",
            f"    type bigBed 5 +",
            f"    bigDataUrl {bigbed_name}",
            f"    subGroups technology={subgroup_val}",
            f"    color {colour}",
            f"    visibility dense",
            f"    priority {10 if is_combined else 20}",
            "",
        ])

    path.write_text("\n".join(lines))
    return path


def _write_description_html(output_dir: pathlib.Path,
                            hub_name: str,
                            bed_files: list[pathlib.Path]) -> pathlib.Path:
    """Write description.html with a brief data description."""
    path = output_dir / "description.html"

    tech_names = sorted({
        _technology_from_bed(f) for f in bed_files
        if _technology_from_bed(f) != "ALL_TECHNOLOGIES"
    })

    total = sum(_count_bed_records(f) for f in bed_files
                if _technology_from_bed(f) != "ALL_TECHNOLOGIES")

    tech_list_html = "\n".join(f"  <li>{t}</li>" for t in tech_names)

    html = f"""\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>{hub_name} Track Hub</title>
</head>
<body>
<h2>{hub_name}: Population-Stratified Clinical Assay Probe Coverage</h2>

<p>
This track hub displays probe-level coverage metrics derived from projecting
clinical assay probes onto the Human Pangenome Reference Consortium (HPRC)
pangenome graph.  Coverage rates are stratified by five continental
populations (AFR, EAS, NFE, SAS, AMR) to identify probes with
population-dependent performance.
</p>

<h3>Data Summary</h3>
<ul>
  <li>Total probe records: {total}</li>
  <li>Reference assembly: GRCh38 (hg38)</li>
  <li>Technologies included:</li>
  <ul>
{tech_list_html}
  </ul>
</ul>

<h3>BED Column Descriptions</h3>
<table border="1" cellpadding="4" cellspacing="0">
  <tr><th>Column</th><th>Description</th></tr>
  <tr><td>chrom</td><td>Chromosome</td></tr>
  <tr><td>chromStart</td><td>0-based start position</td></tr>
  <tr><td>chromEnd</td><td>End position (half-open)</td></tr>
  <tr><td>name</td><td>Probe identifier</td></tr>
  <tr><td>score</td><td>Overall coverage rate &times; 1000 (0-1000)</td></tr>
  <tr><td>product_id</td><td>Assay product identifier</td></tr>
  <tr><td>technology</td><td>Assay technology (FISH, NGS, etc.)</td></tr>
  <tr><td>coverage_AFR</td><td>Coverage rate in African population</td></tr>
  <tr><td>coverage_EAS</td><td>Coverage rate in East Asian population</td></tr>
  <tr><td>coverage_NFE</td><td>Coverage rate in Non-Finnish European population</td></tr>
  <tr><td>coverage_SAS</td><td>Coverage rate in South Asian population</td></tr>
  <tr><td>coverage_AMR</td><td>Coverage rate in American (admixed) population</td></tr>
  <tr><td>disparity_magnitude</td><td>Max inter-population coverage difference</td></tr>
  <tr><td>significant</td><td>1 if disparity is statistically significant, 0 otherwise</td></tr>
</table>

<h3>Methodology</h3>
<p>
Probes are projected onto HPRC haplotype paths via pangenome graph coordinate
liftover.  Per-population coverage rates reflect the fraction of haplotypes
in each population where the probe maps successfully.  Disparity significance
is assessed via permutation testing with correction for multiple comparisons.
</p>

<p>
<strong>Research Use Only.</strong>  Population-level conclusions should not be
extrapolated beyond the HPRC cohort.  Performance variation is attributed to
assay design, not to patient ancestry.
</p>

<h3>Source</h3>
<p>
Generated by the <a href="https://github.com/lauferva/ProbeProjection">ProbeProjection</a> pipeline.
</p>

</body>
</html>
"""
    path.write_text(html)
    return path


# ---------------------------------------------------------------------------
# Post-generation instructions
# ---------------------------------------------------------------------------

def _print_upload_instructions(output_dir: pathlib.Path,
                               bed_files: list[pathlib.Path],
                               hub_url: str) -> None:
    """Print instructions for converting BED to bigBed and uploading."""
    print()
    print("=== Next Steps ===")
    print()
    print("1. Convert BED files to bigBed format using UCSC bedToBigBed:")
    print()
    print("   # Download the hg38 chrom sizes file:")
    print("   wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes")
    print()
    print("   # Create an AutoSQL file for the extra columns (probe_coverage.as):")
    print("   # (A template has been written to the output directory.)")
    print()

    for bed_path in bed_files:
        bb_name = bed_path.stem + ".bb"
        print(f"   bedToBigBed -type=bed5+9 -as=probe_coverage.as \\")
        print(f"       {bed_path} hg38.chrom.sizes {output_dir / 'hg38' / bb_name}")
        print()

    print("2. Upload the track hub directory to a public web server or")
    print("   UCSC HubSpace (https://hgwdev.gi.ucsc.edu/~user/hubspace/).")
    print()
    print("3. Load the hub in the UCSC Genome Browser:")
    print(f"   My Data > Track Hubs > My Hubs > {hub_url}/hub.txt")
    print()


def _write_autosql(output_dir: pathlib.Path) -> pathlib.Path:
    """Write an AutoSQL definition file for bedToBigBed."""
    path = output_dir / "probe_coverage.as"
    content = """\
table probeCoverage
"Population-stratified clinical assay probe coverage from ProbeProjection"
    (
    string chrom;              "Chromosome"
    uint chromStart;           "Start position (0-based)"
    uint chromEnd;             "End position (half-open)"
    string name;               "Probe identifier"
    uint score;                "Overall coverage rate * 1000 (0-1000)"
    string product_id;         "Assay product identifier"
    string technology;         "Assay technology"
    float coverage_AFR;        "Coverage rate in African population"
    float coverage_EAS;        "Coverage rate in East Asian population"
    float coverage_NFE;        "Coverage rate in Non-Finnish European population"
    float coverage_SAS;        "Coverage rate in South Asian population"
    float coverage_AMR;        "Coverage rate in American (admixed) population"
    float disparity_magnitude; "Max inter-population coverage difference"
    uint significant;          "1 if disparity is statistically significant"
    )
"""
    path.write_text(content)
    return path


# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------

def generate_track_hub(bed_dir: pathlib.Path,
                       output_dir: pathlib.Path,
                       hub_name: str,
                       hub_url: str) -> None:
    """Generate UCSC Track Hub files from BED data."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Discover BED files
    bed_files = _discover_bed_files(bed_dir)
    if not bed_files:
        print(f"[INFO] No *.probe_coverage.bed files found in {bed_dir}")
        print("       Run generate_bed_files.py first to produce BED files,")
        print("       or run generate_bed_files.py --demo for example data.")
        print()
        print("       Generating track hub skeleton with no subtracks.")
        print()

    # Write hub configuration files
    hub_path = _write_hub_txt(output_dir, hub_name, hub_url)
    print(f"[OK] {hub_path}")

    genomes_path = _write_genomes_txt(output_dir)
    print(f"[OK] {genomes_path}")

    trackdb_path = _write_track_db(output_dir, bed_files, hub_name, hub_url)
    print(f"[OK] {trackdb_path}")

    desc_path = _write_description_html(output_dir, hub_name, bed_files)
    print(f"[OK] {desc_path}")

    autosql_path = _write_autosql(output_dir)
    print(f"[OK] {autosql_path}")

    # Summary
    print()
    tech_names = sorted({
        _technology_from_bed(f) for f in bed_files
        if _technology_from_bed(f) != "ALL_TECHNOLOGIES"
    })
    print(f"Track hub generated with {len(bed_files)} BED track(s)")
    if tech_names:
        print(f"Technologies: {', '.join(tech_names)}")

    _print_upload_instructions(output_dir, bed_files, hub_url)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate UCSC Track Hub configuration files from "
            "ProbeProjection BED files."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  %(prog)s --bed-dir ../distribution/bed\n"
            "  %(prog)s --hub-url https://example.com/hubs/ProbeProjection\n"
        ),
    )
    parser.add_argument(
        "--bed-dir",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parent.parent / "distribution" / "bed",
        help="Path to BED files (default: ../distribution/bed/)",
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=(
            pathlib.Path(__file__).resolve().parent.parent / "distribution" / "track_hub"
        ),
        help="Path for track hub output (default: ../distribution/track_hub/)",
    )
    parser.add_argument(
        "--hub-name",
        type=str,
        default="ProbeProjection",
        help="Hub name identifier (default: ProbeProjection)",
    )
    parser.add_argument(
        "--hub-url",
        type=str,
        default="https://your-server.example.com/ProbeProjection",
        help="Base URL where hub will be hosted (default: placeholder)",
    )
    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()
    generate_track_hub(args.bed_dir, args.output_dir, args.hub_name, args.hub_url)


if __name__ == "__main__":
    main()

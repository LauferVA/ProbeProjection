# ProbeProjection

Pangenome-based evaluation of clinical assay probe coverage across human haplotype diversity.

## What This Is

ProbeProjection is a dataset and analysis repository that uses the [HPRC pangenome](https://humanpangenome.org/) to evaluate whether commercial clinical laboratory assays — FISH, CMA, MLPA, NGS panels, OGM, and others — provide equitable probe coverage across global human populations.

Clinical assay probes are designed against single-reference genomes (GRCh37/GRCh38), which are predominantly European-ancestry. Structural variants, insertions, deletions, and sequence divergence in non-European haplotypes may cause individual probes to fail silently in certain patient populations. This project quantifies that problem across the full breadth of clinical genomic technologies.

## Why It Matters

A FISH probe for a common cancer translocation might map perfectly to all European-ancestry haplotypes in the HPRC pangenome but fail to map in a subset of African-ancestry haplotypes — not because of biology, but because the probe was designed against a reference genome that doesn't represent that patient's structural variation. The patient's test result could be a false negative due to assay design, not disease status.

This has been discussed in theory (reference genome bias is well-documented) but never systematically quantified across the range of technologies used in clinical cytogenomics and molecular diagnostics. This project does that quantification.

## Relationship to Complement

This repository consumes pipeline infrastructure from the [Complement](https://github.com/lauferva/Complement) project, which provides:

- **Profile catalog** (`complement/reference/profiles/`): 1,628 NormalizedProfile JSON files spanning 16 clinical technologies, each containing standardized probe coordinates.
- **Graph manager** (`complement/graph/manager.py`): LRU-cached per-chromosome `.og` pangenome subgraph loading with odgi/bdsg backend selection.
- **ODGI projector** (`complement/graph/projector.py`): Subprocess wrapper around `odgi position` for projecting genomic coordinates from reference paths to haplotype paths.
- **Haplotype coverage analysis** (`complement/analysis/haplotype_coverage.py`): Projects every probe in an assay to every HPRC haplotype and classifies mapping outcomes (ABSENT, DIVERGENT, DUPLICATED, REARRANGED).
- **Ancestry/diversity analysis** (`complement/analysis/ancestry.py`): Population-stratified chi-squared testing to detect ancestry-correlated probe coverage disparities with Bonferroni correction.
- **Orchestration script** (`Complement/scripts/run_probe_projection_corpus.py`): End-to-end pipeline tying all components together.

ProbeProjection is the output repository: it stores run manifests, aggregated results, analysis notebooks, figures, and specifications produced by running the Complement pipeline against the full probe catalog.

## Pipeline Overview

```
1,628 NormalizedProfile JSONs (16 technologies)
    |
    v
Filter to GRCh38-assembly profiles (625 eligible)
    |
    v
Per profile, per probe, per chromosome:
    Load .og pangenome subgraph (24 chromosomes, ~23 GB total)
    |
    v
    For each of ~272 HPRC haplotype paths:
        odgi position: project probe from GRCh38 reference -> haplotype
        Classify: MAPPED | ABSENT | DIVERGENT | DUPLICATED | REARRANGED
    |
    v
Per-profile outputs:
    - HaplotypeCoverageReport: probe-level mapping success across all haplotypes
    - HaplotypeDiversityReport: chi-squared test per probe, stratified by
      superpopulation (AFR, AMR, EAS, NFE, SAS), Bonferroni-corrected
    |
    v
Cross-profile aggregation:
    - Per-technology summary statistics
    - Global disparity metrics
    - Identification of most-affected assays and populations
```

## Technologies Covered

| Technology | Profiles | Typical Probes/Profile | Description |
|-----------|----------|----------------------|-------------|
| FISH | 319 | 1-4 | Fluorescence in situ hybridization (BAC/fosmid, ~100-200 kb) |
| NGS | 184 | 50-5,000 | Next-generation sequencing capture panels |
| HLA | 38 | Gene-level | Human leukocyte antigen typing |
| PGx | 23 | Gene + star alleles | Pharmacogenomic testing |
| Methylation | 16 | CpG sites | Methylation-specific assays |
| RT-PCR | 11 | Fusion/expression | Real-time PCR for gene fusions |
| CE/FA | 9 | STR markers | Capillary electrophoresis / fragment analysis |
| ddPCR | 6 | Single targets | Digital droplet PCR |
| MLPA | 6 | 40-60 | Multiplex ligation-dependent probe amplification |
| Southern Blot | 6 | Restriction fragments | Southern blot hybridization |
| Long-Read | 4 | SV targets | Nanopore/PacBio structural variant panels |
| OGM | 2 | ~500K+ | Optical genome mapping (DLE-1 labels) |
| CMA | 1 | ~6.85M | Chromosomal microarray (oligonucleotide probes) |

## HPRC Population Representation

| Superpopulation | Samples | Haplotypes | Description |
|----------------|---------|------------|-------------|
| AFR | 53 | 106 | African/African American |
| NFE | 35 | 70 | Non-Finnish European |
| EAS | 24 | 48 | East Asian |
| SAS | 19 | 38 | South Asian |
| AMR | 5 | 10 | Admixed American |
| **Total** | **136** | **272** | |

## Statistical Framework

- **Test**: Chi-squared test of independence (populations x {mapped, unmapped}) per probe
- **Multiple testing correction**: Bonferroni (p_corrected = p_raw x n_probes)
- **Significance gates**: Corrected p < 0.05 AND disparity magnitude > 20%
- **Minimum sample**: 10 haplotypes per superpopulation (smaller groups marked INSUFFICIENT_DATA)
- **Framing**: Disparities attributed to assay design limitations, not patient ancestry
- **Regulatory status**: All outputs are Research Use Only (RUO)

## Repository Structure

```
ProbeProjection/
    README.md               # This file
    results/
        manifest.json       # Run parameters and provenance metadata
    scripts/                # Analysis and visualization scripts (planned)
    specs/                  # Figure and tool specifications (planned)
    analysis/               # Research outputs and literature review (planned)
```

## Current Status

**Blocked on graph data.** The pipeline code is fully implemented and validated (dry run passes, profile loading works, graph loading works). However, the per-chromosome `.og` files currently contain only reference assembly paths (GRCh38, GRCh37, GRCh36, CHM13) and lack the ~272 HPRC sample haplotype paths needed for projection. Every `odgi position` call targeting a haplotype path (e.g., `HG00272#1#chr17`) fails with "path not found."

The graph topology is correct — the `.og` files encode pangenome variation across ~2M+ nodes per chromosome. The sample-specific path annotations were not included when the files were built. The fix is to rebuild the `.og` files from the HPRC `.gbz` (which contains all haplotype paths) or to download pre-built per-chromosome graphs from the HPRC data release.

See `UNBLOCK_STATUS.md` for full technical details on the blocker and three resolution options.

## License

To be determined.

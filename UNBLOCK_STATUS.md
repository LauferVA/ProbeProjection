# Probe Projection Dataset 1 — Status

**Last updated:** 2026-03-30
**Status:** BLOCKED — graph files lack HPRC haplotype paths
**Preprint target:** September 2026

---

## Current Blocker

The 24 per-chromosome `.og` graph files contain ONLY reference assembly paths (GRCh38, GRCh37, GRCh36, CHM13). Zero HPRC sample haplotype paths. Every `odgi position` call targeting a haplotype (e.g., `HG00272#1#chr17`) fails with "path not found."

The graph topology is correct (~2M+ nodes/edges per chromosome). The sample-specific path annotations were not included when the `.og` files were built.

### Fix Options

**Option A (recommended): Install vg + extract from GBZ**
```bash
conda install -n complement -c bioconda vg
vg paths --list hprc-v2.0-mc-grch38.gbz | head -30
# Extract per-chromosome subgraphs with all paths
```

**Option B: Download pre-built HPRC per-chromosome .og files**
- Source: `https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze2/minigraph-cactus/`
- Size: ~100-200GB total

**Option C: Build from GFA**
```bash
odgi build -g chrN.gfa -o chrN.og
# GFA at: ~/genomics_reference_files/hprc/release2/hprc-v2.0-mc-grch38.gfa.gz
```

---

## What Works

| Capability | Status |
|-----------|--------|
| Profile loading (1,563 profiles, 16 technologies) | WORKS |
| Graph loading (24 chromosomes, 23GB) | WORKS |
| HPRC cohort + ancestry registry (136 samples, 5 superpopulations) | WORKS |
| Cross-assembly liftover (GRCh38 <-> GRCh37) | WORKS |
| Dry run (625 eligible GRCh38 profiles) | WORKS |
| Haplotype projection | BLOCKED (no sample paths in .og files) |
| Diversity analysis | BLOCKED (depends on projection) |
| Figures 1, 2 (Panel A), 5 | DONE (data-independent) |
| Figures 2 (B, C), 3, 4 | BLOCKED (require projection results) |
| BED/track hub generation scripts | READY (await pipeline output) |
| Power analysis | DONE |

---

## What Was Completed This Session

### Research (4 agents, ~6 hours total)
1. **Causal chain analysis** — Pipeline is a structural variant screen, not hybridization simulator. FISH detection >90%, CMA <10%, MLPA <5%. Report: `analysis/causal_chain_evidence_review.md`
2. **Expected magnitude** — >95% of probes map fine. Value in specific smoking guns (17q21.31, 8p23.1, CYP2D6, HLA) + reassuring baseline. HLA/PGx: 10-30% disparity.
3. **Competitive landscape** — Dong et al. (Cell Reports 2025) validated for methylation arrays. 16-tech scope unoccupied. Scoop risk MEDIUM-HIGH. Report: `analysis/competitive_landscape_2026_03_30.md`
4. **Regulatory reality** — No FDA/CAP requirement for ancestry-stratified probe validation. ACMG 2023 advisory only. LDT rule vacated. Distribution: BED + UCSC track hub.

### Specifications
- `specs/SPEC_PRESENTATION.md` — Full presentation plan with ranked contributions, publication strategy, explicit non-claims
- `specs/SPEC_FIGURES.md` — Per-figure specs for 5 figures and 3 tables

### Scripts and Figures
- 3 figure generation scripts producing Figs 1, 2, 5 (PNG + PDF)
- 4 utility scripts (BED generation, track hub config, summary stats, power analysis)
- Power analysis output: `results/summaries/power_analysis.{json,tsv}`

### Repository
- GitHub: https://github.com/LauferVA/ProbeProjection
- 13 granular commits on `main` (single clean branch)

---

## Once Unblocked: Execution Plan

1. Dry run with new graphs to confirm haplotype paths found
2. **FISH first** (319 profiles) — highest priority, strongest causal chain
3. **HLA + PGx** (61 profiles) — expected highest disparity rates
4. Remaining technologies — `--parallel 4`
5. CMA (6.85M probes) — `--parallel 2`, `--max-loaded-chromosomes 2`
6. Regenerate Figures 2 (B, C), 3, 4 with real data
7. Generate BED files and UCSC track hub
8. Manuscript draft -> preprint September 2026

---

## Graph Data Inventory

| Resource | Location | Haplotype Paths? |
|----------|----------|------------------|
| `.og` (ref liftover) | `~/genomics_reference_files/complement/graph/hprc-v2.0-mc-grch38-grch37-grch36/og/` (24 chr, 23GB) | NO |
| `.og` (partial) | `~/genomics_reference_files/complement/graph/hprc-v2.0-mc-grch38/` (chr1-8, 21, 22) | NO |
| `.pg` | `~/genomics_reference_files/complement/graph/hprc-v2.0-mc-grch38/pg/` (24 chr, 6.6GB) | NO |
| `.gbz` (full pangenome) | `~/genomics_reference_files/hprc/release2/hprc-v2.0-mc-grch38.gbz` (5GB) | LIKELY YES |
| `.gbz` (v1.1 grch38) | `~/genomics_reference_files/complement/graph/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz` | LIKELY YES |

---

## Pipeline Code Reference

- **Orchestration:** `~/software/Complement/scripts/run_probe_projection_corpus.py`
- **Haplotype coverage:** `complement/analysis/haplotype_coverage.py`
- **Diversity analysis:** `complement/analysis/ancestry.py`
- **Graph manager:** `complement/graph/manager.py`
- **ODGI projector:** `complement/graph/projector.py`
- **Profile models:** `complement/reference/profiles/models.py`
- **Specs:** `Complement/.spec/DATASET_01_EXECUTION_NOTES.md`, `Complement/.spec/SPEC_CROSS_GOAL_DATASET_01_PROBE_PROJECTION_CORPUS.md`

## Key Technical Detail

Pipeline expects PanSN-format haplotype paths: `{sample}#{haplotype}#{chromosome}` (e.g., `HG00097#1#chr9`). Reference paths: `GRCh38#0#chr1`. Expected ~200 paths per autosome.

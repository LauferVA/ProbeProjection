# Probe Projection Dataset 1 — Session Summary

**Date:** 2026-03-30
**Status:** BLOCKED — graph files lack HPRC haplotype paths

---

## Session Work Completed

### 1. Environment Validated
- `odgi` v0.9.4 confirmed installed in `complement` conda env
- `bdsg` Python bindings work (PackedGraph imports clean)
- Profile data restored: 1,563 profiles across 16 technologies
- 24 `.og` chromosome graphs present (23GB total)
- 136 HPRC samples loaded, 5 superpopulations

### 2. Dry Run Passed
- **625 eligible GRCh38 profiles** across 13 technologies
- Technology breakdown: FISH=319, NGS=184, HLA=38, PGX=23, METHYLATION=16, RTPCR=11, CE_FA=9, DDPCR=6, MLPA=6, SOUTHERN_BLOT=6, LONG_READ=4, OGM=2, CMA=1
- 6,854,584 total probes (dominated by CMA_0004 at ~6.85M)
- Pipeline loads all resources, filters profiles, and exits cleanly

### 3. FISH Projection Attempted and Stopped
- Launched 319 FISH profiles with `--parallel 4`
- Pipeline entered Phase 3 (processing) but produced zero completed profiles after ~6 minutes
- Every `odgi position` call failed: `ref path HG00272#1#chr17 not found in graph`
- 381+ path-not-found errors before pipeline was killed

### 4. Root Cause Identified

**The .og graph files contain ONLY reference assembly paths. Zero HPRC sample haplotype paths.**

Verified across multiple chromosomes:
- `chr1.og`: 106 paths — all CHM13, GRCh36, GRCh37, GRCh38 fragments. Zero HG00*/NA* paths.
- `chr17.og`: 21 paths — same pattern.
- `chr22.og`: 57 paths — same pattern.
- `chr7.og`, `chr12.og`, `chrX.og` — same.

The graph *topology* is correct (2.1M nodes, 3M edges for chr22 alone), encoding pangenome variation. But individual sample path annotations (e.g., `HG00272#1#chr17`, `NA20799#2#chr12`) were not included when the `.og` files were built.

**Impact:** `project_probes_to_haplotypes()` calls `odgi position` with target haplotype paths that don't exist. Every lookup fails. Coverage reports show 100% ABSENT, diversity reports show 0% disparity. The pipeline produces structurally valid but scientifically empty results.

---

## What Works vs. What Doesn't

| Capability | Status | Why |
|-----------|--------|-----|
| Profile loading (1,563 profiles) | WORKS | NormalizedProfile JSONs restored from stale dir |
| Graph loading (24 chromosomes) | WORKS | GraphManager finds .og files, odgi available |
| HPRC cohort + ancestry registry | WORKS | 136 samples, 5 superpopulations |
| Cross-assembly liftover (GRCh38 ↔ GRCh37) | WORKS | Reference paths present in graphs |
| Haplotype projection (GRCh38 → HG00272#1#chrN) | FAILS | No sample paths in .og files |
| Diversity analysis (chi-squared per population) | FAILS | Depends on haplotype projection |

---

## Graph Data Inventory

| Resource | Location | Haplotype Paths? |
|----------|----------|------------------|
| `.og` (ref liftover) | `~/genomics_reference_files/complement/graph/hprc-v2.0-mc-grch38-grch37-grch36/og/` (24 chr, 23GB) | NO — reference only |
| `.og` (partial) | `~/genomics_reference_files/complement/graph/hprc-v2.0-mc-grch38/` (chr1-8, 21, 22) | NO — reference only |
| `.pg` | `~/genomics_reference_files/complement/graph/hprc-v2.0-mc-grch38/pg/` (24 chr, 6.6GB) | NO — CHM13 only |
| `.gbz` (full pangenome) | `~/genomics_reference_files/hprc/release2/hprc-v2.0-mc-grch38.gbz` (5GB, symlinked) | LIKELY YES — needs vg |
| `.gbz` (v1.1 grch38) | `~/genomics_reference_files/complement/graph/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz` | LIKELY YES — needs vg |
| `.gbz` (v1.1 chm13) | `~/genomics_reference_files/complement/graph/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.gbz` | LIKELY YES — needs vg |

---

## Fix: Rebuild .og Files with Haplotype Paths (in Complement dir)

The `.og` files need to be rebuilt from the GBZ to include all ~200 paths per chromosome (94 diploid samples × 2 haplotypes + references).

### Option A: Install vg + extract per-chromosome graphs (recommended)
```bash
conda install -n complement -c bioconda vg

# Verify GBZ has haplotype paths
conda activate complement
vg paths --list hprc-v2.0-mc-grch38.gbz | head -30
# Expect: GRCh38#0#chr1, HG00097#1#chr1, HG00097#2#chr1, ...

# Extract per-chromosome subgraphs with all paths
# (exact commands TBD — depends on vg version and graph structure)
```

### Option B: Download pre-built per-chromosome graphs from HPRC
HPRC distributes per-chromosome `.og` files with all haplotype paths:
- Source: `https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze2/minigraph-cactus/`
- Size: ~100-200GB total for all chromosomes with full paths
- These would replace the current reference-only `.og` files

### Option C: Build from GFA
If a per-chromosome `.gfa` with all paths exists (e.g., extracted from the `.gfa.gz`):
```bash
# For each chromosome
odgi build -g chrN.gfa -o chrN.og
```
Note: The compressed GFA is at `~/genomics_reference_files/hprc/release2/hprc-v2.0-mc-grch38.gfa.gz`

---

## Once Graphs Are Rebuilt

Resume the pipeline with the execution plan from `Complement/.spec/DATASET_01_EXECUTION_NOTES.md`:

1. **Dry run** with new graphs to confirm haplotype paths are found
2. **FISH** (319 profiles, ~1-4 probes each) — fast validation
3. **Batch** (FISH + MLPA + ddPCR + RT-PCR + PGx + HLA) — ~400 profiles, `--parallel 4`
4. **NGS** (184 profiles) — `--parallel 2`
5. **Large-scale** (CMA + methylation + OGM) — CMA has 6.85M probes, `--parallel 2`, `--max-loaded-chromosomes 2`

Output goes to `~/software/ProbeProjection/results/`.

---

## Pipeline Code Reference

- **Orchestration script:** `~/software/Complement/scripts/run_probe_projection_corpus.py` (1,147 lines)
- **Haplotype coverage:** `complement/analysis/haplotype_coverage.py` (1,311 lines)
- **Diversity analysis:** `complement/analysis/ancestry.py` (1,598 lines)
- **GraphManager:** `complement/graph/manager.py` (994 lines)
- **ODGI projector:** `complement/graph/projector.py`
- **Profile models:** `complement/reference/profiles/models.py` (504 lines)
- **Execution notes:** `Complement/.spec/DATASET_01_EXECUTION_NOTES.md`
- **Full spec:** `Complement/.spec/SPEC_CROSS_GOAL_DATASET_01_PROBE_PROJECTION_CORPUS.md`

---

## Key Technical Details for Graph Rebuild

The pipeline expects PanSN-format haplotype path names in the `.og` files:
- Format: `{sample}#{haplotype}#{chromosome}` (e.g., `HG00097#1#chr9`, `NA20799#2#chr17`)
- Reference paths: `GRCh38#0#chr1`, `CHM13#0#chr1`, etc.
- Expected: ~200 paths per autosome (94 haplotypes + reference fragments)
- The `odgi position` command translates coordinates: `GRCh38#0#chrN` → `HG00097#1#chrN`

The current `.og` files have the correct graph topology (nodes + edges) but are missing the sample-specific path annotations. Whether the rebuild needs to reconstruct the full graph or can inject paths into existing `.og` files depends on the tooling.

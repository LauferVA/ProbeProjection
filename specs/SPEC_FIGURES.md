# Figure and Table Specifications: ProbeProjection Corpus

**Version:** 1.0
**Date:** 2026-03-30
**Status:** Pre-execution (data-dependent details marked with [DATA])

---

## Figure 1: Pipeline Overview

**Purpose:** Orient the reader to the method.

**Type:** Schematic / flow diagram.

**Content:**
```
[Profile Catalog]          [HPRC Pangenome]        [Ancestry Registry]
 1,628 profiles              24 per-chr .og          136 samples
 16 technologies             ~200M nodes/edges       5 superpopulations
 625 eligible (GRCh38)       272 haplotype paths     AFR/EAS/NFE/SAS/AMR
       |                           |                        |
       v                           v                        v
 Probe Extraction -----> Coordinate Projection -----> Population Stratification
 (chr, start, end)       (odgi position per           (chi-squared / Fisher
  per profile)            probe x haplotype)           exact per probe)
       |                           |                        |
       v                           v                        v
 Per-profile results:    Classification:             Disparity metrics:
 coverage.json           MAPPED / ABSENT             p-value, magnitude,
 diversity.json          (per haplotype)              affected populations
       |
       v
 Cross-profile aggregation -> summary.json + BED files + UCSC track hub
```

**Design notes:**
- Clean, minimal, use muted colors (no neon)
- Technology icons or labels at input stage showing scope (FISH, CMA, MLPA, NGS, OGM, PGx, HLA, ...)
- Single-page width for journal column format
- Software: matplotlib + manual annotation, or Inkscape/Illustrator for final

**Data required:** None (schematic only).

---

## Figure 2: Technology Landscape (MAIN FIGURE)

**Purpose:** The central result figure. Shows scope, baseline, and exceptions by technology.

### Panel A: Catalog Scope

**Type:** Horizontal bar chart or Cleveland dot plot.

**X-axis:** Number of probes (log scale) or number of profiles.
**Y-axis:** Technology (FISH, CMA, NGS, HLA, PGx, MLPA, OGM, ...), ordered by probe count or profile count.
**Annotations:** Profile count label at bar end.

**Data required:** Profile and probe counts per technology. Available from ProfileCatalog metadata (no pipeline run needed).

### Panel B: Projection Results Summary

**Type:** Stacked bar chart or diverging bar chart.

**X-axis:** Fraction of probes (0-100%).
**Y-axis:** Technology (same order as Panel A).
**Segments:**
- Green: probes with universal coverage (all haplotypes mapped)
- Yellow: probes with any haplotype-level failure (not population-specific)
- Red: probes with population-specific disparity (significant after correction AND magnitude > 20%)

**Data required:** [DATA] Per-technology aggregate projection results from summary.json. Requires pipeline execution.

### Panel C: Population Disparity Profile

**Type:** Grouped dot plot with 95% confidence intervals.

**X-axis:** Mean probe coverage rate (0-1).
**Y-axis:** Superpopulation (AFR, EAS, NFE, SAS, AMR).
**Groups:** One group per technology (showing only technologies with detectable disparities: likely FISH, HLA, PGx, possibly NGS).
**AMR:** Displayed with hollow marker and wider CI to visually flag underpowered status.

**Data required:** [DATA] Per-superpopulation mean coverage rates from diversity reports. Requires pipeline execution.

**Design notes:**
- Three panels side by side or stacked vertically
- Consistent technology ordering across panels
- Color palette: colorblind-safe (viridis or similar)
- Journal target: full-page width

---

## Figure 3: Smoking Gun Loci (NARRATIVE ANCHOR)

**Purpose:** Concrete examples of population-specific probe coverage gaps at known structural variant loci. This is the figure that makes the abstract framework tangible.

**Layout:** 2x2 grid of panels. Each panel shows one locus.

### Per-Panel Structure (repeated for each locus)

**Top strip:** Chromosome ideogram with probe location marked (red bar).

**Middle:** Simplified pangenome subgraph at the locus.
- Reference path (GRCh38) shown as horizontal line
- Alternate haplotype paths shown as branches/bubbles
- Structural variant (inversion, deletion) indicated by path divergence
- Color-coded by superpopulation where possible
- This panel requires extracting subgraph topology — NOT currently in pipeline outputs. Requires: `odgi viz` or `odgi layout` to produce a 1D visualization of the local graph, or a simplified hand-drawn schematic based on known SV structure.

**Bottom:** Bar chart of projection success rate by superpopulation.
- X-axis: superpopulation (AFR, EAS, NFE, SAS, AMR)
- Y-axis: fraction of haplotypes where probe maps successfully (0-1)
- Error bars: 95% Wilson confidence intervals
- AMR bar: hatched or outlined to flag low power

### Candidate Loci

**Panel 3A: 17q21.31 (KANSL1/MAPT)**
- SV: H2 inversion, 1.08-1.49 Mbp
- Population: ~20% EUR, ~0% AFR/EAS
- Clinical assay: FISH for Koolen-De Vries syndrome
- Expected pattern: EUR haplotypes show reduced projection success; AFR/EAS show near-100%
- Note: European-ancestry haplotypes are the affected group here

**Panel 3B: 8p23.1 (Defensin cluster)**
- SV: ~4.5 Mb inversion
- Population: ~26% EUR, clinal distribution
- Clinical assay: FISH for 8p deletion/duplication syndromes
- Expected pattern: Similar to 17q21.31; inversion carriers show projection disruption

**Panel 3C: 22q11.2 (DiGeorge/VCFS)**
- SV: 50 structural configurations, 11-fold length variation in LCRA
- Population: Population-specific SD expansion patterns
- Clinical assay: FISH (TUPLE1/ARSA), MLPA, NGS cardiac panels
- Expected pattern: Variable by haplotype; probes in LCR flanks most affected

**Panel 3D: CYP2D6 (22q13.2)**
- SV: Whole-gene deletion (*5), duplications, hybrid genes
- Population: Deletion 6% AFR, 14.5% African American, varies by population
- Clinical assay: PGx star allele panels
- Expected pattern: AFR haplotypes show highest gene-level absence rate

**Data required:** [DATA] Per-probe per-haplotype projection results at these specific loci. Requires pipeline execution. Additionally, subgraph visualization requires `odgi viz` or manual schematic — not part of current pipeline output.

**Design notes:**
- Each panel should be self-contained with a 1-line title (e.g., "17q21.31: FISH probe for Koolen-De Vries syndrome")
- The middle graph schematic can be simplified/stylized; does not need to show every node
- Reference published ideograms or genome browser screenshots for the top strip

---

## Figure 4: HLA/PGx Disparity Landscape

**Purpose:** Show the high-disparity zone where population-specific variation is most pronounced.

**Type:** Heatmap or tile plot.

**Rows:** HLA/PGx gene targets (HLA-A, HLA-B, HLA-C, HLA-DRB1, CYP2D6, CYP2C19, CYP2A6, UGT1A1, ...). Ordered by mean disparity magnitude descending.

**Columns:** Superpopulations (AFR, EAS, NFE, SAS). AMR excluded or shown with stippled cells to flag low power.

**Cell color:** Projection coverage rate (0-1). Diverging palette: dark blue (0% coverage) through white (50%) to dark red (100%). Or sequential: white (100%, no gap) to dark red (0%, complete gap).

**Side annotations:**
- Clinical significance tier or pharmacogenomic actionability level
- Number of clinical profiles that include this target

**Data required:** [DATA] Per-gene per-superpopulation projection coverage rates from HLA and PGx diversity reports. Requires pipeline execution.

**Design notes:**
- Compact format suitable for half-page width
- Row count depends on data; expect 20-50 rows for HLA + PGx combined
- Cluster rows by technology (HLA block, PGx block) with a dividing line

---

## Figure 5 (Supplementary): Statistical Power and Limitations

**Purpose:** Transparency about what the data can and cannot support.

### Panel A: Minimum Detectable Effect Size

**Type:** Line plot.

**X-axis:** Number of haplotypes in the smaller comparison group.
**Y-axis:** Minimum detectable absolute disparity (%) at 80% power, alpha=0.05.
**Lines:** One per significance framework (Bonferroni for 4 probes, Bonferroni for 500 probes, FDR q=0.05).
**Vertical markers:** Actual sample sizes (AFR=106, NFE=70, EAS=48, SAS=38, AMR=10).

**Data required:** Power calculations (can be computed without pipeline data).

### Panel B: Technology Detection Capability

**Type:** Matrix/table visualization.

**Rows:** Technologies (FISH, CMA, MLPA, NGS, OGM).
**Columns:** Failure modes (ABSENT/large SV, DIVERGENT/sequence mismatch, DUPLICATED/multi-mapping, REARRANGED/translocation).
**Cells:** Filled circle = detectable by coordinate projection; empty circle = not detectable; half-filled = partially detectable.

**Data required:** None (based on methodology analysis; already established by research agents).

---

## Table 1: Technology Summary

| Column | Description | Source |
|---|---|---|
| Technology | FISH, CMA, MLPA, NGS, OGM, PGx, HLA, ... | ProfileCatalog |
| # Profiles | Count of eligible GRCh38 profiles | ProfileCatalog |
| # Probes (total) | Sum of probes across profiles | ProfileCatalog |
| Probe Size Range | Min-max probe length | ProfileCatalog |
| % Universal Coverage | Fraction of probes mapping to all haplotypes | [DATA] summary.json |
| % Population-Specific Gap | Fraction with significant disparity | [DATA] summary.json |
| Most Affected Population | Superpopulation with lowest mean coverage | [DATA] summary.json |

**Data required:** First two columns available now. Remaining require pipeline execution.

---

## Table 2: Priority Loci with Population-Specific Gaps

| Column | Description | Source |
|---|---|---|
| Locus | Genomic region (e.g., 17q21.31) | [DATA] diversity reports |
| Chromosome | chr17, chr8, chr22, chr6, ... | [DATA] |
| SV Type | Inversion, deletion, duplication, complex | Known biology + [DATA] |
| Population(s) Affected | Superpopulation(s) with lowest coverage | [DATA] |
| SV Population Frequency | Frequency of the structural variant | Published literature |
| Clinical Assay(s) | Product IDs/names targeting this locus | ProfileCatalog |
| Technology | FISH, PGx, HLA, NGS, ... | ProfileCatalog |
| Gap Rate | Fraction of affected-population haplotypes with projection failure | [DATA] |
| Clinical Context | Condition tested for (e.g., Koolen-De Vries, DiGeorge) | Manual curation |
| Recommended Action | Awareness, confirmatory testing, redesign | Manual curation |

**Ordering:** Ranked by clinical importance x disparity magnitude x confidence.

**Expected rows:** 10-30 loci (based on Tier 1 + Tier 2 smoking gun candidates plus data-driven discoveries).

---

## Table 3: HPRC Cohort Composition

| Superpopulation | Code | Samples | Haplotypes | Statistical Note |
|---|---|---|---|---|
| African/African American | AFR | 53 | 106 | Best powered; ~10% min detectable disparity |
| Non-Finnish European | NFE | 35 | 70 | Well powered; ~12% min detectable |
| East Asian | EAS | 24 | 48 | Moderate power; ~15% min detectable |
| South Asian | SAS | 19 | 38 | Limited power; ~18% min detectable |
| Admixed American | AMR | 5 | 10 | Underpowered; Fisher exact only; ~30% min detectable |
| **Total** | | **136** | **272** | |

**Data required:** Available now from HPRC cohort metadata.

---

## Additional Pipeline Requirements (Not Currently Implemented)

The following are needed to produce the figures above but are not part of the current pipeline output:

1. **Subgraph visualization for Figure 3 middle panels:** Requires running `odgi viz` or `odgi layout` on the .og files at specific loci to produce a visual representation of the pangenome graph structure. Alternative: hand-drawn schematics based on published structural variant descriptions.

2. **BED file generation for distribution:** Requires a post-processing script that reads per-profile diversity reports and writes GA4GH BED v1.0 files with population-stratified coverage columns.

3. **UCSC track hub configuration:** Requires `hub.txt`, `genomes.txt`, and `trackDb.txt` files pointing to bigBed versions of the BED files. `bedToBigBed` conversion needed.

4. **Fisher's exact test fallback:** The current pipeline uses chi-squared only. Needs modification to fall back to Fisher's exact when expected cell counts < 5 or for all AMR comparisons.

5. **FDR correction option:** Benjamini-Hochberg as alternative to Bonferroni for high-probe-count technologies (CMA, OGM).

---

## Production Software

- **Figure generation:** Python (matplotlib + seaborn) for data-driven panels; Inkscape or Adobe Illustrator for schematics and final assembly
- **Statistical analysis:** scipy.stats (chi2_contingency, fisher_exact), statsmodels (multipletests for BH FDR)
- **BED generation:** Custom Python script reading pipeline JSON outputs
- **Track hub:** UCSC kent tools (bedToBigBed), manual hub configuration files
- **Color palette:** Colorblind-safe throughout (viridis, or custom diverging palette tested with Color Oracle)

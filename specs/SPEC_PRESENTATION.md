# Presentation Specification: ProbeProjection Corpus

**Version:** 1.0
**Date:** 2026-03-30
**Status:** Pre-execution (pipeline blocked on graph data; specifications written prospectively)

---

## 1. Framing

### 1.1 One-Sentence Summary

The first systematic projection of clinical assay probe coordinates from 16 diagnostic technologies through the HPRC human pangenome, quantifying where population-specific structural variation overlaps clinical assay targets across 625 profiles and ~272 haplotypes.

### 1.2 What It Is

A coordinate-level structural variant screen of clinical assay probe regions. For each (probe, haplotype) pair, the pipeline determines whether the probe's GRCh38 reference coordinates project to a valid location on the target haplotype path in the HPRC pangenome graph. Failures indicate that a large structural variant (deletion, inversion, rearrangement) overlaps the probe region on that haplotype. Population-stratified statistical testing (chi-squared or Fisher's exact, Bonferroni-corrected) identifies ancestry-correlated patterns in these failures.

### 1.3 What It Is Not

Not a hybridization simulator. The pipeline detects coordinate-level disruptions, not sequence-level mismatches within probe footprints. For short probes (CMA 25-mers, MLPA ~50nt), clinically relevant failure modes (single mismatches affecting hybridization or ligation) are outside the detection scope. This is stated explicitly in the manuscript.

### 1.4 Why Coordinate Projection Is Still Valuable

The pipeline's detection capability aligns with the technologies where coordinate-level events are the dominant clinical failure mode:

- **FISH probes (50kb-1Mb BACs):** Immune to SNP-level variation. Signal loss requires large structural events (deletions, inversions). Coordinate projection captures >90% of clinically relevant failures.
- **HLA/PGx targets:** Extreme structural polymorphism (whole-gene deletions, duplications, complex rearrangements) is the primary source of population-specific variation.
- **OGM (DLE-1 labels):** Pattern-level technology; large SVs are the relevant failure mode.

For CMA and MLPA, coordinate projection captures only large-scale events and explicitly does not detect the sequence-level mismatches that drive most per-probe failures. This is a methodology limitation, not a flaw — it scopes the contribution honestly.

### 1.5 Relationship to Prior Work

**Dong/Heiss et al. (Cell Reports, June 2025):** Projected Illumina methylation array probes against 94 HPRC haplotype assemblies. Identified population-specific "unambiguous probe" sets. This validates the computational approach for a single technology.

**This work extends Dong et al. by:**
- Covering 16 clinical technologies (vs. methylation arrays only)
- Using the full HPRC cohort (136+ samples, 272+ haplotypes vs. 47 samples, 94 haplotypes)
- Applying formal population-stratified statistical testing with multiple testing correction
- Providing a unified cross-technology probe catalog (1,628 profiles)

---

## 2. Contributions

Ranked by combined defensibility and impact.

### C1. Multi-Technology Pangenome Probe Projection (LEAD)

**Claim:** First study to project probes from 13+ clinical technologies through the HPRC pangenome in a unified framework.

**Evidence for novelty:** Literature search (March 2026) found no peer-reviewed multi-technology pangenome probe projection study. The closest precedent (Dong et al. 2025) covers one technology. The literature is siloed by technology — CMA, MLPA, FISH, NGS, and OGM communities operate independently with respect to probe validation.

**Scoping:** "To our knowledge, no peer-reviewed study has simultaneously projected probe coordinates from 13+ clinical technologies through the HPRC pangenome and reported population-stratified coverage statistics."

### C2. FISH Smoking Guns at Known SV Loci

**Claim:** Specific, named clinical FISH assays show demonstrable population-specific coordinate projection failures at known structural variant loci.

**Why this is the narrative anchor:**
- FISH has the strongest causal chain from projection failure to clinical relevance
- 319 profiles (largest technology group) provide statistical breadth
- Tier 1 candidate loci have well-characterized population-specific structural variants:
  - **17q21.31 (KANSL1/MAPT):** H2 inversion at ~20% in Europeans, ~0% in AFR/EAS. 1.08-1.49 Mbp. Clinical FISH for Koolen-De Vries syndrome.
  - **8p23.1 (defensin cluster):** ~4.5 Mb inversion, ~26% European. Clinical FISH for 8p deletion/duplication syndromes.
  - **22q11.2 (DiGeorge/VCFS):** 50 distinct structural configurations. 11-fold length variation. Population-specific SD expansion. Clinical FISH (TUPLE1/ARSA).
  - **CYP2D6 (22q13.2):** Whole-gene deletion frequency varies 3-5x by population. PGx star allele panels.

**Note:** 17q21.31 is a case where European-ancestry haplotypes show the disparity (H2 inversion), not African-ancestry. This is scientifically interesting and avoids the framing that only non-European populations are affected.

### C3. HLA/PGx High-Disparity Zone

**Claim:** HLA (38 profiles) and PGx (23 profiles) show the highest rate of population-specific projection variation because they target the genome's most polymorphic regions.

**Expected magnitude:** 10-30% of targets with population-specific signals. Biologically expected but never quantified in this framework.

### C4. Reassuring Baseline

**Claim:** >95% of probes project equivalently across all superpopulations at the coordinate level. The existing clinical assay landscape is largely robust to population-specific structural variation.

**Why this matters:** This is actionable for clinical lab directors and accreditation. It provides evidence that most existing assay designs are sound, while identifying the specific exceptions that require attention.

### C5. Normalized 16-Technology Probe Catalog

**Claim:** 1,628 NormalizedProfile JSONs with standardized coordinates across 16 technologies. No equivalent public resource exists.

**Standalone value:** This data contribution is independent of the projection results. It enables any downstream analysis that requires unified probe coordinate data across clinical technologies.

### C6. Methodology Characterization

**Not a claim — context:** Coordinate projection's detection capability varies by technology. This is presented as honest methodology characterization: what the pipeline can see, what requires sequence-level analysis, and where future work should focus.

---

## 3. Expected Findings

Based on structural variation data from HPRC, gnomAD-SV, and 1000 Genomes Project long-read SV calls.

| Technology | Expected Disparity Rate | Confidence | Notes |
|---|---|---|---|
| HLA | 20-50% of targets | High | Extreme polymorphism by design |
| PGx | 10-30% of targets | High | CYP2D6, CYP2A6 well-characterized |
| FISH | 5-15% of assays with >=1 probe signal | High | Concentrated at known SV loci |
| OGM | 0.5-2% of labels | Medium | Pattern-level, not probe-level |
| NGS | 1-5 genes per panel | Medium | Concentrated in SV-rich loci |
| CMA | 685-6,850 probes (0.01-0.1%) | Medium | None pass Bonferroni; magnitude-only |
| MLPA | 0-2 probes total | Low | Most failures invisible to pipeline |

---

## 4. Statistical Framework

### 4.1 Primary Tests

- **Chi-squared test of independence:** populations x {mapped, unmapped} per probe. Default for all comparisons where expected cell counts >= 5.
- **Fisher's exact test:** Fallback for any comparison involving AMR (n=10) or any cell with expected count < 5.
- **Bonferroni correction:** p_corrected = p_raw x n_probes. Conservative; appropriate for FISH/MLPA with small probe counts.
- **Benjamini-Hochberg FDR:** Alternative for high-probe-count technologies (CMA, OGM) where Bonferroni is too conservative.

### 4.2 Significance Criteria

- Bonferroni-corrected p < 0.05 AND disparity magnitude > 20% (absolute difference between most- and least-affected populations)
- For CMA (6.85M probes): Bonferroni threshold (7.3e-9) is unachievable with these sample sizes. Rely on magnitude threshold with uncorrected p as supporting evidence.

### 4.3 Power by Superpopulation

| Superpopulation | Haplotypes | Min Detectable Disparity (FISH) | Min Detectable (NGS, 500 probes) |
|---|---|---|---|
| AFR | 106 | ~10% | ~10% |
| NFE | 70 | ~12% | ~12% |
| EAS | 48 | ~15% | ~15% |
| SAS | 38 | ~18% | ~18% |
| AMR | 10 | ~30% (Fisher exact only) | Insufficient |

### 4.4 AMR Handling

AMR (n=10 haplotypes, 5 samples) is underpowered for all but extreme effects. All AMR results are:
- Computed with Fisher's exact test (not chi-squared)
- Flagged in output as LOW_POWER
- Reported separately in tables/figures with visual distinction
- Discussed as a limitation in the manuscript

---

## 5. Publication Strategy

### 5.1 Target Journals

**Primary:** Genome Research
- Strong computational genomics scope
- Has published HPRC ecosystem work
- Methods + resource framing fits well

**Secondary:** American Journal of Human Genetics
- Clinical genetics audience
- Health equity angle resonates
- Higher visibility but higher bar

**Alternative:** Genetics in Medicine
- Reaches clinical lab directors most directly
- ACMG-affiliated, aligns with ACMG 2023 "points to consider" citation

### 5.2 Framing in Abstract

Lead: "We present the first systematic projection of clinical assay probes from 16 technologies through the HPRC human pangenome."

Cite: Dong et al. (2025) as single-technology precedent. ACMG (2023) as professional endorsement of the need.

Findings: "[X]% of probes across [Y] technologies project equivalently to all HPRC superpopulations at the coordinate level. We identify [N] specific loci across [M] clinical assays where population-specific structural variation disrupts probe coverage, including [lead example]."

Close: Data distributed as BED files and UCSC track hub. All outputs RUO.

### 5.3 Explicit Non-Claims

The manuscript will NOT:
- Claim detection of hybridization-level probe failures
- Claim this will change FDA or CAP guidance (contributes to evidence base only)
- Compare the distribution format to ClinVar or gnomAD
- Present AMR results as well-powered
- Claim CMA/MLPA sequence-level probe failures are detected
- Overstate the fraction of probes affected

### 5.4 Timeline

| Phase | Duration | Depends On |
|---|---|---|
| Unblock graph data (rebuild .og from .gbz) | 1-2 weeks | Immediate priority |
| Run FISH profiles (319, fastest validation) | 1-2 weeks | Graph data |
| Run HLA + PGx (61 profiles) | 1 week | Graph data |
| Run remaining technologies | 2-4 weeks | Graph data |
| Analysis and figure generation | 2 weeks | Pipeline results |
| Manuscript draft | 2-3 weeks | Analysis |
| Internal review and revision | 1-2 weeks | Draft |
| **Preprint submission** | **Target: September 2026** | |

---

## 6. Distribution

### 6.1 Supplementary BED Files (Primary)

Per-technology BED files containing:
- Probe coordinates (chr, start, end)
- Probe ID and assay product ID
- Per-superpopulation coverage rate (AFR, EAS, NFE, SAS, AMR)
- Overall coverage rate
- Disparity magnitude
- Statistical significance flag

Format: GA4GH BED v1.0 standard. Usable with bedtools, UCSC Genome Browser, IGV. Zero ongoing maintenance.

### 6.2 UCSC Track Hub (Secondary)

- Hosted on UCSC HubSpace (free, 10GB per user)
- Integrated with existing HPRC assembly hub (132 assemblies)
- Tracks: probe coverage by technology, population-stratified disparity scores
- bigBed format for efficient random access
- Low maintenance: upload once, update with new HPRC releases

### 6.3 Not Building: Standalone Web Tool

A standalone web tool is not feasible without institutional backing. ClinVar (NCBI, 10+ years, federal funding), gnomAD (Broad Institute, 800K+ individuals), and DECIPHER (Wellcome Sanger, Wellcome Trust funding) are products of sustained institutional investment. BED files + UCSC track hub achieve equivalent accessibility with sustainable maintenance.

---

## 7. Regulatory Context (Honest Framing)

### What exists:
- FDA 2018 NGS guidance: recommends "diverse genotypes" in validation specimens. General, not prescriptive.
- ACMG 2023 "points to consider": identifies pangenome references as a strategy to address health inequities. Advisory, not binding.
- ACMG 2025 GRCh38 transition bulletin: only 7% of labs migrated to GRCh38 as of 2021. Pangenome adoption is a decade or more away.
- FDA LDT rule: vacated (March 2025), rescinded (September 2025). Regulatory landscape for LDTs is loosening, not tightening.

### What doesn't exist:
- No FDA guidance requires ancestry-stratified probe validation
- No CAP checklist item requires population-specific probe performance demonstration
- No FDA enforcement action has targeted population-specific assay performance
- No AMP/CAP validation consensus mentions ancestry-specific validation

### Defensible framing:
"This data contributes to the evidence base that professional societies and, eventually, regulatory bodies may draw upon when developing future guidelines. The most immediate utility is to clinical lab directors evaluating assay selection, to manufacturers assessing probe content, and to the genomics equity research community."

---

## 8. Key References

| Reference | Relevance |
|---|---|
| Dong/Heiss et al., Cell Reports 2025 | Methylation array probe projection (single-technology precedent) |
| Andrews et al., Clin Epigenetics 2016 | "Gap hunting" — SNPs under array probes cause population-stratified signal disruption |
| Liao et al., Nature 2023 | HPRC Year 1 draft pangenome (47 assemblies, 119 Mbp added) |
| Garg et al., Nature 2025 | 1000 Genomes long-read SV calls (167,291 SVs, 1,434 highly population-stratified) |
| Matalon et al., GIM 2023 | ACMG "points to consider" on biases in clinical genetics testing |
| ACMG Lab QA Bulletin, GIM Open 2025 | GRCh37-to-GRCh38 transition; 7% lab adoption as of 2021 |
| FDA NGS Germline Guidance, 2018 | "Diverse genotypes" language (general, not prescriptive) |
| Collins et al., Nature 2020 | gnomAD-SV catalog (publication model: vast baseline + targeted findings) |

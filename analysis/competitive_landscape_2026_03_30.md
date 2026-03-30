# Competitive Landscape Analysis: Pangenome-Based Clinical Assay Probe Evaluation

**Date:** 2026-03-30
**Analyst:** Research Agent
**Scope:** HPRC clinical activity, conference abstracts, preprints/publications, industry activity, competing groups
**Assessment period:** 2025-01 through 2026-03

---

## Executive Summary

No group has published or publicly presented work that directly replicates the ProbeProjection approach: systematic projection of clinical assay probe coordinates (FISH, CMA, MLPA, NGS panels, OGM, etc.) through HPRC pangenome graphs to quantify population-stratified coverage disparities across 16 clinical technologies.

However, the competitive environment is tightening on multiple fronts. One published paper (Dong et al., Cell Reports 2025) has already applied an analogous methodology -- mapping array probe sequences against HPRC haplotype assemblies to identify population-specific probe failures -- for methylation arrays specifically. Industry players (Thermo Fisher, Illumina) are actively building pangenome-informed products. The HPRC consortium itself, now in Phase 2 with explicit "community adoption" and "outreach" mandates, is moving toward clinical translation, though it lacks a dedicated clinical working group.

**Overall scoop risk: MEDIUM, trending toward MEDIUM-HIGH.**

The specific niche of projecting clinical cytogenomic probes (FISH, CMA, MLPA) through pangenome graphs remains unoccupied. But the surrounding territory is being colonized rapidly. The window for a first-mover publication is approximately 6-12 months before adjacent work closes in from the methylation side, the variant-calling side, or the industry product side.

---

## 1. HPRC Clinical Activity Summary

### 1.1 Working Groups

The HPRC has nine working groups. None is explicitly focused on clinical applications:

1. Population Sampling and Representation (Karen Miga, Eimear Kenny)
2. Technology and Production (Karen Miga, Bob Fulton)
3. T2T Assembly (Evan Eichler, Karen Miga, Adam Phillippy)
4. Pangenomes (Ira Hall, Heng Li, Benedict Paten)
5. Resource Improvement and Maintenance (Tina Lindsay, Valerie Schneider, Fergal Martin)
6. Resource Sharing and Outreach (Heather Lawson)
7. ELSI (Bob Cook-Deegan, Malia Fullerton, Nanibaa' Garrison)
8. International Outreach (Heather Lawson, Angela Page)
9. Publication (Benedict Paten, Ting Wang)

**Key finding:** There is no "Clinical Applications" or "Clinical Genomics" working group within HPRC. This is significant -- it means clinical probe evaluation is not a formal consortium deliverable, at least not yet.

### 1.2 Phase 2 Mandate

Phase 2 of NHGRI funding began December 2024 and explicitly includes:
- Adding ~200 more individuals to the resource
- Emphasis on "outreach and community adoption"
- Development of "informatics tools that will enable the wider community to use the pangenome reference resource to improve their research"

This mandate increases the likelihood that HPRC-affiliated researchers will pursue clinical translation work, even without a formal clinical working group.

### 1.3 ASHG 2025 Workshop

The HPRC held a workshop at ASHG 2025 titled "The Human Pangenome: Data, Tools, and Workflows" (speakers: Heather Lawson, Benedict Paten, Karen Miga, Eimear Kenny). The workshop taught attendees to:
- Access HPRC data and resources
- Conduct variant analyses using the pangenome
- Map short- and long-read data to the pangenome

This workshop specifically framed the pangenome as relevant to "identifying and predicting the functional outcomes of variants across diverse populations" for "researchers and clinicians." This signals that the HPRC is actively promoting clinical uptake.

### 1.4 HPRC-Affiliated Researchers with Clinical Backgrounds

Several HPRC co-leads have clinical genomics relevance:
- **Eimear Kenny** (Mount Sinai): Founding Director, Institute for Genomic Health. Explicit health equity focus. Strong AI/genomics clinical integration.
- **Ira Hall** (Yale): Co-leads the pangenome working group. Research on disease-associated structural variants, with emphasis on underrepresented populations.
- **Evan Eichler** (U. Washington): Expert on structural variants in clinical contexts, including autism. Published 2026 paper on using pangenome linear references to find missing autism variants.
- **Charles Lee** (Jackson Laboratory): FACMG credential. Investigates structural genomic variation and its role in human disease. Published on clinical cytogenomics.

**Assessment:** While no HPRC PI is working on clinical assay probe evaluation specifically, several have the expertise and infrastructure to pivot to such work on short notice. Eimear Kenny's health equity focus and clinical genomics position at Mount Sinai make her the most likely HPRC-affiliated PI to initiate clinically-oriented pangenome probe work.

---

## 2. Conference Abstract Findings

### 2.1 ASHG 2025

- The HPRC presented a plenary abstract on "Accurate representation of globally diverse human haplotypes in the second release of the human pangenome reference."
- A separate abstract addressed "Leveraging genome assembly and pangenome graphs to investigate an African-origin protective haplotype for Alzheimer's Disease in APOE-e4 Carriers."
- The plenary abstract PDF could not be fully parsed, so additional pangenome+clinical abstracts may exist but were not identified.
- No abstracts were found specifically addressing clinical assay probe evaluation using pangenome data.

### 2.2 ACMG 2026 (March 10-14, Baltimore)

**Critical finding:** Ambry Genetics presented an abstract titled "Use of Pangenome reference improves variant calling in clinical genome sequencing" (authors: Erica Smith, Julian Stone-Farhat, Matthew Schultz, Ali Alhafidh, David Brohawn, Marcus Mahar, Ermanno Florio).

Matthew Schultz, PhD (Director of Bioinformatics, Ambry Genetics) is described as "particularly dedicated to leveraging pangenomes to address the technical challenges of genomic equity." This is the closest commercial entity to directly addressing pangenome-driven diagnostic equity in a clinical laboratory context. However, their focus appears to be on variant calling in sequencing, not on assay probe evaluation.

### 2.3 AMP 2025 (November, Boston)

Bionano Genomics showcased 13 OGM studies at AMP 2025, highlighting "the increasing use of OGM and its capacity to detect novel and clinically relevant structural variants." No abstracts specifically linking pangenome data to OGM probe/label validation were identified.

### 2.4 ESHG 2025 (May, Milan)

Abstracts from the 58th ESHG were published in the European Journal of Human Genetics. No specific pangenome+clinical assay abstracts were identified, though full text search of the supplement was not feasible.

### 2.5 Summary

No conference abstract was identified that directly addresses the ProbeProjection approach. The Ambry Genetics ACMG 2026 abstract is the closest -- same thematic territory (pangenome + clinical genomics + equity) but different methodology (variant calling, not probe projection).

---

## 3. Preprint and Publication Scan

### 3.1 Directly Adjacent Work (HIGH relevance)

#### Dong et al. (2025) -- CLOSEST COMPETING WORK
**"Complete reference genome and pangenome improve genome-wide detection and interpretation of DNA methylation using sequencing and array data"**
- Published: Cell Reports, June 2025
- Authors: Zheng Dong, Joanne Whitehead, Maggie Fu, Julia L MacIsaac, David H Rehkopf, Luis Rosero-Bixby, Michael S Kobor, Keegan Korthauer
- Institutions: University of British Columbia, Stanford, Universidad de Costa Rica
- **What they did:** Mapped Illumina methylation array probes (HM450K, EPIC, EPICv2) against 94 HPRC haplotype-level assemblies to identify "unambiguous" probes -- those that perform consistently across populations.
- **Key results:**
  - Identified 5,858-12,035 additional unambiguous probes compared to T2T-CHM13 alone
  - Created population-specific unambiguous probe sets for AFR, AMR, EAS, EUR, SAS
  - Found population-specific probe numbers vary: HM450K 405,457-415,138; EPIC 734,436-751,048; EPICv2 820,417-832,466
  - EAS-specific unambiguous probes overlapped promoters and gene bodies of 140 cancer driver genes
- **Method:** BLAT alignment of probe sequences against pangenome assemblies (not graph-based projection via odgi, but conceptually analogous)
- **Code/data:** GitHub repository published (functionalepigenomics/Illumina_Infinium_HumanMethylation_BeadChips_Annotation)

**Assessment:** This is the most directly adjacent published work. They did for methylation arrays what ProbeProjection aims to do for 16 clinical technologies. Key differences: (1) they evaluated only methylation arrays, not FISH/CMA/MLPA/NGS/OGM; (2) they used BLAT alignment against individual assemblies, not graph-based coordinate projection via odgi; (3) they did not perform population-stratified chi-squared testing for disparity detection. Their existence validates the approach but also means the "pangenome probe evaluation" concept is no longer entirely novel. ProbeProjection's differentiation lies in the breadth (16 technologies), the graph-based method, and the statistical framework for disparity detection.

#### T2T-CHM13 Reference Bias Paper (2025)
**"T2T-CHM13 reference genome reduces mapping bias and enhances alignment accuracy at disease-associated variants"**
- bioRxiv, December 2025
- Key finding: "Sequence dissimilarities among reference genomes in the proximity of ClinVar annotated variants" suggest "the need for data re-analysis and potential redesign of probes targeting clinically relevant regions."
- This paper explicitly calls for probe redesign but does not perform the systematic evaluation itself.

### 3.2 Closely Related Work (MEDIUM relevance)

#### Heng Li -- Easy Regions for Variant Calling (2025)
- Derived sample-agnostic "easy regions" from pangenome data covering 88.2% of GRCh38, 92.2% of coding regions, 96.3% of ClinVar pathogenic variants.
- Relevant because it identifies which genomic regions are tractable vs. problematic for reference-based analysis. Does not evaluate clinical assay probes specifically.

#### Swave -- Population-level SV Characterization (Nature Genetics 2026)
- Wang et al. "Population-level structural variant characterization using pangenome graphs"
- Deep learning method resolving SVs and their population characteristics from pangenome graphs.
- Applied to healthy and rare-disease cohorts. Does not evaluate clinical assay probes.

#### Pangenome-Based Rare Disease Diagnostics (medRxiv 2025)
- Korean rare disease study: long-read sequencing + pangenome graph workflow identified pathogenic variants in 27.3% of previously unsolved families.
- Demonstrates clinical utility of pangenome for diagnostics but does not evaluate assay probe coverage.

#### Using Linear References from the Pangenome to Discover Missing Autism Variants (Nature Communications 2026)
- Eichler and colleagues. Uses pangenome data for clinical variant discovery in autism.
- Adjacent thematically but does not address probe evaluation.

#### Reference Genome Choice Compromises Population Genetic Analyses (Cell 2025)
- Akopyan et al. Demonstrated that reference choice affects diversity estimates by >30%, FST values, and effective population size.
- Performed in gray foxes but the principle applies directly to human clinical assays.

### 3.3 Review Articles

#### Nyaga et al. (Frontiers in Genetics 2025)
**"Beyond single references: pangenome graphs and the future of genomic medicine"**
- From University of Auckland and University of Otago, New Zealand
- Explicitly notes: "patients of non-European ancestry experience substantially lower diagnostic rates" and ~23% increase in VUS burden
- Notes tension between pangenome complexity and clinical interpretability
- Does NOT discuss clinical assay probes (FISH, CMA, MLPA) -- focuses on sequencing-based diagnostics

### 3.4 Summary

No published paper or preprint was identified that replicates the full ProbeProjection approach (multi-technology clinical assay probe projection through pangenome graphs with population-stratified disparity analysis). The Dong et al. (2025) methylation array paper is the closest and validates the general approach for one technology. The gap for FISH, CMA, MLPA, NGS panels, OGM, and other clinical cytogenomic technologies remains open.

---

## 4. Industry Activity

### 4.1 Thermo Fisher Scientific

**Most active in pangenome-informed array products.**

- **Axiom PangenomiX Array** (launched January 2024): Largest and most ethnically diverse array, combining SNP genotyping, whole-genome CNV detection, fixed copy number discovery, and HLA typing. >800,000 markers selected from 1000 Genomes Project phase 3 across EUR, AFR, AMR, EAS, SAS populations. Already used on nearly half a million samples at a major US biobank.
- **Axiom PangenomePro Array** (launched October 2025): 942,415 markers. Paired with new SwiftArrayStudio Microarray Analyzer and Axiom PharmacoPro Array.

**Assessment:** Thermo Fisher is building pangenome-aware genotyping products but these are research-use-only population genomics arrays, not redesigned clinical diagnostic arrays (FISH, CMA, MLPA). They are not evaluating existing clinical assay probes against pangenome data -- they are designing new research genotyping products. However, their PangenomiX branding signals awareness of the market opportunity in pangenome-informed assay design.

### 4.2 Illumina

- **DRAGEN v4.3** (June 2024): Expanded pangenome reference from 32 to 128 population samples across 26 ancestries. 40% error reduction in variant calling.
- **DRAGEN v4.4** (May 2025): 30% improvement in germline SV calling accuracy.
- **DRAGEN v4.5** (2026 webinar): Enhanced pangenome reference, ML-driven variant calling in hard-to-resolve regions.
- **DesignStudio**: Web-based tool for designing custom array probes and NGS assays.

**Assessment:** Illumina is deeply invested in pangenome-aware variant calling but has not publicly announced pangenome-based evaluation or redesign of their existing clinical array products (CytoSNP, Infinium, etc.) or capture panel designs. Their focus is on improving sequencing analysis, not on evaluating legacy probe-based assays.

### 4.3 Bionano Genomics

- Participated in HPRC as a technology provider (optical maps for scaffolding assemblies).
- Held AMP 2025 symposium with 13 OGM studies.
- Announced Bionano Symposium 2026 (February 23-26) covering hematologic malignancies, solid tumors, constitutional genetic disorders, and gene/cell therapy.
- Maintains a dataset of 279 OGM controls from diverse ethnic backgrounds (24.6% African, 8.9% Admixed American, 9.5% East Asian, 24.6% European, 8.4% South Asian).

**Assessment:** Bionano has diverse population OGM data but has not published pangenome-based evaluation of DLE-1 label site coverage across populations. Their focus is on demonstrating OGM's SV detection capabilities, not on questioning whether OGM's labeling chemistry works equitably across populations. This is a gap ProbeProjection could fill.

### 4.4 Patent Landscape

No patents were identified for pangenome-aware clinical probe design or evaluation. The patent landscape search returned results for CRISPR diagnostics and nanopore sequencing patents but nothing directly relevant to pangenome-based probe validation.

### 4.5 Summary

Industry is building pangenome-informed products (Thermo Fisher PangenomiX/PangenomePro, Illumina DRAGEN pangenome reference) but is not publicly evaluating whether existing clinical assay probes work equitably across pangenome-represented populations. The gap between "new pangenome-aware products" and "evaluation of existing clinical assay portfolios" is exactly where ProbeProjection sits.

---

## 5. Groups to Watch

### Tier 1: Most likely to produce closely related work (estimated 12-18 months)

#### 1. Keegan Korthauer Lab (University of British Columbia)
- **Published:** Dong et al. 2025 (Cell Reports) -- pangenome-based methylation array probe evaluation
- **Why they matter:** Already demonstrated the core methodology for one array type. Natural next step is to extend to clinical diagnostic arrays (CMA, FISH panel design, etc.) or to collaborate with clinical labs.
- **Risk level:** HIGH for methylation arrays (already published). MEDIUM for extension to other technologies.
- **Mitigation:** ProbeProjection covers 16 technologies including FISH, CMA, MLPA, NGS, OGM -- far broader scope. Cite Dong et al. as precedent for one technology while demonstrating the approach across the full clinical landscape.

#### 2. Ambry Genetics / Matthew Schultz
- **Presented:** ACMG 2026 abstract on pangenome + clinical variant calling + genomic equity
- **Why they matter:** Commercial clinical lab with explicit pangenome equity focus, bioinformatics expertise, and access to massive clinical cohorts. Schultz is "dedicated to leveraging pangenomes to address the technical challenges of genomic equity."
- **Risk level:** MEDIUM. Current focus is variant calling, not probe evaluation. But they have the data, infrastructure, and motivation to pivot.
- **Mitigation:** Publish first. Ambry's strength is clinical validation on patient cohorts; ProbeProjection's strength is in silico population-wide evaluation across all technologies.

#### 3. Eimear Kenny Lab (Mount Sinai / Institute for Genomic Health)
- **Why they matter:** HPRC co-lead. Founding Director of Institute for Genomic Health. Explicit health equity mandate. AI/genomics clinical integration. Access to BioMe Biobank.
- **Risk level:** MEDIUM. Has not published on probe evaluation specifically, but has the infrastructure, mandate, and expertise to do so.
- **Mitigation:** Kenny's focus appears to be on genomic medicine implementation (sequencing, interpretation), not on evaluating legacy cytogenomic assays.

### Tier 2: Could produce adjacent work (estimated 18-24 months)

#### 4. Ira Hall Lab (Yale)
- **Why they matter:** Co-leads HPRC pangenome working group. Focuses on SV detection methods. Works on "unbiased genomic analyses" in underrepresented populations. Recent work on rare disease pangenome diagnostics.
- **Risk level:** LOW-MEDIUM. Focus is on SV discovery and genome assembly, not on evaluating clinical assay probe coverage.

#### 5. Benedict Paten Lab (UC Santa Cruz)
- **Why they matter:** Co-leads HPRC pangenome working group. Develops core pangenome graph infrastructure (vg toolkit). Giraffe mapper now supports long reads (published 2025).
- **Risk level:** LOW-MEDIUM. Focus is on tooling and infrastructure, not on clinical assay evaluation. However, a student or postdoc could easily apply vg tools to probe projection.

#### 6. Andrea Guarracino Lab (TGen, Phoenix)
- **Why they matter:** Lead developer of odgi (the exact tool ProbeProjection uses for coordinate projection). Lab opened January 2026. Focus on "cancer pangenomics."
- **Risk level:** LOW-MEDIUM. Guarracino has the most intimate knowledge of odgi position liftover. A cancer pangenomics application using odgi to project clinical panel probes is conceivable.

#### 7. Charles Lee Lab (Jackson Laboratory)
- **Why they matter:** FACMG-credentialed. Works on structural genomic variation and clinical cytogenomics. Contributed to HPRC complex genetic variation paper (Nature 2025).
- **Risk level:** LOW-MEDIUM. Has clinical cytogenomics expertise and pangenome involvement, but published output focuses on SV characterization, not probe evaluation.

### Tier 3: Working in adjacent space, lower immediate risk

#### 8. Evan Eichler Lab (University of Washington)
- Recent 2026 Nature Communications paper on autism variant discovery using pangenome linear references.
- Focus on complex SV biology, not on clinical assay probe evaluation.

#### 9. Tobias Marschall Lab (Heinrich Heine University, Dusseldorf)
- Developer of PanGenie for pangenome-based genotyping. Published Locityper (2025) for challenging gene genotyping.
- Focus on genotyping algorithms, not clinical assay evaluation.

#### 10. Michael Schatz Lab (Johns Hopkins)
- AnVIL project PI. Work on graph genomes for rare disease SV analysis (>200,000 SV alleles detected).
- Focus on variant discovery infrastructure, not probe evaluation.

---

## 6. Scoop Risk Assessment

### Rating: MEDIUM, trending toward MEDIUM-HIGH

### Rationale

**Factors increasing risk:**
1. Dong et al. (2025) has already demonstrated the core approach for methylation arrays, validating the concept and reducing novelty.
2. The HPRC is entering Phase 2 with an explicit community adoption mandate.
3. Ambry Genetics is actively working on pangenome + clinical equity at the commercial lab level.
4. Thermo Fisher is branding products around "pangenome" (PangenomiX, PangenomePro).
5. Multiple review papers (Nyaga et al. 2025, among others) are calling attention to pangenome + diagnostic equity.
6. The ASHG 2025 workshop explicitly targeted clinicians and researchers.
7. The T2T-CHM13 mapping bias paper (2025) explicitly calls for "potential redesign of probes targeting clinically relevant regions."

**Factors reducing risk:**
1. No HPRC clinical working group exists.
2. No one has published on FISH, CMA, MLPA, NGS panel, or OGM probe evaluation using pangenome data.
3. The technical barrier is non-trivial (graph-based coordinate projection, not just BLAT alignment).
4. The clinical cytogenomics community and the pangenome toolchain community do not overlap much.
5. HPRC PIs are focused on genome assembly, SV discovery, and variant calling -- not on evaluating existing clinical assays.
6. Industry focus is on building new products, not retrospectively evaluating existing ones.
7. The project's scope (16 technologies, 1,628 profiles, 272 haplotypes) is substantially broader than anything published or in progress.

### Timeline Recommendation

**Aim to submit a preprint within 6 months (by September 2026).**

Reasoning:
- The methylation array precedent (Dong et al. 2025) means the general concept is now public. Someone will extend it to other technologies.
- ACMG 2027 abstracts will be due in late 2026. If Ambry or another clinical lab pursues pangenome-based probe evaluation, they would likely submit by then.
- A Nature Genetics or Genome Research preprint posted by September 2026 would establish clear priority before the next wave of ASHG/ACMG abstracts.

**Recommended publication strategy:**
1. **Preprint first** (bioRxiv) to establish priority, even with partial results (e.g., FISH + MLPA + NGS panels only).
2. **Full paper** targeting Genome Research, Nature Genetics, American Journal of Human Genetics, or Genome Medicine.
3. **Conference abstracts** for ASHG 2026 (deadline likely June-July 2026) and ACMG 2027.

### Critical Path

The project is currently blocked on graph data (`.og` files lack haplotype paths). Resolving this blocker and running even a subset of technologies (e.g., 319 FISH profiles) would produce publishable results that no one else has. The FISH results alone, showing population-stratified probe coverage disparities for ~320 clinical FISH assays across 272 HPRC haplotypes, would be a first-of-kind publication.

---

## Sources

- [HPRC Working Groups](https://humanpangenome.org/working-groups/)
- [HPRC Data Release 2](https://humanpangenome.org/hprc-data-release-2/)
- [ASHG 2025 Pangenome Workshop](https://www.ashg.org/product/ashg-2025-workshop-the-human-pangenome-data-tools-and-workflows/)
- [ASHG 2025 Updates and Applications of the Human Pangenome Reference](https://www.ashg.org/product/updates-and-applications-of-the-human-pangenome-reference/)
- [Ambry Genetics ACMG 2026](https://www.ambrygen.com/conference/580/acmg)
- [Dong et al. 2025 - Cell Reports - Pangenome methylation array probe evaluation](https://www.cell.com/cell-reports/fulltext/S2211-1247(25)00526-1)
- [Dong et al. GitHub repository](https://github.com/functionalepigenomics/Illumina_Infinium_HumanMethylation_BeadChips_Annotation)
- [Nyaga et al. 2025 - Frontiers in Genetics - Pangenome review](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2025.1679660/full)
- [T2T-CHM13 mapping bias paper (bioRxiv 2025)](https://www.biorxiv.org/content/10.64898/2025.12.17.694618v2)
- [Thermo Fisher PangenomiX Array](https://www.thermofisher.com/us/en/home/life-science/microarray-analysis/applications/predictive-genomics/population-genomics/arrays/axiom-pangenomix.html)
- [Thermo Fisher PangenomiX launch](https://www.businesswire.com/news/home/20240104569788/en/Thermo-Fisher-Scientific-Launches-PangenomiX-Array-to-Support-Population-Scale-Disease-Studies-and-Pharmacogenomics-Research)
- [Illumina DRAGEN v4.3 pangenome reference](https://emea.illumina.com/science/genomics-research/articles/second-gen-multigenome-mapping.html)
- [Illumina DRAGEN v4.5 webinar](https://www.illumina.com/events/webinar/2026/dragen-v4.5.html)
- [Bionano AMP 2025](https://bionanogenomics.gcs-web.com/news-releases/news-release-details/bionano-recaps-advances-optical-genome-mapping-showcased-amp)
- [Bionano Symposium 2026](https://ir.bionanogenomics.com/news-releases/news-release-details/bionano-announces-bionano-symposium-2026-global-experts-showcase/)
- [vg Giraffe long-read mapping 2025](https://www.biorxiv.org/content/10.1101/2025.09.29.678807v1.full)
- [Swave - Population-level SV characterization (Nature Genetics 2026)](https://www.nature.com/articles/s41588-026-02538-6)
- [Pangenome rare disease diagnostics (medRxiv 2025)](https://www.medrxiv.org/content/10.1101/2025.07.08.25330875v1.full)
- [Heng Li - Easy regions for variant calling 2025](https://pubmed.ncbi.nlm.nih.gov/40905373/)
- [Eimear Kenny profile](https://profiles.icahn.mssm.edu/eimear-kenny)
- [Ira Hall profile](https://medicine.yale.edu/profile/ira-hall/)
- [Charles Lee profile](https://www.jax.org/research-and-faculty/faculty/charles-lee)
- [Guarracino Lab](https://guarracinolab.github.io/)
- [Complex genetic variation in nearly complete human genomes (Nature 2025)](https://www.nature.com/articles/s41586-025-09140-6)
- [Using pangenome linear references to discover missing autism variants (Nature Communications 2026)](https://www.nature.com/articles/s41467-026-68378-4)
- [Reference genome choice compromises population genetic analyses (Cell 2025)](https://www.cell.com/cell/abstract/S0092-8674(25)01026-8)
- [WHO genomics equity gaps analysis (December 2025)](https://www.who.int/news/item/21-12-2025-who-publishes-new-global-analysis-revealing-major-equity-gaps-in-human-genomics-research)
- [ACMG CMA technical standard 2021 revision](https://www.sciencedirect.com/science/article/pii/S1098360021051285)
- [T2T-CHM13 improves clinically relevant variant detection in Swedish population](https://pmc.ncbi.nlm.nih.gov/articles/PMC12581843/)

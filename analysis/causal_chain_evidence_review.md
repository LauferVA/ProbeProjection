# Evidence Review: In Silico Pangenome Projection to In Vitro Probe Failure

**Date:** 2026-03-30
**Status:** Complete literature and technical review
**Scope:** Biophysics of probe failure, validation precedents, tool capabilities, gap estimation

---

## Table of Contents

1. [Task 1: Biophysics of Probe Failure by Technology](#task-1)
2. [Task 2: In Silico to In Vitro Validation Precedents](#task-2)
3. [Task 3: What Coordinate-Level Projection Can and Cannot Tell Us](#task-3)
4. [Task 4: Estimate the Detectable vs Undetectable Gap](#task-4)
5. [Deliverable 1: Technology-by-Technology Evidence Table](#deliverable-1)
6. [Deliverable 2: Validation Precedent Summary](#deliverable-2)
7. [Deliverable 3: Tool Comparison](#deliverable-3)
8. [Deliverable 4: Honest Assessment of the Causal Chain](#deliverable-4)

---

<a name="task-1"></a>
## Task 1: Biophysics of Probe Failure by Technology

### 1.1 CMA (Chromosomal Microarray)

#### Affymetrix 25-mer probes (e.g., CytoScan HD)

**Mechanistic basis for failure.** Affymetrix SNP/CMA arrays interrogate each
target with paired perfect-match (PM) and mismatch (MM) 25-mer probes. The MM
probe differs from the PM by a single base substitution at position 13 (the
center). The entire Affymetrix genotyping architecture is built on the
premise that a single centrally-positioned mismatch produces a measurable and
reliable signal reduction in a 25-mer context.

**Quantitative mismatch tolerance.** Published thermodynamic modeling and
empirical studies demonstrate the following for 25-mer probes:

- A single mismatch at the central position reduces hybridization signal by
  approximately 40-60% relative to the perfect match, though the exact
  magnitude depends on the identity of the mismatch pair, nearest-neighbor
  context, and GC content (Naef & Magnasco, 2003; Binder et al., 2009, PLOS
  ONE 4:e7862).
- Mismatches at positions closer to the probe termini have smaller effects
  than central mismatches, but even a single terminal mismatch is detectable.
- For SNP genotyping, the PM/MM signal ratio is the core discriminator:
  a SNP under the probe footprint produces a quantifiable shift in this ratio.
- The Affymetrix CytoScan HD platform contains 750,000 SNP probes and 1.9M
  copy-number probes, all 25-mers. Each SNP is interrogated by six probes
  (three per allele), providing redundancy.

**Practical implication.** A SNP anywhere within a 25-mer CMA probe's binding
site will alter hybridization intensity. For copy-number detection (not SNP
genotyping), the pipeline averages signal across many adjacent probes, so
individual probe failures are tolerated. However, if a population-specific
variant affects a cluster of adjacent probes, it could create a focal
signal dropout that mimics or masks a copy-number change.

**Detectability by coordinate projection.** A SNP within a 25-mer binding
site will NOT be detected by coordinate-level projection. The path through
the graph will exist (the probe coordinates map), but the underlying sequence
will differ at 1-2 positions. This is a clear blind spot. Only sequence-level
alignment of the probe to the haplotype would detect it.

#### Agilent 60-mer probes (e.g., SurePrint G3 CGH+SNP)

**Mechanistic basis for failure.** Agilent CGH arrays use 60-mer
oligonucleotide probes. Longer probes have higher absolute binding energy,
which means they tolerate more mismatches before signal is lost.

**Quantitative mismatch tolerance.** Published data on 50-60 mer probes
show:

- A single mismatch reduces signal by only ~40% of the PM/MM differential
  (i.e., the relative signal of a 1-mismatch probe is still ~0.60 of PM),
  meaning a single SNP under a 60-mer probe is insufficient to fully
  discriminate it from the perfect match (He et al., 2005, NAR 33:e62;
  Kreil et al., BMC Genomics 2008 9:491).
- Three to five mismatches are needed to reduce signal to background levels
  in 50-mer probes, depending on hybridization temperature (42-50 C).
- Agilent probes are designed with in silico validation against the reference
  genome for specificity, but this validation is inherently reference-biased.

**Practical implication.** Individual SNPs under 60-mer probes are unlikely
to cause clinical-grade signal dropout. The risk is lower than for 25-mers
on a per-probe basis. However, small indels (2-5 bp) or clustered SNPs
within the probe footprint could cause failure, particularly if they affect
the central portion of the probe.

**Detectability by coordinate projection.** Same blind spot as 25-mers:
coordinate projection cannot detect sequence-level mismatches. However, the
biological significance is lower because 60-mers are more mismatch-tolerant.

---

### 1.2 MLPA (Multiplex Ligation-Dependent Probe Amplification)

**Mechanistic basis for failure.** MLPA uses pairs of hemiprobes (~25-30 nt
each, total ~50-60 nt) that hybridize adjacently on the target DNA. A
thermostable ligase joins the two hemiprobes ONLY if they are perfectly
base-paired at and immediately around the ligation junction. This is
followed by PCR amplification of the ligated product. The critical
vulnerability is at the ligation junction itself.

**Quantitative mismatch tolerance.** This is the most sensitive technology
to single-nucleotide variants:

- A single mismatch at the ligation site (the terminal nucleotide of either
  hemiprobe) effectively prevents ligation, reducing signal to near zero.
  MRC Holland explicitly states that ligase will not join hybrids with
  mismatches next to the ligation site (Schouten et al., 2002, NAR 30:e57).
- This sensitivity is so extreme that MLPA probes can be deliberately designed
  to detect point mutations by placing the ligation junction at the variant
  site (mutation-specific MLPA).
- SNPs within the hybridization region but away from the ligation junction
  (>5 nt) may reduce hybridization stability but are less likely to prevent
  ligation entirely.
- MRC Holland documents that polymorphisms at or near the ligation site can
  produce apparent single-exon deletions (false positives for copy-number
  loss) because allele dropout from ligation failure mimics deletion.

**Known clinical impact.** MRC Holland explicitly warns that "apparent
deletion could actually consist of a sequence change hampering correct probe
hybridization" and that "apparent single exon deletions should be checked by
an independent method." This is a documented, acknowledged failure mode in
clinical practice.

**Detectability by coordinate projection.** Coordinate projection is almost
entirely blind to this failure mode. MLPA failure is caused by SNPs at
specific nucleotide positions within the ~50 nt probe footprint. The
haplotype path will exist in the graph. Only sequence extraction and
alignment of the probe to the haplotype sequence would identify the risk.
This is arguably the technology where the gap between coordinate projection
and actual probe performance is most dangerous, because a single SNP at the
ligation junction causes complete signal loss.

---

### 1.3 FISH (Fluorescence In Situ Hybridization)

**Mechanistic basis for failure.** FISH probes used in clinical cytogenetics
are typically BAC or fosmid clones spanning 50 kb to 1 Mb. They are labeled
with fluorophores and hybridized to denatured chromosomal DNA under
controlled stringency conditions (typically 37 C hybridization, 73 C
stringency washes in 0.4x SSC/0.3% NP-40). The probe is a complex mixture
of sequences across its entire insert length.

**Quantitative mismatch tolerance.** BAC/fosmid FISH probes are
fundamentally different from short oligonucleotide probes:

- Because the probe spans 100-200 kb, it contains thousands of independent
  hybridization events across its length. Overall signal intensity is the
  aggregate of all these events.
- Published cross-species CGH data shows that at ~90% sequence identity
  (i.e., 10% divergence), hybridization signal begins to degrade
  measurably but is not lost (Drosophila CGH studies; BMC Genomics
  2010 11:271). A roughly linear relationship exists between sequence
  divergence and signal reduction.
- The Tm drops by approximately 10 C when probes hybridize to homeologous
  sequences with several percent divergence.
- Within the human species, sequence divergence between haplotypes is
  typically 0.1-0.5% (SNP density), which is far below the threshold for
  FISH signal loss. Even in highly polymorphic regions (e.g., HLA), local
  divergence rarely exceeds 5%.
- What WILL cause FISH signal loss or splitting: large structural variants
  (deletions >50 kb eliminating part of the probe target), inversions
  relocating the target, or duplications producing extra signals. These are
  exactly what coordinate projection can detect.

**Detectability by coordinate projection.** This is the technology where
coordinate projection is MOST informative. The failure modes that matter
clinically for FISH (large SVs, inversions, translocations, deletions) are
precisely the events that cause a haplotype path to be absent, split, or
rearranged in the pangenome graph. SNP-level divergence is irrelevant to
BAC FISH performance. The pipeline's binary MAPPED/ABSENT classification is
well-matched to the actual biophysics of FISH probe behavior.

---

### 1.4 NGS Capture Probes (Hybridization-Based Target Enrichment)

**Mechanistic basis for failure.** NGS hybridization capture uses 80-120 nt
biotinylated RNA or DNA baits that hybridize in solution to target DNA
fragments. Captured fragments are pulled down with streptavidin beads and
sequenced. Unlike microarray probes (solid-phase), capture occurs in liquid
phase with excess bait, which increases tolerance for mismatches.

**Quantitative mismatch tolerance.**

- Published studies demonstrate that hybridization capture tolerates
  10-20% sequence divergence between bait and target while still achieving
  successful enrichment (Mamanova et al., 2010; Bi et al., 2012; NAR
  2016 44:4504).
- Cross-species capture experiments have successfully enriched targets
  across species separated by >200 million years of evolution (~10%
  sequence divergence in conserved regions), though efficiency degrades
  with divergence.
- There IS a measurable bias: probes designed from one genotype capture
  that genotype more efficiently. Mouse exome capture studies showed that
  C57BL/6J-designed probes captured C57BL/6J alleles more efficiently,
  with bias proportional to sequence divergence.
- RNA baits provide better hybridization stability than DNA baits,
  contributing to higher mismatch tolerance.

**Practical implication.** For human intra-species variation, sequence
divergence at most loci is well within the tolerance of capture probes.
The risk is concentrated at:
  (a) Highly polymorphic regions (HLA, KIR) where local divergence may
      approach or exceed 5%.
  (b) Regions with population-specific structural variants that eliminate
      the target entirely.
  (c) Complex rearrangements that cause capture of the wrong genomic
      region, leading to misinterpretation.

**Detectability by coordinate projection.** Coordinate projection can detect
(b) and (c) -- complete target absence and rearrangements. It cannot detect
the gradual degradation of capture efficiency due to sequence mismatches.
However, given the high mismatch tolerance of capture probes (10-20%), the
undetectable fraction is less clinically relevant than for CMA or MLPA.

---

### 1.5 OGM (Optical Genome Mapping)

**Mechanistic basis for failure.** OGM is fundamentally different from all
other technologies in this analysis. It does not use hybridization probes.
The DLE-1 methyltransferase enzyme covalently labels all CTTAAG hexamer
motifs in ultra-long DNA molecules (>150 kb median) with a fluorophore.
The pattern of label spacing along the molecule is then compared to an
in silico digest of the reference genome to identify structural variants.

**Quantitative mismatch tolerance.** The concept of "mismatch tolerance" does
not apply in the traditional sense:

- DLE-1 labels occur at approximately 14-15 sites per 100 kb in human
  genomic DNA.
- A SNP that destroys a CTTAAG motif (e.g., CTTAAG -> CTTATG) eliminates
  one label. A SNP that creates a new CTTAAG motif adds one label. These
  are rare events because any given SNP has only a small probability of
  affecting the specific 6-bp recognition sequence.
- Individual label gain/loss events are below the resolution of OGM SV
  calling, which requires consistent pattern changes across multiple
  molecules. The minimum SV detection threshold is approximately 500 bp
  for insertions/deletions and 30 kb for inversions and translocations.
- The failure mode for OGM is not probe hybridization failure but rather
  changes in label density patterns that can confuse the alignment algorithm
  or fall below detection thresholds.

**Detectability by coordinate projection.** Coordinate projection is
minimally relevant for OGM. The "probes" in the pipeline's OGM profiles
likely represent label sites or reference intervals, not hybridization
targets. Structural variants detectable by OGM (large indels, inversions,
translocations) would also be detectable by coordinate projection as path
absences or rearrangements. However, the primary clinical concern for OGM
equity is whether label density patterns differ systematically across
populations, which requires motif-level analysis, not coordinate projection.

---

<a name="task-2"></a>
## Task 2: In Silico to In Vitro Validation Precedents

### 2.1 Direct Validation Studies (In Silico Prediction vs. Wet-Lab Outcome)

#### 2.1.1 SARS-CoV-2 Diagnostic PCR Primer/Probe Failures

The strongest published precedent for in silico prediction of probe failure
validated against wet-lab outcomes comes from the COVID-19 diagnostic assay
literature:

- Vanaerschot et al. (2020, Royal Society Open Science 7:200636) and Khan
  & Cheung (2022, Scientific Reports 12:17953) systematically identified
  mutations in SARS-CoV-2 genomes that fell within primer and probe binding
  sites of WHO-recommended RT-qPCR assays.
- In vitro validation confirmed that single mismatches at the 3' end of
  primers could abolish PCR amplification, and even single mismatches in
  TaqMan probe regions reduced sensitivity and caused false-negative results.
- The ThermoFisher TaqPath S-gene target dropout caused by the B.1.1.7
  variant's del69-70 is a real-world example where a variant within a probe
  binding site caused clinical assay failure at population scale, initially
  detected in silico from genome surveillance data.

**Relevance.** This demonstrates the principle that in silico identification
of variants at probe/primer binding sites can predict wet-lab assay failure.
However, these are PCR-based assays (20-30 nt primers/probes), not
hybridization arrays, FISH, or capture. The mismatch intolerance of PCR
systems is well-characterized. The causal chain is shorter and better-
validated for PCR than for the technologies in this project.

#### 2.1.2 Methylation Array Probe Ambiguity: T2T and Pangenome Validation

The most directly relevant published study is:

- Heiss et al. (2025, Cell Reports): "Complete reference genome and pangenome
  improve genome-wide detection and interpretation of DNA methylation using
  sequencing and array data."

This study performed exactly the type of analysis this project proposes, but
for methylation arrays specifically:

- They aligned all probe sequences from Illumina 450K, EPIC, and EPICv2
  arrays to both GRCh38 and T2T-CHM13.
- They identified probes with cross-reactive mappings or mismatches that were
  not detected when using GRCh38 alone.
- Using the HPRC pangenome, they identified "cross-population and
  population-specific unambiguous probes" -- probes whose mapping status
  differs depending on which population haplotypes are considered.
- In EWASs of 24 cancer types, an average of 945 additional differentially
  methylated CpG sites were identified using the new unambiguous probe set
  vs. the GRCh38-based set.

**Critical gap.** This study validates probe AMBIGUITY (cross-reactivity,
multi-mapping) using the pangenome, but does not validate probe FAILURE
(loss of hybridization signal) against wet-lab data. They show that certain
probes are computationally problematic, but they do not show that these
probes actually produce aberrant signals in real experiments. The gap between
"computationally ambiguous" and "experimentally failed" remains open.

#### 2.1.3 Methylation Array "Gap Probes" and SNP Effects

- Andrews et al. (2016, Epigenetics & Chromatin): "'Gap hunting' to
  characterize clustered probe signals in Illumina methylation array data."

This study identified 11,007 "gap probes" on the 450K array that showed
clustered (multimodal) beta-value distributions across 590 samples. The
vast majority (9,199) were attributed to underlying SNPs in the probe
binding site. Gap probes showed population-stratified patterns consistent
with ancestry-informative SNPs.

**Relevance.** This IS a wet-lab observation of probe behavior affected by
sequence variants. The "gap" signal pattern is a direct phenotype of
hybridization disruption caused by SNPs under the probe. However, the study
does not formally validate in silico SNP predictions against the observed
gap patterns -- it works from observed data backwards to infer SNP effects.

#### 2.1.4 In Silico PCR Primer Validation (Microbiology)

- Zagordi et al. (2024, Frontiers in Bioinformatics): "In silico PCR
  analysis: a comprehensive bioinformatics tool for enhancing nucleic acid
  amplification assays."
- AOAC guidelines (2020, JAOAC 103:882): "Recommendations for Developing
  Molecular Assays for Microbial Pathogen Detection Using Modern In Silico
  Approaches."

These establish that in silico primer/probe evaluation is standard practice
in microbiology assay development. Computational screening of primer binding
sites against genome databases is widely accepted as a necessary (but not
sufficient) step. However, wet-lab validation is still required.

### 2.2 What Does NOT Exist in the Published Literature

Based on extensive searching, the following validation studies do NOT appear
to exist:

1. **No study has taken pangenome-predicted probe absence/divergence for
   FISH BAC probes and validated it against clinical FISH performance
   differences across populations.**

2. **No study has computationally predicted MLPA probe failure due to
   population-specific SNPs at ligation junctions and then confirmed allele
   dropout in clinical MLPA testing across populations.**

3. **No study has used pangenome graphs to predict NGS capture efficiency
   differences across populations and validated with actual capture metrics
   (on-target rate, coverage uniformity) stratified by ancestry.**

4. **No study has systematically predicted CMA probe signal intensity
   differences across populations using pangenome data and validated against
   clinical CMA signal profiles.**

5. **No study has used OGM label density predictions from pangenome
   sequences to predict population-specific OGM performance differences.**

### 2.3 The Minimum Viable Validation Experiment

Given the evidence gaps, the smallest feasible validation would be:

**For FISH (highest alignment between pipeline and biophysics):**
- Select 10-20 FISH probes predicted by the pipeline to have
  population-correlated mapping disparities (e.g., probe maps to all NFE
  haplotypes but is absent in >20% of AFR haplotypes).
- Perform FISH on cell lines or clinical specimens from individuals of
  matching ancestry.
- Score signal presence/absence/split and compare to pipeline predictions.
- This is feasible because FISH signal is binary-ish (present/absent/split)
  and the pipeline's coordinate-level prediction is well-matched.

**For MLPA (highest clinical risk from the gap):**
- Identify MLPA probes where the pipeline predicts universal mapping but
  population-specific SNPs exist at ligation junctions (requires sequence
  extraction, not just coordinate projection).
- Test those probes against DNA samples from variant carriers.
- This would require extending the pipeline to sequence-level analysis.

---

<a name="task-3"></a>
## Task 3: What Coordinate-Level Projection Can and Cannot Tell Us

### 3.1 What `odgi position` Actually Computes

Based on the odgi documentation (pangenome.github.io, odgi.readthedocs.io)
and the ODGI paper (Guarracino et al., 2022, Bioinformatics 38:3319):

`odgi position` performs coordinate translation between paths embedded in a
pangenome graph. When the pipeline calls:

```
odgi position -i chr17.og -b input.bed -r HG00272#1#chr17
```

The tool does the following:

1. Takes a BED interval on the source path (e.g., GRCh38#0#chr17, 100000,
   200000).
2. Identifies which graph nodes are traversed by the source path in that
   coordinate range.
3. Performs a breadth-first search (BFS) from those nodes to find where the
   target path (e.g., HG00272#1#chr17) traverses nearby nodes.
4. The BFS is limited to a configurable search radius (default: 10,000 bp).
5. Reports the target path coordinates corresponding to the query interval.

**The critical question: does a successful mapping mean the sequence is
identical?**

In a pangenome graph built by Minigraph-Cactus, the graph topology encodes
sequence variation:

- **Shared nodes**: When two haplotype paths traverse the same node, the
  underlying sequence is IDENTICAL at that node. Nodes in a variation graph
  are labeled with DNA sequences. If path A and path B both pass through
  node 42, they share exactly the sequence stored in node 42.

- **Bubble/snarl regions**: Where haplotypes differ, they traverse different
  nodes within a bubble structure. One path goes through node 42 -> 43 -> 45
  while another goes 42 -> 44 -> 45. Nodes 43 and 44 contain different
  sequences (the variant alleles).

- **What `odgi position` reports**: When it successfully maps a source
  interval to a target path, it means the BFS found the target path within
  the search radius. Crucially:
  - If the source and target paths share the SAME nodes, the sequence is
    identical, and the mapping is reliable.
  - If the target path traverses DIFFERENT nodes in a bubble but re-joins
    downstream, `odgi position` may still report a mapping, but the
    underlying sequence at the divergent nodes is DIFFERENT.
  - The tool does NOT report whether the mapped region traverses the same
    or different nodes. It reports coordinates, not sequence identity.

- **What it does NOT compute**: `odgi position` does not perform sequence
  alignment. It does not compute edit distance, percent identity, or
  mismatch count. It does not extract the sequence at the target coordinates.
  It is a coordinate translation tool, not an alignment tool.

### 3.2 Implications for the Pipeline

The pipeline's `OdgiCliProjector` assigns `mapping_quality = 1.0` and
`mapping_type = "exact"` for ANY successful `odgi position` output. The
code in `projector.py` (lines 478-484, 583-593) shows that if odgi returns
parseable coordinates, the result is classified as mapped with quality 1.0.

The `MissingnessCategory` enum in `haplotype_coverage.py` defines DIVERGENT
as requiring `mapping_quality < threshold`, but the code documentation
explicitly acknowledges (lines 188-196): "The current OdgiCliProjector
produces binary quality values (1.0 for mapped, 0.0 for unmapped). DIVERGENT
will only be classified when the projector is enhanced to produce graded
quality scores."

**This means the pipeline currently has exactly two states:**
- MAPPED (odgi returned coordinates): classified as successful, quality 1.0
- ABSENT (odgi failed or returned nothing): classified as absent, quality 0.0

There is no intermediate state. The DIVERGENT, DUPLICATED, and REARRANGED
categories exist in the code but are unreachable with the current projector.

### 3.3 Alternative Tools for Graded Mapping Quality

| Tool | What It Does | Mismatch Scoring | Pros | Cons |
|------|-------------|------------------|------|------|
| **odgi position** (current) | Coordinate translation between paths via BFS | None -- binary mapped/unmapped | Fast, native to .og format, batch BED input | No sequence-level information |
| **vg giraffe** | Short-read alignment to pangenome graph (GBZ format) | Full scoring: +1 match, -4 mismatch, -6/-1 gap | Produces mapping quality scores, identity metrics, designed for Illumina-length reads (25-150bp) | Requires GBZ index; designed for reads not arbitrary intervals; overhead for probe-level queries |
| **vg map** | Original vg graph aligner | Full dynamic programming alignment | Most flexible, handles arbitrary sequences | Slow for large graphs |
| **GraphAligner** | Long-read to graph alignment | Edit distance, banded DP | Fast, handles long sequences, low memory | Designed for noisy long reads; may need parameter tuning for short exact probes |
| **vg inject + vg surject** | Convert BAM alignments to graph space and back | Inherits aligner scoring | Can project existing alignments | Requires pre-aligned reads |
| **Sequence extraction + BLAST/Tm calculation** | Extract haplotype sequence at projected coordinates, then align probe to extracted sequence | Full thermodynamic modeling possible | Most biologically accurate; can compute Tm, delta-G | Requires sequence extraction step; computationally heavier |

### 3.4 What Would Need to Change in the Pipeline

To go beyond binary mapping and predict hybridization success:

**Option A: Lightweight (sequence extraction + mismatch counting)**
1. When `odgi position` returns coordinates on a target path, extract the
   underlying sequence at those coordinates using `odgi extract` or by
   reading the graph nodes.
2. Align the probe sequence to the extracted target sequence (simple
   pairwise alignment for short probes, or just count mismatches).
3. Report mismatch count, percent identity, and estimated Tm shift.
4. This would enable the DIVERGENT classification for probes that map
   but have sequence differences.

**Option B: Moderate (vg giraffe for short probes)**
1. Convert probe sequences to FASTQ format.
2. Align each probe to the pangenome using vg giraffe (requires GBZ index).
3. Parse mapping quality and identity from GAM/GAF output.
4. This would give per-probe, per-haplotype alignment scores.
5. Well-suited for 25-120 nt probes (CMA, MLPA, NGS baits).

**Option C: Full thermodynamic modeling**
1. Extract haplotype sequences as in Option A.
2. Apply nearest-neighbor thermodynamic models (SantaLucia, 1998) to
   compute predicted Tm and delta-G for each (probe, haplotype) pair.
3. Model hybridization success as a function of the free energy difference
   between the probe-haplotype duplex and the probe-reference duplex.
4. Most biologically accurate but computationally expensive (6.85M probes
   x 272 haplotypes = 1.86B calculations for CMA alone).

---

<a name="task-4"></a>
## Task 4: Estimate the Detectable vs Undetectable Gap

### 4.1 What the Pipeline DETECTS (True Positives for Probe Failure)

The pipeline detects events where a haplotype path is completely absent
through the probe's reference coordinate range. This corresponds to:

- **Large deletions** that eliminate the entire probe target region on a
  haplotype. For FISH probes (100-200 kb), this requires a deletion
  comparable in size to the probe itself.
- **Large inversions** that reorient the region, potentially preventing
  path traversal in the expected orientation.
- **Complex rearrangements** (translocations, large insertions disrupting
  the region) that break the path continuity.
- **Haplotype-specific structural absence** -- regions present in the
  reference but absent from specific haplotypes (e.g., large population-
  specific CNVs).

These are genuine, clinically significant events. The HPRC pangenome
contains 119 Mb of euchromatic polymorphic sequence relative to GRCh38,
with ~90 Mb from structural variation. The HPRC cataloged 167,291 SV sites
including 65,075 deletions, 74,125 insertions, and 25,371 complex sites.
One hundred seventeen loci show evidence of population stratification.

### 4.2 What the Pipeline MISSES (False Negatives)

The pipeline misses events where the haplotype path exists but the sequence
differs from the reference. This includes:

- **SNPs within probe binding sites.** At the human intra-species SNP rate
  (~1 per 1,000 bp average, higher in diverse populations), a 25-mer probe
  has approximately a 2.5% chance of overlapping at least one SNP. For the
  6.85M probes on the CMA, this represents ~170,000 probes potentially
  affected by at least one SNP on at least one haplotype.
- **Small indels (1-50 bp)** within probe footprints. These may be
  represented as bubbles in the graph, and `odgi position` may still map
  through the region despite the sequence change.
- **Clustered variants** in highly polymorphic regions (HLA, KIR, etc.) that
  cause multiple mismatches within a single probe.
- **Population-specific common SNPs** at MLPA ligation junctions -- the
  highest-risk failure mode that is completely invisible to coordinate
  projection.

### 4.3 Estimated Proportions by Technology

| Technology | Probe Size | Failure Mode | Pipeline Detects | Pipeline Misses | Estimated Detection Rate |
|-----------|-----------|-------------|-----------------|----------------|------------------------|
| **FISH** | 100-200 kb | Large SV (deletion, inversion) | Yes -- path absence | SNP-level (irrelevant for BAC FISH) | **>90%** of clinically significant failures |
| **CMA (25-mer)** | 25 bp | SNP under probe | No | Yes -- all SNP-based failures | **<10%** of per-probe failures; possibly more for cluster-level effects from large SVs |
| **CMA (60-mer)** | 60 bp | SNP cluster, small indel | No | Yes -- but less clinically relevant per probe | **<10%** of per-probe failures |
| **MLPA** | ~50 nt | SNP at ligation junction | No | Yes -- the most dangerous failures | **<5%** unless large SV eliminates entire exon |
| **NGS capture** | 80-120 nt | Large SV, high-divergence region | Partially | Gradual efficiency loss | **30-50%** (captures large absences, misses efficiency degradation) |
| **OGM** | N/A (6-bp motif) | Label pattern change from SV | Yes for large SV | Individual motif gain/loss | **>80%** for SV detection; motif-level changes are below resolution |

### 4.4 Is the Detectable Fraction Still Scientifically Interesting?

**Yes, particularly for FISH.** The pipeline's strength aligns precisely
with the technology where coordinate-level projection is most predictive.
FISH is the largest technology category in the profile catalog (319
profiles), and FISH probe failure due to structural variation is exactly the
kind of event that: (a) the pipeline can detect, (b) has clear clinical
consequences (false negative FISH result), and (c) is plausible given known
population-level structural variation in clinically tested regions.

For other technologies, the detectable fraction is smaller but not
negligible:

- **Large structural variants that eliminate entire gene targets** would
  cause failure across ALL probe-based technologies. If a large deletion
  removes the target exon of an MLPA probe, the pipeline detects this
  regardless of the MLPA-specific biophysics.
- **The HPRC contains 117 population-stratified SV loci.** If any of these
  overlap clinical assay target regions, the pipeline can identify them.

The honest framing is: **the pipeline is a structural variant detector
applied to clinical assay coordinates**, not a hybridization simulator. Its
results are strongest for technologies where structural variation is the
primary failure mode (FISH, OGM) and weakest for technologies where
single-nucleotide variation drives failure (MLPA, CMA).

---

<a name="deliverable-1"></a>
## Deliverable 1: Technology-by-Technology Evidence Table

| Technology | Probe Size | Mismatch Tolerance | Primary Failure Mechanism | Coordinate Projection Detects? | Evidence Strength |
|-----------|-----------|-------------------|--------------------------|-------------------------------|------------------|
| **CMA (Affymetrix 25-mer)** | 25 bp | 1 SNP = 40-60% signal reduction; 2-3 SNPs = near background | SNP under probe reduces PM/MM ratio; clustered SNPs cause focal dropout | **No** -- path exists, sequence differs. Only large SVs deleting probe clusters detectable. | Strong empirical evidence for mismatch effects (Binder 2009, Naef 2003). No pangenome validation. |
| **CMA (Agilent 60-mer)** | 60 bp | 1 SNP = ~40% of PM/MM differential retained; 3-5 SNPs for full loss | Similar to 25-mer but more tolerant; small indels more problematic | **No** -- same blind spot but lower clinical risk per probe. | Moderate evidence (Kreil 2008, He 2005). No pangenome validation. |
| **MLPA** | ~50 nt (2x ~25 nt hemiprobes) | 1 SNP at ligation junction = complete signal loss | SNP at ligation site prevents ligase activity; allele dropout mimics deletion | **No** -- path exists, single nucleotide change at junction invisible. Most dangerous blind spot. | Strong evidence from MRC Holland, multiple clinical reports of false positives from SNPs. |
| **FISH** | 50 kb - 1 Mb (BAC/fosmid) | Tolerates ~5-10% divergence; within-species SNPs irrelevant | Large SV (deletion, inversion) removing/rearranging probe target | **Yes** -- this is precisely what path absence/rearrangement means. | Strong biophysical basis. No population-specific validation but mechanism is clear. |
| **NGS capture** | 80-120 nt (RNA/DNA baits) | 10-20% divergence tolerated | Gradual efficiency loss with divergence; complete failure only with target deletion | **Partially** -- detects complete target absence but not efficiency gradation. | Strong cross-species evidence (Bi 2012, NAR 2016). Capture bias with genotype documented (mouse exome). |
| **OGM** | N/A (DLE-1 labels CTTAAG) | N/A -- not a hybridization technology | Label pattern changes from large SVs; individual motif gain/loss below resolution | **Mostly** -- large SVs detectable; motif-level changes not relevant to coordinate projection. | Moderate evidence. DLE-1 mechanism well-documented. Population-specific label density not studied. |

---

<a name="deliverable-2"></a>
## Deliverable 2: Validation Precedent Summary

### What Exists

1. **SARS-CoV-2 primer/probe failure prediction.** In silico identification
   of binding site variants validated by in vitro RT-qPCR sensitivity loss.
   Demonstrates the principle for short PCR-based probes. (Multiple studies,
   2020-2022.)

2. **Methylation array probe ambiguity from pangenome/T2T.** Heiss et al.
   (2025, Cell Reports) aligned probe sequences to T2T-CHM13 and the HPRC
   pangenome, identifying population-specific unambiguous probe sets.
   Closest precedent to this project but validates computational mapping,
   not wet-lab hybridization outcomes.

3. **Methylation array "gap probes" from SNPs.** Andrews et al. (2016)
   showed that SNPs under 450K array probes produce characteristic
   multimodal beta-value distributions observable in real data, with
   population-stratified patterns. This is observational wet-lab evidence
   consistent with in silico predictions but not a formal validation study.

4. **Reference genome bias in variant calling.** Wang et al. (2022) and
   the HPRC consortium have documented that non-European genomes have
   higher VUS burden, lower diagnostic rates, and more alignment artifacts
   when analyzed against GRCh38. This establishes the upstream premise
   (reference bias exists) but does not validate the downstream claim
   (probe-level hybridization failure).

### What Does NOT Exist

1. No published validation of pangenome-predicted FISH probe failure against
   clinical FISH signal in diverse populations.
2. No published validation of pangenome-predicted MLPA allele dropout
   against clinical MLPA results stratified by ancestry.
3. No published validation of pangenome-predicted NGS capture bias against
   on-target rates stratified by ancestry.
4. No published validation of pangenome-predicted CMA signal dropout against
   clinical CMA intensity profiles stratified by ancestry.
5. No systematic cross-technology validation study linking in silico
   pangenome analysis to in vitro assay performance.

### Minimum Viable Validation

**Tier 1 (most feasible, strongest alignment with pipeline):**
- Select 10-20 FISH probes with predicted population-correlated path
  absence.
- FISH on HapMap/1000G cell lines from matching populations.
- Score signal presence/absence. Binary outcome matches pipeline prediction.
- Estimated cost: moderate (cell lines, FISH reagents, microscopy time).

**Tier 2 (harder but higher impact):**
- For MLPA: Identify probes where pipeline predicts universal mapping but
  sequence-level analysis reveals common SNPs at ligation junctions in
  specific populations.
- Test MLPA on DNA from SNP carriers vs. non-carriers.
- Score for allele dropout (apparent deletion in carrier).
- Requires extending the pipeline to sequence-level analysis first.

**Tier 3 (comprehensive but expensive):**
- For CMA: Run clinical CMA on DNA from HPRC samples (or matching
  populations). Compare per-probe signal intensities to pipeline predictions.
- For NGS: Run capture sequencing on diverse samples. Compare on-target
  rates and coverage uniformity to pipeline-predicted probe mapping status.

---

<a name="deliverable-3"></a>
## Deliverable 3: Tool Comparison

### odgi position (Current Pipeline Tool)

- **What it does:** Translates coordinates between paths in an odgi (.og)
  pangenome graph using BFS within a configurable search radius (default
  10 kb).
- **Output:** Target path coordinates for the query interval, or failure
  if the target path is not reachable.
- **Sequence information:** None. Reports coordinates only.
- **Strengths:** Native to .og format; fast BED batch mode; well-matched
  to the existing pipeline architecture; sufficient for detecting large
  structural events (path absence, rearrangement).
- **Weaknesses:** Binary output (mapped/unmapped) with no graded quality.
  Cannot distinguish "mapped to identical sequence" from "mapped to nearby
  divergent sequence." The BFS walk may bridge small bubbles without
  reporting them.
- **Verdict for probe prediction:** Adequate for FISH-scale analysis where
  the relevant failure mode is structural. Inadequate for CMA/MLPA/NGS
  where sequence-level mismatches matter.

### vg giraffe (Recommended Alternative for Short Probes)

- **What it does:** Aligns short sequences (designed for Illumina reads,
  25-300 bp) to a pangenome graph in GBZ format.
- **Output:** Alignments in GAM/GAF format with mapping quality, alignment
  score (+1 match, -4 mismatch, -6/-1 gap), and identity metrics.
- **Sequence information:** Full alignment with mismatch positions.
- **Strengths:** Designed for exactly the sequence lengths of CMA (25-60
  bp), MLPA (~50 nt), and NGS capture probes (80-120 nt). Produces
  graded quality scores. Maps to specific haplotype paths.
- **Weaknesses:** Requires GBZ index (the pipeline already has a .gbz
  file). Designed for read mapping, not coordinate projection -- would
  require reformulating probes as query sequences rather than coordinate
  intervals. Overhead per query may be higher than odgi position.
- **Verdict:** The most promising tool for adding sequence-level probe
  quality assessment. Would enable detection of mismatches, gaps, and
  identity scores for each (probe, haplotype) pair.

### GraphAligner

- **What it does:** Aligns long sequences to genome graphs using seeded
  and banded dynamic programming.
- **Output:** Alignments with edit distance.
- **Sequence information:** Full alignment.
- **Strengths:** Fast, low memory. Handles sequences of any length.
  Could align BAC-scale probe sequences (though this is unusual).
- **Weaknesses:** Optimized for noisy long reads, not short exact probes.
  May need parameter tuning for 25-60 nt sequences with low error rates.
- **Verdict:** Less appropriate than vg giraffe for the short-probe use
  case. Could be useful if the project needs to analyze the full probe
  sequence rather than just the binding site.

### Sequence Extraction + Pairwise Alignment

- **What it does:** Uses `odgi extract` or `odgi paths` to extract the
  haplotype sequence at the projected coordinates, then performs standard
  pairwise alignment (Needleman-Wunsch, BLAST, or simple mismatch count).
- **Output:** Percent identity, mismatch count, gap count, optionally Tm
  prediction using nearest-neighbor thermodynamic models.
- **Sequence information:** Complete.
- **Strengths:** Most biologically accurate. Can feed directly into
  thermodynamic modeling (delta-G, Tm prediction) for each (probe,
  haplotype) pair. Technology-specific: can apply different thresholds
  for FISH (irrelevant), CMA 25-mer (1 mismatch matters), MLPA (junction
  position matters), NGS (10-20% divergence tolerated).
- **Weaknesses:** Two-step process (project + extract + align).
  Computationally heavier. For CMA (6.85M probes x 272 haplotypes),
  this is ~1.86 billion extractions.
- **Verdict:** The gold standard for accuracy. Could be applied selectively
  to probes flagged as "mapped" to determine whether the mapping is
  sequence-identical or divergent. A hybrid approach (odgi position for
  coarse filtering, then sequence extraction for flagged probes) would
  balance accuracy and computational cost.

---

<a name="deliverable-4"></a>
## Deliverable 4: Honest Assessment of the Causal Chain

### The Causal Chain

The project's implicit causal chain is:

```
Probe fails to map in pangenome graph
  => Probe target is structurally absent on that haplotype
    => Probe will fail to hybridize in vitro
      => Clinical assay produces false negative
        => Patient harm (missed diagnosis)
```

### Assessment by Technology

#### FISH: STRONG (7/10)

The causal chain for FISH is the strongest because:
- The failure mode (large structural absence) IS what coordinate projection
  detects.
- BAC FISH probes are insensitive to SNP-level variation, so the
  undetectable gap is small.
- FISH is the largest technology category (319 profiles).
- Clinical FISH false negatives from probe design issues are plausible and
  consequential (e.g., missed BCR-ABL fusion).

Weaknesses:
- No published validation against wet-lab FISH in diverse populations.
- The magnitude of the effect is unknown (how many FISH probes actually fail
  on real haplotypes?).
- Pipeline is blocked on graph data; no results yet.

#### OGM: MODERATE (5/10)

The causal chain for OGM is moderate because:
- OGM does not use hybridization probes, so the "probe failure" framing
  is partially mismatched.
- Coordinate projection can detect structural variants that would also be
  detected by OGM, creating a circular argument (the pipeline detects what
  OGM itself would detect).
- The equity question for OGM is about label density patterns, not probe
  hybridization.

#### NGS Capture: MODERATE (4/10)

The causal chain for NGS capture is moderate because:
- Coordinate projection catches complete target absences (real but rare).
- The main concern (capture efficiency degradation from mismatches) is not
  detectable by the current pipeline.
- NGS capture's high mismatch tolerance (10-20%) means that within-species
  variation rarely causes outright failure.
- The clinical impact of slightly reduced capture efficiency (lower coverage,
  not absent coverage) is less severe than for other technologies.

#### CMA: WEAK (2/10)

The causal chain for CMA is weak because:
- The primary failure mode (SNP under 25-mer probe) is invisible to
  coordinate projection.
- CMA uses massive probe redundancy (6.85M probes), so individual probe
  failures are mitigated by averaging across adjacent probes.
- The clinically relevant failure would be a region where MANY adjacent
  probes fail simultaneously, which would require a large structural event
  -- bringing us back to what coordinate projection can already detect.
- For CMA, the real question is whether large structural absence creates
  focal probe dropouts, not whether individual SNPs affect individual
  probes. The pipeline CAN address this question.

#### MLPA: WEAK-TO-ABSENT (1/10)

The causal chain for MLPA is the weakest because:
- The primary failure mode (SNP at ligation junction) is completely
  invisible to coordinate projection.
- This is the most clinically dangerous failure mode (complete allele
  dropout from a single SNP, mimicking a pathogenic deletion).
- MLPA has only 6 profiles in the catalog, limiting the scope of analysis.
- The pipeline would need to be fundamentally extended (to sequence-level
  alignment with ligation-junction annotation) to address MLPA.

### Overall Assessment

**The pipeline, as designed, is a structural variant screen of clinical
assay probe coordinates.** Its strength is detecting cases where
population-specific structural variation eliminates or rearranges the target
of a clinical probe. This is a scientifically interesting and clinically
relevant question, particularly for FISH.

**The pipeline is NOT a hybridization simulator.** It cannot predict
whether a probe that maps (coordinates exist on the haplotype) will
actually hybridize successfully. This distinction must be stated clearly
in any publication.

**The causal chain has two segments:**

1. **"Path absent in graph" => "Probe target structurally absent"**: This is
   a valid inference with high confidence. If the haplotype path does not
   traverse the probe region, the target DNA is structurally different. This
   segment of the chain is solid.

2. **"Probe target structurally absent" => "Clinical assay fails"**: This
   varies by technology. For FISH (where the entire probe target must be
   present for signal), it is strong. For CMA (where probe redundancy
   compensates), it is weaker. For MLPA and other technologies, structural
   absence of the entire target is not the primary concern -- SNP-level
   variation is.

**The project's scientific contribution is real but bounded.** It answers
the question: "Are there population-correlated structural variants in the
human pangenome that overlap clinical assay probe target regions?" This is
a valuable and novel contribution. The additional claim that "these overlaps
predict clinical assay failure" is supported for FISH but requires
qualification for other technologies.

---

## Key References

### Probe Biophysics
- Binder et al. (2009). Mismatch and G-Stack Modulated Probe Signals on SNP Microarrays. PLOS ONE 4:e7862.
- Kreil et al. (2008). Design and analysis of mismatch probes for long oligonucleotide microarrays. BMC Genomics 9:491.
- Schouten et al. (2002). Relative quantification of 40 nucleic acid sequences by multiplex ligation-dependent probe amplification. NAR 30:e57.
- SantaLucia (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. PNAS 95:1460.
- Bi et al. (2012). Sequence capture by hybridization to explore modern and ancient genomic diversity. NAR 44:4504.

### Validation Precedents
- Heiss et al. (2025). Complete reference genome and pangenome improve genome-wide detection and interpretation of DNA methylation. Cell Reports.
- Andrews et al. (2016). "Gap hunting" to characterize clustered probe signals in Illumina methylation array data. Epigenetics & Chromatin.
- Vanaerschot et al. (2020). Presence of mismatches between diagnostic PCR assays and coronavirus SARS-CoV-2 genome. Royal Society Open Science 7:200636.
- Khan & Cheung (2022). Identification of mutations in SARS-CoV-2 PCR primer regions. Scientific Reports 12:17953.

### Tools and Pangenome
- Guarracino et al. (2022). ODGI: understanding pangenome graphs. Bioinformatics 38:3319.
- Siren et al. (2021). Pangenomics enables genotyping of known structural variants in 5202 diverse genomes. Science.
- Hickey et al. (2023). Pangenome graph construction from genome alignments with Minigraph-Cactus. Nature Biotechnology.
- Rautiainen & Marschall (2020). GraphAligner: rapid and versatile sequence-to-graph alignment. Genome Biology 21:253.

### Reference Bias and Population Genomics
- Liao et al. (2023). A draft human pangenome reference. Nature 617:312.
- Ebler et al. (2022). Pangenome-based genome inference. Nature Methods 19:1516.
- Manso et al. (2025). Structural variation in 1,019 diverse humans based on long-read sequencing. Nature.

### Source URLs
- [odgi position documentation](https://odgi.readthedocs.io/en/latest/rst/commands/odgi_position.html)
- [odgi position (pangenome.github.io)](https://pangenome.github.io/odgi.github.io/rst/commands/odgi_position.html)
- [ODGI paper (Bioinformatics)](https://academic.oup.com/bioinformatics/article/38/13/3319/6585331)
- [vg Giraffe wiki](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe)
- [GraphAligner paper](https://link.springer.com/article/10.1186/s13059-020-02157-2)
- [Heiss et al. 2025 (Cell Reports)](https://www.sciencedirect.com/science/article/pii/S2211124725005261)
- [Andrews et al. 2016 (gap hunting)](https://link.springer.com/article/10.1186/s13072-016-0107-z)
- [MRC Holland MLPA technique](https://www.mrcholland.com/technology/mlpa/technique)
- [Binder et al. 2009 (PLOS ONE)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0007862)
- [Cross-species capture (NAR 2016)](https://academic.oup.com/nar/article/44/10/4504/2516864)
- [SARS-CoV-2 primer mismatches](https://royalsocietypublishing.org/doi/10.1098/rsos.200636)
- [Pangenome and methylation (bioRxiv)](https://www.biorxiv.org/content/10.1101/2024.10.07.617116v1.full)
- [HPRC SV in 1019 genomes (Nature 2025)](https://www.nature.com/articles/s41586-025-09290-7)
- [Bionano OGM overview](https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/optical-genome-mapping)
- [OGM in routine diagnostics (Genes 2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8701374/)
- [Minigraph-Cactus (Nature Biotechnology)](https://www.nature.com/articles/s41587-023-01793-w)
- [Pangenome reference bias (Frontiers in Genetics 2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12492951/)
- [HPRC draft pangenome (Nature 2023)](https://www.nature.com/articles/s41586-023-05896-x)

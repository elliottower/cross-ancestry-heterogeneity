# Pre-registration: Vector-Valued Sheaf Q Extension

**Freeze date: 2026-07-17**
**Frozen before any vector Q computation is run on real data.**
**Frozen before the author has viewed any results from the analysis scripts.**

## Overview

We extend the scalar sheaf Q validation (71 MR pairs, 84.5% accuracy) to
vector-valued stalks. On vector stalks the sheaf-Cochran equivalence breaks:
the multivariate Q_V captures cross-component covariance structure that d
independent scalar Q tests cannot. The specific failure mode is *effect vector
rotation* — per-trait effect magnitudes are preserved across ancestries but the
joint direction of the effect vector changes.

## Data sources

1. **Pan-UK Biobank multi-ancestry GWAS** (pan.ukbb.broadinstitute.org):
   per-ancestry summary statistics for blood cell traits across EUR, AFR, CSA,
   EAS, AMR, MID. Used for multi-trait vector Q at known loci.

2. **CAMeRa multi-ancestry MVMR** (github.com/yoonsucho/camera_analysis):
   per-ancestry multi-exposure MR estimates for cardiometabolic outcomes.
   Already used for LDL→stroke scalar analysis.

3. **JASS multi-ancestry joint analysis** (PMC11022331 / biorxiv
   2023.06.23.546248): 19 hematological + glycemic traits jointly analyzed
   across 5 ancestries. If per-ancestry summary statistics are publicly
   available; otherwise substitute Pan-UKB blood traits.

**Primary data source**: Pan-UKB blood traits (source 1), as JASS
per-ancestry summary statistics are not currently publicly available.
CAMeRa (source 2) is used for MVMR-specific prediction H_V4 only.

## Predictions

### H_V1: ACKR1/Duffy effect vector rotation

At the ACKR1/Duffy locus (rs2814778 or LD proxy r² > 0.8 in EUR), the
multivariate sheaf Q_V computed on >= 10 blood cell traits across >= 3
ancestry groups will be significant at p < 0.001, AND the maximum pairwise
rotation angle θ_max between ancestry effect vectors will exceed 30°.

**Rationale**: The Duffy-null allele (fixed in sub-Saharan African populations
due to malaria selection) eliminates ACKR1-mediated chemokine scavenging,
producing a coordinated multi-trait shift (neutrophil count, monocyte count,
basophil count, WBC count) that is qualitatively different from the EUR/EAS
effect profile — a rotation, not a uniform scaling.

### H_V2: Rotation-only loci exist

Among genome-wide significant loci (p < 5e-8 in any ancestry) tested across
>= 3 ancestries on >= 5 blood cell traits, at least 5% of loci with
significant vector Q_V (at Bonferroni-corrected threshold) will be
"rotation-only": vector Q_V significant but NO individual scalar Q_j
significant after Bonferroni correction for d traits.

**Rationale**: Pleiotropic loci affecting multiple blood traits through shared
regulatory pathways can have ancestry-specific pathway weighting (due to LD,
allele frequency, or gene-environment interaction) that rotates the effect
vector without changing per-component magnitudes enough for any single scalar
test to detect.

### H_V3: Vector Q exceeds sum of scalar Q

Across all loci where vector Q_V is significant, the median ratio
Q_V / sum(Q_j^scalar, j=1..d) will exceed 1.0, indicating that the vector
test captures covariance-mediated heterogeneity beyond the sum of per-trait
heterogeneity.

**Rationale**: When cross-trait covariance is non-trivial (correlated blood
traits sharing regulatory pathways), the off-diagonal blocks of the block
covariance matrix contribute to Q_V. If heterogeneity is distributed across
components in a correlated pattern (rotation), Q_V exceeds the sum of marginal
Q tests.

### H_V4: CAMeRa MVMR vector-scalar divergence

For the lipid triad (LDL + HDL + triglycerides) → CAD MVMR model across >= 3
ancestries, vector Q_V will be significant at p < 0.05 while at least one
individual exposure's scalar Q is non-significant. The looser threshold
(p < 0.05 vs p < 0.001 for H_V1) reflects the lower statistical power of
MVMR estimates, which have fewer ancestries and higher variance than
single-trait GWAS.

**Rationale**: The relative contribution of LDL, HDL, and triglycerides to CAD
risk differs across ancestries due to population-specific dietary patterns,
metabolic regulation, and LD structure at lipid-associated loci. The joint
effect structure (pathway weighting) varies even when individual exposure
effects are comparable.

### H_V5: ACKR1 selective scalar detection

At the ACKR1 locus, scalar Q on neutrophil count and monocyte count will
individually reach significance (p < 0.05) across ancestries, but scalar Q on
>= 5 of the tested traits will be non-significant (p > 0.05).

**Rationale**: ACKR1/Duffy affects a specific subset of blood cell lineages
(myeloid-derived cells). Traits not on the ACKR1 pathway (platelet count,
red cell indices, glycemic traits) should show no ancestry-specific
heterogeneity at this locus. This creates a mixed pattern where some scalar
tests fire and others don't — the vector Q_V integrates the coordinated
multi-trait signal more powerfully than any single scalar test.

### H_V6: Anti-correlated cancellation is rare

The "component-only" category (scalar Q significant on >= 1 trait but vector
Q_V non-significant) will account for < 2% of all heterogeneous loci.

**Rationale**: Anti-correlated heterogeneity that cancels in the joint test
requires a specific and biologically unusual pattern: one trait's
cross-ancestry heterogeneity must exactly oppose another's in the covariance-
weighted sense. This is possible in principle but should be rare in practice.

## Analysis plan

1. Download per-ancestry summary statistics for d >= 10 blood cell traits at
   genome-wide significant loci across >= 3 ancestries.
2. At each locus, compute:
   a. d independent scalar Cochran's Q tests (chi²_{K-1} each)
   b. Multivariate sheaf Q_V (chi²_{d(K-1)})
   c. Pairwise rotation angles θ_{ij} between ancestry effect vectors
3. Classify each locus as: concordant-significant, concordant-null,
   rotation-only, or component-only.
4. Test H_V1-H_V6.
5. For CAMeRa MVMR: extract multi-exposure per-ancestry estimates, compute
   scalar Q per exposure and vector Q_V on the joint model.

## Outcome

Fill TODO placeholders in paper_A_v8_eje.tex Section 5 (Vector-valued sheaf Q)
with the real numbers from this analysis.

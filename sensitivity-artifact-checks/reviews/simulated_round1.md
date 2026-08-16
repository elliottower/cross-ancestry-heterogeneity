# Simulated referee round 1 — paper_v8 + supplementary_v5 (2026-08-16)

Three independent referee agents (statistical methods; applied multi-ancestry GWAS;
transparency/reporting), each reading the full manuscript and supplement as a Genetic
Epidemiology submission. All three returned MAJOR REVISION. Full reports preserved in the
session record; consolidated findings below.

## Convergent (all three referees)

1. Ascertainment controls (non-EUR-only; discovery-traits-excluded) were run on the
   24-trait panel; abstract and Discussion attach them to the 14-trait corrected estimate.
   Rerun on the 14-trait panel; add a base-panel column to Table 2.
2. H_R3 covariance-stability gate scored on the full-panel share (5.3 pts, held) while the
   primary-panel AFR shift is 29.2 pts (66.4% -> 37.2%, 506 -> 148 significant). Verdict
   must be re-scored on the primary quantity; majority claim conditioned on the covariance
   model wherever stated.
3. ACKR1 contradiction: pipeline's own region locus chr1_159175494_C_T is multivariate-only
   (theta_max 102.5) while §3.6 presents hand-queried rs2814778 (chr1:159174683, GRCh37)
   as concordant-heterogeneous. Must be reconciled explicitly; it is itself a demonstration
   of index-choice sensitivity.

## Major, individual

- R1: Table 3 Q_V p-value computed at df=125 (pilot 25-trait config) against a 24-trait
  caption (df=120, p ~ 3e-57). 14-trait panel spectrum never reported; WBC ~ sum of four
  retained leukocyte counts is a surviving near-identity; consider dropping WBC or adding
  basophil count. Diagonal statistic used inferentially despite failing the calibration
  gate (0.139) — recalibrate by moment matching before quoting 231/10.0%. Tail calibration
  needed at BH-relevant alphas (0.01, 0.001), not only 0.05. Rotation mixes panels
  (24-trait vectors at 14-trait-selected loci) and is scale-dependent. Multiplicity regimes
  differ between arms (BH across loci vs within-locus Bonferroni); report a matched regime.
  Correlation-matrix uncertainty (0.094 element-wise between EUR windows) unpropagated;
  genome-wide LD-pruned estimate should replace the single 5-Mb window.
- R2: Parametric bootstrap is circular (same Sigma generates and tests); empirical
  calibration at random null variants is feasible and decisive. Per-ancestry correlation
  matrices feasible (AFR estimable from 15,967 variants) and should be primary. IRNT
  per-ancestry standardization can manufacture the rotation signature via SD ratios;
  transformation never stated; diagonal-rescaling absorption test needed. Non-EUR 89.8%
  rejection with only small groups remaining is itself evidence of covariance misfit.
  No positive control exists for the multivariate-only category; spike-in simulations
  under real Sigma structure needed. Provide per-locus supplementary data. Missing
  literature: MR-MEGA, Popcorn, S-LDXR, Hou et al., CPASSOC/metaUSAT.
- R3: Methods §2.8 sentence ("verdict rule frozen before the effect vectors were
  assembled") is false given post-outcome amendment 2 — main text must carry the
  disclosure. Verdict vocabulary (M*, bands, G2, unresolvable) undefined for readers;
  registrations must be reproduced verbatim in the supplement and deposited with a DOI at
  submission. Multivariate-only ignores per-trait evidence on the ten excluded traits:
  check the 336 against Bonferroni-significant per-trait Q among all 24 traits and
  re-derive the headline. Pilot loci and selection rule never listed; H_V4
  "not evaluated" needs a dated deviation rationale; H_GW2 failure needs one main-text
  sentence. Marginal-only (89 loci, 15% of union) undiscussed. Exploratory labels missing
  in main text for the 16.3% and 37.2% figures.

## Post-round verification (session)

Confirmed against genomewide_vector_q_v2_results.json: chr1_159175494_C_T (811 bp from
rs2814778 at GRCh37 1:159174683) is rotation_only with zero significant scalars;
chr1_159892088_G_A (~717 kb away, a separate clumped locus) is concordant_significant
with 5 scalars. Table 3's rs2814778 numbers are a genuine direct query of that variant
(WBC scalar Q reproduces to all printed digits against the current release).

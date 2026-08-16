# Provenance notebook

Dated entries for corrections and findings that the papers themselves do not narrate.
Each entry names the record that settled it.

## 2026-08-16 — corrections found while preparing paper_v7 / supplementary_v5

- **Diagonal-statistic claim.** paper_v6 §3.3 stated the diagonal Q_V is "also significant
  at 95% of loci." Recomputed from `analysis/results/genomewide_vector_q_v2_results.json`
  (`sum_scalar_Q` vs chi-square at each locus's stored df): 20.7% at nominal alpha on the
  2,567-locus set; 10.0% under BH FDR 0.05 on the 2,303-locus evaluation set. The
  bootstrap additionally shows the diagonal reference is anti-conservative under
  correlated traits (rejection 0.139 at nominal 0.05). v7 reports the correct numbers.
- **Null-variant count.** paper_v6 said the correlation matrix used 8,412 chr22 EUR
  variants; the stored estimate (`analysis/data/phenotypic_correlations_estimated.json`)
  records 9,117 shared variants. v7/supplement report 9,117.
- **Chromosome-1 stability claim.** paper_v6 supplement claimed mean absolute element-wise
  difference < 0.02 between chr22 and chr1 correlation estimates. Recomputed with the
  chr1 estimate produced under the frozen robustness registration
  (`sensitivity-artifact-checks/results/corr_chr1_eur.json`): 0.094. The functional
  sensitivity (full-panel multivariate-only share shifts 0.2 points) is what v7 reports.
- **Supplement S2 prediction list.** The old supplement described "six pre-registered
  directional predictions" (ACKR1, HBB region, SH2B3, CEBPA region, JAK2 region, ABO)
  with all six confirmed. No registration containing these per-locus predictions exists
  in the repository or its git history; the six loci trace to the hand-picked
  `TARGET_LOCI` dictionary in `analysis/scripts/download_panukb_blood_traits.py`. The
  actual frozen registration of 2026-07-17 (`analysis/PREREGISTRATION.md`) contains
  H_V1–H_V6, whose outcomes (5 held, 1 not evaluated) are now reported in supplementary
  Note S2, together with the script-registered H_GW1–H_GW6.
- **Unreported registered failure.** H_GW2 (all seven pilot loci recovered by the
  genome-wide scan) failed: no locus within 500 kb of rs9349379 (PHACTR1) entered the
  discovery set, in the July data or the current release. First reported in
  supplementary_v5.
- **Trait-exclusion story.** The old supplement's Table S1 said immature reticulocyte
  fraction was the excluded trait and marked basophil percentage as used. The stored
  correlation metadata shows the 24-trait panel excludes basophil percentage (fewer than
  100 usable null variants in the chr22 window); immature reticulocyte fraction was not
  part of the 25-trait extraction set. Corrected in supplementary_v5 Table S1.
- **Data release.** Pan-UKB re-released flat files between 2026-07 and 2026-08-16; 262 of
  2,567 loci are absent from the current release (uniform variant mask), values at
  surviving variants unchanged. Recorded in the robustness registration chain
  (amendment 1, commit 10fbfba26adc); the papers report the current release only.

Robustness suite outcomes and the SURVIVES WITH CAVEAT verdict are recorded in
`sensitivity-artifact-checks/` (registration chain: 4955185ed63b, 10fbfba26adc,
c4eb0b684a38) and its results directory.

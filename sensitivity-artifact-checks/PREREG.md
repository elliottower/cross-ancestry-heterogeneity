# Pre-registration: Robustness of the multivariate-only heterogeneity claim

**Status:** FROZEN at `4955185ed63b`
**Plan sha256:** `4594b705b27fd880c363a28b5443f67be43bab6e9f496015912153048ef52f94`
**Frozen:** 2026-08-16

Extension of the 2026-07-17 registration (`analysis/PREREGISTRATION.md`). That
registration governs discovery (H_V1–H_V6, H_GW1–H_GW6); its outcomes exist on disk and
will be reported separately. Nothing here modifies it. This registration governs one
question: does the v2 headline — 88.5% of significant loci multivariate-only — survive
correction for the artifact channels identified after submission?

## Foreknowledge

Stated as fact; the predictions below are informed, not blind, and this registration binds
interpretation rather than discovery. Known before drafting: the full v2 results (95.9% of
2,567 loci significant; 88.5% multivariate-only; Q_V/ΣQ ratio median 1.79, IQR 1.54–2.08);
the diagonal statistic recomputed 2026-08-15 from the stored v2 file is significant at
20.7% of loci overall and 15.3% of headline loci, against the manuscript's erroneous claim
of 95%; a v1 run under identity correlation found 22.2% significant; a separate session
reported condition number ≈3,553, smallest eigenvalue ≈0.00116, effective dimensionality
11–16 (unverified; will be recomputed); the panel's arithmetic redundancies
(MCV = Hct/RBC, MCH = Hb/RBC, MCHC = Hb/Hct, plateletcrit ≈ PLT×MPV, differential
percentages = counts/WBC) are known from field definitions; the power simulations and
ACKR1 results are known. No fine-mapping resource has been opened: the BCX2 and Vuckovic
credible-set files named in H_R6 are known only from their documentation, and H_R6's gate
covers the possibility that their coverage of the 2,567 loci is insufficient.

## Definitions

- **Significant**: Benjamini–Hochberg FDR 0.05 across all 2,567 loci, per test variant.
- **Multivariate-only share (M)**: among significant loci, the fraction with no per-trait
  Q significant at Bonferroni α/d, d the variant's trait count. Baseline M = 88.5%.
- **Non-redundant panel**: from the paper's 24 traits, drop hematocrit, MCV, MCH, MCHC,
  plateletcrit, the five differential percentages, and reticulocyte percentage; keep the
  rest.
- **Truncation rank k**: smallest number of leading eigenvalues of the EUR correlation
  matrix whose sum reaches 95% of the trace; truncated Q_V uses df = k(K−1).
- **Discovery traits**: the traits genome-wide significant (p < 5×10⁻⁸) in EUR at that
  locus — the ascertainment criterion.
- **M\***: M under H_R1's non-redundant panel — the corrected headline rate. Two
  registered bands: **minimal correction** (M\* ≥ 75%) and **substantial correction**
  (50% ≤ M\* < 75%).

## Predictions

**H_R1 — conditioning (carries the design).** On the non-redundant panel with full
estimated covariance, M > 50%.
*Rationale*: if the dominant-pattern claim needs the arithmetically redundant traits, it is
a property of the matrix, not the biology.

**H_R2 — truncation consistency.** Under rank-k truncated Q_V on the full panel, M > 50%.
*Rationale*: the same claim through an independent conditioning treatment; agreement of
H_R1 and H_R2 rules out that either fix manufactured the result.

**H_R3 — covariance-estimate stability.** Re-estimating the correlation matrix from
chromosome-1 EUR null variants and from AFR null variants each changes full-panel M by
less than 10 percentage points.
*Rationale*: sensitivity of M to sampling error in the plug-in matrix, the channel a
same-matrix simulation cannot detect.

**H_R4 — ascertainment, non-EUR-only.** With EUR excluded (K = 5), M > 50% among
significant loci. Gated: evaluated only if ≥ 50 loci are significant in this
configuration; otherwise reported as unresolvable, not failed.
*Rationale*: loci were selected on EUR significance; excluding EUR makes the tested data
disjoint from selection, removing winner's curse entirely at the cost of power.

**H_R5 — ascertainment, discovery traits excluded.** Dropping each locus's discovery
traits from its vector, M > 50% among significant loci.
*Rationale*: winner's-curse inflation lives only on the ascertained traits; if the signal
persists on the unascertained remainder, selection bias does not explain it.

**H_R6 — LD tagging, credible-set subset.** Among loci whose index variant is a member of
a 95% credible set of size ≤ 20 for at least one discovery trait in the BCX2 trans-ancestry
fine-mapping (fallback: Vuckovic 2020 UK Biobank credible sets; exact GRCh37
chr:pos:ref:alt match), M > 50%, computed with H_R1's non-redundant-panel statistic
restricted to the qualifying loci. Gated: evaluated only if ≥ 100 loci qualify; otherwise
reported as unresolvable and the LD channel is deferred to a later registration.
*Rationale*: differential tagging predicts heterogeneity concentrated where the index
variant merely tags the causal variant; if M persists where the index variant is likely
causal, tagging does not explain it.

**H_R7 — magnitude (strengthening outcome).** Among loci significant under H_R1's panel,
the median maximum pairwise rotation angle θ_max exceeds the 95th percentile of the
bootstrap-null distribution of that median.
*Rationale*: converts significance into magnitude; a claim that survives H_R1–H_R6 and
also exceeds null-expected rotation is stronger than the submitted version.

## Gates

**G1 — fidelity. All predictions are void if recomputed full-panel Q_V fails to match the
stored v2 Q_V within 1% at ≥ 99% of loci** (guards against a Pan-UKB data-version
mismatch in the tabix recovery).

**G2 — calibration, evaluated separately per Q_V variant. A variant is void if its
parametric-bootstrap rejection under the homogeneous null falls outside [0.03, 0.07] at
nominal 0.05.** A hypothesis whose variant is void is unresolvable; if H_R1's variant is
void, H_R2 becomes the design-carrying hypothesis and M\* is computed from it.

## Overall verdict

What "passing robustness" means, committed now:

| verdict | condition | consequence for the paper |
|---|---|---|
| SURVIVES, STRENGTHENED | H_R1, H_R2, H_R5 hold; H_R3 holds; H_R4 and H_R6 each hold or unresolvable; H_R7 holds | resubmit with corrected rate and rotation distribution as the headline |
| SURVIVES | as above but H_R7 fails | resubmit with corrected rate; magnitude reported descriptively |
| SURVIVES WITH CAVEAT | H_R1 and H_R2 hold; exactly one of H_R3/H_R5/H_R6 fails, or H_R4 fails outright | corrected rate as headline; the failed channel stated in Limitations; methods venue |
| DOES NOT SURVIVE | H_R1 fails, or H_R2 fails, or two or more of H_R3/H_R4/H_R5/H_R6 fail | the paper becomes the corrected-estimate paper: uncorrected 88.5% reported as the naive figure, the corrected value as the finding |

Every surviving verdict additionally carries the M\* band. Minimal correction (M\* ≥ 75%):
the corrected rate replaces 88.5% and the artifact analyses are reported as sensitivity
checks. Substantial correction (50% ≤ M\* < 75%): the artifact decomposition becomes a
co-result — the naive and corrected rates are reported together wherever the headline
appears.

Both terminal rows produce a submittable paper; no outcome is shelved.

## Analysis plan

1. Recover β̂ and SE per ancestry per trait at the 2,567 v2 loci by remote tabix against
   Pan-UKB; one resumable append-only shard per trait file; skip existing shards.
2. Verify G1. Recompute condition number, eigenvalues, and k.
3. Compute Q_V variants: full plug-in, non-redundant panel, rank-k truncated, diagonal,
   non-EUR-only, discovery-traits-excluded; per-trait Q throughout; BH per variant.
4. Parametric bootstrap: β̂_i ~ N(β̂_w, Σ_i) per locus for G2 and the H_R6 null;
   misspecification stress for H_R3's matrices.
5. Download BCX2 95% credible sets (Lettre Lab; fallback Vuckovic 2020); annotate each
   locus's index variant by exact GRCh37 chr:pos:ref:alt match; evaluate H_R6's gate.
6. θ_max and pairwise rotation angles on the v2 configuration, reusing the pilot code in
   `run_vector_analysis.py`.
7. Evaluate H_R1–H_R7, apply the verdict table and M\* band, correct the manuscript's
   diagonal sentence with the recomputed value.

## Exploratory (labeled so in any writeup)

Eigendirection decomposition of per-locus Q_V; conditional winner's-curse debiasing of
EUR betas; pairwise rotation matrices; reporting of the 2026-07-17 registration's
outcomes, which discharges that registration and is not part of this one.

---

## Log

Append only. Never edit above the line.

The last column is what distinguishes an amendment from a deviation, so you do not have to
decide which word to use: `nothing run`, `no results seen`, `results not opened`, `results seen`.

```
2026-08-15  created                              nothing run
2026-08-15  rewritten in house format (H_R*, verdict table) before freeze   nothing run
2026-08-15  external review: added M* bands and H_R6 LD-tagging check; strengthening outcome renumbered H_R7   nothing run
2026-08-16  frozen at 4955185ed63b                nothing run
2026-08-16  erratum: analysis plan step 4 'H_R6 null' should read 'H_R7 null' (renumbering typo found post-freeze)  nothing run
2026-08-16  reviewer note: credible-set membership in H_R6 is a soft proxy for index-variant causality; residual tagging ambiguity to be stated in Limitations  nothing run
```

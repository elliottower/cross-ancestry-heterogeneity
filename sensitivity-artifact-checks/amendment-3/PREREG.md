# Amendment 3: revision analyses on the primary panel

**Status:** FROZEN at `67bf3060f5fd`
**Plan sha256:** `709a5435bae0eeb1d55673c4647f4b897e119f8dde77c4a307cfedbb1a5f3694`
**Frozen:** 2026-08-16

Amendment to the chain 4955185ed63b → 10fbfba26adc → c4eb0b684a38. A three-referee
review of the v8 manuscript (record: `sensitivity-artifact-checks/reviews/simulated_round1.md`)
found that the ascertainment controls (H_R4, H_R5) and the covariance-stability gate
(H_R3) were evaluated on the full 24-trait panel while the paper's primary estimate is
the 14-trait panel, and that the multivariate-only definition ignores per-trait evidence
on the ten excluded traits. This amendment re-scores the verdict on the primary panel and
registers the analyses that have not yet been run.

## Foreknowledge

Stated as fact. All prior outcomes are known. Two primary-panel quantities relevant to
this amendment are already computed and are therefore foreknowledge, not predictions:
under the AFR-estimated correlation matrix the primary-panel share is 37.2% (shift 29.2
points — the re-scored stability check **fails** the registered 10-point bound), and
under the chromosome-1 EUR matrix it is 64.9% (shift 1.5 points, within bound). The
following have not been computed: the non-EUR and discovery-traits-excluded analyses on
the 14-trait panel; the all-24-trait per-trait check of the 336 multivariate-only loci;
rotation angles on 14-trait vectors; the moment-matched recalibration of the diagonal
statistic; the eigenvalue spectrum of the 14-trait correlation matrix.

## Registered analyses and criteria

- **H_R3′ — covariance stability, primary panel.** Both alternative-matrix arms shift the
  14-trait share by < 10 points. **Outcome already known: fails** (AFR arm 29.2 points).
  Recorded here so the verdict is scored on the quantity the paper certifies.
- **H_R4′ — non-EUR only, primary panel.** 14-trait Q_V with EUR excluded; M > 50% among
  BH-significant loci. Gated: evaluated only if ≥ 50 loci are significant; otherwise
  unresolvable.
- **H_R5′ — discovery traits excluded, primary panel.** Drop each locus's discovery traits
  from the 14-trait vector (d′ ≥ 2 required); M > 50% among significant loci.
- **H_R7′ — rotation, primary panel.** θ_max computed on 14-trait vectors among
  H_R1-significant loci; observed median exceeds the 95th percentile of the
  homogeneity-bootstrap null medians.
- **H_R8 — practitioner-relevant multivariate-only share.** M₂₄ = the share of
  H_R1-significant loci with no Bonferroni-significant per-trait Q among all 24 measured
  traits (α/24, per-locus). The per-trait procedure is identical to the one already used
  in the 24-trait classification: scalar Cochran's Q against χ²_{K−1} on the locus's
  recorded ancestry set, Bonferroni-corrected within locus. Criterion: M₂₄ > 50%. M₂₄
  replaces the panel-internal share as the headline number in the manuscript regardless
  of outcome.

Descriptive, no criteria: eigenvalue spectrum and condition number of the 14-trait
matrix; moment-matched diagonal recalibration (replacing the miscalibrated 231-locus
count wherever quoted). Exploratory, labeled so: a 13-trait panel without white blood
cell count (the surviving WBC ≈ sum-of-differential near-identity).

## Verdict re-scoring

The frozen verdict table is re-evaluated with H_R3′, H_R4′, H_R5′, H_R7′ in place of
their originals, H_R6 unchanged (already computed with the 14-trait statistic), and the
dominance requirement read as H_R1 **and** H_R8. The channels are H_R3′, H_R4′, H_R5′,
and H_R6, and only channels count toward failure tallies. H_R7′ is not a channel: as in
the parent table, its outcome separates SURVIVES, STRENGTHENED from SURVIVES and cannot
contribute to DOES NOT SURVIVE; with H_R3′ failed the strengthened row is unreachable, so
H_R7′ is reported but cannot change the verdict. Known consequence, disclosed now: with
H_R3′ failed, the best available verdict is SURVIVES WITH CAVEAT (the caveat row's
"exactly one channel fails"), reachable only if H_R5′, H_R6, and H_R8 hold and H_R4′
holds or is unresolvable. If H_R5′ or H_R8 fails, or a second channel fails, the verdict
is DOES NOT SURVIVE and the manuscript moves to the pre-committed co-finding framing: the
conditioning artifact and the diagnostic as the durable contribution, with the corrected
share reported as a specification-dependent range (16.3%–66.4% on current numbers).

**Disclosure.** This amendment was drafted with all prior outcomes known, including that
H_R3′ fails. The open quantities are H_R4′, H_R5′, H_R7′, and H_R8. No sham prediction is
made about quantities already observed.

Manuscript corrections accompanying this amendment but outside it (errors, not analyses):
the Table 3 degrees of freedom (120, not 125), base-panel labels on every Table 2 row,
the Methods registration sentence (must disclose the post-outcome amendment 2 in the main
text), the ACKR1 reconciliation (the pipeline's region locus chr1_159175494 is
multivariate-only while hand-queried rs2814778 is concordant-heterogeneous — to be
presented explicitly as index-choice sensitivity), and verdict-vocabulary definitions in
the supplement.

---

## Log

Append only. Never edit above the line.

The last column is what distinguishes an amendment from a deviation, so you do not have to
decide which word to use: `nothing run`, `no results seen`, `results not opened`, `results seen`.

```
2026-08-16  created                              results seen
2026-08-16  frozen at 67bf3060f5fd                nothing run
2026-08-16  Tier-1 outcomes: H_R4' holds (221 sig >= gate, M=70.6%, calibration 0.051); H_R5' holds (147 sig, M=67.4%); H_R7' holds (117.0 vs null95 111.0); H_R8 holds (M24 = 333/506 = 65.8%, becomes headline). Channels: H_R3' fails (known), H_R4'/H_R5'/H_R6 hold -> exactly one channel fails -> verdict SURVIVES WITH CAVEAT re-earned on primary panel. Descriptive: 14-trait spectrum cond=189, lambda_min=0.017; diagonal recalibrated to 53 sig (2.3%). Exploratory 13-trait (no WBC): M=65.4%, cond=22.6.  results seen
```

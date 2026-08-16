# Amendment 1: evaluation set under the current Pan-UKB release

**Status:** FROZEN at `10fbfba26adc`
**Plan sha256:** `dae8e549a02acd54ab730c2623c7774ad323f9f7c9d69f644b92b97c6980bdf0`
**Frozen:** 2026-08-16

Amendment to the registration frozen at `4955185ed63b` (sha256 `4594b705b27fd880…`,
`sensitivity-artifact-checks/PREREG.md`). The parent's G1 gate failed on coverage and its
outcome is logged there; this amendment redefines the evaluation set and the fidelity gate
for the current data release. Everything not restated here — H_R1 through H_R7, gate G2,
the verdict table, the M\* bands, the analysis plan — carries over from the parent
verbatim, applied to the amended evaluation set.

## What changed and why

Pan-UK Biobank has re-released its summary-statistic flat files since the July analysis.
The current release lacks 262 of the 2,567 registered loci, uniformly across trait files
(a revised variant mask, not per-file gaps); two further loci lack one of their recorded
ancestry entries. Values at surviving variants are unchanged: 99.87% of the 2,303
rebuildable loci reproduce the stored v2 Q_V within 1%, and a spot check reproduces a
stored scalar Q to all printed digits. Three loci differ by more than 1% (worst relative
difference 0.47); they remain in the evaluation set because the registered fidelity bar
is a set-level criterion and it is met.

The paper will state the release it uses and report numbers from that release. The
version history lives in this registration chain and the repository log, not in the
paper.

## Amended definitions

- **Evaluation set**: the 2,303 loci present in the current release with all 24 panel
  traits and all v2-recorded ancestry entries. All H_R hypotheses, both G2 arms, the
  verdict table, and the M\* bands are evaluated on this set.
- **G1a — coverage (descriptive)**: report the count and identity of registered loci
  absent from the current release. No threshold; recorded as a data note.
- **G1b — value fidelity (gate)**: recomputed full-panel Q_V must match the stored v2
  Q_V within 1% at ≥ 99% of evaluation-set loci. Already evaluated: 99.87%, PASS.

## Foreknowledge added since the parent freeze

Stated as fact. The G1 numbers above are known. The intersection baseline is known:
restricted to the 2,303 evaluation loci, the stored v2 statistics give 2,203 significant
under Benjamini–Hochberg FDR 0.05, of which 1,951 are multivariate-only — M = 88.56%,
against 88.5% on the full registered set. The eigen recompute confirms condition number
3,553, smallest eigenvalue 0.00116, truncation rank k = 12 under the registered 95%-trace
rule. Both H_R3 correlation matrices are estimated and complete on the full panel (chr1
EUR: 8,639 shared null variants; chr22 AFR: 15,967); their contents have not been
compared to the chr22 EUR matrix. BCX2 trans-ancestry credible sets are downloaded; no
locus has been matched against them. No H_R hypothesis outcome has been computed.

## Unchanged commitments, restated for clarity

The 50% dominance thresholds, the 10-point stability bound, the H_R4 and H_R6 gates, the
G2 calibration band [0.03, 0.07] per variant, the H_R7 rotation criterion, the verdict
table, and the M\* bands (minimal ≥ 75%; substantial 50–75%) apply exactly as frozen in
the parent, with M\* computed on the evaluation set.

---

## Log

Append only. Never edit above the line.

The last column is what distinguishes an amendment from a deviation, so you do not have to
decide which word to use: `nothing run`, `no results seen`, `results not opened`, `results seen`.

```
2026-08-16  created                              results seen
2026-08-16  frozen at 10fbfba26adc                nothing run
```

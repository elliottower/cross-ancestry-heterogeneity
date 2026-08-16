# Amendment 4: resolving the resolvable limitations

**Status:** FROZEN at `9fbed73f34cf`
**Plan sha256:** `b7cc3db59a7e90719fa28a7488e2384ddec9386b59a6558b8dcfcda6262e2989`
**Frozen:** 2026-08-16

Amendment to the chain 4955185ed63b → 10fbfba26adc → c4eb0b684a38 → 67bf3060f5fd. Three
of the manuscript's stated or planned limitations are resolvable with data in hand: the
dependence of the estimate on a single EUR-derived correlation matrix, the circularity of
the parametric-bootstrap calibration, and the possibility that per-ancestry trait
standardization manufactures the distributed-heterogeneity signature. This amendment
registers one check for each. **It is the final registration before submission; after
these checks, the manuscript is written and submitted regardless of outcome, and only
error corrections follow.**

## Foreknowledge

All prior outcomes are known, including: primary-panel M = 66.4%, M₂₄ = 65.8%; the
full-matrix AFR substitution gives 37.2%; per-ancestry null-variant counts from the H_R3
estimation (AFR: 15,967 shared; chr1 EUR: 8,639). Nothing has been computed for:
per-ancestry-weighted Q_V, empirical calibration at null variants, or scale-absorption.

## Registered checks and criteria

- **H_R9 — empirical calibration at null variants.** Sample one variant per 1-Mb genomic
  window from the trait files, excluding any variant with EUR $-\log_{10}p \geq 4$ for
  any of the 24 traits (target $n \geq 2{,}000$ testable variants with complete data,
  same inclusion rules as the locus analysis). Compute the 14-trait and 24-trait $Q_V$
  at these variants against their $\chi^2$ references. Criterion: 14-trait rejection
  rate at nominal 0.05 lies in $[0.03, 0.07]$. Also reported, no criterion: rejection at
  0.01 and 0.001, and the 24-trait rate. **Gate: if H_R9 fails, all $\chi^2$-based
  prevalence estimates are demoted to descriptive and the manuscript moves to the
  co-finding framing.** H_R9 is computed under the single-EUR-matrix specification and
  stands as the general diagnostic regardless of the other checks; for the H_R10
  primary specification, H_R10's own calibration criterion supersedes it.
- **H_R10 — per-ancestry covariance as the primary specification.** Estimate $R_i$ for
  each ancestry from chromosome-22 null variants by the standing protocol. Shrinkage
  rule, pinned: $\tilde R_i = \lambda_i R_i + (1-\lambda_i) R_{\mathrm{EUR}}$ with
  $\lambda_i = n_i/(n_i + 5{,}000)$, where $n_i$ is that ancestry's shared-null-variant
  count; any ancestry with $n_i < 1{,}000$ uses $R_{\mathrm{EUR}}$. Recompute the
  14-trait analysis with $\Sigma_i = D_i \tilde R_i D_i$, including its bootstrap
  calibration and BH significance. Criteria: the analysis is calibrated
  ($[0.03, 0.07]$), and its $M_{24}$ analog exceeds 50\%. **Commitment: this
  specification becomes the manuscript's primary analysis regardless of outcome; the
  single-EUR-matrix analysis becomes the sensitivity.** If H_R10 holds, the
  correlation-structure caveat narrows to the shrinkage estimation of small-group
  matrices; if it fails, the verdict drops to DOES NOT SURVIVE and the co-finding
  framing applies with the share reported as a range.
- **H_R11 — scale-absorption.** For each non-EUR ancestry $i$ and trait $j$, fit one
  scale factor $s_{ij}$ by inverse-variance-weighted regression of $\beta_{ij\ell}$ on
  the EUR estimate $\beta_{\mathrm{EUR},j\ell}$ across all 2{,}303 loci. Rescale
  ($\beta_{ij\ell}/s_{ij}$, $\mathrm{se}_{ij\ell}/s_{ij}$) and recompute the H_R10
  primary analysis. Criteria: rescaled $M_{24}$ analog exceeds 50\%. If it holds, the
  standardization concern is a resolved check; if it fails, per-ancestry scaling
  explains the pattern and is reported as the finding, under the co-finding framing.
  Scope note: with EUR as the regression reference, this check detects
  non-EUR-relative-to-EUR scale mismatch only; a scale artifact shared by all groups,
  including EUR, is outside its reach and remains covered by the standing limitation.

Descriptive, no criteria: the fraction of aggregate $Q_V$ mass at significant loci
absorbed by the fitted rescaling; the fitted $s_{ij}$ themselves (deviations from 1
indicate scale mismatch or ascertainment effects, and the discovery-trait entries are
expected to be inflated by winner's curse — noted so their deviation is not read as
scale mismatch).

## Manuscript consequence

If H_R9, H_R10, and H_R11 all hold: verdict remains SURVIVES WITH CAVEAT with a narrowed
caveat; headline numbers come from the H_R10 primary specification; the three former
limitations are reported as checks in Results. If any fails: DOES NOT SURVIVE; the
conditioning artifact, the calibration diagnostic, and the specification-dependence
finding become the paper, with the multivariate-only share reported as a range across
specifications. Both branches are written and submitted without further registration.

**Disclosure.** Drafted with all prior outcomes known; the three registered quantities
are unobserved. The shrinkage constant (5,000) and the null-variant sampling rule were
fixed at drafting, before any estimation.

---

## Log

Append only. Never edit above the line.

The last column is what distinguishes an amendment from a deviation, so you do not have to
decide which word to use: `nothing run`, `no results seen`, `results not opened`, `results seen`.

```
2026-08-16  created                              results seen
2026-08-16  frozen at 9fbed73f34cf                nothing run
```

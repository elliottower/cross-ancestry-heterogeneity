"""H_R10 and H_R11 (amendment 4, frozen 9fbed73f34cf).

H_R10: 14-trait primary analysis with per-ancestry shrunk correlation matrices
(lambda_i = n_i/(n_i + 5000); n_i < 1000 falls back to R_EUR), with bootstrap
calibration under the same per-ancestry model. H_R11: inverse-variance-weighted
per-ancestry-trait scale factors fitted against EUR, data rescaled, H_R10
analysis recomputed. Criteria per the frozen amendment.

    uv run python sensitivity-artifact-checks/scripts/per_ancestry_primary.py
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats

from robustness_analysis import (ALPHA, DROP_CODES, GWS_NEGLOG10, bh_significant,
                                 build_matrices, load_shards, log, scalar_q)

SHRINK_N0 = 5_000
MIN_SHARED = 1_000
G2_BAND = (0.03, 0.07)
ANCS = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]


def load_matrices(codes):
    base = json.load(open("analysis/data/phenotypic_correlations_estimated.json"))
    assert base["trait_codes"] == codes
    R_eur = np.array(base["correlation_matrix"])
    paths = {"AFR": "corr_chr22_afr.json", "CSA": "corr_chr22_csa.json",
             "EAS": "corr_chr22_eas.json", "MID": "corr_chr22_mid.json",
             "AMR": "corr_chr22_amr.json"}
    R = {"EUR": R_eur}
    lam = {"EUR": 1.0}
    for anc, fn in paths.items():
        d = json.load(open(Path("sensitivity-artifact-checks/data") / fn))
        assert d["trait_codes"] == codes and d["complete"]
        n = d["n_shared_variants"]
        if n < MIN_SHARED:
            R[anc], lam[anc] = R_eur, 0.0
        else:
            l = n / (n + SHRINK_N0)
            R[anc] = l * np.array(d["correlation_matrix"]) + (1 - l) * R_eur
            lam[anc] = round(l, 3)
    return R, lam


def analyze(loci, codes, nonred, R_nr_inv, rescale=None):
    """Per-ancestry-covariance 14-trait Q_V. rescale: {anc: s-vector over 24}."""
    entries = []
    for lid, (betas, ses, ancs, disc, sp24) in loci.items():
        b = betas.copy()
        s = ses.copy()
        if rescale:
            for k, a in enumerate(ancs):
                if a in rescale:
                    b[k] = b[k] / rescale[a]
                    s[k] = s[k] / rescale[a]
        b = b[:, nonred]
        s = s[:, nonred]
        K = len(ancs)
        W = [R_nr_inv[a] * np.outer(1 / s[k], 1 / s[k]) for k, a in enumerate(ancs)]
        Wsi = np.linalg.inv(np.sum(W, axis=0))
        bbar = Wsi @ np.sum([W[k] @ b[k] for k in range(K)], axis=0)
        Q = float(sum((b[k] - bbar) @ W[k] @ (b[k] - bbar) for k in range(K)))
        p = float(stats.chi2.sf(Q, len(nonred) * (K - 1)))
        entries.append({"locus_id": lid, "p": p, "Q": Q, "sp24": sp24})
    sig = bh_significant([e["p"] for e in entries])
    n_sig = int(sig.sum())
    n_mo24 = sum(1 for e, g in zip(entries, sig)
                 if g and not any(p < ALPHA / 24 for p in e["sp24"]))
    return entries, sig, {"n_fdr_significant": n_sig, "n_multivariate_only_24": n_mo24,
                          "M24": round(n_mo24 / n_sig, 4) if n_sig else None}


def calibrate(loci, nonred, R_nr, R_nr_inv, n_reps, rng):
    chol = {a: np.linalg.cholesky(R_nr[a]) for a in R_nr}
    n_rej = n = 0
    for lid, (betas, ses, ancs, _, _) in loci.items():
        s = ses[:, nonred]
        K = len(ancs)
        W = [R_nr_inv[a] * np.outer(1 / s[k], 1 / s[k]) for k, a in enumerate(ancs)]
        Wsi = np.linalg.inv(np.sum(W, axis=0))
        crit = stats.chi2.isf(ALPHA, len(nonred) * (K - 1))
        for _ in range(n_reps):
            e = np.stack([(rng.standard_normal(len(nonred)) @ chol[a].T) * s[k]
                          for k, a in enumerate(ancs)])
            bbar = Wsi @ np.sum([W[k] @ e[k] for k in range(K)], axis=0)
            Q = sum((e[k] - bbar) @ W[k] @ (e[k] - bbar) for k in range(K))
            n_rej += Q > crit
            n += 1
    rate = n_rej / n
    return {"n_draws": int(n), "reject_rate": round(float(rate), 5),
            "pass": bool(G2_BAND[0] <= rate <= G2_BAND[1])}


def main():
    rng = np.random.default_rng()
    corr = json.load(open("analysis/data/phenotypic_correlations_estimated.json"))
    codes = corr["trait_codes"]
    nonred = [j for j, c in enumerate(codes) if c not in DROP_CODES]
    R, lam = load_matrices(codes)
    R_nr = {a: R[a][np.ix_(nonred, nonred)] for a in R}
    R_nr_inv = {a: np.linalg.inv(R_nr[a]) for a in R}
    log(f"shrinkage lambdas: {lam}")

    v2 = {r["locus_id"]: r for r in json.load(
        open("analysis/results/genomewide_vector_q_v2_results.json"))["results"]}
    data = load_shards(Path("sensitivity-artifact-checks/data/shards"), codes)
    loci = {}
    for lid, rec in v2.items():
        built = build_matrices(data.get(lid, {}), codes, rec["ancestries"])
        if built is None:
            continue
        betas, ses = built
        ancs = rec["ancestries"]
        disc = [j for j, c in enumerate(codes)
                if (e := data[lid][c].get("EUR")) and e["neglog10_pval"] is not None
                and e["neglog10_pval"] >= GWS_NEGLOG10]
        sp24 = [float(stats.chi2.sf(scalar_q(betas[:, j], ses[:, j]), len(ancs) - 1))
                for j in range(len(codes))]
        loci[lid] = (betas, ses, ancs, disc, sp24)
    log(f"{len(loci)} loci")

    # --- H_R10
    cal = calibrate(dict(list(loci.items())[:600]), nonred, R_nr, R_nr_inv, 40, rng)
    entries, sig, summ = analyze(loci, codes, nonred, R_nr_inv)
    h10_ok = cal["pass"] and summ["M24"] is not None and summ["M24"] > 0.5
    h10 = {**summ, "calibration": cal, "shrinkage": lam,
           "outcome": ("holds" if h10_ok else
                       "fails" if cal["pass"] else "fails (calibration)")}
    log(f"H_R10: sig={summ['n_fdr_significant']} M24={summ['M24']} "
        f"cal={cal['reject_rate']} -> {h10['outcome']}")

    # --- H_R11: fit s_ij against EUR across all loci
    s_fit = {}
    for a in ANCS:
        if a == "EUR":
            continue
        s_vec = np.ones(len(codes))
        for j in range(len(codes)):
            num = den = 0.0
            for lid, (betas, ses, ancs, _, _) in loci.items():
                if a not in ancs or "EUR" not in ancs:
                    continue
                k, ke = ancs.index(a), ancs.index("EUR")
                w = 1.0 / ses[k, j] ** 2
                num += w * betas[ke, j] * betas[k, j]
                den += w * betas[ke, j] ** 2
            s_vec[j] = num / den if den > 0 else 1.0
        s_fit[a] = s_vec
    s_summary = {a: {"mean": round(float(np.mean(v)), 3),
                     "min": round(float(np.min(v)), 3),
                     "max": round(float(np.max(v)), 3)} for a, v in s_fit.items()}
    log(f"fitted scale factors: {s_summary}")

    entries_r, sig_r, summ_r = analyze(loci, codes, nonred, R_nr_inv, rescale=s_fit)
    h11_ok = summ_r["M24"] is not None and summ_r["M24"] > 0.5
    q_before = sum(e["Q"] for e, g in zip(entries, sig) if g)
    sig_ids = {e["locus_id"] for e, g in zip(entries, sig) if g}
    q_after = sum(e["Q"] for e in entries_r if e["locus_id"] in sig_ids)
    absorbed = 1 - q_after / q_before if q_before else float("nan")
    h11 = {**summ_r, "outcome": "holds" if h11_ok else "fails",
           "absorbed_fraction_at_h10_sig": round(float(absorbed), 4),
           "scale_factors": s_summary}
    log(f"H_R11: sig={summ_r['n_fdr_significant']} M24={summ_r['M24']} "
        f"absorbed={absorbed:.3f} -> {h11['outcome']}")

    out = {"metadata": {"timestamp": datetime.now(timezone.utc).isoformat(),
                        "amendment": "9fbed73f34cf"},
           "H_R10": h10, "H_R11": h11}
    with open("sensitivity-artifact-checks/results/per_ancestry_primary.json", "w") as f:
        json.dump(out, f, indent=1)
    log("saved sensitivity-artifact-checks/results/per_ancestry_primary.json")


if __name__ == "__main__":
    main()

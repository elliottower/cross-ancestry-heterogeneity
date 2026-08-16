"""Parametric bootstrap: G2 calibration per Q_V variant and the H_R7 rotation null.

Step 4 of the frozen registration. Under exact homogeneity the common effect
cancels from every weighted-deviation quadratic form, so calibration draws are
mean-zero: eps_k ~ N(0, Sigma_k), Sigma_k = D_k rho D_k. The H_R7 null keeps the
pooled mean per locus because rotation angles depend on the vector itself.

    uv run python sensitivity-artifact-checks/scripts/bootstrap_calibration.py
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats

from robustness_analysis import (
    DROP_CODES, GWS_NEGLOG10, build_matrices, load_shards, log,
)

ALPHA = 0.05
G2_BAND = (0.03, 0.07)


def weights_for(ses_k, r_inv):
    Dinv = 1.0 / ses_k
    return r_inv * np.outer(Dinv, Dinv)


def calibrate_variant(loci, r_gen, r_inv_test, df_per_stratum, drop_eur, keep_idx,
                      n_reps, rng):
    """Pooled bootstrap rejection rate at nominal ALPHA for one variant.

    loci: list of (betas, ses, ancestries, disc_idx).
    r_gen: generating correlation (on the variant's trait set).
    r_inv_test: testing inverse-correlation (weight matrix core).
    keep_idx: 'disc' to drop per-locus discovery traits, else a fixed index list.
    """
    L_gen = np.linalg.cholesky(r_gen) if not isinstance(keep_idx, str) else None
    n_reject = n_draws = 0
    for betas, ses, ancs, disc_idx in loci:
        if keep_idx == "disc":
            cols = [j for j in range(betas.shape[1]) if j not in disc_idx]
            if len(cols) < 2:
                continue
            r_sub = r_gen[np.ix_(cols, cols)]
            L = np.linalg.cholesky(r_sub)
            r_inv = np.linalg.inv(r_sub)
            d_df = len(cols)
        else:
            cols = keep_idx
            L = L_gen
            r_inv = r_inv_test
            d_df = df_per_stratum
        s = ses[:, cols]
        rows = [i for i, a in enumerate(ancs) if a != "EUR"] if drop_eur else list(range(len(ancs)))
        if len(rows) < 3:
            continue
        s = s[rows]
        K, d = s.shape
        df = d_df * (K - 1)
        W = [weights_for(s[k], r_inv) for k in range(K)]
        W_sum_inv = np.linalg.inv(np.sum(W, axis=0))
        crit = stats.chi2.isf(ALPHA, df)
        # draws: eps (n_reps, K, d) with eps_k = D_k L z
        z = rng.standard_normal((n_reps, K, d))
        eps = np.einsum("dj,rkj->rkd", L, z) * s[None, :, :]
        for r in range(n_reps):
            e = eps[r]
            bbar = W_sum_inv @ np.sum([W[k] @ e[k] for k in range(K)], axis=0)
            Q = sum((e[k] - bbar) @ W[k] @ (e[k] - bbar) for k in range(K))
            n_reject += Q > crit
            n_draws += 1
    rate = n_reject / n_draws if n_draws else float("nan")
    return {"n_draws": int(n_draws), "reject_rate": round(float(rate), 5),
            "pass": bool(G2_BAND[0] <= rate <= G2_BAND[1]) if n_draws else None}


def theta_max_of(betas):
    K = betas.shape[0]
    tmax = 0.0
    for i in range(K):
        for j in range(i + 1, K):
            ni, nj = np.linalg.norm(betas[i]), np.linalg.norm(betas[j])
            if ni < 1e-12 or nj < 1e-12:
                continue
            t = np.degrees(np.arccos(np.clip(betas[i] @ betas[j] / (ni * nj), -1, 1)))
            tmax = max(tmax, t)
    return tmax


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", type=Path, default=Path("sensitivity-artifact-checks/data/shards"))
    ap.add_argument("--corr", type=Path,
                    default=Path("analysis/data/phenotypic_correlations_estimated.json"))
    ap.add_argument("--v2", type=Path,
                    default=Path("analysis/results/genomewide_vector_q_v2_results.json"))
    ap.add_argument("--robustness", type=Path,
                    default=Path("sensitivity-artifact-checks/results/robustness_results.json"))
    ap.add_argument("--out", type=Path,
                    default=Path("sensitivity-artifact-checks/results/bootstrap_results.json"))
    ap.add_argument("--reps-calib", type=int, default=40)
    ap.add_argument("--reps-theta", type=int, default=500)
    args = ap.parse_args()

    rng = np.random.default_rng()

    with open(args.corr) as f:
        corr = json.load(f)
    trait_codes = corr["trait_codes"]
    rho = np.array(corr["correlation_matrix"])
    with open(args.v2) as f:
        v2 = {r["locus_id"]: r for r in json.load(f)["results"]}
    data = load_shards(args.shards, trait_codes)

    loci = []
    for lid, v2rec in v2.items():
        built = build_matrices(data.get(lid, {}), trait_codes, v2rec["ancestries"])
        if built is None:
            continue
        betas, ses = built
        disc_idx = []
        for j, code in enumerate(trait_codes):
            e = data[lid][code].get("EUR")
            if e and e["neglog10_pval"] is not None and e["neglog10_pval"] >= GWS_NEGLOG10:
                disc_idx.append(j)
        loci.append((lid, betas, ses, v2rec["ancestries"], disc_idx))
    log(f"{len(loci)} loci loaded for bootstrap")

    d = len(trait_codes)
    evals, evecs = np.linalg.eigh(rho)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    k_rank = int(np.searchsorted(np.cumsum(evals), 0.95 * d) + 1)
    rho_inv = np.linalg.inv(rho)
    rho_trunc_inv = evecs[:, :k_rank] @ np.diag(1.0 / evals[:k_rank]) @ evecs[:, :k_rank].T
    nonred_idx = [j for j, c in enumerate(trait_codes) if c not in DROP_CODES]
    rho_nr = rho[np.ix_(nonred_idx, nonred_idx)]

    plain = [(b, s, a, di) for _, b, s, a, di in loci]
    all_idx = list(range(d))
    g2 = {}
    for name, (r_gen, r_inv_test, dfs, drop_eur, keep) in {
        "full": (rho, rho_inv, d, False, all_idx),
        "nonredundant": (rho_nr, np.linalg.inv(rho_nr), len(nonred_idx), False, nonred_idx),
        "rank_truncated": (rho, rho_trunc_inv, k_rank, False, all_idx),
        "non_eur": (rho, rho_inv, d, True, all_idx),
        "no_discovery": (rho, None, None, False, "disc"),
        "diagonal": (rho, np.eye(d), d, False, all_idx),
    }.items():
        log(f"G2 calibrating {name}...")
        g2[name] = calibrate_variant(plain, r_gen, r_inv_test, dfs, drop_eur, keep,
                                     args.reps_calib, rng)
        log(f"  {name}: reject={g2[name]['reject_rate']} pass={g2[name]['pass']}")

    # --- H_R7 null: median theta_max over the H_R1-significant locus set
    with open(args.robustness) as f:
        rob = json.load(f)
    h_r1_sig = {r["locus_id"] for r in rob["variant_records"]["nonredundant"]
                if r["fdr_significant"]}
    subset = [(lid, b, s, a) for lid, b, s, a, _ in loci if lid in h_r1_sig]
    log(f"H_R7: {len(subset)} H_R1-significant loci")

    observed = float(np.median([theta_max_of(b) for _, b, s, a in subset]))
    L_full = np.linalg.cholesky(rho)
    null_medians = []
    pre = []
    for _, b, s, a in subset:
        W = [weights_for(s[k], rho_inv) for k in range(len(a))]
        W_sum_inv = np.linalg.inv(np.sum(W, axis=0))
        bbar = W_sum_inv @ np.sum([W[k] @ b[k] for k in range(len(a))], axis=0)
        pre.append((bbar, s))
    for rep in range(args.reps_theta):
        tms = []
        for bbar, s in pre:
            K = s.shape[0]
            z = rng.standard_normal((K, d))
            eps = (z @ L_full.T) * s
            tms.append(theta_max_of(bbar[None, :] + eps))
        null_medians.append(float(np.median(tms)))
    null_95 = float(np.percentile(null_medians, 95))
    h_r7 = {
        "observed_median_theta_max": round(observed, 2),
        "null_median_mean": round(float(np.mean(null_medians)), 2),
        "null_95th_percentile": round(null_95, 2),
        "n_loci": len(subset),
        "n_reps": args.reps_theta,
        "pass": bool(observed > null_95),
    }
    log(f"H_R7: observed={observed:.1f} null95={null_95:.1f} pass={h_r7['pass']}")

    out = {
        "metadata": {"timestamp": datetime.now(timezone.utc).isoformat(),
                     "reps_calib": args.reps_calib, "reps_theta": args.reps_theta,
                     "k_rank": k_rank},
        "G2": g2,
        "H_R7": h_r7,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    log(f"saved {args.out}")


if __name__ == "__main__":
    main()

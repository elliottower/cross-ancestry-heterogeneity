"""Amendment 3 (frozen 67bf3060f5fd): Tier-1 revision analyses on the primary panel.

Computes H_R4', H_R5', H_R7', H_R8, the descriptive items (14-trait spectrum,
moment-matched diagonal recalibration), the exploratory 13-trait panel, and
G2-style calibration for the two new variants.

    uv run python sensitivity-artifact-checks/scripts/tier1_revision.py
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats

from robustness_analysis import (
    ALPHA, DROP_CODES, GWS_NEGLOG10, bh_significant, build_matrices,
    load_shards, log, scalar_q,
)

G2_BAND = (0.03, 0.07)
H_R4_GATE = 50


def qv(betas, ses, r_inv):
    K, d = betas.shape
    W = [r_inv * np.outer(1 / ses[k], 1 / ses[k]) for k in range(K)]
    Wsum_inv = np.linalg.inv(np.sum(W, axis=0))
    bbar = Wsum_inv @ np.sum([W[k] @ betas[k] for k in range(K)], axis=0)
    Q = float(sum((betas[k] - bbar) @ W[k] @ (betas[k] - bbar) for k in range(K)))
    return Q, bbar


def classify_variant(entries):
    """entries: list of dicts with p, scalar_p, d. Returns summary with M."""
    sig = bh_significant([e["p"] for e in entries])
    n_sig = int(sig.sum())
    n_mo = sum(1 for e, s in zip(entries, sig)
               if s and not any(sp < ALPHA / e["d"] for sp in e["scalar_p"]))
    return {"n_loci": len(entries), "n_fdr_significant": n_sig,
            "n_multivariate_only": n_mo,
            "M": round(n_mo / n_sig, 4) if n_sig else None}, sig


def theta_max_of(b):
    K = b.shape[0]
    t = 0.0
    for i in range(K):
        for j in range(i + 1, K):
            ni, nj = np.linalg.norm(b[i]), np.linalg.norm(b[j])
            if ni < 1e-12 or nj < 1e-12:
                continue
            t = max(t, float(np.degrees(np.arccos(np.clip(b[i] @ b[j] / (ni * nj), -1, 1)))))
    return t


def calibrate(loci_rows, r_gen_chol, r_inv, d, n_reps, rng):
    n_rej = n = 0
    for ses_sub, K in loci_rows:
        W = [r_inv * np.outer(1 / ses_sub[k], 1 / ses_sub[k]) for k in range(K)]
        Wsum_inv = np.linalg.inv(np.sum(W, axis=0))
        crit = stats.chi2.isf(ALPHA, d * (K - 1))
        z = rng.standard_normal((n_reps, K, d))
        eps = np.einsum("dj,rkj->rkd", r_gen_chol, z) * ses_sub[None]
        for r in range(n_reps):
            e = eps[r]
            bbar = Wsum_inv @ np.sum([W[k] @ e[k] for k in range(K)], axis=0)
            Q = sum((e[k] - bbar) @ W[k] @ (e[k] - bbar) for k in range(K))
            n_rej += Q > crit
            n += 1
    rate = n_rej / n if n else float("nan")
    return {"n_draws": int(n), "reject_rate": round(float(rate), 5),
            "pass": bool(G2_BAND[0] <= rate <= G2_BAND[1]) if n else None}


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
                    default=Path("sensitivity-artifact-checks/results/tier1_revision.json"))
    ap.add_argument("--reps-theta", type=int, default=500)
    args = ap.parse_args()

    rng = np.random.default_rng()
    corr = json.load(open(args.corr))
    codes = corr["trait_codes"]
    rho = np.array(corr["correlation_matrix"])
    nonred_idx = [j for j, c in enumerate(codes) if c not in DROP_CODES]
    rho_nr = rho[np.ix_(nonred_idx, nonred_idx)]
    rho_nr_inv = np.linalg.inv(rho_nr)
    L_nr = np.linalg.cholesky(rho_nr)

    v2 = {r["locus_id"]: r for r in json.load(open(args.v2))["results"]}
    data = load_shards(args.shards, codes)
    rob = json.load(open(args.robustness))
    h_r1_sig = {r["locus_id"] for r in rob["variant_records"]["nonredundant"]
                if r["fdr_significant"]}
    full_scalar = {r["locus_id"]: r["scalar_p"] for r in rob["variant_records"]["full"]
                   } if "scalar_p" in rob["variant_records"]["full"][0] else None

    loci = {}
    for lid, rec in v2.items():
        built = build_matrices(data.get(lid, {}), codes, rec["ancestries"])
        if built:
            loci[lid] = (*built, rec["ancestries"])
    log(f"{len(loci)} loci loaded")

    out = {"metadata": {"timestamp": datetime.now(timezone.utc).isoformat(),
                        "amendment": "67bf3060f5fd"}}

    # --- descriptive: 14-trait spectrum
    ev = np.linalg.eigvalsh(rho_nr)
    out["spectrum_14"] = {"condition_number": float(ev.max() / ev.min()),
                          "lambda_min": float(ev.min()),
                          "n_below_0.05": int((ev < 0.05).sum()),
                          "eigenvalues": [round(float(x), 5) for x in sorted(ev)]}
    log(f"14-trait spectrum: cond={out['spectrum_14']['condition_number']:.1f} "
        f"lambda_min={out['spectrum_14']['lambda_min']:.4f}")

    # --- H_R8: all-24 per-trait check of H_R1-significant loci (from stored records)
    scal24 = {r["locus_id"]: r for r in rob["variant_records"]["full"]}
    n_mo24 = 0
    for lid in h_r1_sig:
        rec = scal24[lid]
        # stored records had scalar_p popped during summarize; recompute from matrices
        betas, ses, ancs = loci[lid]
        sp = [float(stats.chi2.sf(scalar_q(betas[:, j], ses[:, j]), len(ancs) - 1))
              for j in range(len(codes))]
        if not any(p < ALPHA / 24 for p in sp):
            n_mo24 += 1
    m24 = n_mo24 / len(h_r1_sig)
    out["H_R8"] = {"n_h_r1_significant": len(h_r1_sig), "n_multivariate_only_24": n_mo24,
                   "M24": round(m24, 4), "outcome": "holds" if m24 > 0.5 else "fails"}
    log(f"H_R8: M24 = {n_mo24}/{len(h_r1_sig)} = {m24:.4f} -> {out['H_R8']['outcome']}")

    # --- H_R4': non-EUR, 14-trait
    entries = []
    rows_for_cal = []
    for lid, (betas, ses, ancs) in loci.items():
        keep = [i for i, a in enumerate(ancs) if a != "EUR"]
        if len(keep) < 3:
            continue
        b = betas[np.ix_(keep, nonred_idx)]
        s = ses[np.ix_(keep, nonred_idx)]
        K = len(keep)
        Q, _ = qv(b, s, rho_nr_inv)
        p = float(stats.chi2.sf(Q, len(nonred_idx) * (K - 1)))
        sp = [float(stats.chi2.sf(scalar_q(b[:, j], s[:, j]), K - 1))
              for j in range(len(nonred_idx))]
        entries.append({"locus_id": lid, "p": p, "scalar_p": sp, "d": len(nonred_idx)})
        rows_for_cal.append((s, K))
    summ, _ = classify_variant(entries)
    cal = calibrate(rows_for_cal[:600], L_nr, rho_nr_inv, len(nonred_idx), 40, rng)
    gate = summ["n_fdr_significant"] >= H_R4_GATE
    outcome = ("unresolvable (G2 void)" if cal["pass"] is False
               else "unresolvable (gate)" if not gate
               else "holds" if summ["M"] and summ["M"] > 0.5 else "fails")
    out["H_R4p"] = {**summ, "calibration": cal, "outcome": outcome}
    log(f"H_R4': sig={summ['n_fdr_significant']} M={summ['M']} cal={cal['reject_rate']} -> {outcome}")

    # --- H_R5': discovery traits excluded, within 14
    entries = []
    for lid, (betas, ses, ancs) in loci.items():
        disc = {j for j, c in enumerate(codes)
                if (e := data[lid][c].get("EUR")) and e["neglog10_pval"] is not None
                and e["neglog10_pval"] >= GWS_NEGLOG10}
        cols = [j for j in nonred_idx if j not in disc]
        if len(cols) < 2:
            continue
        r_sub = rho[np.ix_(cols, cols)]
        b = betas[:, cols]
        s = ses[:, cols]
        K = len(ancs)
        Q, _ = qv(b, s, np.linalg.inv(r_sub))
        p = float(stats.chi2.sf(Q, len(cols) * (K - 1)))
        sp = [float(stats.chi2.sf(scalar_q(b[:, j], s[:, j]), K - 1)) for j in range(len(cols))]
        entries.append({"locus_id": lid, "p": p, "scalar_p": sp, "d": len(cols)})
    summ5, _ = classify_variant(entries)
    out["H_R5p"] = {**summ5,
                    "outcome": "holds" if summ5["M"] and summ5["M"] > 0.5 else "fails"}
    log(f"H_R5': sig={summ5['n_fdr_significant']} M={summ5['M']} -> {out['H_R5p']['outcome']}")

    # --- H_R7': rotation on 14-trait vectors among H_R1-significant loci
    subset = [(lid, *loci[lid]) for lid in h_r1_sig if lid in loci]
    obs = float(np.median([theta_max_of(b[:, nonred_idx]) for _, b, s, a in subset]))
    pre = []
    for _, b, s, a in subset:
        s_nr = s[:, nonred_idx]
        K = len(a)
        W = [rho_nr_inv * np.outer(1 / s_nr[k], 1 / s_nr[k]) for k in range(K)]
        Wsum_inv = np.linalg.inv(np.sum(W, axis=0))
        bbar = Wsum_inv @ np.sum([W[k] @ b[k, nonred_idx] for k in range(K)], axis=0)
        pre.append((bbar, s_nr))
    null_med = []
    for _ in range(args.reps_theta):
        tms = [theta_max_of(bbar[None, :] + (rng.standard_normal(s.shape) @ L_nr.T) * s)
               for bbar, s in pre]
        null_med.append(float(np.median(tms)))
    n95 = float(np.percentile(null_med, 95))
    out["H_R7p"] = {"observed_median": round(obs, 2), "null_95th": round(n95, 2),
                    "null_mean": round(float(np.mean(null_med)), 2),
                    "n_loci": len(subset),
                    "outcome": "holds" if obs > n95 else "fails"}
    log(f"H_R7': obs={obs:.1f} null95={n95:.1f} -> {out['H_R7p']['outcome']}")

    # --- descriptive: moment-matched diagonal recalibration (full 24-trait diagonal)
    L_full = np.linalg.cholesky(rho)
    by_k = {}
    diag_records = {r["locus_id"]: r for r in rob["variant_records"]["diagonal"]}
    for lid, (betas, ses, ancs) in loci.items():
        K = len(ancs)
        z = rng.standard_normal((25, K, len(codes)))
        eps = np.einsum("dj,rkj->rkd", L_full, z) * ses[None]
        for r in range(25):
            e = eps[r]
            Q = sum(scalar_q(e[:, j] , ses[:, j]) for j in range(len(codes)))
            by_k.setdefault(K, []).append(Q)
    params = {K: (np.var(v) / (2 * np.mean(v)), np.mean(v) ** 2 / (np.var(v) / 2))
              for K, v in by_k.items()}
    ps = []
    lids = []
    for lid, (betas, ses, ancs) in loci.items():
        c, nu = params[len(ancs)]
        ps.append(float(stats.chi2.sf(diag_records[lid]["Q"] / c, nu)))
        lids.append(lid)
    sig = bh_significant(ps)
    out["diagonal_recalibrated"] = {"n_fdr_significant": int(sig.sum()),
                                     "pct": round(100 * sig.sum() / len(ps), 1)}
    log(f"diagonal recalibrated: {int(sig.sum())} significant ({out['diagonal_recalibrated']['pct']}%)")

    # --- exploratory: 13-trait panel without WBC
    idx13 = [j for j in nonred_idx if codes[j] != "30000"]
    r13 = rho[np.ix_(idx13, idx13)]
    r13_inv = np.linalg.inv(r13)
    entries = []
    for lid, (betas, ses, ancs) in loci.items():
        b = betas[:, idx13]
        s = ses[:, idx13]
        K = len(ancs)
        Q, _ = qv(b, s, r13_inv)
        p = float(stats.chi2.sf(Q, len(idx13) * (K - 1)))
        sp = [float(stats.chi2.sf(scalar_q(b[:, j], s[:, j]), K - 1)) for j in range(len(idx13))]
        entries.append({"locus_id": lid, "p": p, "scalar_p": sp, "d": len(idx13)})
    summ13, _ = classify_variant(entries)
    ev13 = np.linalg.eigvalsh(r13)
    out["exploratory_13"] = {**summ13, "condition_number": float(ev13.max() / ev13.min())}
    log(f"13-trait (no WBC): sig={summ13['n_fdr_significant']} M={summ13['M']} "
        f"cond={out['exploratory_13']['condition_number']:.1f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    log(f"saved {args.out}")


if __name__ == "__main__":
    main()

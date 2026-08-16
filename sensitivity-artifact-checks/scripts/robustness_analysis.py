"""Robustness analysis for the multivariate-only heterogeneity claim.

Steps 2-3 and 6-7 of the frozen registration (sensitivity-artifact-checks/PREREG.md):
G1 fidelity gate, eigen recompute, Q_V variants (H_R1/H_R2/H_R4/H_R5 + diagonal),
per-variant BH and M, rotation angles. Bootstrap (G2/H_R3/H_R7) runs separately.

    uv run python sensitivity-artifact-checks/scripts/robustness_analysis.py
"""

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats

GWS_NEGLOG10 = 7.30103  # -log10(5e-8)
ALPHA = 0.05
FDR_Q = 0.05
TRACE_FRACTION = 0.95

# Non-redundant panel rule from the registration, keyed by UKB field code
# (stored trait names for 30090/30110 are swapped relative to UKB fields;
# 30090 is Platelet crit, 30110 is Platelet distribution width).
DROP_CODES = {
    "30030",  # Haematocrit
    "30040",  # Mean corpuscular volume
    "30050",  # Mean corpuscular haemoglobin
    "30060",  # Mean corpuscular haemoglobin concentration
    "30090",  # Platelet crit
    "30180", "30190", "30200", "30210",  # differential percentages in panel
    "30240",  # Reticulocyte percentage
}


def log(msg):
    print(f"[{datetime.now(timezone.utc).isoformat(timespec='seconds')}] {msg}", flush=True)


def load_shards(shard_dir, trait_codes):
    """Return locus_id -> {code -> {anc -> {beta, se, neglog10_pval}}}."""
    data = {}
    for code in trait_codes:
        path = shard_dir / f"{code}.jsonl"
        if not path.exists():
            raise FileNotFoundError(f"missing shard {path}")
        with open(path) as f:
            for line in f:
                rec = json.loads(line)
                data.setdefault(rec["locus_id"], {})[code] = rec["ancestries"]
    return data


def build_matrices(locus_traits, trait_codes, ancestries):
    """(K, d) beta and SE matrices in trait_codes order; None if incomplete."""
    K, d = len(ancestries), len(trait_codes)
    betas = np.zeros((K, d))
    ses = np.zeros((K, d))
    for j, code in enumerate(trait_codes):
        per_anc = locus_traits.get(code)
        if per_anc is None:
            return None
        for i, anc in enumerate(ancestries):
            entry = per_anc.get(anc)
            if entry is None:
                return None
            betas[i, j] = entry["beta"]
            ses[i, j] = entry["se"]
    return betas, ses


def q_v(betas, ses, weight_builder):
    """Generalized multivariate Q with per-stratum weight matrices."""
    K, d = betas.shape
    W = [weight_builder(ses[k]) for k in range(K)]
    W_sum = np.sum(W, axis=0)
    beta_bar = np.linalg.solve(W_sum, np.sum([W[k] @ betas[k] for k in range(K)], axis=0))
    Q = float(sum((betas[k] - beta_bar) @ W[k] @ (betas[k] - beta_bar) for k in range(K)))
    return Q, beta_bar


def full_weight(rho_inv):
    def build(se):
        Dinv = 1.0 / se
        return rho_inv * np.outer(Dinv, Dinv)
    return build


def scalar_q(betas_j, ses_j):
    w = 1.0 / ses_j**2
    bbar = np.sum(w * betas_j) / np.sum(w)
    return float(np.sum(w * (betas_j - bbar) ** 2))


def variant_stats(betas, ses, rho_inv, k_df):
    """Q_V, df, p, and per-trait scalar Q p-values for one variant at one locus."""
    K, d = betas.shape
    Q, _ = q_v(betas, ses, full_weight(rho_inv))
    df = k_df * (K - 1)
    p = float(stats.chi2.sf(Q, df))
    scalar_p = [float(stats.chi2.sf(scalar_q(betas[:, j], ses[:, j]), K - 1)) for j in range(d)]
    return Q, df, p, scalar_p


def classify(vector_fdr_sig, scalar_p, d):
    any_scalar = any(p < ALPHA / d for p in scalar_p)
    if vector_fdr_sig and not any_scalar:
        return "multivariate_only"
    if vector_fdr_sig and any_scalar:
        return "concordant_significant"
    if any_scalar:
        return "component_only"
    return "concordant_null"


def bh_significant(pvals, q=FDR_Q):
    """Boolean array: Benjamini-Hochberg at level q."""
    p = np.asarray(pvals)
    n = len(p)
    order = np.argsort(p)
    thresh = q * (np.arange(1, n + 1)) / n
    passed = p[order] <= thresh
    k = np.max(np.nonzero(passed)[0]) + 1 if passed.any() else 0
    sig = np.zeros(n, dtype=bool)
    sig[order[:k]] = True
    return sig


def summarize_variant(records, name):
    ps = [r["p"] for r in records]
    sig = bh_significant(ps)
    n_sig = int(sig.sum())
    n_multi_only = 0
    for r, s in zip(records, sig):
        r["fdr_significant"] = bool(s)
        r["category"] = classify(s, r.pop("scalar_p"), r["d"])
        if r["category"] == "multivariate_only":
            n_multi_only += 1
    M = n_multi_only / n_sig if n_sig else float("nan")
    return {
        "variant": name,
        "n_loci": len(records),
        "n_fdr_significant": n_sig,
        "n_multivariate_only": n_multi_only,
        "M": round(M, 4) if n_sig else None,
    }


def rotations(betas, ancestries):
    K = betas.shape[0]
    theta_max, pair = 0.0, None
    pairwise = {}
    for i in range(K):
        for j in range(i + 1, K):
            ni, nj = np.linalg.norm(betas[i]), np.linalg.norm(betas[j])
            if ni < 1e-12 or nj < 1e-12:
                continue
            theta = float(np.degrees(np.arccos(
                np.clip(betas[i] @ betas[j] / (ni * nj), -1, 1))))
            pairwise[f"{ancestries[i]}-{ancestries[j]}"] = round(theta, 2)
            if theta > theta_max:
                theta_max, pair = theta, f"{ancestries[i]}-{ancestries[j]}"
    return theta_max, pair, pairwise


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", type=Path, default=Path("sensitivity-artifact-checks/data/shards"))
    ap.add_argument("--corr", type=Path,
                    default=Path("analysis/data/phenotypic_correlations_estimated.json"))
    ap.add_argument("--v2", type=Path,
                    default=Path("analysis/results/genomewide_vector_q_v2_results.json"))
    ap.add_argument("--out", type=Path,
                    default=Path("sensitivity-artifact-checks/results/robustness_results.json"))
    ap.add_argument("--no-g1", action="store_true",
                    help="Skip the G1 gate (H_R3 reruns with alternative correlation matrices)")
    args = ap.parse_args()

    with open(args.corr) as f:
        corr = json.load(f)
    trait_codes = corr["trait_codes"]
    rho = np.array(corr["correlation_matrix"])
    with open(args.v2) as f:
        v2 = {r["locus_id"]: r for r in json.load(f)["results"]}

    log(f"panel: {len(trait_codes)} traits; v2 loci: {len(v2)}")
    data = load_shards(args.shards, trait_codes)
    log(f"recovered loci: {len(data)}")

    # --- eigen recompute and registered truncation rank
    evals, evecs = np.linalg.eigh(rho)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    k_rank = int(np.searchsorted(np.cumsum(evals), TRACE_FRACTION * rho.shape[0]) + 1)
    eigen_report = {
        "condition_number": float(evals[0] / evals[-1]),
        "lambda_min": float(evals[-1]),
        "n_eigenvalues_below_0.05": int((evals < 0.05).sum()),
        "truncation_rank_k": k_rank,
    }
    log(f"eigen: cond={eigen_report['condition_number']:.0f} "
        f"lambda_min={eigen_report['lambda_min']:.5f} k={k_rank}")

    rho_inv = np.linalg.inv(rho)
    rho_trunc_inv = evecs[:, :k_rank] @ np.diag(1.0 / evals[:k_rank]) @ evecs[:, :k_rank].T

    # --- G1 fidelity gate
    g1_total = g1_ok = 0
    g1_worst = 0.0
    missing = []
    matrices = {}
    for lid, v2rec in v2.items():
        locus_traits = data.get(lid)
        built = build_matrices(locus_traits, trait_codes, v2rec["ancestries"]) if locus_traits else None
        if built is None:
            missing.append(lid)
            continue
        betas, ses = built
        matrices[lid] = (betas, ses, v2rec["ancestries"])
        if args.no_g1:
            continue
        Q = variant_stats(betas, ses, rho_inv, len(trait_codes))[0]
        g1_total += 1
        rel = abs(Q - v2rec["vector_Q_V"]) / max(v2rec["vector_Q_V"], 1e-12)
        g1_worst = max(g1_worst, rel)
        if rel <= 0.01:
            g1_ok += 1
    # Amendment 1 (frozen 10fbfba26adc): G1b is value fidelity on the evaluation
    # set; coverage is a descriptive data note (G1a), not a gate.
    g1_frac = g1_ok / g1_total if g1_total else 0.0
    g1_pass = args.no_g1 or g1_frac >= 0.99
    g1_report = {
        "n_v2_loci": len(v2),
        "n_rebuilt": g1_total,
        "n_missing": len(missing),
        "frac_within_1pct": round(g1_frac, 5),
        "worst_rel_diff": round(g1_worst, 5),
        "pass": bool(g1_pass),
    }
    log(f"G1: rebuilt {g1_total}/{len(v2)}, within 1%: {g1_frac:.4f}, pass={g1_pass}")

    out = {
        "metadata": {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "registration": "sensitivity-artifact-checks/PREREG.md",
        },
        "eigen": eigen_report,
        "G1": g1_report,
        "missing_loci": missing[:50],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if not g1_pass:
        with open(args.out, "w") as f:
            json.dump(out, f, indent=1)
        log("G1 FAILED — all predictions void per registration. Stopping.")
        sys.exit(2)

    # --- variant definitions
    nonred_idx = [j for j, c in enumerate(trait_codes) if c not in DROP_CODES]
    rho_nr = rho[np.ix_(nonred_idx, nonred_idx)]
    rho_nr_inv = np.linalg.inv(rho_nr)
    log(f"non-redundant panel: {len(nonred_idx)} traits")

    variant_records = {v: [] for v in
                       ["full", "nonredundant", "rank_truncated", "non_eur", "no_discovery", "diagonal"]}
    locus_out = []

    for lid, (betas, ses, ancs) in matrices.items():
        d = len(trait_codes)
        eur_neglog10 = {}
        for j, code in enumerate(trait_codes):
            e = data[lid][code].get("EUR")
            eur_neglog10[code] = e["neglog10_pval"] if e and e["neglog10_pval"] is not None else 0.0
        disc_idx = [j for j, c in enumerate(trait_codes) if eur_neglog10[c] >= GWS_NEGLOG10]

        theta_max, pair, pairwise = rotations(betas, ancs)
        rec_common = {"locus_id": lid, "K": len(ancs),
                      "theta_max": round(theta_max, 2), "theta_max_pair": pair,
                      "n_discovery_traits": len(disc_idx)}

        def add(variant, b, s, r_inv, k_df, anc_note=None):
            if b.shape[0] < 3 or b.shape[1] < 2:
                return
            Q, df, p, scalar_p = variant_stats(b, s, r_inv, k_df)
            variant_records[variant].append(
                {"locus_id": lid, "Q": round(Q, 3), "df": df, "p": p,
                 "d": b.shape[1], "K": b.shape[0], "scalar_p": scalar_p})

        add("full", betas, ses, rho_inv, d)
        add("nonredundant", betas[:, nonred_idx], ses[:, nonred_idx], rho_nr_inv, len(nonred_idx))
        add("rank_truncated", betas, ses, rho_trunc_inv, k_rank)
        add("diagonal", betas, ses, np.eye(d), d)
        if "EUR" in ancs:
            keep = [i for i, a in enumerate(ancs) if a != "EUR"]
            add("non_eur", betas[keep], ses[keep], rho_inv, d)
        keep_t = [j for j in range(d) if j not in disc_idx]
        if len(keep_t) >= 2:
            r_sub = rho[np.ix_(keep_t, keep_t)]
            add("no_discovery", betas[:, keep_t], ses[:, keep_t],
                np.linalg.inv(r_sub), len(keep_t))

        locus_out.append(rec_common)

    summaries = {}
    for name, records in variant_records.items():
        if records:
            summaries[name] = summarize_variant(records, name)
            log(f"{name}: n={summaries[name]['n_loci']} "
                f"sig={summaries[name]['n_fdr_significant']} M={summaries[name]['M']}")

    out["variant_summaries"] = summaries
    out["variant_records"] = {k: v for k, v in variant_records.items()}
    out["locus_rotations"] = locus_out
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    log(f"saved {args.out}")


if __name__ == "__main__":
    main()

"""Re-run ACKR1/Duffy positive control with the 24-trait genome-wide methodology.

Matches the genome-wide v2 analysis exactly: same 24 traits, same estimated
correlation matrix, same Q_V computation. Queries rs2814778 (chr1:159174683)
directly rather than using the clumped variant.
"""

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from cross_ancestry_heterogeneity.vector_q import (
    scalar_cochran_q,
    vector_sheaf_q_correlated,
    vector_sheaf_q_diagonal,
)

PANUKB_DATA = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files"
PANUKB_IDX = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files_tabix"
ANCESTRIES = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]

ACKR1_CHROM = "1"
ACKR1_POS = 159174683

COLUMNS = [
    "chr", "pos", "ref", "alt",
    "af_meta_hq", "beta_meta_hq", "se_meta_hq",
    "neglog10_pval_meta_hq", "neglog10_pval_heterogeneity_hq",
    "af_meta", "beta_meta", "se_meta",
    "neglog10_pval_meta", "neglog10_pval_heterogeneity",
    "af_AFR", "af_AMR", "af_CSA", "af_EAS", "af_EUR", "af_MID",
    "beta_AFR", "beta_AMR", "beta_CSA", "beta_EAS", "beta_EUR", "beta_MID",
    "se_AFR", "se_AMR", "se_CSA", "se_EAS", "se_EUR", "se_MID",
    "neglog10_pval_AFR", "neglog10_pval_AMR", "neglog10_pval_CSA",
    "neglog10_pval_EAS", "neglog10_pval_EUR", "neglog10_pval_MID",
    "low_confidence_AFR", "low_confidence_AMR", "low_confidence_CSA",
    "low_confidence_EAS", "low_confidence_EUR", "low_confidence_MID",
]


def panukb_url(pheno_code):
    data = f"{PANUKB_DATA}/continuous-{pheno_code}-both_sexes-irnt.tsv.bgz"
    idx = f"{PANUKB_IDX}/continuous-{pheno_code}-both_sexes-irnt.tsv.bgz.tbi"
    return f"{data}##idx##{idx}"


def tabix_query(url, chrom, pos):
    region = f"{chrom}:{pos}-{pos}"
    result = subprocess.run(
        ["tabix", url, region],
        capture_output=True, text=True, timeout=60,
    )
    if result.returncode != 0 or not result.stdout.strip():
        return None
    for line in result.stdout.strip().split("\n"):
        fields = line.split("\t")
        if len(fields) >= len(COLUMNS) and int(fields[1]) == pos:
            return dict(zip(COLUMNS, fields))
    return None


def main():
    corr_path = REPO / "analysis" / "data" / "phenotypic_correlations_estimated.json"
    if not corr_path.exists():
        print(f"ERROR: Correlation matrix not found at {corr_path}")
        print("Run estimate_pheno_corr.py first.")
        sys.exit(1)

    with open(corr_path) as f:
        corr_data = json.load(f)

    corr_matrix = np.array(corr_data["correlation_matrix"])
    trait_codes = corr_data["trait_codes"]
    trait_names = corr_data["trait_names"]
    d = len(trait_codes)

    print(f"Loaded {d}x{d} correlation matrix ({d} traits)")
    print(f"Traits: {trait_names}")
    print(f"\nFetching rs2814778 (chr{ACKR1_CHROM}:{ACKR1_POS}) across {d} traits...")

    trait_data = {}
    for code, name in zip(trait_codes, trait_names):
        url = panukb_url(code)
        try:
            row = tabix_query(url, ACKR1_CHROM, ACKR1_POS)
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT: {name}")
            continue
        if row is None:
            print(f"  MISSING: {name}")
            continue

        anc_data = {}
        for anc in ANCESTRIES:
            try:
                beta = float(row[f"beta_{anc}"])
                se = float(row[f"se_{anc}"])
            except (ValueError, KeyError):
                continue
            if np.isnan(beta) or np.isnan(se) or se <= 0:
                continue
            anc_data[anc] = {"beta": beta, "se": se}

        if anc_data:
            trait_data[code] = {"name": name, "ancestries": anc_data}
            print(f"  {name}: {len(anc_data)} ancestries")

    available_codes = [c for c in trait_codes if c in trait_data]
    available_names = [trait_data[c]["name"] for c in available_codes]
    d_avail = len(available_codes)
    print(f"\nAvailable traits: {d_avail}/{d}")

    all_anc = set()
    for tc in available_codes:
        all_anc.update(trait_data[tc]["ancestries"].keys())

    good_anc = []
    for anc in sorted(all_anc):
        n_avail = sum(1 for tc in available_codes if anc in trait_data[tc]["ancestries"])
        if n_avail == d_avail:
            good_anc.append(anc)

    K = len(good_anc)
    print(f"Ancestries with complete data: {K} ({good_anc})")

    betas = np.zeros((K, d_avail))
    ses = np.zeros((K, d_avail))
    for j, tc in enumerate(available_codes):
        for i, anc in enumerate(good_anc):
            ad = trait_data[tc]["ancestries"][anc]
            betas[i, j] = ad["beta"]
            ses[i, j] = ad["se"]

    trait_idx = [trait_codes.index(c) for c in available_codes]
    rho = corr_matrix[np.ix_(trait_idx, trait_idx)]

    vec_corr = vector_sheaf_q_correlated(betas, ses, rho)
    vec_diag = vector_sheaf_q_diagonal(betas, ses)

    alpha = 0.05
    alpha_bonf = alpha / d_avail

    print(f"\n{'='*60}")
    print(f"ACKR1/Duffy (rs2814778) — 24-trait reanalysis")
    print(f"{'='*60}")
    print(f"K = {K} ancestries, d = {d_avail} traits")
    print(f"Multivariate Q_V (full cov): {vec_corr['Q_V']:.1f}, df={vec_corr['df']}, p={vec_corr['p']:.2e}")
    print(f"Multivariate Q_V (diagonal): {vec_diag['Q_V']:.1f}, df={vec_diag['df']}, p={vec_diag['p']:.2e}")
    print(f"Ratio full/diagonal: {vec_corr['Q_V']/vec_diag['Q_V']:.2f}")
    print(f"\nBonferroni threshold: {alpha_bonf:.4f}")

    scalar_results = {}
    sig_traits = []
    print(f"\nPer-trait scalar Q results:")
    for j, (tc, name) in enumerate(zip(available_codes, available_names)):
        r = scalar_cochran_q(betas[:, j], ses[:, j])
        scalar_results[name] = {
            "Q": round(r["Q"], 1),
            "df": r["df"],
            "p": float(r["p"]),
            "significant": bool(r["p"] < alpha_bonf),
        }
        sig = "***" if r["p"] < alpha_bonf else ""
        print(f"  {name:45s}: Q={r['Q']:8.1f}, p={r['p']:.2e} {sig}")
        if r["p"] < alpha_bonf:
            sig_traits.append(name)

    n_sig = len(sig_traits)
    vector_sig = vec_corr["p"] < alpha

    if vector_sig and n_sig > 0:
        category = "concordant_heterogeneous"
    elif vector_sig and n_sig == 0:
        category = "multivariate_only"
    elif not vector_sig and n_sig > 0:
        category = "marginal_only"
    else:
        category = "concordant_null"

    print(f"\nSignificant traits ({n_sig}): {sig_traits}")
    print(f"Category: {category}")

    results = {
        "locus": "ACKR1/Duffy (rs2814778)",
        "chrom": ACKR1_CHROM,
        "pos": ACKR1_POS,
        "K": K,
        "d": d_avail,
        "ancestries": good_anc,
        "traits": available_names,
        "vector_Q_V_full": round(float(vec_corr["Q_V"]), 3),
        "vector_Q_V_diagonal": round(float(vec_diag["Q_V"]), 3),
        "vector_df": int(vec_corr["df"]),
        "vector_p_full": float(vec_corr["p"]),
        "vector_p_diagonal": float(vec_diag["p"]),
        "ratio_full_over_diagonal": round(float(vec_corr["Q_V"]) / float(vec_diag["Q_V"]), 3),
        "n_scalar_significant": n_sig,
        "significant_traits": sig_traits,
        "bonferroni_threshold": alpha_bonf,
        "category": category,
        "scalar_results": scalar_results,
        "timestamp": datetime.now().isoformat(),
        "method": "24-trait genome-wide methodology with estimated EUR correlation matrix",
    }

    out_path = REPO / "analysis" / "results" / "ackr1_24trait_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()

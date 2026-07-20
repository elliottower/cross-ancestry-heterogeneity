"""Run vector Q analysis on Pan-UKB blood trait data at target loci.

Tests pre-registered predictions H_V1 through H_V6.
Outputs results to experiments/batch6_vector_extension/results/.

Usage:
    uv run python experiments/batch6_vector_extension/scripts/run_vector_analysis.py
"""

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vector_q import (
    format_p,
    max_rotation,
    rotation_angle,
    scalar_cochran_q,
    vector_sheaf_q_correlated,
    vector_sheaf_q_diagonal,
)

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "data"
RESULTS_DIR = REPO / "results"


def load_data():
    locus_path = DATA_DIR / "panukb_locus_data.json"
    if not locus_path.exists():
        raise FileNotFoundError(f"Run fetch_panukb.py first: {locus_path}")
    with open(locus_path) as f:
        locus_data = json.load(f)

    corr_path = DATA_DIR / "phenotypic_correlations_estimated.json"
    if corr_path.exists():
        with open(corr_path) as f:
            corr_data = json.load(f)
        return locus_data, corr_data
    return locus_data, None


def build_rho_submatrix(locus_traits, corr_data):
    if corr_data is None:
        return np.eye(len(locus_traits))

    corr_names = corr_data["trait_names"]
    rho_full = np.array(corr_data["correlation_matrix"])

    d = len(locus_traits)
    rho = np.eye(d)
    for i, t_i in enumerate(locus_traits):
        for j, t_j in enumerate(locus_traits):
            if i == j:
                continue
            if t_i in corr_names and t_j in corr_names:
                ii = corr_names.index(t_i)
                jj = corr_names.index(t_j)
                rho[i, j] = rho_full[ii][jj]
    return rho


def analyze_locus(rsid, locus, corr_data, alpha=0.05):
    beta_matrix = np.array(locus["betas"])
    se_matrix = np.array(locus["ses"])
    K, d = beta_matrix.shape
    ancestries = locus["ancestries"]
    traits = locus["traits"]

    rho = build_rho_submatrix(traits, corr_data)

    vec_corr = vector_sheaf_q_correlated(beta_matrix, se_matrix, rho)
    vec_diag = vector_sheaf_q_diagonal(beta_matrix, se_matrix)

    scalar_results = {}
    for j, trait in enumerate(traits):
        r = scalar_cochran_q(beta_matrix[:, j], se_matrix[:, j])
        scalar_results[trait] = r

    rot = max_rotation(beta_matrix)
    pairwise_rotations = {}
    for i in range(K):
        for j2 in range(i + 1, K):
            theta = rotation_angle(beta_matrix[i], beta_matrix[j2])
            key = f"{ancestries[i]}-{ancestries[j2]}"
            pairwise_rotations[key] = round(theta, 2)

    alpha_bonf = alpha / d
    any_scalar_sig = any(r["p"] < alpha_bonf for r in scalar_results.values())
    n_scalar_sig = sum(1 for r in scalar_results.values() if r["p"] < alpha_bonf)
    vector_sig = vec_corr["p"] < alpha

    if vector_sig and any_scalar_sig:
        category = "concordant_significant"
    elif not vector_sig and not any_scalar_sig:
        category = "concordant_null"
    elif vector_sig and not any_scalar_sig:
        category = "rotation_only"
    else:
        category = "component_only"

    return {
        "rsid": rsid,
        "gene": locus.get("gene", ""),
        "description": locus.get("description", ""),
        "K": K,
        "d": d,
        "ancestries": ancestries,
        "traits": traits,
        "category": category,
        "vector_Q_V": round(vec_corr["Q_V"], 3),
        "vector_df": vec_corr["df"],
        "vector_p": vec_corr["p"],
        "vector_Q_V_diagonal": round(vec_diag["Q_V"], 3),
        "sum_scalar_Q": round(vec_diag["sum_scalar_Q"], 3),
        "ratio_QV_over_sumQ": round(
            vec_corr["Q_V"] / max(vec_diag["sum_scalar_Q"], 1e-12), 3
        ),
        "n_scalar_significant": n_scalar_sig,
        "scalar_results": {
            trait: {
                "Q": round(r["Q"], 3),
                "p": r["p"],
                "I2": round(r["I2"], 3),
                "significant": r["p"] < alpha_bonf,
            }
            for trait, r in scalar_results.items()
        },
        "theta_max": round(rot["theta_max"], 2),
        "theta_max_pair": f"{ancestries[rot['pair'][0]]}-{ancestries[rot['pair'][1]]}",
        "pairwise_rotations": pairwise_rotations,
    }


def test_predictions(results):
    tests = {}

    ackr1 = next((r for r in results if r["rsid"] == "rs2814778"), None)
    if ackr1:
        tests["H_V1"] = {
            "prediction": "ACKR1: vector Q_V significant (p<0.001) AND theta_max > 30 degrees",
            "Q_V": ackr1["vector_Q_V"],
            "Q_V_p": ackr1["vector_p"],
            "theta_max": ackr1["theta_max"],
            "theta_max_pair": ackr1["theta_max_pair"],
            "pass": ackr1["vector_p"] < 0.001 and ackr1["theta_max"] > 30,
        }
    else:
        tests["H_V1"] = {"prediction": "ACKR1 not in data", "pass": None}

    n_total = len(results)
    bonf_alpha = 0.05 / n_total if n_total > 0 else 0.05
    het_loci = [r for r in results if r["vector_p"] < bonf_alpha]
    rotation_only = [r for r in het_loci if r["category"] == "rotation_only"]
    pct = len(rotation_only) / max(len(het_loci), 1) * 100
    tests["H_V2"] = {
        "prediction": ">= 5% of heterogeneous loci are rotation-only",
        "n_het_loci": len(het_loci),
        "n_rotation_only": len(rotation_only),
        "rotation_only_names": [f"{r['rsid']}({r['gene']})" for r in rotation_only],
        "pct": round(pct, 1),
        "pass": pct >= 5 if het_loci else None,
    }

    het_ratios = [r["ratio_QV_over_sumQ"] for r in results if r["vector_p"] < 0.05]
    if het_ratios:
        median_ratio = float(np.median(het_ratios))
        tests["H_V3"] = {
            "prediction": "Median Q_V/sum(Q_scalar) > 1.0 among significant loci",
            "median_ratio": round(median_ratio, 3),
            "all_ratios": [round(r, 3) for r in sorted(het_ratios)],
            "n_significant": len(het_ratios),
            "pass": median_ratio > 1.0,
        }
    else:
        tests["H_V3"] = {"prediction": "No significant vector Q loci", "pass": None}

    if ackr1:
        alpha_bonf = 0.05 / ackr1["d"]
        scalar_ps = {t: r["p"] for t, r in ackr1["scalar_results"].items()}
        sig_traits = [t for t, p in scalar_ps.items() if p < alpha_bonf]
        nonsig_traits = [t for t, p in scalar_ps.items() if p >= alpha_bonf]
        tests["H_V5"] = {
            "prediction": "ACKR1: neutrophil+monocyte scalar Q sig, >= 5 traits non-sig",
            "n_scalar_sig": len(sig_traits),
            "sig_traits": sig_traits,
            "n_scalar_nonsig": len(nonsig_traits),
            "pass": len(nonsig_traits) >= 5,
        }
    else:
        tests["H_V5"] = {"prediction": "ACKR1 not in data", "pass": None}

    component_only = [r for r in results if r["category"] == "component_only"]
    het_total = [
        r for r in results
        if r["category"] in ("concordant_significant", "rotation_only", "component_only")
    ]
    pct_co = len(component_only) / max(len(het_total), 1) * 100
    tests["H_V6"] = {
        "prediction": "Component-only < 2% of heterogeneous loci",
        "n_component_only": len(component_only),
        "component_only_names": [f"{r['rsid']}({r['gene']})" for r in component_only],
        "n_het_total": len(het_total),
        "pct": round(pct_co, 1),
        "pass": pct_co < 2 if het_total else None,
    }

    return tests


def print_ackr1_detail(result):
    print(f"\n{'=' * 70}")
    print("ACKR1/DUFFY LOCUS DETAIL (for paper Section 5.3)")
    print(f"{'=' * 70}")
    print(f"Variant: {result['rsid']}")
    print(f"Ancestries ({result['K']}): {', '.join(result['ancestries'])}")
    print(f"Traits: {result['d']} blood cell traits")
    print(f"\nVector Q_V = {result['vector_Q_V']:.1f} (df={result['vector_df']}, "
          f"p={format_p(result['vector_p'])})")
    print(f"Vector Q_V (diagonal) = {result['vector_Q_V_diagonal']:.1f}")
    print(f"Sum scalar Q = {result['sum_scalar_Q']:.1f}")
    print(f"Ratio Q_V / sum(Q) = {result['ratio_QV_over_sumQ']:.2f}")
    print(f"Max rotation = {result['theta_max']:.1f} degrees ({result['theta_max_pair']})")

    print(f"\nPairwise rotation angles:")
    for pair, theta in sorted(result['pairwise_rotations'].items()):
        print(f"  {pair}: {theta:.1f} degrees")

    alpha_bonf = 0.05 / result['d']
    print(f"\nPer-trait scalar Q (Bonferroni alpha = {alpha_bonf:.4f}):")
    for trait, r in sorted(result['scalar_results'].items(),
                           key=lambda x: x[1]['Q'], reverse=True):
        sig = "***" if r['p'] < 0.001 else "**" if r['p'] < 0.01 else "*" if r['p'] < alpha_bonf else ""
        print(f"  {trait:45s} Q={r['Q']:8.1f}  p={format_p(r['p']):>12s}  I2={r['I2']:.2f}  {sig}")


def main():
    print("Vector-Valued Sheaf Q Analysis — Pan-UKB Blood Traits")
    print("=" * 55)
    print(f"Timestamp: {datetime.now().isoformat()}")
    print()

    RESULTS_DIR.mkdir(exist_ok=True)

    locus_data, corr_data = load_data()
    loci = locus_data["loci"]
    n_corr = corr_data["n_shared_variants"] if corr_data else 0
    print(f"Loaded {len(loci)} loci")
    print(f"Phenotypic correlations: {'estimated from ' + str(n_corr) + ' null variants' if corr_data else 'DIAGONAL (identity)'}")
    print()

    results = []
    for rsid in tqdm(loci, desc="Analyzing loci"):
        result = analyze_locus(rsid, loci[rsid], corr_data)
        results.append(result)

    print(f"\n{'=' * 60}")
    print("LOCUS CLASSIFICATION SUMMARY")
    print(f"{'=' * 60}")
    categories = {}
    for r in results:
        categories[r["category"]] = categories.get(r["category"], 0) + 1
    for cat, n in sorted(categories.items()):
        print(f"  {cat:30s}: {n:3d} ({n / len(results) * 100:.1f}%)")
    print(f"  {'TOTAL':30s}: {len(results):3d}")

    print(f"\n{'=' * 60}")
    print("PER-LOCUS RESULTS")
    print(f"{'=' * 60}")
    for r in sorted(results, key=lambda x: x["vector_p"]):
        sig = "***" if r["vector_p"] < 0.001 else "**" if r["vector_p"] < 0.01 else "*" if r["vector_p"] < 0.05 else ""
        print(f"\n{r['rsid']} ({r['gene']}): {r['description']}")
        print(f"  Category: {r['category']}")
        print(f"  Vector Q_V = {r['vector_Q_V']:.1f}, p = {format_p(r['vector_p'])} "
              f"(df={r['vector_df']}) {sig}")
        print(f"  Q_V(diagonal) = {r['vector_Q_V_diagonal']:.1f}, "
              f"sum(Q_scalar) = {r['sum_scalar_Q']:.1f}, "
              f"ratio = {r['ratio_QV_over_sumQ']:.2f}")
        print(f"  Max rotation: {r['theta_max']:.1f} degrees ({r['theta_max_pair']})")
        print(f"  Scalar significant: {r['n_scalar_significant']}/{r['d']} traits")

    ackr1 = next((r for r in results if r["rsid"] == "rs2814778"), None)
    if ackr1:
        print_ackr1_detail(ackr1)

    print(f"\n{'=' * 60}")
    print("PRE-REGISTERED PREDICTION TESTS")
    print(f"{'=' * 60}")
    tests = test_predictions(results)
    for hid, test in sorted(tests.items()):
        status = "PASS" if test["pass"] else "FAIL" if test["pass"] is False else "N/A"
        marker = "[+]" if test["pass"] else "[-]" if test["pass"] is False else "[?]"
        print(f"\n  {marker} {hid}: {test['prediction']}")
        for k, v in test.items():
            if k not in ("prediction", "pass"):
                print(f"      {k}: {v}")
        print(f"      Result: {status}")

    output = {
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "n_loci": len(results),
            "source": locus_data["metadata"]["source"],
            "phenotypic_correlations": corr_data["method"] if corr_data else "diagonal",
            "n_corr_variants": n_corr,
        },
        "results": [
            {k: (float(v) if isinstance(v, (np.floating, np.integer)) else v)
             for k, v in r.items()}
            for r in results
        ],
        "prediction_tests": {
            k: {kk: (float(vv) if isinstance(vv, (np.floating, np.integer)) else vv)
                for kk, vv in v.items()}
            for k, v in tests.items()
        },
        "summary": {"categories": categories},
    }
    output_path = RESULTS_DIR / "vector_q_results.json"
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nResults saved to {output_path}")

    n_pass = sum(1 for t in tests.values() if t["pass"] is True)
    n_fail = sum(1 for t in tests.values() if t["pass"] is False)
    n_na = sum(1 for t in tests.values() if t["pass"] is None)
    print(f"\nPrediction summary: {n_pass} PASS, {n_fail} FAIL, {n_na} N/A out of {len(tests)}")


if __name__ == "__main__":
    main()

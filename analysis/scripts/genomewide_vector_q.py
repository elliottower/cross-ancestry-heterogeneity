"""Genome-wide vector Q scan on Pan-UKB pleiotropic blood-trait loci.

Replaces the 7-locus hand-picked analysis with a systematic genome-wide screen:
1. Scan all 27 blood trait EUR summary stats for GWS variants (p < 5e-8)
2. Identify pleiotropic loci: GWS in >= 2 traits
3. Distance-clump at 500kb (keep strongest EUR signal)
4. Fetch per-ancestry data via tabix for all 6 ancestries x 27 traits
5. Run vector Q_V and 25 scalar Q tests at each locus
6. Classify: concordant_significant / rotation_only / component_only / concordant_null

Pre-registration: frozen at the commit containing this file.

Pre-registered predictions (alpha = 0.05):
- H_GW1: ACKR1 (rs2814778) will be recovered and classified as concordant_significant
         or rotation_only, with theta_max > 30 degrees.
- H_GW2: All 7 original hand-picked loci will be recovered in the genome-wide scan.
- H_GW3: >= 10% of loci with significant vector Q_V are rotation-only.
- H_GW4: Median Q_V / sum(Q_scalar) > 1.0 among rotation-only loci.
- H_GW5: Total pleiotropic loci (GWS in >= 2 traits, after clumping) > 50.
- H_GW6: Component-only loci (scalar sig but vector non-sig) < 5% of all
         heterogeneous loci.

Two-stage pipeline:
  Stage 1 (Modal, ~30 min): Download EUR-only columns, discover pleiotropic loci,
           clump, save locus list.
  Stage 2 (local, ~2 hr): Tabix per-ancestry data at discovered loci, run vector Q.

Usage:
    # Stage 1: discover loci (Modal)
    modal run modal_genomewide_discover.py

    # Stage 2: analyze (local, after downloading locus list)
    uv run python genomewide_vector_q.py --stage analyze --loci data/genomewide_loci.json
"""

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from tqdm import tqdm

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "data"
RESULTS_DIR = REPO / "results"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vector_q import (
    format_p,
    max_rotation,
    scalar_cochran_q,
    vector_sheaf_q_correlated,
    vector_sheaf_q_diagonal,
)

ANCESTRIES = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]

PANUKB_DATA = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files"
PANUKB_IDX = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files_tabix"

BLOOD_TRAITS = [
    ("30000", "White blood cell count"),
    ("30010", "Red blood cell count"),
    ("30020", "Haemoglobin concentration"),
    ("30030", "Haematocrit"),
    ("30040", "Mean corpuscular volume"),
    ("30050", "Mean corpuscular haemoglobin"),
    ("30060", "Mean corpuscular haemoglobin concentration"),
    ("30070", "Red blood cell distribution width"),
    ("30080", "Platelet count"),
    ("30090", "Platelet distribution width"),
    ("30100", "Mean platelet volume"),
    ("30110", "Platelet crit"),
    ("30120", "Lymphocyte count"),
    ("30130", "Monocyte count"),
    ("30140", "Neutrophil count"),
    ("30150", "Eosinophil count"),
    ("30180", "Lymphocyte percentage"),
    ("30190", "Monocyte percentage"),
    ("30200", "Neutrophil percentage"),
    ("30210", "Eosinophil percentage"),
    ("30220", "Basophil percentage"),
    ("30240", "Reticulocyte percentage"),
    ("30250", "Reticulocyte count"),
    ("30260", "Mean reticulocyte volume"),
    ("30270", "Mean sphered cell volume"),
    ("30280", "Immature reticulocyte fraction"),
    ("30290", "High light scatter reticulocyte percentage"),
]

GWS_THRESHOLD = 7.301  # -log10(5e-8)
MIN_TRAITS_PLEIOTROPIC = 2
CLUMP_WINDOW_BP = 500_000
MIN_ANCESTRIES = 3

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

ORIGINAL_7_LOCI = {
    "rs2814778", "rs1354034", "rs3184504", "rs2476601",
    "rs855791", "rs9349379", "rs11065987",
}


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


def parse_ancestry_data(row):
    out = {}
    for anc in ANCESTRIES:
        try:
            beta = float(row[f"beta_{anc}"])
            se = float(row[f"se_{anc}"])
        except (ValueError, KeyError):
            continue
        if np.isnan(beta) or np.isnan(se) or se <= 0:
            continue
        out[anc] = {"beta": beta, "se": se}
    return out


def fetch_locus_data(loci, traits):
    """Fetch per-ancestry data for all loci via tabix."""
    results = {}
    total = len(loci) * len(traits)
    pbar = tqdm(total=total, desc="Fetching per-ancestry data")

    for locus in loci:
        chrom = locus["chr"]
        pos = locus["pos"]
        locus_id = locus["id"]
        trait_data = {}

        for pheno_code, pheno_name in traits:
            pbar.update(1)
            pbar.set_postfix_str(f"{locus_id}/{pheno_name[:15]}")
            url = panukb_url(pheno_code)
            try:
                row = tabix_query(url, chrom, pos)
            except subprocess.TimeoutExpired:
                continue
            if row is None:
                continue
            anc_data = parse_ancestry_data(row)
            if anc_data:
                trait_data[pheno_code] = {
                    "name": pheno_name,
                    "ancestries": anc_data,
                }

        if len(trait_data) >= MIN_TRAITS_PLEIOTROPIC:
            results[locus_id] = {**locus, "traits": trait_data}

    pbar.close()
    return results


def build_matrices(locus_data, traits):
    """Build (K, d) beta and SE matrices for a locus."""
    trait_codes = [t[0] for t in traits]
    available = [c for c in trait_codes if c in locus_data["traits"]]
    if len(available) < 2:
        return None

    all_anc = set()
    for tc in available:
        all_anc.update(locus_data["traits"][tc]["ancestries"].keys())

    good_anc = []
    for anc in sorted(all_anc):
        n_avail = sum(
            1 for tc in available
            if anc in locus_data["traits"][tc]["ancestries"]
        )
        if n_avail == len(available):
            good_anc.append(anc)

    if len(good_anc) < MIN_ANCESTRIES:
        return None

    K = len(good_anc)
    d = len(available)
    betas = np.zeros((K, d))
    ses = np.zeros((K, d))
    trait_names = []

    for j, tc in enumerate(available):
        trait_names.append(locus_data["traits"][tc]["name"])
        for i, anc in enumerate(good_anc):
            ad = locus_data["traits"][tc]["ancestries"][anc]
            betas[i, j] = ad["beta"]
            ses[i, j] = ad["se"]

    return {
        "betas": betas,
        "ses": ses,
        "ancestries": good_anc,
        "traits": trait_names,
        "trait_codes": available,
        "K": K,
        "d": d,
    }


def analyze_locus(locus_id, matrices, corr_matrix, alpha=0.05):
    betas = matrices["betas"]
    ses = matrices["ses"]
    K, d = matrices["K"], matrices["d"]
    ancestries = matrices["ancestries"]
    traits = matrices["traits"]

    rho = corr_matrix if corr_matrix is not None else np.eye(d)

    vec_corr = vector_sheaf_q_correlated(betas, ses, rho)
    vec_diag = vector_sheaf_q_diagonal(betas, ses)

    scalar_results = {}
    for j, trait in enumerate(traits):
        r = scalar_cochran_q(betas[:, j], ses[:, j])
        scalar_results[trait] = r

    rot = max_rotation(betas)
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
        "locus_id": locus_id,
        "K": K,
        "d": d,
        "ancestries": ancestries,
        "category": category,
        "vector_Q_V": round(float(vec_corr["Q_V"]), 3),
        "vector_df": int(vec_corr["df"]),
        "vector_p": float(vec_corr["p"]),
        "vector_Q_V_diagonal": round(float(vec_diag["Q_V"]), 3),
        "sum_scalar_Q": round(float(vec_diag["sum_scalar_Q"]), 3),
        "ratio_QV_over_sumQ": round(
            float(vec_corr["Q_V"]) / max(float(vec_diag["sum_scalar_Q"]), 1e-12), 3
        ),
        "n_scalar_significant": n_scalar_sig,
        "theta_max": round(float(rot["theta_max"]), 2),
        "theta_max_pair": f"{ancestries[rot['pair'][0]]}-{ancestries[rot['pair'][1]]}",
    }


def test_predictions(results, locus_metadata):
    tests = {}

    ackr1 = [r for r in results if "rs2814778" in r["locus_id"] or "ACKR1" in r.get("locus_id", "")]
    if not ackr1:
        ackr1 = [r for r in results if locus_metadata.get(r["locus_id"], {}).get("chr") == "1"
                 and abs(locus_metadata.get(r["locus_id"], {}).get("pos", 0) - 159174683) < CLUMP_WINDOW_BP]
    ackr1_r = ackr1[0] if ackr1 else None

    tests["H_GW1"] = {
        "prediction": "ACKR1 recovered, classified concordant/rotation-only, theta_max > 30",
        "found": ackr1_r is not None,
        "category": ackr1_r["category"] if ackr1_r else None,
        "theta_max": ackr1_r["theta_max"] if ackr1_r else None,
        "pass": (ackr1_r is not None
                 and ackr1_r["category"] in ("concordant_significant", "rotation_only")
                 and ackr1_r["theta_max"] > 30) if ackr1_r else False,
    }

    recovered = set()
    for r in results:
        lid = r["locus_id"]
        meta = locus_metadata.get(lid, {})
        rsid = meta.get("rsid", lid)
        if rsid in ORIGINAL_7_LOCI:
            recovered.add(rsid)
    tests["H_GW2"] = {
        "prediction": "All 7 original loci recovered",
        "recovered": sorted(recovered),
        "missing": sorted(ORIGINAL_7_LOCI - recovered),
        "n_recovered": len(recovered),
        "pass": len(recovered) == 7,
    }

    het_loci = [r for r in results if r["vector_p"] < 0.05]
    rotation_only = [r for r in het_loci if r["category"] == "rotation_only"]
    pct_ro = len(rotation_only) / max(len(het_loci), 1) * 100
    tests["H_GW3"] = {
        "prediction": ">= 10% of vector-sig loci are rotation-only",
        "n_het": len(het_loci),
        "n_rotation_only": len(rotation_only),
        "pct": round(pct_ro, 1),
        "pass": pct_ro >= 10 if het_loci else None,
    }

    ro_ratios = [r["ratio_QV_over_sumQ"] for r in rotation_only]
    if ro_ratios:
        med = float(np.median(ro_ratios))
        tests["H_GW4"] = {
            "prediction": "Median Q_V/sum(Q_scalar) > 1.0 among rotation-only",
            "median_ratio": round(med, 3),
            "n": len(ro_ratios),
            "pass": med > 1.0,
        }
    else:
        tests["H_GW4"] = {"prediction": "No rotation-only loci", "pass": None}

    tests["H_GW5"] = {
        "prediction": "Total pleiotropic loci > 50",
        "n_loci": len(results),
        "pass": len(results) > 50,
    }

    component_only = [r for r in results if r["category"] == "component_only"]
    pct_co = len(component_only) / max(len(het_loci), 1) * 100
    tests["H_GW6"] = {
        "prediction": "Component-only < 5% of heterogeneous loci",
        "n_component_only": len(component_only),
        "n_het": len(het_loci),
        "pct": round(pct_co, 1),
        "pass": pct_co < 5 if het_loci else None,
    }

    return tests


def run_analysis(loci_path, corr_path=None):
    """Stage 2: analyze discovered loci."""
    with open(loci_path) as f:
        loci_data = json.load(f)

    loci = loci_data["loci"]
    print(f"Loaded {len(loci)} pleiotropic loci from {loci_path}")

    corr_matrix = None
    if corr_path and Path(corr_path).exists():
        with open(corr_path) as f:
            cd = json.load(f)
        corr_matrix = np.array(cd["correlation_matrix"])
        print(f"Loaded phenotypic correlations ({corr_matrix.shape[0]} traits)")

    print(f"\nFetching per-ancestry data for {len(loci)} loci x {len(BLOOD_TRAITS)} traits...")
    locus_full = fetch_locus_data(list(loci.values()), BLOOD_TRAITS)
    print(f"Got data for {len(locus_full)} loci")

    results = []
    locus_metadata = {}
    for lid, ldata in tqdm(locus_full.items(), desc="Analyzing"):
        matrices = build_matrices(ldata, BLOOD_TRAITS)
        if matrices is None:
            continue
        result = analyze_locus(lid, matrices, corr_matrix)
        results.append(result)
        locus_metadata[lid] = ldata

    print(f"\n{'='*60}")
    print(f"GENOME-WIDE VECTOR Q RESULTS ({len(results)} loci)")
    print(f"{'='*60}")

    categories = {}
    for r in results:
        categories[r["category"]] = categories.get(r["category"], 0) + 1
    for cat, n in sorted(categories.items()):
        print(f"  {cat:30s}: {n:3d} ({n/len(results)*100:.1f}%)")
    print(f"  {'TOTAL':30s}: {len(results):3d}")

    tests = test_predictions(results, {r["locus_id"]: locus_metadata.get(r["locus_id"], {}) for r in results})

    print(f"\n{'='*60}")
    print("PRE-REGISTERED PREDICTIONS")
    print(f"{'='*60}")
    for hid, test in sorted(tests.items()):
        status = "PASS" if test["pass"] else "FAIL" if test["pass"] is False else "N/A"
        marker = "[+]" if test["pass"] else "[-]" if test["pass"] is False else "[?]"
        print(f"\n  {marker} {hid}: {test['prediction']}")
        for k, v in test.items():
            if k not in ("prediction", "pass"):
                print(f"      {k}: {v}")

    RESULTS_DIR.mkdir(exist_ok=True)
    output = {
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "source": "genome-wide Pan-UKB pleiotropic blood-trait loci",
            "n_loci_analyzed": len(results),
            "n_loci_discovered": len(loci),
        },
        "results": results,
        "prediction_tests": {
            k: {kk: (float(vv) if isinstance(vv, (np.floating, np.integer)) else vv)
                for kk, vv in v.items()}
            for k, v in tests.items()
        },
        "summary": {"categories": categories},
    }
    out_path = RESULTS_DIR / "genomewide_vector_q_results.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")

    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--loci", type=Path, default=DATA_DIR / "genomewide_loci.json")
    parser.add_argument("--corr", type=Path, default=DATA_DIR / "phenotypic_correlations_estimated.json")
    args = parser.parse_args()
    run_analysis(args.loci, args.corr)

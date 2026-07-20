"""Estimate phenotypic correlations between blood traits from GWAS z-scores.

Under the null (no genetic effect), z-scores at common variants are correlated
across traits in proportion to the phenotypic correlation * sample overlap.
Since all Pan-UKB traits use the same UK Biobank samples (100% overlap within
each ancestry), the z-score correlation at null loci directly estimates the
phenotypic correlation.

We sample ~500 common variants from chr22 (small, fast to query) for each
trait, compute z = beta/se, and estimate the Pearson correlation matrix.

Usage:
    uv run python experiments/batch6_vector_extension/scripts/estimate_pheno_corr.py
"""

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from tqdm import tqdm

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "data"

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
]

ANCESTRIES = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]

# chr22 region: 22:20000000-25000000 (5Mb, should have ~1000+ common variants)
SAMPLE_REGION = "22:20000000-25000000"


def panukb_url(pheno_code):
    data = f"{PANUKB_DATA}/continuous-{pheno_code}-both_sexes-irnt.tsv.bgz"
    idx = f"{PANUKB_IDX}/continuous-{pheno_code}-both_sexes-irnt.tsv.bgz.tbi"
    return f"{data}##idx##{idx}"


def tabix_region(url, region):
    result = subprocess.run(
        ["tabix", url, region],
        capture_output=True, text=True, timeout=120,
    )
    if result.returncode != 0 or not result.stdout.strip():
        return []
    return result.stdout.strip().split("\n")


def extract_zscores(lines, ancestry="EUR"):
    """Extract z-scores (beta/se) for a given ancestry from tabix output rows."""
    # Column indices (0-based):
    # beta_AFR=20, beta_AMR=21, beta_CSA=22, beta_EAS=23, beta_EUR=24, beta_MID=25
    # se_AFR=26, se_AMR=27, se_CSA=28, se_EAS=29, se_EUR=30, se_MID=31
    # af_EUR=18, low_confidence_EUR=40
    anc_idx = ANCESTRIES.index(ancestry)
    beta_col = 20 + anc_idx
    se_col = 26 + anc_idx
    af_col = 14 + anc_idx
    lc_col = 38 + anc_idx

    positions = []
    zscores = []
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 44:
            continue
        try:
            af = float(fields[af_col])
            beta = float(fields[beta_col])
            se = float(fields[se_col])
            lc = fields[lc_col].lower() == "true"
        except (ValueError, IndexError):
            continue
        if lc or se <= 0 or np.isnan(beta) or np.isnan(se):
            continue
        # Filter to common variants (MAF > 5%) to reduce noise
        maf = min(af, 1 - af)
        if maf < 0.05:
            continue
        z = beta / se
        # Filter out genome-wide significant hits (these are NOT null)
        if abs(z) > 4:
            continue
        pos = int(fields[1])
        positions.append(pos)
        zscores.append(z)

    return positions, zscores


def main():
    print("Estimating phenotypic correlations from EUR z-scores at null loci")
    print("=" * 65)
    print(f"Region: {SAMPLE_REGION}")
    print(f"Traits: {len(BLOOD_TRAITS)}")
    print(f"Ancestry: EUR (largest sample)")
    print(f"Timestamp: {datetime.now().isoformat()}")
    print()

    try:
        subprocess.run(["tabix", "--version"], capture_output=True, timeout=5, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        print("ERROR: tabix not found")
        sys.exit(1)

    # Fetch z-scores for each trait
    trait_positions = {}
    trait_zscores = {}

    for pheno_code, pheno_name in tqdm(BLOOD_TRAITS, desc="Fetching chr22 data"):
        url = panukb_url(pheno_code)
        try:
            lines = tabix_region(url, SAMPLE_REGION)
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT: {pheno_name}")
            continue
        positions, zscores = extract_zscores(lines, ancestry="EUR")
        trait_positions[pheno_code] = positions
        trait_zscores[pheno_code] = dict(zip(positions, zscores))
        print(f"  {pheno_name}: {len(positions)} common null variants")

    # Drop traits with too few variants
    min_variants = 100
    good_traits = [(c, n) for c, n in BLOOD_TRAITS if len(trait_zscores.get(c, {})) >= min_variants]
    dropped = [n for c, n in BLOOD_TRAITS if len(trait_zscores.get(c, {})) < min_variants]
    if dropped:
        print(f"\nDropped {len(dropped)} traits with <{min_variants} variants: {dropped}")
    print(f"Using {len(good_traits)} traits for correlation estimation")

    # Find positions shared across all good traits
    all_positions = None
    for pheno_code, _ in good_traits:
        pos_set = set(trait_zscores[pheno_code].keys())
        if all_positions is None:
            all_positions = pos_set
        else:
            all_positions = all_positions & pos_set

    shared = sorted(all_positions)
    print(f"Shared positions across all traits: {len(shared)}")

    if len(shared) < 50:
        print("ERROR: Too few shared positions to estimate correlations")
        sys.exit(1)

    # Build z-score matrix: (n_variants, n_traits)
    n = len(shared)
    d = len(good_traits)
    Z = np.zeros((n, d))
    for j, (pheno_code, _) in enumerate(good_traits):
        for i, pos in enumerate(shared):
            Z[i, j] = trait_zscores[pheno_code][pos]

    # Compute Pearson correlation
    rho = np.corrcoef(Z.T)
    print(f"Correlation matrix: {rho.shape}")
    print(f"Mean off-diagonal |r|: {np.mean(np.abs(rho[np.triu_indices(d, k=1)])):.3f}")
    print(f"Max off-diagonal |r|: {np.max(np.abs(rho[np.triu_indices(d, k=1)])):.3f}")

    # Print top correlations
    print("\nTop 15 trait-trait correlations:")
    pairs = []
    for i in range(d):
        for j in range(i + 1, d):
            pairs.append((abs(rho[i, j]), rho[i, j], good_traits[i][1], good_traits[j][1]))
    pairs.sort(reverse=True)
    for _, r, t1, t2 in pairs[:15]:
        print(f"  {r:+.3f}  {t1} <-> {t2}")

    # Save
    DATA_DIR.mkdir(exist_ok=True)
    output = {
        "method": "z-score correlation at common null variants (|z|<4, MAF>5%) on chr22",
        "ancestry": "EUR",
        "region": SAMPLE_REGION,
        "n_shared_variants": len(shared),
        "n_traits": d,
        "trait_codes": [t[0] for t in good_traits],
        "trait_names": [t[1] for t in good_traits],
        "correlation_matrix": rho.tolist(),
        "timestamp": datetime.now().isoformat(),
    }
    out_path = DATA_DIR / "phenotypic_correlations_estimated.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()

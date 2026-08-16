"""Re-estimate the cross-trait correlation matrix for H_R3.

Same protocol as analysis/scripts/estimate_pheno_corr.py (z-score correlation at
common null variants: MAF > 5%, |z| < 4, not low-confidence), parameterized by
region and ancestry, fixed to the registered 24-trait panel. If any panel trait
yields fewer than 100 usable variants the matrix is incomplete and the arm is
unresolvable; the output records this rather than silently shrinking the panel.

    uv run python sensitivity-artifact-checks/scripts/estimate_corr_variant.py \
        --region 1:20000000-25000000 --ancestry EUR --out .../corr_chr1_eur.json
"""

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PANUKB_DATA = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files"
PANUKB_IDX = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files_tabix"
ANCESTRIES = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]
MIN_VARIANTS = 100


def log(msg):
    print(f"[{datetime.now(timezone.utc).isoformat(timespec='seconds')}] {msg}", flush=True)


def panukb_url(code):
    return (f"{PANUKB_DATA}/continuous-{code}-both_sexes-irnt.tsv.bgz"
            f"##idx##{PANUKB_IDX}/continuous-{code}-both_sexes-irnt.tsv.bgz.tbi")


def zscores_for(code, region, ancestry):
    result = subprocess.run(["tabix", panukb_url(code), region],
                            capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        log(f"{code}: tabix failed: {result.stderr.strip()[:200]}")
        return {}
    out = {}
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 34:
            continue
        anc_start = len(fields) - 30
        i = ANCESTRIES.index(ancestry)
        try:
            af = float(fields[anc_start + i])
            beta = float(fields[anc_start + 6 + i])
            se = float(fields[anc_start + 12 + i])
        except ValueError:
            continue
        if fields[anc_start + 24 + i].lower() == "true" or se <= 0:
            continue
        if min(af, 1 - af) < 0.05:
            continue
        z = beta / se
        if abs(z) > 4:
            continue
        out[int(fields[1])] = z
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", required=True)
    ap.add_argument("--ancestry", required=True, choices=ANCESTRIES)
    ap.add_argument("--corr", type=Path,
                    default=Path("analysis/data/phenotypic_correlations_estimated.json"))
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    with open(args.corr) as f:
        base = json.load(f)
    trait_codes = base["trait_codes"]

    per_trait = {}
    for n, code in enumerate(trait_codes, 1):
        per_trait[code] = zscores_for(code, args.region, args.ancestry)
        log(f"[{n}/{len(trait_codes)}] {code}: {len(per_trait[code])} null variants")

    short = [c for c in trait_codes if len(per_trait[c]) < MIN_VARIANTS]
    complete = not short
    shared = None
    for code in trait_codes:
        s = set(per_trait[code])
        shared = s if shared is None else shared & s
    shared = sorted(shared or [])
    log(f"shared positions: {len(shared)}; traits below {MIN_VARIANTS}: {short}")

    out = {
        "method": f"z-score correlation at common null variants (|z|<4, MAF>5%) on {args.region}",
        "ancestry": args.ancestry,
        "region": args.region,
        "n_shared_variants": len(shared),
        "n_traits": len(trait_codes),
        "trait_codes": trait_codes,
        "traits_below_minimum": short,
        "complete": complete and len(shared) >= 50,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    if out["complete"]:
        Z = np.array([[per_trait[c][p] for c in trait_codes] for p in shared])
        out["correlation_matrix"] = np.corrcoef(Z.T).tolist()
    else:
        log("INCOMPLETE: matrix not estimable on the full panel under the protocol")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    log(f"saved {args.out} (complete={out['complete']})")


if __name__ == "__main__":
    main()

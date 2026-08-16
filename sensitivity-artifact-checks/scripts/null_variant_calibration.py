"""H_R9 (amendment 4, frozen 9fbed73f34cf): empirical calibration at null variants.

Samples one variant per 1-Mb genomic window (first variant in the window's opening
50 kb of the 30000 file), pulls all 24 traits at the sampled positions, excludes
variants with EUR -log10 p >= 4 for any trait, and computes 14- and 24-trait Q_V
against their chi-square references. Checkpointed per stage and per trait shard.

    uv run python sensitivity-artifact-checks/scripts/null_variant_calibration.py
"""

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats

from robustness_analysis import DROP_CODES, build_matrices, load_shards, log
from recover_vectors import panukb_url, parse_row

GRCH37_MB = {"1": 249, "2": 243, "3": 198, "4": 191, "5": 180, "6": 171, "7": 159,
             "8": 146, "9": 141, "10": 135, "11": 135, "12": 133, "13": 115,
             "14": 107, "15": 102, "16": 90, "17": 81, "18": 78, "19": 59,
             "20": 63, "21": 48, "22": 51, "X": 155}
EUR_NULL_MAX_NEGLOG10 = 4.0
ALPHA_LEVELS = (0.05, 0.01, 0.001)
BAND = (0.03, 0.07)

BASE = Path("sensitivity-artifact-checks/data/null_calibration")


def stage_a_positions():
    """One candidate position per 1-Mb window: first 30000-file variant in the
    window's opening 50 kb."""
    pos_file = BASE / "positions.json"
    if pos_file.exists():
        return json.load(open(pos_file))
    bed = BASE / "enum_regions.bed"
    with open(bed, "w") as f:
        for chrom, mb in GRCH37_MB.items():
            for w in range(1, mb):
                f.write(f"{chrom}\t{w * 1_000_000}\t{w * 1_000_000 + 50_000}\n")
    log("stage A: enumerating candidate positions (one tabix -R call)")
    result = subprocess.run(["tabix", "-R", str(bed), panukb_url("30000")],
                            capture_output=True, text=True, timeout=7200)
    if result.returncode != 0:
        raise RuntimeError(result.stderr[:300])
    chosen = {}
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 34:
            continue
        chrom, pos = fields[0], int(fields[1])
        win = (chrom, pos // 1_000_000)
        if win not in chosen:
            chosen[win] = {"chrom": chrom, "pos": pos, "ref": fields[2], "alt": fields[3]}
    positions = sorted(chosen.values(), key=lambda v: (v["chrom"], v["pos"]))
    json.dump(positions, open(pos_file, "w"))
    log(f"stage A: {len(positions)} candidate positions")
    return positions


def stage_b_shards(positions, trait_codes):
    bed = BASE / "sample_regions.bed"
    with open(bed, "w") as f:
        for v in positions:
            f.write(f"{v['chrom']}\t{v['pos'] - 1}\t{v['pos']}\n")
    targets = {}
    for v in positions:
        targets.setdefault((v["chrom"], v["pos"]), []).append(v)
    for i, code in enumerate(trait_codes, 1):
        shard = BASE / f"{code}.jsonl"
        done = BASE / f"{code}.done"
        if done.exists():
            continue
        log(f"stage B: trait {i}/{len(trait_codes)} ({code})")
        shard.unlink(missing_ok=True)
        result = subprocess.run(["tabix", "-R", str(bed), panukb_url(code)],
                                capture_output=True, text=True, timeout=7200)
        if result.returncode != 0:
            log(f"{code}: tabix failed: {result.stderr[:200]}")
            continue
        n = 0
        with open(shard, "w") as f:
            for line in result.stdout.splitlines():
                fields = line.split("\t")
                if len(fields) < 34:
                    continue
                key = (fields[0], int(fields[1]))
                for t in targets.get(key, []):
                    if fields[2] == t["ref"] and fields[3] == t["alt"]:
                        anc = parse_row(fields, len(fields) - 30)
                        vid = f"{t['chrom']}_{t['pos']}_{t['ref']}_{t['alt']}"
                        f.write(json.dumps({"locus_id": vid, "ancestries": anc}) + "\n")
                        n += 1
        done.write_text(json.dumps({"n": n}))
        log(f"{code}: {n} variants")


def main():
    BASE.mkdir(parents=True, exist_ok=True)
    corr = json.load(open("analysis/data/phenotypic_correlations_estimated.json"))
    codes = corr["trait_codes"]
    rho = np.array(corr["correlation_matrix"])
    nonred = [j for j, c in enumerate(codes) if c not in DROP_CODES]
    rho_nr = rho[np.ix_(nonred, nonred)]
    inv24, inv14 = np.linalg.inv(rho), np.linalg.inv(rho_nr)

    positions = stage_a_positions()
    stage_b_shards(positions, codes)

    data = load_shards(BASE, codes)
    log(f"stage C: {len(data)} variants with any data")

    rates = {a: [0, 0] for a in ALPHA_LEVELS}
    rates14 = {a: [0, 0] for a in ALPHA_LEVELS}
    n_tested = 0
    for vid, traits in data.items():
        eur_max = 0.0
        for c in codes:
            e = traits.get(c, {}).get("EUR")
            if e is None or e.get("neglog10_pval") is None:
                eur_max = 99.0
                break
            eur_max = max(eur_max, e["neglog10_pval"])
        if eur_max >= EUR_NULL_MAX_NEGLOG10:
            continue
        ancs = sorted(a for a in ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]
                      if all(a in traits.get(c, {}) for c in codes))
        built = build_matrices(traits, codes, ancs) if len(ancs) >= 3 else None
        if built is None:
            continue
        betas, ses = built
        K = len(ancs)
        n_tested += 1
        for r_inv, d, store in ((inv24, 24, rates), (inv14, 14, rates14)):
            b = betas[:, nonred] if d == 14 else betas
            s = ses[:, nonred] if d == 14 else ses
            W = [r_inv * np.outer(1 / s[k], 1 / s[k]) for k in range(K)]
            Wsi = np.linalg.inv(np.sum(W, axis=0))
            bbar = Wsi @ np.sum([W[k] @ b[k] for k in range(K)], axis=0)
            Q = float(sum((b[k] - bbar) @ W[k] @ (b[k] - bbar) for k in range(K)))
            p = float(stats.chi2.sf(Q, d * (K - 1)))
            for a in ALPHA_LEVELS:
                store[a][0] += p < a
                store[a][1] += 1

    out = {"metadata": {"timestamp": datetime.now(timezone.utc).isoformat(),
                        "amendment": "9fbed73f34cf"},
           "n_candidate": len(positions), "n_tested": n_tested}
    for name, store in (("panel24", rates), ("panel14", rates14)):
        out[name] = {str(a): round(c / max(n, 1), 5) for a, (c, n) in store.items()}
    r14 = out["panel14"]["0.05"]
    out["H_R9"] = {"reject_rate_14_at_0.05": r14,
                   "outcome": "holds" if BAND[0] <= r14 <= BAND[1] else "fails"}
    with open("sensitivity-artifact-checks/results/null_calibration.json", "w") as f:
        json.dump(out, f, indent=1)
    log(f"H_R9: n={n_tested} reject14@.05={r14} -> {out['H_R9']['outcome']} "
        f"(24-trait: {out['panel24']['0.05']}; tails 14: "
        f"{out['panel14']['0.01']}@.01 {out['panel14']['0.001']}@.001)")


if __name__ == "__main__":
    main()

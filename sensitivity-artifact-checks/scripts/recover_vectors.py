"""Recover per-locus per-ancestry effect vectors from Pan-UKB via remote tabix.

Step 1 of the frozen robustness registration (sensitivity-artifact-checks/PREREG.md).
Queries the 24-trait panel at the 2,567 v2 loci with one batched `tabix -R` call
per trait, writing one append-only JSONL shard per trait. A trait with a .done
marker is skipped on restart, so the job resumes cleanly.

Run from the repository root so cached .tbi index files are reused:

    uv run python sensitivity-artifact-checks/scripts/recover_vectors.py
"""

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

PANUKB_DATA = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files"
PANUKB_IDX = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files_tabix"
ANCESTRIES = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]
TABIX_TIMEOUT_S = 3600


def log(msg):
    print(f"[{datetime.now(timezone.utc).isoformat(timespec='seconds')}] {msg}", flush=True)


def panukb_url(pheno_code):
    data = f"{PANUKB_DATA}/continuous-{pheno_code}-both_sexes-irnt.tsv.bgz"
    idx = f"{PANUKB_IDX}/continuous-{pheno_code}-both_sexes-irnt.tsv.bgz.tbi"
    return f"{data}##idx##{idx}"


def load_loci(results_path):
    """Locus list from the stored v2 results: locus_id 'chr10_10037366_C_T'."""
    with open(results_path) as f:
        records = json.load(f)["results"]
    loci = {}
    for r in records:
        chrom, pos, ref, alt = r["locus_id"].rsplit("_", 3)
        chrom = chrom.removeprefix("chr")
        loci.setdefault((chrom, int(pos)), []).append(
            {"ref": ref, "alt": alt, "locus_id": r["locus_id"]}
        )
    return loci


CHROM_ORDER = {str(c): c for c in range(1, 23)} | {"X": 23, "Y": 24}


def write_regions_bed(loci, path):
    keys = sorted(loci, key=lambda k: (CHROM_ORDER.get(k[0], 99), k[1]))
    with open(path, "w") as f:
        for chrom, pos in keys:
            f.write(f"{chrom}\t{pos - 1}\t{pos}\n")


def parse_row(fields, anc_start):
    """Ancestry blocks are the final 30 columns: af, beta, se, neglog10p,
    low_confidence, each x 6 ancestries (AFR, AMR, CSA, EAS, EUR, MID)."""

    def safe_float(x):
        if x in ("NA", "", "true", "false"):
            return None
        try:
            return float(x)
        except ValueError:
            return None

    out = {}
    for i, anc in enumerate(ANCESTRIES):
        beta = safe_float(fields[anc_start + 6 + i])
        se = safe_float(fields[anc_start + 12 + i])
        if beta is None or se is None or se <= 0:
            continue
        out[anc] = {
            "af": safe_float(fields[anc_start + i]),
            "beta": beta,
            "se": se,
            "neglog10_pval": safe_float(fields[anc_start + 18 + i]),
            "low_confidence": fields[anc_start + 24 + i],
        }
    return out


def recover_trait(pheno_code, loci, bed_path, out_dir):
    shard = out_dir / f"{pheno_code}.jsonl"
    done = out_dir / f"{pheno_code}.done"
    if done.exists():
        log(f"{pheno_code}: done marker present, skipping")
        return

    shard.unlink(missing_ok=True)  # partial shard from an interrupted run
    t0 = time.time()
    result = subprocess.run(
        ["tabix", "-R", str(bed_path), panukb_url(pheno_code)],
        capture_output=True, text=True, timeout=TABIX_TIMEOUT_S,
    )
    if result.returncode != 0:
        log(f"{pheno_code}: TABIX FAILED: {result.stderr.strip()[:300]}")
        return

    n_rows = n_matched = 0
    with open(shard, "w") as f:
        for line in result.stdout.splitlines():
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 34:
                continue
            n_rows += 1
            key = (fields[0].removeprefix("chr"), int(fields[1]))
            for target in loci.get(key, []):
                if fields[2] == target["ref"] and fields[3] == target["alt"]:
                    anc = parse_row(fields, len(fields) - 30)
                    f.write(json.dumps(
                        {"locus_id": target["locus_id"], "ancestries": anc}
                    ) + "\n")
                    n_matched += 1

    done.write_text(json.dumps({
        "pheno_code": pheno_code,
        "n_rows_returned": n_rows,
        "n_loci_matched": n_matched,
        "elapsed_s": round(time.time() - t0, 1),
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }))
    log(f"{pheno_code}: {n_matched} loci matched ({n_rows} rows) in {time.time()-t0:.0f}s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path,
                    default=Path("analysis/results/genomewide_vector_q_v2_results.json"))
    ap.add_argument("--corr", type=Path,
                    default=Path("analysis/data/phenotypic_correlations_estimated.json"))
    ap.add_argument("--out", type=Path,
                    default=Path("sensitivity-artifact-checks/data/shards"))
    ap.add_argument("--traits", nargs="*", default=None,
                    help="Subset of phenocodes (default: the 24 panel codes from --corr)")
    ap.add_argument("--max-loci", type=int, default=None, help="Truncate locus list (testing)")
    args = ap.parse_args()

    with open(args.corr) as f:
        trait_codes = json.load(f)["trait_codes"]
    if args.traits:
        trait_codes = args.traits

    loci = load_loci(args.results)
    if args.max_loci:
        loci = dict(list(loci.items())[: args.max_loci])
    n_ids = sum(len(v) for v in loci.values())
    log(f"{len(loci)} positions / {n_ids} locus ids; {len(trait_codes)} traits")

    args.out.mkdir(parents=True, exist_ok=True)
    bed_path = args.out / "regions.bed"
    write_regions_bed(loci, bed_path)

    for i, code in enumerate(trait_codes, 1):
        log(f"--- trait {i}/{len(trait_codes)}: {code}")
        recover_trait(code, loci, bed_path, args.out)

    n_done = sum(1 for c in trait_codes if (args.out / f"{c}.done").exists())
    log(f"finished: {n_done}/{len(trait_codes)} trait shards complete")
    if n_done < len(trait_codes):
        sys.exit(1)


if __name__ == "__main__":
    main()

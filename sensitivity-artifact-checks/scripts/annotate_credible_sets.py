"""H_R6 credible-set annotation: BCX2 trans-ancestry CS95 membership.

A locus qualifies if its index variant is a member of a 95% credible set of
size <= 20 for at least one of its EUR-discovery traits, in the BCX2
trans-ancestry fine-mapping (GRCh37; alleles matched as an unordered pair,
since BCX2 MarkerNames list alleles lexicographically). Gate: >= 100
qualifying loci, else H_R6 is unresolvable.

    uv run python sensitivity-artifact-checks/scripts/annotate_credible_sets.py
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from robustness_analysis import GWS_NEGLOG10, load_shards, log

CS_MAX_SIZE = 20
GATE_MIN_LOCI = 100

# UKB field code -> BCX2 trait tag (traits without a BCX2 counterpart are
# unmappable; loci whose discovery traits are all unmappable cannot qualify).
BCX2_MAP = {
    "30000": "WBC", "30010": "RBC", "30020": "HGB", "30030": "HCT",
    "30040": "MCV", "30050": "MCH", "30060": "MCHC", "30070": "RDW",
    "30080": "PLT", "30100": "MPV", "30120": "LYM", "30130": "MON",
    "30140": "NEU", "30150": "EOS",
}


def load_cs95(bcx2_dir, tag):
    """(chrom, pos, frozenset(alleles)) -> smallest credible-set size."""
    path = bcx2_dir / f"{tag}_Trans_Credible_Sets_VEP.txt"
    out = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        i_marker = header.index("MarkerName")
        i_ncc = header.index("N.CC")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            marker = fields[i_marker]
            try:
                ncc = int(fields[i_ncc])
            except ValueError:
                continue
            chrom_pos, _, alleles = marker.partition("_")
            chrom, _, pos = chrom_pos.partition(":")
            a1, _, a2 = alleles.partition("_")
            key = (chrom, int(pos), frozenset((a1, a2)))
            out[key] = min(out.get(key, 10**9), ncc)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", type=Path, default=Path("sensitivity-artifact-checks/data/shards"))
    ap.add_argument("--corr", type=Path,
                    default=Path("analysis/data/phenotypic_correlations_estimated.json"))
    ap.add_argument("--v2", type=Path,
                    default=Path("analysis/results/genomewide_vector_q_v2_results.json"))
    ap.add_argument("--bcx2", type=Path, default=Path("sensitivity-artifact-checks/data/bcx2"))
    ap.add_argument("--out", type=Path,
                    default=Path("sensitivity-artifact-checks/results/credible_set_annotation.json"))
    args = ap.parse_args()

    with open(args.corr) as f:
        trait_codes = json.load(f)["trait_codes"]
    with open(args.v2) as f:
        v2_ids = [r["locus_id"] for r in json.load(f)["results"]]
    data = load_shards(args.shards, trait_codes)

    cs = {}
    for code, tag in BCX2_MAP.items():
        cs[code] = load_cs95(args.bcx2, tag)
        log(f"{tag}: {len(cs[code])} CS95 variants")

    qualifying = []
    n_no_mappable = 0
    for lid in v2_ids:
        locus = data.get(lid)
        if locus is None:
            continue
        chrom, pos, ref, alt = lid.rsplit("_", 3)
        key = (chrom.removeprefix("chr"), int(pos), frozenset((ref, alt)))
        disc = []
        for code in trait_codes:
            e = locus.get(code, {}).get("EUR")
            if e and e.get("neglog10_pval") is not None and e["neglog10_pval"] >= GWS_NEGLOG10:
                disc.append(code)
        mappable = [c for c in disc if c in BCX2_MAP]
        if not mappable:
            n_no_mappable += 1
            continue
        hits = [c for c in mappable if cs[c].get(key, 10**9) <= CS_MAX_SIZE]
        if hits:
            qualifying.append({"locus_id": lid, "traits": hits})

    gate_pass = len(qualifying) >= GATE_MIN_LOCI
    out = {
        "metadata": {"timestamp": datetime.now(timezone.utc).isoformat(),
                     "cs_max_size": CS_MAX_SIZE, "gate_min_loci": GATE_MIN_LOCI},
        "n_loci_checked": sum(1 for lid in v2_ids if lid in data),
        "n_no_mappable_discovery_trait": n_no_mappable,
        "n_qualifying": len(qualifying),
        "gate_pass": bool(gate_pass),
        "qualifying": qualifying,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    log(f"qualifying: {len(qualifying)} (gate >= {GATE_MIN_LOCI}: {gate_pass}); saved {args.out}")


if __name__ == "__main__":
    main()

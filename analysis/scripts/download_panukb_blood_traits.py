"""Download Pan-UKB per-ancestry blood cell trait data at target loci via tabix.

Requires htslib (tabix) installed locally.
Uses remote tabix queries — no large file downloads needed.

Usage:
    uv run python experiments/batch6_vector_extension/scripts/download_panukb_blood_traits.py
"""

import subprocess
import json
import csv
import sys
from pathlib import Path
from datetime import datetime, timezone

TABIX = "/usr/local/opt/htslib/bin/tabix"
BASE_URL = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files"

ANCESTRIES = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]

BLOOD_TRAITS = {
    "30000": "White blood cell count",
    "30010": "Red blood cell count",
    "30020": "Haemoglobin concentration",
    "30030": "Haematocrit percentage",
    "30040": "Mean corpuscular volume",
    "30050": "Mean corpuscular haemoglobin",
    "30060": "Mean corpuscular haemoglobin concentration",
    "30070": "Red blood cell distribution width",
    "30080": "Platelet count",
    "30090": "Platelet crit",
    "30100": "Mean platelet volume",
    "30120": "Lymphocyte count",
    "30130": "Monocyte count",
    "30140": "Neutrophil count",
    "30150": "Eosinophil count",
    "30180": "Lymphocyte percentage",
    "30190": "Monocyte percentage",
    "30200": "Neutrophil percentage",
    "30210": "Eosinophil percentage",
    "30220": "Basophil percentage",
    "30240": "Reticulocyte percentage",
    "30250": "Reticulocyte count",
    "30260": "Mean reticulocyte volume",
    "30270": "Mean sphered cell volume",
    "30280": "Immature reticulocyte fraction",
}

TARGET_LOCI = {
    "ACKR1_rs2814778": {"chr": "1", "pos": 159174683, "gene": "ACKR1/Duffy",
                         "note": "Duffy-null, malaria selection in AFR"},
    "SH2B3_rs3184504": {"chr": "12", "pos": 111884608, "gene": "SH2B3/ATXN2",
                         "note": "Multi-lineage blood trait locus"},
    "ABO_rs505922": {"chr": "9", "pos": 136131322, "gene": "ABO",
                      "note": "Blood group, platelet/VWF effects"},
    "JAK2_rs3780367": {"chr": "9", "pos": 5073770, "gene": "JAK2",
                        "note": "Myeloproliferative signaling"},
    "HFE_rs1800562": {"chr": "6", "pos": 26093141, "gene": "HFE",
                       "note": "Hereditary hemochromatosis, RBC traits"},
    "IL6R_rs2228145": {"chr": "1", "pos": 154426264, "gene": "IL6R",
                        "note": "Inflammation, WBC traits"},
    "IKZF1_rs4917014": {"chr": "7", "pos": 50344378, "gene": "IKZF1",
                          "note": "Lymphoid differentiation"},
    "TMPRSS6_rs855791": {"chr": "22", "pos": 37462936, "gene": "TMPRSS6",
                          "note": "Iron metabolism, RBC traits"},
    "MPL_rs6141": {"chr": "1", "pos": 43814027, "gene": "MPL",
                    "note": "Thrombopoietin receptor, platelet traits"},
    "GATA2_rs2335052": {"chr": "3", "pos": 128202534, "gene": "GATA2",
                         "note": "Hematopoietic transcription factor"},
    "CEBPA_rs2239633": {"chr": "19", "pos": 33792622, "gene": "CEBPA",
                         "note": "Myeloid differentiation"},
    "EPO_rs1617640": {"chr": "7", "pos": 100318423, "gene": "EPO",
                       "note": "Erythropoietin, RBC traits"},
    "CSF3R_rs3917932": {"chr": "1", "pos": 36932196, "gene": "CSF3R",
                         "note": "G-CSF receptor, neutrophil traits"},
    "GP1BA_rs2243093": {"chr": "17", "pos": 4836381, "gene": "GP1BA",
                         "note": "Platelet glycoprotein, platelet traits"},
    "ITGA2B_rs5911": {"chr": "17", "pos": 42449532, "gene": "ITGA2B",
                       "note": "Platelet integrin"},
    "RHD_rs590787": {"chr": "1", "pos": 25598832, "gene": "RHD",
                      "note": "Rh blood group, RBC traits"},
    "GATA1_rs11798027": {"chr": "X", "pos": 48649534, "gene": "GATA1",
                          "note": "Erythroid/megakaryocyte TF"},
    "HBS1L_MYB_rs9399137": {"chr": "6", "pos": 135418916, "gene": "HBS1L-MYB",
                              "note": "Multi-lineage locus, fetal hemoglobin"},
    "DARC_region_proxy": {"chr": "1", "pos": 159175354, "gene": "ACKR1 region",
                           "note": "Nearby variant in ACKR1 region for sensitivity"},
}


def parse_panukb_row(fields):
    """Parse a Pan-UKB TSV row into structured dict.

    Two layouts exist depending on whether meta_hq block is present:
    44 cols: 4 variant + 5 meta + 5 meta_hq + 6×5 ancestry
    39 cols: 4 variant + 5 meta + 6×5 ancestry
    """
    n_pops = 6
    n_cols = len(fields)
    n_meta = n_cols - 4 - (5 * n_pops)
    anc_start = 4 + n_meta

    def safe_float(x):
        if x == "NA" or x == "" or x in ("true", "false"):
            return None
        return float(x)

    result = {
        "chr": fields[0],
        "pos": int(fields[1]),
        "ref": fields[2],
        "alt": fields[3],
        "meta": {
            "af": safe_float(fields[4]),
            "beta": safe_float(fields[5]),
            "se": safe_float(fields[6]),
            "neglog10_pval": safe_float(fields[7]),
        },
    }

    for i, anc in enumerate(ANCESTRIES):
        af_idx = anc_start + i
        beta_idx = anc_start + n_pops + i
        se_idx = anc_start + 2 * n_pops + i
        pval_idx = anc_start + 3 * n_pops + i
        lc_idx = anc_start + 4 * n_pops + i
        result[anc] = {
            "af": safe_float(fields[af_idx]) if af_idx < n_cols else None,
            "beta": safe_float(fields[beta_idx]) if beta_idx < n_cols else None,
            "se": safe_float(fields[se_idx]) if se_idx < n_cols else None,
            "neglog10_pval": safe_float(fields[pval_idx]) if pval_idx < n_cols else None,
            "low_confidence": fields[lc_idx] if lc_idx < n_cols else None,
        }

    return result


def tabix_query(phenocode, region):
    """Query a Pan-UKB phenotype file for a genomic region."""
    filename = f"continuous-{phenocode}-both_sexes-irnt.tsv.bgz"
    url = f"{BASE_URL}/{filename}"
    cmd = [TABIX, url, region]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    if result.returncode != 0:
        print(f"  WARNING: tabix failed for {phenocode} at {region}: {result.stderr.strip()}")
        return []
    rows = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        rows.append(parse_panukb_row(fields))
    return rows


def download_locus_data(locus_name, locus_info, traits):
    """Download all blood trait data at a single locus."""
    chrom = locus_info["chr"]
    pos = locus_info["pos"]
    region = f"{chrom}:{pos}-{pos}"

    print(f"\n{'='*60}")
    print(f"Locus: {locus_name} ({locus_info['gene']}, {chrom}:{pos})")
    print(f"  {locus_info['note']}")
    print(f"{'='*60}")

    locus_data = {"locus": locus_name, "info": locus_info, "traits": {}}

    for phenocode, trait_name in traits.items():
        rows = tabix_query(phenocode, region)
        if not rows:
            window = f"{chrom}:{pos-500}-{pos+500}"
            rows = tabix_query(phenocode, window)
            if rows:
                closest = min(rows, key=lambda r: abs(r["pos"] - pos))
                print(f"  {trait_name}: exact position not found, using closest at {closest['pos']} (delta={abs(closest['pos']-pos)}bp)")
                rows = [closest]
            else:
                print(f"  {trait_name}: NO DATA in +/-500bp window")
                continue

        row = rows[0]
        locus_data["traits"][phenocode] = {
            "name": trait_name,
            "variant": f"{row['chr']}:{row['pos']}:{row['ref']}:{row['alt']}",
            "meta_neglog10_pval": row["meta"]["neglog10_pval"],
            "ancestries": {},
        }

        for anc in ANCESTRIES:
            anc_data = row[anc]
            if anc_data["beta"] is not None:
                locus_data["traits"][phenocode]["ancestries"][anc] = {
                    "af": anc_data["af"],
                    "beta": anc_data["beta"],
                    "se": anc_data["se"],
                    "neglog10_pval": anc_data["neglog10_pval"],
                    "low_confidence": anc_data["low_confidence"],
                }

        n_anc = len(locus_data["traits"][phenocode]["ancestries"])
        meta_p = row["meta"]["neglog10_pval"]
        print(f"  {trait_name}: {n_anc} ancestries, meta -log10(p) = {meta_p:.1f}" if meta_p else f"  {trait_name}: {n_anc} ancestries")

    return locus_data


def main():
    out_dir = Path(__file__).resolve().parent.parent / "data"
    out_dir.mkdir(exist_ok=True)

    all_results = {
        "download_timestamp": datetime.now(timezone.utc).isoformat(),
        "source": "Pan-UK Biobank (pan.ukbb.broadinstitute.org)",
        "genome_build": "GRCh37",
        "n_traits": len(BLOOD_TRAITS),
        "n_loci": len(TARGET_LOCI),
        "ancestries": ANCESTRIES,
        "loci": {},
    }

    for locus_name, locus_info in TARGET_LOCI.items():
        locus_data = download_locus_data(locus_name, locus_info, BLOOD_TRAITS)
        all_results["loci"][locus_name] = locus_data

    out_file = out_dir / "panukb_blood_traits_at_loci.json"
    with open(out_file, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved to {out_file}")

    write_flat_tsv(all_results, out_dir / "panukb_blood_traits_flat.tsv")
    print(f"Flat TSV: {out_dir / 'panukb_blood_traits_flat.tsv'}")


def write_flat_tsv(results, out_path):
    """Write a flat TSV for easy inspection."""
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "locus", "gene", "phenocode", "trait", "variant",
            "ancestry", "af", "beta", "se", "neglog10_pval", "low_confidence",
        ])
        for locus_name, locus_data in results["loci"].items():
            gene = locus_data["info"]["gene"]
            for phenocode, trait_data in locus_data["traits"].items():
                for anc, anc_data in trait_data["ancestries"].items():
                    writer.writerow([
                        locus_name, gene, phenocode, trait_data["name"],
                        trait_data["variant"], anc,
                        anc_data["af"], anc_data["beta"], anc_data["se"],
                        anc_data["neglog10_pval"], anc_data["low_confidence"],
                    ])


if __name__ == "__main__":
    main()

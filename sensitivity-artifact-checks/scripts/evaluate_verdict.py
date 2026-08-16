"""Evaluate H_R1-H_R7, apply the registered verdict table and M* band.

Final step of the frozen registration (parent 4955185ed63b, amendment 10fbfba26adc).
Reads the results of every prior stage; computes nothing new except H_R6's
subset M (BH within the credible-set-qualifying loci on the non-redundant
statistic, per the registered wording).

    uv run python sensitivity-artifact-checks/scripts/evaluate_verdict.py
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from robustness_analysis import ALPHA, bh_significant, log

H3_MAX_SHIFT = 10.0  # percentage points
H4_GATE = 50
DOMINANCE = 0.50


def main():
    ap = argparse.ArgumentParser()
    base = Path("sensitivity-artifact-checks/results")
    ap.add_argument("--robustness", type=Path, default=base / "robustness_results.json")
    ap.add_argument("--chr1", type=Path, default=base / "robustness_chr1.json")
    ap.add_argument("--afr", type=Path, default=base / "robustness_afr.json")
    ap.add_argument("--credible", type=Path, default=base / "credible_set_annotation.json")
    ap.add_argument("--bootstrap", type=Path, default=base / "bootstrap_results.json")
    ap.add_argument("--out", type=Path, default=base / "verdict.json")
    args = ap.parse_args()

    rob = json.load(open(args.robustness))
    boot = json.load(open(args.bootstrap))
    cred = json.load(open(args.credible))
    S = rob["variant_summaries"]
    g2 = boot["G2"]

    def variant_ok(name):
        v = g2.get(name, {})
        return v.get("pass") is True

    hypotheses = {}

    # H_R1 / H_R2 / H_R4 / H_R5: dominance on their variants, void if G2 voids the variant
    for hid, variant in [("H_R1", "nonredundant"), ("H_R2", "rank_truncated"),
                         ("H_R4", "non_eur"), ("H_R5", "no_discovery")]:
        s = S[variant]
        entry = {"variant": variant, "M": s["M"],
                 "n_significant": s["n_fdr_significant"],
                 "g2_calibrated": variant_ok(variant)}
        if not variant_ok(variant):
            entry["outcome"] = "unresolvable (G2 void)"
        elif hid == "H_R4" and s["n_fdr_significant"] < H4_GATE:
            entry["outcome"] = "unresolvable (gate: <50 significant loci)"
        else:
            entry["outcome"] = "holds" if (s["M"] or 0) > DOMINANCE else "fails"
        hypotheses[hid] = entry

    # H_R3: both alternative-matrix arms shift full-panel M by < 10 points
    m_base = S["full"]["M"] * 100
    arms = {}
    for tag, path in [("chr1_EUR", args.chr1), ("chr22_AFR", args.afr)]:
        m_alt = json.load(open(path))["variant_summaries"]["full"]["M"] * 100
        arms[tag] = {"M": round(m_alt, 2), "shift_points": round(abs(m_alt - m_base), 2)}
    h3_ok = all(a["shift_points"] < H3_MAX_SHIFT for a in arms.values())
    hypotheses["H_R3"] = {"baseline_M_pct": round(m_base, 2), "arms": arms,
                          "g2_calibrated": variant_ok("full"),
                          "outcome": ("unresolvable (G2 void)" if not variant_ok("full")
                                      else "holds" if h3_ok else "fails")}

    # H_R6: BH within qualifying subset on the non-redundant statistic
    qual = {q["locus_id"] for q in cred["qualifying"]}
    nr = [r for r in rob["variant_records"]["nonredundant"] if r["locus_id"] in qual]
    if not cred["gate_pass"]:
        hypotheses["H_R6"] = {"outcome": "unresolvable (coverage gate)",
                              "n_qualifying": cred["n_qualifying"]}
    else:
        sig = bh_significant([r["p"] for r in nr])
        n_sig = int(sig.sum())
        n_mo = sum(1 for r, s in zip(nr, sig) if s and r["category"] == "multivariate_only")
        m6 = n_mo / n_sig if n_sig else None
        hypotheses["H_R6"] = {"n_qualifying": len(nr), "n_significant": n_sig,
                              "M": round(m6, 4) if m6 is not None else None,
                              "g2_calibrated": variant_ok("nonredundant"),
                              "outcome": ("unresolvable (G2 void)" if not variant_ok("nonredundant")
                                          else "holds" if (m6 or 0) > DOMINANCE else "fails")}

    # H_R7: from bootstrap
    h7 = boot["H_R7"]
    hypotheses["H_R7"] = {**h7, "outcome": "holds" if h7["pass"] else "fails"}

    # --- verdict table (parent registration, carried by amendment 1)
    o = {h: hypotheses[h]["outcome"] for h in hypotheses}
    holds = {h: o[h] == "holds" for h in o}
    unres = {h: o[h].startswith("unresolvable") for h in o}

    def hold_or_unres(h):
        return holds[h] or unres[h]

    channel_fails = sum(1 for h in ["H_R3", "H_R4", "H_R5", "H_R6"] if o[h] == "fails")
    if not (holds["H_R1"] or unres["H_R1"]) or not (holds["H_R2"] or unres["H_R2"]) \
            or channel_fails >= 2 or o["H_R1"] == "fails" or o["H_R2"] == "fails":
        verdict = "DOES NOT SURVIVE"
    elif holds["H_R1"] and holds["H_R2"] and holds["H_R5"] and holds["H_R3"] \
            and hold_or_unres("H_R4") and hold_or_unres("H_R6"):
        verdict = "SURVIVES, STRENGTHENED" if holds["H_R7"] else "SURVIVES"
    elif holds["H_R1"] and unres["H_R2"] and holds["H_R5"] and holds["H_R3"] \
            and hold_or_unres("H_R4") and hold_or_unres("H_R6"):
        # Amendment 2 (frozen): unresolvable H_R2 with contrary exploratory
        # recalibration assigns SURVIVES WITH CAVEAT (correlation-structure dependence).
        verdict = "SURVIVES WITH CAVEAT"
    elif holds["H_R1"] and holds["H_R2"] and (
            channel_fails == 1 or o["H_R4"] == "fails"):
        verdict = "SURVIVES WITH CAVEAT"
    else:
        verdict = "NOT CLASSIFIED BY TABLE (registration gap)"

    m_star = S["nonredundant"]["M"]
    if not variant_ok("nonredundant") and variant_ok("rank_truncated"):
        m_star = S["rank_truncated"]["M"]  # G2 fallback per parent registration
    band = None
    if m_star is not None and "SURVIVES" in verdict:
        band = "minimal correction" if m_star >= 0.75 else \
               "substantial correction" if m_star >= 0.50 else "below bands"

    out = {
        "metadata": {"timestamp": datetime.now(timezone.utc).isoformat(),
                     "parent": "4955185ed63b", "amendment": "10fbfba26adc"},
        "hypotheses": hypotheses,
        "verdict": verdict,
        "M_star": m_star,
        "band": band,
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    log(json.dumps({h: o[h] for h in sorted(o)}, indent=1))
    log(f"VERDICT: {verdict}  M*={m_star}  band={band}")
    log(f"saved {args.out}")


if __name__ == "__main__":
    main()

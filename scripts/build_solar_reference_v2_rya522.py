#!/usr/bin/env python3
"""
scripts/build_solar_reference_v2_rya522.py  (RYA-522)
=====================================================
Regenerate the solar gold reference v2 candidate ENTIRELY from the phase_c verdict
channel (RYA-521 single authoritative source) — NOT by hand-editing v1's C cell.
Produces two artifacts (candidates; the write-once freeze is a separate, Ryan-
ratified step via promote_solar_reference):

  data/reference/solar/solar_abundances_v2_candidate.csv   — verdict-sourced rows
  docs/audit/solar_gold_v2_ratification_rya522.md          — the diff table (step-4
      ratification artifact): v2 | scale | v1 | Δ(v2-v1) | Asplund2021 | Δ(v2-Asplund)

Asplund column = the in-repo cited reference SOLAR_ASPLUND2021 (Asplund, Amarsi &
Grevesse 2021, A&A 653 A141), never typed from memory. Scale is stated per row so a
documented 1D-NLTE-vs-3D offset is not misread as a discrepancy.

Usage: python scripts/build_solar_reference_v2_rya522.py --verdict <path>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import SOLAR_ASPLUND2021, TARGET_ELEMENTS  # noqa: E402

ASPLUND_CITE = "Asplund, Amarsi & Grevesse 2021, A&A 653, A141"


def _scale_and_note(el: str, ch: str, verdict: str, a, asp):
    ch = ch or ""
    d = (a - asp) if (a is not None and asp is not None) else None
    if el == "Fe":
        return "1D-NLTE (Fe I)", "our 1D-NLTE runs ~+0.05 above Asplund 3D-true 7.46 (RYA-336) — documented offset, NOT a discrepancy"
    if "synthesis" in ch:
        base = "synthesis"
        return base, ("CH G-band + C I (RYA-237)" if el == "C" else
                      "O I 777 + [O I] 6300 (RYA-237)" if el == "O" else
                      "HFS-resolved synthesis (RYA-411/466/473)")
    if "kittpeak" in ch:
        n = {"N": "N I red multiplets; +0.37 = owed NLTE (RYA-369)",
             "P": "near-IR multiplet, gf-limited (RYA-460)",
             "K": "K I 7699 + K NLTE grid (RYA-462)",
             "Co": "blue-edge, SNR-limited — value not fully trusted (RYA-460)",
             "Sc": "blue-edge HFS single line (RYA-460)"}.get(el, "Kitt Peak atlas (RYA-460)")
        return "atlas 1D", n
    # EW-curation path
    if el == "Li":
        return "EW (upper limit)", "CN-blended, carried as UPPER LIMIT (RYA-103) — a clean low value is a red flag"
    if verdict == "CURATION-OWED" and d is not None and d > 0.15:
        return "EW 1D-LTE/NLTE", f"gf-limited residual floor ({d:+.2f}); NOT an Asplund disagreement (RYA-399)"
    return "EW 1D-LTE/NLTE", "curated EW; low-confidence" if verdict == "CURATION-OWED" else ""


def load_v1(path):
    df = pd.read_csv(path, comment="#")
    out = {}
    for _, r in df.iterrows():
        el = str(r["element"])
        a = r.get("A_X_nlte")
        a = a if pd.notna(a) else r.get("A_X")
        if el not in out and pd.notna(a):
            out[el] = float(a)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verdict", required=True, help="phase_c verdict json (authoritative source)")
    ap.add_argument("--v1", default=str(ROOT / "data/reference/solar/solar_abundances_v1.csv"))
    args = ap.parse_args()

    V = {r["element"]: r for r in json.loads(Path(args.verdict).read_text())["verdicts"]}
    v1 = load_v1(args.v1)

    cand_rows, table = [], []
    for el in TARGET_ELEMENTS:
        r = V.get(el, {})
        a = r.get("A_measured")
        a = float(a) if a is not None else None
        ch, verdict = r.get("channel", ""), r.get("verdict", "")
        asp = SOLAR_ASPLUND2021.get(el)
        scale, note = _scale_and_note(el, ch, verdict, a, asp)
        v1v = v1.get(el)
        cand_rows.append({"element": el, "A_authoritative": a, "verdict": verdict,
                          "channel": ch, "scale": scale, "n_lines": r.get("n_lines"),
                          "source": "phase_c_verdict (RYA-521)"})
        table.append({
            "Element": el,
            "v2 (verdict)": f"{a:.3f}" if a is not None else "owed",
            "method/scale": scale,
            "v1 (old)": f"{v1v:.3f}" if v1v is not None else "—",
            "Δ(v2−v1)": f"{a - v1v:+.3f}" if (a is not None and v1v is not None) else "—",
            f"Asplund 2021": f"{asp:.2f}" if asp is not None else "—",
            "Δ(v2−Asp)": f"{a - asp:+.3f}" if (a is not None and asp is not None) else "—",
            "verdict": verdict, "note": note,
        })

    cand = pd.DataFrame(cand_rows)
    cand_path = ROOT / "data/reference/solar/solar_abundances_v2_candidate.csv"
    cand.to_csv(cand_path, index=False)

    tdf = pd.DataFrame(table)
    cols = list(tdf.columns)
    md_table = ["| " + " | ".join(cols) + " |",
                "|" + "|".join("---" for _ in cols) + "|"]
    md_table += ["| " + " | ".join(str(r[c]) for c in cols) + " |" for _, r in tdf.iterrows()]
    md = ["# Solar gold reference v2 — ratification diff table (RYA-522)", "",
          f"Source: the phase_c **verdict** channel (RYA-521), regenerated — not hand-edited. "
          f"Asplund column = `SOLAR_ASPLUND2021` ({ASPLUND_CITE}).", "",
          "**Scales differ:** our values are 1D-NLTE / synthesis on our stack; Asplund 2021 is "
          "3D-NLTE photospheric. The `note` states each row's scale so documented offsets "
          "(e.g. Fe I +0.05, RYA-336) are not misread as disagreement.", "",
          *md_table, "",
          f"**Headline:** C {v1.get('C')} → {float(V['C']['A_measured']):.3f} "
          f"(the RYA-520 saturated-C I-5380 artifact corrected; −1.77 dex).",
          "", "Freeze is GATED on Ryan's explicit ratification of this table (RYA-522 step 4)."]
    doc = ROOT / "docs/audit/solar_gold_v2_ratification_rya522.md"
    doc.parent.mkdir(parents=True, exist_ok=True)
    doc.write_text("\n".join(md) + "\n")
    print(f"wrote {cand_path.relative_to(ROOT)} ({len(cand)} elements)")
    print(f"wrote {doc.relative_to(ROOT)}")
    print(f"C: {v1.get('C')} -> {float(V['C']['A_measured']):.3f}")


if __name__ == "__main__":
    main()

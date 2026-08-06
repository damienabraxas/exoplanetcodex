#!/usr/bin/env python3
"""
scripts/build_solar_reference_v2_rya522.py  (RYA-522)
=====================================================
Regenerate the solar gold reference v2 candidate ENTIRELY from the phase_c verdict
channel (RYA-521 single authoritative source) — NOT by hand-editing v1's C cell —
and TIER it by row-confidence per Ryan's ratification (RYA-522 comment 2026-07-05).

Tiers (ratified): a value is frozen ONLY if we would stake a differential on it.
  gold        — C, O, K, Mn, Fe, Sc            (normal method/scale spread of Asplund)
  gf_floor    — Cr, Si                          (characterized ~+0.4 floor; Si<=0.5)
  upper_limit — Li                              (RYA-103; tagged, not a point value)
  owed        — everything else                 (NO frozen value; suspect -> held, not
                                                  immortalised — the C=10.26 lesson)
S is owed (offset +0.63 > ~0.5). Sr is the +2.1-tail saturation suspect (owed).

Artifacts (candidates; the write-once freeze is promote_solar_reference --apply):
  data/reference/solar/solar_abundances_v2_candidate.csv
  docs/audit/solar_gold_v2_ratification_rya522.md
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import SOLAR_ASPLUND2021, TARGET_ELEMENTS  # noqa: E402
from pipeline.provenance_honesty import (  # noqa: E402  RYA-653 shared tripwire
    assert_blank_cause_is_honest)

ASPLUND_CITE = "Asplund, Amarsi & Grevesse 2021, A&A 653, A141"
HONESTY_SITE = "gold reference builder (RYA-522)"

# Ratified tiers (RYA-522, Ryan 2026-07-05). Everything not listed -> owed.
TIER_GOLD = {"C", "O", "K", "Mn", "Fe", "Sc"}
TIER_GF_FLOOR = {"Cr", "Si"}
TIER_UPPER_LIMIT = {"Li"}
# Dominant ion for the frozen rows (cosmetic except Fe I, which the anchor check needs).
ION = {"Sc": "II"}


def confidence_of(el: str) -> str:
    if el in TIER_GOLD:
        return "gold"
    if el in TIER_GF_FLOOR:
        return "gf_floor"
    if el in TIER_UPPER_LIMIT:
        return "upper_limit"
    return "owed"


# Per-element synthesis notes. These are the elements whose synthesis channel this
# builder was written against; anything else carries the VERDICT's own provenance
# rather than one of these citations (RYA-653 — the fallthrough used to stamp
# "HFS-resolved synthesis (RYA-411/466/473)" on any synthesis row, which is a
# citation this stage never established: it is wrong for Ba, whose synthesis is
# RYA-559).
SYNTHESIS_NOTES = {
    "C": "CH G-band + C I (RYA-237)",
    "O": "O I 777 + [O I] 6300 (RYA-237)",
    "Mn": "HFS-resolved synthesis (RYA-411/466/473)",
    "Cu": "HFS-resolved synthesis (RYA-411/466/473)",
    "V": "HFS-resolved synthesis (RYA-411/466/473)",
}
KITTPEAK_NOTES = {"N": "N I red multiplets; +0.37 owed NLTE (RYA-369)",
                  "P": "near-IR multiplet, gf-limited (RYA-460)",
                  "K": "K I 7699 + K NLTE grid (RYA-462)",
                  "Co": "blue-edge, SNR-limited — not trusted (RYA-460)",
                  "Sc": "blue-edge HFS single line (RYA-460)"}


def _scale_and_note(el, ch, verdict, a, asp, conf, provenance=""):
    """Return (method_scale, note) for one row.

    RYA-653: every branch here must state a cause this builder can stand behind.
    Where it has none of its own it carries the verdict channel's — attributed,
    not invented. The old tail fabricated "no independent-gf line survives the
    graded cull" for EVERY value-less row, which is how gold v2's Ba row came to
    blame a cull for an element RYA-559 has measured.
    """
    ch = ch or ""
    d = (a - asp) if (a is not None and asp is not None) else None
    if el == "Fe":
        return "1D-NLTE (Fe I)", "our 1D-NLTE runs ~+0.05 above Asplund 3D-true 7.46 (RYA-336) — documented offset, NOT a discrepancy"
    if el == "Li":
        return "EW (upper limit)", "CN-blended, carried as UPPER LIMIT (RYA-103) — a clean low value is a red flag"
    if "synthesis" in ch:
        return "synthesis", (SYNTHESIS_NOTES.get(el) or provenance or ch)
    if "kittpeak" in ch:
        return "atlas 1D", KITTPEAK_NOTES.get(el, "Kitt Peak atlas (RYA-460)")
    if conf == "gf_floor" and d is not None:
        return "EW 1D-LTE/NLTE", f"characterized gf-scale floor ({d:+.2f}); 3D not the lever (RYA-398/399)"
    if el == "Sr" and d is not None:
        return "EW (SUSPECT)", f"+{d:.2f} — NOT a gf-floor; saturated-line-on-flat-COG signature (the RYA-520 disease) → saturation-trace owed"
    if d is not None:
        return "EW 1D-LTE/NLTE", f"LOW_CONFIDENCE / thin graded pool ({d:+.2f})"
    # No delta => this stage has established NO cause of its own. Carry the
    # verdict's stated cause verbatim (the verdict's own RYA-596 tripwire has
    # already vetted it), and never fabricate one.
    return "EW 1D-LTE/NLTE", (ch or provenance or
                              "no cause stated by the verdict channel")


def _rel(path):
    p = Path(path).resolve()
    try:
        return p.relative_to(ROOT)
    except ValueError:
        return p


def _headline(v1, V):
    """C's re-derivation, stated from the two artifacts — never from memory
    (the old text carried a hardcoded '−1.77 dex' that no rebuild could keep true)."""
    c_old, c_new = v1.get("C"), V.get("C", {}).get("A_measured")
    if c_old is None or c_new is None:
        return "**Headline:** C not present in both the v1 reference and the verdict channel."
    return (f"**Headline:** C {c_old:.3f} → {float(c_new):.3f} "
            f"(RYA-520 saturated-C I-5380 fix, {float(c_new) - c_old:+.2f} dex).")


def _sr_line(v1, V):
    """Sr's standing, read off the verdict channel rather than frozen in prose."""
    sr = V.get("Sr", {})
    a, asp = sr.get("A_measured"), SOLAR_ASPLUND2021.get("Sr")
    if a is not None and asp is not None:
        return (f"**+2.1-tail suspect:** Sr (Δ{float(a) - asp:+.2f}) — "
                "saturated-line-on-flat-COG signature (RYA-520 class) → held owed, "
                "routed to a saturation-trace ticket.")
    return (f"**Sr:** no value in the verdict channel — {sr.get('channel', 'element absent')}. "
            "The +2.1 saturation-trace remains owed (RYA-520 class).")


def diff_vs_frozen(cand: pd.DataFrame, frozen: pd.DataFrame, frozen_version: str):
    """Per-cell diff of a rebuilt candidate against the live FROZEN reference.

    This is the promote-time question stated honestly: "what changes if this
    candidate is frozen?" — derived from the two artifacts, never narrated.
    """
    o = frozen.set_index("element")
    lines = [f"| element | field | frozen {frozen_version} | rebuilt candidate |",
             "|---|---|---|---|"]
    n_changed = 0
    for _, row in cand.iterrows():
        el = row["element"]
        for col in cand.columns:
            if col == "element":
                continue
            old = o.loc[el, col] if el in o.index and col in o.columns else "<absent>"
            new = row[col]
            if _cell(old) == _cell(new):
                continue
            n_changed += 1
            lines.append(f"| {el} | {col} | {_cell(old)} | {_cell(new)} |")
    return lines, n_changed


def _cell(v):
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return "—"
    return str(v).replace("|", "\\|")


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


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--verdict", required=True)
    ap.add_argument("--v1", default=str(ROOT / "data/reference/solar/solar_abundances_v1.csv"))
    # Output paths are overridable so the honesty tests can build against a
    # scratch dir without touching the real candidate artifact (RYA-653).
    ap.add_argument("--out-csv", default=str(ROOT / "data/reference/solar/solar_abundances_v2_candidate.csv"))
    ap.add_argument("--out-md", default=str(ROOT / "docs/audit/solar_gold_v2_ratification_rya522.md"))
    ap.add_argument("--diff-md", default=None,
                    help="also write a per-cell diff of the rebuilt candidate vs the "
                         "live FROZEN reference (RYA-653 deliverable; freezing stays "
                         "a separate, ratified act)")
    args = ap.parse_args(argv)

    V = {r["element"]: r for r in json.loads(Path(args.verdict).read_text())["verdicts"]}
    v1 = load_v1(args.v1)

    cand_rows, table = [], []
    for el in TARGET_ELEMENTS:
        r = V.get(el, {})
        a_verdict = r.get("A_measured")
        a_verdict = float(a_verdict) if a_verdict is not None else None
        conf = confidence_of(el)
        # OWED tier freezes NO authoritative value (held), even if the verdict has one.
        a_frozen = a_verdict if conf != "owed" else None
        ch, verdict = r.get("channel", ""), r.get("verdict", "")
        asp = SOLAR_ASPLUND2021.get(el)
        scale, note = _scale_and_note(el, ch, verdict, a_verdict, asp, conf,
                                      provenance=r.get("provenance", ""))
        # An `owed` row that HAS a verdict value is held, not absent — say so on
        # the row, so the withheld value is visible instead of the row reading as
        # "nothing was ever measured" (RYA-653; the Ba/RYA-559 shape).
        if conf == "owed" and a_verdict is not None:
            note = (f"{note} — A(X) {a_verdict:.3f} HELD at tier 'owed' (RYA-522), "
                    f"not frozen")
        n_lines = r.get("n_lines")
        # RYA-653 tripwire (shared with the phase_c verdict): neither the cause we
        # write nor the cause we inherit may claim a zero-survivor cull on a row
        # that carries survivors or a measured value.
        for claim in (ch, note):
            assert_blank_cause_is_honest(el, claim, n_lines, a_measured=a_verdict,
                                         site=HONESTY_SITE)
        v1v = v1.get(el)
        cand_rows.append({
            "element": el, "ion": ION.get(el, "I"),
            "A_X": a_frozen if a_frozen is not None else np.nan,
            "A_X_nlte": a_frozen if a_frozen is not None else np.nan,
            "confidence": conf, "verdict": verdict, "method_scale": scale,
            "asplund2021": asp, "n_lines": n_lines,
            "source": "phase_c_verdict (RYA-521)", "note": note,
        })
        val = (f"{a_frozen:.3f}" if a_frozen is not None
               else (f"[{a_verdict:.3f} held]" if a_verdict is not None else "owed"))
        table.append({
            "Element": el, "conf": conf, "v2 (frozen)": val, "method/scale": scale,
            "v1 (old)": f"{v1v:.3f}" if v1v is not None else "—",
            "Δ(v2−v1)": f"{a_frozen - v1v:+.3f}" if (a_frozen is not None and v1v is not None) else "—",
            "Asplund 2021": f"{asp:.2f}" if asp is not None else "—",
            "Δ(v2−Asp)": f"{a_verdict - asp:+.3f}" if (a_verdict is not None and asp is not None) else "—",
            "note": note,
        })

    cand = pd.DataFrame(cand_rows)
    cand_path = Path(args.out_csv)
    cand_path.parent.mkdir(parents=True, exist_ok=True)
    cand.to_csv(cand_path, index=False)

    cols = list(table[0].keys())
    mdt = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    mdt += ["| " + " | ".join(str(row[c]) for c in cols) + " |" for row in table]
    n = {t: sum(1 for e in TARGET_ELEMENTS if confidence_of(e) == t)
         for t in ("gold", "gf_floor", "upper_limit", "owed")}
    md = ["# Solar gold reference v2 — TIERED ratification table (RYA-522)", "",
          f"Verdict-sourced (RYA-521), tiered by row-confidence per Ryan's ratification. "
          f"Asplund = `SOLAR_ASPLUND2021` ({ASPLUND_CITE}).", "",
          f"**Tiers:** gold={n['gold']} · gf_floor={n['gf_floor']} · upper_limit={n['upper_limit']} "
          f"· owed(held, no frozen value)={n['owed']}.", "",
          "Scales differ (ours 1D-NLTE/synth vs Asplund 3D-NLTE) — `note` states each row's scale so "
          "documented offsets are not misread as disagreement. `owed` rows freeze NO value (the C=10.26 "
          "lesson: suspect → held, not immortalised).", "",
          *mdt, "",
          _headline(v1, V), _sr_line(v1, V), ""]
    doc = Path(args.out_md)
    doc.parent.mkdir(parents=True, exist_ok=True)
    doc.write_text("\n".join(md) + "\n")
    print(f"wrote {_rel(cand_path)} ({len(cand)} rows; "
          f"gold={n['gold']} gf_floor={n['gf_floor']} upper_limit={n['upper_limit']} owed={n['owed']})")
    print(f"wrote {_rel(doc)}")

    if args.diff_md:
        from pipeline.data_namespace import read_solar_reference
        frozen, version = read_solar_reference()
        rows, n_changed = diff_vs_frozen(cand, frozen, version)
        dpath = Path(args.diff_md)
        dpath.parent.mkdir(parents=True, exist_ok=True)
        dpath.write_text("\n".join([
            f"# Rebuilt solar gold candidate vs frozen {version} — per-cell diff", "",
            f"Candidate: `{_rel(cand_path)}` (rebuilt from `{_rel(args.verdict)}`).",
            f"Frozen reference: `{version}` via `pipeline.data_namespace.read_solar_reference()` "
            "(the CURRENT pointer). **This diff is a proposal, not a freeze** — promoting it is "
            "`scripts/promote_solar_reference.py --apply` and a ratification decision (RYA-527).", "",
            f"**{n_changed} cell(s) differ.**", "", *rows, ""]) + "\n")
        print(f"wrote {_rel(dpath)} ({n_changed} cells differ vs frozen {version})")


if __name__ == "__main__":
    main()

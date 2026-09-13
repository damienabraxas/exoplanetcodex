#!/usr/bin/env python3
"""RYA-1171 — does any stored gf grade claim BETTER precision than the source it cites?

RYA-1160's EP-aware control validated the 5 pre-existing CNO graded rows against NIST ASD:
every log gf reproduced to <= 0.001 dex, but the store graded the O I 777 triplet **A+
(<=2%)** where NIST publishes **A (<=3%)**. The values were right and the claimed precision
was not -- so every gf_sigma_dex derived from that grade was too small, by a factor of
log10(1.03)/log10(1.02) = 1.49, on the most-used oxygen diagnostic in the programme.

This sweep re-runs that comparison and widens it: every grade the store holds is checked
against every INDEPENDENT record of NIST's published grade that we hold, and any row
claiming a tighter accuracy class than its source is reported.

⚠️ IT NEVER INVENTS A GRADE. Where RYA-1160 refused a match -- [O I] 6300.304 has two
EP-aware candidates -- there is no NIST grade to compare against, and the row is reported
as UNRESOLVED rather than assigned one. A refusal is the correct answer to an ambiguity
(RYA-1072), and silently replacing it with a guess is the defect, not the fix.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.gf_grades import NIST_ACC_PCT, nist_sigma_dex  # noqa: E402

CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"
CONTROL = ROOT / "data" / "audit" / "rya1160_cno_nist_gf" / "control.csv"
#: Independent records of a published NIST grade that we hold locally.
NIST_RECORDS = (
    ROOT / "data" / "linelists" / "nist_reference.csv",
    ROOT / "data" / "linelists" / "nist_crosscheck.csv",
)

#: Wavelength/EP tolerances of the EP-aware match (RYA-1037's dual key).
TOL_A, TOL_EP = 0.005, 0.002


def _rows(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


def tighter(a: str, b: str) -> bool:
    """Is grade `a` a TIGHTER accuracy claim than `b`? Unknown classes never compare."""
    pa, pb = NIST_ACC_PCT.get(a.strip()), NIST_ACC_PCT.get(b.strip())
    return pa is not None and pb is not None and pa < pb


def store_grades() -> list[dict]:
    return [r for r in _rows(CANONICAL) if (r["nist_grade"] or "").strip()]


def sweep() -> tuple[list[str], list[str]]:
    """(over-claims, notes). An over-claim is a hard finding; a note is reported only."""
    over, suspect, notes = [], [], []
    store = store_grades()
    by_key = {}
    for r in store:
        by_key.setdefault((r["species"].strip(),), []).append(r)

    # ── 1. the RYA-1160 control, re-run from its own recorded NIST columns ──
    for c in _rows(CONTROL):
        ours, theirs = c["ours_nist_grade"].strip(), c["nist_grade"].strip()
        if c["match"] != "UNIQUE":
            notes.append(f"{c['species']} {c['wavelength_air_A']}: NIST match "
                         f"{c['match']} — no published grade to compare (correctly refused)")
            continue
        # re-read OUR grade from the live store rather than trusting the control's copy
        live = [r for r in store
                if r["species"].strip() == c["species"].strip()
                and abs(float(r["wavelength_air_A"]) - float(c["wavelength_air_A"])) <= TOL_A
                and abs(float(r["excitation_potential_eV"]) - float(c["EP_eV"])) <= TOL_EP]
        if len(live) != 1:
            notes.append(f"{c['species']} {c['wavelength_air_A']}: {len(live)} store rows "
                         f"on the EP-aware key — not comparable")
            continue
        now = (live[0]["nist_grade"] or "").strip()
        if tighter(now, theirs):
            over.append(f"{c['species']} {c['wavelength_air_A']}: store {now} "
                        f"({NIST_ACC_PCT[now]}%) is TIGHTER than NIST {theirs} "
                        f"({NIST_ACC_PCT[theirs]}%)")
        elif now != ours:
            notes.append(f"{c['species']} {c['wavelength_air_A']}: store grade corrected "
                         f"{ours} -> {now}, now agrees with NIST {theirs}")

    # ── 2. every local NIST record vs the live store, EP-aware ──
    for rec in NIST_RECORDS:
        for r in _rows(rec):
            g = (r.get("nist_grade") or "").strip()
            if not g:
                continue
            sp = f"{r['element']} {r['ion']}".strip()
            live = [s for s in store
                    if s["species"].strip() == sp
                    and abs(float(s["wavelength_air_A"]) - float(r["wavelength_air_A"])) <= TOL_A
                    and abs(float(s["excitation_potential_eV"])
                            - float(r["excitation_potential_eV"])) <= TOL_EP]
            if len(live) != 1:
                continue
            now = (live[0]["nist_grade"] or "").strip()
            if now != g:
                who = rec.name
                if tighter(now, g):
                    suspect.append(
                        f"{sp} {r['wavelength_air_A']}: store {now} "
                        f"({NIST_ACC_PCT[now]}%) is TIGHTER than {who} {g} "
                        f"({NIST_ACC_PCT[g]}%) — sigma understated "
                        f"{nist_sigma_dex(g)/nist_sigma_dex(now):.2f}x if the record is right")
                else:
                    notes.append(f"{sp} {r['wavelength_air_A']}: store {now} vs {who} {g} "
                                 f"(the RECORD is tighter — not a store over-claim)")
    return over, suspect, notes


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    over, suspect, notes = sweep()
    print(f"=== RYA-1171 grade over-claim sweep ===")
    print(f"\n  CONFIRMED over-claims vs the RYA-1160 authoritative NIST pull: {len(over)}")
    for o in over:
        print(f"    🔴 {o}")
    if not over:
        print("    none — no stored grade claims tighter accuracy than NIST publishes")

    print(f"\n  SUSPECTED over-claims vs a LOCAL NIST record: {len(suspect)}")
    print(f"    (nist_crosscheck/nist_reference are our own files, not an authority pull;")
    print(f"     confirming these needs a NIST ASD pull for those species — RYA-1160 was CNO only)")
    for x in suspect:
        print(f"    ⚠️ {x}")
    if not suspect:
        print("    none")

    print(f"\n  {len(notes)} other note(s):")
    for n in notes:
        print(f"    {n}")

    store = store_grades()
    ap_rows = [r for r in store if (r["nist_grade"] or "").strip() == "A+"]
    print(f"\n  A+ rows remaining in canonical_gf: {len(ap_rows)}")
    for r in ap_rows:
        print(f"    {r['species']:6s} {r['wavelength_air_A']:>10s}  {r['loggf_reference'][:44]}")
    print(f"\n  sigma: A+ -> {nist_sigma_dex('A+'):.6f} dex,  A -> {nist_sigma_dex('A'):.6f} dex "
          f"({nist_sigma_dex('A')/nist_sigma_dex('A+'):.2f}x)")
    return 1 if (a.check and over) else 0


if __name__ == "__main__":
    sys.exit(main())

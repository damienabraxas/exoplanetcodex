#!/usr/bin/env python3
"""
RYA-853 — the Fe II arbiter lines claim a NIST grade NIST never gave them.
==========================================================================
`canonical_gf` labels Fe II 6149.246 and 6247.557 `NIST ASD v5.11 grade B`. Live NIST ASD,
queried in AIR and EP-matched on both sides, uniquely returns

    6149.246 -> log gf -2.854, Acc. E      (we store -2.724, +0.130 above)
    6247.557 -> log gf -2.444, Acc. D      (we store -2.329, +0.115 above)

Wrong on BOTH axes, and the grade is the load-bearing half: `B` sits in NIST_GRADE_HIGH
("lab, <=10%, trusted") while the true grades sit in NIST_GRADE_CULL. A fabricated B
publishes 0.041 dex on a line whose source says 0.176-0.301.

WHERE IT CAME FROM — three files, one defect, propagating downhill:

    nist_reference.csv / nist_crosscheck.csv   the value AND the grade originate here.
              |                                 These are hand-maintained extracts; the
              |                                 header of nist_reference.csv already warns
              |                                 its rows "should not be assumed current",
              |                                 and RYA-853 scope 1 measured 70% of its
              |                                 grades wrong.
              |  build_linelist.crosscheck_nist()
              v
    linelist_solar.csv                         VALD3 value (-2.72) carrying the stamped B
              |  migrate_gf_single_source.py
              v
    canonical_gf.csv                           adjudicated to the extract: -2.724 + grade B

WHY THE VALUE IS HELD (Ryan's call, 2026-08-28). Every alternative is worse:
  * NIST covers these lines at grade E/D -- and RYA-853 scope 3 measured our Fe II scale
    AGAINST pure-lab Den Hartog 2019 and found NIST is the LOW one (NIST - DH = +0.021 on
    gf, i.e. NIST gf high / abundance low). Adopting NIST's value would move the ionization
    arbiter ~0.13 dex onto the scale we just measured as low, on grade E/D data.
  * Den Hartog 2019 is pure lab but STOPS AT 4584 A. These are 6149 / 6247.
  * Melendez & Barbuy 2009 has both, flagged `S` -- solar-fitted, FIREWALLED by RYA-161.
So there is no adjudicable source. The value stands as what it is (VALD3 scale), the false
NIST claim is removed, and the absence of a primary lab gf here is recorded as a DECLARED
GAP rather than papered over (RYA-833).

This script changes NO log_gf anywhere. It removes a claim.

Usage:
    python3 scripts/rya853_fe2_provenance_fix.py            # dry run
    python3 scripts/rya853_fe2_provenance_fix.py --apply
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
SOLAR = ROOT / "data" / "linelists" / "linelist_solar.csv"
NIST_REF = ROOT / "data" / "linelists" / "nist_reference.csv"
NIST_XC = ROOT / "data" / "linelists" / "nist_crosscheck.csv"
OUT = ROOT / "data" / "results" / "rya853"

#: (wavelength, EP) of the two lines, and what NIST actually says. EP is carried because a
#: wavelength is not a line identifier -- that is the defect this ticket is about.
TARGETS = [
    {"wl": 6149.246, "ep": 3.8892, "ours": -2.724, "nist_loggf": -2.853872, "nist_acc": "E"},
    {"wl": 6247.557, "ep": 3.8916, "ours": -2.329, "nist_loggf": -2.443697, "nist_acc": "D"},
]

WAVE_TOL_A, EP_TOL_EV = 0.05, 0.05

HONEST_REF = ("VALD3-scale value via nist_reference.csv; NOT NIST "
              "(ASD gives {nist:.3f} acc {acc}) -- RYA-853")


def _split_csv(line: str) -> list[str]:
    """csv.reader on one line, so a quoted field containing a comma survives."""
    import csv, io
    return next(csv.reader(io.StringIO(line)))


def _join_csv(fields: list[str]) -> str:
    import csv, io
    buf = io.StringIO()
    csv.writer(buf, lineterminator="").writerow(fields)
    return buf.getvalue()


def edit_rows(path: Path, match, edits: dict, apply: bool, changes: list) -> int:
    """Rewrite ONLY the lines `match` selects, in place, leaving every other byte alone.

    🔴 NOT a pandas round-trip. The first cut of this script read each file with pandas and
    wrote it back, which reformatted 46 unrelated canonical_gf rows (-1.2500 -> -1.25) and
    rewrote both extracts wholesale. A provenance fix that silently reformats rows it was
    never asked to touch is the same class of defect it is fixing. Line surgery only.
    """
    # newline='' — no universal-newline translation. nist_crosscheck.csv has MIXED
    # endings (21 CRLF, 5 LF); read_text() silently normalises them and the write-back
    # then rewrites 19 lines this script was never asked to touch.
    with open(path, newline='') as fh:
        raw = fh.readlines()
    header_idx = next(i for i, l in enumerate(raw) if not l.lstrip().startswith("#"))
    cols = _split_csv(raw[header_idx].rstrip("\n"))
    n = 0
    for i in range(header_idx + 1, len(raw)):
        term = "\r\n" if raw[i].endswith("\r\n") else ("\n" if raw[i].endswith("\n") else "")
        line = raw[i][:len(raw[i]) - len(term)]
        if not line.strip():
            continue
        f = _split_csv(line)
        if len(f) != len(cols):
            continue
        row = dict(zip(cols, f))
        hit = match(row)
        if not hit:
            continue
        n += 1
        for field, value in edits(row).items():
            if field not in cols:
                continue
            j = cols.index(field)
            changes.append({"file": path.name,
                            "wavelength_air_A": row.get("wavelength_air_A"),
                            "field": field, "before": f[j], "after": value})
            if apply:
                f[j] = value
        if apply:
            raw[i] = _join_csv(f) + term
    if apply and n:
        with open(path, "w", newline='') as fh:
            fh.write("".join(raw))
    return n


def _near(row, wl, ep, species=None) -> bool:
    """EP-guarded on BOTH sides — the defect this ticket is about is a wavelength-only
    match, so this matcher must not be one."""
    try:
        w = float(row["wavelength_air_A"]); e = float(row["excitation_potential_eV"])
    except (KeyError, ValueError, TypeError):
        return False
    if abs(w - wl) > WAVE_TOL_A or abs(e - ep) > EP_TOL_EV:
        return False
    if species is not None:
        return str(row.get("species", "")).strip() == species
    return str(row.get("element", "")).strip() == "Fe" and \
        str(row.get("ion", "")).strip() == "II"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    changes: list = []

    for t in TARGETS:
        ref = HONEST_REF.format(nist=t["nist_loggf"], acc=t["nist_acc"])
        note = (f"RYA-853: grade REMOVED — ASD gives {t['nist_loggf']:.3f} acc "
                f"{t['nist_acc']}, not the B stored here; log_gf is VALD3-scale and HELD "
                f"(docs/science/rya853_fe2_provenance.md)")

        n = edit_rows(CANON, lambda r, t=t: _near(r, t["wl"], t["ep"], species="Fe II"),
                      lambda r, ref=ref: {"loggf_reference": ref, "nist_grade": "",
                                          "adjudication_status": "held_rya853"},
                      a.apply, changes)
        if n != 1:
            raise SystemExit(f"canonical_gf: {n} rows matched {t['wl']} — refusing")

        if edit_rows(SOLAR, lambda r, t=t: _near(r, t["wl"], t["ep"]),
                     lambda r: {"nist_grade": ""}, a.apply, changes) != 1:
            raise SystemExit(f"linelist_solar: not exactly one row at {t['wl']}")

        for path in (NIST_REF, NIST_XC):
            edit_rows(path, lambda r, t=t: _near(r, t["wl"], t["ep"]),
                      lambda r, note=note: {"nist_grade": "", "notes": note},
                      a.apply, changes)

    df = pd.DataFrame(changes)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "rya853_fe2_provenance_changes.csv", index=False)
    (OUT / "rya853_fe2_declared_gap.json").write_text(json.dumps({
        "ticket": "RYA-853",
        "gap": "no primary laboratory gf covers Fe II 6149.246 or 6247.557",
        "checked_not_assumed": {
            "DenHartog2019": "pure lab, but its table stops at 4584 A — these are redward",
            "MelendezBarbuy2009": "covers both, flagged S (solar-fitted) — FIREWALLED (RYA-161)",
            "NIST_ASD": "covers both at acc E / D; values -2.854 / -2.444, "
                        "0.130 / 0.115 dex below ours",
        },
        "value_disposition": "HELD — no adjudicable source; the claim was removed, not the number",
        "so_the_sigma_is": ("uncited: gf_sigma_dex stays NaN and the band budget charges the "
                            "ungraded blanket, which is the honest state for a line with no "
                            "graded source"),
        "owner": "RYA-953 (Fe II has no primary-lab gf table)",
        "targets": [{k: t[k] for k in ("wl", "ep", "ours", "nist_loggf", "nist_acc")}
                    for t in TARGETS],
    }, indent=2) + "\n")

    print(f"=== RYA-853 Fe II provenance fix ({'APPLY' if a.apply else 'DRY RUN'}) ===")
    for f, g in df.groupby("file"):
        print(f"  {f:<24} {len(g)} field(s) on {g.wavelength_air_A.nunique()} row(s)")
    print(f"\n  log_gf touched: {int((df.field == 'log_gf').sum())} — the value is HELD")
    print(f"\n[out] {OUT}/rya853_fe2_provenance_changes.csv")


if __name__ == "__main__":
    main()

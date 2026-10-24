#!/usr/bin/env python3
"""RYA-1268 -- adopt every remaining element's published-set gf through the orchestrator's own
prepare -> apply (pipeline.element_prepare), WITHOUT culls (culls stay Ryan's yes/no) and
WITHOUT any adoption that would move a line by more than HOLD_MOVE_DEX ("super bad dex should
not be included", Ryan 2026-10-23) or land on a row another ticket adjudicated -- those are
listed in held_set_adoptions.csv.
Missing published lines are added by apply's own published_line_add plan (clash-free only).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import element_prepare as ep  # noqa: E402

HOLD_MOVE_DEX = 0.5
# finished work: a row another ticket adjudicated is held for review, never re-adopted
# (rya1268_adopt_lab_gf.OPEN_STATUS; its own adoptions stay open to the published set)
OPEN_STATUS = {"", "single_source", "pending", "agreed"}
GEM = ROOT / "data" / "reference" / "gf_error_model.csv"
REFLINES = ROOT / "data" / "reference" / "agss21_cited_sets" / "agss21_reference_lines.csv"
WISC = ("Wisconsin A = BF x 1/tau: LIF lifetimes '+/-5% (+/-0.02 dex)' per upper level, independent "
        "level to level; calibration error '0.001%/cm-1 x wavenumber difference' is per line -- "
        "the DH14 model, shared scale 2%")
LAB_ROWS = {  # lab_source_tag: (group, correlation, scale, basis, citation)
    **{k: ("WISCONSIN_" + k.upper(), "scale_only", 0.0086, WISC, c) for k, c in [
        ("wood2013_ti2", "Wood+2013 ApJS 208 27 Sect. 3"), ("lawler2014_v1", "Lawler+2014 ApJS 215 20"),
        ("wood2014_v2", "Wood+2014 ApJS 214 18 Sect. 3"), ("lawler2017_cr2", "Lawler+2017 ApJS 228 10"),
        ("lawler2015", "Lawler+2015 ApJS 220 13"), ("wood2014", "Wood+2014 ApJS 211 20 Sect. 3"),
        ("denhartog2011", "Den Hartog+2011 ApJS 194 35"), ("sobeck2007", "Sobeck+2007 ApJ 667 1267")]},
    **{k: ("NIST_ASD", "scale_only", 0.0, "lab gf (via VALD code) priced by its NIST ASD class on the same gf: "
                                         "the NIST-C+ model", c) for k, c in [
        ("lawler2019_sc", "Lawler+2019 ApJS 241 21 via VALD LHSNWC + NIST ASD"),
        ("lawler2013", "Lawler+2013 ApJS 205 11 via VALD LGWSC + NIST ASD"),
        ("lawler2001_eu2", "Lawler+2001 ApJ 563 1075 via VALD LWHS + NIST ASD")]},
    # one stated percentage for every line of a paper is a shared-scale statement: charged fully
    # correlated across that paper's lines (the conservative reading, as MG_AGSS21 / AL_AGSS21)
    **{k: (k.upper(), "fully_correlated", "", b, c) for k, b, c in [
        ("grevesse2015_Heavies", "Zr II: 'uncertainties of the mean gf-values we adopt are of order 5-10%'",
         "Grevesse+2015 A&A 573 A27 Sect. Zr"),
        ("gallagher2020_ba", "per-line gf error printed in Table 1", "Gallagher+2020 A&A 634 A55 Table 1"),
        ("storm2024_y_eu", "Eu II Lawler+2001 'estimated at 12%'; Y II Palmeri+2017 '~10%'",
         "Storm+2024 A&A 683 A200 Sect. 2")]},
}


def ensure_gf_error_rows() -> None:
    """Every gf source this ticket adopts gets a correlation class (unlisted = UNREVIEWED).
    A set tag takes its lab paper's class when every sigma in the set comes from that one
    paper; a mixed set is charged fully correlated (conservative, as MG_AGSS21 / AL_AGSS21)."""
    raw = GEM.read_bytes().decode()
    nl = "\r\n" if "\r\n" in raw else "\n"
    have = {l.split(",")[0] for l in raw.splitlines()}
    rows = dict(LAB_ROWS)
    ref = pd.read_csv(REFLINES)
    for name, g in ref.groupby("set_name"):
        papers = set(g.lab_paper.dropna()) if g.gf_sigma_basis.str.startswith("lab").all() else set()
        if len(papers) == 1:
            p = papers.pop()
            rows[name] = (LAB_ROWS[p][0], "scale_only", 0.0086, f"every sigma in the set is {p}'s A +/- dA; {WISC}",
                          LAB_ROWS[p][4])
        else:
            rows[name] = (name, "fully_correlated", "", "per-line sigma from mixed sources (lab table / NIST class / "
                          "stated in text); charged fully correlated (conservative)", g.citation.iloc[0])
    add = []
    for k, (grp, cor, sc, basis, cite) in rows.items():
        if k not in have:
            q = lambda x: '"' + str(x).replace('"', "'") + '"'
            add.append(",".join([k, grp, cor, str(sc), q(basis), q(cite), "RYA-1268"]))
    if add:
        GEM.write_bytes((raw.rstrip("\r\n") + nl + nl.join(add) + nl).encode())
ELEMENTS = ["Li", "Na", "K", "Ca", "P", "S", "Sc", "Ti", "V", "Cr", "Mn", "Co", "Ni", "Cu", "Zn",
            "Sr", "Y", "Zr", "Ba", "Eu"]


def main() -> None:
    ensure_gf_error_rows()
    if "--report" in sys.argv:  # re-run prepare only: the post-apply state of every report
        for el in ELEMENTS:
            try:
                print(el, "READY" if ep.prepare("solar", el)["ready"] else "NOT READY", flush=True)
            except Exception as exc:
                print(el, f"PREPARE FAILED: {type(exc).__name__}: {exc}", flush=True)
        return
    held, summary = [], {}
    for el in (sys.argv[1:] or ELEMENTS):
        try:
            rep = ep.prepare("solar", el)
        except Exception as exc:  # report and move on; one element must not stop the rest
            summary[el] = f"PREPARE FAILED: {type(exc).__name__}: {exc}"
            print(el, summary[el], flush=True)
            continue
        keep = []
        status = pd.read_csv(ep.CANONICAL, usecols=["line_id", "adjudication_status"], dtype=str,
                             keep_default_na=False).set_index("line_id").adjudication_status
        for a in rep["adopt"]:
            d = a.get("d_loggf")
            st = status.get(a["line_id"], "")
            why = ("gf move > %.1f dex" % HOLD_MOVE_DEX if d is not None and abs(d) > HOLD_MOVE_DEX
                   else f"row already adjudicated ({st})"
                   if st not in OPEN_STATUS and not st.startswith("adopted_rya1268_") else None)
            if why:
                held.append({"element": el, "reason": why, **{k: a.get(k) for k in (
                    "set", "species", "wavelength_A", "canonical_loggf", "authors_loggf", "d_loggf")}})
            else:
                keep.append(a)
        rep["adopt"], rep["cull_candidates"] = keep, []
        done = ep.apply("solar", el, rep)
        summary[el] = done
        print(el, done, flush=True)
    out = ROOT / "data" / "audit" / "rya1268_line_coverage"
    pd.DataFrame(held).to_csv(out / "held_set_adoptions.csv", index=False)
    (out / "apply_sets_summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()

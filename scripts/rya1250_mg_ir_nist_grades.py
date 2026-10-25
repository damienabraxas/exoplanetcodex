#!/usr/bin/env python3
"""RYA-1250: NIST ASD accuracy grades for the solar Mg IR rows of canonical_gf (Y/J, H, K).

The Codex / Deep pools take a line when its gf carries a published uncertainty
(`gf_grades.is_gf_graded`: lab, a stored sigma, or NIST A..C). The IR Mg rows were seeded
from VALD (RYA-762 / RYA-1053) and never matched to NIST, so K offered 2 graded lines where
NIST grades ~13 features (21458.86 B+, 24566.30 A, 22808 B+, ...; Nandakumar+2023 IGRINS use
21458.87 for exactly that reason).

A row takes a grade only where NIST's value IS the row's value:
  * single row: |log gf(row) - log gf(NIST)| <= GF_TOL, same lower level (EP_TOL);
  * a row holding a cluster TOTAL (hfs_n_components > 1): the cluster's dominant synthesis
    component matches NIST to GF_TOL and carries >= DOMINANT of the cluster's gf -- the
    19723/19733 precedent (`rya1250_mg_reference_lines.RESTORE_TOTAL`).
log_gf is never changed; a row already LAB / PUBLISHED-SET / NIST-graded is left alone.

    python scripts/rya1250_mg_ir_nist_grades.py            # report
    python scripts/rya1250_mg_ir_nist_grades.py --apply    # write canonical_gf + snapshot

Needs the network (astroquery.nist). The NIST rows used are frozen in
data/reference/nist/mg_ir_rya1250.csv, which --offline re-applies without querying.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.gf_grades import GRADED_NIST_CLASSES  # noqa: E402

CANON = ROOT / "data/linelists/canonical_gf.csv"
SNAP = ROOT / "data/reference/nist/mg_ir_rya1250.csv"
#: (lo, hi, synthesis list) -- the IR lists the solar Mg cells read.
BANDS = [(9200, 13000, "ispec_ir_9200_13000"), (15000, 17500, "ispec_h_15007_17494"),
         (19450, 24850, "ispec_k_19452_24846")]
W_TOL, EP_TOL, GF_TOL, DOMINANT = 0.15, 0.01, 0.03, 0.5
KEEP_TIERS = {"LAB", "PUBLISHED-SET", "NIST-C+", "PR2024-CALC"}


def _num(x) -> float:
    try:
        return float(str(x).strip().strip("[]()+?*"))
    except ValueError:
        return float("nan")


def query_nist() -> pd.DataFrame:
    from astroquery.nist import Nist
    import astropy.units as u
    out = []
    for sp in ("Mg I", "Mg II"):
        for lo, hi, _ in BANDS:
            t = Nist.query(lo * u.AA, hi * u.AA, linename=sp, wavelength_type="vac+air")
            for r in t:
                acc = str(r["Acc."]).strip()
                if acc not in GRADED_NIST_CLASSES:
                    continue
                w = _num(r["Observed"])
                w = w if np.isfinite(w) else _num(r["Ritz"])
                gi = _num(str(r["gi   gk"]).split("-")[0])
                f = _num(r["fik"])
                ei = _num(str(r["Ei           Ek"]).split("-")[0])
                if not (np.isfinite(w) and np.isfinite(gi) and np.isfinite(f) and f > 0):
                    continue
                out.append(dict(species=sp, wavelength_air_A=round(w, 3), ei_eV=ei,
                                log_gf=round(float(np.log10(gi * f)), 3), acc=acc,
                                lower=" ".join(str(r["Lower level"]).split()),
                                upper=" ".join(str(r["Upper level"]).split()),
                                tp=str(r["TP"]).strip()))
    return pd.DataFrame(out)


def _components(lst: str, species: str) -> pd.DataFrame:
    d = pd.read_csv(ROOT / "data/linelists" / lst / "atomic_lines.tsv", sep="\t",
                    usecols=["element", "wave_A", "loggf", "lower_state_eV"], low_memory=False)
    return d[d.element == ("Mg 1" if species == "Mg I" else "Mg 2")]


def plan(nist: pd.DataFrame, canon: pd.DataFrame) -> list[dict]:
    out = []
    for lo, hi, lst in BANDS:
        for sp in ("Mg I", "Mg II"):
            comps = _components(lst, sp)
            rows = canon[(canon.species == sp) & canon.wavelength_air_A.between(lo, hi)]
            for _, r in rows.iterrows():
                if str(r.gf_tier) in KEEP_TIERS or str(r.nist_grade) in GRADED_NIST_CLASSES:
                    continue
                w, ep, g = float(r.wavelength_air_A), float(r.excitation_potential_eV), float(r.log_gf)
                cand = nist[(nist.species == sp) & ((nist.wavelength_air_A - w).abs() <= W_TOL)
                            & ((nist.ei_eV - ep).abs() <= EP_TOL)]
                if cand.empty:
                    continue
                how, n = None, None
                single = cand[(cand.log_gf - g).abs() <= GF_TOL]
                if len(single):
                    n, how = single.iloc[(single.wavelength_air_A - w).abs().argmin()], "the line"
                elif int(r.hfs_n_components or 1) > 1:
                    cl = comps[((comps.wave_A - w).abs() <= 0.10)
                               & ((comps.lower_state_eV - ep).abs() <= EP_TOL)]
                    if len(cl):
                        dom = cl.loc[cl.loggf.idxmax()]
                        share = 10 ** dom.loggf / (10 ** cl.loggf).sum()
                        m = cand[(cand.log_gf - dom.loggf).abs() <= GF_TOL]
                        if share >= DOMINANT and len(m):
                            n = m.iloc[(m.wavelength_air_A - dom.wave_A).abs().argmin()]
                            how = f"the dominant component ({share:.0%} of the line)"
                if n is None:
                    continue
                out.append(dict(line_id=r.line_id, species=sp, wavelength_air_A=w, log_gf=g,
                                nist_w=n.wavelength_air_A, nist_log_gf=n.log_gf, acc=n.acc,
                                tp=n.tp, how=how, prior_tier=r.gf_tier))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--offline", action="store_true", help="use the frozen NIST snapshot")
    a = ap.parse_args()
    nist = pd.read_csv(SNAP) if a.offline else query_nist()
    canon = pd.read_csv(CANON, low_memory=False)
    p = plan(nist, canon)
    for x in p:
        print(f"{x['species']:5} {x['wavelength_air_A']:10.3f} {x['log_gf']:+.3f}  NIST {x['acc']:3} "
              f"{x['nist_log_gf']:+.3f} @ {x['nist_w']:.3f}  ({x['how']}; was {x['prior_tier']})")
    print(f"{len(p)} rows to grade")
    if not a.apply:
        return 0
    if not a.offline:
        SNAP.parent.mkdir(parents=True, exist_ok=True)
        nist.to_csv(SNAP, index=False)
    with CANON.open(newline="") as fh:
        rows = list(csv.reader(fh))
    head, body = rows[0], rows[1:]
    col = {k: i for i, k in enumerate(head)}
    by_id = {r[0]: r for r in body}
    for x in p:
        r = by_id[x["line_id"]]
        prior = r[col["loggf_reference"]].split(" | prior: ")[-1]
        r[col["nist_grade"]] = x["acc"]
        r[col["gf_tier"]] = "NIST-C+"
        r[col["adjudication_status"]] = "adopted_rya1250_nist_ir_grade"
        r[col["loggf_reference"]] = (f"NIST ASD {x['tp']} (grade of {x['how']}; gf unchanged, "
                                     f"RYA-1250 IR sweep)" + (f" | prior: {prior}" if prior else ""))
    with CANON.open("w", newline="") as fh:
        csv.writer(fh, lineterminator="\n").writerows([head, *body])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

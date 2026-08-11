#!/usr/bin/env python3
"""RYA-773 step 1 — the LEVEL-REACH check for Al: does the model atom already know the
levels of the lines we measure but the departure grid does not serve?

    python3 scripts/rya773_al_level_reach.py

WHY THIS RUNS FIRST (the RYA-763 move, one element over)
--------------------------------------------------------
`Al_Amarsi2020_PySME.csv` is a per-LINE extract: it carries departures for 6696/6698 and
nothing else. The underlying Nordlander & Lind (2017) atom behind the Amarsi-2020 grid is
per-LEVEL. Departure coefficients are a property of levels, not wavelengths — so whether
extending the extract to Al's clean doublet is cheap or expensive is decided entirely by
whether the atom already models the levels those transitions connect:

  LEVELS PRESENT  -> the departures already exist; extending is a per-line delta run.
  LEVELS ABSENT   -> the atom must be rebuilt/enlarged offline. An atom-upgrade decision,
                     not something to undertake silently.

WHAT IT MATCHES ON, AND WHY NOT THE ENERGY ALONE
------------------------------------------------
RYA-763's lesson was that a level INDEX is a false coordinate across decks — key on the
LABEL. Here the label is (configuration, term, J, energy), and energy alone is NOT
sufficient for Al: the atom carries `3s2.5g 2G` only 5.5 meV above `3s2.5f 2F*`, and
`3s2.6g 2G` 3.4 meV above `3s2.6f 2F*`. A nearest-energy match at any tolerance wider
than ~3 meV is ambiguous between them.

The E1 selection rule breaks the tie for free and is not a tolerance at all: from a lower
`3s2.3d 2D` level (l=2) an electric-dipole transition must land on l=1 or l=3, so 2F* is
allowed and 2G (l=4) is forbidden. This script therefore requires BOTH — an energy match
AND selection-rule consistency — and reports the margin between the accepted level and
the nearest rejected one, so the tolerance is measured rather than assumed.

INPUTS, all read from files (no hardcoded wavelengths — RYA-773 guardrail)
    data/measured/band_ew/AlI_*_ew.csv            the lines we actually measure
    data/linelists/linelist_solar.csv             EP + loggf for those lines
    data/nlte_grids/amarsi_galah/label_Al.txt     the atom's levels (shipped NIST table)
    data/nlte_grids/Al_Amarsi2020_PySME.csv       what the extract currently serves
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.wavelength_util import air_to_vac   # noqa: E402

# hc in eV.A — from scipy's CODATA constants, never a typed-in 12398.4.
import scipy.constants as _sc                     # noqa: E402
HC_EV_A = _sc.h * _sc.c / _sc.e * 1e10

ELEMENT, ION = "Al", "I"
LABELS = ROOT / "data" / "nlte_grids" / "amarsi_galah" / f"label_{ELEMENT}.txt"
EXTRACT = ROOT / "data" / "nlte_grids" / f"{ELEMENT}_Amarsi2020_PySME.csv"
LINELIST = ROOT / "data" / "linelists" / "linelist_solar.csv"
MEASURED = str(ROOT / "data" / "measured" / "band_ew" / f"{ELEMENT}{ION}_*_ew.csv")

# Orbital angular momentum of the outermost electron, read off the configuration's last
# subshell (…3d, …6f, …5g). Used only for the E1 rule, so a config we cannot parse yields
# None and the rule abstains rather than guessing.
_L = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4, "h": 5, "i": 6}


def orbital_l(conf: str):
    """l of the last subshell in a configuration string, or None if unparseable."""
    tail = conf.split(".")[-1].strip()
    for ch in tail:
        if ch in _L:
            return _L[ch]
    return None


def read_levels(path: Path = LABELS) -> pd.DataFrame:
    """The atom's level table: index, species, configuration, term, J, energy (eV).

    Whitespace-delimited, and the term is NOT a single token: a Rydberg series level is
    written `3s2.nd   y 2D`, where the `y` is part of the TERM, not the configuration.
    The `.grd` itself encodes exactly that split (conf='3s2.nd', term='y 2D'), so a
    parser that put the series letter on the configuration would produce labels that
    disagree with the grid's own — and label agreement is the entire basis of the
    identification. The row is therefore parsed from both ends, with the series letter
    claimed by the term.
    """
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        p = line.split()
        term_start = -3
        if len(p[-4]) == 1 and p[-4].isalpha() and p[-4].islower():
            term_start = -4                       # series designation belongs to the term
        rows.append(dict(idx=int(p[0]), species=f"{p[1]} {p[2]}",
                         conf=" ".join(p[3:term_start]),
                         term=" ".join(p[term_start:-2]),
                         J=float(p[-2]), E_eV=float(p[-1])))
    df = pd.DataFrame(rows)
    df["l"] = df.conf.map(orbital_l)
    df["odd"] = df.term.str.endswith("*")
    return df


def e1_allowed(lo, up) -> bool:
    """Electric-dipole selection rules between two level rows: parity must flip, |dJ| <= 1
    with J=0 -> J=0 forbidden, and dl = +/-1 where both l are known."""
    if bool(lo.odd) == bool(up.odd):
        return False
    dJ = abs(float(lo.J) - float(up.J))
    if dJ > 1.0 + 1e-9 or (float(lo.J) == 0.0 and float(up.J) == 0.0):
        return False
    if lo.l is not None and up.l is not None and abs(int(lo.l) - int(up.l)) != 1:
        return False
    return True


def match_level(levels: pd.DataFrame, energy_eV: float, *, lower=None, tol_eV: float):
    """Nearest level to `energy_eV`, restricted to E1-allowed partners when `lower` is
    given. Returns (level, resid_eV, margin_eV) — margin is the energy gap to the nearest
    level this match BEAT, which is the number that says whether the tolerance is safe.
    Returns (None, nan, nan) when nothing is within tolerance.
    """
    cand = levels.copy()
    if lower is not None:
        cand = cand[cand.apply(lambda r: e1_allowed(lower, r), axis=1)]
    if not len(cand):
        return None, float("nan"), float("nan")
    d = (cand.E_eV - energy_eV).abs()
    i = int(d.idxmin())
    best, resid = cand.loc[i], float(d.loc[i])
    if resid > tol_eV:
        return None, resid, float("nan")
    # Margin = distance to the nearest level of a DIFFERENT (configuration, term).
    #
    # Deliberately not "the nearest other level": the fine-structure siblings of the
    # accepted term (e.g. 6f 2F* J=2.5 and J=3.5, which this atom lists at the SAME
    # energy to 7 decimals) sit ~0 eV away and would report a margin of 0.00 meV. That
    # is not an ambiguity — both siblings are the same term with the same departures,
    # and where both are E1-allowed the line genuinely has two components. The
    # identification that has to be safe is the TERM, so the margin is measured to the
    # nearest competing term, which for Al is the 5g/6g pair a few meV up.
    other = levels[(levels.conf != best.conf) | (levels.term != best.term)]
    margin = float((other.E_eV - energy_eV).abs().min()) if len(other) else float("nan")
    return best, resid, margin


def measured_lines() -> pd.DataFrame:
    """Every Al line the pipeline has measured, with its in-aggregate disposition."""
    files = sorted(glob.glob(MEASURED))
    if not files:
        raise SystemExit(f"no measured Al EW files matching {MEASURED}")
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    d = d[(d.element == ELEMENT) & (d.ion == ION)]
    return d.sort_values("wavelength_air_A").reset_index(drop=True)


def linelist_entries(waves, tol_A: float = 0.02) -> dict:
    """EP / loggf / vdW for each measured wavelength, from the project's own line list.

    A feature can carry several catalogue components (7836.134 has two J-components from
    one lower level). All are returned, strongest first — the strongest sets the level
    identification and the rest are reported, because a second component sharing the
    lower level is evidence FOR the identification, not noise.
    """
    ll = pd.read_csv(LINELIST, comment="#", low_memory=False)
    ll = ll[(ll.element == ELEMENT) & (ll.ion == ION)]
    out = {}
    for w in waves:
        m = ll[(ll.wavelength_air_A - w).abs() <= tol_A]
        out[w] = m.sort_values("log_gf", ascending=False).reset_index(drop=True)
    return out


def served_waves(tol_A: float = 0.15) -> np.ndarray:
    """Wavelengths the current departure extract serves (the registry matches at 0.15 A —
    `rya716_al_products.py`'s ENGINE-A UNCOVERED disposition uses the same window)."""
    if not EXTRACT.exists():
        return np.array([])
    return np.unique(pd.read_csv(EXTRACT).wave_A.to_numpy())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tol-eV", type=float, default=0.006,
                    help="energy tolerance for a level identification. Default 0.006 eV "
                         "is deliberately WIDER than the 3.4 meV 6f/6g near-degeneracy, "
                         "to prove the selection rule (not the tolerance) is what "
                         "resolves it.")
    ap.add_argument("--json", metavar="PATH", help="also write the table as JSON")
    a = ap.parse_args()

    levels = read_levels()
    n_i = int((levels.species == f"{ELEMENT} 1").sum())
    print(f"atom: {len(levels)} levels in {LABELS.name} "
          f"({n_i} {ELEMENT} I, {len(levels) - n_i} {ELEMENT} II)  "
          f"E {levels.E_eV.min():.3f}-{levels.E_eV.max():.3f} eV")

    meas = measured_lines()
    ent = linelist_entries(meas.wavelength_air_A.to_numpy())
    served = served_waves()
    print(f"measured: {len(meas)} Al I lines ({int(meas.in_aggregate.sum())} in-aggregate)")
    print(f"extract:  {EXTRACT.name} serves {len(served)} wavelengths {list(served)}\n")

    print(f"{'line':>10} {'agg':>4} {'served':>7}  {'lower level':<26} "
          f"{'upper level':<26} {'resid_ueV':>10} {'margin_meV':>11}  verdict")
    rows = []
    for _, m in meas.iterrows():
        w = float(m.wavelength_air_A)
        cat = ent[w]
        is_served = bool(len(served) and np.min(np.abs(served - w)) <= 0.15)
        if not len(cat):
            print(f"{w:10.3f} {str(bool(m.in_aggregate)):>4} {str(is_served):>7}  "
                  f"{'-- not in the line list --':<26}")
            rows.append(dict(wave_A=w, verdict="NO-CATALOGUE-ENTRY"))
            continue
        c0 = cat.iloc[0]
        elow = float(c0.excitation_potential_eV)
        lvac = float(air_to_vac(np.array([w]))[0])
        eup = elow + HC_EV_A / lvac

        lo, lo_res, _ = match_level(levels, elow, tol_eV=a.tol_eV)
        if lo is None:
            print(f"{w:10.3f} {str(bool(m.in_aggregate)):>4} {str(is_served):>7}  "
                  f"LOWER ABSENT (nearest {lo_res * 1e3:.1f} meV)")
            rows.append(dict(wave_A=w, verdict="LOWER-LEVEL-ABSENT"))
            continue
        up, up_res, up_margin = match_level(levels, eup, lower=lo, tol_eV=a.tol_eV)
        if up is None:
            print(f"{w:10.3f} {str(bool(m.in_aggregate)):>4} {str(is_served):>7}  "
                  f"{lo.conf + ' ' + lo.term + f' J={lo.J}':<26} "
                  f"UPPER ABSENT (nearest allowed {up_res * 1e3:.1f} meV)")
            rows.append(dict(wave_A=w, verdict="UPPER-LEVEL-ABSENT",
                             upper_nearest_meV=up_res * 1e3))
            continue

        verdict = "PRESENT" + ("" if is_served else "  <- reachable, NOT served")
        lo_s = f"{lo.idx:3d} {lo.conf} {lo.term} J={lo.J}"
        up_s = f"{up.idx:3d} {up.conf} {up.term} J={up.J}"
        print(f"{w:10.3f} {str(bool(m.in_aggregate)):>4} {str(is_served):>7}  "
              f"{lo_s:<26} {up_s:<26} {up_res * 1e6:10.1f} {up_margin * 1e3:11.2f}  "
              f"{verdict}")
        rows.append(dict(
            wave_A=w, in_aggregate=bool(m.in_aggregate), served=is_served,
            n_catalogue_components=len(cat), loggf=float(c0.log_gf),
            damping_vdW=float(c0.damping_vdW),
            lower_idx=int(lo.idx), lower=f"{lo.conf} {lo.term}", lower_J=float(lo.J),
            upper_idx=int(up.idx), upper=f"{up.conf} {up.term}", upper_J=float(up.J),
            upper_resid_ueV=up_res * 1e6, upper_margin_meV=up_margin * 1e3,
            verdict="PRESENT"))

    d = pd.DataFrame(rows)
    present = d[d.verdict == "PRESENT"] if "verdict" in d else d
    gap = present[(~present.served.fillna(False)) & present.in_aggregate.fillna(False)] \
        if len(present) else present

    print(f"\n{'=' * 78}")
    print(f"LEVEL-REACH VERDICT: {len(present)} of {len(d)} measured lines resolve to "
          f"levels the atom already carries.")
    if len(present):
        print(f"  worst upper-level residual : {present.upper_resid_ueV.max():.1f} ueV")
        print(f"  smallest margin to the next level : "
              f"{present.upper_margin_meV.min():.2f} meV  "
              f"(what the E1 selection rule, not the tolerance, resolves)")
    print(f"  IN-AGGREGATE AND REACHABLE BUT UNSERVED: {len(gap)} line(s) "
          f"{sorted(gap.wave_A.tolist()) if len(gap) else ''}")
    print("  => " + ("CHEAP: extending the extract is a per-line delta run, no new atom."
                     if len(gap) else
                     "nothing to extend — every in-aggregate line is already served."))

    if a.json:
        Path(a.json).write_text(json.dumps(rows, indent=2, default=float))
        print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()

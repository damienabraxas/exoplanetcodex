#!/usr/bin/env python3
"""RYA-1233 -- extend the Si I ENGINE-A extract to every SOLAR-GRADED Si I line. SIRIUS ONLY.

    python3 scripts/rya1233_extend_si_nlte.py --solar-only   # non-regression + new lines, solar
    python3 scripts/rya1233_extend_si_nlte.py --write        # all 11 nodes, append to the CSV

WHAT THIS FIXES
---------------
`Si_Amarsi2020_PySME.csv` (RYA-410) serves seven optical lines 5622-6254 A. The published
Si products measure Asplund's set (5645-7226 A), Deshmukh 2022's IR lines and Elgueta 2026's
solar-graded IR lines, so ENGINE-A covered 2 of Asplund's 7 VIS lines, nothing red-optical
and nothing IR: the "NLTE" headline rested on 2 lines (+/-0.102) while the 7-line LTE
product read +/-0.067. Same story as Al (RYA-773): a per-LINE extract of a per-LEVEL grid.

SAME METHOD AS RYA-410, ONE CODE PATH
-------------------------------------
pipeline.pysme_nlte.nlte_delta: PySME synthesis NLTE vs LTE on the Amarsi et al. 2020 GALAH
departure grid (nlte_Si_scatt_pysme.grd, Zenodo 3982506), delta = A_NLTE - A_LTE by EW
matching. Grid level labels: lower by nearest energy, upper by the LS SELECTION RULES
(pysme_nlte.labels_by_selection_rules) -- RYA-410's nearest-energy auto_labels gave four Si
lines an unreachable upper level (6125 J1->5, 6741 J3->1, 7034 J2->4, 7226 J1->3), which
is the case Amarsi & Asplund 2017 resolve by the same rules. Same 11 nodes as RYA-410.

Atomic data: log gf and EP from canonical_gf -- the values the published measurement
synthesised with (Garz+0.097 optical, Pehlivan Rhodin 2024 IR); E_up = EP + hc/lambda.
vdW from the synthesis linelist covering the line (GES 4200-9200 A, ispec_ir / ispec_h
beyond), with RYA-410's rule (a non-positive value -> 0, PySME's default).

It does NOT touch the 77 existing rows. Non-regression: the 7 banked lines are re-derived
at the solar node by this script's own atomic-data path and compared (--solar-only).
Si II 6371 is not here: these are Si I level populations, and Amarsi & Asplund 2017 put
its 3D NLTE correction at zero.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import pysme_nlte as pn          # noqa: E402

CSV = ROOT / "data" / "nlte_grids" / "Si_Amarsi2020_PySME.csv"
PROV = ROOT / "data" / "nlte_grids" / "Si_Amarsi2020_PySME.prov.json"
ELEMENT = "Si"
SETS = ("si_agss21_SiI", "si_deshmukh2022_SiI", "si_elgueta2026_SiI")
GES = Path(os.environ.get("ISPEC_DIR", "/mnt/codex-data/engines/ispec_src")) / \
    "input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv"
IR_LISTS = (ROOT / "data/linelists/ispec_ir_9200_13000/atomic_lines.tsv",
            ROOT / "data/linelists/ispec_h_15007_17494/atomic_lines.tsv")
NODES = [(5772, 4.44, 0.0), (5100, 4.44, 0.0), (6200, 4.44, 0.0), (5772, 4.0, 0.0),
         (5772, 4.7, 0.0), (5772, 4.44, -0.3), (5172, 4.43, 0.31), (5100, 4.5, 0.5),
         (6000, 4.3, -0.2), (5772, 4.44, 0.6), (5400, 4.5, 0.5)]
SOLAR = {"teff": 5772, "logg": 4.44, "feh": 0.0, "vmic": 1.0}


def graded_lines() -> list[float]:
    """Every Si I line any solar-graded set measures (the solar graded file where one exists)."""
    from pipeline.run_descriptor import star_graded_csv
    reg = {r["csv"]: r for r in csv.DictReader(open(ROOT / "data/reference/line_sets/REGISTRY.csv"))}
    waves = set()
    for s in SETS:
        row = reg[f"data/reference/line_sets/{s}.csv"]
        for r in csv.DictReader(open(ROOT / star_graded_csv(row, "solar"))):
            waves.add(round(float(r["wavelength_air_A"]), 3))
    return sorted(waves)


_LL: dict = {}


def _synth_rows(wl: float) -> pd.DataFrame:
    path = GES if wl < 9200 else (IR_LISTS[0] if wl < 15000 else IR_LISTS[1])
    if path not in _LL:
        d = pd.read_csv(path, sep="\t", low_memory=False,
                        usecols=["element", "wave_A", "lower_state_eV", "waals"])
        _LL[path] = d[d.element == "Si 1"]
    return _LL[path]


def atomic(wl: float):
    """(log gf, EP, E_up, vdW) -- gf/EP from canonical_gf, vdW from the synthesis list."""
    cg = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False,
                     usecols=["species", "wavelength_air_A", "excitation_potential_eV", "log_gf"])
    m = cg[(cg.species == "Si I") & ((cg.wavelength_air_A - wl).abs() < 0.01)]
    if len(m) != 1:
        raise ValueError(f"Si I {wl}: {len(m)} canonical rows within 0.01 A")
    loggf, ep = float(m.log_gf.iloc[0]), float(m.excitation_potential_eV.iloc[0])
    ll = _synth_rows(wl)
    s = ll[((ll.wave_A - wl).abs() < 0.05) & ((ll.lower_state_eV - ep).abs() < 0.03)]
    vw = float(s.waals.iloc[0]) if len(s) else 0.0
    return loggf, ep, ep + 12398.42 / wl, (vw if vw > 0 else 0.0)


def build_lines(waves):
    out = []
    for wl in waves:
        loggf, ep, eup, vw = atomic(wl)
        tl, tu, jl, ju, note = pn.labels_by_selection_rules(ELEMENT, ep, eup)
        out.append((wl, loggf, ep, jl, eup, ju, tl, tu, vw))
        print(f"  {wl:10.3f} loggf={loggf:+.3f} EP={ep:.3f} {tl} -> {tu} "
              f"(J{jl:.0f}->{ju:.0f}) vdW={vw:.1f}" + (f"  [{note}]" if note else ""),
              flush=True)
    return out


def derive(star: dict, lines) -> dict:
    """delta per line at one node, derived in groups a few A wide (RYA-773)."""
    out = {}
    lines = sorted(lines, key=lambda l: l[0])
    grp = [lines[0]]
    for l in lines[1:]:
        if l[0] - grp[0][0] > 30.0:
            out.update(pn.nlte_delta(ELEMENT, star=star, lines=grp)["per_line"])
            grp = [l]
        else:
            grp.append(l)
    out.update(pn.nlte_delta(ELEMENT, star=star, lines=grp)["per_line"])
    return out


def _node(args):
    (te, lg, fe), lines = args
    from config.constants import SOLAR_ASPLUND2021
    pn._A_SUN.setdefault(ELEMENT, SOLAR_ASPLUND2021[ELEMENT])
    pn.NLTE_LINES[ELEMENT] = list(lines)
    d = derive({"teff": te, "logg": lg, "feh": fe, "vmic": 1.0}, lines)
    print(f"  node {te}/{lg}/{fe:+.2f} done", flush=True)
    return d


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--solar-only", action="store_true")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--jobs", type=int, default=3)
    a = ap.parse_args()
    from config.constants import SOLAR_ASPLUND2021
    pn._A_SUN.setdefault(ELEMENT, SOLAR_ASPLUND2021[ELEMENT])
    old = pd.read_csv(CSV)
    banked = sorted(old.wave_A.unique())
    new = [w for w in graded_lines() if not any(abs(w - b) < 0.01 for b in banked)]
    print(f"{len(banked)} banked lines, {len(new)} new solar-graded Si I lines to derive")
    new_lines = build_lines(new)
    pn.NLTE_LINES[ELEMENT] = list(new_lines)

    if a.solar_only:
        print("\nNON-REGRESSION -- the banked lines, re-derived by this atomic-data path:")
        bl = build_lines(banked)
        d = derive(SOLAR, bl)
        sol = old[(old.teff_K == 5772) & (old.logg == 4.44) & (old.feh == 0.0)]
        for w in banked:
            ref = float(sol[sol.wave_A == w].delta_nlte.iloc[0])
            print(f"  {w:9.3f} banked {ref:+.4f} re-derived {d[w]:+.4f} diff {d[w] - ref:+.4f}")
        print("\nNEW (solar node):")
        dn = derive(SOLAR, new_lines)
        for w in sorted(dn):
            print(f"  {w:9.3f} delta {dn[w]:+.4f}")
        return 0

    if not a.write:
        raise SystemExit("pass --solar-only or --write")
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(a.jobs) as ex:
        results = list(ex.map(_node, [(n, new_lines) for n in NODES]))
    rows = []
    for (te, lg, fe), d in zip(NODES, results):       # node order, then wavelength
        for w in sorted(d):
            rows.append(dict(element=ELEMENT, ion=1, wave_A=round(w, 3), teff_K=te, logg=lg,
                             feh=fe, delta_nlte=round(float(d[w]), 4)))
    out = pd.concat([old, pd.DataFrame(rows)], ignore_index=True)
    out.to_csv(CSV, index=False)
    print(f"wrote {CSV.name}: {len(old)} kept + {len(rows)} new = {len(out)} rows")
    prov = json.loads(PROV.read_text())
    prov.setdefault("extensions", []).append({
        "ticket": "RYA-1233", "date": str(date.today()),
        "added_waves": sorted({r["wave_A"] for r in rows}),
        "why": ("ENGINE-A served 2 of Asplund's 7 VIS lines and no red-optical or IR line; "
                "the published Si products measure Asplund 2021, Deshmukh 2022 and "
                "Elgueta 2026 (solar-graded) lines."),
        "atomic_data": ("log gf + EP from canonical_gf (the values the measurement used); "
                        "vdW from the covering synthesis linelist (RYA-410 rule); labels by "
                        "pysme_nlte.auto_labels"),
        "existing_rows": "untouched; banked lines re-derived at solar for non-regression",
        "script": "scripts/rya1233_extend_si_nlte.py"})
    PROV.write_text(json.dumps(prov, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

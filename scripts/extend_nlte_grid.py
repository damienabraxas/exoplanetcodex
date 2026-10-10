#!/usr/bin/env python3
"""Extend an element's ENGINE-A (Amarsi 2020 / PySME) table to the lines prepare says it lacks.
SIRIUS ONLY (PySME: /mnt/codex-data/venv_pysme).

    python scripts/extend_nlte_grid.py --element Al --solar-only   # non-regression + new deltas
    python scripts/extend_nlte_grid.py --element Al --write [--jobs 3]

RYA-1233. The Si extension (scripts/rya1233_extend_si_nlte.py) hard-coded Si and its three
line-set files, so Al would have needed a copy. The line list here is the prepare report's
own `nlte_missing` (data/results/orchestrator/prepare/<star>_<El>.json, check D): every
gf-graded neutral line of the element, not culled for the star, with no row in the table.
Run `run_pipeline.py --prepare` first; re-run it after --write and check D reads 0.

Atomic data exactly as the Si extension (RYA-410 rule): log gf + EP from canonical_gf -- the
values the measurement uses -- vdW from the covering synthesis list, level labels by LS
selection rules (pysme_nlte.labels_by_selection_rules). Existing rows are never touched;
--solar-only re-derives the banked lines at the solar node and prints the difference.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import pysme_nlte as pn                     # noqa: E402
from pipeline.element_prepare import SYNTH_LISTS, OUT     # noqa: E402

GES = Path(os.environ.get("ISPEC_DIR", "/mnt/codex-data/engines/ispec_src")) / \
    "input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv"
#: The node set the Si extension used (solar, Teff/logg/[Fe/H] steps, the benchmark stars).
NODES = [(5772, 4.44, 0.0), (5100, 4.44, 0.0), (6200, 4.44, 0.0), (5772, 4.0, 0.0),
         (5772, 4.7, 0.0), (5772, 4.44, -0.3), (5172, 4.43, 0.31), (5100, 4.5, 0.5),
         (6000, 4.3, -0.2), (5772, 4.44, 0.6), (5400, 4.5, 0.5)]
SOLAR = {"teff": 5772, "logg": 4.44, "feh": 0.0, "vmic": 1.0}
_LL: dict = {}
ELEMENT = ""


def _synth_rows(wl: float) -> pd.DataFrame:
    path = GES
    for lo, f in SYNTH_LISTS:
        if wl >= lo:
            path = f
    if path not in _LL:
        d = pd.read_csv(path, sep="\t", low_memory=False,
                        usecols=["element", "wave_A", "lower_state_eV", "waals"])
        _LL[path] = d[d.element == f"{ELEMENT} 1"]
    return _LL[path]


def atomic(wl: float):
    """(log gf, EP, E_up, vdW) -- gf/EP from canonical_gf, vdW from the synthesis list."""
    cg = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False,
                     usecols=["species", "wavelength_air_A", "excitation_potential_eV", "log_gf"])
    m = cg[(cg.species == f"{ELEMENT} I") & ((cg.wavelength_air_A - wl).abs() < 0.01)]
    if len(m) != 1:
        raise ValueError(f"{ELEMENT} I {wl}: {len(m)} canonical rows within 0.01 A")
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
              f"(J{jl:.1f}->{ju:.1f}) vdW={vw:.1f}" + (f"  [{note}]" if note else ""), flush=True)
    return out


def derive(star: dict, lines) -> dict:
    """delta per line at one node, derived in groups <= 30 A wide (RYA-773)."""
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
    (te, lg, fe), lines, el = args
    global ELEMENT
    ELEMENT = el
    from config.constants import SOLAR_ASPLUND2021
    pn._A_SUN.setdefault(el, SOLAR_ASPLUND2021[el])
    pn.NLTE_LINES[el] = list(lines)
    d = derive({"teff": te, "logg": lg, "feh": fe, "vmic": 1.0}, lines)
    print(f"  node {te}/{lg}/{fe:+.2f} done", flush=True)
    return d


def main(argv=None) -> int:
    global ELEMENT
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--element", required=True)
    ap.add_argument("--star", default="solar")
    ap.add_argument("--solar-only", action="store_true")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--jobs", type=int, default=3)
    a = ap.parse_args(argv)
    ELEMENT = a.element
    csv_path = ROOT / "data" / "nlte_grids" / f"{ELEMENT}_Amarsi2020_PySME.csv"
    prov_path = csv_path.with_suffix(".prov.json")
    rep = json.loads((OUT / f"{a.star}_{ELEMENT}.json").read_text())
    from config.constants import SOLAR_ASPLUND2021
    pn._A_SUN.setdefault(ELEMENT, SOLAR_ASPLUND2021[ELEMENT])
    old = pd.read_csv(csv_path)
    banked = sorted(old.wave_A.unique())
    new = sorted({round(m["wavelength_A"], 3) for m in rep["nlte_missing"]
                  if m["species"] == f"{ELEMENT} I"})
    print(f"{len(banked)} banked lines, {len(new)} lines prepare says the table lacks")
    if not new:
        return 0
    new_lines = build_lines(new)

    if a.solar_only:
        print("\nNON-REGRESSION -- the banked lines, re-derived by this atomic-data path:")
        bl = build_lines(banked)
        pn.NLTE_LINES[ELEMENT] = list(bl)
        d = derive(SOLAR, bl)
        sol = old[(old.teff_K == 5772) & (old.logg == 4.44) & (old.feh == 0.0)]
        for w in banked:
            ref = sol[sol.wave_A == w].delta_nlte
            print(f"  {w:9.3f} banked {float(ref.iloc[0]) if len(ref) else float('nan'):+.4f} "
                  f"re-derived {d[w]:+.4f}")
        print("\nNEW (solar node):")
        pn.NLTE_LINES[ELEMENT] = list(new_lines)
        for w, v in sorted(derive(SOLAR, new_lines).items()):
            print(f"  {w:9.3f} delta {v:+.4f}")
        return 0

    if not a.write:
        raise SystemExit("pass --solar-only or --write")
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(a.jobs) as ex:
        results = list(ex.map(_node, [(n, new_lines, ELEMENT) for n in NODES]))
    # A line with no finite delta at any node is not written (pysme_nlte returns NaN when the
    # NLTE EW falls outside the LTE curve of growth); prepare's check D keeps naming it.
    rows = [dict(element=ELEMENT, ion=1, wave_A=round(w, 3), teff_K=te, logg=lg, feh=fe,
                 delta_nlte=round(float(d[w]), 4))
            for (te, lg, fe), d in zip(NODES, results) for w in sorted(d)
            if all(np.isfinite(dd.get(w, np.nan)) for dd in results)]
    out = pd.concat([old, pd.DataFrame(rows)], ignore_index=True)
    out.to_csv(csv_path, index=False)
    print(f"wrote {csv_path.name}: {len(old)} kept + {len(rows)} new = {len(out)} rows")
    prov = json.loads(prov_path.read_text()) if prov_path.exists() else {}
    prov.setdefault("extensions", []).append({
        "ticket": "RYA-1233", "date": str(date.today()),
        "added_waves": sorted({r["wave_A"] for r in rows}),
        "why": f"prepare check D: {len(new)} gf-graded {ELEMENT} I lines had no ENGINE-A row",
        "atomic_data": ("log gf + EP from canonical_gf (the values the measurement used); vdW "
                        "from the covering synthesis list (RYA-410 rule); labels by "
                        "pysme_nlte.labels_by_selection_rules"),
        "existing_rows": "untouched", "script": "scripts/extend_nlte_grid.py"})
    prov_path.write_text(json.dumps(prov, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Extend an element's ENGINE-A (Amarsi 2020 / PySME) table to the lines prepare says it lacks.
SIRIUS ONLY (PySME: /mnt/codex-data/venv_pysme).

    python scripts/extend_nlte_grid.py --element Al --ion I --solar-only   # non-regression + new deltas
    python scripts/extend_nlte_grid.py --element Al --write [--jobs 3]

RYA-1233. The Si extension (scripts/rya1233_extend_si_nlte.py) hard-coded Si and its three
line-set files, so Al would have needed a copy. The line list here is the prepare report's
own `nlte_missing` (data/results/orchestrator/prepare/<star>_<El>.json, check D): every
gf-graded line in the selected ion of the element, not culled for the star, with no row in the table.
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
from pipeline.species import parse_ion, _INT2ROMAN
from pipeline import pysme_nlte as pn                     # noqa: E402
from pipeline.element_prepare import (SYNTH_LISTS, OUT, engine_a_table,  # noqa: E402
                                      engine_a_provenance)

GES = Path(os.environ.get("ISPEC_DIR", "/mnt/codex-data/engines/ispec_src")) / \
    "input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv"
#: The near-UV band's list (RYA-759; generated, not tracked). GES starts at 4200 A, so without
#: this every near-UV line looked unidentified and was held (RYA-1315: 211 Ti II lines).
NEARUV = (3000.0, 3780.0, ROOT / "data/linelists/ispec_nearuv_3000_3780/atomic_lines.tsv")
#: The node set the Si extension used (solar, Teff/logg/[Fe/H] steps, the benchmark stars).
NODES = [(5772, 4.44, 0.0), (5100, 4.44, 0.0), (6200, 4.44, 0.0), (5772, 4.0, 0.0),
         (5772, 4.7, 0.0), (5772, 4.44, -0.3), (5172, 4.43, 0.31), (5100, 4.5, 0.5),
         (6000, 4.3, -0.2), (5772, 4.44, 0.6), (5400, 4.5, 0.5)]
SOLAR = {"teff": 5772, "logg": 4.44, "feh": 0.0, "vmic": 1.0}
_LL: dict = {}
ELEMENT = ""
ION = 1


def _synth_rows(wl: float) -> pd.DataFrame:
    path = NEARUV[2] if NEARUV[0] <= wl < NEARUV[1] else GES
    for lo, f in SYNTH_LISTS:
        if wl >= lo:
            path = f
    key = (path, ELEMENT, ION)
    if key not in _LL:
        d = pd.read_csv(path, sep="\t", low_memory=False,
                        usecols=["element", "wave_A", "lower_state_eV", "waals", "rad", "stark", "upper_state_eV", "lower_j", "upper_j"])
        _LL[key] = d[d.element == f"{ELEMENT} {ION}"]
    return _LL[key]


def atomic(wl: float):
    """(log gf, EP, E_up, vdW) -- gf/EP from canonical_gf, vdW from the synthesis list."""
    cg = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False,
                     usecols=["species", "wavelength_air_A", "excitation_potential_eV", "log_gf"])
    m = cg[(cg.species == f"{ELEMENT} {_INT2ROMAN[ION]}") & ((cg.wavelength_air_A - wl).abs() < 0.01)]
    if len(m) != 1:
        raise ValueError(f"{ELEMENT} {_INT2ROMAN[ION]} {wl}: {len(m)} canonical rows within 0.01 A")
    loggf, ep = float(m.log_gf.iloc[0]), float(m.excitation_potential_eV.iloc[0])
    ll = _synth_rows(wl)
    s = ll[((ll.wave_A - wl).abs() < 0.05) & ((ll.lower_state_eV - ep).abs() < 0.03)]
    if ION != 1:
        if len(s) != 1:
            raise ValueError(f"{ELEMENT} {ION} {wl}: need one identified synthesis component, got {len(s)}")
        row = s.iloc[0]
        values = [row.rad, row.stark, row.waals, row.upper_state_eV, row.lower_j, row.upper_j]
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{ELEMENT} {ION} {wl}: incomplete synthesis atomic data")
        return loggf, ep, float(row.upper_state_eV), float(row.waals), {
            'gamrad': float(row.rad), 'gamqst': float(row.stark)}, float(row.lower_j), float(row.upper_j)
    vw = float(s.waals.iloc[0]) if len(s) else 0.0
    return loggf, ep, ep + 12398.42 / wl, (vw if vw > 0 else 0.0)


#: wavelength -> why its levels could not be resolved (RYA-1315). One unlabelled line used
#: to abort every other line in the batch; it is now held, printed, and written to provenance.
UNLABELLED: dict = {}


def build_lines(waves):
    out = []
    for wl in waves:
        try:
            data = atomic(wl)
            loggf, ep, eup, vw = data[:4]
            tl, tu, jl, ju, note = pn.labels_by_selection_rules(ELEMENT, ep, eup, ion=ION)
        except ValueError as exc:
            UNLABELLED[round(wl, 3)] = str(exc)
            print(f"  {wl:10.3f} HELD: {exc}", flush=True)
            continue
        if ION != 1 and (jl != data[5] or ju != data[6]):
            UNLABELLED[round(wl, 3)] = f"{ELEMENT} {ION} {wl}: grid J values disagree with identified component"
            print(f"  {wl:10.3f} HELD: {UNLABELLED[round(wl, 3)]}", flush=True)
            continue
        line = (wl, loggf, ep, jl, eup, ju, tl, tu, vw)
        out.append(line + (None, data[4]) if ION != 1 else line)
        print(f"  {wl:10.3f} loggf={loggf:+.3f} EP={ep:.3f} {tl} -> {tu} "
              f"(J{jl:.1f}->{ju:.1f}) vdW={vw:.1f}" + (f"  [{note}]" if note else ""), flush=True)
    return out


def derive(star: dict, lines) -> dict:
    """delta per line at one node, derived in groups <= 30 A wide (RYA-773)."""
    out = {}
    lines = sorted(lines, key=lambda l: l[0])
    if not lines:
        return out
    grp = [lines[0]]
    for l in lines[1:]:
        if l[0] - grp[0][0] > 30.0:
            out.update(pn.nlte_delta(ELEMENT, star=star, lines=grp, ion=ION)["per_line"])
            grp = [l]
        else:
            grp.append(l)
    out.update(pn.nlte_delta(ELEMENT, star=star, lines=grp, ion=ION)["per_line"])
    return out


def _node(args):
    (te, lg, fe), lines, el, stage = args
    global ELEMENT, ION
    ELEMENT, ION = el, stage
    from config.constants import SOLAR_ASPLUND2021
    pn._A_SUN.setdefault(el, SOLAR_ASPLUND2021[el])
    pn.NLTE_LINES[el] = list(lines)
    d = derive({"teff": te, "logg": lg, "feh": fe, "vmic": 1.0}, lines)
    print(f"  node {te}/{lg}/{fe:+.2f} done", flush=True)
    return d


def main(argv=None) -> int:
    global ELEMENT, ION
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--element", required=True)
    ap.add_argument("--ion", type=parse_ion, default=1, choices=range(1, 9))
    ap.add_argument("--star", default="solar")
    ap.add_argument("--solar-only", action="store_true")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--jobs", type=int, default=3)
    a = ap.parse_args(argv)
    ELEMENT, ION = a.element, a.ion
    csv_path = engine_a_table(ELEMENT, root=ROOT)
    prov_path = engine_a_provenance(csv_path)
    rep = json.loads((OUT / f"{a.star}_{ELEMENT}.json").read_text())
    from config.constants import SOLAR_ASPLUND2021
    pn._A_SUN.setdefault(ELEMENT, SOLAR_ASPLUND2021[ELEMENT])
    old = pd.read_csv(csv_path)
    banked = sorted(old.loc[old.ion.map(parse_ion) == ION, "wave_A"].unique())
    new = sorted({round(m["wavelength_A"], 3) for m in rep["nlte_missing"]
                  if m["species"] == f"{ELEMENT} {_INT2ROMAN[ION]}"})
    print(f"{len(banked)} banked lines, {len(new)} lines prepare says the table lacks")
    if not new:
        return 0
    new_lines = build_lines(new)

    if a.solar_only:
        print("\nNON-REGRESSION -- the banked lines, re-derived by this atomic-data path:")
        bl = build_lines(banked)
        pn.NLTE_LINES[ELEMENT] = list(bl)
        d = derive(SOLAR, bl)
        sol = old[(old.teff_K == 5772) & (old.logg == 4.44) & (old.feh == 0.0) & (old.ion.map(parse_ion) == ION)]
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
        results = list(ex.map(_node, [(n, new_lines, ELEMENT, ION) for n in NODES]))
    # A line with no finite delta at any node is not written (pysme_nlte returns NaN when the
    # NLTE EW is unbracketed or the LTE curve is non-monotonic); provenance records the hold.
    # The synthesis now confirms departure application explicitly. An exactly zero
    # delta is not evidence that the levels are absent. Preserve the existing no-zero
    # admission policy (RYA-413), recording the observed result without inventing a cause.
    zero = sorted(w for w in results[0]
                  if all(np.isfinite(dd.get(w, np.nan)) and dd[w] == 0.0 for dd in results))
    rows = [dict(element=ELEMENT, ion=ION, wave_A=round(w, 3), teff_K=te, logg=lg, feh=fe,
                 delta_nlte=round(float(d[w]), 4))
            for (te, lg, fe), d in zip(NODES, results) for w in sorted(d)
            if w not in zero and all(np.isfinite(dd.get(w, np.nan)) for dd in results)]
    out = pd.concat([old, pd.DataFrame(rows)], ignore_index=True)
    out.to_csv(csv_path, index=False)
    print(f"wrote {csv_path.name}: {len(old)} kept + {len(rows)} new = {len(out)} rows")
    prov = json.loads(prov_path.read_text()) if prov_path.exists() else {}
    prov.setdefault("extensions", []).append({
        "ticket": os.environ.get("CODEX_TICKET", "RYA-1233"), "date": str(date.today()),
        "ion": ION, "added_waves": sorted({r["wave_A"] for r in rows}),
        "why": f"prepare check D: {len(new)} gf-graded {ELEMENT} {_INT2ROMAN[ION]} lines had no ENGINE-A row",
        "atomic_data": ("log gf + EP from canonical_gf (the values the measurement used); vdW "
                        "from the covering synthesis list (RYA-410 rule); labels by "
                        "pysme_nlte.labels_by_selection_rules"),
        "existing_rows": "untouched", "script": "scripts/extend_nlte_grid.py"})
    nd = {round(w, 3) for w in zero} | {round(w, 3) for w in new
                                        if not any(abs(r["wave_A"] - w) < 0.01 for r in rows)}
    if nd:
        bucket = prov.setdefault("not_derivable_by_ion", {}).setdefault(str(ION), {})
        for w in sorted(nd):
            bucket[f"{w:.3f}"] = (
                "identically 0 at every node despite confirmed departure application; "
                "held under the no-placeholder-zero admission policy; cause unestablished" if w in {round(z, 3) for z in zero} else
                UNLABELLED[w] if w in UNLABELLED else
                "Unbracketed NLTE EW or non-monotonic LTE curve at one or more nodes (no value)")
    prov_path.write_text(json.dumps(prov, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

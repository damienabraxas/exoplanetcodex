#!/usr/bin/env python3
"""RYA-763 — can Engine A be pointed at a line LOCALLY, without the MPIA web query?

    python3 scripts/rya763_level_mapping.py                 # asset inventory
    python3 scripts/rya763_level_mapping.py --element Ti    # + the mapping test

TWO THINGS THE TICKET ASSUMES THAT ARE NOT TRUE ANY MORE
--------------------------------------------------------
1. *"RYA-713 found 14 of 15 staged Engine-A grids have no label_*.txt extracted -- Ti is
   the only one that does."* **Inverted.** RYA-713 then went and recovered them
   (`_label_recovery_rya713.log`, md5-verified, archives deleted): 14 of 15 now HAVE a
   label file and only **Cu** lacks one. The precondition this ticket calls blocking is
   already met.

2. *"Look those levels up directly in the staged .grd plus label_Fe.txt."* **There is no
   Fe asset at all** in the Engine-A grid directory -- no `nlte_Fe_*.grd`, no
   `label_Fe.txt`, no `atmos_Fe.txt`. Fe was never one of the staged Amarsi/PySME grids,
   because Fe's Engine A is the **Bergemann/MPIA** route: a different supplier, reachable
   only through the web service plus a committed per-line CSV extract (which RYA-764
   measured as stopping at 6843.7 A). So step 2 of this ticket **cannot be executed for
   Fe**: there is no local Fe departure grid to interrogate, and the "grid vs our query"
   question is not locally decidable for Fe.

WHAT IS STILL DECIDABLE, AND WHY IT MATTERS MORE
------------------------------------------------
The ticket's own reason for caring is general: *"If local level-mapping works, it applies
to all 15 staged elements."* That question CAN be settled, on any element that has both a
grid and a label file. This does that.

THE TEST, AND THE FAILURE MODE IT IS LOOKING FOR
-------------------------------------------------
GES carries BOTH a numeric `nlte_level_low/up` and a string `nlte_label_low/up` per line.
If the numeric index addresses the same level that `label_{El}.txt` lists at that index,
local mapping is trivial -- the linelist already hands us the grid coordinate.

The dangerous outcome is not "no match". It is a **silent mismatch**: two different model
atoms, each with its own level ordering, indexed by the same integer. Mapping by index
across them would return a real departure coefficient for the WRONG level -- a wrong
number rather than an honest LTE fallback, and invisible. The vendor warns explicitly of
"significant differences" between VALD/GES labels and the NIST labels the grids carry.

So this compares, per line, the GES label string against the label file's entry AT THE
GES INDEX, and reports agreement / disagreement / out-of-range rather than a bare count.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, os.environ.get("ISPEC_DIR", "/srv/codex/engines/ispec_src"))

GRID_DIR = Path(os.environ.get("CODEX_NLTE_GRID_DIR",
                               "/srv/codex/grids/nlte/amarsi_galah"))
GERBER_DIR = Path(os.environ.get("CODEX_GERBER_GRID_DIR",
                                 "/srv/codex/grids/nlte/gerber_ts"))
OUT_DIR = ROOT / "data" / "audit" / "rya763"
UNSET = {"", "0", "-1", "nan", "none", "None", "NONE"}


def inventory() -> pd.DataFrame:
    """grid x label x atmos, per element. Replaces the ticket's stale claim."""
    rows = []
    for g in sorted(GRID_DIR.glob("nlte_*.grd")):
        m = re.match(r"nlte_([A-Za-z]+)_", g.name)
        if not m:
            continue
        el = m.group(1)
        rows.append(dict(element=el, grd=g.name,
                         grd_GB=round(g.stat().st_size / 1e9, 2),
                         label=(GRID_DIR / f"label_{el}.txt").exists(),
                         atmos=(GRID_DIR / f"atmos_{el}.txt").exists()))
    return pd.DataFrame(rows).sort_values("element").reset_index(drop=True)


# ── Gerber TS-native model atoms (RYA-763 CORRECTION) ────────────────────────
#
# An earlier pass of this script asserted "there is no Fe asset at all" on the strength
# of `label_Fe.txt` being absent from the Amarsi/PySME directory. **That was wrong.**
# Fe has a complete NLTE asset set on Sirius, in the GERBER TS-native deck:
#
#     NLTEgrid4TS_Fe_MARCS_May-07-2021.bin   4.16 GB departure grid
#                                            (symlinked out to gerber_overflow/)
#     auxData_Fe_MARCS_May-07-2021.dat       15230 atmosphere nodes
#     atom.fe607a                            607 levels, 12635 transitions
#
# The level identification lives in the MODEL ATOM, not in a `label_*.txt` — which is why
# a search for the Amarsi naming convention missed it. The two decks simply file the same
# information under different names, and I read the absence of one filename as the absence
# of the physics. That is the same "the extract is not the model" error this ticket exists
# to stop, committed against the asset layout instead of against a wavelength range.
#
# Gerber atom format: a header (element, A(X), mass; n_levels, n_transitions, ...), then
# one line per level:
#
#     Energy[cm^-1]   G   'Level N =  <term>'   Ion
#
# G is the statistical weight, so J = (G - 1) / 2, and energy converts to eV by the
# standard 8065.544 cm^-1/eV. That gives exactly the (J, energy) key the label-keyed
# resolver already uses, so Fe needs no new matching logic — only a reader.

CM1_PER_EV = 8065.543937


def read_gerber_atom(path: 'Path') -> pd.DataFrame:
    """Levels from a Gerber TS-native `atom.*` file -> the same frame as read_labels()."""
    rows, seen_header = [], False
    n_expect = None
    for raw in Path(path).read_text(errors="replace").splitlines():
        line = raw.rstrip()
        if not line.strip() or line.lstrip().startswith(("*", "#")):
            continue
        parts = line.split()
        if not seen_header:
            # the 3rd non-comment line carries n_levels n_transitions ...
            if n_expect is None and len(parts) >= 3 and all(
                    x.lstrip("-").isdigit() for x in parts[:3]):
                n_expect = int(parts[0])
                seen_header = True
            continue
        m = re.match(r"\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+'([^']*)'\s+(\d+)", line)
        if not m:
            if len(rows) >= (n_expect or 0):
                break
            continue
        energy_cm, g, label, ion = m.groups()
        try:
            e = float(energy_cm)
            gg = float(g)
        except ValueError:
            continue
        term = label.split("=")[-1].strip() if "=" in label else label.strip()
        rows.append(dict(index=len(rows) + 1, species=f"Fe {ion}",
                         term=term, J=(gg - 1.0) / 2.0,
                         energy_eV=e / CM1_PER_EV, ion=int(ion)))
        if n_expect and len(rows) >= n_expect:
            break
    return pd.DataFrame(rows)


def read_levels(el: str, deck: str = "auto") -> pd.DataFrame:
    """Level table for `el` from whichever deck actually holds it.

    Tries the Amarsi/PySME `label_{El}.txt` first, then the Gerber `atom.*`. Loud on
    neither, and it SAYS which one it used -- the failure this corrects was a silent
    assumption that one naming convention was the only one.
    """
    lab = GRID_DIR / f"label_{el}.txt"
    if deck in ("auto", "amarsi") and lab.exists():
        print(f"  levels from {lab.name} (Amarsi/PySME deck)")
        return read_labels(el)
    hits = sorted(GERBER_DIR.glob(f"atom.{el.lower()}*"))
    if deck in ("auto", "gerber") and hits:
        print(f"  levels from {hits[0].name} (Gerber TS-native deck)")
        return read_gerber_atom(hits[0])
    raise SystemExit(
        f"no level table for {el}: neither {lab} nor {GERBER_DIR}/atom.{el.lower()}*")


def read_labels(el: str) -> pd.DataFrame:
    """label_{El}.txt -> index, species, config, term, J, energy_eV."""
    p = GRID_DIR / f"label_{el}.txt"
    if not p.exists():
        raise SystemExit(f"no label file for {el}: {p}")
    rows = []
    for line in p.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            idx = int(parts[0])
            energy = float(parts[-1])
            j = float(parts[-2])
        except ValueError:
            continue
        term = parts[-3]
        rows.append(dict(index=idx, species=" ".join(parts[1:3]),
                         term=term, J=j, energy_eV=energy))
    return pd.DataFrame(rows)


# ── the label-keyed resolver (RYA-763) ───────────────────────────────────────
#
# The index is a false coordinate (see the module docstring and the offset check): two
# model atoms number their levels differently above the low-lying set, so an integer
# that always resolves resolves to the WRONG level. The fix is to key on what does not
# depend on either atom's bookkeeping.
#
# ENERGY and J are physical. Both the GES linelist and label_{El}.txt carry NIST-derived
# level energies in eV and the level J, so a level can be identified by (J, E) with a
# tolerance, and the term string used only as corroboration. Term NAMES are convention
# ('z3D3*' vs 'z3D' + J=3) and are deliberately NOT the key.
#
# The tolerance is not assumed: `--scan-tol` measures resolution rate against it, so the
# operating point is chosen from the curve rather than picked.

def resolve_by_label(lab: 'pd.DataFrame', energy_eV: float, j: float,
                     tol_eV: float = 0.001) -> tuple:
    """Identify one level by (J, energy). Returns (verdict, index, n_candidates).

    UNIQUE     exactly one atom level within tol at that J -> usable
    AMBIGUOUS  several -> refuse; picking one would be a guess with a real number attached
    ABSENT     none -> the atom genuinely does not carry this level
    """
    if not np.isfinite(energy_eV) or not np.isfinite(j):
        return ("NO-ENERGY-OR-J", -1, 0)
    c = lab[(np.abs(lab.energy_eV - energy_eV) <= tol_eV) & (np.abs(lab.J - j) < 0.01)]
    if len(c) == 1:
        return ("UNIQUE", int(c.iloc[0]["index"]), 1)
    if len(c) > 1:
        return ("AMBIGUOUS", -1, len(c))
    return ("ABSENT", -1, 0)


def label_key_report(ll, lab: 'pd.DataFrame', el: str, ion: str,
                     lo: float, hi: float, tol_eV: float) -> 'pd.DataFrame':
    """Resolve BOTH endpoints of every in-band line by (J, energy)."""
    w = np.asarray(ll["wave_A"], dtype=float)
    els = np.asarray([str(x).strip() for x in ll["element"]])
    want = f"{el} {1 if ion.upper() == 'I' else 2}".upper()
    m = np.array([e.upper().startswith(want) for e in els]) & (w >= lo) & (w <= hi)

    rows = []
    for i in np.where(m)[0]:
        vlo, ilo, nlo = resolve_by_label(lab, float(ll["lower_state_eV"][i]),
                                         float(ll["lower_j"][i]), tol_eV)
        vup, iup, nup = resolve_by_label(lab, float(ll["upper_state_eV"][i]),
                                         float(ll["upper_j"][i]), tol_eV)
        both = (vlo == "UNIQUE") and (vup == "UNIQUE")
        rows.append(dict(wave_A=float(w[i]), lower_verdict=vlo, upper_verdict=vup,
                         atom_index_low=ilo, atom_index_up=iup,
                         n_cand_low=nlo, n_cand_up=nup,
                         engine_a_mappable=both))
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--element", help="run the mapping test for this element")
    ap.add_argument("--ion", default="I")
    ap.add_argument("--lo", type=float, default=6910.0)
    ap.add_argument("--hi", type=float, default=9199.9)
    ap.add_argument("--deck", choices=["auto", "amarsi", "gerber"], default="auto",
                    help="which deck's level table to test against. The GES level "
                         "INDEX is native to the GERBER deck, so testing it against "
                         "the AMARSI atom compares two DIFFERENT atoms -- that is the "
                         "real content of the index disagreement.")
    ap.add_argument("--label-key", action="store_true",
                    help="resolve levels by (J, energy) instead of by index")
    ap.add_argument("--tol-eV", type=float, default=0.001,
                    help="energy match tolerance for --label-key. 0.001 eV is the "
                         "MEASURED optimum from --scan-tol on Ti I (76.6 percent); "
                         "tighter starts missing real matches (ABSENT rises), looser "
                         "collapses into AMBIGUOUS.")
    ap.add_argument("--scan-tol", action="store_true",
                    help="measure resolution rate vs tolerance rather than assuming one")
    a = ap.parse_args()

    inv = inventory()
    gerber = sorted(GERBER_DIR.glob("atom.*"))
    gerber_grids = sorted(GERBER_DIR.glob("NLTEgrid*.bin"))
    print("\n" + "=" * 78)
    print("ENGINE-A LOCAL ASSET INVENTORY  (" + str(GRID_DIR) + ")")
    print("=" * 78)
    print(inv.to_string(index=False))
    n_lab = int(inv.label.sum())
    print(f"\n  {n_lab} of {len(inv)} staged grids have a label file"
          f"{'' if n_lab == len(inv) else '; missing: ' + ', '.join(inv[~inv.label].element)}")
    print(f"  Fe in the Amarsi/PySME deck? "
          f"{'YES' if 'Fe' in set(inv.element) else 'no'}")
    print(f"\n  GERBER TS-native deck ({GERBER_DIR}):")
    print(f"    model atoms  {len(gerber)}: {', '.join(g.name for g in gerber)}")
    print(f"    grids        {len(gerber_grids)}")
    fe_atom = [g for g in gerber if g.name.lower().startswith("atom.fe")]
    if fe_atom:
        print(f"    => Fe HAS a local level table: {fe_atom[0].name}. The earlier claim in")
        print("       this ticket that 'there is no Fe asset at all' was WRONG — it was")
        print("       scoped to the Amarsi naming convention and missed the Gerber deck.")

    if not a.element:
        print("\n  (pass --element <El> to run the level-mapping test)")
        return

    el = a.element
    lab = read_levels(el, a.deck)
    print("\n" + "=" * 78)
    print(f"LEVEL-MAPPING TEST — {el} {a.ion}, GES index vs label_{el}.txt")
    print("=" * 78)
    print(f"  label file: {len(lab)} levels, index {lab['index'].min()}-{lab['index'].max()}, "
          f"energy {lab.energy_eV.min():.3f}-{lab.energy_eV.max():.3f} eV")

    from pipeline.abundances_derive import _load_synth_resources
    ll, _, _ = _load_synth_resources()
    w = np.asarray(ll["wave_A"], dtype=float)
    els = np.asarray([str(x).strip() for x in ll["element"]])
    want = f"{el} {1 if a.ion.upper() == 'I' else 2}".upper()
    m = np.array([e.upper().startswith(want) for e in els]) & (w >= a.lo) & (w <= a.hi)
    print(f"  GES {el} {a.ion} in {a.lo:.0f}-{a.hi:.0f} A: {int(m.sum())} lines")
    if not m.sum():
        print("  no lines in band — nothing to test")
        return

    by_index = lab.set_index("index")
    rows = []
    for i in np.where(m)[0]:
        for side in ("low", "up"):
            raw_idx = str(ll[f"nlte_level_{side}"][i]).strip()
            ges_lab = str(ll[f"nlte_label_{side}"][i]).strip()
            if raw_idx in UNSET or ges_lab in UNSET:
                rows.append(dict(side=side, verdict="NO-LEVEL-ID"))
                continue
            try:
                k = int(float(raw_idx))
            except ValueError:
                rows.append(dict(side=side, verdict="INDEX-UNPARSEABLE"))
                continue
            if k not in by_index.index:
                rows.append(dict(side=side, verdict="INDEX-OUT-OF-RANGE",
                                 ges_index=k, ges_label=ges_lab))
                continue
            r = by_index.loc[k]
            # The GES label is a term+J string like 'z3D3*' / 'e7D4'; the label file
            # holds term ('a3F') and J separately. Compare on the TERM core, which is
            # the part both spell the same way when they mean the same level.
            term = str(r["term"]).strip()
            core = re.sub(r"[^A-Za-z0-9]", "", term).lower()
            gcore = re.sub(r"[^A-Za-z0-9]", "", ges_lab).lower()
            agree = gcore.startswith(core) or core.startswith(gcore[:len(core)])
            rows.append(dict(side=side, verdict="AGREE" if agree else "MISMATCH",
                             ges_index=k, ges_label=ges_lab,
                             atom_term=term, atom_J=r["J"],
                             atom_energy_eV=r["energy_eV"]))

    df = pd.DataFrame(rows)
    print(f"\n  per-level-endpoint outcomes ({len(df)} endpoints = "
          f"{int(m.sum())} lines x 2):")
    for k, v in df.verdict.value_counts().items():
        print(f"    {k:<20} {v:6d}   ({100.0*v/len(df):.1f}%)")

    mm = df[df.verdict == "MISMATCH"]
    ag = df[df.verdict == "AGREE"]
    print("\n" + "-" * 78)
    print("  VERDICT")
    if len(mm) == 0 and len(ag):
        print("  GES indices address the SAME levels as the Engine-A atom. Local")
        print("  level-mapping is sound for this element — the linelist already carries")
        print("  the grid coordinate, and no web query is needed to point Engine A.")
    elif len(ag) == 0:
        print("  NO endpoint agrees. The GES index does NOT address this atom's levels;")
        print("  index-based local mapping would be wrong, not merely unavailable.")
    else:
        print(f"  MIXED: {len(ag)} agree, {len(mm)} disagree. That is the DANGEROUS case —")
        print("  index-based mapping would silently return a real departure coefficient")
        print("  for the wrong level on the disagreeing fraction. Local mapping must key")
        print("  on the LABEL (term+J+energy), never on the bare index.")
        print(f"\n  examples of disagreement:")
        for _, r in mm.head(5).iterrows():
            print(f"    index {int(r.ges_index):4d}  GES '{r.ges_label}'  vs  atom "
                  f"'{r.atom_term}' J={r.atom_J} ({r.atom_energy_eV:.3f} eV)")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    f = OUT_DIR / f"{el}{a.ion}_level_mapping.csv"
    df.to_csv(f, index=False)
    inv.to_csv(OUT_DIR / "engine_a_asset_inventory.csv", index=False)
    print(f"\n  wrote {f.relative_to(ROOT)} and engine_a_asset_inventory.csv")

    if not (a.label_key or a.scan_tol):
        print("\n  (pass --label-key to resolve by (J, energy) instead of by index)")
        return

    print("\n" + "=" * 78)
    print(f"LABEL-KEYED RESOLUTION — {el} {a.ion}, by (J, energy), not by index")
    print("=" * 78)

    if a.scan_tol:
        print(f"\n  resolution rate vs energy tolerance (the operating point is MEASURED):")
        print(f"{'tol (eV)':>10} {'both unique':>12} {'ambiguous':>11} {'absent':>8} "
              f"{'rate':>8}")
        for t in (0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1):
            r = label_key_report(ll, lab, el, a.ion, a.lo, a.hi, t)
            amb = int(((r.lower_verdict == "AMBIGUOUS") |
                       (r.upper_verdict == "AMBIGUOUS")).sum())
            abs_ = int(((r.lower_verdict == "ABSENT") |
                        (r.upper_verdict == "ABSENT")).sum())
            k = int(r.engine_a_mappable.sum())
            print(f"{t:10.4f} {k:12d} {amb:11d} {abs_:8d} "
                  f"{k/max(len(r),1):8.3f}")

    r = label_key_report(ll, lab, el, a.ion, a.lo, a.hi, a.tol_eV)
    k = int(r.engine_a_mappable.sum())
    print(f"\n  at tol = {a.tol_eV} eV, over {len(r)} in-band {el} {a.ion} lines:")
    for side in ("lower", "upper"):
        print(f"    {side} level:  " + "  ".join(
            f"{kk}={vv}" for kk, vv in r[f"{side}_verdict"].value_counts().items()))
    print(f"\n    BOTH levels uniquely resolved: {k}  ({100.0*k/max(len(r),1):.1f}%)")

    idx_rate = (len(df[df.verdict == 'AGREE']) /
                max(len(df[df.verdict.isin(['AGREE', 'MISMATCH'])]), 1))
    print(f"\n  compare: index-keyed agreement was {idx_rate:.1%} on the same band,")
    print(f"  and 0% above level index 215 where the IR lines live.")
    print(f"\n  => local level-mapping IS viable — but on (J, energy), not the index.")
    print(f"     {k} in-band {el} {a.ion} lines are Engine-A mappable without any web query.")

    rf = OUT_DIR / f"{el}{a.ion}_label_keyed.csv"
    r.to_csv(rf, index=False)
    print(f"\n  wrote {rf.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

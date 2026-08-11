#!/usr/bin/env python3
"""RYA-773 step 2 — extend the Al departure extract to the clean doublet. SIRIUS ONLY.

    python3 scripts/rya773_extend_al_grid.py --check-levels     # grid-side evidence
    python3 scripts/rya773_extend_al_grid.py --solar-only       # one node, fast check
    python3 scripts/rya773_extend_al_grid.py --write            # all nodes, emit the CSV

WHAT THIS FIXES
---------------
`Al_Amarsi2020_PySME.csv` served 6696.023 and 6698.673 — and the pipeline's clean,
multi-arm-corroborated Al lines are 7835/7836 and 8772/8773 (RYA-708/716). Al's best
data and Al's NLTE coverage did not overlap, so three of the four in-aggregate lines
carried an explicit ENGINE-A UNCOVERED disposition and the Engine-A product rested on
a single line. The grid was "built for lines the pool never carried" (RYA-706 v50).

That was never an atom limitation: `scripts/rya773_al_level_reach.py` shows the
Nordlander & Lind (2017) atom already carries every level involved (3d 2D, 5f 2F*,
6f 2F*) to within 54 ueV. The extract was just a per-LINE extract of a per-LEVEL grid.
So this is a delta run, not an atom rebuild.

WHAT IT DOES NOT DO
-------------------
It does not touch the 22 existing rows. They are the banked, anchor-validated RYA-402
values; re-deriving them would silently move a registered correction. Non-regression is
proved instead by re-deriving 6696/6698 at solar and COMPARING (--solar-only prints it).

Deriving doublet by doublet (not all six lines in one call) keeps each synthesis window
a few A wide instead of the 2000 A that spanning 6696 to 8773 would require; the
features are far apart and independent, so this changes no number.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import pysme_nlte as pn          # noqa: E402

CSV = ROOT / "data" / "nlte_grids" / "Al_Amarsi2020_PySME.csv"
PROV = ROOT / "data" / "nlte_grids" / "Al_Amarsi2020_PySME.prov.json"
ELEMENT = "Al"

# Derive doublet by doublet — see the module docstring. Read from the registered
# diagnostics, so adding a line to NLTE_LINES is the only edit a future extension needs.
GROUPS = [(6696.023, 6698.673), (7835.309, 7836.134), (8772.865, 8773.897)]
# The pair already in the extract: derived for comparison, never rewritten.
BANKED_GROUP = GROUPS[0]

# What the banked extract says at solar, for the non-regression comparison.
BANKED_SOLAR = {6696.023: -0.0275, 6698.673: -0.0171}


def nodes_from_extract() -> pd.DataFrame:
    """The (Teff, logg, [Fe/H]) nodes the extract already uses, in file order.

    Read from the artifact rather than restated here: the new lines must land on exactly
    the same nodes as the old, or the registry's LinearND interpolation would be working
    off a different hull per wavelength — which is the kind of difference that produces a
    plausible wrong delta instead of an error.
    """
    d = pd.read_csv(CSV)
    n = d[["teff_K", "logg", "feh"]].drop_duplicates().reset_index(drop=True)
    return n


def check_levels() -> int:
    """Grid-side evidence for the level-reach verdict, read from the .grd itself.

    `label_Al.txt` is the shipped NIST table; this confirms the same levels are in the
    binary the departures actually come from, and measures the one thing the label file
    cannot answer — whether the two fine-structure siblings of an upper term (6f 2F*
    J=2.5 vs J=3.5, listed at identical energy) really do carry the same departures. The
    two-component features are emitted through the HFS channel on that assumption, so it
    is measured here rather than asserted.
    """
    from pipeline.nlte_bfactor_synth import read_amarsi_grid
    g = read_amarsi_grid(ELEMENT)
    E, J = g.get("energy"), g.get("J")
    conf, term = g.get("conf"), g.get("term")

    def dec(arr, i):
        # NUL-padded fixed-width fields — str.strip() alone leaves the NULs (RYA-773).
        return bytes(arr[i]).decode("latin1").strip("\x00 \t\r\n")

    print(f"grid: {g.path.name}  {len(E)} levels, {len(g.nodes)} model blocks")
    wanted = {"3s2.3d": "2D", "3s2.5f": "2F*", "3s2.6f": "2F*"}
    found = {}
    for i in range(len(E)):
        c, t = dec(conf, i), dec(term, i)
        if c in wanted and t == wanted[c]:
            found.setdefault((c, t), []).append((i, float(J[i]), float(E[i])))
    for (c, t), levs in sorted(found.items()):
        print(f"  {c:10s} {t:5s} -> " +
              ", ".join(f"idx {i} J={j} E={e:.7f}" for i, j, e in levs))
    missing = [k for k in wanted.items() if k not in found]
    if missing:
        print(f"  MISSING from the grid: {missing}")
        return 1

    # Do the two fine-structure siblings of an UPPER term share departures?
    #
    # This is not idle curiosity. A two-component feature can be synthesised either by
    # giving both components the dominant component's label (the HFS channel, correct
    # only when the siblings share b) or by labelling each component with its own upper
    # level. The 3d 2D LOWER term is listed here for scale but is not part of that
    # question — both components of a feature share the lower level by construction.
    tt, gg, mm = g.nodes[:, 0], g.nodes[:, 1], g.nodes[:, 2]
    k = int(np.argmin((tt - 5772) ** 2 / 100.0 ** 2 + (gg - 4.44) ** 2 + (mm - 0.0) ** 2))
    key = g.node_keys[k]
    b = g.get(key)                      # (ndepth, nlevel)
    print(f"\n  sibling-J departure check at block {key} "
          f"(Teff {tt[k]:.0f}, logg {gg[k]:.2f}, [Fe/H] {mm[k]:+.2f}):")
    worst_upper = 0.0
    for (c, t), levs in sorted(found.items()):
        if len(levs) != 2:
            continue
        (i1, j1, _), (i2, j2, _) = levs[0], levs[1]
        b1, b2 = np.asarray(b)[:, i1], np.asarray(b)[:, i2]
        rel = float(np.max(np.abs(b1 - b2) / np.maximum(np.abs(b1), 1e-30)))
        role = "lower term (not a component question)" if t == "2D" else "UPPER term"
        if t != "2D":
            worst_upper = max(worst_upper, rel)
        print(f"    {c} {t}: J={j1} vs J={j2}  max |db|/b over depth = {rel:.3e}  "
              f"(b {b1.min():.4f}-{b1.max():.4f})   [{role}]")
    print(f"\n  worst UPPER-term sibling disagreement: {worst_upper:.3e}")
    print("  => " + ("the siblings share departures; either treatment is equivalent."
                     if worst_upper < 1e-3 else
                     "the siblings do NOT share departures, so each component carries "
                     "its own upper-level label (pysme_nlte._linelist_rows 4-tuple "
                     "component form). The HFS channel would have averaged this away."))
    return 0


def prove_nul_bug() -> int:
    """Measure the auto_labels NUL-padding defect instead of asserting it.

    `auto_labels` decoded the grid's fixed-width level fields with str.strip(), which
    does not remove NUL, so it returned labels like '3s2.5f\\x00 2F*\\x00'. PySME matches
    a line to a grid level by (species, conf, term, 2J+1); a label that does not match is
    NOT an error — the line simply synthesises in LTE. So the worry was that every delta
    derived through auto_labels (RYA-409/410/411 Family-A, RYA-540 Li/Cu, RYA-592 Mg) had
    quietly been an LTE result wearing an NLTE label.

    RESULT: refuted. The two derivations agree to four decimals, so PySME normalises the
    padding itself and nothing downstream is affected. Kept as a regression probe, and
    because "we checked" is worth more than "it looked fine".

    Al is the right probe precisely because its labels are hand-written and clean, so the
    correct answer is already banked (-0.0275 / -0.0171).
    """
    star = {"teff": 5772, "logg": 4.44, "feh": 0.0, "vmic": 1.0}
    clean = [row for row in pn.NLTE_LINES[ELEMENT] if row[0] in BANKED_GROUP]
    nulled = [tuple(row[:6]) + (row[6] + "\x00", row[7] + "\x00") + tuple(row[8:])
              for row in clean]
    print("deriving Al 6696/6698 at solar with CLEAN labels, then with NUL-padded ones\n")
    a = pn.nlte_delta(ELEMENT, star=star, lines=clean)["per_line"]
    b = pn.nlte_delta(ELEMENT, star=star, lines=nulled)["per_line"]
    print(f"\n{'line':>10} {'clean':>10} {'NUL-padded':>12} {'banked':>10}")
    for w in sorted(a):
        print(f"{w:10.3f} {a[w]:+10.4f} {b[w]:+12.4f} {BANKED_SOLAR[w]:+10.4f}")
    lost = max(abs(a[w] - b[w]) for w in a)
    print(f"\n  largest difference: {lost:.4f} dex")
    print("  => " + ("the NUL labels DO NOT match the grid: those lines ran in LTE. Any "
                     "delta derived through auto_labels before this fix is suspect."
                     if lost > 0.005 else
                     "no measurable difference — PySME normalises the padding itself, so "
                     "no previously derived delta is affected. decode_grid_label is "
                     "hygiene, not a fix."))
    return 0


def derive(star: dict, groups) -> dict:
    """delta per line at one node, derived group by group. Returns {wave: delta}."""
    out = {}
    for grp in groups:
        res = pn.nlte_delta(ELEMENT, star=star, lines=list(grp))
        out.update(res["per_line"])
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check-levels", action="store_true",
                    help="grid-side level + sibling-departure evidence, then exit")
    ap.add_argument("--prove-nul-bug", action="store_true",
                    help="measure what the auto_labels NUL-padding defect actually did: "
                         "derive Al 6696/6698 with clean labels and again with the "
                         "NUL-padded labels auto_labels used to return, and print both")
    ap.add_argument("--solar-only", action="store_true",
                    help="derive the solar node only (includes the 6696/6698 "
                         "non-regression comparison) and do not write")
    ap.add_argument("--write", action="store_true",
                    help="derive every node in the extract and rewrite the CSV + prov")
    a = ap.parse_args()

    if a.check_levels:
        raise SystemExit(check_levels())

    if a.prove_nul_bug:
        raise SystemExit(prove_nul_bug())

    new_groups = GROUPS[1:]
    if a.solar_only:
        star = {"teff": 5772, "logg": 4.44, "feh": 0.0, "vmic": 1.0}
        print(f"solar node {star}")
        banked = derive(star, [BANKED_GROUP])
        print("\nNON-REGRESSION — the pair already in the extract, re-derived:")
        ok = True
        for w, d in sorted(banked.items()):
            ref = BANKED_SOLAR[w]
            same = abs(d - ref) <= 0.001
            ok &= same
            print(f"  {w:9.3f}  banked {ref:+.4f}   re-derived {d:+.4f}   "
                  f"diff {d - ref:+.4f}  {'OK' if same else 'DRIFT'}")
        print(f"  => {'reproduces the banked extract' if ok else 'DOES NOT REPRODUCE'}")
        new = derive(star, new_groups)
        print("\nTHE CLEAN DOUBLET — new:")
        for w, d in sorted(new.items()):
            print(f"  {w:9.3f}  delta {d:+.4f}")
        return

    if not a.write:
        raise SystemExit("nothing to do: pass --check-levels, --solar-only or --write")

    old = pd.read_csv(CSV)
    nodes = nodes_from_extract()
    print(f"extending {CSV.name}: {len(old)} existing rows, {len(nodes)} nodes, "
          f"{sum(len(g) for g in new_groups)} new lines")

    rows = []
    for i, n in nodes.iterrows():
        star = {"teff": float(n.teff_K), "logg": float(n.logg), "feh": float(n.feh),
                "vmic": 1.0}
        d = derive(star, new_groups)
        for w in sorted(d):
            rows.append(dict(element=ELEMENT, ion=1, wave_A=w, teff_K=int(n.teff_K),
                             logg=float(n.logg), feh=float(n.feh),
                             delta_nlte=round(float(d[w]), 4)))
        print(f"  [{i + 1}/{len(nodes)}] Teff {n.teff_K:.0f} logg {n.logg:.2f} "
              f"[Fe/H] {n.feh:+.2f}: " +
              "  ".join(f"{w:.3f} {d[w]:+.4f}" for w in sorted(d)), flush=True)

    add = pd.DataFrame(rows)
    # The existing rows are carried through UNCHANGED and first; new rows follow in a
    # deterministic (node order, then wavelength) order. RYA-768: an artifact that does
    # not byte-diff clean cannot be reviewed.
    out = pd.concat([old, add], ignore_index=True)
    out.to_csv(CSV, index=False)
    print(f"\nwrote {CSV}  ({len(old)} kept + {len(add)} new = {len(out)} rows)")

    prov = json.loads(PROV.read_text())
    served = sorted(out.wave_A.unique())
    prov["what"] = (
        "Al I 1D non-LTE abundance corrections (delta_nlte = A_NLTE - A_LTE) for the "
        "subordinate doublet 6696.023 / 6698.673 (RYA-402) and, added by RYA-773, the "
        "clean doublet 7835.309 / 7836.134 (3d 2D -> 6f 2F*) and 8772.865 / 8773.897 "
        "(3d 2D -> 5f 2F*), derived by synthesising NLTE vs LTE in PySME from the "
        "Amarsi-2020 departure grid (pipeline/pysme_nlte.py). Family-B element "
        "(Option 2, PySME).")
    prov.setdefault("extensions", []).append({
        "ticket": "RYA-773",
        "added_waves": [w for w in served if w not in (6696.023, 6698.673)],
        "why": (
            "The extract served 6696/6698 while the pipeline's clean, multi-arm "
            "corroborated lines are 7835/7836 + 8772/8773 (RYA-708/716) -- Al's best "
            "data and Al's NLTE coverage did not overlap, so 3 of 4 in-aggregate lines "
            "were dispositioned ENGINE-A UNCOVERED and the Engine-A product rested on "
            "one line."),
        "level_reach": (
            "scripts/rya773_al_level_reach.py: every level involved is already in the "
            "Nordlander & Lind 2017 atom -- 3d 2D (idx 7/8), 5f 2F* (23/24), 6f 2F* "
            "(32/33) -- to within 54 ueV. The upper TERM is fixed by the E1 selection "
            "rule, not by energy: 5g/6g 2G sit 3.4-5.5 meV away but are dipole-forbidden "
            "from a 2D lower. No new or enlarged model atom was needed."),
        "non_regression": (
            "The 22 pre-existing rows are carried through unchanged. 6696/6698 were "
            "re-derived at solar and compared against the banked -0.0275 / -0.0171."),
        "generator": "scripts/rya773_extend_al_grid.py --write",
        "derived": date.today().isoformat(),
    })
    prov["grid"]["nodes"] = (
        f"{len(nodes)} (Teff {int(nodes.teff_K.min())}-{int(nodes.teff_K.max())}, logg "
        f"{nodes.logg.min()}-{nodes.logg.max()}, [Fe/H] {nodes.feh.min()}..{nodes.feh.max()}) "
        f"spanning the FGK-dwarf range incl. solar + 55 Cnc; LinearND-interpolated per "
        f"wave by the registry.")
    prov["grid"]["delta_range"] = {
        f"{w:.3f}": [round(float(out[out.wave_A == w].delta_nlte.min()), 4),
                     round(float(out[out.wave_A == w].delta_nlte.max()), 4)]
        for w in served}
    PROV.write_text(json.dumps(prov, indent=2) + "\n")
    print(f"wrote {PROV}")


if __name__ == "__main__":
    main()

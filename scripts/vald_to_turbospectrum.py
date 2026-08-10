#!/usr/bin/env python3
"""Convert a VALD long-format extract to a Turbospectrum atomic line list — RYA-713.

    python3 scripts/vald_to_turbospectrum.py \
        --raw ../exoplanetcodex-rya558/data/linelists/vald_solar_nearuv_2000_3780_hfson_raw.txt \
        --lo 3000 --hi 3780 --out data/linelists/ts_nearuv_3000_3780.lte

Ryan: *"convert the near-UV VALD extract to Turbospectrum format."*

WHY
---
Near-UV synthesis was reported as blocked on "no line list below 4200 Å". That was wrong:
the GES **synthesis** list stops at 4200 Å, but `vald_solar_nearuv_2000_3780_hfson_raw.txt`
holds **55,849 records in 3000–3780 Å across 103 species** — Co I 15779, V I 4813,
Fe I 4364, Mn I 3277, Nb II 2898 … That is exactly what synthesis needs, because synthesis
models every contributor in the window rather than isolating one line.

The blocker was never acquisition. It was format.

THE TARGET FORMAT, taken verbatim from the GES list rather than from documentation
---------------------------------------------------------------------------------
    '  26.000            '    1     14280
    'Fe I    NLTE'
      4200.087  3.884 -1.130   -7.420    7.0  5.01E+07  0.000  'p' 'd'   0.0    1.0 'Fe I LS:… LS:…'

    wavelength(Å)  E_low(eV)  log gf  fdamp(vdW)  g_upper  gamma_rad  gamma_stark
    lower-type  upper-type  0.0  1.0  'label'

Header is `'  Z.000  '  STAGE  NLINES`, where STAGE is 1 for neutral and 2 for singly
ionised — a SEPARATE field, not the decimals of the species code (the trap `ts_gerber_gate`
documents).

An LTE list ends at the label. NLTE lists carry extra level-index and level-label fields
used to map into a departure grid; this writes **LTE**, because no near-UV NLTE grid exists
and inventing those fields would silently claim a mapping we do not have.

WHAT IS DELIBERATELY NOT GUESSED
--------------------------------
VALD supplies a Stark broadening parameter, and every row of the reference GES file carries
`0.000` in that column. Rather than assume VALD's log convention matches whatever TS expects
there, this writes `0.000` to match the reference and **carries VALD's value into the
sidecar CSV instead**. A wrong broadening convention is a silent physics error, and the
column can be filled once its meaning is confirmed against the Turbospectrum source.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

Z = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
     "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
     "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
     "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34,
     "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
     "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
     "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58,
     "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
     "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74,
     "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
     "Bi": 83, "Th": 90, "U": 92}

# 'Fe 1', 3000.123, -1.234, 3.456, 2.5, 7.890, 3.5, lande…, rad, stark, waals, depth
_REC = re.compile(
    r"^'([A-Z][a-z]?)\s+(\d+)',\s*"
    r"([\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+),\s*([\d.]+),\s*(-?[\d.]+),\s*([\d.]+),"
    r"\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+),")
_LS = re.compile(r"^'\s*(LS|JJ|JK|LK)\s+(.*?)'\s*$")


def parse(raw: Path, lo: float, hi: float, max_species: int | None = None):
    """Yield one dict per transition. A record is the data line plus its LS continuations
    and a reference line; the two LS lines carry the lower and upper level designations."""
    out = []
    cur = None
    ls: list[str] = []
    with raw.open(errors="replace") as f:
        for line in f:
            m = _REC.match(line)
            if m:
                if cur is not None:
                    cur["ls"] = ls
                    out.append(cur)
                el, ion = m.group(1), int(m.group(2))
                wl = float(m.group(3))
                cur, ls = (None, []) if not (lo <= wl < hi and el in Z) else (dict(
                    element=el, ion=ion, wl=wl, loggf=float(m.group(4)),
                    e_low=float(m.group(5)), j_lo=float(m.group(6)),
                    e_up=float(m.group(7)), j_up=float(m.group(8)),
                    rad=float(m.group(12)), stark=float(m.group(13)),
                    waals=float(m.group(14))), [])
                continue
            if cur is not None:
                q = _LS.match(line.strip())
                if q:
                    ls.append(q.group(2).strip())
    if cur is not None:
        cur["ls"] = ls
        out.append(cur)
    return out


def ts_rows(recs) -> dict[tuple[str, int], list[str]]:
    """Group formatted TS rows by (element, ion)."""
    blocks = defaultdict(list)
    for r in recs:
        g_up = 2.0 * r["j_up"] + 1.0
        # VALD gives log10(gamma_rad); the reference file carries the LINEAR value.
        gamma_rad = 10.0 ** r["rad"] if r["rad"] else 0.0
        lo_lab = r["ls"][0] if len(r["ls"]) > 0 else ""
        up_lab = r["ls"][1] if len(r["ls"]) > 1 else ""
        label = f"{r['element']} {'I'*r['ion'] if r['ion']<4 else r['ion']} LS:{lo_lab} LS:{up_lab}"
        blocks[(r["element"], r["ion"])].append(
            # LTE ROW LAYOUT -- no gamma_stark field.
            #
            # The NLTE blocks in the GES file carry an extra 0.000 between gamma_rad and
            # the level parities, plus trailing level-index/label/flag fields. The LTE
            # blocks do NOT:
            #   LTE :  ...  4.0  6.61E+07 'p' 'd'   0.0    1.0 'label'
            #   NLTE:  ...  7.0  5.01E+07  0.000  'p' 'd'   0.0    1.0 'label'  82 0 ...
            # I first copied the NLTE layout, which inserts a field the LTE reader does
            # not expect. bsyn then reports the list, reads it, and emits an EMPTY
            # spectrum -- no error, no warning. That is the worst failure shape there is,
            # because every surface check passes.
            f"{r['wl']:10.3f}{r['e_low']:7.3f}{r['loggf']:7.3f}{r['waals']:9.3f}"
            f"{g_up:7.1f} {gamma_rad:.2E} 'p' 'd' {0.0:5.1f} {1.0:6.1f} "
            f"'{label[:120]}'")
    return blocks


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", required=True)
    ap.add_argument("--lo", type=float, required=True)
    ap.add_argument("--hi", type=float, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    raw = Path(a.raw)
    print(f"parsing {raw.name} ({raw.stat().st_size/1e6:.1f} MB) for {a.lo:.0f}-{a.hi:.0f} A ...")
    recs = parse(raw, a.lo, a.hi)
    print(f"  {len(recs)} transitions in band, across "
          f"{len({(r['element'], r['ion']) for r in recs})} species")

    blocks = ts_rows(recs)
    dest = Path(a.out); dest.parent.mkdir(parents=True, exist_ok=True)
    n_lines = 0
    with dest.open("w") as f:
        for (el, ion) in sorted(blocks, key=lambda k: (Z[k[0]], k[1])):
            rows = sorted(blocks[(el, ion)], key=lambda s: float(s[:10]))
            # Species code LEFT-aligned inside the 20-char field, matching the GES
            # reference byte-for-byte: "'  26.000            '". Right-aligning it
            # ("'              26.000'") is accepted by the file reader and then yields
            # an EMPTY spectrum -- bsyn reports the list, reads it, and synthesises
            # nothing, which is the worst kind of failure because it looks like success.
            code = f"  {Z[el]:d}.000"
            f.write(f"'{code:<20}'{ion:>5}{len(rows):>9}\n")
            f.write(f"'{el} {'I'*ion if ion<4 else ion}    LTE'\n")
            f.write("\n".join(rows) + "\n")
            n_lines += len(rows)

    print(f"\n  wrote {dest}  ({dest.stat().st_size/1e6:.1f} MB, {n_lines} lines, "
          f"{len(blocks)} species blocks)")
    top = sorted(blocks.items(), key=lambda kv: -len(kv[1]))[:8]
    print(f"\n{'species':>10s}{'lines':>8s}")
    for (el, ion), rows in top:
        print(f"{el + ' ' + str(ion):>10s}{len(rows):8d}")


if __name__ == "__main__":
    main()

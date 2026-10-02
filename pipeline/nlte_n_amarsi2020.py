"""N I departure legs for the SUN, from Amarsi et al. 2020's own per-line table. RYA-1230.

WHY A TABLE AND NOT A GRID
--------------------------
C I and O I have Amarsi, Nissen & Skuladottir 2019's 4D grids (`pipeline.nlte_cno`). N I
has no such grid: Amarsi, Grevesse, Grumer, Asplund, Barklem & Collet 2020 (A&A 636,
A120) computed the SOLAR N I lines only. Their Table 3 gives, per line, the abundance
from six models on one set of disk-centre equivalent widths -- so the difference of two
columns is a published, per-line, same-EW model correction:

    1D leg  : A(1D non-LTE) - A(1D LTE)        (MARCS)            -0.008 .. -0.011
    3D leg  : A(3D non-LTE) - A(1D LTE)        (STAGGER vs MARCS) -0.037 .. -0.053

The existing N ENGINE-A source (the Amarsi 2020 departure grid via PySME) serves only
8216 and 8683 of the four in-band lines; this table serves all four, and 10108.90.

SOLAR ONLY, BY CONSTRUCTION. These are the Sun's numbers. Any other star is refused: a
per-line solar correction transplanted to another star is the borrowed-number error.

The transcription is RYA-1220's `ni_five_line_reference.csv` (Tables 1-3), read here rather
than re-typed (point at the SSOT, never copy from it) and checked against the paper's
Table 3 text in RYA-1220's `source_text/amarsi2020.txt` by RYA-1230.

⚠️ WHAT THE 3D LEG IS AND IS NOT. It is Amarsi's 3D-NLTE minus Amarsi's 1D-LTE MARCS, added
to OUR 1D-LTE fit. Our 1D model (ATLAS9) differs from their MARCS by ~0.01 dex on these
lines (RYA-1230 ladder, median MARCS step -0.010), which is inside the leg's own spread.
The product says it is a CITED per-line correction, never an in-transfer 3D-NLTE result.
"""
from __future__ import annotations

import csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TABLE = REPO / "data/output/rya1220/audit_v1/ni_five_line_reference.csv"
CITATION = ("Amarsi, Grevesse, Grumer, Asplund, Barklem & Collet 2020, A&A 636, A120, "
            "Table 3 (solar, disk-centre EWs)")
#: Half AGSS21/Amarsi's printed resolution (nm to 2 dp), the RYA-1214 set tolerance.
MATCH_TOL_A = 0.05

_COLUMNS = {"1D": ("A_1d_nlte", "A_1d_lte"), "3D": ("A_3d_nlte", "A_1d_lte")}


class NotSolarError(ValueError):
    """The Amarsi 2020 N I table is solar-only."""


def _rows() -> list[dict]:
    with TABLE.open() as fh:
        return list(csv.DictReader(fh))


def source(leg: str) -> str:
    what = {"1D": "1D non-LTE - 1D LTE (MARCS)",
            "3D": "3D non-LTE (STAGGER) - 1D LTE (MARCS)"}[leg]
    return (f"{CITATION} -- {what}; SOLAR per-line additive correction, same EWs in both "
            f"columns (match <= {MATCH_TOL_A} A)")


def deltas(wavelengths_A, leg: str, *, star: str) -> dict[float, float]:
    """Per-line correction keyed by the caller's own wavelength. Unmatched lines are absent."""
    if star != "solar":
        raise NotSolarError(f"Amarsi 2020 Table 3 is the Sun's; refusing to apply it to {star!r}")
    hi, lo = _COLUMNS[leg]
    table = [(float(r["wavelength_air_A"]), float(r[hi]) - float(r[lo])) for r in _rows()]
    out: dict[float, float] = {}
    for w in wavelengths_A:
        w = float(w)
        best = min(table, key=lambda t: abs(t[0] - w))
        if abs(best[0] - w) <= MATCH_TOL_A:
            out[w] = round(best[1], 4)
    return out

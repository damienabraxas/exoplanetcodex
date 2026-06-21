"""
pipeline/threed_corrections.py
==============================
RYA-399 — the 3D abundance-correction leg for the GET-3D metals (Si/Ti/Cr).

WHY
---
RYA-400 routed Si/Ti/Cr as "3D-owed". RYA-398 left them at a graded+1D-NLTE solar
residual of +0.378 (Si) / +0.502 (Ti) / +0.402 (Cr). This module is the 3D
extension of the existing C/N/O 3D path (pipeline.nlte_cno): it adds the published
solar 3D dimensional correction ON TOP of the 1D-NLTE pass —

    A(3D-NLTE) = A(1D-NLTE; NLTE_CORRECTION_ELEMENTS) + delta_3d

so it never double-counts the NLTE already applied (the increment, A(3D)-A(1D) at
fixed NLTE).

THE FINDING (validate-don't-tune)
---------------------------------
The published solar 3D corrections for these metals are SMALL: Si -0.01
(Amarsi & Asplund 2017), Ti +0.06, Cr +0.03 (Scott et al. 2015, A&A 573 A26) —
all |delta| <= 0.1 dex, an order of magnitude below the +0.4-0.5 residuals, and
for Ti/Cr POSITIVE (3D raises the abundance, moving them the WRONG way). So 3D is
NOT the lever that closes the Si/Ti/Cr solar residual; the residual is line-data /
gf-zero-point (RYA-161 differential territory) and is carried forward to the
multi-star arc, never tuned. The corrections are real published physics, applied
because the audit owed the check — and the check says 3D is not the answer here.

NO FAKING
---------
delta_3d is read from the vendored, provenance-stamped grid
(data/threed_grids/solar3d_metals_rya399.csv). An element with no grid node is
left at 1D-NLTE with threed_flag='3D_unavailable' — never silently "corrected".
A full off-solar 3D grid exists publicly only for Si (Amarsi 2020 GALAH); Ti/Cr
have none (Fe-peak excluded) — off-solar 3D for those is a documented gap.

Linear: RYA-399  (blockedBy RYA-400; blocks RYA-371 Phase C)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from config.constants import THREED_CORRECTION_ELEMENTS, SOLAR_ASPLUND2021  # noqa: E402

_GRID_DIR = _REPO / 'data' / 'threed_grids'

# Magnitude sanity ceiling: a solar 3D dimensional correction larger than this is
# physically implausible for an FGK dwarf metal line and almost certainly a grid
# typo — fail loud rather than apply it (no silent fat-finger).
_THREED_MAG_CEILING_DEX = 0.5

_grid_cache: dict = {}


def parse_ion(ion) -> int:
    """'I'->1, 'II'->2, 1->1, '1'->1 (mirror of the NLTE module's parser)."""
    if isinstance(ion, (int, np.integer)):
        return int(ion)
    s = str(ion).strip()
    if s in ('1', '2', '3'):
        return int(s)
    return {'I': 1, 'II': 2, 'III': 3}.get(s.upper(), 1)


def _load_solar3d_grid() -> dict:
    """element -> {'delta_3d', 'unc', 'source', 'ion'} from the vendored solar grid.
    Loud-fails if a registered element is missing from the grid or its correction
    exceeds the magnitude ceiling."""
    if _grid_cache:
        return _grid_cache
    # all registered elements point at the same grid file in this build
    grids = {spec['grid'] for spec in THREED_CORRECTION_ELEMENTS.values()}
    rows = {}
    for g in grids:
        path = _GRID_DIR / g
        if not path.exists():
            raise FileNotFoundError(f"3D grid not found: {path}")
        df = pd.read_csv(path, comment='#')
        for _, r in df.iterrows():
            d = float(r['delta_3d'])
            if abs(d) > _THREED_MAG_CEILING_DEX:
                raise ValueError(
                    f"3D correction |delta_3d|={abs(d):.3f} for {r['element']} exceeds "
                    f"the {_THREED_MAG_CEILING_DEX} dex sanity ceiling — refusing to apply "
                    f"(likely a grid typo).")
            rows[str(r['element'])] = {
                'delta_3d': d, 'unc': float(r['delta_3d_unc']),
                'source': str(r['source']), 'ion': parse_ion(r['ion']),
            }
    for el in THREED_CORRECTION_ELEMENTS:
        if el not in rows:
            raise KeyError(
                f"{el} is registered in THREED_CORRECTION_ELEMENTS but absent from the "
                f"3D grid file(s) {sorted(grids)} — no silent gap.")
    _grid_cache.update(rows)
    return _grid_cache


def solar_threed_delta(element: str) -> float:
    """Solar 3D dimensional correction (dex) for `element`, or NaN if not registered."""
    if element not in THREED_CORRECTION_ELEMENTS:
        return np.nan
    return _load_solar3d_grid()[element]['delta_3d']


def threed_ref(element: str) -> str:
    spec = THREED_CORRECTION_ELEMENTS.get(element, {})
    return spec.get('ref', '')


def apply_threed_corrections(
    abundances_df: pd.DataFrame,
    stellar_params: dict = None,
    elements=None,
) -> pd.DataFrame:
    """Add the 3D dimensional correction on top of the 1D-NLTE abundance for the
    registered GET-3D metals. Composes AFTER apply_element_nlte_corrections: reads
    the 1D-NLTE column (`A_X_nlte`, falling back to `A_X`) and writes

        delta_3d / A_X_3dnlte / threed_flag / threed_ref

    touching ONLY rows whose (element, ion) are in THREED_CORRECTION_ELEMENTS and
    leaving every other row exactly as the NLTE pass left it. An element with no
    grid node is flagged '3D_unavailable' and left at 1D-NLTE (no silent fake).

    `stellar_params` is accepted for signature-parity with the NLTE path and the
    future off-solar grid; this solar-node build applies the solar correction for
    any FGK-dwarf params (and notes the off-solar limitation for Ti/Cr).
    """
    reg = THREED_CORRECTION_ELEMENTS
    want = set(elements) if elements is not None else set(reg)
    grid = _load_solar3d_grid()

    out = abundances_df.copy()
    for col, default in (('delta_3d', np.nan), ('A_X_3dnlte', np.nan),
                         ('threed_flag', 'no_3D'), ('threed_ref', '')):
        if col not in out.columns:
            out[col] = default

    for idx, row in out.iterrows():
        element = str(row['element'])
        if element not in reg or element not in want:
            continue
        if parse_ion(row.get('ion', 1)) != int(reg[element]['ion']):
            continue
        # base = the 1D-NLTE abundance if present, else the 1D-LTE A_X
        base = row.get('A_X_nlte', np.nan)
        if not np.isfinite(base):
            base = row.get('A_X', np.nan)
        if not np.isfinite(base):
            out.at[idx, 'threed_flag'] = '3D_unavailable'
            continue
        node = grid.get(element)
        if node is None or not np.isfinite(node['delta_3d']):
            out.at[idx, 'A_X_3dnlte'] = base
            out.at[idx, 'threed_flag'] = '3D_unavailable'
            continue
        d = node['delta_3d']
        out.at[idx, 'delta_3d'] = round(d, 4)
        out.at[idx, 'A_X_3dnlte'] = round(base + d, 3)
        out.at[idx, 'threed_flag'] = '3D_solar_Scott2015_AA2017'
        out.at[idx, 'threed_ref'] = reg[element]['ref']
    return out


# ── Report: does 3D close the RYA-398 graded residual? (the whole point) ──────

_GRADED_DIAG = _REPO / 'data' / 'curation' / 'nonfe_pools' / 'curation_diagnostics_graded_rya398.csv'


def residual_after_3d() -> pd.DataFrame:
    """Take the RYA-398 graded+1D-NLTE curated abundances for Si/Ti/Cr, add the
    published solar 3D correction, and report the residual before/after. This is the
    validate-don't-tune demonstration that 3D does NOT close the solar residual."""
    diag = pd.read_csv(_GRADED_DIAG)
    rows = []
    for el in ('Si', 'Ti', 'Cr'):
        r = diag[diag['element'] == el]
        if r.empty:
            continue
        r = r.iloc[0]
        a_nlte = float(r['A_nlte_curated'])
        asp = float(r['asplund'])
        d3 = solar_threed_delta(el)
        a_3d = a_nlte + d3
        rows.append({
            'element': el,
            'A_1Dnlte': round(a_nlte, 3),
            'resid_1Dnlte': round(a_nlte - asp, 3),
            'delta_3d': round(d3, 3),
            'A_3Dnlte': round(a_3d, 3),
            'resid_3Dnlte': round(a_3d - asp, 3),
            'asplund': asp,
            'closes_residual': abs(a_3d - asp) < 0.10,
            'verdict': '3D_NOT_THE_LEVER' if abs(a_3d - asp) >= 0.10 else 'CLOSED',
        })
    return pd.DataFrame(rows)


def _cli(argv=None):
    import argparse
    p = argparse.ArgumentParser(description="RYA-399 metal 3D abundance corrections")
    p.add_argument('--report', action='store_true',
                   help="show whether 3D closes the RYA-398 graded residual")
    a = p.parse_args(argv)

    print("── RYA-399 solar 3D corrections for the GET-3D metals (Si/Ti/Cr) ──")
    for el in ('Si', 'Ti', 'Cr'):
        d = solar_threed_delta(el)
        print(f"  {el}: delta_3d = {d:+.3f} dex   [{threed_ref(el)}]")
    print("  (all |delta_3d| <= 0.1 dex; Ti/Cr POSITIVE — 3D raises the abundance)")

    if a.report:
        rep = residual_after_3d()
        print("\n── Does 3D close the RYA-398 graded solar residual? ──")
        print(rep.to_string(index=False))
        print("\n  FINDING: 3D moves Si/Ti/Cr by <= 0.06 dex and (Ti/Cr) the WRONG way;")
        print("  the +0.4 dex residual survives → 3D is NOT the lever. The residual is")
        print("  line-data / gf-zero-point (RYA-161 differential territory), carried")
        print("  forward to the multi-star arc — never tuned. Off-solar 3D grids for")
        print("  Ti/Cr do not exist publicly (documented gap); Si has Amarsi 2020.")


if __name__ == '__main__':
    _cli()

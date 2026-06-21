"""
pipeline/pysme_nlte.py
======================
RYA-402 Option 2 — derive NLTE abundance corrections from the Amarsi PySME
departure grids by synthesising NLTE vs LTE in PySME (the native consumer of the
`.grd` files), for the Family-B elements (Al, K, Cu, S) whose grids exist only in
SME form (no Turbospectrum model atom).

VALIDATED against the Na anchor (RYA-402 Step 3, the guard): PySME reproduces the
INSPECT Lind-2011 solar Na correction to within tolerance —

    Na I 5682.633  delta -0.121
    Na I 5688.205  delta -0.138
    median        -0.129   vs anchor -0.107 +/- 0.03   -> PASS

so the machinery (grid read + line<->level match + NLTE synthesis + delta
extraction, with the correct sign) is trusted. This is validate-don't-tune: the
anchor is REPRODUCED, never fitted.

Method
------
For each diagnostic line at fixed stellar params: synthesise the line in NLTE (the
grid's departures applied) -> EW_NLTE, and on an LTE curve of growth EW_LTE(A);
the NLTE correction is delta = A(NLTE) - A(LTE) for a fixed EW, i.e.
  delta = A_used - A*,   where EW_LTE(A*) = EW_NLTE   (note the sign).

Gotchas baked in (each was a real failure mode while validating):
  * NLTE needs the VALD3 LONG line format (short silently runs LTE).
  * lines match grid levels by (species, configuration, term, 2J+1) — the long
    line's term_lower/upper must read 'conf term' exactly as the grid encodes them.
  * PySME's NLTE LFS resolves the grid via a file URI, which breaks on a path with
    spaces — we symlink the grid into a space-free temp dir.
  * the ABO van der Waals broadening (gamvw, e.g. 1955.327) matters: these lines
    are saturated and delta is damping-sensitive.

Requires `pysme-astro` (pip). Heavy (synthesis) + needs the MARCS atmosphere grid
(auto-downloaded once) — so this is a manual/offline derivation tool, not a CI test.
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

import numpy as np

# Diagnostic lines per element, with the grid level labels for NLTE matching.
# (wl_A, loggf, Elow_eV, Jlow, Eup_eV, Jup, term_lower, term_upper, gamvw_ABO)
# Na is the VALIDATION anchor (Lind 2011 INSPECT, solar delta -0.107).
NLTE_LINES = {
    'Na': [
        (5682.633, -0.706, 2.102, 0.5, 4.284, 1.5, '2p6.3p 2P*', '2p6.4d 2D', 1955.327),
        (5688.205, -0.404, 2.104, 1.5, 4.283, 2.5, '2p6.3p 2P*', '2p6.4d 2D', 1955.327),
    ],
    # Al I subordinate doublet 4s 2S -> 5p 2P* (the clean Na-analog: subordinate,
    # weak, not the saturated resonance regime). gamvw=0 -> Unsold (GES waals=0).
    'Al': [
        (6696.023, -1.569, 3.143, 0.5, 4.994, 1.5, '3s2.4s 2S', '3s2.5p 2P*', 0.0),
        (6698.673, -1.870, 3.143, 0.5, 4.993, 0.5, '3s2.4s 2S', '3s2.5p 2P*', 0.0),
    ],
    # K / Cu / S diagnostic lines + their grid level labels are added as each
    # element is derived (read the labels from PySMEGrid.get('conf'/'term'/'J'/'energy')).
}

# Solar A(X) reference (Asplund 2021) for the COG zero point.
_A_SUN = {'Na': 6.24, 'Al': 6.43, 'K': 5.07, 'Cu': 4.18, 'S': 7.12}

_GRID_FILENAME = {
    'Na': 'nlte_Na_scatt_pysme.grd', 'Al': 'nlte_Al_scatt_pysme.grd',
    'K': 'nlte_K_scatt_pysme.grd', 'Cu': 'nlte_Cu_caliskan_Oct2024_pysme.grd',
    'S': 'nlte_S_ama51_Sep2024_pysme.grd',
}
_REPO = Path(__file__).resolve().parents[1]
_GRID_DIR = _REPO / 'data' / 'nlte_grids' / 'amarsi_galah'


def _spacefree_grid(element: str) -> str:
    """PySME resolves the NLTE grid via a file URI, which breaks on a path with
    spaces. Symlink the vendored grid into a space-free temp dir and return that."""
    src = _GRID_DIR / _GRID_FILENAME[element]
    if not src.exists():
        raise FileNotFoundError(f"PySME grid for {element} not intaken: {src}")
    d = Path(tempfile.gettempdir()) / 'rya402_pysme_grids'
    d.mkdir(exist_ok=True)
    link = d / _GRID_FILENAME[element]
    if not link.exists():
        os.symlink(src, link)
    return str(link)


def _synth_ew(element, offset, nlte, star, lines, grid_path):
    """One PySME synthesis; returns {wl: EW_mA}. Lazy PySME import (fail-loud)."""
    import pandas as pd
    from pysme.sme import SME_Structure
    from pysme.abund import Abund
    from pysme.linelist.linelist import LineList
    from pysme.synthesize import synthesize_spectrum

    rows = []
    for wl, gf, elo, jlo, eup, jup, tl, tu, vw in lines:
        rows.append(dict(species=f'{element} 1', wlcent=wl, excit=elo, gflog=gf,
                         gamrad=7.8, gamqst=0.0, gamvw=vw, atom_number=_Z(element),
                         ionization=1, lande_lower=0.0, lande_upper=0.0, lande=0.0,
                         j_lo=jlo, j_up=jup, e_upp=eup, term_lower=tl, term_upper=tu,
                         error=0.0, depth=0.6, reference='RYA-402'))
    sme = SME_Structure()
    sme.teff, sme.logg, sme.monh = star['teff'], star['logg'], star['feh']
    sme.vmic, sme.vmac, sme.vsini = star.get('vmic', 1.0), 0.0, 0.0
    ab = Abund.solar(); ab[element] = _A_SUN[element] + offset; sme.abund = ab
    sme.linelist = LineList(pd.DataFrame(rows), lineformat='long')
    wmin = min(l[0] for l in lines) - 2.0
    wmax = max(l[0] for l in lines) + 2.0
    sme.wave = np.linspace(wmin, wmax, int((wmax - wmin) * 220))
    sme.atmo.source = 'marcs2012.sav'; sme.atmo.method = 'grid'; sme.atmo.geom = 'PP'
    if nlte:
        sme.nlte.set_nlte(element, grid_path)
    sme = synthesize_spectrum(sme)
    w = np.asarray(sme.wave[0]); f = np.asarray(sme.synth[0])

    def ew(c, hw=0.8):
        m = (w > c - hw) & (w < c + hw)
        return float(np.trapz(1 - f[m], w[m]) * 1000.0)
    return {l[0]: ew(l[0]) for l in lines}


def _Z(el):
    return {'Na': 11, 'Al': 13, 'K': 19, 'Cu': 29, 'S': 16}[el]


def nlte_delta(element: str, star: dict = None, offs=(-0.2, -0.1, 0.0, 0.1, 0.2)) -> dict:
    """Per-line NLTE abundance correction delta = A(NLTE) - A(LTE) via PySME, plus
    the median. Raises if the element has no registered diagnostic lines/grid."""
    if element not in NLTE_LINES:
        raise KeyError(f"No NLTE diagnostic lines registered for {element} "
                       f"(have {list(NLTE_LINES)}). Add them from the grid level labels.")
    star = star or {'teff': 5772, 'logg': 4.44, 'feh': 0.0, 'vmic': 1.0}
    lines = NLTE_LINES[element]
    grid = _spacefree_grid(element)
    ew_nlte = _synth_ew(element, 0.0, True, star, lines, grid)
    cog = {l[0]: [] for l in lines}
    for off in offs:
        ew_lte = _synth_ew(element, off, False, star, lines, grid)
        for l in lines:
            cog[l[0]].append(ew_lte[l[0]])
    A = _A_SUN[element] + np.array(offs)
    per_line = {}
    for l in lines:
        wl = l[0]
        a_star = float(np.interp(ew_nlte[wl], np.array(cog[wl]), A))
        per_line[wl] = _A_SUN[element] - a_star      # A(NLTE)=A_sun minus A_LTE(=a_star)
    return {'element': element, 'per_line': per_line,
            'delta_median': float(np.median(list(per_line.values())))}


# Published solar anchors for the Step-3 guard — each element validated against its
# OWN published delta (reproduce, never fit). (anchor, tol, reference).
#   Na: Lind 2011 INSPECT solar median -0.107 (the strong subordinate doublet).
#   Al: Nordlander & Lind 2017 (A&A 607, A75) — subordinate optical lines have
#       "small negative abundance corrections of at most -0.04 dex on the lower MS";
#       so the solar 6696/6698 correction is near-zero-to-slightly-negative. Anchor
#       the band there (NOT the +0.2 resonance 3961 line).
_ANCHOR = {
    'Na': (-0.107, 0.03, 'Lind et al. 2011 INSPECT (Na I 5682/5688)'),
    'Al': (-0.02, 0.04, 'Nordlander & Lind 2017, A&A 607 A75 (subordinate Al I; <=-0.04 dex on the lower MS)'),
}


def validate(element: str = 'Na') -> dict:
    """The RYA-402 Step-3 guard, generalised per element: reproduce the element's OWN
    published solar delta. `passed` True iff within tolerance. NEVER tunes the anchor —
    if it cannot be reproduced, `passed` is False and the caller STOPS."""
    if element not in _ANCHOR:
        raise KeyError(f"No published anchor for {element}; have {list(_ANCHOR)}")
    res = nlte_delta(element)
    anchor, tol, ref = _ANCHOR[element]
    res['anchor'] = anchor
    res['anchor_tol'] = tol
    res['passed'] = abs(res['delta_median'] - anchor) <= tol
    res['ref'] = ref
    return res


def validate_na() -> dict:
    """Back-compat alias — the Na guard."""
    return validate('Na')


if __name__ == '__main__':
    import json, sys
    el = sys.argv[1] if len(sys.argv) > 1 else 'Na'
    print(json.dumps(validate(el), indent=2, default=float))

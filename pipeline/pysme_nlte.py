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
    # Cu I 5782 (3d9.4s2 2D -> 3d10.4p 2P*), the standard NLTE-studied line, used to
    # validate the machinery + the HFS x NLTE interaction. RYA-402 finding: Cu NLTE is
    # SMALL (single-component delta +0.001; HFS-resolved 10-component +0.003 -> HFS x
    # NLTE interaction ~+0.002 dex, NEGLIGIBLE for this weak line; would only bite a
    # saturated HFS line). Consistent with Shi et al. 2014 (small positive optical Cu).
    # NOT registered for production: our MEASURED Cu lines are inadequate (4767 EW=198
    # mA at EP 5.7 is unphysical -> junk; 4704 upper is a high-Rydberg state with no
    # clean grid level) and the standard lines (5105/5220/5782) are unmeasured. Cu's
    # blocker is line quality (RYA-395 curation), not the NLTE machinery.
    'Cu': [
        (5782.130, -1.488, 1.642, 1.5, 3.786, 0.5, '3d9.4s2 2D', '3d10.4p 2P*', -7.79),
    ],
    # S I multiplet-8 high-excitation subordinate lines 4p 5P -> 5d 5D* (6748/6757),
    # the standard optical S diagnostics (we measure them, clean). Small negative NLTE
    # (high-excitation -> both levels ~thermalised): solar median -0.016. gamvw=0 Unsold.
    'S': [
        (6748.680, -0.639, 7.868, 2.0, 9.704, 3.0, '3s2.3p3.(4S*).4p 5P', '3s2.3p3.(4S*).5d 5D*', 0.0),
        (6757.150, -0.240, 7.870, 3.0, 9.704, 4.0, '3s2.3p3.(4S*).4p 5P', '3s2.3p3.(4S*).5d 5D*', 0.0),
    ],
    # K I resonance doublet 7665/7699 (4s 2S -> 4p 2P*, GROUND state lower) -- the boss
    # fight: SATURATED + LARGE negative NLTE (line much stronger in NLTE) + sits in the
    # O2 telluric A-band. ABO van der Waals (487/485). Needs the wide-window/bracket
    # derivation (_DERIV_OPTS['K']). Solar: 7665 -0.27, 7699 -0.31, median -0.29.
    'K': [
        (7664.899, 0.149, 0.0, 0.5, 1.617, 1.5, '3p6.4s 2S', '3p6.4p 2P*', 487.23),
        (7698.964, -0.154, 0.0, 0.5, 1.610, 0.5, '3p6.4s 2S', '3p6.4p 2P*', 485.23),
    ],
}

# Solar A(X) reference (Asplund 2021) for the COG zero point.
_A_SUN = {'Na': 6.24, 'Al': 6.43, 'K': 5.07, 'Cu': 4.18, 'S': 7.12}

_GRID_FILENAME = {
    'Na': 'nlte_Na_scatt_pysme.grd', 'Al': 'nlte_Al_scatt_pysme.grd',
    'K': 'nlte_K_scatt_pysme.grd', 'Cu': 'nlte_Cu_caliskan_Oct2024_pysme.grd',
    'S': 'nlte_S_ama51_Sep2024_pysme.grd',
    # RYA-409 Part B re-source (v3 Amarsi-2020 grids, [Fe/H] -> +1):
    'Mg': 'nlte_Mg_scatt_pysme.grd', 'Si': 'nlte_Si_scatt_pysme.grd',
    'Ca': 'nlte_Ca_scatt_pysme.grd', 'Mn': 'nlte_Mn_scatt_pysme.grd',
}
_REPO = Path(__file__).resolve().parents[1]
_GRID_DIR = _REPO / 'data' / 'nlte_grids' / 'amarsi_galah'


def auto_labels(element: str, elow_eV: float, eup_eV: float, tol: float = 0.06):
    """Resolve the grid level labels for a line by NEAREST ENERGY — so Family-A
    re-derivations (Mg/Si/Ca/Mn) don't need hand-mapping like Al/S did. Returns
    (term_lower 'conf term', term_upper 'conf term', j_lo, j_up). Raises if either
    level is >tol eV from any grid level (a wrong-level match is exactly where a
    silent NLTE error hides). RYA-409."""
    from pipeline.nlte_bfactor_synth import read_amarsi_grid
    g = read_amarsi_grid(element)
    E = g.get('energy'); J = g.get('J')
    conf = g.get('conf'); term = g.get('term')

    def _dec(arr, i):
        return bytes(arr[i]).decode('latin1').strip()

    def _match(e):
        i = int(np.argmin(np.abs(E - e)))
        if abs(float(E[i]) - e) > tol:
            raise ValueError(f"{element}: no grid level within {tol} eV of {e:.3f} eV "
                             f"(nearest {float(E[i]):.3f}).")
        return i

    il, iu = _match(elow_eV), _match(eup_eV)
    return (f"{_dec(conf, il)} {_dec(term, il)}", f"{_dec(conf, iu)} {_dec(term, iu)}",
            float(J[il]), float(J[iu]))


def _spacefree_grid(element: str) -> str:
    """PySME resolves the NLTE grid via a file URI, which breaks on a path with
    spaces. Symlink the vendored grid into a space-free temp dir and return that."""
    src = _GRID_DIR / _GRID_FILENAME[element]
    if not src.exists():
        raise FileNotFoundError(f"PySME grid for {element} not intaken: {src}")
    d = Path(tempfile.gettempdir()) / 'rya402_pysme_grids'
    d.mkdir(exist_ok=True)
    link = d / _GRID_FILENAME[element]
    # Replace a stale/broken link (e.g. left by a previous run whose grid was since
    # freed to manage disk) — os.path.lexists catches a broken symlink that
    # link.exists() reports absent, which would make os.symlink raise FileExists.
    if os.path.lexists(link):
        if not link.exists() or os.path.realpath(link) != str(src.resolve()):
            os.unlink(link)
    if not os.path.lexists(link):
        os.symlink(src, link)
    return str(link)


# Per-element derivation options. Saturated resonance lines (K 7665/7699) need a
# WIDE EW window (broad damping wings) and a WIDE abundance bracket (the COG is flat,
# so EW barely moves with A) — the standard narrow settings under-integrate the wing
# and fail to bracket. (ew_hw in A, offs = abundance offsets for the LTE COG.)
_DERIV_OPTS = {
    'K': {'ew_hw': 2.0, 'offs': (-0.4, -0.2, 0.0, 0.2, 0.4)},
}
_DEFAULT_OPTS = {'ew_hw': 0.8, 'offs': (-0.2, -0.1, 0.0, 0.1, 0.2)}


def _linelist_rows(element, lines):
    """Build PySME long-format line rows from diagnostic lines. HFS-resolved (RYA-411):
    a line's optional 10th element is a list of (wl, gflog) hyperfine components — each is
    emitted as its own row sharing the feature's lower/upper NLTE level labels (HFS splits
    the level by ~ueV, so the departure coefficients are identical). Pure logic (no PySME),
    so the HFS expansion is unit-testable without the gitignored grid."""
    rows = []
    for line in lines:
        wl, gf, elo, jlo, eup, jup, tl, tu, vw = line[:9]
        comps = line[9] if len(line) > 9 and line[9] else [(wl, gf)]
        for cwl, cgf in comps:
            rows.append(dict(species=f'{element} 1', wlcent=cwl, excit=elo, gflog=cgf,
                             gamrad=7.8, gamqst=0.0, gamvw=vw, atom_number=_Z(element),
                             ionization=1, lande_lower=0.0, lande_upper=0.0, lande=0.0,
                             j_lo=jlo, j_up=jup, e_upp=eup, term_lower=tl, term_upper=tu,
                             error=0.0, depth=0.6, reference='RYA-411'))
    return rows


def _synth_ew(element, offset, nlte, star, lines, grid_path, ew_hw=0.8):
    """One PySME synthesis; returns {feature_wl: EW_mA}. Lazy PySME import (fail-loud).

    HFS-resolved (RYA-411): a diagnostic 'line' may carry an optional 10th element — a
    list of (wl, gflog) hyperfine components. When present, each component is emitted as
    its OWN line (sharing the feature's lower/upper NLTE level labels, EP, J, vdW — HFS
    splits the level by ~ueV so the departure coefficients are identical), so the feature
    DESATURATES correctly instead of being a single over-saturated gf-summed line (which
    suppressed the Mn NLTE delta in RYA-410). EW is still integrated per FEATURE center."""
    import pandas as pd
    from pysme.sme import SME_Structure
    from pysme.abund import Abund
    from pysme.linelist.linelist import LineList
    from pysme.synthesize import synthesize_spectrum

    rows = _linelist_rows(element, lines)
    sme = SME_Structure()
    sme.teff, sme.logg, sme.monh = star['teff'], star['logg'], star['feh']
    sme.vmic, sme.vmac, sme.vsini = star.get('vmic', 1.0), 0.0, 0.0
    ab = Abund.solar(); ab[element] = _A_SUN[element] + offset; sme.abund = ab
    # PySME requires the line list ascending in wavelength.
    sme.linelist = LineList(pd.DataFrame(rows).sort_values('wlcent').reset_index(drop=True),
                            lineformat='long')
    wmin = min(l[0] for l in lines) - 2.0
    wmax = max(l[0] for l in lines) + 2.0
    sme.wave = np.linspace(wmin, wmax, int((wmax - wmin) * 220))
    sme.atmo.source = 'marcs2012.sav'; sme.atmo.method = 'grid'; sme.atmo.geom = 'PP'
    if nlte:
        sme.nlte.set_nlte(element, grid_path)
    sme = synthesize_spectrum(sme)
    w = np.asarray(sme.wave[0]); f = np.asarray(sme.synth[0])

    def ew(c, hw=ew_hw):
        m = (w > c - hw) & (w < c + hw)
        return float(np.trapz(1 - f[m], w[m]) * 1000.0)
    return {l[0]: ew(l[0]) for l in lines}


def _Z(el):
    return {'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'S': 16, 'K': 19,
            'Ca': 20, 'Mn': 25, 'Cu': 29}[el]


def nlte_delta(element: str, star: dict = None, offs=None) -> dict:
    """Per-line NLTE abundance correction delta = A(NLTE) - A(LTE) via PySME, plus
    the median. Uses per-element derivation options (_DERIV_OPTS) — wide EW window +
    bracket for saturated lines (K). Raises if the element has no diagnostic lines."""
    if element not in NLTE_LINES:
        raise KeyError(f"No NLTE diagnostic lines registered for {element} "
                       f"(have {list(NLTE_LINES)}). Add them from the grid level labels.")
    opts = _DERIV_OPTS.get(element, _DEFAULT_OPTS)
    offs = offs if offs is not None else opts['offs']
    ew_hw = opts['ew_hw']
    star = star or {'teff': 5772, 'logg': 4.44, 'feh': 0.0, 'vmic': 1.0}
    lines = NLTE_LINES[element]
    grid = _spacefree_grid(element)
    ew_nlte = _synth_ew(element, 0.0, True, star, lines, grid, ew_hw=ew_hw)
    cog = {l[0]: [] for l in lines}
    for off in offs:
        ew_lte = _synth_ew(element, off, False, star, lines, grid, ew_hw=ew_hw)
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
    'Cu': (0.01, 0.05, 'Shi et al. 2014 (small positive ~+0.02 for optical Cu I in the Sun; approximate band)'),
    'S':  (-0.04, 0.07, 'Amarsi et al. 2025 (A&A 703 A35, grid source) / Takeda 2005: optical high-excitation S I (6757 mult-8) small negative, 0 to ~-0.1'),
    'K':  (-0.27, 0.10, 'Reggiani et al. 2019 (A&A 627 A177) / Andrievsky et al. 2006: K I 7665/7699 resonance, severe negative NLTE (line stronger in NLTE); solar ~-0.2..-0.3'),
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

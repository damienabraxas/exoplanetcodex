"""
pipeline/cno_synthesis.py
=========================
Region-aware C / N / O synthesis engine — Turbospectrum flux-fit (RYA-237).

C/N/O is the flagship science (the C/O ratio drives the rocky-planet-composition
thesis and the 55 Cnc controversy). C, N and O are coupled through molecular
equilibrium (CO / CN / CH), so we **synthesize** rather than invert EWs:
Turbospectrum's equation of state solves the molecular partial pressures
internally and self-consistently at each set of A(C)/A(N)/A(O). There is no
hand-coded "CO correction" — the iteration is just re-fit + re-synthesize until
A(C)/A(N)/A(O) are self-consistent (the EOS recomputes equilibrium each call).

This module builds the **region-aware engine** and validates it on the
**HARPS-VIS arm** (the data we have now: solar Dumusque HARPS + Procyon HARPS).
The engine is built region-aware — per-instrument LSF, a telluric-correction
gate for IR arms, and per-region product output — so the UV (STIS) and IR
(CO-band) arms plug in without rework as their data clears the campaign gates
(RYA-351 / RYA-162 / RYA-119). Those arms are sequenced by data readiness, not
deferred.

Plugs into (already on main):
  * synth-v2 flux-fitting core — `abundances_derive._synth_flux_at_abund` /
    `_fit_synth_flux` (RYA-285/287); we reuse its resources + broadening.
  * single-source gf via `gf_resolver` (RYA-353) — the GES linelist loaded by
    `_load_synth_resources()` is already rescaled to canonical gf.
  * per-star broadening (RYA-288) — `_resolve_broadening`, fail-loud (no solar
    default leak into Procyon).
  * the sign-corrected NLTE module (RYA-339) — wired as a pluggable per-arm hook.

Molecular bands: iSpec's Turbospectrum wrapper auto-includes
`input/linelists/turbospectrum/molecules/*.bsyn` (CH Masseron, CN Brooke+Sneden,
C2, OH, NH, CO; RYA-236, verified at the iSpec tool path) whenever
`use_molecules=True`. Setting C/N/O as fixed abundances feeds babsma's molecular
equilibrium, so the band depths respond to A(C)/A(N)/A(O) — that coupling is the
whole point of synthesizing CNO.

NLTE — VIS arm is LTE-correct by physics (run scope, Ryan 2026-06-19; this is
correct treatment, NOT a silent fallback — see `VIS_NLTE_POLICY`):
  * [O I] 6300 — forbidden, LTE-insensitive.
  * CH G-band / C2 Swan / CN red — molecular bands (no molecular NLTE grid).
  * C I 5052/5380 — optical C I, ΔNLTE < 0.03 dex (Alexeeva & Mashonkina 2015);
    stamp `cI_vis_lte_assumed` so it is revisited when RYA-359 lands and the
    red-arm O I 777 cross-check runs.
The 339+359 post-hoc / in-synthesis NLTE path is exercised when the
red-optical / IR / UV arms land (O I 777 triplet, IR/FUV C I) — that is where
corrections are large and the grid is mandatory. The NLTE application is kept
pluggable per arm (`nlte_backend`) so those arms slot in without rework.

Linear issue: RYA-237   (Procyon-VIS shakedown: RYA-348)
"""

import os
import sys
import json
import time
import argparse
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.interpolate import interp1d

from config.constants import (
    PATHS, SOLAR_ASPLUND2021, HARPS_R, ISPEC_DIR, get_star_params,
)
# Reuse the synth-v2 core (RYA-285/287/288/353) — same atmosphere interpolation,
# linelist (canonical-gf rescaled), isotopes, observed-spectrum loader and
# per-star broadening resolver the Fe synthesis path uses. No parallel machinery.
from pipeline.abundances_derive import (
    ispec,
    _load_atmosphere,
    _load_synth_resources,
    _load_observed_spectrum,
    _resolve_broadening,
    _ISPEC_SOLAR_ABUND_FILE,
)
from pipeline.gf_resolver import resolve as resolve_gf   # RYA-365: canonical Ni gf assert
from pipeline import nlte_cno   # RYA-359: vendored Amarsi 2019 C I / O I 3D-NLTE grid (Phase-A C I correction)

# Molecular line lists (RYA-236) — iSpec globs this dir when use_molecules=True.
_MOLECULES_DIR = ISPEC_DIR / 'input' / 'linelists' / 'turbospectrum' / 'molecules'

# χ²ᵣ fit-quality gate per band (the established synth convention, RYA-342).
from config.constants import SYNTH_CHI2_GATE

_C_KMS = 299792.458


# ── Region awareness ──────────────────────────────────────────────────────────
# Built now, lit per arm as data clears (RYA-351/162/119). Each arm carries its
# instrument LSF, a telluric-clearance requirement, and the NLTE backend the arm
# routes through. Products are emitted PER REGION (never coadded across
# resolutions — combined at the abundance layer, RYA-282 presentation decision).

@dataclass(frozen=True)
class RegionConfig:
    name: str                       # 'vis' | 'red' | 'ir' | 'uv'
    instrument: str                 # 'HARPS' | 'UVES' | 'CRIRES+' | 'STIS' ...
    R: float                        # instrumental resolving power (LSF)
    wave_min_A: float
    wave_max_A: float
    telluric_correction_required: bool
    nlte_backend: str               # 'lte_by_design' (VIS) | 'amarsi_grid' | 'in_synthesis'
    notes: str = ''


HARPS_VIS = RegionConfig(
    name='vis', instrument='HARPS', R=float(HARPS_R),
    wave_min_A=3780.0, wave_max_A=6910.0,
    telluric_correction_required=False,   # optical HARPS; the chosen VIS windows
                                          # are free of significant telluric bands
    nlte_backend='lte_by_design',         # VIS CNO lines are LTE-correct (see module docstring)
    notes='HARPS 3780-6910 Angstrom; solar Dumusque + Procyon. RYA-237 deliverable.',
)

REGIONS = {'vis': HARPS_VIS}


# ── Diagnostic registry (HARPS-VIS, wavelength-correct) ───────────────────────
# Windows are AIR wavelengths in Angstrom, the fit sub-windows for each band.
# `depends_on` lists elements whose abundance must be fixed before this diagnostic
# is fit (the equilibrium coupling — CN needs A(C); [O I] is tied to A(C) via CO).

@dataclass(frozen=True)
class Diagnostic:
    key: str
    element: str                    # 'C' | 'N' | 'O'
    kind: str                       # 'molecular_band' | 'atomic' | 'forbidden_blend'
    windows_A: tuple                # tuple of (lo, hi) Angstrom fit sub-windows
    use_molecules: bool
    role: str                       # 'primary' | 'cross_check'
    nlte_flag: str
    nlte_ref: str
    depends_on: tuple = ()
    pinned_blends: tuple = ()        # elements pinned from EW/canonical (e.g. Ni for [O I])
    reference: str = ''


VIS_DIAGNOSTICS = (
    # ── Carbon ────────────────────────────────────────────────────────────────
    Diagnostic(
        key='CH_Gband', element='C', kind='molecular_band',
        windows_A=((4303.5, 4306.5), (4310.0, 4313.0)),
        use_molecules=True, role='primary',
        nlte_flag='lte_molecular_band',
        nlte_ref='molecular band — no NLTE grid (LTE-by-design)',
        reference='CH A-X G-band 4290-4315 (Masseron+2014/2022); primary solar/Procyon A(C)',
    ),
    Diagnostic(
        key='CI_5052', element='C', kind='atomic',
        windows_A=((5051.3, 5053.0),),
        use_molecules=False, role='cross_check',
        nlte_flag='cI_vis_lte_assumed',
        nlte_ref='Alexeeva & Mashonkina 2015: optical C I ΔNLTE < 0.03 dex',
        reference='C I 5052.17 atomic; LTE cross-check vs CH (revisit when RYA-359 lands)',
    ),
    Diagnostic(
        key='CI_5380', element='C', kind='atomic',
        windows_A=((5379.3, 5381.3),),
        use_molecules=False, role='cross_check',
        nlte_flag='cI_vis_lte_assumed',
        nlte_ref='Alexeeva & Mashonkina 2015: optical C I ΔNLTE < 0.03 dex',
        reference='C I 5380.34 atomic; LTE cross-check vs CH',
    ),
    Diagnostic(
        key='C2_Swan', element='C', kind='molecular_band',
        windows_A=((5160.0, 5166.0),),
        use_molecules=True, role='cross_check',
        nlte_flag='lte_molecular_band',
        nlte_ref='molecular band — no NLTE grid (LTE-by-design)',
        reference='C2 Swan (0,0) bandhead 5165 (metal-rich cross-check, 55 Cnc)',
    ),
    # ── Nitrogen ──────────────────────────────────────────────────────────────
    Diagnostic(
        key='CN_red', element='N', kind='molecular_band',
        windows_A=((6125.0, 6130.0), (6195.0, 6200.0)),
        use_molecules=True, role='primary', depends_on=('C',),
        nlte_flag='lte_molecular_band',
        nlte_ref='molecular band — no NLTE grid (LTE-by-design)',
        reference='CN red A-X 6000-6200 (Brooke+Sneden 2014); primary N, given A(C). '
                  'Also deblends Li 6707 (RYA-103).',
    ),
    # ── Oxygen ────────────────────────────────────────────────────────────────
    Diagnostic(
        key='OI_6300', element='O', kind='forbidden_blend',
        windows_A=((6299.5, 6301.0),),
        use_molecules=True, role='primary', depends_on=('C',),
        pinned_blends=('Ni',),
        nlte_flag='lte_forbidden_insensitive',
        nlte_ref='[O I] forbidden — LTE-insensitive',
        reference='[O I] 6300.30 + Ni I 6300.34 joint synthesis (A(Ni) pinned). '
                  'Ni I 6300.34 gf resolves via gf_resolver = Johansson+2003 '
                  'log gf -2.11 (RYA-365 adjudication; the [O I]-blend lab authority).',
    ),
)


# ── VIS NLTE policy — LTE-by-design, pluggable per arm ────────────────────────
# This is the `lte_by_design` backend. It applies ZERO correction and stamps the
# physics-justified flag per diagnostic. The red/IR/UV arms pass a different
# backend (Amarsi post-hoc grid RYA-359, or in-synthesis departures RYA-361/363)
# with the SAME signature: (diagnostic, A_lte, params) -> (A_nlte, delta, flag, ref).

def vis_lte_backend(diag: 'Diagnostic', a_lte: float, params: dict) -> tuple:
    """LTE-by-design VIS NLTE backend. delta = 0; flag justified per diagnostic.

    Correct treatment, NOT a silent fallback (Ryan handoff 2026-06-19): the VIS
    CNO lines are LTE or near-LTE by their physics. Does NOT loud-fail on an
    absent C/O grid — the grid (RYA-359) is mandatory only for the red/IR/UV arms.
    """
    return float(a_lte), 0.0, diag.nlte_flag, diag.nlte_ref


def amarsi_grid_backend(diag: 'Diagnostic', a_lte: float, params: dict) -> tuple:
    """Placeholder for the red/IR-arm Amarsi C I / O I post-hoc grid (RYA-359).

    Loud-fails until the grid lands — the red-arm O I 777 / FUV C I corrections
    are large and negative; routing them through LTE silently is exactly the
    failure this guards. Wired here so the per-arm NLTE path is real and pluggable;
    NOT used by the VIS deliverable.
    """
    raise NotImplementedError(
        f"Amarsi C/O NLTE grid backend (RYA-359) not yet available for diagnostic "
        f"'{diag.key}'. The red-optical/IR/UV arms require it (O I 777, FUV C I — "
        f"large negative corrections). The VIS arm is LTE-by-design (use "
        f"vis_lte_backend); do not route VIS through here."
    )


NLTE_BACKENDS = {'lte_by_design': vis_lte_backend, 'amarsi_grid': amarsi_grid_backend}


# ── Phase-A cited correction layer (RYA-371) ──────────────────────────────────
# The VIS synthesis MEASURES in 1D-LTE; Phase A applies the CITED 3D/NLTE
# correction on top, per diagnostic, each tagged with its source. VALIDATE-DON'T-
# TUNE: every value here is published / vendored, NEVER fitted to the Asplund
# anchors. Three kinds, by what the literature actually provides for the line:
#   * vendored grid delta — C I 5052/5380: Amarsi 2019 3D-NLTE (pipeline.nlte_cno,
#     RYA-359). A real interpolated Delta = A(3D-NLTE) - A(1D-LTE).
#   * cited 3D anchor     — [O I] 6300: Caffau et al. 2015 (A&A 579 A88, POSP III)
#     full-3D solar A(O)=8.73 with OUR EXACT atomic data (Ni I -2.11 Johansson
#     2003 + [O I] -9.717 Storey & Zeippen 2000; RYA-367). NO hardcoded 3D-1D
#     offset — Caffau 2015 publishes the absolute, not a grid node (RYA-367 rule).
#   * 3D-offset-owed      — CH/CN/C2 molecular bands: no vendored solar 3D grid →
#     reported 1D-LTE, flagged owed (honest, not silently called LTE).

# Cited per-line full-3D solar anchors: {key: (A_3d, unc, flag, source)}. The
# published absolute for the exact line + our atomic data, surfaced as the
# reconciled value — a citation, not a fit.
CITED_3D_ANCHORS = {
    'OI_6300': (
        8.73, 0.05, '3d_lte_caffau2015',
        'Caffau et al. 2015, A&A 579 A88 (POSP III): [O I] 630 nm CO5BOLD full-3D, '
        'Ni I -2.11 (Johansson 2003) + [O I] -9.717 (Storey & Zeippen 2000) = our '
        'atomic data; in gate 8.69+/-0.05. RYA-367 (no hardcoded 3D-1D offset).'),
}
# Atomic C I lines that carry a vendored Amarsi-2019 3D-NLTE grid delta.
_CI_GRID_KEYS = {'CI_5052': 5052.17, 'CI_5380': 5380.34}


def apply_cited_corrections(per_band, params, region) -> list:
    """Attach the cited Phase-A correction to each VIS diagnostic. Returns a list of
    records {key, element, role, a_lte, kind, a_corr, delta, flag, source}. Cited /
    vendored values only — never fitted (RYA-371 validate-don't-tune)."""
    teff = float(params['teff_K']); logg = float(params['logg'])
    feh = float(params['feh']); vmic = float(params['vturb_kms'])
    out = []
    for r in per_band:
        key = r.get('key'); a_lte = r.get('A_X')
        rec = {'key': key, 'element': r.get('element'), 'role': r.get('role'),
               'a_lte': a_lte, 'kind': None, 'a_corr': a_lte, 'delta': 0.0,
               'flag': r.get('nlte_flag'), 'source': r.get('nlte_ref')}
        if not (a_lte is not None and np.isfinite(a_lte)):
            out.append(rec); continue
        if key in CITED_3D_ANCHORS:                       # cited full-3D anchor
            a3d, unc, flag, src = CITED_3D_ANCHORS[key]
            rec.update(kind='cited_3d_anchor', a_corr=a3d, delta=round(a3d - a_lte, 3),
                       flag=flag, source=src, unc=unc)
        elif key in _CI_GRID_KEYS:                        # vendored Amarsi-2019 grid delta
            try:
                label = nlte_cno.resolve_line('CI', _CI_GRID_KEYS[key])
                delta = nlte_cno.cno_nlte_delta('CI', label, teff, logg, feh, vmic, a_lte)
                if np.isfinite(delta):
                    nlte_cno.assert_cno_sign('CI', label, delta)
                    rec.update(kind='amarsi2019_grid', a_corr=round(a_lte + delta, 3),
                               delta=round(delta, 3), flag='3d_nlte_amarsi2019',
                               source=f'Amarsi, Nissen & Skuladottir 2019 A&A 630 A104, '
                                      f'C I {label} 3D-NLTE leg ({nlte_cno.select_leg(teff)})')
                else:                                     # outside 4D hull → flag, no silent LTE
                    rec.update(kind='grid_out_of_hull',
                               flag='amarsi2019_out_of_hull', source='Amarsi 2019 grid: query outside 4D hull')
            except Exception as exc:                       # noqa: BLE001 — surface, never fake
                rec.update(kind='grid_error', flag='amarsi2019_error', source=f'grid error: {exc}')
        elif r.get('nlte_flag') == 'lte_molecular_band':  # molecular band, no vendored 3D grid
            rec.update(kind='3d_offset_owed', flag='lte_molecular_band_3d_offset_owed',
                       source='molecular band; no vendored solar 3D-LTE offset grid → '
                              'reported 1D-LTE, 3D offset OWED (Asplund 2005b CH/C2 3D-1D '
                              '0.00..-0.15; not applied — would be uncited)')
        out.append(rec)
    return out


# ── Abundance state + low-level synthesis ─────────────────────────────────────

_ISPEC_SCALE_OFFSET = 12.036   # A(X) = log(N/Ntot)_iSpec + 12.036 (matches abundances_derive)


def _atom_codes(elements, chem_elements, solar_abund) -> dict:
    """{element: iSpec atom code} via create_free_abundances_structure."""
    codes = {}
    for el in elements:
        s = ispec.create_free_abundances_structure([el], chem_elements, solar_abund)
        codes[el] = int(s['code'][0])
    return codes


def _solar_A(elements, chem_elements, solar_abund) -> dict:
    """{element: A(X) on the iSpec internal (Asplund 2009) scale}."""
    out = {}
    for el in elements:
        s = ispec.create_free_abundances_structure([el], chem_elements, solar_abund)
        out[el] = float(s['Abund'][0]) + _ISPEC_SCALE_OFFSET
    return out


def _fixed_ab(state: dict, codes: dict) -> np.recarray:
    """Build a multi-element fixed-abundance recarray on the iSpec SPECTRUM scale.

    state: {element: A(X)} for every element to pin (C, N, O, Ni). All listed
    elements are passed to Turbospectrum as fixed abundances; the abundance under
    fit is just the one we vary between synth calls. babsma uses these to solve the
    molecular equilibrium, so CH/CN/CO band strengths track the set C/N/O.
    """
    elems = list(state.keys())
    fa = np.recarray(len(elems), dtype=[('code', int), ('Abund', float),
                                        ('element', '|U30')])
    for i, el in enumerate(elems):
        fa['code'][i]    = codes[el]
        fa['Abund'][i]   = float(state[el]) - _ISPEC_SCALE_OFFSET
        fa['element'][i] = el
    return fa


def _synth_window(sw_nm, atm, params, ll, iso, sab, fixed_ab,
                  broadening, use_molecules, tmp_dir) -> np.ndarray:
    """Normalized Turbospectrum flux over sw_nm at the given fixed composition.

    Asserts code='turbospectrum' (no SPECTRUM fallback) and threads the
    region-aware instrumental R + per-star vmac/vsini broadening through.
    """
    R, vmac, vsini = broadening
    return ispec.generate_spectrum(
        sw_nm, atm,
        float(params['teff_K']), float(params['logg']), float(params['feh']), 0.0,
        ll, iso, sab, fixed_ab,
        microturbulence_vel=float(params['vturb_kms']),
        macroturbulence=vmac, vsini=vsini, R=R,
        verbose=0, code='turbospectrum',
        use_molecules=use_molecules, tmp_dir=tmp_dir,
    )


# ── Window flux fit (single free element) ─────────────────────────────────────
# Fit A(free_el) by minimizing reduced χ² between observed normalized flux and
# the broadened synthetic flux over the diagnostic's sub-windows. All other CNO
# elements (and pinned blends) are held fixed at `state` — the EOS recomputes
# molecular equilibrium each eval. Single free parameter (the abundance), per the
# RYA-287 convention. σ is the constant 0.01-flux model-adequacy floor (RYA-287):
# its scale shifts χ²ᵣ magnitude but NOT the χ² minimum location (the fitted A).

_SIGMA_FLUX = 0.01
_WSTEP_NM = 0.0002          # fine synthesis grid (0.002 Angstrom)


def _fit_element(obs_w_nm, obs_f, atm, params, free_el, state, codes,
                 windows_A, use_molecules, broadening, a_lo, a_hi,
                 ll, iso, sab, tmp_dir) -> dict:
    """Minimize χ²ᵣ over A(free_el) across the diagnostic windows.

    Returns {A_X, red_chi2, sigma_fit, n_pix, n_eval, status}. No silent fallback:
    a synthesis error → status='failed', A_X=nan (never substituted).
    """
    # Pre-slice observed pixels per window (rest-frame), build synthesis grids.
    segs = []
    for (lo_A, hi_A) in windows_A:
        lo_nm, hi_nm = lo_A / 10.0, hi_A / 10.0
        m = (obs_w_nm >= lo_nm) & (obs_w_nm <= hi_nm)
        if m.sum() < 5:
            continue
        sw = np.arange(lo_nm, hi_nm + _WSTEP_NM * 0.5, _WSTEP_NM)
        segs.append((obs_w_nm[m], obs_f[m], sw))
    if not segs:
        return {'A_X': np.nan, 'red_chi2': np.nan, 'sigma_fit': np.nan,
                'n_pix': 0, 'n_eval': 0, 'status': 'failed',
                'reason': 'no observed pixels in windows'}

    n_pix = int(sum(len(ow) for ow, _, _ in segs))
    fail = [None]
    n_eval = [0]
    local = dict(state)

    def chi2(a_x):
        n_eval[0] += 1
        local[free_el] = float(a_x)
        fa = _fixed_ab(local, codes)
        tot = 0.0
        for ow, of, sw in segs:
            try:
                sf = _synth_window(sw, atm, params, ll, iso, sab, fa,
                                   broadening, use_molecules, tmp_dir)
            except Exception as exc:          # noqa: BLE001 — surface, never fake
                fail[0] = str(exc)
                return 1e30
            sf_i = interp1d(sw, sf, bounds_error=False, fill_value=1.0)(ow)
            r = (of - sf_i) / _SIGMA_FLUX
            tot += float(np.nansum(r * r))
        return tot

    res = minimize_scalar(chi2, bounds=(a_lo, a_hi), method='bounded',
                          options={'xatol': 1e-3})
    if fail[0] is not None:
        return {'A_X': np.nan, 'red_chi2': np.nan, 'sigma_fit': np.nan,
                'n_pix': n_pix, 'n_eval': int(n_eval[0]), 'status': 'failed',
                'reason': f'synthesis error: {fail[0]}'}

    a_best = float(res.x)
    dof = max(n_pix - 1, 1)
    chi2_min = chi2(a_best)
    red_chi2 = chi2_min / dof
    # 1σ fit uncertainty from the χ² parabola: Δχ² = 1 above the minimum. Probe a
    # small step and invert the local curvature (σ_A = step / sqrt(Δχ²_step)).
    step = 0.05
    a_hi_probe = min(a_best + step, a_hi)
    dchi = max(chi2(a_hi_probe) - chi2_min, 1e-6)
    sigma_fit = float(abs(a_hi_probe - a_best) / np.sqrt(dchi))
    sigma_fit = float(np.clip(sigma_fit, 0.0, 1.0))
    edge = min(abs(a_best - a_lo), abs(a_best - a_hi)) < 1e-2
    return {'A_X': round(a_best, 3), 'red_chi2': round(float(red_chi2), 3),
            'sigma_fit': round(sigma_fit, 3), 'n_pix': n_pix,
            'n_eval': int(n_eval[0]),
            'status': 'edge_pinned' if edge else 'ok'}


# ── Preflight assertions (no silent fallback) ─────────────────────────────────

def _molecules_cover(windows_A) -> bool:
    """True iff at least one .bsyn molecular file spans each requested window.

    Filenames encode nm ranges, e.g. 12CH_400-450.bsyn. iSpec globs them at synth
    time; we fail loud here if a molecular band is requested but uncovered.
    """
    import re
    files = list(_MOLECULES_DIR.glob('*.bsyn'))
    spans = []
    for f in files:
        m = re.match(r'.*_(\d+)-(\d+)\.bsyn', f.name)
        if m:
            spans.append((float(m.group(1)), float(m.group(2))))
    for (lo_A, hi_A) in windows_A:
        lo_nm, hi_nm = lo_A / 10.0, hi_A / 10.0
        if not any(s_lo <= lo_nm and hi_nm <= s_hi for s_lo, s_hi in spans):
            return False
    return True


def preflight(region: RegionConfig, star_id: str, diagnostics) -> dict:
    """Assert the CRITICAL no-silent-fallback invariants; return broadening tuple.

    - Turbospectrum only (RT code is hardcoded to 'turbospectrum' in _synth_window).
    - gf via gf_resolver (canonical table present; the GES linelist is rescaled in
      _load_synth_resources — RYA-353).
    - per-star broadening (RYA-288, fail-loud — never a solar default for Procyon).
    - molecular lists cover every molecular band (RYA-236, at the iSpec path).
    - telluric clearance for arms that require it (IR; VIS does not).
    """
    print(f"  [preflight] region={region.name} ({region.instrument}, R={region.R:.0f})")

    # gf_resolver / canonical gf table
    canon = PATHS.get('canonical_gf') if hasattr(PATHS, 'get') else None
    canon = canon or (Path(PATHS['linelist_solar']).parent / 'canonical_gf.csv')
    if not Path(str(canon)).exists():
        raise FileNotFoundError(
            f"Canonical gf table not found ({canon}); the synth path must resolve gf "
            f"via gf_resolver (RYA-353). Refusing to run with unresolved/duplicated gf.")
    print(f"  [preflight] gf_resolver canonical table present: {Path(str(canon)).name}")

    # per-star broadening (fail-loud)
    R_star, vmac, vsini, fit_vmac = _resolve_broadening(star_id)
    if fit_vmac:
        raise NotImplementedError(
            f"[{star_id}] vmac='fit' (RYA-309 §3.3) — RT vmac fit not wired in CNO "
            f"engine; refusing to use the init guess as a fixed vmac (RYA-288).")
    broadening = (float(region.R), float(vmac), float(vsini))   # region LSF + per-star
    print(f"  [preflight] broadening (per-star, RYA-288): R={region.R:.0f}, "
          f"vmac={vmac} km/s, vsini={vsini} km/s")

    # molecular coverage for molecular bands
    for d in diagnostics:
        if d.use_molecules and not _molecules_cover(d.windows_A):
            raise FileNotFoundError(
                f"No molecular .bsyn list covers {d.key} windows {d.windows_A} at "
                f"{_MOLECULES_DIR}. RYA-236 lists must be present (verify at iSpec "
                f"tool path; secured by RYA-360). Refusing to synthesize the band "
                f"without its molecular opacity.")
    print(f"  [preflight] molecular lists cover all molecular bands "
          f"({_MOLECULES_DIR.name}/*.bsyn)")

    # telluric clearance gate (real per-arm hook; VIS not required)
    if region.telluric_correction_required:
        raise RuntimeError(
            f"Region {region.name} requires telluric-corrected input (IR arm: "
            f"cr2res+molecfit / APERO+Wapiti, RYA-351). No clearance flag supplied — "
            f"refusing to fit CNO over uncorrected telluric bands.")
    print(f"  [preflight] telluric gate: not required for {region.name} (optical)")

    # NLTE backend resolves
    if region.nlte_backend not in NLTE_BACKENDS:
        raise KeyError(f"Unknown NLTE backend '{region.nlte_backend}'")
    print(f"  [preflight] NLTE backend = {region.nlte_backend}")
    return broadening


# ── Engine ────────────────────────────────────────────────────────────────────

@dataclass
class CNOResult:
    star_id: str
    region: str
    abundances: dict = field(default_factory=dict)       # element -> A(X)
    per_band: list = field(default_factory=list)         # list of band dicts
    iterations: int = 0
    converged: bool = False
    flags: list = field(default_factory=list)
    provenance: dict = field(default_factory=dict)
    uncertainty: dict = field(default_factory=dict)      # element -> {stat, sys, tot}
    phase_a_corrections: list = field(default_factory=list)  # RYA-371 cited 3D/NLTE per diagnostic


def _seed_abundances(star_id, params, codes, solar_A_ispec, feh) -> dict:
    """Initial A(C/N/O/Ni). C from the solar anchor scaled by [Fe/H]; N, O likewise;
    Ni pinned (canonical/EW). The CH/CN/[O I] fits refine C/N/O from here."""
    seed = {}
    for el in ('C', 'N', 'O', 'Ni'):
        seed[el] = float(SOLAR_ASPLUND2021[el]) + (feh if star_id != 'solar' else 0.0)
    return seed


def run_cno(star_id: str, region_name: str = 'vis', *,
            params_override: dict = None, max_iter: int = 5,
            with_systematics: bool = True, out_dir: Path = None,
            tmp_dir: str = '/tmp/ispec_cno') -> CNOResult:
    """Region-aware C/N/O synthesis for `star_id` over `region_name`.

    Flow: seed → fit CH (A(C)) → fit CN given A(C) → fit [O I]+Ni given A(C) →
    re-synthesize (EOS recomputes equilibrium) and re-fit until A(C)/A(N)/A(O)
    converge (Δ < 0.01 dex, max `max_iter`); flag co_equilibrium_not_converged
    otherwise. C I 5052/5380 (LTE) and C2 Swan are fit once as cross-checks.
    """
    region = REGIONS[region_name]
    diagnostics = VIS_DIAGNOSTICS
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)     # RYA-344: TS tmp_dir must exist
    out_dir = Path(out_dir) if out_dir else (Path(PATHS['solar_ew']).parent.parent /
                                             'audit' / 'cno_synthesis')
    out_dir.mkdir(parents=True, exist_ok=True)

    rec = get_star_params(star_id)
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    if params_override:
        params.update(params_override)
    feh = params['feh']

    print(f"\n{'='*72}\n  C/N/O synthesis — {star_id} / {region.name} "
          f"(Teff={params['teff_K']:.0f} logg={params['logg']:.2f} "
          f"[Fe/H]={feh:+.2f} xi={params['vturb_kms']:.2f})\n{'='*72}")

    broadening = preflight(region, star_id, diagnostics)
    nlte_backend = NLTE_BACKENDS[region.nlte_backend]

    atm = _load_atmosphere(params['teff_K'], params['logg'], feh, params['vturb_kms'])
    ll, iso, chem = _load_synth_resources()
    sab = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    obs_w, obs_f = _load_observed_spectrum(star_id)

    codes = _atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)
    solar_A_ispec = _solar_A(('C', 'N', 'O', 'Ni'), chem, sab)
    state = _seed_abundances(star_id, params, codes, solar_A_ispec, feh)
    print(f"  seed A: " + "  ".join(f"{e}={state[e]:.2f}" for e in ('C', 'N', 'O', 'Ni')))

    result = CNOResult(star_id=star_id, region=region.name)
    by_key = {d.key: d for d in diagnostics}

    def _fit(diag, a_center):
        a_lo, a_hi = a_center - 1.2, a_center + 1.2
        t0 = time.time()
        r = _fit_element(obs_w, obs_f, atm, params, diag.element, state, codes,
                         diag.windows_A, diag.use_molecules, broadening,
                         a_lo, a_hi, ll, iso, sab, tmp_dir)
        r['key'] = diag.key
        r['element'] = diag.element
        r['role'] = diag.role
        r['wall_s'] = round(time.time() - t0, 1)
        a_nlte, delta, flag, ref = nlte_backend(diag, r['A_X'], params)
        r.update({'A_X_nlte': a_nlte, 'nlte_delta': delta,
                  'nlte_flag': flag, 'nlte_ref': ref})
        return r

    # ── CNO equilibrium iteration: CH → CN(|C) → [O I](|C) ────────────────────
    primary = {'C': by_key['CH_Gband'], 'N': by_key['CN_red'], 'O': by_key['OI_6300']}
    converged = False
    last = {}
    for it in range(1, max_iter + 1):
        print(f"\n  ── iteration {it} ──")
        deltas = {}
        for el in ('C', 'N', 'O'):
            d = primary[el]
            center = state[el]
            r = _fit(d, center)
            last[d.key] = r
            if np.isfinite(r['A_X']):
                deltas[el] = abs(r['A_X'] - state[el])
                state[el] = r['A_X']     # pin for the next element (EOS coupling)
            else:
                deltas[el] = np.nan
            print(f"    {d.key:10s} A({el})={r['A_X']}  χ²ᵣ={r['red_chi2']}  "
                  f"σfit={r['sigma_fit']}  [{r['status']}]  ({r['wall_s']}s)")
        result.iterations = it
        finite = [v for v in deltas.values() if np.isfinite(v)]
        if finite and max(finite) < 0.01:
            converged = True
            print(f"  converged after {it} iter (max Δ={max(finite):.4f} dex)")
            break
        print(f"  Δ this iter: " + "  ".join(f"{e}={deltas[e]:.3f}" for e in deltas))
    result.converged = converged
    if not converged:
        result.flags.append('co_equilibrium_not_converged')

    # ── Carbon cross-checks (LTE): C I 5052/5380 atomic + C2 Swan ─────────────
    print(f"\n  ── carbon cross-checks (LTE) ──")
    for key in ('CI_5052', 'CI_5380', 'C2_Swan'):
        d = by_key[key]
        r = _fit(d, state['C'])
        last[d.key] = r
        print(f"    {d.key:10s} A(C)={r['A_X']}  χ²ᵣ={r['red_chi2']}  "
              f"[{r['status']}]  flag={r['nlte_flag']}  ({r['wall_s']}s)")

    result.per_band = [last[d.key] for d in diagnostics if d.key in last]
    result.abundances = {'C': state['C'], 'N': state['N'], 'O': state['O']}

    # ── Phase-A cited correction layer (RYA-371): 1D-LTE → cited 3D/NLTE ───────
    corrections = apply_cited_corrections(result.per_band, params, region)
    result.phase_a_corrections = corrections
    print(f"\n  ── Phase-A cited corrections ({region.instrument} arm) — "
          f"validate-don't-tune (cited/vendored only) ──")
    print(f"    {'diagnostic':10s} {'el':2s} {'1D-LTE':>7s} {'corr':>7s} "
          f"{'recon':>7s}  kind / source")
    for c in corrections:
        al = f"{c['a_lte']:.3f}" if isinstance(c['a_lte'], float) and np.isfinite(c['a_lte']) else '  N/A '
        ac = f"{c['a_corr']:.3f}" if isinstance(c['a_corr'], float) and np.isfinite(c['a_corr']) else '  N/A '
        dl = f"{c['delta']:+.3f}" if np.isfinite(c.get('delta', np.nan)) else '   -- '
        print(f"    {c['key']:10s} {c['element']:2s} {al:>7s} {dl:>7s} {ac:>7s}  "
              f"[{c['kind']}] {(c['source'] or '')[:64]}")

    # ── Uncertainty budget (Type A statistical + Type B systematic) ───────────
    result.uncertainty = _uncertainty_budget(
        result, last, by_key, star_id, params, rec, state, codes, obs_w, obs_f,
        atm, broadening, ll, iso, sab, tmp_dir, with_systematics)

    # ── C/O ratio ─────────────────────────────────────────────────────────────
    aC, aO = state['C'], state['O']
    co = float(10 ** (aC - aO)) if np.isfinite(aC) and np.isfinite(aO) else np.nan
    result.abundances['C/O'] = round(co, 3)

    # ── Provenance ─────────────────────────────────────────────────────────────
    result.provenance = {
        'engine': 'pipeline.cno_synthesis (RYA-237)',
        'rt_code': 'turbospectrum',
        'region': region.name, 'instrument': region.instrument, 'R_LSF': region.R,
        'atomic_linelist': 'GESv6_atom_hfs_iso (canonical-gf via gf_resolver, RYA-353)',
        'molecular_lists': f'{_MOLECULES_DIR.name}/*.bsyn (RYA-236: CH/CN/C2/CO/OH/NH)',
        'broadening': {'R': broadening[0], 'vmac': broadening[1], 'vsini': broadening[2],
                       'source': rec.get('source', ''), 'rule': 'per-star RYA-288'},
        'nlte': {'backend': region.nlte_backend,
                 'policy': 'VIS synthesis is 1D-LTE; Phase-A cited corrections applied '
                           'on top (RYA-371): C I 5052/5380 Amarsi-2019 3D-NLTE grid; '
                           '[O I] 6300 cited Caffau-2015 full-3D anchor 8.73 (RYA-367, '
                           'no hardcoded offset); CH/CN/C2 molecular 3D offset OWED'},
        'phase_a_corrections': result.phase_a_corrections,
        'solar_reference': 'Asplund 2021 (A&A 653, A141) via SOLAR_ASPLUND2021',
        'params': params,
        'caveats': ['Ni I 6300.34 gf in the [O I] blend resolves via gf_resolver '
                    '= Johansson+2003 log gf -2.11 (RYA-365). Broader store-#2 '
                    'per-star-master reroute remains the RYA-353 follow-on umbrella.'],
    }

    _write_product(result, out_dir)
    return result


# ── Uncertainty budget ────────────────────────────────────────────────────────

def _uncertainty_budget(result, last, by_key, star_id, params, rec, state, codes,
                        obs_w, obs_f, atm, broadening, ll, iso, sab, tmp_dir,
                        with_systematics) -> dict:
    """Type A (statistical) + Type B (systematic stellar-param sensitivity).

    Statistical:
      C — std across the carbon diagnostics (CH primary + C I 5052/5380 + C2) / √N
          (multi-diagnostic scatter is the honest random term).
      N, O — single primary band → the χ² fit 1σ (sigma_fit).
    Systematic (opt-in, primary band per element): refit the primary band at
      Teff+e_teff, logg+e_logg, [Fe/H]+e_feh, ξ±0.1; quadrature of the ΔA. A
      perturbation whose parameter error is negligible (e.g. solar e_teff=1,
      e_logg≈0) is skipped — keeps the solar validate cheap, real for Procyon.
    """
    budget = {}
    # Statistical
    c_vals = [last[k]['A_X'] for k in ('CH_Gband', 'CI_5052', 'CI_5380', 'C2_Swan')
              if k in last and np.isfinite(last[k]['A_X'])]
    stat = {}
    if len(c_vals) >= 2:
        stat['C'] = float(np.std(c_vals, ddof=1) / np.sqrt(len(c_vals)))
    else:
        stat['C'] = float(last.get('CH_Gband', {}).get('sigma_fit', np.nan))
    stat['N'] = float(last.get('CN_red', {}).get('sigma_fit', np.nan))
    stat['O'] = float(last.get('OI_6300', {}).get('sigma_fit', np.nan))

    # Systematic
    sys_budget = {e: 0.0 for e in ('C', 'N', 'O')}
    if with_systematics:
        primary = {'C': by_key['CH_Gband'], 'N': by_key['CN_red'], 'O': by_key['OI_6300']}
        e_teff = float(rec.get('e_teff', 0.0))
        e_logg = float(rec.get('e_logg', 0.0))
        e_feh = float(rec.get('e_feh', 0.0))
        perturbs = []
        if e_teff >= 5.0:
            perturbs.append(('teff_K', e_teff))
        if e_logg >= 0.02:
            perturbs.append(('logg', e_logg))
        if e_feh >= 0.02:
            perturbs.append(('feh', e_feh))
        perturbs.append(('vturb_kms', 0.1))   # always include a ξ sensitivity term
        print(f"\n  ── systematics: {len(perturbs)} perturbation(s) × 3 primary bands ──")
        for el in ('C', 'N', 'O'):
            d = primary[el]
            base = state[el]
            terms = []
            for pkey, pstep in perturbs:
                pp = dict(params); pp[pkey] = pp[pkey] + pstep
                try:
                    atm_p = _load_atmosphere(pp['teff_K'], pp['logg'], pp['feh'],
                                             pp['vturb_kms'])
                except Exception:
                    continue
                r = _fit_element(obs_w, obs_f, atm_p, pp, el, state, codes,
                                 d.windows_A, d.use_molecules, broadening,
                                 base - 1.0, base + 1.0, ll, iso, sab, tmp_dir)
                if np.isfinite(r['A_X']):
                    terms.append(abs(r['A_X'] - base))
            sys_budget[el] = float(np.sqrt(np.sum(np.square(terms)))) if terms else 0.0
            print(f"    {el}: σ_sys={sys_budget[el]:.3f} dex ({len(terms)} terms)")

    for el in ('C', 'N', 'O'):
        s_a = stat[el] if np.isfinite(stat[el]) else 0.0
        s_b = sys_budget[el]
        budget[el] = {'stat': round(s_a, 3), 'sys': round(s_b, 3),
                      'tot': round(float(np.sqrt(s_a ** 2 + s_b ** 2)), 3)}
    return budget


# ── Output ────────────────────────────────────────────────────────────────────

def _write_product(result: CNOResult, out_dir: Path) -> None:
    rows = []
    for el in ('C', 'N', 'O'):
        unc = result.uncertainty.get(el, {})
        rows.append({
            'element': el, 'A_X': result.abundances.get(el),
            'sigma_stat': unc.get('stat'), 'sigma_sys': unc.get('sys'),
            'sigma_tot': unc.get('tot'),
        })
    rows.append({'element': 'C/O', 'A_X': result.abundances.get('C/O'),
                 'sigma_stat': None, 'sigma_sys': None, 'sigma_tot': None})
    prod = pd.DataFrame(rows)
    band = pd.DataFrame(result.per_band)
    base = out_dir / f'{result.star_id}_{result.region}_cno'
    prod.to_csv(f'{base}_product.csv', index=False)
    band.to_csv(f'{base}_per_band.csv', index=False)
    with open(f'{base}_provenance.json', 'w') as fh:
        json.dump({'abundances': result.abundances, 'iterations': result.iterations,
                   'converged': result.converged, 'flags': result.flags,
                   'uncertainty': result.uncertainty,
                   'provenance': result.provenance}, fh, indent=2, default=str)
    print(f"\n  [out] {base}_product.csv / _per_band.csv / _provenance.json")


# ── Validation (solar VIS gates) ──────────────────────────────────────────────

# Acceptance anchors — Asplund 2021 solar VIS reference (RYA-162 target table).
SOLAR_VIS_GATES = {
    'C': (8.46, 0.05), 'N': (7.83, 0.07), 'O': (8.69, 0.05), 'C/O': (0.59, 0.08),
}


def validate_solar(result: CNOResult) -> bool:
    """Print the solar-VIS gate table; CH vs C I agreement. Returns all-pass bool."""
    print(f"\n{'='*72}\n  SOLAR-VIS ACCEPTANCE GATES (Asplund 2021)\n{'='*72}")
    print(f"  {'qty':6s} {'derived':>10s} {'σ_tot':>7s}  {'target':>8s} {'±':>5s}  result")
    all_pass = True
    for q, (tgt, terr) in SOLAR_VIS_GATES.items():
        val = result.abundances.get(q)
        if val is None or not np.isfinite(val):
            print(f"  {q:6s} {'nan':>10s}  — gate INDETERMINATE"); all_pass = False; continue
        sig = result.uncertainty.get(q, {}).get('tot', 0.0) if q != 'C/O' else 0.0
        within = abs(val - tgt) <= (terr + (sig or 0.0) + 1e-9)
        all_pass = all_pass and within
        print(f"  {q:6s} {val:>10.3f} {sig or 0.0:>7.3f}  {tgt:>8.2f} {terr:>5.2f}  "
              f"{'PASS' if within else 'FAIL'} (Δ={val-tgt:+.3f})")

    # CH vs C I agreement
    bands = {b['key']: b for b in result.per_band}
    ch = bands.get('CH_Gband', {}).get('A_X')
    ci = [bands[k]['A_X'] for k in ('CI_5052', 'CI_5380')
          if k in bands and np.isfinite(bands[k].get('A_X', np.nan))]
    if ch is not None and np.isfinite(ch) and ci:
        ci_mean = float(np.mean(ci))
        agree = abs(ch - ci_mean) <= 0.15
        print(f"\n  CH G-band A(C)={ch:.3f}  vs  C I 5052/5380 A(C)={ci_mean:.3f}  "
              f"Δ={ch-ci_mean:+.3f}  {'AGREE' if agree else 'DISAGREE'} (±0.15)")
    print(f"\n  CO-equilibrium iterations: {result.iterations}  "
          f"(converged={result.converged}); flags={result.flags or 'none'}")
    print(f"  Overall: {'ALL GATES PASS' if all_pass else 'GATES NOT ALL MET'}")
    return all_pass


# ── [O I] 6300 blend partition (RYA-365) ─────────────────────────────────────
# The Ni I 6300.34 gf feeding the [O I] 6300 joint fit must resolve through
# gf_resolver (canonical = Johansson+2003 -2.11). This diagnostic "sees the
# blend": it decomposes total / Ni / [O I] contribution + inferred A(O), BEFORE
# (the stale NIST-grade-B gf -2.310 that RYA-353/354 had seeded) vs AFTER (the
# canonical Johansson-2003 -2.11). The before/after is in-memory only — the live
# linelist already carries the canonical value; we transiently override the one
# Ni row to show the decomposition. No hardcoded gf in the live path.

_NI_6300_WL = 6300.342          # Å (air)
_NI_6300_EP = 4.266             # eV
_STALE_NI_GF = -2.310           # pre-RYA-365 canonical (NIST ASD grade B) — diagnostic only
_OI_PART_WIN_A = (6300.00, 6300.55)   # tight window over the [O I]+Ni blend core for EW


def _ni6300_idx(ll) -> int:
    """Row index of Ni I 6300.34 in the synth linelist (raises if absent)."""
    w = ll['wave_A'].astype(float)
    ep = ll['lower_state_eV'].astype(float)
    for i in np.where((np.abs(w - _NI_6300_WL) <= 0.02) & (np.abs(ep - _NI_6300_EP) <= 0.05))[0]:
        if str(ll['element'][i]).strip().startswith('Ni') and int(ll['ion'][i]) == 1:
            return int(i)
    raise RuntimeError(
        "Ni I 6300.34 absent from the synth linelist — cannot run the [O I] blend "
        "partition (RYA-365). The joint fit requires the Ni blend partner.")


def _ew_mA(sw_nm, flux) -> float:
    """Equivalent width in mÅ of (1 - normalized flux) over sw_nm (a nm grid)."""
    return float(np.trapz(1.0 - np.asarray(flux), np.asarray(sw_nm)) * 1.0e4)


def oi_blend_partition(star_id: str = 'solar', *, tmp_dir: str = '/tmp/ispec_cno',
                       a_C: float = None) -> dict:
    """Decompose the [O I] 6300.30 + Ni I 6300.34 blend and fit A(O), BEFORE
    (stale Ni gf -2.310) vs AFTER (canonical gf_resolver = Johansson+2003 -2.11).

    C/N pinned at the solar anchor (CO coupling to [O I] is negligible at solar);
    A(Ni) pinned from the canonical solar Ni abundance (Asplund 2021, sourced),
    not a literal. Returns {'before': {...}, 'after': {...}, 'canon_gf': ...}.
    """
    region = REGIONS['vis']
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    rec = get_star_params(star_id)
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    feh = params['feh']
    off = (feh if star_id != 'solar' else 0.0)

    print(f"\n{'='*72}\n  [O I] 6300 BLEND PARTITION — {star_id} (RYA-365)\n{'='*72}")
    broadening = preflight(region, star_id, [d for d in VIS_DIAGNOSTICS if d.key == 'OI_6300'])

    atm = _load_atmosphere(params['teff_K'], params['logg'], feh, params['vturb_kms'])
    ll, iso, chem = _load_synth_resources()        # already canonical-gf (RYA-353)
    sab = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    obs_w, obs_f = _load_observed_spectrum(star_id)
    codes = _atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)

    # CRITICAL: the live Ni gf MUST equal the gf_resolver canonical (no hardcoded copy).
    ni_i = _ni6300_idx(ll)
    canon_gf = float(resolve_gf((28, 1), _NI_6300_WL, _NI_6300_EP))
    live_gf = float(ll['loggf'][ni_i])
    if abs(live_gf - canon_gf) > 1e-4:
        raise AssertionError(
            f"Ni I 6300.34 synth gf {live_gf:+.4f} != gf_resolver canonical "
            f"{canon_gf:+.4f} — the [O I] path is NOT resolving via gf_resolver "
            f"(RYA-365 invariant violated).")
    print(f"  Ni I 6300.34 gf via gf_resolver = {canon_gf:+.3f} "
          f"(Johansson+2003); live synth row matches ✓")

    a_Ni = float(SOLAR_ASPLUND2021['Ni']) + off      # pinned, canonical solar Ni
    a_C0 = float(a_C if a_C is not None else SOLAR_ASPLUND2021['C']) + (0.0 if a_C is not None else off)
    a_N0 = float(SOLAR_ASPLUND2021['N']) + off
    print(f"  pinned: A(Ni)={a_Ni:.2f} (Asplund2021 solar Ni)  A(C)={a_C0:.2f}  A(N)={a_N0:.2f}")

    win = (_OI_PART_WIN_A,)
    fit_win = ((6299.5, 6301.0),)        # same fit window run_cno uses for OI_6300
    sw = np.arange(_OI_PART_WIN_A[0] / 10.0, _OI_PART_WIN_A[1] / 10.0 + _WSTEP_NM * 0.5,
                   _WSTEP_NM)

    out = {'canon_gf': canon_gf}
    for tag, ni_gf in (('before', _STALE_NI_GF), ('after', canon_gf)):
        ll_v = ll.copy()
        ll_v['loggf'][ni_i] = ni_gf
        seed_O = float(SOLAR_ASPLUND2021['O']) + off
        rfit = _fit_element(obs_w, obs_f, atm, params, 'O',
                            {'C': a_C0, 'N': a_N0, 'O': seed_O, 'Ni': a_Ni}, codes,
                            fit_win, True, broadening,
                            seed_O - 1.2, seed_O + 1.2, ll_v, iso, sab, tmp_dir)
        a_O = rfit['A_X']
        st_full = {'C': a_C0, 'N': a_N0, 'O': a_O, 'Ni': a_Ni}
        st_noNi = dict(st_full, Ni=a_Ni - 6.0)       # remove Ni (6 dex weaker → ~0)
        st_noO = dict(st_full, O=a_O - 6.0)          # remove [O I]
        f_full = _synth_window(sw, atm, params, ll_v, iso, sab, _fixed_ab(st_full, codes),
                               broadening, True, tmp_dir)
        f_noNi = _synth_window(sw, atm, params, ll_v, iso, sab, _fixed_ab(st_noNi, codes),
                               broadening, True, tmp_dir)
        f_noO = _synth_window(sw, atm, params, ll_v, iso, sab, _fixed_ab(st_noO, codes),
                              broadening, True, tmp_dir)
        ew_tot = _ew_mA(sw, f_full)
        ew_ni = ew_tot - _ew_mA(sw, f_noNi)
        ew_o = ew_tot - _ew_mA(sw, f_noO)
        depth = float(1.0 - np.min(f_full))
        out[tag] = {'ni_gf': ni_gf, 'A_O': a_O, 'red_chi2': rfit['red_chi2'],
                    'ew_total_mA': round(ew_tot, 2), 'ew_Ni_mA': round(ew_ni, 2),
                    'ew_OI_mA': round(ew_o, 2), 'core_depth': round(depth, 3),
                    'ni_frac': round(ew_ni / ew_tot, 3) if ew_tot else float('nan')}
    return out


def print_oi_partition(part: dict) -> bool:
    """Print the before/after blend-partition table; return whether O clears the gate."""
    tgt, terr = SOLAR_VIS_GATES['O']
    print(f"\n  [O I] 6300 blend partition (EW over {_OI_PART_WIN_A[0]}-{_OI_PART_WIN_A[1]} Å):")
    print(f"  {'case':<26}{'Ni gf':>8}{'EW_tot':>9}{'EW_Ni':>8}{'EW_[OI]':>9}"
          f"{'Ni frac':>9}{'A(O)':>8}{'χ²ᵣ':>8}")
    for tag, label in (('before', 'before (NIST B)'), ('after', 'after (Johansson03)')):
        p = part[tag]
        print(f"  {label:<26}{p['ni_gf']:>+8.3f}{p['ew_total_mA']:>9.2f}"
              f"{p['ew_Ni_mA']:>8.2f}{p['ew_OI_mA']:>9.2f}{p['ni_frac']:>9.3f}"
              f"{p['A_O']:>8.3f}{p['red_chi2']:>8.3f}")
    aO = part['after']['A_O']
    within = abs(aO - tgt) <= terr + 1e-9
    print(f"\n  A(O)☉: before={part['before']['A_O']:.3f}  →  after={aO:.3f}  "
          f"(target {tgt:.2f} ± {terr:.2f}, Δ={aO - tgt:+.3f})")
    print(f"  GATE: {'PASS — solar O cleared' if within else 'FAIL — RCA finding, do NOT force (see RYA-354)'}")
    return within


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(description='Region-aware C/N/O synthesis (RYA-237)')
    ap.add_argument('--star', default='solar')
    ap.add_argument('--region', default='vis', choices=sorted(REGIONS))
    ap.add_argument('--species', default=None,
                    help='restrict to a channel, e.g. "O I" → the [O I] 6300 blend '
                         'partition (Ni vs [O I], before/after canonical gf). '
                         'Omit for the full C/N/O run.')
    ap.add_argument('--validate', action='store_true',
                    help='print the solar-VIS acceptance gate table')
    ap.add_argument('--no-systematics', action='store_true',
                    help='skip the Type-B stellar-parameter sensitivity refits')
    ap.add_argument('--max-iter', type=int, default=5)
    ap.add_argument('--out', default=None)
    args = ap.parse_args(argv)

    # --species "O I" → the focused [O I] 6300 blend-partition diagnostic (RYA-365).
    if args.species and args.species.replace(' ', '').upper() in ('OI', 'O'):
        if args.region != 'vis':
            raise SystemExit("[O I] 6300 blend partition is a VIS diagnostic (--region vis).")
        part = oi_blend_partition(args.star, tmp_dir='/tmp/ispec_cno')
        if args.validate:
            print_oi_partition(part)
        return part

    result = run_cno(args.star, args.region, max_iter=args.max_iter,
                     with_systematics=not args.no_systematics,
                     out_dir=Path(args.out) if args.out else None)
    if args.validate and args.star == 'solar':
        validate_solar(result)
    elif args.validate:
        print(f"\n  --validate: gate table is defined for solar only; "
              f"{args.star} compared against GBS downstream (RYA-348).")
    return result


if __name__ == '__main__':
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        main()

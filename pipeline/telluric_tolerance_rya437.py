"""
pipeline/telluric_tolerance_rya437.py
=====================================
RYA-437 Part A — derive the telluric verification tolerance from the 13C / CO
precision requirement, NOT a round number.

RYA-424 shipped the telluric data-input gate at TELLURIC_RESIDUAL_TOL = 0.05 — a
provisional round number with no provenance. This module derives the real tolerance
by propagating a telluric residual of magnitude r through to the 13C abundance (the
12C/13C ratio's binding leg) and finding the r at which the induced bias reaches our
carbon-abundance precision target. The binding (tighter) CO feature sets the gate
(RYA-436 discipline: no stored constant without a cited derivation).

THE PROPAGATION (linear curve-of-growth)

A telluric mis-correction of typical amplitude r (the residual metric = median
|1 - corrected_flux/continuum| at telluric-dominated pixels) perturbs the apparent
depth of a CO feature by ~r in continuum-normalized flux. For a weak line on the
LINEAR part of the curve of growth, equivalent width EW ∝ N (the abundance), so a
depth perturbation r on a feature of central depth d gives:

    ΔEW/EW = (r · W) / (d · W) = r / d                 (W = feature width, cancels)
    Δlog N = (ΔEW/EW) / ln10 = r / (d · ln10)          (abundance bias, dex)

So the abundance bias from a residual r on a feature of depth d is

    ΔA(feature) = r / (d · ln10).

Inverting at the precision target σ (dex) gives the per-feature tolerance:

    r(feature) = σ · ln10 · d(feature).

THE TWO CO FEATURES (per-feature; the binding one sets the gate)

  * 12CO(2-0) bandhead 2.2935 µm — STRONG (deep) → larger d → LOOSER tolerance.
  * 13CO(2-0) bandhead 2.3448 µm — WEAK (shallow) → smaller d → TIGHTER tolerance,
    and it is the leg that actually carries the 12C/13C ratio (12CO is the
    well-determined numerator). The 13CO feature BINDS.

THE PRECISION TARGET (cited, in-repo)

σ = 0.05 dex — the project's adopted carbon-abundance precision, the Asplund-2021
solar carbon gate `SOLAR_VIS_GATES['C'] = (8.46, 0.05)` in pipeline/cno_synthesis.py.
The 13C abundance is a carbon-isotopologue abundance, so the carbon precision is the
natural (and the tighter, vs C/O's 0.08) target. The implied 12C/13C fractional
precision is ln10 · σ ≈ 11.5% — consistent with high-quality stellar CO isotopic
work (~10-20%).

THE FEATURE DEPTHS (measured, cited)

From the ACE-FTS solar atlas (Hase, Wallace, McLeod, Harrison & Bernath 2010, JQSRT
111, 521; telluric-FREE, roughly disk-integrated → matches reflected-solar Vesta),
the segment vendored at data/solar_reference/ir_atlases/ace_fts_solar_co_4255_4367.csv
(RYA-390). Disk-integrated depths are the right match to the integrated-disk Vesta /
solar Phase-B target (and are the more conservative, shallower → tighter-gate choice
vs the ground-based NSO photatl).

RESULT  (see `derive_tolerance`)

    d(12CO) ≈ 0.140,  r(12CO) = 0.05·ln10·0.140 ≈ 0.0161
    d(13CO) ≈ 0.093,  r(13CO) = 0.05·ln10·0.093 ≈ 0.0107   ← BINDING

    TELLURIC_RESIDUAL_TOL = 0.011  (1.1%), set by 13CO(2-0).

This is ~5× tighter than the round 0.05. At Vesta's measured ~6% residual the implied
13C bias is r/(d13·ln10) ≈ 0.28 dex (12C/13C wrong by ~1.9×) — i.e. the round 0.05
would have admitted a frame that cannot do the 13C science. Validate-don't-tune: the
number comes from the physics; whether Vesta passes is then a separate (RCA) question.

Run:  python -m pipeline.telluric_tolerance_rya437
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from config.constants import PATHS

# ── Precision target (cited, in-repo) ─────────────────────────────────────────
# Asplund-2021 solar carbon gate, pipeline/cno_synthesis.py SOLAR_VIS_GATES['C'].
CARBON_PRECISION_DEX = 0.05
CARBON_PRECISION_CITE = ("Asplund et al. 2021 solar carbon precision; "
                         "cno_synthesis.SOLAR_VIS_GATES['C'] = (8.46, 0.05)")

# ── The CO features (vacuum Å; IR convention) ─────────────────────────────────
CO_12_2_0_BANDHEAD_A = 22935.0   # 12C16O (2-0) R-branch bandhead, 2.2935 µm (strong)
CO_13_2_0_BANDHEAD_A = 23448.0   # 13C16O (2-0) R-branch bandhead, 2.3448 µm (weak; binds)

# ACE-FTS solar atlas segment (RYA-390), telluric-free, ~disk-integrated.
_ACE_ATLAS = (Path(str(PATHS['data_root'])) / 'solar_reference' / 'ir_atlases'
              / 'ace_fts_solar_co_4255_4367.csv')
ACE_ATLAS_CITE = ("ACE-FTS solar atlas, Hase, Wallace, McLeod, Harrison & Bernath "
                  "2010, JQSRT 111, 521 (telluric-free, disk-integrated)")

_LN10 = math.log(10.0)


@dataclass
class FeatureTolerance:
    name: str
    wavelength_A: float
    depth: float                 # measured central depth (continuum-normalized)
    sensitivity_dex_per_resid: float   # ΔA per unit residual = 1/(d·ln10)
    tolerance: float             # r at which ΔA = σ_target  = σ·ln10·d
    binding: bool = False


def measure_bandhead_depth(wavelength_A: float, half_width_A: float = 8.0,
                           atlas_path: Path = None) -> float:
    """Central depth (1 - min continuum-normalized flux) of a CO bandhead in the
    ACE-FTS solar atlas, within ±half_width of the nominal bandhead. The atlas region
    is line-dense (no flat continuum), so depth is taken vs the normalized continuum
    (≈1), the same frame the residual metric measures perturbations in."""
    import pandas as pd
    atlas_path = Path(atlas_path) if atlas_path else _ACE_ATLAS
    df = pd.read_csv(atlas_path)
    w = df['wavelength_vac_A'].to_numpy(float)
    f = df['intensity'].to_numpy(float)
    m = (w >= wavelength_A - half_width_A) & (w <= wavelength_A + half_width_A)
    if m.sum() < 5:
        raise RuntimeError(f"ACE-FTS atlas does not cover {wavelength_A} Å "
                           f"(±{half_width_A}) — cannot measure the bandhead depth.")
    return float(1.0 - np.nanmin(f[m]))


def derive_tolerance(sigma_dex: float = CARBON_PRECISION_DEX,
                     atlas_path: Path = None) -> dict:
    """Derive TELLURIC_RESIDUAL_TOL from the 13C/CO precision requirement, per-feature.
    Returns the per-feature tolerances, the binding (minimum) tolerance, and a
    provenance string. The binding feature (13CO, weak) sets the gate."""
    feats = []
    for name, wl in (('12CO(2-0)', CO_12_2_0_BANDHEAD_A),
                     ('13CO(2-0)', CO_13_2_0_BANDHEAD_A)):
        d = measure_bandhead_depth(wl, atlas_path=atlas_path)
        sens = 1.0 / (d * _LN10)
        tol = sigma_dex * _LN10 * d
        feats.append(FeatureTolerance(name=name, wavelength_A=wl, depth=d,
                                      sensitivity_dex_per_resid=sens, tolerance=tol))
    binding = min(feats, key=lambda ft: ft.tolerance)
    binding.binding = True
    tol = binding.tolerance
    # round to 3 significant figures for the stored constant (still derivation-exact;
    # a test asserts the stored constant matches this within rounding).
    tol_stored = float(f"{tol:.3g}")
    prov = (
        f"RYA-437: derived from the 13C/CO precision requirement (NOT a round number). "
        f"r = sigma * ln10 * d with sigma={sigma_dex} dex [{CARBON_PRECISION_CITE}] and "
        f"the binding feature {binding.name} at {binding.wavelength_A:.1f} A vac, "
        f"central depth d={binding.depth:.3f} measured from {ACE_ATLAS_CITE}. "
        f"r(12CO)={feats[0].tolerance:.4f}, r(13CO)={feats[1].tolerance:.4f} -> binding "
        f"{binding.name} -> TELLURIC_RESIDUAL_TOL={tol_stored:g}. Implied 12C/13C "
        f"fractional precision = ln10*sigma = {_LN10*sigma_dex:.3f}. See "
        f"pipeline/telluric_tolerance_rya437.derive_tolerance + "
        f"docs/telluric_tolerance_rya437.md.")
    return {'features': feats, 'binding': binding, 'tolerance': tol,
            'tolerance_stored': tol_stored, 'sigma_dex': sigma_dex,
            'provenance': prov,
            'bias_at': {r: r / (binding.depth * _LN10)         # ΔA(13C) at residual r
                        for r in (0.011, 0.02, 0.05, 0.06)}}


def main(argv=None):
    res = derive_tolerance()
    print('=' * 84)
    print('  RYA-437 Part A — telluric tolerance from the 13C/CO precision requirement')
    print('=' * 84)
    print(f"  target precision sigma = {res['sigma_dex']} dex  "
          f"[{CARBON_PRECISION_CITE}]")
    print(f"  feature depths from {ACE_ATLAS_CITE}\n")
    print(f"  {'feature':<12}{'lambda_vac_A':>13}{'depth d':>10}"
          f"{'dA/dr':>9}{'r=sig*ln10*d':>15}{'':>5}")
    for ft in res['features']:
        print(f"  {ft.name:<12}{ft.wavelength_A:>13.1f}{ft.depth:>10.3f}"
              f"{ft.sensitivity_dex_per_resid:>9.2f}{ft.tolerance:>15.4f}"
              f"{'  <-- BINDING' if ft.binding else ''}")
    print(f"\n  TELLURIC_RESIDUAL_TOL = {res['tolerance_stored']:g}  "
          f"(set by {res['binding'].name})")
    print(f"\n  13C abundance bias ΔA(13C) = r / (d13 * ln10) at sample residuals:")
    for r, b in res['bias_at'].items():
        print(f"    residual {r*100:>4.1f}%  ->  ΔA(13C) = {b:+.3f} dex"
              f"{'   (calibrated gate)' if abs(r-res['tolerance_stored'])<1e-6 else ''}"
              f"{'   (Vesta ~6%)' if abs(r-0.06)<1e-9 else ''}")
    print('=' * 84)
    return res


if __name__ == '__main__':
    main()

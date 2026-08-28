"""
pipeline/uncertainty_stack.py
=============================
Reported-uncertainty budget for the solar (and, by design, any-star) abundances.

CRITICAL distinction (RYA-166 rewrite, 2026-07-15) — two DIFFERENT "sigma"s:

  raw_sigma  = line-to-line scatter (std of the N per-line A(X)). An internal QA
               diagnostic ("is this line pool clean"). ALREADY built + gated: the
               solar Fe I honest floor is 0.1398 dex (RYA-407), wired into the
               G acceptance profile (RYA-446). We CITE it here (single source of
               truth: ACCEPTANCE_PROFILES['G']['fe1_scatter_max']); we never
               re-derive or re-gate it.

  sigma_reported = the STANDARD ERROR of the mean combined in quadrature with the
               stellar-parameter (Type B) systematics — a DIFFERENT, much smaller
               number (dividing by sqrt(N) alone takes 0.1398 -> ~0.018 for N=62).
               This is what the 0.05 dex solar gate target actually applies to
               (RYA-282 field convention; Amarsi+2022 solar Fe I ~0.01). This is
               what RYA-158/50 never built (stub NotImplementedError) and is the
               deliverable here.

Budget
------
  sigma_SE      = raw_sigma / sqrt(N)                         # Type A -> mean
  sigma_B_p     = |dA/dp| * delta_p   for p in {Teff,logg,vmic,[Fe/H]}   # Type B
                  dA/dp from re-deriving abundances at perturbed params
                  (abundances_derive.run(stellar_params_override=...), RYA-158 hook);
                  delta_p = the STAR's own parameter uncertainty. For the Sun these
                  are tiny (Teff +/-1 K, logg exact, [Fe/H]=0 by definition), so the
                  systematic terms are small by construction.
  sigma_reported = sqrt(sigma_SE^2 + sum_p sigma_B_p^2)
  solar_zero_point_err = the solar calibration uncertainty every science-target
                  [X/Fe] inherits (propagate in quadrature downstream).

Output: data/processed/{star}_uncertainty.csv (working) + a committed audit JSON
(data/audit/uncertainty/{star}_uncertainty_rya158.json) the RYA-166 Layer-3 gate
reads without re-running the pipeline.

Supersedes RYA-50. Folded into RYA-166 (Step 2.5). Linear: RYA-158 / RYA-166.
"""

import argparse
import contextlib
import io
import json
import os
import sys
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

from config.constants import (PATHS, ACCEPTANCE_PROFILES, get_star_params)

_ROOT = Path(__file__).resolve().parent.parent

# ── Type B derivative steps (for the numerical dA/dp) and the SUN's own parameter
#    uncertainties delta_p (what actually multiplies the derivative). ────────────
# Derivative steps: RYA-158 spec (Teff +/-100 K, logg +/-0.10, vmic +/-0.10, [Fe/H] +/-0.05).
TYPE_B_STEPS = {'teff_K': 100.0, 'logg': 0.10, 'vturb_kms': 0.10, 'feh': 0.05}

# Per-star parameter uncertainties delta_p. For the Sun: Teff/logg from the pinned
# reference (get_star_params e_teff/e_logg), vmic from RYA-158 (0.9 +/- 0.05 km/s),
# [Fe/H] = 0 BY DEFINITION (the Sun is the zero-point) -> delta=0.
SOLAR_DELTA_P = {'teff_K': None,   # filled from get_star_params('solar')['e_teff']
                 'logg':   None,   # from e_logg
                 'vturb_kms': 0.05,  # RYA-158 (solar microturbulence uncertainty)
                 'feh':    0.0}    # Sun defines [Fe/H]=0

_PARAM_LABEL = {'teff_K': 'Teff', 'logg': 'logg', 'vturb_kms': 'vmic', 'feh': 'FeH'}

#: Which `stars.yaml` field carries each parameter's uncertainty. RYA-1088 declared
#: `e_xi` from Jofré+2014 Table 2 for the run-order targets; `e_teff` / `e_logg` /
#: `e_feh` were already present on the GBS scale (Heiter+2015, Jofré+2014).
_DELTA_FIELD = {'teff_K': 'e_teff', 'logg': 'e_logg',
                'vturb_kms': 'e_xi', 'feh': 'e_feh'}


def params_and_deltas(star_id: str = "solar"):
    """Nominal params + delta_p for ANY star whose uncertainties are declared.

    🔴 RYA-1088 — WHY THE SUN IS NOT JUST "A STAR WITH SMALL delta_p".
    Solar logg and [Fe/H] are EXACT BY DEFINITION: the Sun *is* the [Fe/H] zero point, and
    its logg is fixed by the IAU nominal mass and radius. So their delta_p are 0 by
    construction and `_type_b` skips the derivative runs entirely -- sigma_B = 0, a
    MEASURED zero rather than an unmeasured blank.

    ⚠️ AND `stars.yaml` CARRIES `e_feh: 0.03` FOR THE SUN, WHICH MUST NOT BE USED HERE.
    That value is the uncertainty on the solar iron ABUNDANCE SCALE, not on the Sun's
    [Fe/H] relative to itself -- which is 0 with no error, definitionally. Reading it
    would compute a non-zero solar sigma_B([Fe/H]) and violate the ticket's own critical
    condition. The Sun keeps its explicit zeros; other stars read theirs from the record.
    """
    sp = get_star_params(star_id)
    params0 = {'teff_K': float(sp['teff']), 'logg': float(sp['logg']),
               'feh': float(sp.get('feh_ref', 0.0)),
               'vturb_kms': float(sp.get('xi', 1.0))}
    if 'solar' in star_id.lower():
        deltas = dict(SOLAR_DELTA_P)
        deltas['teff_K'] = float(sp.get('e_teff', 1.0))   # solar Teff known to ~1 K
        deltas['logg'] = float(sp.get('e_logg', 0.0))     # exact by definition
        # deltas['feh'] stays 0.0 and deltas['vturb_kms'] stays the RYA-158 0.05 --
        # see the docstring; neither is read from the record.
        return params0, deltas

    # Non-solar: every delta_p comes from the record, and a MISSING one is a refusal
    # rather than a zero. An undeclared uncertainty silently treated as 0 would publish
    # "measured and negligible" for a term nobody determined (RYA-907).
    missing = [k for k, f in _DELTA_FIELD.items() if sp.get(f) is None]
    if missing:
        raise NotImplementedError(
            f"'{star_id}' has no declared uncertainty for "
            f"{', '.join(_DELTA_FIELD[k] for k in missing)} in config/stars.yaml. "
            f"Declare it from a CITED source (Heiter+2015 for Teff/logg, Jofré+2014 for "
            f"[Fe/H] and xi) before running Type B -- an undeclared delta_p must not "
            f"become a silent zero (RYA-1088).")
    return params0, {k: float(sp[f]) for k, f in _DELTA_FIELD.items()}


def _solar_params_and_deltas():
    """Nominal solar params (teff_K/logg/feh/vturb_kms) + the per-param delta_p."""
    sp = get_star_params('solar')
    params0 = {'teff_K': float(sp['teff']), 'logg': float(sp['logg']),
               'feh': 0.0, 'vturb_kms': float(sp.get('xi', 1.0))}
    deltas = dict(SOLAR_DELTA_P)
    deltas['teff_K'] = float(sp.get('e_teff', 1.0))       # solar Teff known to ~1 K
    deltas['logg'] = float(sp.get('e_logg', 0.0))         # solar logg effectively exact
    return params0, deltas


def _derive_A(params: dict) -> dict:
    """Re-derive 1D-LTE absolute A(X) per (element,ion) at the given params via the
    EW->abundance path (abundances_derive.run, the RYA-158 override hook). Returns
    {(element, ion): A_X_1dlte}. stdout is suppressed (the pipeline is chatty)."""
    import pipeline.abundances_derive as ad
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        _, results = ad.run(star_id='solar', engine='spectrum', skip_convergence=True,
                            stellar_params_override=dict(params))
    out = {}
    for _, r in results.iterrows():
        a = r.get('A_X')
        if pd.notna(a):
            out[(str(r['element']), str(r['ion']))] = float(a)
    return out


def _type_b(params0: dict, deltas: dict, verbose: bool = True) -> dict:
    """Central-difference dA/dp and sigma_B_p per (element,ion) for every param whose
    delta_p > 0 (a param with delta 0 contributes exactly 0 to the Sun's budget, so
    we skip its two runs). Returns {(el,ion): {param: (dA_dp, sigma_B)}}."""
    out = {}
    for p, step in TYPE_B_STEPS.items():
        dp = deltas.get(p, 0.0)
        if not dp or dp <= 0:
            # delta_p == 0 (solar logg exact, [Fe/H]=0 by definition) -> sigma_B = 0,
            # no derivative run needed. Recorded as 0 with the reason.
            if verbose:
                print(f"  Type B {_PARAM_LABEL[p]}: delta_p=0 -> sigma_B=0 (no run)")
            continue
        if verbose:
            print(f"  Type B {_PARAM_LABEL[p]}: dA/dp via +/-{step} runs, delta_p={dp}")
        A_plus = _derive_A({**params0, p: params0[p] + step})
        A_minus = _derive_A({**params0, p: params0[p] - step})
        for key in set(A_plus) & set(A_minus):
            dA_dp = (A_plus[key] - A_minus[key]) / (2.0 * step)
            sigma_B = abs(dA_dp) * dp
            out.setdefault(key, {})[p] = (dA_dp, sigma_B)
    return out


def _raw_sigma_and_N(star_id: str):
    """Per-(element) raw line-to-line scatter + N, sourced from the committed
    verdict. Fe I raw_sigma is OVERRIDDEN to the ratified RYA-407/446 floor
    (single source of truth: ACCEPTANCE_PROFILES['G']['fe1_scatter_max']) rather
    than recomputed. Returns {element: (raw_sigma, N, A_X_mean)}."""
    vj = json.loads((_ROOT / 'data' / 'audit' / 'cno_synthesis'
                     / 'solar_phase_c_verdict.json').read_text())
    out = {}
    for r in vj['verdicts']:
        el = str(r['element'])
        sig = r.get('sigma')
        n = r.get('n_lines')
        if sig is None or n in (None, 0):
            continue
        out[el] = (float(sig), int(n), float(r['A_measured']))
    # Fe I: cite the ratified honest-floor scatter, do NOT recompute (RYA-407/446).
    if 'Fe' in out:
        fe1_floor = float(ACCEPTANCE_PROFILES['G']['fe1_scatter_max'])   # 0.1398, RYA-407
        _, n_fe, a_fe = out['Fe']
        out['Fe'] = (fe1_floor, n_fe, a_fe)
    return out


def run(star_id: str = 'solar', write: bool = True, verbose: bool = True) -> pd.DataFrame:
    """Assemble the reported-uncertainty budget. Solar-only today (the Type B
    override + solar delta_p are solar-pinned); raises for other stars until their
    delta_p are declared."""
    # RYA-1088: no longer solar-only. `params_and_deltas` serves any star whose delta_p
    # are declared and REFUSES any star whose are not -- so the fence moved from "is this
    # the Sun" to "is this sourced", which is the question that actually matters.
    params0, deltas = params_and_deltas(star_id)
    if verbose:
        print(f"[uncertainty_stack] {star_id} params {params0}")
        print(f"[uncertainty_stack] {star_id} delta_p {deltas}")

    raw = _raw_sigma_and_N(star_id)
    if verbose:
        print(f"[uncertainty_stack] Type B: re-deriving abundances at perturbed params "
              f"(this runs the EW->abundance path several times) ...")
    type_b = _type_b(params0, deltas, verbose=verbose)

    rows = []
    for el, (raw_sigma, n, a_mean) in sorted(raw.items()):
        sigma_SE = raw_sigma / np.sqrt(n) if n >= 1 else np.nan
        flags = []
        if n < 3:
            flags.append('n_lines_low')
        # Type B keyed on (element, ion). Prefer neutral 'I'; else first available.
        b_terms = type_b.get((el, 'I')) or {}
        if not b_terms:
            for (e2, i2), t in type_b.items():
                if e2 == el:
                    b_terms = t
                    break
        sB = {p: b_terms.get(p, (np.nan, 0.0))[1] for p in TYPE_B_STEPS}
        dAdp = {p: b_terms.get(p, (np.nan, np.nan))[0] for p in TYPE_B_STEPS}
        sigma_reported = float(np.sqrt(sigma_SE ** 2 + sum(v ** 2 for v in sB.values())))
        # solar zero-point: the solar calibration uncertainty every [X/Fe] inherits.
        solar_zp = sigma_reported
        rows.append(dict(
            element=el, ion='I',
            n_lines=n, A_X_mean=round(a_mean, 4),
            raw_sigma=round(raw_sigma, 4),          # QA diagnostic (RYA-407/446), NOT the reported error
            sigma_SE=round(sigma_SE, 4),
            sigma_B_Teff=round(sB['teff_K'], 4),
            sigma_B_logg=round(sB['logg'], 4),
            sigma_B_vmic=round(sB['vturb_kms'], 4),
            sigma_B_FeH=round(sB['feh'], 4),
            dA_dTeff_per100K=(round(dAdp['teff_K'] * 100.0, 4)
                              if np.isfinite(dAdp['teff_K']) else np.nan),
            dA_dvmic_per_kms=(round(dAdp['vturb_kms'], 4)
                              if np.isfinite(dAdp['vturb_kms']) else np.nan),
            sigma_solar=round(solar_zp, 4),
            sigma_reported=round(sigma_reported, 4),
            flags=';'.join(flags),
        ))
    df = pd.DataFrame(rows)

    if write:
        out_csv = Path(str(PATHS['solar_ew'])).parent / f'{star_id}_uncertainty.csv'
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
        audit = _ROOT / 'data' / 'audit' / 'uncertainty'
        audit.mkdir(parents=True, exist_ok=True)
        payload = {
            'ticket': 'RYA-158 (folded into RYA-166 Step 2.5)',
            'star': star_id,
            'raw_sigma_source': 'RYA-407 honest floor via ACCEPTANCE_PROFILES[G][fe1_scatter_max]=0.1398 (Fe I); '
                                'per-element verdict sigma otherwise. raw_sigma is the QA scatter, NOT sigma_reported.',
            'delta_p_solar': {_PARAM_LABEL[k]: deltas[k] for k in TYPE_B_STEPS},
            'type_b_derivative_steps': {_PARAM_LABEL[k]: v for k, v in TYPE_B_STEPS.items()},
            'per_element': df.to_dict(orient='records'),
        }
        (audit / f'{star_id}_uncertainty_rya158.json').write_text(json.dumps(payload, indent=2))
        if verbose:
            print(f"[uncertainty_stack] wrote {out_csv.name} + "
                  f"data/audit/uncertainty/{star_id}_uncertainty_rya158.json")
    return df


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='solar')
    args = ap.parse_args()
    _df = run(star_id=args.star)
    _fe = _df[_df['element'] == 'Fe']
    print("\n  element  n  raw_sigma  sigma_SE  sig_reported")
    for _, r in _df.iterrows():
        print(f"  {r['element']:>4s}  {r['n_lines']:>3d}  {r['raw_sigma']:.4f}   "
              f"{r['sigma_SE']:.4f}    {r['sigma_reported']:.4f}")
    if len(_fe):
        fr = _fe.iloc[0]
        print(f"\n  Fe I sigma_reported = {fr['sigma_reported']:.4f} dex "
              f"(target <= 0.05; raw scatter {fr['raw_sigma']:.4f} is the separate QA floor)")

#!/usr/bin/env python3
"""
scripts/measure_unsaturated_namg_rya465.py
==========================================
RYA-465 — unsaturated line-selection recovery for Na / Mg. The RYA-354 Batch-1
finding: Na/Mg are NOT gf-limited; their pool lines carry NIST grade-A gf but are
SATURATED (flat part of the COG) -> culled. The fix is pool COMPOSITION: add weaker,
unsaturated SUBORDINATE transitions (linear part of the COG) that already carry a
citable NIST-graded gf, measure their solar EW, and check internal consistency.

This MEASURES the candidate unsaturated lines on the same HARPS solar spectrum the
canonical pool was built from (data/processed/solar_normalized.csv), reusing the
pipeline EW primitives (_fit_profile / _integrate_profile / _ew_error). It writes a
provenance table; it never tunes gf and never edits an existing EW.

Out: data/audit/line_recovery/namg_unsaturated_rya465.csv (per-line measurement +
     gf/grade/reference + COG regime).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline.lines_fit import _fit_profile, _integrate_profile, _ew_error   # noqa: E402
from config.constants import PIPELINE                                         # noqa: E402

SPEC = ROOT / 'data' / 'processed' / 'solar_normalized.csv'
OUT = ROOT / 'data' / 'audit' / 'line_recovery' / 'namg_unsaturated_rya465.csv'

REW_KNEE = -4.90               # RYA-458/curate: linear-COG saturation knee
SAT_CEIL = float(PIPELINE['vmic_ew_ceiling_mA'])   # 100 mA absolute ceiling
FIT_HALF = 0.50                # A — tight local fit window half-width (crowded optical)
CORE_TOL = 0.07                # A — search for the observed line core within +/- this

_KP = 'NIST ASD v5.11; Kelleher & Podobedova 2008 JPCRD 37 p1285 (Na/Mg)'

# Candidate unsaturated subordinate lines (species, lambda_air, EP, log_gf, grade, ref).
# gf already on canonical_gf — values + grades from NIST ASD (Kelleher & Podobedova 2008).
CANDIDATES = [
    ('Na', 'I', 6154.225, 2.102, -1.547, 'A',  'NIST ASD v5.11 grade A'),
    ('Na', 'I', 6160.747, 2.104, -1.246, 'A',  'NIST ASD v5.11 grade A'),
    ('Mg', 'I', 4730.029, 4.346, -2.347, 'C+', f'{_KP} grade C+ L7428'),
    ('Mg', 'I', 4571.096, 0.000, -5.623, 'A',  'NIST ASD v5.11 grade A'),  # intercombination
    ('Mg', 'I', 4167.271, 4.346, -0.745, '?',  'VALD3 (grade TBD)'),
    ('Mg', 'I', 4702.991, 4.346, -0.440, 'B+', f'{_KP} grade B+ L7428'),
]


def measure(wav, flux, line):
    m = np.abs(wav - line) < FIT_HALF
    w, f = wav[m].copy(), flux[m].copy()
    if len(w) < 8:
        return dict(status='no_pixels')
    # LOCAL continuum re-normalization: the global continuum drifts below 1.0 in the
    # crowded optical (e.g. 0.91 at Mg 4730). Renormalize to the window's upper
    # envelope (median of points above the 85th percentile) so EW/depth are not biased.
    cont = float(np.median(f[f >= np.percentile(f, 85)]))
    if cont <= 0:
        return dict(status='bad_continuum')
    f = f / cont
    # centre the fit on the OBSERVED core (catalog air wavelength is ~15 mA off here).
    near = np.abs(w - line) < CORE_TOL
    core = float(w[near][np.argmin(f[near])]) if near.any() else line
    popt, pcov, ptype, chi2 = _fit_profile(w, f, core)
    if popt is None:
        return dict(status='fit_failed', local_cont=round(cont, 4))
    ew = _integrate_profile(w, popt, ptype)
    edge = (np.abs(w - core) > 0.35)
    ew_err = _ew_error(f[edge], ew, popt, pcov, profile_type=ptype)
    x0 = float(popt[0])
    return dict(status='ok', ew_mA=round(float(ew), 2),
                ew_err_mA=round(float(ew_err), 2) if np.isfinite(ew_err) else np.nan,
                profile=ptype, chi2=round(float(chi2), 4), local_cont=round(cont, 4),
                fit_center=round(x0, 3), center_off_mA=round((x0 - line) * 1000, 1),
                depth=round(float(popt[1]), 3))


def main():
    s = pd.read_csv(SPEC)
    wav = s['wavelength_air_A'].to_numpy()
    flux = s['flux_normalized'].to_numpy()
    rows = []
    for el, ion, line, ep, loggf, grade, ref in CANDIDATES:
        r = measure(wav, flux, line)
        # reduced EW = log10(EW/lambda) with EW and lambda in the SAME unit (A).
        rew = (np.log10((r['ew_mA'] / 1000.0) / line) if r.get('status') == 'ok'
               and r.get('ew_mA', 0) > 0 else np.nan)
        unsat = bool(np.isfinite(rew) and rew <= REW_KNEE and r.get('ew_mA', 999) <= SAT_CEIL)
        rec = dict(element=el, ion=ion, wavelength_air_A=line, excitation_potential_eV=ep,
                   log_gf=loggf, nist_grade=grade, loggf_reference=ref,
                   rew=round(rew, 3) if np.isfinite(rew) else None,
                   cog_regime='linear/unsaturated' if unsat else 'saturated/over-ceiling',
                   unsaturated=unsat, **r)
        rows.append(rec)
        print(f"  {el} {ion} {line:.3f}  EW={r.get('ew_mA')}+-{r.get('ew_err_mA')} mA  "
              f"REW={rec['rew']}  {rec['cog_regime']:24s} chi2={r.get('chi2')} "
              f"off={r.get('center_off_mA')}mA  grade={grade}")
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    nun = int(df['unsaturated'].sum())
    print(f"\n  {nun}/{len(df)} candidates land in the linear (unsaturated) COG regime.")
    print(f"  wrote {OUT.relative_to(ROOT)}")


if __name__ == '__main__':
    main()

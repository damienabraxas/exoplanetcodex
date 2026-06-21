"""
pipeline/co_validation_rya390.py
================================
RYA-390 Part B — three-way validation of RYA-373's telluric-conditioned CRIRES+
CO spectrum against the IR reference atlases intaked in Part A.

The three checks (the issue's Part B design):
  1. TELLURIC-REMOVAL  conditioned CRIRES CO  vs  ACE-FTS (telluric-free) —
     did we recover the true solar CO? Self-consistency precedent: the solar CO
     should match the solar reference at 1.0× (the ¹³CO/¹⁸O solar-twin check).
  2. TELLURIC-MODEL     molecfit telluric transmission  vs  Wallace ratio —
     did the telluric MODEL match reality? (needs molecfit `mtrans`, see below.)
  3. CROSS-INSTRUMENT   conditioned CO  vs  NSO photatl (independent reduction).

Plus a RESIDUAL-TELLURIC diagnostic: correlate the conditioned CO against the
telluric reference at v≈0 (tellurics are stationary in the topocentric frame). A
high correlation means residual telluric still dominates the "corrected" product —
the decisive check that the removal actually worked.

The conditioned CO is TOPOCENTRIC and RV-INSUFFICIENT (RYA-373 `RVSTATUS`), so the
reflected-solar velocity is unknown. The solar checks therefore MEASURE the velocity
by cross-correlating against the rest-frame atlas, then score recovery in the aligned
overlap. WAVELENGTHS ARE VACUUM throughout (CRIRES IDP + molecfit `WAVELENGTH_FRAME=VAC`
→ the conditioned `wave_A` is vacuum; the atlases carry `wavelength_vac_A`). Mixing in
air would inject a ~6.3 Å ≈ 82 km/s offset at 2.3 µm.

GAP it closes: RYA-373 stamped `GAP1 = FTS solar IR atlas (RYA-162) absent` on the
conditioned product — Part A's ACE/photatl atlases ARE that missing reference.

LIMITATION (telluric-MODEL check 2): RYA-373 persists only the corrected flux, not the
molecfit transmission `mtrans`. `check_telluric_model` is built and ready but reports
BLOCKED-pending-mtrans on the current products; we instead cross-check the two
INDEPENDENT telluric references (Wallace vs photatl atmospheric) so the reference side
is validated for when RYA-373 persists `mtrans`.

Usage:  python -m pipeline.co_validation_rya390
Out:    data/audit/crires_co_conditioned/rya390_co_validation.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

ROOT = Path(__file__).resolve().parent.parent
ATLAS = ROOT / 'data' / 'solar_reference' / 'ir_atlases'
COND = ROOT / 'data' / 'audit' / 'crires_co_conditioned'

C_KMS = 299792.458
R_CRIRES = 86000.0          # CRIRES+ K-band nominal resolving power (atlas → this R)
V_MAX_KMS = 120.0           # reflected-solar topocentric search range (RYA-370: −17..+33)
V_STEP_KMS = 0.5
V_TELLURIC_KMS = 10.0       # tellurics are stationary (topocentric) → search near 0

# Conditioned-CO snapshots: committed CSV (the repo policy ignores *.fits, so RYA-373's
# 56 KB FITS products are snapshotted to CSV — wave_vac_A,flux_norm,err + a header-comment
# provenance block carrying the original FITS cards: RESTFRM/SPECSYS/RVSTATUS/GAP1/...).
CONDITIONED = {
    'K2192': 'vesta_crires_K_CO_K2192_topocent_PROVISIONAL.csv',
    'K2217': 'vesta_crires_K_CO_K2217_topocent_PROVISIONAL.csv',
}


def _read_csv_provenance(path: Path) -> dict:
    """Pull the `# KEY = value` provenance comments from a snapshot CSV header."""
    meta = {}
    for line in path.read_text().splitlines():
        if not line.startswith('#'):
            break
        if '=' in line:
            k, _, v = line[1:].partition('=')
            meta[k.strip()] = v.strip()
    return meta


# ── loaders (VACUUM wavelength throughout) ─────────────────────────────────────

def load_conditioned_co(setting: str) -> dict:
    """Conditioned CRIRES CO order: VACUUM wavelength + continuum-normalized flux,
    gap/bad pixels (wave≈0, flux≤0, non-finite) dropped, ascending λ."""
    path = COND / CONDITIONED[setting]
    meta = _read_csv_provenance(path)
    df = pd.read_csv(path, comment='#')
    w = df['wave_vac_A'].to_numpy(float)        # CRIRES/molecfit VAC
    f = df['flux_norm'].to_numpy(float)
    good = (w > 1000) & np.isfinite(f) & (f > 0)
    w, f = w[good], f[good]
    o = np.argsort(w)
    return {'setting': setting, 'wave_A': w[o], 'flux': f[o],
            'restfrm': meta.get('RESTFRM') == 'True', 'specsys': meta.get('SPECSYS'),
            'rvstatus': meta.get('RVSTATUS', ''), 'provisional': meta.get('PROVIS') == 'True'}


def load_atlas(name: str, flux_col: str) -> dict:
    """An atlas CO segment: VACUUM wavelength + a chosen flux column, ascending λ."""
    files = {'ace': 'ace_fts_solar_co_4255_4367.csv',
             'photatl': 'nso_photatl_co_4255_4367.csv',
             'wallace': 'wallace_telluric_co_ratio.csv'}
    df = pd.read_csv(ATLAS / files[name]).sort_values('wavelength_vac_A')
    return {'name': name, 'wave_A': df['wavelength_vac_A'].to_numpy(float),
            'flux': df[flux_col].to_numpy(float)}


# ── core ───────────────────────────────────────────────────────────────────────

def _smooth_to_crires(wave: np.ndarray, flux: np.ndarray, R: float = R_CRIRES) -> np.ndarray:
    """Gaussian-smooth a high-res atlas to CRIRES resolution (on its own grid)."""
    dl = np.median(np.diff(wave))
    if dl <= 0:
        return flux
    sigma_px = (np.median(wave) / R / 2.3548) / dl
    return gaussian_filter1d(flux, max(sigma_px, 0.5))


def _overlap(w1, w2):
    return max(w1.min(), w2.min()), min(w1.max(), w2.max())


def measure_velocity(w_obs, f_obs, w_ref, f_ref,
                     vmax=V_MAX_KMS, dv=V_STEP_KMS, smooth=True) -> dict:
    """Cross-correlate observed (topocentric) vs a reference over a velocity grid; the
    peak is the velocity that shifts the reference INTO the observed frame. Reference is
    smoothed to CRIRES R first, then resampled at w_obs/(shifted)."""
    lo, hi = _overlap(w_obs, w_ref)
    m_obs = (w_obs >= lo) & (w_obs <= hi)
    wo, fo = w_obs[m_obs], f_obs[m_obs] - np.nanmean(f_obs[m_obs])
    fr_src = _smooth_to_crires(w_ref, f_ref) if smooth else f_ref
    vs = np.arange(-vmax, vmax + dv, dv)
    corr = np.full(vs.size, np.nan)
    for i, v in enumerate(vs):
        fr = np.interp(wo, w_ref * (1.0 + v / C_KMS), fr_src, left=np.nan, right=np.nan)
        fr = fr - np.nanmean(fr)
        ok = np.isfinite(fo) & np.isfinite(fr)
        if ok.sum() < 50:
            continue
        denom = np.sqrt(np.sum(fo[ok] ** 2) * np.sum(fr[ok] ** 2))
        corr[i] = np.sum(fo[ok] * fr[ok]) / denom if denom > 0 else np.nan
    if np.all(np.isnan(corr)):
        return {'v_kms': np.nan, 'peak_xcorr': np.nan, 'n_overlap_px': int(m_obs.sum()),
                'overlap_vac_A': [round(lo, 2), round(hi, 2)]}
    j = int(np.nanargmax(corr))
    return {'v_kms': float(vs[j]), 'peak_xcorr': float(corr[j]),
            'n_overlap_px': int(m_obs.sum()), 'overlap_vac_A': [round(lo, 2), round(hi, 2)]}


def _aligned_residual(w_obs, f_obs, w_ref, f_ref, v_kms) -> dict:
    """Resample the velocity-aligned, CRIRES-smoothed reference onto the observed grid;
    score recovery: RMS(obs−ref) and absorption-depth (1−flux) correlation."""
    lo, hi = _overlap(w_obs, w_ref)
    m = (w_obs >= lo) & (w_obs <= hi)
    wo, fo = w_obs[m], f_obs[m]
    fr = np.interp(wo, w_ref * (1.0 + v_kms / C_KMS), _smooth_to_crires(w_ref, f_ref),
                   left=np.nan, right=np.nan)
    ok = np.isfinite(fo) & np.isfinite(fr)
    resid_rms = float(np.sqrt(np.nanmean((fo[ok] - fr[ok]) ** 2)))
    do, dr = 1.0 - fo[ok], 1.0 - fr[ok]
    denom = np.sqrt(np.sum((do - do.mean()) ** 2) * np.sum((dr - dr.mean()) ** 2))
    depth_corr = float(np.sum((do - do.mean()) * (dr - dr.mean())) / denom) if denom > 0 else np.nan
    return {'resid_rms': resid_rms, 'depth_corr': depth_corr, 'n_px': int(ok.sum())}


# ── the three checks + the residual-telluric diagnostic ───────────────────────

def _solar_check(co: dict, ref: dict, label: str) -> dict:
    vel = measure_velocity(co['wave_A'], co['flux'], ref['wave_A'], ref['flux'])
    rec = (_aligned_residual(co['wave_A'], co['flux'], ref['wave_A'], ref['flux'], vel['v_kms'])
           if np.isfinite(vel['v_kms']) else {'resid_rms': np.nan, 'depth_corr': np.nan, 'n_px': 0})
    return {'reference': label, **vel, **rec}


def check_telluric_removal(co: dict, ace: dict) -> dict:
    """Check 1: conditioned CO vs ACE-FTS (telluric-free solar)."""
    return _solar_check(co, ace, 'ACE-FTS (telluric-free solar)')


def check_cross_instrument(co: dict, photatl_solar: dict) -> dict:
    """Check 3: conditioned CO vs photatl solar column (independent ground reduction)."""
    return _solar_check(co, photatl_solar, 'NSO photatl (terrestrial solar)')


def diagnose_residual_telluric(co: dict, telluric_ref: dict) -> dict:
    """Decisive diagnostic: correlate the conditioned CO against the TELLURIC reference
    near v≈0 (tellurics are stationary in the topocentric frame). High xcorr ⇒ residual
    telluric still dominates the 'corrected' product."""
    vel = measure_velocity(co['wave_A'], co['flux'], telluric_ref['wave_A'],
                           telluric_ref['flux'], vmax=V_TELLURIC_KMS)
    return {'reference': f"{telluric_ref['name']} (telluric)", **vel}


def check_telluric_model(mtrans_wave_A, mtrans, wallace: dict) -> dict:
    """Check 2: molecfit telluric transmission vs the Wallace telluric ratio. Ready for
    when RYA-373 persists `mtrans`; mtrans is None on the current products → BLOCKED."""
    if mtrans is None or mtrans_wave_A is None:
        return {'reference': 'Wallace telluric near-IR atlas',
                'status': 'BLOCKED — molecfit transmission (mtrans) not persisted by '
                          'RYA-373 (only corrected flux is saved). RECOMMEND RYA-373 add '
                          'the BEST_FIT_MODEL.mtrans column to the conditioned FITS; then '
                          'this check runs directly.', 'verdict': 'BLOCKED'}
    wgrid = np.asarray(mtrans_wave_A)
    lo, hi = _overlap(wgrid, wallace['wave_A'])
    m = (wgrid >= lo) & (wgrid <= hi)
    wr = np.interp(wgrid[m], wallace['wave_A'], _smooth_to_crires(wallace['wave_A'], wallace['flux']))
    mt = np.asarray(mtrans)[m]
    ok = np.isfinite(wr) & np.isfinite(mt)
    rms = float(np.sqrt(np.nanmean((mt[ok] - wr[ok]) ** 2)))
    return {'reference': 'Wallace telluric near-IR atlas', 'status': 'RAN', 'resid_rms': rms,
            'n_px': int(ok.sum()),
            'verdict': 'PASS — telluric model matches' if rms <= 0.05 else 'FAIL — model–reality mismatch'}


def crosscheck_telluric_references(wallace: dict, photatl_atm: dict) -> dict:
    """Validate the telluric REFERENCE side: do the two INDEPENDENT telluric atlases
    (Wallace ratio vs photatl atmospheric) agree in their overlap? Stationary → no shift."""
    lo, hi = _overlap(wallace['wave_A'], photatl_atm['wave_A'])
    if hi <= lo:
        return {'status': 'NO OVERLAP'}
    m = (wallace['wave_A'] >= lo) & (wallace['wave_A'] <= hi)
    pa = np.interp(wallace['wave_A'][m], photatl_atm['wave_A'], photatl_atm['flux'])
    wl = wallace['flux'][m]
    ok = np.isfinite(pa) & np.isfinite(wl)
    do, dr = 1 - wl[ok], 1 - pa[ok]
    denom = np.sqrt(np.sum((do - do.mean()) ** 2) * np.sum((dr - dr.mean()) ** 2))
    corr = float(np.sum((do - do.mean()) * (dr - dr.mean())) / denom) if denom > 0 else np.nan
    return {'status': 'RAN', 'overlap_vac_A': [round(lo, 2), round(hi, 2)],
            'telluric_depth_corr': corr, 'n_px': int(ok.sum()),
            'verdict': 'PASS — independent telluric refs agree' if corr >= 0.6
                       else 'WEAK — telluric refs differ (coverage/epoch)'}


def _setting_verdict(solar_best, tell_corr) -> str:
    """Synthesize the per-setting verdict from solar recovery vs residual telluric."""
    if np.isfinite(tell_corr) and tell_corr >= 0.6 and (not np.isfinite(solar_best) or solar_best < 0.45):
        return ('TELLURIC-DOMINATED — residual telluric still dominates the corrected '
                f'product (xcorr {tell_corr:.2f} with the telluric atlas at v≈0); solar CO '
                'NOT recovered. Telluric removal incomplete — re-run RYA-373.')
    if np.isfinite(solar_best) and solar_best >= 0.7:
        return 'PASS — solar CO recovered'
    if np.isfinite(solar_best) and solar_best >= 0.45:
        return 'PARTIAL — solar CO pattern present but imperfect'
    return 'INCONCLUSIVE'


def run() -> dict:
    ace = load_atlas('ace', 'intensity')
    ph_solar = load_atlas('photatl', 'solar')
    ph_atm = load_atlas('photatl', 'atmospheric')
    wallace = load_atlas('wallace', 'telluric_ratio')

    report = {'ticket': 'RYA-390 Part B', 'co_band_vac_A': [22899.0, 23502.0],
              'wavelength_convention': 'VACUUM (CRIRES/molecfit VAC; atlas wavelength_vac_A)',
              'note_gap1': 'closes RYA-373 GAP1 (FTS solar IR atlas absent)', 'settings': {}}
    for setting in CONDITIONED:
        co = load_conditioned_co(setting)
        c1 = check_telluric_removal(co, ace)
        c3 = check_cross_instrument(co, ph_solar)
        rt = diagnose_residual_telluric(co, ph_atm)
        solar_best = np.nanmax([c1['depth_corr'], c3['depth_corr']])
        report['settings'][setting] = {
            'frame': co['specsys'], 'restfrm': co['restfrm'],
            'rvstatus': co['rvstatus'][:60], 'provisional': co['provisional'],
            'check1_telluric_removal_vs_ACE': c1,
            'check3_cross_instrument_vs_photatl': c3,
            'check2_telluric_model_vs_Wallace': check_telluric_model(None, None, wallace),
            'diagnostic_residual_telluric_vs_photatl_atm': rt,
            'verdict': _setting_verdict(solar_best, rt['peak_xcorr']),
        }
    report['telluric_reference_crosscheck'] = crosscheck_telluric_references(wallace, ph_atm)

    out = COND / 'rya390_co_validation.json'
    out.write_text(json.dumps(report, indent=2))

    print("\n========= RYA-390 Part B — three-way CO validation (VACUUM λ) =========")
    print("  conditioned CO is PROVISIONAL / topocentric / RV-insufficient; solar checks")
    print("  MEASURE the reflected-solar velocity by cross-correlation, then score recovery.\n")
    for s, r in report['settings'].items():
        c1, c3, rt = (r['check1_telluric_removal_vs_ACE'], r['check3_cross_instrument_vs_photatl'],
                      r['diagnostic_residual_telluric_vs_photatl_atm'])
        print(f"  [{s}]  frame={r['frame']}")
        print(f"    1 solar  vs ACE     : v={c1['v_kms']:+6.1f}  xcorr={c1['peak_xcorr']:.3f}  "
              f"depth_corr={c1['depth_corr']:.3f}  rms={c1['resid_rms']:.3f}")
        print(f"    3 solar  vs photatl : v={c3['v_kms']:+6.1f}  xcorr={c3['peak_xcorr']:.3f}  "
              f"depth_corr={c3['depth_corr']:.3f}  rms={c3['resid_rms']:.3f}")
        print(f"    ! residual telluric : v={rt['v_kms']:+6.1f}  xcorr={rt['peak_xcorr']:.3f}  (v≈0 ⇒ telluric)")
        print(f"    2 telluric-model    : {r['check2_telluric_model_vs_Wallace']['verdict']}")
        print(f"    → {r['verdict']}\n")
    tc = report['telluric_reference_crosscheck']
    print(f"  telluric-ref cross-check (Wallace vs photatl-atm): {tc.get('verdict', tc['status'])}"
          + (f"  (depth_corr={tc['telluric_depth_corr']:.3f})" if tc.get('telluric_depth_corr') is not None else ""))
    print(f"\n  [out] {out}")
    return report


if __name__ == '__main__':
    run()

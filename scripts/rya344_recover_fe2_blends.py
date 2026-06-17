"""
scripts/rya344_recover_fe2_blends.py
====================================
RYA-344 — recover the 5 χ²ᵣ-failing solar Fe II lines by adding the missing
in-window blend absent from the synth line list (RYA-343 found Fe II is 63–78%
of each feature but the model under-models it 0.65–0.79× obs → a missing absorber,
not the Fe II line's own strength).

Two phases, driven by the current state of the synth line list (atomic_lines.tsv):

  --identify  : for each of the 5 windows, fit the Fe II abundance in flux space
                (same machinery as RYA-287/342), compute the residual (obs − synth)
                at the best fit, locate the missing-flux centroid, and search
                linelist_solar.csv for candidate absorbers in that sub-window NOT
                already in the synth line list. Ranks by central_depth × proximity.
                READ-ONLY. Writes data/results/rya344_candidates.csv.

  --verify    : re-fit the 5 windows against the CURRENT synth line list and report
                χ²ᵣ + Fe II abundance per line. Run before AND after the additions
                to get the before/after deltas. Writes rya344_verify.csv.

Source of absorbers = linelist_solar.csv (our VALD3 extraction; loggf_source carried
inline). The synth line-list edit itself (append to atomic_lines.tsv) is done by a
separate explicit step once a candidate is chosen — this script never writes the
line list. NO guessing: if no linelist_solar candidate explains a residual, the
window is reported NOT-recovered for manual VALD extraction.
"""
import sys, argparse
from pathlib import Path
import numpy as np, pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))

from config.constants import STAR_PARAMS, ISPEC_DIR
from pipeline.species import species_key
import pipeline.abundances_derive as ad

_OUT = _REPO / 'data' / 'results'; _OUT.mkdir(parents=True, exist_ok=True)
_LINELIST_SOLAR = _REPO / 'data' / 'linelists' / 'linelist_solar.csv'
_TMP = '/tmp/ispec_rya344'
Path(_TMP).mkdir(parents=True, exist_ok=True)   # Turbospectrum needs the tmp dir to exist

# RYA-343 recoverable set: (wave_air_A, observed EW mÅ, RYA-342 best-fit a_synth)
TARGETS = [
    (5234.623, 120.7, 7.494),
    (6084.102,  28.5, 7.634),
    (6247.557,  62.9, 7.555),
    (5991.371,  44.3, 7.586),
    (6456.380,  85.4, 7.685),
]
CORE_HW_A   = 0.18    # Å — search the line core ± this for the missing absorber
CLEAN_POOL_A = 7.50   # RYA-342 clean-pool Fe II consensus (n=5) → abundance-consistency ref


def _solar_context():
    """Set up the shared solar flux-synthesis context (params/atmos/linelist/broadening)."""
    sp = STAR_PARAMS['solar']
    teff, logg, vturb = float(sp['teff']), float(sp['logg']), float(sp['xi'])
    feh = float(sp['feh_ref'])
    obs_wave_nm, obs_flux = ad._load_observed_spectrum('solar')
    atmosphere = ad._load_atmosphere(teff, logg, feh, vturb)
    import ispec
    linelist, isotopes, chem = ad._load_synth_resources()
    solar_abund = ispec.read_solar_abundances(ad._ISPEC_SOLAR_ABUND_FILE)
    tmp = ispec.create_free_abundances_structure(['Fe'], chem, solar_abund)
    a_code = int(tmp['code'][0]); a_solar = float(tmp['Abund'][0]) + 12.036
    R, vmac, vsini, _ = ad._resolve_broadening('solar')
    return dict(teff=teff, logg=logg, feh=feh, vturb=vturb, obs_wave_nm=obs_wave_nm,
                obs_flux=obs_flux, atmosphere=atmosphere, linelist=linelist,
                isotopes=isotopes, solar_abund=solar_abund, a_code=a_code,
                a_solar=a_solar, R=R, vmac=vmac, vsini=vsini)


def _fit_line(ctx, wave_A, ew_obs):
    """Flux-space Fe II fit for one window → (a_synth, red_chi2, residual dict)."""
    wave_nm = wave_A / 10.0
    wbase, wtop = ad._wingwide_window_nm(wave_nm, ew_obs)
    a_lo, a_hi = max(ctx['a_solar'] - 3.0, 1.0), ctx['a_solar'] + 5.0
    r = ad._fit_synth_flux(
        ctx['obs_wave_nm'], ctx['obs_flux'], ctx['atmosphere'],
        ctx['teff'], ctx['logg'], ctx['feh'], ctx['vturb'],
        ctx['linelist'], ctx['isotopes'], ctx['solar_abund'], 'Fe', ctx['a_code'],
        wbase, wtop, a_lo, a_hi, ctx['R'], ctx['vmac'], ctx['vsini'],
        tmp_dir=_TMP)
    # residual (obs − synth) at the best fit, on the observed grid in-window
    resid = None
    if np.isfinite(r['A_X']):
        sw = np.arange(wbase, wtop + 1e-4, 0.0002)
        sf = ad._synth_flux_at_abund(
            sw, ctx['atmosphere'], ctx['teff'], ctx['logg'], ctx['feh'], ctx['vturb'],
            ctx['linelist'], ctx['isotopes'], ctx['solar_abund'], 'Fe', ctx['a_code'],
            float(r['A_X']), R=ctx['R'], macroturbulence=ctx['vmac'], vsini=ctx['vsini'],
            tmp_dir=_TMP)
        from scipy.interpolate import interp1d
        mask = (ctx['obs_wave_nm'] >= wbase) & (ctx['obs_wave_nm'] <= wtop)
        ow = ctx['obs_wave_nm'][mask] * 10.0          # Å
        of = ctx['obs_flux'][mask]
        sfi = interp1d(sw * 10.0, sf, bounds_error=False, fill_value=1.0)(ow)
        d = of - sfi                                   # <0 where obs deeper (missing absorption)
        resid = dict(wave_A=ow, dflux=d)
    return float(r['A_X']), float(r['red_chi2']), r['status'], resid


def _missing_centroid(resid, wave_A):
    """Flux-weighted centroid of the missing absorption (obs deeper than model)
    within the core window, and its peak depth."""
    w, d = resid['wave_A'], resid['dflux']
    core = (np.abs(w - wave_A) <= CORE_HW_A) & (d < 0)
    if not core.any():
        return np.nan, 0.0
    wt = -d[core]
    return float(np.sum(w[core] * wt) / np.sum(wt)), float(wt.max())


def _in_synth_list(ctx, sp_key, wave_A, tol=0.05):
    """Is a line of this species already in atomic_lines.tsv within tol Å?"""
    ll = ctx['linelist']
    el = np.array([str(e) for e in ll['element']])
    wv = np.array(ll['wave_A'], dtype=float)
    for i in np.where(np.abs(wv - wave_A) <= tol)[0]:
        try:
            if species_key(el[i]) == sp_key:
                return True, float(wv[i])
        except ValueError:
            continue
    return False, np.nan


def identify(ctx):
    sol = pd.read_csv(_LINELIST_SOLAR)
    rows = []
    for wave_A, ew_obs, _a342 in TARGETS:
        a, chi, status, resid = _fit_line(ctx, wave_A, ew_obs)
        cen, depth = (_missing_centroid(resid, wave_A) if resid else (np.nan, 0.0))
        print(f"\n── {wave_A:.3f} Å  a_synth={a:.3f}  χ²ᵣ={chi:.1f}  "
              f"missing-flux centroid={cen:.3f} Å (peak {depth:.3f})")
        # candidate absorbers in the core, from our VALD3 extraction, NOT the Fe II target
        win = sol[(sol.wavelength_air_A >= wave_A - CORE_HW_A) &
                  (sol.wavelength_air_A <= wave_A + CORE_HW_A)].copy()
        win = win[~((win.element == 'Fe') & (win.ion == 'II') &
                    ((win.wavelength_air_A - wave_A).abs() < 0.05))]   # exclude target
        cand_rows = []
        for _, lr in win.iterrows():
            sp = species_key(str(lr.element), str(lr.ion))
            present, w_in = _in_synth_list(ctx, sp, float(lr.wavelength_air_A))
            dist = abs(float(lr.wavelength_air_A) - (cen if np.isfinite(cen) else wave_A))
            cd = float(lr.central_depth) if pd.notna(lr.central_depth) else 0.0
            score = cd / (1.0 + 10.0 * dist)          # depth, penalized by distance to residual
            cand_rows.append(dict(
                target=wave_A, cand_wave=round(float(lr.wavelength_air_A), 3),
                species=f"{lr.element} {lr.ion}", central_depth=round(cd, 3),
                log_gf=lr.log_gf, ep_eV=lr.excitation_potential_eV,
                loggf_source=lr.loggf_source, in_synth_list=present,
                synth_wave=round(w_in, 3) if present else np.nan,
                dist_to_resid=round(dist, 3), score=round(score, 4)))
        _cols = ['target', 'cand_wave', 'species', 'central_depth', 'log_gf',
                 'ep_eV', 'loggf_source', 'in_synth_list', 'synth_wave',
                 'dist_to_resid', 'score']
        cdf = (pd.DataFrame(cand_rows, columns=_cols).sort_values('score', ascending=False)
               if cand_rows else pd.DataFrame(columns=_cols))
        # the recovery candidates = strong lines MISSING from the synth list
        missing = cdf[(~cdf.in_synth_list.astype(bool)) & (cdf.central_depth > 0.01)]
        print("  Top candidates (missing from synth list, by score):")
        if len(missing):
            print(missing.head(5)[['cand_wave', 'species', 'central_depth',
                                    'log_gf', 'loggf_source', 'dist_to_resid',
                                    'score']].to_string(index=False))
        else:
            print("  *** NONE missing — no linelist_solar absorber explains the residual ***")
            print("  → FLAG for manual VALD3 extraction (do not guess).")
        cdf['fit_a_synth'] = round(a, 3); cdf['fit_chi2r'] = round(chi, 1)
        cdf['resid_centroid'] = round(cen, 3) if np.isfinite(cen) else np.nan
        rows.append(cdf)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(_OUT / 'rya344_candidates.csv', index=False)
    print(f"\nWrote {(_OUT / 'rya344_candidates.csv').relative_to(_REPO)}")


def verify(ctx):
    rows = []
    for wave_A, ew_obs, a342 in TARGETS:
        a, chi, status, _ = _fit_line(ctx, wave_A, ew_obs)
        recovered = (chi < 10.0) and (abs(a - CLEAN_POOL_A) <= 0.10)
        rows.append(dict(wave_A=wave_A, a_synth=round(a, 3), red_chi2=round(chi, 1),
                         status=status, abund_consistent=abs(a - CLEAN_POOL_A) <= 0.10,
                         chi2_pass=chi < 10.0, recovered=recovered))
        print(f"  {wave_A:.3f}  a_synth={a:.3f}  χ²ᵣ={chi:.1f}  "
              f"recovered={'Y' if recovered else 'N'}")
    out = pd.DataFrame(rows)
    out.to_csv(_OUT / 'rya344_verify.csv', index=False)
    print(f"\nRecovered {int(out.recovered.sum())}/5 | clean pool 5 → {5 + int(out.recovered.sum())}")
    print(f"Wrote {(_OUT / 'rya344_verify.csv').relative_to(_REPO)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--identify', action='store_true')
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    ctx = _solar_context()
    if args.identify:
        identify(ctx)
    if args.verify:
        verify(ctx)
    if not (args.identify or args.verify):
        ap.error("pass --identify and/or --verify")


if __name__ == '__main__':
    main()

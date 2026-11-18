"""Direct fit of solar lines with OpenMULTI3D spectra: profile fit and EW inversion, LTE + NLTE.

Inputs: the wide run (run_openmulti3d_cog.py with a background line list, full spectra per
abundance), the O-only curve of growth (same driver, no line list), the spec of each, and the
ORCH holdings. Observed flux comes through the ORCH loader (measure_band_ew.load_window_ex)
and the continuum through the ORCH standing rule (local_continuum.fit_model_guided with the
RYA-1350 neighbourhood check), with the OpenMULTI3D synthesis as the reference model.

Profile route: per line, chi^2 over +/- half_width_A; free A (linear interpolation of the
model flux between grid abundances), Gaussian broadening (macroturbulence + instrument, km/s)
and wavelength shift.
EW route: EW_obs over the COG integration window minus the synthetic blend EW (wide run minus
O-only run at the profile-route abundance), inverted on the O-only LTE and NLTE COGs.

Run with the tumult venv (CODEX_OPENMULTI3D_PY) and PYTHONPATH=<engine>:<engine>/HOP:<repo>.
"""
import argparse, hashlib, json, os, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'scripts'))
C_KMS = 299792.458


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_grid(wide, spec):
    import tumult as tm
    cache = wide / 'spectra_cache.npz'
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        return d['lam'], d['grid'].item()
    grid, lam = {}, None
    for atom in spec['atoms']:
        for kind in ['nlte', 'lte']:
            rows = []
            for A in spec['abundances']:
                run = tm.read(str(wide / 'runs' / f'{atom}_A{A:.3f}' / 'data' / 'run'))
                l = np.asarray(run.lam, float)
                lam = l if lam is None else lam
                if l.shape != lam.shape or np.max(abs(l - lam)) > 1e-6:
                    raise ValueError('wavelength axis differs between runs')
                f = run.flux if kind == 'nlte' else run.flux_lte
                rows.append(np.asarray(f, float) / np.asarray(run.flux_cnt, float))
            grid[(atom, kind)] = np.array(rows)
    np.savez_compressed(cache, lam=lam, grid=np.array(grid, dtype=object))
    return lam, grid


def model_at(lam, rows, abund, A):
    A = float(np.clip(A, abund[0], abund[-1]))
    i = min(np.searchsorted(abund, A) - 1, len(abund) - 2); i = max(i, 0)
    t = (A - abund[i]) / (abund[i + 1] - abund[i])
    return (1 - t) * rows[i] + t * rows[i + 1]


def broaden(lam, flux, fwhm_kms):
    """Gaussian in velocity on a log-uniform resample."""
    from scipy.ndimage import gaussian_filter1d
    if fwhm_kms <= 0:
        return flux
    dv = C_KMS * np.median(np.diff(lam) / lam[:-1])
    ll = np.exp(np.arange(np.log(lam[0]), np.log(lam[-1]), dv / C_KMS))
    f = np.interp(ll, lam, flux)
    f = gaussian_filter1d(f, fwhm_kms / 2.35482 / dv, mode='nearest')
    return np.interp(lam, ll, f)


def ew_mA(w, depth, lo, hi):
    m = (w >= lo) & (w <= hi)
    return 1e3 * float(np.trapezoid(depth[m], w[m]))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--wide', type=Path, required=True)
    p.add_argument('--wide-spec', type=Path, required=True)
    p.add_argument('--cog', type=Path, required=True)
    p.add_argument('--holdings', nargs='+', required=True, help='instrument:holding:resolving_power')
    p.add_argument('--half-width-A', type=float, default=1.1)
    p.add_argument('--ew-half-A', type=float, default=0.6)
    p.add_argument('--cog-spec', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    from scipy.optimize import least_squares
    import measure_band_ew as H
    from pipeline import local_continuum as LC

    spec = json.loads(a.wide_spec.read_text())
    cog_spec = json.loads(a.cog_spec.read_text())
    abund = np.array(spec['abundances'], float)
    cabund = np.array(cog_spec['abundances'], float)
    lam, grid = load_grid(a.wide, spec)
    clam, cgrid = load_grid(a.cog, cog_spec)
    wide_lo, wide_hi = spec['windows_air_A'][0]
    centre_all = 0.5 * (min(l['wavelength_air_A'] for l in spec['lines']) + max(l['wavelength_air_A'] for l in spec['lines']))
    rows, plots = [], {}
    for hs in a.holdings:
        inst, hold, R = hs.split(':')
        win = H.load_window_ex(inst, centre_all, 0.5 * (wide_hi - wide_lo) - 2.0, holding=hold)
        if not win.holding.pre_normalised:
            raise ValueError(f'{hold} is not pre-normalised; this fit does not place a band continuum')
        W, F = np.asarray(win.wave, float), np.asarray(win.flux, float)
        for atom in spec['atoms']:
            for kind in ['nlte', 'lte']:
                rows_k, crows_k = grid[(atom, kind)], cgrid[(atom, kind)]
                for line in spec['lines']:
                    c0 = line['wavelength_air_A']
                    env = (W >= c0 - LC.ENV_HALF_WIDTH_A) & (W <= c0 + LC.ENV_HALF_WIDTH_A)
                    legs = {}
                    # (a) RYA-1232 nominal: synthesis reference (here the OM synthesis + SIU list)
                    ref = broaden(lam, model_at(lam, rows_k, abund, 8.80), 4.0)
                    check = LC.neighbourhood_check(W, F, c0, lam, ref, exclude_half_width_A=a.half_width_A)
                    legs['synth_ref'] = LC.fit_model_guided(W[env], F[env], c0, lam, ref, exclude_half_width_A=a.half_width_A,
                                                            reference='OpenMULTI3D synthesis', corroborate=check)[0]
                    # (b) RYA-587 contref leg: the IAG atlas as reference
                    LC.CONT_REF = 'iag'
                    iw, iff = LC.solar_atlas_reference(c0 - LC.NEIGH_A - 6, c0 + LC.NEIGH_A + 6, float(R))
                    LC.CONT_REF = 'synthesis'
                    check = LC.neighbourhood_check(W, F, c0, iw, iff, exclude_half_width_A=a.half_width_A)
                    legs['iag_ref'] = LC.fit_model_guided(W[env], F[env], c0, iw, iff, exclude_half_width_A=a.half_width_A,
                                                          reference='IAG solar atlas', corroborate=check)[0]
                    # (c) local upper envelope (Amarsi+2020 local maxima analogue)
                    legs['envelope'] = LC.fit(W, F, c0)[0]
                    for leg, rec in legs.items():
                        level = rec.level_at_centre if rec.applied else 1.0
                        slope = rec.slope_per_A if rec.applied else 0.0
                        Fn = F / (level + slope * (W - c0))
                        m = np.abs(W - c0) <= a.half_width_A
                        wv, fv = W[m], Fn[m]
                        sub = (lam >= c0 - a.half_width_A - 3) & (lam <= c0 + a.half_width_A + 3)

                        def resid(x):
                            A, fwhm, shift = x
                            mod = broaden(lam[sub], model_at(lam, rows_k, abund, A)[sub], fwhm)
                            return np.interp(wv - shift, lam[sub], mod) - fv
                        best = None
                        for A0 in [8.7, 8.9]:
                            r = least_squares(resid, [A0, 3.5, 0.0], bounds=([abund[0], 0.0, -0.08], [abund[-1], 12.0, 0.08]),
                                              x_scale=[0.05, 0.5, 0.005], diff_step=[1e-3, 1e-2, 1e-2])
                            best = r if best is None or r.cost < best.cost else best
                        A_fit, fwhm, shift = best.x
                        rms = float(np.sqrt(np.mean(best.fun ** 2)))
                        # EW route: +/- ew_half around the (shifted) centre on obs, wide synthesis
                        # (blend) and O-only synthesis, all with the fitted broadening.
                        lo, hi = c0 + shift - a.ew_half_A, c0 + shift + a.ew_half_A
                        ew_obs = ew_mA(W, 1 - Fn, lo, hi)
                        csub = (clam >= lo - 3) & (clam <= hi + 3)
                        ew_wide = ew_mA(lam[sub], 1 - broaden(lam[sub], model_at(lam, rows_k, abund, A_fit)[sub], fwhm), lo, hi)
                        cog_ew = np.array([ew_mA(clam[csub], 1 - broaden(clam[csub], crows_k[i][csub], fwhm), lo, hi) for i in range(len(cabund))])
                        ew_o_only = float(np.interp(A_fit, cabund, cog_ew))
                        blend = ew_wide - ew_o_only
                        ew_o = ew_obs - blend
                        if not cog_ew[0] < ew_o < cog_ew[-1]:
                            A_ew = float('nan')
                        else:
                            A_ew = float(np.interp(np.log10(ew_o), np.log10(cog_ew), cabund))
                        rows.append(dict(holding=hold, atom=atom, treatment=f'1D-{kind.upper()}', continuum_leg=leg, line=line['name'],
                                         wavelength_air_A=c0, A_profile=round(A_fit, 4), fwhm_kms=round(fwhm, 3),
                                         shift_A=round(shift, 4), rms=round(rms, 5),
                                         at_grid_edge=bool(A_fit <= abund[0] + 1e-3 or A_fit >= abund[-1] - 1e-3),
                                         A_ew=round(A_ew, 4), ew_obs_mA=round(ew_obs, 2), blend_ew_mA=round(blend, 2),
                                         continuum_level=round(rec.level_at_centre, 5), continuum_applied=bool(rec.applied),
                                         continuum_reason=rec.reason[:120]))
                        if atom == list(spec['atoms'])[0] and kind == 'nlte' and leg == 'envelope':
                            mod = broaden(lam[sub], model_at(lam, rows_k, abund, A_fit)[sub], fwhm)
                            plots[(hold, line['name'])] = (wv, fv, np.interp(wv - shift, lam[sub], mod))
                        print(hold, atom, kind, leg, line['name'], f'A_prof={A_fit:.3f} A_ew={A_ew:.3f} ew={ew_obs:.1f}-{blend:.1f} fwhm={fwhm:.2f} rms={rms:.4f} cont={rec.level_at_centre:.4f}/{rec.applied}', flush=True)
    t = pd.DataFrame(rows)
    a.output.mkdir(parents=True, exist_ok=False)
    t.to_csv(a.output / 'fits.csv', index=False)
    summ = (t.groupby(['continuum_leg', 'holding', 'atom', 'treatment'])[['A_profile', 'A_ew']].agg(['mean', 'std', 'count']).round(4))
    summ.to_csv(a.output / 'summary.csv')
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    holds = sorted({k[0] for k in plots}); names = [l['name'] for l in spec['lines']]
    fig, ax = plt.subplots(len(holds), len(names), figsize=(4 * len(names), 2.8 * len(holds)), squeeze=False)
    for i, h in enumerate(holds):
        for j, n in enumerate(names):
            wv, fv, mv = plots[(h, n)]
            ax[i, j].plot(wv, fv, 'k.', ms=2, label='observed'); ax[i, j].plot(wv, mv, 'r-', lw=1, label='OM NLTE fit')
            ax[i, j].plot(wv, fv - mv + 0.4, 'b-', lw=0.6, label='resid+0.4'); ax[i, j].set_title(f'{h} {n}', fontsize=8)
    ax[0, 0].legend(fontsize=6); fig.tight_layout(); fig.savefig(a.output / 'profiles.png', dpi=120)
    (a.output / 'report.json').write_text(json.dumps(dict(
        schema='codex.openmulti3d_line_fit/1', generator_sha256=sha(__file__), wide_spec_sha256=sha(a.wide_spec),
        wide_report=json.loads((a.wide / 'report.json').read_text()), cog_report=json.loads((a.cog / 'report.json').read_text()),
        holdings=a.holdings, half_width_A=a.half_width_A,
        continuum_legs=dict(synth_ref='RYA-1232 nominal: fit_model_guided + neighbourhood_check, OM synthesis A=8.80 4 km/s', iag_ref='RYA-587 contref leg: IAG atlas reference', envelope='local_continuum.fit upper envelope'), ew_half_A=a.ew_half_A,
        profile_model='OM flux/cont, linear in A between grid points, Gaussian FWHM free (km/s), shift free'), indent=2) + '\n')
    print(summ.to_string())


if __name__ == '__main__':
    main()

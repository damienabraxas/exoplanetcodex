"""RYA-1220 diagnostic N I profile reproduction on Sirius; no publication.

Loads the existing CNO engine from --engine-root, recording source hashes.
All synthesized windows span the line centre plus/minus 6 A; only the central mask
contributes to the fit. Wide-line wings beyond that buffer remain a limitation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--engine-root', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--holding', default='solar_iag')
    ap.add_argument('--instrument', default='iag_fts_solar_atlas')
    ap.add_argument('--line', type=float, required=True)
    ap.add_argument('--smoke', action='store_true')
    ap.add_argument('--quick', action='store_true', help='Nominal narrow/wide masks only')
    ap.add_argument('--code', default='turbospectrum',
                    choices=('turbospectrum','spectrum','moog','moog-scat','synthe','sme'),
                    help='iSpec synthesis backend; molecular opacity is enabled only for Turbospectrum')
    ap.add_argument('--amarsi-reference-10108', action='store_true',
                    help='Diagnostic-only Table 1 transition overlay; never modifies canonical gf')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(args.engine_root.resolve()))
    sys.path.insert(0, str(args.engine_root.resolve() / 'scripts'))
    import numpy as np
    from scipy.optimize import minimize_scalar
    from pipeline import cno_synthesis as c
    from config.constants import get_star_params
    from config.synth_bands import SYNTH_BANDS
    from pipeline.band_policy import resolve
    from pipeline.nearuv_synth import gf_provenance
    from measure_band_ew import load_window_ex

    wave = args.line
    band = resolve(wave).name
    config = SYNTH_BANDS[band]
    llpath = config.linelist
    gf = gf_provenance(wave-6, wave+6)
    ll, iso, chem = c._load_synth_resources(str(llpath), apply_canonical_gf=gf['apply_canonical_gf'])
    overlay = None
    if args.amarsi_reference_10108:
        if abs(wave-10108.90) > 1e-6:
            raise ValueError('Reference overlay is restricted to 10108.90 A')
        if np.any((np.abs(ll['wave_A']-wave)<.05) & (ll['element']=='N 1')):
            raise ValueError('Existing N I line must be adjudicated, not duplicated')
        # Published Table 1, not a solar-adjusted gf. These are the fields the
        # installed iSpec Turbospectrum writer consumes (ispec.lines).
        row = np.zeros(1, dtype=ll.dtype)
        values = dict(element='N 1', wave_A=wave, wave_nm=wave/10, loggf=.444,
                      lower_state_eV=11.753, upper_state_eV=12.979,
                      lower_j=1.5, upper_j=2.5, upper_g=6., rad=7.773,
                      turbospectrum_rad=10**7.773, turbospectrum_fdamp=750.2671,
                      waals=750.2671, lower_orbital_type='p', upper_orbital_type='d',
                      molecule='F', ion=1, turbospectrum_species='7.000000',
                      turbospectrum_support='T', nlte='F',
                      nlte_label_low='none', nlte_label_up='none')
        for name, value in values.items():
            row[name] = value
        ll = np.concatenate((ll, row))
        overlay = {'source': 'Amarsi2020 Table 1; DOI 10.1051/0004-6361/202037890',
                   'lower': '3p 4D* J=3/2', 'upper': '3d 4F J=5/2',
                   'fields': values, 'canonical_changed': False,
                   'limitation': 'LTE diagnostic; TS writer rounds ABO alpha to 0.267; no NLTE level indices'}
    targets = ll[(np.abs(ll['wave_A']-wave)<.05) & (ll['element']=='N 1')]
    if not len(targets):
        raise RuntimeError(f'No N I transition within 0.05 A of {wave}; abundance fit refused')
    sab = c.ispec.read_solar_abundances(c._ISPEC_SOLAR_ABUND_FILE)
    codes = c._atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)
    sp = get_star_params('solar')
    params = dict(teff_K=float(sp['teff']), logg=float(sp['logg']),
                  feh=float(sp['feh_ref']), vturb_kms=float(sp.get('xi', 1.)))
    state = c._seed_abundances('solar', params, codes,
                              c._solar_A(('C', 'N', 'O', 'Ni'), chem, sab), params['feh'])
    atm = c._load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    win = load_window_ex(args.instrument, wave, 6., holding=args.holding)
    ow, of = np.asarray(win.wave)/10, np.asarray(win.flux)
    valid = np.isfinite(ow) & np.isfinite(of)
    ow, of = ow[valid], of[valid]
    _, vmac, vsini, fit_vmac = c._resolve_broadening('solar')
    if fit_vmac:
        raise RuntimeError('vmac needs fitting; fixed guess refused')
    # Same instrument-catalog R used by the established band-product route.
    with (args.engine_root/'data/catalog/instrument_catalog.csv').open() as f:
        instruments = list(csv.DictReader(f))
    inst = [r for r in instruments if r['instrument_id'] == args.instrument]
    if len(inst) != 1:
        raise RuntimeError('instrument resolution cannot be resolved uniquely')
    resolution = float(inst[0]['resolving_power_max'])
    broad = (float(resolution), float(vmac), float(vsini))
    sw = np.arange((wave-6)/10, (wave+6)/10, c._WSTEP_NM)
    tmp = args.out / 'tmp'; tmp.mkdir()
    cache = {}

    def synth(a, carbon=None, oxygen=None, xi=None):
        key = (round(float(a), 7), carbon, oxygen, xi)
        if key not in cache:
            st = {**state, 'N': float(a)}
            if carbon is not None: st['C'] = carbon
            if oxygen is not None: st['O'] = oxygen
            pp = dict(params)
            if xi is not None: pp['vturb_kms'] = xi
            R, vmac, vsini = broad
            sf = c.ispec.generate_spectrum(
                sw, atm, float(pp['teff_K']), float(pp['logg']), float(pp['feh']), 0.0,
                ll, iso, sab, c._fixed_ab(st, codes),
                microturbulence_vel=float(pp['vturb_kms']), macroturbulence=vmac,
                vsini=vsini, R=R, verbose=0, code=args.code,
                use_molecules=(args.code == 'turbospectrum'), tmp_dir=str(tmp))
            if not np.all(np.isfinite(sf)):
                raise RuntimeError('nonfinite synthesis')
            cache[key] = sf
        return cache[key]

    if args.smoke:
        of = np.interp(ow, sw, synth(8.1))
    rows = []
    for half, cont, dc, do, dx in (
        (.25, 0., 0., 0., 0.), (1.1, 0., 0., 0., 0.),
        (.25, -.001, 0., 0., 0.), (.25, .001, 0., 0., 0.),
        (.25, 0., -.1, 0., 0.), (.25, 0., .1, 0., 0.),
        (.25, 0., 0., -.1, 0.), (.25, 0., 0., .1, 0.),
        (.25, 0., 0., 0., -.1), (.25, 0., 0., 0., .1),
    ):
        mask = np.abs(ow*10-wave) <= half
        if mask.sum() < 5: raise RuntimeError('insufficient pixels')
        def objective(a):
            sf = synth(a, state['C']+dc, state['O']+do, params['vturb_kms']+dx)
            return float(np.sum(((of[mask]/(1+cont)-np.interp(ow[mask], sw, sf))/.01)**2))
        fit = minimize_scalar(objective, bounds=(6.3, 9.5), method='bounded',
                              options={'xatol': .005})
        if not fit.success or not np.isfinite(fit.fun) or min(fit.x-6.3, 9.5-fit.x) < .01:
            raise RuntimeError(f'Unconverged or boundary abundance fit: {fit}')
        row = dict(wavelength_air_A=wave, half_width_A=half, continuum_fraction=cont,
                   delta_C=dc, delta_O=do, delta_xi=dx, A_N=float(fit.x),
                   red_chi2=float(fit.fun/(mask.sum()-1)), n_pixels=int(mask.sum()),
                   status='diagnostic_only', geometry='disk_integrated_flux',
                   scale='1D-LTE', synthetic_observation=args.smoke)
        rows.append(row)
        (args.out/'fits.json').write_text(json.dumps(rows, indent=2)+'\n')
        print(json.dumps(row), flush=True)
        if args.smoke: break
        if args.quick and len(rows) == 2: break
    best = rows[0]['A_N']
    # Preserve the full best-fit profile; no atomic-only deblend is claimed.
    np.savetxt(args.out/'profiles.csv', np.column_stack((sw*10, synth(best))),
               delimiter=',', header='wavelength_air_A,model_flux', comments='')
    if overlay:
        selected = (np.abs(ll['wave_A']-wave)<.05) & (ll['element']=='N 1')
        without = c._synth_window(sw, atm, params, ll[~selected], iso, sab,
            c._fixed_ab({**state, 'N': best}, codes), broad, True, str(tmp))
        np.savetxt(args.out/'atomic_contribution.csv',
                   np.column_stack((sw*10, synth(best), without)), delimiter=',',
                   header='wavelength_air_A,full_flux,only_target_NI_removed_flux', comments='')
    np.savetxt(args.out/'observed.csv', np.column_stack((ow*10,of)), delimiter=',',
               header='wavelength_air_A,flux', comments='')
    sources = {}
    for p in [args.engine_root/'pipeline/cno_synthesis.py',
              args.engine_root/'pipeline/abundances_derive.py', llpath]:
        sources[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    provenance = dict(engine_root=str(args.engine_root), source_sha256=sources,
                      holding=str(win.holding), params=params, background=state,
                      broadening=broad, gf=gf, engine=args.code,
                      reference_overlay=overlay,
                      limitation='1D LTE flux; fixed broadening; no empirical CN subtraction; '
                                 'continuum perturbation steps are probes, not adopted errors')
    (args.out/'provenance.json').write_text(json.dumps(provenance,indent=2,default=str)+'\n')
    if args.smoke and abs(best-8.1) > .025:
        raise RuntimeError(f'synthetic recovery failed: {best}')


if __name__ == '__main__':
    main()

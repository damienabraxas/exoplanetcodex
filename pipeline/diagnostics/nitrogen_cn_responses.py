"""Sirius-only exact-window CN refits; measured probes, never a complete budget.

Uses the established CNO fitter and its original IAG CN window pool. Outputs
are exclusive-create and checkpoint after each fit. C/O steps are numerical
probes, distinct from adopted abundance errors. Stellar steps come from the
shared params_and_deltas source, including the solar xi allowance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import math
from pathlib import Path
from pipeline import molecular_identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--engine-root', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--group', choices=['coupling', 'stellar', 'continuum'], required=True)
    ap.add_argument('--prior', type=Path, required=True)
    ap.add_argument('--region', default='nir_cn_iag', help='declared CNO region to probe')
    ap.add_argument('--element', default='N', choices=['C', 'N', 'O'])
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(args.engine_root.resolve()))
    from pipeline import cno_synthesis as c
    from pipeline.uncertainty_stack import params_and_deltas
    import numpy as np

    region = c.REGIONS[args.region]
    nd = [d for d in c.REGION_DIAGNOSTICS[args.region] if d.element == args.element and d.role == 'primary']
    if len(nd) != 1:
        raise RuntimeError(f'expected one primary N diagnostic for {args.region}, got {[d.key for d in nd]}')
    diag = nd[0]
    broad = c.preflight(region, 'solar', (diag,))
    params, deltas = params_and_deltas('solar')
    llpath, label, gf = c.region_atomic_linelist(region, (diag,))
    if llpath is None:
        ll, iso, chem = c._load_synth_resources()
    else:
        ll, iso, chem = c._load_synth_resources(str(llpath), apply_canonical_gf=(gf or {}).get('apply_canonical_gf', False))
    sab = c.ispec.read_solar_abundances(c._ISPEC_SOLAR_ABUND_FILE)
    codes = c._atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)
    state = c._seed_abundances('solar', params, codes,
                             c._solar_A(('C', 'N', 'O', 'Ni'), chem, sab), params['feh'])
    # Existing product composition, not a new solar calibration.
    prior_path = args.prior
    prior = json.loads(prior_path.read_text())
    state.update({k: prior['abundances'][k] for k in ('C', 'N', 'O')})
    ow, of = c._load_region_spectrum('solar', region)
    mask = np.zeros(len(ow), dtype=bool)
    for lo, hi in diag.windows_A:
        m = (ow*10 >= lo) & (ow*10 <= hi)
        if m.sum() < 5 or not np.isfinite(of[m]).all():
            raise RuntimeError('Missing/nonfinite original CN window; changed pool refused')
        mask |= m
    ow, of = ow[mask], of[mask]
    # Fe-style local continuum diagnostic: estimate the residual placement scale
    # from top-percentile anchor pixels around the declared windows. This does not
    # renormalize the source product; it measures the observed local deviation from
    # unity and uses that measured scale for sensitivity probes.
    import numpy as np
    anchor_resid = []
    for lo, hi in diag.windows_A:
        center = 0.5 * (lo + hi); half = 0.5 * (hi - lo)
        anchor_mask = (ow * 10 > center - half - 1.2) & (ow * 10 < center + half + 1.2)
        edge = anchor_mask & ((ow * 10 < center - 0.7 * half) | (ow * 10 > center + 0.7 * half))
        if edge.sum() >= 8:
            vals = of[edge]
            top = vals[vals >= np.percentile(vals, 70)]
            if len(top) >= 4:
                anchor_resid.extend((top - 1.0).tolist())
    continuum_amp = float(1.4826 * np.median(np.abs(np.asarray(anchor_resid) - np.median(anchor_resid)))) if anchor_resid else 0.001
    continuum_amp = max(continuum_amp, 1e-5)
    # Molecular identity is DECLARED per diagnostic, never defaulted.  The previous
    # form here was `'NH' if diag.key == 'NH_AX' else '12C14N'`, which stamped OH_H
    # (1.53-1.69 um) and CO_K (2.31-2.46 um) as CN A-X (0-0) -- a band with no
    # transitions in either window.  pool is hashed into indicator_id, so that wrong
    # label propagated into the identity the uncertainty contract joins on.
    pool = molecular_identity.pool_for(diag.key, diag.windows_A)
    identity = f'{diag.key}:molecular:{args.element}:windows_sha256:' + hashlib.sha256(
        json.dumps(pool, sort_keys=True).encode()).hexdigest()
    sources = [Path(c.__file__), args.engine_root/'pipeline/uncertainty_stack.py',
               prior_path, Path(__file__)]
    if llpath is not None:
        sources.insert(2, Path(llpath))
    provenance = {'engine_root': str(args.engine_root), 'group': args.group,
                  'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                  'params': params, 'adopted_parameter_deltas': deltas, 'background': state,
                  'broadening': broad, 'holding': c.holding_for_region(region),
                  'pool': pool, 'indicator_id': identity,
                  'n_pixels': len(ow), 'geometry': 'disk_integrated_flux',
                  'scale': '1D-LTE', 'statistical_correlation_model': 'UNRESOLVED',
                  'C_O_covariance': None, 'continuum_anchor_amp': continuum_amp,
                  'continuum_anchor_n': len(anchor_resid), 'status': 'diagnostic_only'}
    (args.out/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    jobs = [('nominal', {}, {})]
    if args.group == 'coupling':
        for el in ('C', 'O'):
            for sign in (-1, 1):
                jobs.append((f'{el}_{sign:+d}', {}, {el: state[el]+sign*.1}))
    elif args.group == 'stellar':
        for p, step in deltas.items():
            if step > 0:
                for sign in (-1, 1):
                    jobs.append((f'{p}_{sign:+d}', {p: params[p]+sign*step}, {}))
    else:
        jobs.extend([('continuum_level_-1', {}, {'_continuum_fraction': -continuum_amp}),
                     ('continuum_level_+1', {}, {'_continuum_fraction': continuum_amp}),
                     ('continuum_slope_-1', {}, {'_continuum_slope': -continuum_amp}),
                     ('continuum_slope_+1', {}, {'_continuum_slope': continuum_amp})])
    rows = []
    def json_safe(value):
        if isinstance(value, float) and not math.isfinite(value):
            return None
        if isinstance(value, dict):
            return {k: json_safe(v) for k, v in value.items()}
        if isinstance(value, list):
            return [json_safe(v) for v in value]
        return value
    for name, overrides, composition in jobs:
        pp, st = {**params, **overrides}, {**state, **composition}
        continuum_fraction = float(st.pop('_continuum_fraction', 0.0))
        continuum_slope = float(st.pop('_continuum_slope', 0.0))
        if continuum_slope:
            x = (ow * 10 - float(np.mean([a for w in diag.windows_A for a in w]))) / 10.0
            scale = 1.0 + continuum_slope * x / max(np.ptp(ow * 10), 1.0)
            fit_flux = of / scale
        else:
            fit_flux = of / (1.0 + continuum_fraction)
        atm = c._load_atmosphere(pp['teff_K'], pp['logg'], pp['feh'], pp['vturb_kms'])
        tmp = args.out/name; tmp.mkdir()
        started = time.monotonic()
        fit = c._fit_element(ow, fit_flux, atm, pp, args.element, st, codes, diag.windows_A,
                             True, broad, state[args.element]-1.2, state[args.element]+1.2,
                             ll, iso, sab, str(tmp))
        row = {'name': name, 'params': pp, 'composition': st,
               'indicator_id': identity, 'fit': fit, 'wall_s': time.monotonic()-started}
        rows.append(row)
        (args.out/'fits.json').write_text(json.dumps(json_safe(rows), indent=2, allow_nan=False)+'\n')
        print(json.dumps(row), flush=True)
        if fit['status'] != 'ok' or not fit.get('constrained') or fit['n_pix'] != len(ow):
            raise RuntimeError(f'Invalid or changed-pool CN fit: {name}')


if __name__ == '__main__':
    main()

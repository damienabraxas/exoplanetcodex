"""Finite O777 response geometry from immutable existing profiles; no inversion."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
AXES = ('temperature', 'gravity', 'xi')


def direction_geometry(columns):
    """Unit-column Euclidean geometry; deliberately has no noise/parameter metric."""
    columns = np.asarray(columns, dtype=float)
    if columns.ndim != 2 or columns.shape[1] != 3 or not np.isfinite(columns).all():
        raise ValueError('Require three finite response columns')
    norm = np.linalg.norm(columns, axis=0)
    if np.any(norm == 0): raise ValueError('Zero response direction')
    unit = columns/norm
    singular = np.linalg.svd(unit, compute_uv=False)
    return unit, dict(column_norms=norm.tolist(), signed_cosines=(unit.T@unit).tolist(),
                      singular_values=singular.tolist(), smallest_to_largest_singular_ratio=float(singular[-1]/singular[0]))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True); args = p.parse_args()
    from pipeline.bride_reference import fingerprint
    import matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    roots = [ROOT/'data/results/rya1281'/s for s in ('thermal_slice_v1', 'xi_slice_v1')]
    receipts = [json.loads((d/'analysis/report.json').read_text()) for d in roots]
    t, x = receipts
    for key in ('observation', 'policy', 'broadening', 'ispec_effects'):
        if t[key] != x[key]: raise ValueError('Mismatched shared observational conditioning: '+key)
    for i, name in enumerate(('solar5772_xi090', 'solar5772_xi110')):
        if fingerprint(roots[1]/'runs'/name/'report.json') != t['response_directions']['prior_xi_reports'][i]:
            raise ValueError('Changed source-matched xi response input')
    arrays = []
    for d, r in zip(roots, receipts):
        if fingerprint(d/'analysis/profiles.csv') != r['profiles']: raise ValueError('Changed archived profiles')
        arrays.append(np.genfromtxt(d/'analysis/profiles.csv', delimiter=',', names=True))
    a, b = arrays; wave = a['wavelength_air_A']
    np.testing.assert_array_equal(wave, b['wavelength_air_A'])
    np.testing.assert_array_equal(a['shipped_observation'], b['shipped_observation'])
    mask = np.zeros(len(wave), dtype=bool)
    centres = (7771.944, 7774.166, 7775.388)
    for centre in centres: mask |= (wave >= centre-.5) & (wave <= centre+.5)
    if not np.any(mask) or not np.allclose(np.diff(wave), .005, rtol=0, atol=1e-11):
        raise ValueError('Changed held profile grid/apertures')
    area_columns = np.array([t['response_directions']['area_change_vectors_mA'][k] for k in AXES]).T
    _, area_geometry = direction_geometry(area_columns)
    report = dict(schema='codex.bride_profile_geometry/1', scientific_admission='HOLD', generator=fingerprint(__file__),
        source_reports=[fingerprint(d/'analysis/report.json') for d in roots],
        axis_order=AXES, intervals=t['response_directions']['intervals'], area_geometry=area_geometry,
        pixel_count=int(mask.sum()), wavelength_step_air_A=.005, aperture_centres_air_A=centres, aperture_half_width_A=.5, variants=[],
        interpretation='Unit directions discard response amplitude and parameter units. Fixed-grid unweighted finite-response geometry, not a noise-weighted Jacobian, Fisher matrix, parameter covariance, inverse fit, rank threshold or measured identifiability. Broadening/continuum/atom/opacity uncertainties remain open; no favorable variant selected.')
    fig, panels = plt.subplots(2, 3, figsize=(13, 7), sharex='col')
    profiles = [wave]; labels = ['wavelength_air_A']
    for row, variant in enumerate(('canonical_broadened', 'plus_gravitational_shift')):
        np.testing.assert_array_equal(a['thermal_nominal_'+variant], b['solar5772_xi100_'+variant])
        full = np.column_stack((a['thermal_T5822_'+variant]-a['thermal_T5672_'+variant],
                                a['thermal_g4538_'+variant]-a['thermal_g4338_'+variant],
                                b['solar5772_xi110_'+variant]-b['solar5772_xi090_'+variant]))
        unit, geometry = direction_geometry(full[mask])
        report['variants'].append(dict(name=variant, geometry=geometry, nominal_array_identity=True))
        scaled = full/np.array(geometry['column_norms'])
        for column, name in enumerate(AXES):
            profiles.append(scaled[:, column]); labels.append(variant+'_'+name+'_unit_response')
            for panel in panels[row]: panel.plot(wave, scaled[:, column], label=name)
        for panel, centre in zip(panels[row], centres):
            panel.set_xlim(centre-.5, centre+.5); panel.grid(alpha=.2); panel.legend(fontsize=8)
            panel.set_xlabel('AIR wavelength (Angstrom)')
        panels[row, 0].set_ylabel(variant+'\nUnit response (Euclidean norm)')
    fig.suptitle('Finite temperature/gravity/xi response shapes in held O777 profiles\nBoth fixed broadening/RV variants; geometry only, no noise model or fit')
    fig.tight_layout(); args.output.mkdir(parents=True, exist_ok=False)
    fig.savefig(args.output/'profile_geometry.png', dpi=150); plt.close(fig)
    np.savetxt(args.output/'directions.csv', np.column_stack(profiles), delimiter=',', header=','.join(labels), comments='')
    report.update(figure=fingerprint(args.output/'profile_geometry.png'), directions=fingerprint(args.output/'directions.csv'))
    (args.output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(area=area_geometry, variants=report['variants'])))


if __name__ == '__main__': main()

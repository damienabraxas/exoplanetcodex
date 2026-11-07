"""Decompose fixed O777 residuals and measure existing local-continuum controls.

No fit, abundance change or baseline replacement. Model-guided controls require
complete +/-5 A model coverage; truncated coverage is explicitly withheld.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def integral(wave, flux, lo, hi):
    if lo < wave[0] or hi > wave[-1]:
        raise ValueError('Integral outside coverage')
    x = np.r_[lo, wave[(wave > lo) & (wave < hi)], hi]
    y = np.interp(x, wave, flux)
    return float(np.sum((y[1:] + y[:-1]) * np.diff(x) / 2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison', type=Path, required=True)
    parser.add_argument('--continuum-reference', choices=('synthesis','iag'), default='synthesis')
    parser.add_argument('--wide-spectrum', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    import os
    os.environ['CODEX_CONT_REF'] = args.continuum_reference
    from pipeline import local_continuum as lc
    from pipeline.telluric_policy import gate_holding
    from pipeline.prenormalised_guard import assert_not_renormalising
    if lc.ENV_PERCENTILE != 95 or lc.MODEL_Q != 90 or lc.CONT_REF != args.continuum_reference:
        raise ValueError('Nominal continuum configuration required')
    holding = 'solar_kpno_kurucz2005_corrected'
    allowed, reason = gate_holding(holding, 'kpno_solar_atlas')
    if not allowed:
        raise ValueError(reason)
    assert_not_renormalising(holding, pre_normalised=True, local_window_envelope=True)
    comparison = json.loads(args.comparison.read_text())
    observation = require_codex_path('data.kurucz2005_residual')
    if fingerprint(observation) != comparison['observation']:
        raise ValueError('Changed observation')
    wave, observed = read_kurucz2005_residual(observation, 7765, 7782)
    case = next(c for c in comparison['cases'] if c['name'] == 'solar5772_bound_Fe_sampled')
    profiles = args.comparison.parent / (case['name'] + '_profiles.csv')
    if fingerprint(profiles) != case['profiles']:
        raise ValueError('Changed profiles')
    data = np.genfromtxt(profiles, delimiter=',', names=True)
    wide_receipt = None
    if args.wide_spectrum is not None:
        wide_receipt = json.loads((args.wide_spectrum/'report.json').read_text())
        if (wide_receipt['source_manifest'] != case['source_manifest']
                or wide_receipt['source_flux'] != case['flux']
                or fingerprint(args.wide_spectrum/'flux.csv') != wide_receipt['flux']):
            raise ValueError('Wide spectrum source differs')
        sys.path.insert(0, str(require_codex_path('engines.ispec')))
        import ispec
        from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
        raw = np.genfromtxt(args.wide_spectrum/'flux.csv', delimiter=',', names=True)
        data = np.zeros(len(raw), dtype=[(k,float) for k in ('wavelength_air_A','canonical_broadened','plus_gravitational_shift')])
        data['wavelength_air_A'] = raw['wavelength_air_A']
        segments = np.array([(raw['wavelength_air_A'][0]/10, raw['wavelength_air_A'][-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
        b = comparison['broadening']
        for name, rv in [('canonical_broadened',0.),('plus_gravitational_shift',SOLAR_GRAVITATIONAL_REDSHIFT)]:
            data[name] = ispec.apply_post_fundamental_effects(raw['wavelength_air_A']/10, raw['normalized_flux'].copy(),segments,
                macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(rv,))
    mw = data['wavelength_air_A']
    result = dict(schema='codex.bride_kitt_residual_audit/1', scientific_admission='HOLD',
        generator=fingerprint(__file__), comparison=fingerprint(args.comparison),
        observation=fingerprint(observation), profiles=fingerprint(profiles),
        wide_spectrum=wide_receipt, continuum_code=fingerprint(lc.__file__), holding=holding, policy=reason,
        caveats='Fixed diagnostic partitions are not isolated oxygen EWs or blend identifications. No errors/covariance or optimization. Local controls never replace immutable baseline. Model-guided anchors remain incomplete-opacity diagnostics.', lines=[])
    result['continuum_reference'] = args.continuum_reference
    if args.continuum_reference == 'iag':
        allowed_iag, why_iag = gate_holding('solar_iag','iag_fts_solar_atlas')
        if not allowed_iag: raise ValueError(why_iag)
        result['reference_reader'] = fingerprint(ROOT/'scripts/measure_band_ew.py')
        result['reference_atlas'] = dict(source=fingerprint(require_codex_path('data.solar_iag_baker2020')),policy=why_iag,
            method='Existing solar_atlas_reference; AIR reader and Gaussian catalog-R broadening')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3,1,figsize=(11,9))
    for line_index, centre in enumerate((7771.944, 7774.166, 7775.388)):
        record, envelope = lc.fit(wave, observed, centre)
        variants = [('shipped', observed)]
        if envelope is not None:
            variants.append(('observed_p95_envelope', observed / envelope))
        line = dict(centre_air_A=centre, envelope=asdict(record), variants=[])
        reference_wave, reference_flux = mw, data['plus_gravitational_shift']
        if args.continuum_reference == 'iag':
            ref = lc.solar_atlas_reference(centre-5,centre+5,comparison['broadening']['R'])
            if ref is None: raise ValueError('IAG reference unavailable; no silent synthesis substitution')
            reference_wave, reference_flux = ref
        coverage = reference_wave[0] <= centre - 5 and reference_wave[-1] >= centre + 5
        line['model_guided_complete_coverage'] = bool(coverage)
        if coverage:
            # Avoid iSpec's zeroed end samples; all nominal anchors are interior here.
            keep = (reference_wave >= centre - 5.01) & (reference_wave <= centre + 5.01)
            rec, cont = lc.fit_model_guided(wave, observed, centre, reference_wave[keep],
                reference_flux[keep], exclude_half_width_A=.5, reference=args.continuum_reference)
            line['model_guided'] = asdict(rec)
            mod = np.interp(wave, reference_wave[keep], reference_flux[keep])
            eligible = (abs(wave-centre) <= 5) & (abs(wave-centre) > .5)
            threshold = np.percentile(mod[eligible], lc.MODEL_Q)
            anchors = wave[eligible & (mod >= threshold)]
            line['model_guided_anchor_geometry'] = dict(count=int(len(anchors)),
                min_air_A=float(anchors.min()), max_air_A=float(anchors.max()),
                blue_count=int(sum(anchors < centre)), red_count=int(sum(anchors > centre)),
                note='Pre-clipping selection, not final robust-fit anchors')
            if cont is not None:
                variants.append(('model_guided_'+args.continuum_reference+'_diagnostic', observed / cont))
        else:
            line['model_guided'] = dict(applied=False, reason='WITHHELD: saved model lacks complete +/-5 A neighborhood; no truncated-envelope substitution')
        for label, obs in variants:
            row = dict(name=label, observed_area_mA=1000*integral(wave, 1-obs, centre-.5, centre+.5), residuals=[])
            for model_name in ('canonical_broadened', 'plus_gravitational_shift'):
                model = data[model_name]
                zones = []
                for lo, hi, name in ((-.5,-.15,'blue_wing'),(-.15,.15,'core'),(.15,.5,'red_wing')):
                    observed_area = 1000*integral(wave, 1-obs, centre+lo, centre+hi)
                    model_area = 1000*integral(mw, 1-model, centre+lo, centre+hi)
                    zones.append(dict(zone=name, observed_area_mA=observed_area, model_area_mA=model_area,
                        missing_area_mA=observed_area-model_area))
                row['residuals'].append(dict(model=model_name,zones=zones))
            line['variants'].append(row)
        view = abs(wave-centre) <= .6
        for label, obs in variants:
            axes[line_index].plot(wave[view]-centre,obs[view],label=label,lw=1)
        model_view = abs(mw-centre) <= .6
        axes[line_index].plot(mw[model_view]-centre,data['plus_gravitational_shift'][model_view],label='Fixed O8.70; canonical gravitational shift',color='black',lw=1)
        axes[line_index].set_title(f'{centre:.3f} AIR A; reference {args.continuum_reference}')
        axes[line_index].legend(fontsize=7);axes[line_index].grid(alpha=.2)
        axes[line_index].set_ylabel('Normalized flux')
        result['lines'].append(line)
    # Candidate inventory is wavelength selection only: no strengths/completeness inferred.
    import csv
    canonical = ROOT/'data/linelists/canonical_gf.csv'
    with canonical.open() as stream:
        result['canonical_atomic_candidates'] = [r for r in csv.DictReader(stream) if 7763 <= float(r['wavelength_air_A']) <= 7783]
    result['canonical_atomic_source'] = fingerprint(canonical)
    result['candidate_scope'] = 'Wavelength-only inventory, not installed blend opacity or molecular completeness'
    args.output.mkdir(parents=True, exist_ok=False)
    axes[-1].set_xlabel('Offset from line centre (AIR Angstrom)')
    fig.suptitle('Separate local-continuum controls; no fit or scientific adoption')
    fig.tight_layout();fig.savefig(args.output/'continuum_controls.png',dpi=140);plt.close(fig)
    result['figure'] = fingerprint(args.output/'continuum_controls.png')
    (args.output/'report.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    for line in result['lines']:
        print(line['centre_air_A'], 'envelope', line['envelope']['level_at_centre'],
              [(v['name'],round(v['observed_area_mA'],5)) for v in line['variants']],
              'full_model_coverage', line['model_guided_complete_coverage'])


if __name__ == '__main__':
    main()

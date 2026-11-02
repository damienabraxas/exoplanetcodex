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
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline import local_continuum as lc
    from pipeline.telluric_policy import gate_holding
    from pipeline.prenormalised_guard import assert_not_renormalising
    if lc.ENV_PERCENTILE != 95 or lc.MODEL_Q != 90 or lc.CONT_REF != 'synthesis':
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
    mw = data['wavelength_air_A']
    result = dict(schema='codex.bride_kitt_residual_audit/1', scientific_admission='HOLD',
        generator=fingerprint(__file__), comparison=fingerprint(args.comparison),
        observation=fingerprint(observation), profiles=fingerprint(profiles),
        continuum_code=fingerprint(lc.__file__), holding=holding, policy=reason,
        caveats='Fixed diagnostic partitions are not isolated oxygen EWs or blend identifications. No errors/covariance or optimization. Local controls never replace immutable baseline. Oxygen-only anchors are incomplete-opacity diagnostics.', lines=[])
    for centre in (7771.944, 7774.166, 7775.388):
        record, envelope = lc.fit(wave, observed, centre)
        variants = [('shipped', observed)]
        if envelope is not None:
            variants.append(('observed_p95_envelope', observed / envelope))
        line = dict(centre_air_A=centre, envelope=asdict(record), variants=[])
        coverage = mw[0] <= centre - 5 and mw[-1] >= centre + 5
        line['model_guided_complete_coverage'] = bool(coverage)
        if coverage:
            # Avoid iSpec's zeroed end samples; all nominal anchors are interior here.
            keep = (mw >= centre - 5) & (mw <= centre + 5)
            rec, cont = lc.fit_model_guided(wave, observed, centre, mw[keep],
                data['plus_gravitational_shift'][keep], exclude_half_width_A=.5)
            line['model_guided'] = asdict(rec)
            mod = np.interp(wave, mw[keep], data['plus_gravitational_shift'][keep])
            eligible = (abs(wave-centre) <= 5) & (abs(wave-centre) > .5)
            threshold = np.percentile(mod[eligible], lc.MODEL_Q)
            anchors = wave[eligible & (mod >= threshold)]
            line['model_guided_anchor_geometry'] = dict(count=int(len(anchors)),
                min_air_A=float(anchors.min()), max_air_A=float(anchors.max()),
                blue_count=int(sum(anchors < centre)), red_count=int(sum(anchors > centre)),
                note='Pre-clipping selection, not final robust-fit anchors')
            if cont is not None:
                variants.append(('model_guided_incomplete_opacity_diagnostic', observed / cont))
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
        result['lines'].append(line)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/'report.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    for line in result['lines']:
        print(line['centre_air_A'], 'envelope', line['envelope']['level_at_centre'],
              [(v['name'],round(v['observed_area_mA'],5)) for v in line['variants']],
              'full_model_coverage', line['model_guided_complete_coverage'])


if __name__ == '__main__':
    main()

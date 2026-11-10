"""Read-only Baker IAG quality-flag coverage for the frozen O777 controls."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    from astropy.io import fits
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.uv_conditioning import vac_to_air
    atlas=require_codex_path('data.solar_iag_baker2020')
    library=require_codex_path('refs.library')
    result=dict(schema='codex.bride_iag_quality/1',scientific_admission='HOLD',
        generator=fingerprint(Path(__file__)),atlas=fingerprint(atlas),
        helper=fingerprint(ROOT/'pipeline/local_continuum.py'),
        reader=fingerprint(ROOT/'scripts/measure_band_ew.py'),
        sources={name:fingerprint(library/name) for name in
                 ('Baker_2020_ApJS_247_24.pdf','aa27530-15.pdf')},
        source_interpretation='Baker2020 p7: flag0 robust; flags1/2/3 mark >10percent/>25percent/saturated telluric regions replaced by spline. No mask or flux changes.',
        frame='Native vacuum wavenumber cm^-1 -> vacuum A -> AIR once, same reader conversion.',
        caveats='Flag0 is not proof of zero telluric residual or a covariance model. No second telluric correction. Catalog Gaussian-R smoothing is not a calibrated FTS line-spread function.',
        windows=[])
    with fits.open(atlas) as hdus:
        data=hdus[1].data
        wave=vac_to_air(1e8/np.asarray(data['v'],float))
        flags=np.asarray(data['flags'],float)
        if not np.isfinite(flags).all():raise ValueError('Nonfinite IAG quality flags')
        # Delivered flags are floating point and include fractional values outside
        # this task's windows. Preserve that discrepancy; never round categories.
        result['global_flag_provenance']=dict(minimum=float(flags.min()),maximum=float(flags.max()),
            distinct_values=len(np.unique(flags)),
            noncategorical_count=int(np.count_nonzero(~np.isin(flags,[0,1,2,3]))),
            disposition='Fractional delivered flags outside O777 require source-format clarification; no rounding or global categorical interpretation.')
        for centre in (7771.944,7774.166,7775.388):
            for half,label in ((.5,'O_aperture'),(5.,'continuum_window'),(6.,'full_helper_input')):
                mask=abs(wave-centre)<=half
                if not mask.any():raise ValueError('Missing IAG coverage')
                values,counts=np.unique(flags[mask],return_counts=True)
                result['windows'].append(dict(centre_air_A=centre,half_width_A=half,scope=label,
                    count=int(mask.sum()),flag_value_counts={repr(float(v)):int(n) for v,n in zip(values,counts)},
                    flagged_count=int(np.count_nonzero(flags[mask]))))
    # Fail rather than overwrite any prior diagnostic.
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()

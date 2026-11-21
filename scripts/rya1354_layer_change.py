"""Native n3-to-frozen-n4 E1 energy/strength decomposition; no observed input."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from rya1354_native_parts import load_parts
from rya1354_lifetime_summary import decay
from rya1354_frozen_summary import summarize

ROOT=Path(__file__).resolve().parents[1]


def decompose(previous,current):
    def components(model):
        return {r['lower_J']:r for r in model['channels'] if r['destination']=='to_ground'
                and r['operator']=='E1' and r['lower_ordinal']==1}
    old,new=components(previous),components(current)
    if set(old)!={0,1,2} or set(new)!=set(old):
        raise ValueError('Physical ground multiplet components differ')
    result=[]
    for j in sorted(old):
        a,b=old[j],new[j]
        energy=(b['native_energy_cm1']/a['native_energy_cm1'])**3
        gauges={}
        for gauge in ('B','C'):
            strength=b['native_S_'+gauge+'_au']/a['native_S_'+gauge+'_au']
            rate=b['A_'+gauge+'_s1']/a['A_'+gauge+'_s1']
            if not math.isclose(rate,energy*strength,rel_tol=5e-5):
                raise ValueError('Native E1 energy/strength identity differs')
            gauges[gauge]=dict(A_ratio=rate,native_strength_ratio=strength)
        result.append(dict(lower_J=j,native_energy_cubed_ratio=energy,gauges=gauges))
    return result


def report(base):
    files,_=load_parts(base/'lifetime_relax_v4')
    previous=decay(files,'runs/lifetime_relax_n3_v3/')
    current=summarize(base)['decay']
    return dict(schema='rya1354.layer-change-diagnostic.v13',components=decompose(previous,current),
        previous_model='Core-adapted n3 full2s-SD/DC',current_model='Frozen-prior-layer n4 full2s-SD/DC',
        lifetimes_ns={g:dict(n3=previous['model'][g]['model_one_photon_lifetime_ns'],
                            frozen_n4=current['model'][g]['model_one_photon_lifetime_ns']) for g in ('B','C')},
        scope='THEORY_TO_THEORY_E1_DECOMPOSITION_ONLY_NO_OBSERVED_ENERGY_SUBSTITUTION',
        scientific_admission='HOLD',adoption_allowed=False,uncertainty=None,covariance=None,calibrated_sigma=None,
        caveat='Changes include enlarged CSF space and optimized outer-layer functions; not a single-cause attribution or calibrated convergence error.',
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.write_text(json.dumps(report(args.base),indent=2)+'\n')

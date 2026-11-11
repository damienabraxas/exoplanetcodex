"""Independent replay of the declared frozen-prior-layer rate sensitivity."""
import argparse
import hashlib
import json
from pathlib import Path

from rya1354_native_parts import load_parts
from rya1354_frozen_layer import validate, PARENT, STATES, STATUS
from rya1354_lifetime_summary import decay

ROOT=Path(__file__).resolve().parents[1]
RUN='runs/lifetime_frozen_prior_n4_v1/'


def summarize(base):
    source,_=load_parts(base/'from_n3_terminal_v12')
    target,_=load_parts(base/'lifetime_aux_v4b')
    audit=validate(source,target)
    files,manifest=load_parts(base/'frozen_complete_v12')
    receipt=json.loads(files[RUN+'receipt.json'])
    if receipt['status']!=STATUS or receipt['source_audit']!=audit:
        raise ValueError('Frozen model status/source audit differs')
    for state in STATES:
        for name in ('state.c','state.w','state.m','state.lsj.lbl','isodata'):
            if files[RUN+state+'/'+name]!=source[PARENT+'outer_n4_'+state+'_verify/'+name]:
                raise ValueError('Frozen model selected state differs')
    for name,lower,operators,same in [('to_ground','ground','E1,M2,E3',False),('to_5S','quintet','M1,E2,M3',True)]:
        path=RUN+name+'/'
        for tag,state in [('lower',lower),('upper','triplet')]:
            for ext in ('c','w','m'):
                if files[path+tag+'.'+ext]!=files[RUN+state+'/state.'+ext]:
                    raise ValueError('Native transition input differs')
            if files[path+tag+'.lsj.lbl']!=files[RUN+state+'/state.lsj.lbl']:
                raise ValueError('Native transition label differs')
        if files[path+'isodata']!=files[RUN+lower+'/isodata']:
            raise ValueError('Native transition nucleus differs')
        controls={'rbiotransform':'y\nn\nlower\nupper\ny\n',
                  'rtransition':'y\nn\nlower\nupper\n'+operators+'\n'+('n\nn\n' if same else '')}
        for exe,expected in controls.items():
            if files[path+exe+'.stdin'].decode()!=expected:
                raise ValueError('Native rate control differs')
    result=decay(files,RUN,ci=False)
    return dict(schema='rya1354.frozen-prior-layer-sensitivity.v12',source_audit=audit,
        model=receipt['model'],decay=result,parts=manifest['parts'],
        scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        fully_released_qualified=False,accepted_lifetime=None,uncertainty=None,covariance=None,
        calibrated_sigma=None,generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interpretation='Model sensitivity only. Full-release failure retained. Gauge agreement and absolute accuracy are not calibrated. No grade or reference adoption.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    args.out.write_text(json.dumps(summarize(args.base),indent=2)+'\n')

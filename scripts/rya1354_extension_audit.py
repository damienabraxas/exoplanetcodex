"""Independent n5 frozen-extension replay; incomplete native states are refused."""
import argparse
import hashlib
import json
from pathlib import Path

from rya1354_native_parts import load_parts
from rya1354_frozen_extend import controls
from rya1354_frozen_layer import STATUS,STATES
from rya1354_frozen_summary import summarize as source_summary
from rya1354_lifetime_summary import scf_state,decay
from rya1354_triplet_correlation import require_triplet_quintet_order
from rya1354_radial_controls import initialization_input
from rya1354_wave_probe import decode
from rya1354_mr_inventory import inventory

ROOT=Path(__file__).resolve().parents[1]
RUN='runs/lifetime_frozen_prior_n5_v1/'
CORE='runs/lifetime_core_n5_v1/'
SOURCE='runs/lifetime_frozen_prior_n4_v1/'


def stage(files,target,source,state,verify):
    case='outer_n5_'+state+('_verify' if verify else '')
    path=RUN+case+'/'
    # Read residuals before considering any promoted-state aliases or labels.
    from rya1354_two_method_summary import scf_diagnostics
    d=scf_diagnostics(files[path+'rmcdhf.stdout'].decode(),files[path+'rmcdhf.stderr'].decode())
    if d['iteration_limit_hit'] or not d['all_reported_functions_below_ACC']:
        raise ValueError('Unqualified native n5 orbital residuals')
    diagnostic=scf_state(files,path.rstrip('/'),state)
    if files[path+'rmcdhf.stdin'].decode()!=controls(state,verify):
        raise ValueError('Native n5 optimization controls differ')
    for name in ('isodata','rcsf.inp','rwfn.inp'):
        if files[path+name]!=files[path+'fixed_load/'+name]:
            raise ValueError('Native n5 loading inputs differ')
    if files[path+'fixed_load/rmcdhf.stdin'].decode()!=controls(state,True,True):
        raise ValueError('Native n5 loading controls differ')
    if 'All subshell radial wavefunctions are fixed' not in (files[path+'fixed_load/rmcdhf.stdout']+files[path+'fixed_load/rmcdhf.stderr']).decode():
        raise ValueError('Native n5 all-fixed loading marker absent')
    for selected,native in (('state.c','rcsf.inp'),('state.w','rwfn.out'),('state.m','rmix.out')):
        if files[path+selected]!=files[path+native]:
            raise ValueError('Native n5 selected state differs')
    for name in ('state.c','isodata'):
        if files[path+name]!=target[CORE+'CV_2s_DC_'+state+'/'+name]:
            raise ValueError('Native n5 target CSF/nucleus differs')
    inventory(files[path+'state.c'].decode(),'even' if state=='ground' else 'odd',True)
    a,b=decode(files[path+'fixed_load/rwfn.out']),decode(files[path+'state.w'])
    fixed=[k for k in a if k[0]<=4]
    expected=decode(target[CORE+'CV_2s_DC_'+state+'/state.w'])
    if set(a)!=set(b) or set(a)!=set(expected) or len(a)!=25 or len(fixed)!=16 or any(a[k]!=b[k] for k in fixed):
        raise ValueError('Native n5 frozen radial records differ')
    seed=RUN+('outer_n5_'+state if verify else 'initial_'+state)+'/'
    if files[path+'rwfn.inp']!=files[seed+'state.w']:
        raise ValueError('Native n5 initializer/fresh identity differs')
    initial=RUN+'initial_'+state+'/'
    for name,selected in (('state.c','state.c'),('rcsf.inp','state.c'),('isodata','isodata')):
        if files[initial+name]!=target[CORE+'CV_2s_DC_'+state+'/'+selected]:
            raise ValueError('Native n5 initialization CSF/nucleus differs')
    if files[initial+'previous.w']!=source[SOURCE+state+'/state.w']:
        raise ValueError('Native n5 prior layer source differs')
    if files[initial+'rwfnestimate.stdin'].decode()!=initialization_input(3,439,.03125):
        raise ValueError('Native n5 initialization controls differ')
    if files[initial+'state.w']!=files[initial+'rwfn.inp']:
        raise ValueError('Native n5 initialized wave differs')
    return dict(diagnostics=diagnostic,exact_frozen_n_le_4_records=16)


def summarize(base,bundle):
    prior=source_summary(base)
    source,_=load_parts(base/'frozen_complete_v12')
    target,_=load_parts(base/'lifetime_core_n5_v4')
    files,manifest=load_parts(bundle)
    receipt=json.loads(files[RUN+'receipt.json'])
    if receipt['status']!=STATUS or receipt['n']!=5:
        raise ValueError('N5 extension has no qualified terminal model')
    if (receipt['source']['source_audit']!=prior['source_audit']
            or receipt['source']['receipt_sha256']!=hashlib.sha256(source[SOURCE+'receipt.json']).hexdigest()
            or receipt['target_core_receipt_sha256']!=hashlib.sha256(target[CORE+'receipt.json']).hexdigest()):
        raise ValueError('N5 source/target receipt differs')
    states={}
    for state in STATES:
        for verify in (False,True):
            name='outer_n5_'+state+('_verify' if verify else '')
            states[name]=stage(files,target,source,state,verify)
    require_triplet_quintet_order(files[RUN+'outer_n5_triplet_verify/state.lsj.lbl'].decode(),
                                 files[RUN+'outer_n5_quintet_verify/state.lsj.lbl'].decode())
    for name,lower,operators,same in [('to_ground','ground','E1,M2,E3',False),('to_5S','quintet','M1,E2,M3',True)]:
        path=RUN+name+'/'
        for tag,state in [('lower',lower),('upper','triplet')]:
            parent=RUN+'outer_n5_'+state+'_verify/'
            for ext in ('c','w','m'):
                if files[path+tag+'.'+ext]!=files[parent+'state.'+ext]:
                    raise ValueError('N5 transition input differs')
            if files[path+tag+'.lsj.lbl']!=files[parent+'state.lsj.lbl']:
                raise ValueError('N5 transition physical labels differ')
        if files[path+'isodata']!=files[RUN+'outer_n5_'+lower+'_verify/isodata']:
            raise ValueError('N5 transition nucleus differs')
        if (files[path+'rbiotransform.stdin'].decode()!='y\nn\nlower\nupper\ny\n'
                or files[path+'rtransition.stdin'].decode()!='y\nn\nlower\nupper\n'+operators+'\n'+('n\nn\n' if same else '')):
            raise ValueError('N5 transition controls differ')
    result=decay(files,RUN)
    return dict(schema='rya1354.independent-n5-extension-audit.v14',states=states,decay=result,
        previous_n4=prior['decay']['model'],parts=manifest['parts'],scientific_admission='HOLD',
        exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,
        uncertainty=None,covariance=None,calibrated_sigma=None,fully_released_qualified=False,
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        caveat='Layerwise frozen-prior-layer sensitivity; no full relaxation, calibrated accuracy or gauge selection. Historical n3 component limitation retained.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--bundle',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.write_text(json.dumps(summarize(args.base,args.bundle),indent=2)+'\n')

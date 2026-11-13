"""Independent matched qualified-n5 DC/Breit replay; not accuracy calibration."""
import argparse
import hashlib
import json
from pathlib import Path

from rya1354_native_parts import load_parts
from rya1354_extension_audit import summarize as source_summary
from rya1354_qualified_ci import controls,PHYSICS
from rya1354_lifetime_ladder import STATES,identity
from rya1354_lifetime_summary import decay
from rya1354_triplet_correlation import require_triplet_quintet_order
from rya1354_layer_change import decompose

ROOT=Path(__file__).resolve().parents[1]
RUN='runs/lifetime_qualified_ci_n5_v1/'
SOURCE='runs/lifetime_frozen_prior_n5_v1/'


def summarize(base,bundle):
    files,manifest=load_parts(bundle)
    if RUN+'receipt.json' not in files:
        raise ValueError('No terminal qualified CI receipt; startup evidence is not a result')
    prior=source_summary(base,base/'n5_complete_v15');prior['parts']={}
    source,_=load_parts(base/'n5_complete_v15')
    receipt=json.loads(files[RUN+'receipt.json'])
    if (receipt['status']!='BALANCED_CORE_LIFETIMES_COMPLETED_REVIEW_REQUIRED'
            or receipt['n']!=5 or receipt['cases']!=list(PHYSICS)):
        raise ValueError('Qualified CI terminal scope differs')
    if receipt['source_audit']!=prior or receipt['source_receipt_sha256']!=hashlib.sha256(source[SOURCE+'receipt.json']).hexdigest():
        raise ValueError('Qualified CI source provenance differs')
    cases={};warnings={}
    for case in PHYSICS:
        root=RUN+case+'/'
        for state in STATES:
            path=root+state+'/';seed=SOURCE+'outer_n5_'+state+'_verify/'
            for name in ('state.c','state.w','isodata'):
                if files[path+name]!=source[seed+name]:
                    raise ValueError('Qualified CI matched source input differs')
            if files[path+'source.lsj.lbl']!=source[seed+'state.lsj.lbl'] or files[path+'rcsf.inp']!=files[path+'state.c']:
                raise ValueError('Qualified CI source label/CSF alias differs')
            identity(files[path+'state.lsj.lbl'].decode(),state)
            if not files[path+'state.cm'] or files[path+'rci.stdin'].decode()!=controls(state,case):
                raise ValueError('Qualified CI Hamiltonian/root/mixing differs')
            output=(files[path+'rci.stdout']+files[path+'rci.stderr']).decode()
            warnings[case+'_'+state]=[line for line in output.splitlines()
                if any(word in line.lower() for word in ('warning','mass of the nucleus','infinite nuclear'))]
        require_triplet_quintet_order(files[root+'triplet/state.lsj.lbl'].decode(),files[root+'quintet/state.lsj.lbl'].decode())
        for name,lower,operators,same in [('to_ground','ground','E1,M2,E3',False),('to_5S','quintet','M1,E2,M3',True)]:
            path=root+name+'/'
            for tag,state in [('lower',lower),('upper','triplet')]:
                for ext in ('c','w','cm'):
                    if files[path+tag+'.'+ext]!=files[root+state+'/state.'+ext]:
                        raise ValueError('Qualified CI transition state differs')
                if files[path+tag+'.lsj.lbl']!=files[root+state+'/state.lsj.lbl']:
                    raise ValueError('Qualified CI transition labels differ')
            if files[path+'isodata']!=files[root+lower+'/isodata']:
                raise ValueError('Qualified CI transition nucleus differs')
            if (files[path+'rbiotransform.stdin'].decode()!='y\ny\nlower\nupper\ny\n'
                    or files[path+'rtransition.stdin'].decode()!='y\ny\nlower\nupper\n'+operators+'\n'+('n\nn\n' if same else '')):
                raise ValueError('Qualified CI transition controls differ')
        cases[case]=decay(files,root,ci=True)
    return dict(schema='rya1354.qualified-ci-audit.v16',cases=cases,
        DC_versus_SCF_E1=decompose(prior['decay'],cases['DC']),
        Breit_versus_DC_E1=decompose(cases['DC'],cases['Breit']),
        source_SCF=prior['decay']['model'],warnings=warnings,parts=manifest['parts'],
        lifetime_fractional_changes={g:dict(DC_versus_SCF=cases['DC']['model'][g]['model_one_photon_lifetime_ns']/prior['decay']['model'][g]['model_one_photon_lifetime_ns']-1,
            Breit_versus_DC=cases['Breit']['model'][g]['model_one_photon_lifetime_ns']/cases['DC']['model'][g]['model_one_photon_lifetime_ns']-1) for g in ('B','C')},
        scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interpretation='Matched fixed qualified-source same-implementation Hamiltonian/representation sensitivity. Native RCI loading/interpolation retained. No independently implemented relativistic check, calibrated accuracy, experimental fitting or gauge selection. Historical n3 component limitation and failed joint/core release retained.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--bundle',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.write_text(json.dumps(summarize(args.base,args.bundle),indent=2)+'\n')

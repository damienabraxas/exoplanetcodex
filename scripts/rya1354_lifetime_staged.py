"""Fixed2s preconditioning, then full2s release; no lifetime from intermediate states."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import re
import shutil
import sys

from rya1354_lifetime_relax import validate_core_layer,scf_controls
from rya1354_lifetime_ladder import STATES,identity
from rya1354_solver_smoke import run,sha
from rya1354_two_method_summary import scf_diagnostics
from rya1354_triplet_correlation import transitions,require_triplet_quintet_order
from rya1354_wave_probe import decode
from rya1354_radial_controls import initialization_input
from rya1354_mr_inventory import inventory


def require_frozen(source,target,freeze_2s,fixed_max_n=None):
    """Compare against native fixed-load output, including load normalization."""
    a,b=decode(source.read_bytes()),decode(target.read_bytes())
    if set(a)!=set(b):
        raise ValueError('Orbital inventory changed within matched staged model')
    keys=tuple(k for k in a if k[0]<=fixed_max_n) if fixed_max_n is not None else (((1,-1),(2,-1)) if freeze_2s else ((1,-1),))
    for key in keys:
        if a[key]!=b[key]:
            raise ValueError('Declared fixed radial function changed: '+str(key))


def validate_relaxed_n3(root):
    """Verify historical native stopping, fresh inputs, grid, states and n3 inventory."""
    receipt=json.loads((root/'receipt.json').read_text())
    if receipt['status']!='CORE_RELAXED_MODEL_LIFETIME_COMPLETED_REVIEW_REQUIRED':
        raise ValueError('n3 source campaign is not qualified')
    if receipt.get('n') is not None and receipt['n']!=3:
        raise ValueError('n3 source declared a different layer')
    expected={(1,-1),(2,-1),(2,1),(2,-2),(3,-1),(3,1),(3,-2),(3,2),(3,-3)}
    for state in STATES:
        for suffix in ('','_verify'):
            stage=root/(state+suffix)
            for name in ('state.c','state.w','isodata','state.lsj.lbl','rmcdhf.stdout','rmcdhf.stderr','rmcdhf.sum','rwfn.inp'):
                path=stage/name;record=receipt['files'][str(path.relative_to(root))]
                if sha(path)!=record['sha256'] or path.stat().st_size!=record['bytes']:
                    raise ValueError('n3 source artifact differs: '+str(path))
            identity((stage/'state.lsj.lbl').read_text(),state)
            diagnostic=scf_diagnostics((stage/'rmcdhf.stdout').read_text(),(stage/'rmcdhf.stderr').read_text())
            summary=(stage/'rmcdhf.sum').read_text()
            spacing=re.search(r'^ H\s*=\s*([\d.D+\-]+)',summary,re.M)
            points=re.search(r'^ N\s*=\s*(\d+)',summary,re.M)
            if (diagnostic['iteration_limit_hit'] or not diagnostic['all_reported_functions_below_ACC']
                    or spacing is None or points is None or float(spacing[1].replace('D','E'))!=.03125 or int(points[1])!=439):
                raise ValueError('n3 source native criterion/grid differs')
            if set(decode((stage/'state.w').read_bytes()))!=expected:
                raise ValueError('n3 source orbital inventory differs')
            text=(stage/'state.c').read_text();inventory(text,'even' if state=='ground' else 'odd',True)
            if not re.search(r'2s\s*\(\s*1\)',text):
                raise ValueError('n3 source lacks declared2s-hole configurations')
            if any(int(n)>3 for n in re.findall(r'(\d+)[spdfgh]-?\s*\(',text)):
                raise ValueError('n3 source CSF contains higher layer')
        if (root/state/'state.w').read_bytes()!=(root/(state+'_verify')/'rwfn.inp').read_bytes():
            raise ValueError('n3 fresh restart did not use optimized wavefunction')
    return sha(root/'receipt.json')


def prepare_n4_sources(args,receipt):
    sources={}
    for state in STATES:
        source=args.n3_core_relaxed/(state+'_verify');target=args.core/('CV_2s_DC_'+state)
        if (source/'isodata').read_bytes()!=(target/'isodata').read_bytes():
            raise ValueError('n3/n4 nucleus identity differs')
        out=args.out/('initial_'+state);out.mkdir()
        shutil.copyfile(target/'state.c',out/'rcsf.inp');shutil.copyfile(target/'state.c',out/'state.c')
        shutil.copyfile(target/'isodata',out/'isodata');shutil.copyfile(source/'state.w',out/'previous.w')
        receipt['steps'].append(run(args.grasp_bin/'rwfnestimate',initialization_input(3,439,.03125),out,1800))
        text=(out/'rwfnestimate.stdout').read_text()+(out/'rwfnestimate.stderr').read_text()
        if receipt['steps'][-1]['semantic_errors'] or 'All required subshell radial wavefunctions' not in text:
            raise ValueError('n4 native initializer refused')
        if set(decode((out/'rwfn.inp').read_bytes()))!=set(decode((target/'state.w').read_bytes())):
            raise ValueError('Initialized n4 orbital inventory differs from target')
        shutil.copyfile(out/'rwfn.inp',out/'state.w');sources[state]=out
    return sources


def load_baseline(args,out,state,receipt):
    """Same input/grid, all orbitals fixed: reproduce native load normalization."""
    baseline=out/'fixed_load';baseline.mkdir()
    for name in ('isodata','rcsf.inp','rwfn.inp'):
        shutil.copyfile(out/name,baseline/name)
    receipt['steps'].append(run(args.grasp_bin/'rangular','y\n',baseline,1800))
    controls=scf_controls(state,4,1,args.damping,True,radial_iteration_floor=24,
                          both_2p_spectroscopic=True,all_fixed=True)
    receipt['steps'].append(run(args.strict_bin/'rmcdhf',controls,baseline,43200))
    text=(baseline/'rmcdhf.stdout').read_text()+(baseline/'rmcdhf.stderr').read_text()
    if ('All subshell radial wavefunctions are fixed' not in text or
            any(s['semantic_errors'] for s in receipt['steps']) or not (baseline/'rwfn.out').is_file()):
        raise ValueError('Native all-fixed loading control failed')
    if set(decode((baseline/'rwfn.out').read_bytes()))!=set(decode((out/'rwfn.inp').read_bytes())):
        raise ValueError('Native loading control changed orbital inventory')
    return baseline/'rwfn.out'


def optimize(args,phase,state,seed,receipt):
    freeze=phase=='fixed_2s';outer=phase=='outer_n4';selected=None
    for verify in (False,True):
        out=args.out/(phase+'_'+state+('_verify' if verify else ''));out.mkdir()
        source=selected if verify else seed
        for name in ('isodata','state.c','state.w'):
            shutil.copyfile(source/name,out/name)
        shutil.copyfile(out/'state.c',out/'rcsf.inp');shutil.copyfile(out/'state.w',out/'rwfn.inp')
        receipt['steps'].append(run(args.grasp_bin/'rangular','y\n',out,1800))
        baseline=load_baseline(args,out,state,receipt)
        controls=scf_controls(state,4,args.cycles,args.damping,verify,radial_iteration_floor=24,
                              both_2p_spectroscopic=True,freeze_2s=freeze,optimization_nsic=args.optimization_nsic,outer_layer_only=outer)
        receipt['steps'].append(run(args.strict_bin/'rmcdhf',controls,out,43200))
        if any(step['semantic_errors'] for step in receipt['steps']):
            raise ValueError('Native refusal; no staged-state promotion')
        diagnostics=scf_diagnostics((out/'rmcdhf.stdout').read_text(),(out/'rmcdhf.stderr').read_text())
        if diagnostics['iteration_limit_hit'] or not diagnostics['all_reported_functions_below_ACC']:
            raise ValueError('Staged orbital criterion unmet')
        require_frozen(baseline,out/'rwfn.out',freeze,fixed_max_n=3 if outer else None)
        if (out/'rcsf.inp').read_bytes()!=(seed/'state.c').read_bytes() or (out/'isodata').read_bytes()!=(seed/'isodata').read_bytes():
            raise ValueError('Staged CSF/nucleus identity changed')
        for native,extension in (('rwfn.out','w'),('rmix.out','m')):
            shutil.copyfile(out/native,out/('state.'+extension))
        receipt['steps'].append(run(args.grasp_bin/'jj2lsj','state\nn\nn\ny\n',out,600))
        identity((out/'state.lsj.lbl').read_text(),state)
        receipt['diagnostics'][out.name]=diagnostics
        receipt['last_qualified_stage']=out.name
        (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        selected=out
    return selected


def optimize_stages(args,receipt):
    previous=(prepare_n4_sources(args,receipt) if getattr(args,'n3_core_relaxed',None)
              else {state:args.core/('CV_2s_DC_'+state) for state in STATES})
    for phase in receipt['phases']:
        selected={state:optimize(args,phase,state,previous[state],receipt) for state in STATES}
        require_triplet_quintet_order((selected['triplet']/'state.lsj.lbl').read_text(),(selected['quintet']/'state.lsj.lbl').read_text())
        receipt['last_completed_phase']=phase;previous=selected
        (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return selected


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--strict-bin',type=Path,required=True)
    ap.add_argument('--core',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--cycles',type=int,default=1500)
    ap.add_argument('--damping',type=float,default=.5)
    ap.add_argument('--optimization-nsic',type=int,default=20,help='Extra least-self-consistent orbital updates after every outer cycle; verification remains NSIC0')
    ap.add_argument('--n3-core-relaxed',type=Path,help='Explicit qualified n3 relaxed source; initialize n4 natively, optimize outer layer before full release')
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    if args.cycles<1 or not 0<args.damping<1 or args.optimization_nsic<0:
        raise ValueError('Invalid staged numerical controls')
    for attr in ('grasp_bin','strict_bin','core','out'):
        setattr(args,attr,getattr(args,attr).resolve())
    initializer=validate_core_layer(args.core,4)
    relaxed_source=None
    if args.n3_core_relaxed:
        args.n3_core_relaxed=args.n3_core_relaxed.resolve()
        relaxed_source=validate_relaxed_n3(args.n3_core_relaxed)
    args.out.mkdir(parents=True,exist_ok=False)
    drivers=args.out/'drivers';drivers.mkdir();shutil.copyfile(__file__,drivers/'driver.py')
    for name in ('rya1354_lifetime_relax','rya1354_lifetime_ladder','rya1354_solver_smoke',
                 'rya1354_two_method_summary','rya1354_triplet_correlation','rya1354_wave_probe',
                 'rya1354_radial_controls','rya1354_mr_inventory'):
        shutil.copyfile(sys.modules[name].__file__,drivers/(name+'.py'))
    receipt=dict(schema='rya1354.staged-core-relaxation.v11',status='RUNNING',n=4,
        scientific_admission='HOLD',adoption_allowed=False,started_utc=datetime.now(timezone.utc).isoformat(),
        core_receipt_sha256=initializer,steps=[],diagnostics={},cycles=args.cycles,damping=args.damping,
        optimization_nsic=args.optimization_nsic,frozen_record_baseline='NATIVE_ALL_FIXED_SAME_INPUT_GRID_LOADING_CONTROL',
        n3_core_relaxed_receipt_sha256=relaxed_source,
        phases=['outer_n4','released_2s'] if relaxed_source else ['fixed_2s','released_2s'],accepted_lifetime=None,
        model='Same n4 full2s-SD/DC CSFs. All3states first vary2p*/3*/4* with1s/2s fixed, strict/fresh verification; then vary2*/3*/4* with1s fixed, strict/fresh verification. Explicit2p-/2p spectroscopic, NSOLV24, original439grid/ACCY. No observed/failed-wave inputs; no lifetime from preconditioner.')
    if relaxed_source:
        receipt['model']='Same n4 full2s-SD/DC CSFs; qualified historical core-relaxed n3 wavefunctions plus native hydrogenic missing-n4 initialization. All3states optimize4* withn<=3 fixed/verified, then2*/3*/4* with1s fixed/verified. Strict original439grid/ACCY/freshNSIC0 and same-input native loaded baselines. Explicit2p pair, NSOLV24, no observed/failed-wave inputs, no intermediate lifetime. Historical n3 spectroscopic-component limitation retained until final requalification.'
    try:
        selected=optimize_stages(args,receipt)
        for name,lower,operators,same in [('to_ground',selected['ground'],'E1,M2,E3',False),
                                         ('to_5S',selected['quintet'],'M1,E2,M3',True)]:
            receipt['steps']+=transitions(args.grasp_bin,lower,selected['triplet'],args.out/name,operators,same)
        receipt['status']='CORE_RELAXED_MODEL_LIFETIME_COMPLETED_REVIEW_REQUIRED'
        # Terminal native completion still requires independent offline lifetime audit.
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()

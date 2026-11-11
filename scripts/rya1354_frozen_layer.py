"""Separate frozen-prior-layer sensitivity; never admits the failed full release."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import shutil

from rya1354_native_parts import load_parts
from rya1354_lifetime_summary import scf_state
from rya1354_lifetime_relax import scf_controls, validate_core_layer
from rya1354_triplet_correlation import transitions, require_triplet_quintet_order
from rya1354_wave_probe import decode

PARENT = 'runs/lifetime_relax_n4_from_relaxed_n3_v7/'
TARGET = 'runs/lifetime_core_n4_v1/'
STATES = ('ground', 'quintet', 'triplet')
STATUS = 'FROZEN_PRIOR_LAYER_LIFETIME_COMPLETED_REVIEW_REQUIRED'


def digest(value):
    return hashlib.sha256(value).hexdigest()


def validate(files, target):
    """Replay six closed states from original bytes, independently of parent verdict."""
    raw = files[PARENT+'receipt.json']
    receipt = json.loads(raw)
    if (receipt['status'] != 'FAILED_RETAINED' or receipt['n'] != 4
            or receipt['last_completed_phase'] != 'outer_n4'
            or receipt['last_qualified_stage'] != 'outer_n4_triplet_verify'
            or receipt['phases'] != ['outer_n4', 'released_2s']
            or receipt['optimization_nsic'] != 0):
        raise ValueError('Unexpected source campaign scope')
    if digest(target[TARGET+'receipt.json']) != receipt['core_receipt_sha256']:
        raise ValueError('Target core receipt differs')
    # Every supplied parent artifact must retain its terminal receipt identity.
    for name, record in receipt['files'].items():
        if PARENT+name not in files:
            if Path(name).name.startswith('mcp.'):
                continue
            raise ValueError('Missing source artifact: '+name)
        value = files[PARENT+name]
        if digest(value) != record['sha256'] or len(value) != record['bytes']:
            raise ValueError('Source artifact hash differs: '+name)
    states = {}
    for state in STATES:
        for verify in (False, True):
            case = 'outer_n4_'+state+('_verify' if verify else '')
            path = PARENT+case
            diagnostic = scf_state(files, path, state)
            expected = scf_controls(state,4,1500,.5,verify,radial_iteration_floor=24,
                both_2p_spectroscopic=True,optimization_nsic=0,outer_layer_only=True)
            if files[path+'/rmcdhf.stdin'].decode() != expected:
                raise ValueError('Outer-layer control differs')
            for name in ('isodata','rcsf.inp','rwfn.inp'):
                if files[path+'/'+name] != files[path+'/fixed_load/'+name]:
                    raise ValueError('Loading baseline input differs')
            baseline = scf_controls(state,4,1,.5,True,radial_iteration_floor=24,
                both_2p_spectroscopic=True,all_fixed=True)
            if files[path+'/fixed_load/rmcdhf.stdin'].decode() != baseline:
                raise ValueError('All-fixed loading control differs')
            for adopted, native in (('state.c','rcsf.inp'),('state.w','rwfn.out'),('state.m','rmix.out')):
                if files[path+'/'+adopted] != files[path+'/'+native]:
                    raise ValueError('Selected state differs from native output')
            for name in ('state.c','isodata'):
                if files[path+'/'+name] != target[TARGET+'CV_2s_DC_'+state+'/'+name]:
                    raise ValueError('Target configuration/nucleus differs')
            a = decode(files[path+'/fixed_load/rwfn.out'])
            b = decode(files[path+'/state.w'])
            fixed = [k for k in a if k[0] <= 3]
            if set(a) != set(b) or len(a) != 16 or len(fixed) != 9 or any(a[k] != b[k] for k in fixed):
                raise ValueError('Frozen earlier-layer records changed')
            seed = PARENT+('outer_n4_'+state if verify else 'initial_'+state)+'/state.w'
            if files[path+'/rwfn.inp'] != files[seed]:
                raise ValueError('Initializer/fresh restart differs')
            states[case] = dict(diagnostics=diagnostic,exact_frozen_n_le_3_records=9)
    require_triplet_quintet_order(files[PARENT+'outer_n4_triplet_verify/state.lsj.lbl'].decode(),
                                 files[PARENT+'outer_n4_quintet_verify/state.lsj.lbl'].decode())
    return dict(schema='rya1354.frozen-layer-source-audit.v12',states=states,
        parent_receipt_sha256=digest(raw),parent_status=receipt['status'],parent_error=receipt['error'],
        fully_released_qualified=False,scientific_admission='HOLD',adoption_allowed=False,
        exposure='EXPOSED_NONBLIND',accepted_lifetime=None,uncertainty=None,covariance=None,
        calibrated_sigma=None,generator_sha256=digest(Path(__file__).read_bytes()))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,required=True)
    ap.add_argument('--core',type=Path,required=True)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args = ap.parse_args()
    if platform.node().lower().split('.')[0] != 'sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    validate_core_layer(args.core,4)
    source = {PARENT+str(p.relative_to(args.source)):p.read_bytes()
              for p in args.source.rglob('*') if p.is_file() and not p.name.startswith('mcp.')}
    target = {TARGET+str(p.relative_to(args.core)):p.read_bytes()
              for p in [args.core/'receipt.json', *[args.core/('CV_2s_DC_'+s)/n
                        for s in STATES for n in ('state.c','isodata')]]}
    audit = validate(source,target)
    args.out.mkdir(parents=True,exist_ok=False)
    drivers = args.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*.py'):
        shutil.copyfile(p,drivers/p.name)
    receipt = dict(status='RUNNING',source_audit=audit,steps=[],
        started_utc=datetime.now(timezone.utc).isoformat(),scientific_admission='HOLD',
        adoption_allowed=False,accepted_lifetime=None,
        model='Separate full2s-SD/DC n4 sensitivity: qualified core-adapted n<=3 fixed; only n4 optimized. Not full core relaxation. Historical n3 component constraint limitation retained. No observed-energy substitution.')
    try:
        selected = {}
        for state in STATES:
            dest = args.out/state;dest.mkdir();selected[state] = dest
            for name in ('state.c','state.w','state.m','state.lsj.lbl','isodata'):
                (dest/name).write_bytes(source[PARENT+'outer_n4_'+state+'_verify/'+name])
        for name, lower, operators, same in [('to_ground','ground','E1,M2,E3',False),('to_5S','quintet','M1,E2,M3',True)]:
            receipt['steps'] += transitions(args.grasp_bin,selected[lower],selected['triplet'],args.out/name,operators,same)
        receipt['status'] = STATUS
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc))
        raise
    finally:
        receipt['files'] = {str(p.relative_to(args.out)):dict(sha256=digest(p.read_bytes()),bytes=p.stat().st_size)
                            for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__ == '__main__':
    main()

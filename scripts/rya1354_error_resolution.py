"""Replay convergence and physical-model error controls, without an accuracy claim."""
import hashlib
import json
import re
from pathlib import Path

from rya1354_strict_scf import strict_source
from rya1354_two_method_summary import capsule_files,scf_diagnostics,verify_artifacts
from rya1354_forbidden_summary import parse as forbidden
from rya1354_lifetime_validation import parse_table,model_lifetimes

ROOT=Path(__file__).resolve().parents[1]


def qualified_scf(files,root):
    receipt=json.loads(files[root+'/receipt.json'])
    if receipt['status']!='NATIVE_STOP_COMPLETED_REVIEW_REQUIRED':
        raise ValueError('Strict SCF continuation did not complete')
    labels=files[root+'/state.lsj.lbl'].decode()
    blocks=re.findall(r'^\s*\d+\s+(\d+)\s+([+-])\s+[^\n]+\n([^\n]+)',labels,re.M)
    expected=([('2','3s_5S','-')] if receipt.get('state_kind')=='quintet3s' else
              [('0','_3P','+'),('0','_1S','+'),('1','_3P','+'),('2','_3P','+'),('2','_1D','+')])
    if len(blocks)!=len(expected) or any(j!=ej or p!=ep or not line.rstrip().endswith(term)
                                       for (j,p,line),(ej,term,ep) in zip(blocks,expected)):
        raise ValueError('Strict SCF state identity changed')
    d=scf_diagnostics(files[root+'/rmcdhf.stdout'].decode(),files[root+'/rmcdhf.stderr'].decode(),
                      function_accy=receipt.get('native_ACC',.0625**6))
    if d['iteration_limit_hit'] or not d['all_reported_functions_below_ACC']:
        raise ValueError('Native orbital-change criterion is not satisfied')
    d['effective_stopping_rule']='Isolated strict variant requires native SCNSTY<=ACCY; energy-only stop removed'
    d['measure']='Native CONSIS: maximum radial-function change multiplied by sqrt(UCF occupation); not a differential-equation residual or uncertainty'
    return d


def key(row):
    return row['operator'],row['upper'],row['lower']


def rate_comparison(rows,baseline):
    old={key(row):row for row in baseline}
    return [dict(operator=r['operator'],upper=r['upper'],lower=r['lower'],
                 A_B_s1=r['A_B_s1'],native_energy_cm1=r['native_energy_cm1'],
                 fractional_A_change=None if old[key(r)]['A_B_s1']==0 else r['A_B_s1']/old[key(r)]['A_B_s1']-1)
            for r in rows if r['upper']=='1D2' and r['lower'] in ('3P1','3P2') and r['operator'] in ('M1','E2')]


def summarize(destination=None):
    base=ROOT/'data/audit/rya1354'
    paths=[base/name for name in ('strict_controls_n34_v1.tar.gz','strict_lifetime_balance_v2.tar.gz',
                                  'strict_extensions_v3.tar.gz','strict_grid_correction_v4.tar.gz')]
    parts=[base/f'ground_n5_hamiltonian_v1_parts.{case}.tar.gz' for case in ('DC','Breit_state','Breit_multipoles')]
    files={}
    for path in paths:
        files.update(capsule_files(path))
    for part in parts:
        files.update(capsule_files(part,verify=False))
    verify_artifacts(files)
    paths+=parts
    build='toolchains/grasp-strict-scf-v1'
    build_receipt=json.loads(files[build+'/receipt.json'])
    if build_receipt['upstream_dirty']:
        raise ValueError('Strict diagnostic was not built from the clean pinned source')
    original=files[build+'/scf.original.f90'].decode()
    if strict_source(original)!=files[build+'/scf.f90'].decode():
        raise ValueError('Strict solver contains undeclared source changes')
    control='runs/strict_control_'
    for variant in ('original_v1','variant_v1'):
        record=json.loads(files[control+variant+'/receipt.json'])
        if record['status']!='FAILED_RETAINED' or record['steps'][0]['semantic_errors']!=['MAXIMUM ITERATIONS IN SCF EXCEEDED']:
            raise ValueError('One-cycle control did not retain the expected iteration-limit failure')
    for name in ('rwfn.out','rmix.out'):
        if files[control+'original_v1/'+name]!=files[control+'variant_v1/'+name]:
            raise ValueError('Solver update differs before stopping decisions diverge')
    for variant,expected in [('original_v1',build_receipt['original_executable_sha256']),
                             ('variant_v1',build_receipt['files']['rmcdhf']['sha256'])]:
        if json.loads(files[control+variant+'/receipt.json'])['steps'][0]['sha256']!=expected:
            raise ValueError('Control executable differs from the declared build')
    old=json.loads((base/'source_curation_v2/forbidden_convergence_v2.json').read_text())
    numerical=[]
    for n in (3,4,5):
        root=f'runs/strict_n{n}_v1' if n<5 else 'runs/strict_n5_verify_v1'
        multipoles='strict_multipoles_n34_v1' if n<5 else 'strict_n5_multipoles_v1'
        channels=forbidden(files[f'runs/{multipoles}/n{n}_multipoles/lower.upper.t'].decode())
        numerical.append(dict(n=n,native_SCF_diagnostics=qualified_scf(files,root),channels=channels,
                              versus_energy_stopped=rate_comparison(channels,old['layers'][n-3]['channels'])))
    # The first custom-grid driver assumed ACCY was recomputed. Native GETSCD
    # initializes it before grid input; reject that metadata instead of laundering
    # its one-cycle completion as satisfaction of the finer requested criterion.
    try:
        qualified_scf(files,'runs/strict_n3_halfgrid_v1')
    except ValueError as exc:
        if 'native ACCY prompt' not in str(exc):
            raise
    else:
        raise ValueError('Expected custom-grid criterion mismatch absent')
    fine=qualified_scf(files,'runs/strict_n3_halfgrid_v2')
    grid_summary=files['runs/strict_n3_halfgrid_v2/rmcdhf.sum'].decode()
    grid_h=float(re.search(r'^ H\s*=\s*([\d.D+\-]+)',grid_summary,re.M)[1].replace('D','E'))
    grid_n=int(re.search(r'^ N\s*=\s*(\d+)',grid_summary,re.M)[1])
    if grid_h!=.03125 or grid_n!=439:
        raise ValueError('Requested finer grid not established in native summary')
    fine_channels=forbidden(files['runs/strict_halfgrid_multipoles_v2/n3_multipoles/lower.upper.t'].decode())
    core={}
    dc=None
    for case in ('DC','Breit_static_limit','CV_2s_one_hole','CV_2s_two_holes','CV_2s_two_holes_Breit_static_limit'):
        run='ground_core_v2' if case=='CV_2s_two_holes_Breit_static_limit' else 'ground_core_v1'
        channels=forbidden(files[f'runs/{run}/{case}_multipoles/lower.upper.ct'].decode())
        if case=='DC':
            dc=channels
        core[case]=dict(channels=channels,versus_DC=rate_comparison(channels,dc))
    # Never count the reversed, upward-to-quintet v1 table as a triplet decay.
    try:
        parse_table(files['runs/triplet_correlation_v1/to_5S/lower.upper.t'].decode(),'to_5S')
    except ValueError:
        pass
    else:
        raise ValueError('Expected inverted physical ordering was not detected')
    rows=[]
    for name in ('to_ground','to_5S'):
        rows+=parse_table(files[f'runs/triplet_correlation_v2/{name}/lower.upper.t'].decode(),name)
    lifetime=model_lifetimes(rows)
    parent='runs/triplet_correlation_v2/triplet_3s'
    td=scf_diagnostics(files[parent+'/rmcdhf.stdout'].decode(),files[parent+'/rmcdhf.stderr'].decode())
    if td['iteration_limit_hit'] or not td['all_reported_functions_below_ACC']:
        raise ValueError('Correlated triplet fails orbital criterion')
    previous=base/'source_curation_v2/lifetime_validation.json'
    oldlife=json.loads(previous.read_text())
    experiment=oldlife['experiment']
    result=dict(schema='rya1354.error-resolution.v3',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
                prediction_origin='THIS_UNPUBLISHED_THEORETICAL_CALCULATION',candidate_doi=None,
                capsules={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                software_control=dict(one_cycle_orbitals_and_mixing_byte_identical=True,
                                      declared_patch='Remove alternative energy-only stop, exactly one statement',
                                      upstream_build_unchanged=True),
                numerical_continuations=numerical,
                n5_acceleration=qualified_scf(files,'runs/strict_n5_damp05_v1'),
                interrupted_runs=json.loads(files['numerical_interruptions_v3.json']),
                finer_grid=dict(diagnostics=fine,channels=fine_channels,
                                native_H=grid_h,native_N=grid_n,
                                versus_original_grid=rate_comparison(fine_channels,numerical[0]['channels']),
                                rejected_v1='Native ACCY is initialized before custom grid input and not automatically recomputed. Original driver/receipt incorrectly declared H**6 after changing H; retained/rejected. v2 explicitly revises ACCY; diagnostics now compare declared against native prompt.'),
                numerical_limitation='Same-layer continuations leave inherited lower-layer functions frozen; no complete rebuilt ladder, radial grid or physical convergence claim',
                ground_core_controls=core,
                ground_core_warning='Fixed n3 orbitals and RCI interpolation grid; DC numerical control precedes physical cases. No independent-error or laboratory-gf interpretation. Stationary-mass divide-by-zero and underflow/denormal warnings remain explicit.',
                lifetime=dict(model=lifetime,channels=rows,triplet_SCF_diagnostics=td,
                              quintet_SCF_diagnostics=qualified_scf(files,'runs/strict_quintet_n3_v1'),
                              experiment=experiment,
                              prior_lifetime_summary_sha256=hashlib.sha256(previous.read_bytes()).hexdigest(),
                              residual_ns={g:lifetime[g]['model_one_photon_lifetime_ns']-experiment['value_ns'] for g in ('B','C')},
                              rejected_v1='Triplet correlated while destination quintet uncorrelated; native levels invert and to_5S table is reverse-direction. Raw receipt completion does not qualify lifetime. Retained and rejected; v2 uses matching n3 VV-SD destination.',
                              limitation='Balanced n3 valence model still freezes n<3; omitted core, larger active spaces, independent method, grid and higher-order multiphoton contributions unqualified'),
                uncertainty=None,covariance=None,calibrated_sigma=None,chi_square=None,adoption_allowed=False)
    rows4=[]
    for name in ('to_ground','to_5S'):
        rows4+=parse_table(files[f'runs/triplet_correlation_n4_v1/{name}/lower.upper.t'].decode(),name)
    model4=model_lifetimes(rows4)
    parent4='runs/triplet_correlation_n4_v1/triplet_3s'
    d4=scf_diagnostics(files[parent4+'/rmcdhf.stdout'].decode(),files[parent4+'/rmcdhf.stderr'].decode())
    if d4['iteration_limit_hit'] or not d4['all_reported_functions_below_ACC']:
        raise ValueError('n4 triplet fails orbital criterion')
    result['lifetime_n4']=dict(model=model4,channels=rows4,triplet_SCF_diagnostics=d4,
                               quintet_SCF_diagnostics=qualified_scf(files,'runs/strict_quintet_n4_v1'),
                               residual_ns={g:model4[g]['model_one_photon_lifetime_ns']-experiment['value_ns'] for g in ('B','C')},
                               limitation='Native stopping qualified for varied n4 functions; inherited frozen earlier layers differ slightly from a full rebuilt strict ladder. Gauges/active-space residuals are not uncertainties.')
    hamiltonian={}
    dc5=None
    for case in ('DC','Breit_static_limit'):
        rows5=forbidden(files[f'runs/ground_n5_hamiltonian_v1/{case}_multipoles/lower.upper.ct'].decode())
        if case=='DC':
            dc5=rows5
        hamiltonian[case]=dict(channels=rows5,versus_DC=rate_comparison(rows5,dc5),
                              versus_strict_SCF=rate_comparison(rows5,numerical[-1]['channels']))
    result['ground_n5_hamiltonian']=hamiltonian
    (destination or base/'source_curation_v2/error_resolution_v3.json').write_text(json.dumps(result,indent=2)+'\n')
    print('lifetime ns:',{g:v['model_one_photon_lifetime_ns'] for g,v in lifetime.items()})


if __name__=='__main__':
    summarize()

"""Independent radial metadata replay and energy/amplitude sensitivity accounting."""
import ast
import hashlib
import json
import math
from pathlib import Path
import re

from rya1354_native_parts import load_parts
from rya1354_fac_potential_audit import verify_directory
from rya1354_fac_radial_operator import pairs,raw_rows,review_files
from rya1354_fac_radial_angular_remaining_metadata import hold


def quotient(a,b):
    if b==0:return None
    value=a/b
    if not math.isfinite(value):raise ValueError('Nonfinite diagnostic ratio')
    return value


def factorization(first,second):
    """Pinned FAC formula has A proportional to E**(2|m|+1) times amplitude**2."""
    if first['multipole']!=second['multipole']:raise ValueError('Operators differ')
    power=2*abs(first['multipole'])+1
    energy=(second['energy_au']/first['energy_au'])**power
    amp=quotient(second['amplitude'],first['amplitude'])
    actual=quotient(second['A_s1'],first['A_s1'])
    predicted=None if amp is None else energy*amp*amp
    return dict(energy_power=power,energy_factor=energy,amplitude_squared_factor=None if amp is None else amp*amp,
        predicted_rate_ratio=predicted,actual_rate_ratio=actual,
        floating_factorization_residual=None if actual is None or predicted is None else actual-predicted,
        scope='Algebraic decomposition of exposed model outputs, not an uncertainty or a substituted rate')


def channels(entry,suffix):
    result={}
    for line in entry['raw'][suffix].splitlines():
        if not line.strip() or line.startswith('#'):continue
        f=line.split()
        if len(f)!=15:raise ValueError('Parent raw channel columns differ')
        up,low,ju,jl,nele,m,g,mode,status,nz,nmk=map(int,f[:11]);e,s,gf,a=map(float,f[11:])
        lower=entry['levels'][str(low)]
        if (up,ju,jl,nele,g,mode,status,nmk)!=(entry['upper'],2,lower['two_J'],8,2 if suffix=='B' else 1,0,0,0):
            raise ValueError('Parent raw channel context differs')
        if not all(math.isfinite(x) for x in (e,s,gf,a)) or e<=0 or min(gf,a)<0:raise ValueError('Parent raw channel is invalid')
        key=f"{lower['dominant']}|J2={jl}|parity={lower['parity']}|m={m}"
        if key in result:raise ValueError('Ambiguous nominal physical channel join')
        result[key]=dict(lower=low,lower_identity=lower,multipole=m,energy_au=e,amplitude=s,gf=gf,A_s1=a)
    if len(result)!=13:raise ValueError('Parent full channel inventory differs')
    total=math.fsum(v['A_s1'] for v in result.values())
    return dict(total_A_s1=total,channels={k:dict(**v,rate_fraction=v['A_s1']/total) for k,v in result.items()})


def paired_gauges(b,c):
    if set(b)!=set(c):raise ValueError('Radial gauge pair inventory differs')
    result=[]
    for pair,br in b.items():
        cr=c[pair]
        if any(br[k]!=cr[k] for k in ('used_au','eps1_au','eps2_au','radial_NR_length')):
            raise ValueError('Gauge switch changed orbital/frequency/length comparator')
        ratio=quotient(cr['radial_FR'],br['radial_FR']);gap=abs(br['eps1_au']-br['eps2_au'])
        result.append(dict(n1=pair[0][0],kappa1=pair[0][1],n2=pair[1][0],kappa2=pair[1][1],
            photon_energy_au=br['used_au'],orbital_gap_au=gap,B=br['radial_FR'],C=cr['radial_FR'],
            NR_length=br['radial_NR_length'],signed_C_over_B=ratio,abs_C_over_B=None if ratio is None else abs(ratio),
            orbital_gap_over_photon=gap/br['used_au'],
            abs_gauge_ratio_times_photon_over_gap=None if ratio is None else abs(ratio)*br['used_au']/gap,
            B_over_NR_length=quotient(br['radial_FR'],br['radial_NR_length'])))
    return result


def summarize(base):
    digest=lambda data:hashlib.sha256(data).hexdigest()
    for name in ('fac_radial_operator_failed_v42','fac_radial_operator_complete_v43'):verify_directory(base/name)
    failed,failed_manifest=load_parts(base/'fac_radial_operator_failed_v42')
    files,manifest=load_parts(base/'fac_radial_operator_complete_v43')
    oldroot='runs/fac_radial_operator_v42/';root='runs/fac_radial_operator_v43/'
    old=json.loads(failed[oldroot+'receipt.json']);receipt=json.loads(files[root+'receipt.json']);hold(old);hold(receipt)
    index_raw=files[root+'index.json'];index=json.loads(index_raw);hold(index)
    plan=json.loads((base.parents[2]/'config/atomic/rya1354/lifetime_plan.v43.json').read_text())
    previous_manifest=json.loads((base/'fac_radial_rates_complete_v41/manifest.json').read_text())
    if (receipt['schema']!='rya1354.radial-operator.v43' or index['schema']!='rya1354.radial-operator-input.v43'
        or digest(index_raw)!=plan['index_sha256'] or digest(index_raw)!=receipt['index_sha256']
        or receipt['status']!='FAC_RADIAL_OPERATOR_COMPLETED_REVIEW_REQUIRED' or old['status']!='FAILED_RETAINED'
        or receipt['retained_failure_receipt_sha256']!=digest(failed[oldroot+'receipt.json'])
        or index['retained_failure_receipt_sha256']!=digest(failed[oldroot+'receipt.json'])
        or index['parent_receipt_sha256']!=previous_manifest['files']['runs/fac_radial_rates_v41/receipt.json']['sha256']
        or index['parent_index_sha256']!=previous_manifest['files']['runs/fac_radial_rates_v41/index.json']['sha256']):
        raise ValueError('Diagnostic parent/index/status chain differs')
    radial=files[root+'source_documents/faclib/radial.c'].decode()
    floor=float(re.search(r'static double _awmin = ([0-9.e+-]+);',radial)[1])
    if floor!=1e-4 or index['photon_grid_min']!=plan['photon_grid_min']:
        raise ValueError('Declared/native photon-grid minimum differs')
    if 'if (0 == strcmp(s, "radial:awmin")) {\n    _awmin = dp;\n    return;' not in radial:
        raise ValueError('Photon-only source option semantics differ')
    if receipt['probe']['libfac_sha256']!=index['probe']['libfac_sha256'] or receipt['probe']['upstream_binary_sha256']!=index['probe']['upstream_binary_sha256']:
        raise ValueError('Original library/executable differs')
    for n in receipt['canary']['files']:
        if files[root+'canary_original/'+n]!=files[root+'canary_diagnostic/'+n]:raise ValueError('Native frontend control differs')
    oldindex=json.loads(failed[oldroot+'index.json'])
    expected_controls={'ground_coherent10_B_own','ground_coherent10_B_paired'}
    if set(old['cases'])!=expected_controls or receipt['previous_successful_controls_exact']!={n:True for n in expected_controls}:
        raise ValueError('Retained successful control scope differs')
    for label in expected_controls:
        p=oldroot+label+'/';subset={n[len(p):]:data for n,data in failed.items() if n.startswith(p)}
        entry=oldindex['cases']['ground_coherent10'];energy=entry['energy_au'] if label.endswith('_own') else entry['paired_energy_au']
        replay=review_files(subset,entry,oldindex,2,energy)
        if replay!=old['cases'][label] or replay!=receipt['cases'][label] or failed[p+'radial.tsv']!=files[root+label+'/radial.tsv']:
            raise ValueError('Prior successful native outputs changed')
    result=dict(schema='rya1354.radial-operator-audit.v43',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
        adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        probes={},gauges={},channel_sensitivity={},potential_sensitivity={},common_energy_radial_comparison={})
    parsed={}
    for case,entry in index['cases'].items():
        result['gauges'][case]={};result['channel_sensitivity'][case]={}
        for kind,energy in [('own',entry['energy_au']),('paired',entry['paired_energy_au']),('orbital_gap',0.0)]:
            for suffix,gauge in [('B',2),('C',1)]:
                label=case+'_'+suffix+'_'+kind;p=root+label+'/'
                subset={n[len(p):]:data for n,data in files.items() if n.startswith(p)}
                reviewed=review_files(subset,entry,index,gauge,energy)
                if reviewed!=receipt['cases'][label]:raise ValueError('Independent native diagnostic metadata differs')
                parsed[label]=raw_rows(subset['radial.tsv'].decode(),gauge,energy,index['constants'])
                result['probes'][label]=dict(pairs=reviewed['pairs'],physical_orbitals_unchanged=True,
                    protected_potential_sha256=reviewed['potential']['protected_sha256'])
            result['gauges'][case][kind]=paired_gauges(parsed[case+'_B_'+kind],parsed[case+'_C_'+kind])
        for suffix in ('B','C'):
            path='runs/fac_radial_rates_v41/'+case+'_'+suffix+'/raw.tsv'
            if digest(entry['raw'][suffix].encode())!=previous_manifest['files'][path]['sha256']:
                raise ValueError('Original rate channels differ')
        b=channels(entry,'B');c=channels(entry,'C')
        if set(b['channels'])!=set(c['channels']):raise ValueError('Gauge channel identity join differs')
        result['channel_sensitivity'][case]=dict(B=b,C=c,
            gauge_factorization={k:factorization(v,c['channels'][k]) for k,v in b['channels'].items()})
    expected={case+'_'+g+'_'+kind for case in index['cases'] for g in ('B','C') for kind in ('own','paired','orbital_gap')}
    if set(receipt['cases'])!=expected:raise ValueError('Incomplete cold-probe scope')
    for scale in (10,100):
        ground='ground_coherent'+str(scale);equal='equal_coherent'+str(scale)
        result['potential_sensitivity'][str(scale)]={}
        result['common_energy_radial_comparison'][str(scale)]={}
        for g in ('B','C'):
            a=result['channel_sensitivity'][ground][g]['channels'];b=result['channel_sensitivity'][equal][g]['channels']
            if set(a)!=set(b):raise ValueError('Nominal physical channel join differs across potentials')
            result['potential_sensitivity'][str(scale)][g]={k:factorization(v,b[k]) for k,v in a.items()}
            common=[]
            for energy_owner,ground_kind,equal_kind in [('ground','own','paired'),('equal','paired','own')]:
                x=parsed[ground+'_'+g+'_'+ground_kind];y=parsed[equal+'_'+g+'_'+equal_kind]
                for pair in pairs():
                    if x[pair]['used_au']!=y[pair]['used_au']:raise ValueError('Common photon-energy comparison differs')
                    common.append(dict(energy_owner=energy_owner,n1=pair[0][0],kappa1=pair[0][1],n2=pair[1][0],kappa2=pair[1][1],
                        photon_energy_au=x[pair]['used_au'],equal_over_ground_radial=quotient(y[pair]['radial_FR'],x[pair]['radial_FR'])))
            result['common_energy_radial_comparison'][str(scale)][g]=common
    result['retained_failure']=dict(status=old['status'],qualified_cases=len(old['cases']),native_default_awmin=floor,
        source_inferred_reason='Requested orbital-gap photon point below unchanged default floor; exact-grid guard rejects the clamp',
        corrected_awmin=index['photon_grid_min'],previous_successful_controls_exact=True)
    return result


if __name__=='__main__':print(json.dumps(summarize(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))

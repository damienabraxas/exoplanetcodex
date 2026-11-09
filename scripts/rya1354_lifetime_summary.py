"""Strict native replay of balanced lifetime, core and independent E1 controls."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

from rya1354_two_method_summary import capsule_files,verify_artifacts,scf_diagnostics,NUM
from rya1354_lifetime_ladder import identity
from rya1354_lifetime_validation import parse_table,model_lifetimes
from rya1354_triplet_correlation import require_triplet_quintet_order
from rya1354_strict_mchf import strict_source as strict_mchf_source
from rya1354_lifetime_mchf import validate_hf_seed
from rya1354_wave_probe import decode as decode_wave,difference as wave_difference

ROOT=Path(__file__).resolve().parents[1]


def load_parts(directory):
    manifest=json.loads((directory/'manifest.json').read_text())
    if manifest.get('schema')!='rya1354.native-evidence-parts.v1' or not manifest['parts']:
        raise ValueError('Unsupported or empty native evidence manifest')
    files={}
    for name,record in manifest['parts'].items():
        if Path(name).name!=name:
            raise ValueError('Part filename is not local to the declared bundle')
        part=directory/name
        if hashlib.sha256(part.read_bytes()).hexdigest()!=record['sha256'] or part.stat().st_size!=record['bytes']:
            raise ValueError('Part hash or size differs: '+name)
        for key,value in capsule_files(part,verify=False).items():
            if key in files and files[key]!=value:
                raise ValueError('Conflicting native artifact across parts')
            files[key]=value
    if set(files)!=set(manifest['files']):
        raise ValueError('Missing/extra native artifact in bundle')
    for name,value in files.items():
        record=manifest['files'][name]
        if hashlib.sha256(value).hexdigest()!=record['sha256'] or len(value)!=record['bytes']:
            raise ValueError('Native member differs from bundle manifest')
    verify_artifacts(files,external_restart_caches=manifest.get('external_restart_caches'))
    return files,manifest


def scf_state(files,path,state):
    identity(files[path+'/state.lsj.lbl'].decode(),state)
    diagnostic=scf_diagnostics(files[path+'/rmcdhf.stdout'].decode(),files[path+'/rmcdhf.stderr'].decode())
    if diagnostic['iteration_limit_hit'] or not diagnostic['all_reported_functions_below_ACC']:
        raise ValueError('Unqualified state promoted: '+path)
    summary=files[path+'/rmcdhf.sum'].decode()
    h=float(re.search(r'^ H\s*=\s*([\d.D+\-]+)',summary,re.M)[1].replace('D','E'))
    n=int(re.search(r'^ N\s*=\s*(\d+)',summary,re.M)[1])
    if h!=.03125 or n!=439:
        raise ValueError('Consistent finer grid not established in native summary')
    diagnostic.update(native_H=h,native_N=n,effective_stop='Declared strict native orbital criterion; no energy-only stop')
    return diagnostic


def decay(files,root,*,ci=False):
    rows=[]
    for name in ('to_ground','to_5S'):
        rows+=parse_table(files[root+name+'/lower.upper.'+('ct' if ci else 't')].decode(),name,include_strengths=True)
    model=model_lifetimes(rows)
    ground=[r for r in rows if r['destination']=='to_ground' and r['operator']=='E1' and r['lower_ordinal']==1]
    if {r['lower_J'] for r in ground}!={0,1,2} or len(ground)!=3:
        raise ValueError('Ground3P E1 component inventory incomplete')
    strongest,=[r for r in ground if r['lower_J']==2]
    return dict(model=model,channels=rows,dominant_3P_E1_A_s1={g:sum(r[f'A_{g}_s1'] for r in ground) for g in ('B','C')},
                ground_J2_E1=strongest)


def change(current,previous):
    a=current['ground_J2_E1'];b=previous['ground_J2_E1']
    energy=(a['native_energy_cm1']/b['native_energy_cm1'])**3
    strength=a['native_S_B_au']/b['native_S_B_au']
    ratio=a['A_B_s1']/b['A_B_s1']
    # Printed native values have six significant figures; this checks their E1
    # dimensional identity, not an accuracy threshold or experimental energy fit.
    if not math.isclose(ratio,energy*strength,rel_tol=5e-5):
        raise ValueError('Native E1 rate/energy/strength decomposition inconsistent')
    return dict(A_B_ratio=ratio,native_energy_cubed_ratio=energy,native_S_B_ratio=strength,
                scope='THEORY_TO_THEORY_DIAGNOSTIC_NO_OBSERVED_ENERGY_SUBSTITUTION')


def mchf_rate(text):
    if '2p(4)_3P' not in text or '2p(3)4S3.3s_3S' not in text:
        raise ValueError('Independent LS parent/term identity absent')
    conv=lambda x:float(x.replace('D','E'))
    energies=re.findall(r'([\d.]+) CM-1\s+([\d.]+) ANGS\(VAC\)',text)
    length=re.findall(rf'E1\s+length:\s+S =\s+({NUM})\s+GF =\s+({NUM})\s+AKI =\s+({NUM})',text)
    velocity=re.findall(rf'velocity:\s+S =\s+({NUM})\s+GF =\s+({NUM})\s+AKI =\s+({NUM})',text)
    if len(energies)!=1 or len(length)!=1 or len(velocity)!=1:
        raise ValueError('Independent dominant E1 table missing/ambiguous')
    e,w=map(conv,energies[0]);sl,gl,al=map(conv,length[0]);sv,gv,av=map(conv,velocity[0])
    if not all(math.isfinite(x) and x>0 for x in (e,w,sl,gl,al,sv,gv,av)):
        raise ValueError('Invalid independent native rate')
    return dict(native_energy_cm1=e,native_vacuum_A=w,S_length_au=sl,S_velocity_au=sv,gf_length=gl,gf_velocity=gv,
                partial_E1_A_length_s1=al,partial_E1_A_velocity_s1=av,
                inverse_partial_rate_ns=dict(length=1e9/al,velocity=1e9/av),
                scope='LS_DOMINANT_3S_TO_3P_E1_ONLY_NOT_FULL_RELATIVISTIC_LIFETIME',printed_air_coordinate_used=False)


def mchf_stop(text):
    markers=re.findall(r'RYA1354_ORBITAL_STOP\s+([\d.EeDd+\-]+)\s+([\d.EeDd+\-]+)',text)
    if len(markers)!=1:
        raise ValueError('Independent native orbital stop missing/ambiguous')
    dp,z2=(float(x.replace('D','E')) for x in markers[0])
    if not math.isfinite(dp) or not math.isfinite(z2) or not 0<=dp<=z2:
        raise ValueError('Independent orbital criterion not satisfied')
    return dict(weighted_change=dp,native_threshold=z2,orbital_stop=True)


def summarize(destination=None):
    base=ROOT/'data/audit/rya1354'
    files={};bundles={};external={}
    for name in ('lifetime_aux_v4b','lifetime_ladder_v4','lifetime_core_n5_v4','lifetime_relax_v4'):
        loaded,manifest=load_parts(base/name)
        for key,value in loaded.items():
            # The frozen plan is copied at each subsequent explicit extension.
            if key in files and files[key]!=value and key!='lifetime_plan.v4.json':
                raise ValueError('Conflicting source evidence across bundles: '+key)
            files[key]=value
        external.update(manifest.get('external_restart_caches',{}))
        bundles[name]=dict(manifest_sha256=hashlib.sha256((base/name/'manifest.json').read_bytes()).hexdigest(),parts=manifest['parts'])
    root='runs/lifetime_ladder_v3'
    receipt=json.loads(files[root+'/receipt.json'])
    if receipt['status']!='BALANCED_LIFETIME_LADDER_COMPLETED_REVIEW_REQUIRED' or len(receipt['qualified_states'])!=24:
        raise ValueError('Balanced lifetime sequence incomplete')
    levels=[];previous=None
    for n in (0,3,4,5):
        diagnostics={}
        for state in ('ground','quintet','triplet'):
            for verify in (False,True):
                path=f'{root}/n{n}_{state}'+('_verify' if verify else '')
                diagnostics[state+('_verify' if verify else '')]=scf_state(files,path,state)
                if verify:
                    original=f'{root}/n{n}_{state}/state.w'
                    if files[path+'/rwfn.inp']!=files[original]:
                        raise ValueError('Fresh restart did not use the declared state')
                elif n:
                    parent=0 if n==3 else n-1
                    if files[path+'/previous.w']!=files[f'{root}/n{parent}_{state}_verify/state.w']:
                        raise ValueError('Layer inherits different/unqualified earlier orbitals')
        require_triplet_quintet_order(files[f'{root}/n{n}_triplet_verify/state.lsj.lbl'].decode(),files[f'{root}/n{n}_quintet_verify/state.lsj.lbl'].decode())
        row=decay(files,f'{root}/n{n}_')
        row.update(n=n,SCF=diagnostics)
        if previous:
            row['change_from_previous_layer']=change(row,previous)
        levels.append(row);previous=row
    core={}
    for n in (3,4,5):
        r=f'runs/lifetime_core_n{n}_v1'
        if json.loads(files[r+'/receipt.json'])['status']!='BALANCED_CORE_LIFETIMES_COMPLETED_REVIEW_REQUIRED':
            raise ValueError('Balanced core lifetime cases incomplete')
        cases={}
        for case in ('DC','CV_2s_DC','CV_2s_Breit'):
            for state in ('ground','quintet','triplet'):
                identity(files[f'{r}/{case}_{state}/state.lsj.lbl'].decode(),state)
                seed=f'runs/lifetime_ladder_v3/n{n}_{state}_verify'
                if files[f'{r}/{case}_{state}/state.w']!=files[seed+'/state.w']:
                    raise ValueError('Fixed-orbital control changed the declared state orbitals')
                if files[f'{r}/{case}_{state}/isodata']!=files[seed+'/isodata']:
                    raise ValueError('Core control changed nuclear input identity')
                if case=='DC' and files[f'{r}/{case}_{state}/state.c']!=files[seed+'/state.c']:
                    raise ValueError('Same-space DC control changed the declared CSFs')
                if case=='CV_2s_Breit' and files[f'{r}/{case}_{state}/state.c']!=files[f'{r}/CV_2s_DC_{state}/state.c']:
                    raise ValueError('Core/Breit control changed more than Hamiltonian')
            require_triplet_quintet_order(files[f'{r}/{case}_triplet/state.lsj.lbl'].decode(),files[f'{r}/{case}_quintet/state.lsj.lbl'].decode())
            row=decay(files,f'{r}/{case}_',ci=True)
            row['change_from_VV_DC']=change(row,row if case=='DC' else cases['DC'])
            cases[case]=row
        core[str(n)]=cases
    build='toolchains/atsp-strict-scf-v2'
    build_receipt=json.loads(files[build+'/receipt.json'])
    if build_receipt['upstream_dirty']:
        raise ValueError('Independent variant did not start from clean pinned source')
    if strict_mchf_source(files[build+'/scf.original.f'].decode())!=files[build+'/scf.f'].decode():
        raise ValueError('Independent solver has undeclared source changes')
    control='validation/mchf_stop_control_v1'
    for name in ('wfn.out','1De.l'):
        if files[f'{control}/original/{name}']!=files[f'{control}/strict/{name}']:
            raise ValueError('Independent solver one-cycle software control differs')
    independent={}
    for version in (2,3,4):
        r=f'runs/lifetime_mchf_v{version}'
        for term in ('3P','3S'):
            path=f'{r}/reference_{term}/'
            validate_hf_seed(files[path+'HF.stdout'].decode()+files[path+'HF.stderr'].decode(),term)
        independent[str(version)]=[]
        for n in ((3,4) if version==2 else (3,4,5)):
            row=mchf_rate(files[f'{r}/n{n}_E1/lower.upper.ls'].decode())
            row.update(n=n,SCF={term:mchf_stop(files[f'{r}/n{n}_{term}/mchf.stdout'].decode()) for term in ('3P','3S')})
            if version==4:
                row['reference_SCF']={term:mchf_stop(files[f'{r}/qualified_reference_{term}/mchf.stdout'].decode()) for term in ('3P','3S')}
                grasp=next(item for item in levels if item['n']==n)
                row['relative_rate_difference_vs_GRASP_dominant_E1']={
                    'length_vs_B':row['partial_E1_A_length_s1']/grasp['dominant_3P_E1_A_s1']['B']-1,
                    'velocity_vs_C':row['partial_E1_A_velocity_s1']/grasp['dominant_3P_E1_A_s1']['C']-1}
            independent[str(version)].append(row)
    primary=json.loads((base/'source_curation_v2/lifetime_validation.json').read_text())['experiment']
    for row in [*levels,*(case for cases in core.values() for case in cases.values())]:
        row['primary_nonblind_residual_ns']={g:m['model_one_photon_lifetime_ns']-primary['value_ns'] for g,m in row['model'].items()}
    result=dict(schema='rya1354.balanced-lifetime-study.v4',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
                prediction_origin='THIS_UNPUBLISHED_THEORETICAL_CALCULATION',candidate_doi=None,bundles=bundles,
                balanced_layers=levels,core_controls=core,independent_partial_E1=independent,primary_experiment=primary,
                rejected_attempts=dict(coarse_reference_v1_v2='Native strict1s criterion exhausted500 cycles; finer consistent grid required, no layer promoted',
                                       mchf_v1='Wrong HF fixed-width closed-shell layout creates invalid s orbital and NaN/failed equations; partial seeds and downstream outputs rejected',
                                       mchf_v2_n5='No native orbital-stop marker after500 cycles; maximum function-change floor above requested tighter criterion. Preserve, do not score n5'),
                cache_policy='Native mcp.* and explicitly hashed rci.res restart matrix caches retained on Sirius. Offline verifies all archived result material and matching external-cache metadata; it does not reread external cache bytes.',
                external_restart_caches=external,
                limitations='VV frozen earlier layers, fixed-orbital2s expansion, 1s inactive in correlation. Larger MR/active-space and independent relativistic convergence remain owed. Independent LS rate covers dominant E1 only; ground3P term versus GRASP five-state EOL, native Rydberg/mass/grid conventions differ. No pure relativity isolation, complete independent lifetime or gauge-based uncertainty claim.',
                uncertainty=None,covariance=None,calibrated_sigma=None,chi_square=None,adoption_allowed=False)
    result['core_relaxation_attempts']={}
    for version in (1,2,3):
        root=f'runs/lifetime_relax_n3_v{version}'
        record=json.loads(files[root+'/receipt.json'])
        entry=dict(status=record['status'],error=record.get('error'),diagnostics=record['diagnostics'])
        if record['status']=='CORE_RELAXED_MODEL_LIFETIME_COMPLETED_REVIEW_REQUIRED':
            for state in ('ground','quintet','triplet'):
                for suffix in ('','_verify'):
                    scf_state(files,root+'/'+state+suffix,state)
            entry['decay']=decay(files,root+'/')
            entry['decay']['primary_nonblind_residual_ns']={g:m['model_one_photon_lifetime_ns']-primary['value_ns'] for g,m in entry['decay']['model'].items()}
            entry['change_from_fixed_orbital_core']=change(entry['decay'],core['3']['CV_2s_DC'])
        else:
            entry['lifetime_accepted']=False
        result['core_relaxation_attempts'][str(version)]=entry
    result['shared_odd_models']={}
    for version in (1,2):
        root=f'runs/lifetime_shared_odd_n3_v{version}'
        shared=json.loads(files[root+'/receipt.json'])
        entry=dict(status=shared['status'],error=shared.get('error'),diagnostics=shared['diagnostics'])
        if shared['status']=='SHARED_ODD_LIFETIME_COMPLETED_REVIEW_REQUIRED':
            for name in ('odd','odd_verify'):
                d=scf_diagnostics(files[f'{root}/{name}/rmcdhf.stdout'].decode(),files[f'{root}/{name}/rmcdhf.stderr'].decode())
                if d['iteration_limit_hit'] or not d['all_reported_functions_below_ACC']:
                    raise ValueError('Shared odd state not qualified')
            for state in ('ground','quintet','triplet'):
                identity(files[f'{root}/CI_{state}/state.lsj.lbl'].decode(),state)
            for state in ('quintet','triplet'):
                if files[f'{root}/CI_{state}/state.w']!=files[root+'/odd_verify/state.w']:
                    raise ValueError('Shared odd CI blocks do not use the same qualified orbitals')
            entry['decay']=decay(files,root+'/',ci=True)
            entry['decay']['primary_nonblind_residual_ns']={g:m['model_one_photon_lifetime_ns']-primary['value_ns'] for g,m in entry['decay']['model'].items()}
            entry['change_from_fixed_orbital_core']=change(entry['decay'],core['3']['CV_2s_DC'])
        else:
            entry['lifetime_accepted']=False
        result['shared_odd_models'][str(version)]=entry
    probe='runs/lifetime_relax_wave_probe_v1'
    first=files[probe+'/snapshot001.w'];second=files[probe+'/snapshot002.w']
    recorded=json.loads(files[probe+'/snapshot002.json'])
    replay=wave_difference(decode_wave(first),decode_wave(second))
    if replay!=recorded['versus_previous'] or hashlib.sha256(second).hexdigest()!=recorded['sha256']:
        raise ValueError('Saved-wave diagnostic does not replay')
    result['saved_wave_diagnostic']=dict(scope=recorded['scope'],timing=recorded['timing'],versus_previous=replay)
    result['generator_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (destination or base/'source_curation_v2/lifetime_summary_v4.json').write_text(json.dumps(result,indent=2)+'\n')
    print('balanced ns',[(r['n'],{g:m['model_one_photon_lifetime_ns'] for g,m in r['model'].items()}) for r in levels])
    print('core ns',{n:{case:{g:m['model_one_photon_lifetime_ns'] for g,m in r['model'].items()} for case,r in cases.items()} for n,cases in core.items()})


def plot_summary(source,destination):
    """Plot native model values; gauge pairs are not error bars."""
    import matplotlib.pyplot as plt
    d=json.loads(source.read_text())
    primary=d['primary_experiment'];value=primary['value_ns'];reported=primary['reported_uncertainty_ns']
    with plt.rc_context({'font.size':9,'axes.spines.top':False,'axes.spines.right':False}):
        fig,axes=plt.subplots(1,2,figsize=(11,4.6),sharey=True)
        for ax in axes:
            ax.axhspan(value-reported,value+reported,color='#e2e8f0',label='Reported experimental interval')
            ax.axhline(value,color='#64748b',linewidth=1)
            ax.grid(axis='y',alpha=.2)
        x=[3,4,5]
        vv={row['n']:row for row in d['balanced_layers']}
        colors={'B':'#1d4ed8','C':'#b45309'}
        for gauge in ('B','C'):
            field='model_one_photon_lifetime_ns'
            axes[0].plot(x,[vv[n]['model'][gauge][field] for n in x],'-o',color=colors[gauge],label=f'Valence SD, {gauge}')
            axes[0].plot(x,[d['core_controls'][str(n)]['CV_2s_DC']['model'][gauge][field] for n in x],'--s',color=colors[gauge],label=f'Plus 2s CI, {gauge}')
        axes[0].set(xticks=x,xlabel='Largest active principal quantum number',ylabel='Conditional one-photon lifetime (ns)',title='Balanced sequence and fixed-orbital core expansion')
        models=[('Valence',vv[3]),('2s CI',d['core_controls']['3']['CV_2s_DC'])]
        relaxed=d['core_relaxation_attempts']['3']
        shared=d['shared_odd_models']['2']
        if 'decay' in relaxed:models.append(('Relaxed\norbitals',relaxed['decay']))
        if 'decay' in shared:models.append(('Shared odd\norbitals',shared['decay']))
        for gauge in ('B','C'):
            axes[1].plot(range(len(models)),[row['model'][gauge]['model_one_photon_lifetime_ns'] for _,row in models],'o-',color=colors[gauge],label=gauge)
        axes[1].set(xticks=range(len(models)),xticklabels=[label for label,_ in models],title='Separate n=3 orbital-model controls')
        axes[0].legend(loc='upper right',fontsize=8,frameon=False)
        fig.suptitle('O I 3s triplet lifetime: exploratory, exposed nonblind models',fontsize=12)
        fig.text(.5,.025,'B/C are distinct gauges, not uncertainty bounds. Shading: Smith1971 reported interval, not a calibrated model acceptance gate.',ha='center',fontsize=8)
        fig.tight_layout(rect=(0,.07,1,.94))
        fig.savefig(destination,dpi=180,metadata={'Description':'Native results from '+source.name+'; scientific HOLD. Gauge agreement does not calibrate uncertainty.'})
        plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot',action='store_true')
    args=parser.parse_args()
    summarize()
    if args.plot:
        base=ROOT/'data/audit/rya1354'
        plot_summary(base/'source_curation_v2/lifetime_summary_v4.json',base/'lifetime_convergence_v4.png')

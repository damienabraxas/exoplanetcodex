"""Full two-stage LTE shape and inverse-input audit; no RT/SE or input changes."""
import argparse,csv,json,re,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint,load_plan_engine
    from pipeline.multi_atom import read_multi_atom,to_lightweaver
    from pipeline.oxygen_lte_contract import normalized_weights,native_lte_weights,held_lw_weights
    from pipeline.thermal_profile_contract import collision_summary
    from lightweaver import atomic_set,constants as c
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['plan','replay','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();plan=json.loads(a.plan.read_text());replay=json.loads((a.replay/'report.json').read_text());checks={};start=time.perf_counter()
    def pin(path,expected=None):
        path=Path(path);actual=fingerprint(path)
        if expected is not None and actual!=expected:raise ValueError('Changed bound input '+str(path))
        checks[str(path.resolve())]=actual
        return path
    if not replay['computational_completed']or plan['purpose']!='full-source-O41-HOP-fixed-population-rates':raise ValueError('Require completed full-source replay')
    native=require_codex_path(plan['native_analysis_resource']);run=require_codex_path(plan['native_run_resource'])
    pin(a.plan);pin(a.replay/'report.json');pin(a.replay/'rates.npz')
    for filename,key in [('arrays.npz','native_arrays'),('report.json','native_report')]:
        path=pin(native/filename,plan[key])
        if replay['inputs'].get(str(path))!=fingerprint(path):raise ValueError('Wrong replay native state')
    nr=json.loads((native/'report.json').read_text());_,app,_,engine=load_plan_engine(dict(engine='openmulti3d',profile_variant=nr['profile_variant']))
    if nr['engine_receipt']!=engine:raise ValueError('Changed qualified reference engine')
    run_report=pin(run/'report.json',nr['native_report']);rd=json.loads(run_report.read_text())
    for name in ['native/atoms/O.txt','native/params.nml']:
        entry=next(x for x in rd['output_files']if x['path']==name);pin(run/name,entry['receipt'])
    params=(run/'native/params.nml').read_bytes().replace(b'\x00',b'').decode()
    if re.search(r'\bDEBYE_SHIELD\s*=\s*T\s*,',params)is None or re.search(r'\bLEVEL_PARTF\s*=\s*F\s*,',params)is None:raise ValueError('Unqualified shielding/partition mode')
    source=pin(app/'HOP/const_mod.f90');pin(app/'HOP/composition_mod.f90');pin(app/'SRC/atom_mod.f90');constants={}
    for key in ['hh','bk','ee','em','ec','pi']:
        match=re.search(r'\b'+key+r'\s*=\s*([0-9.+\-dDeE]+)',source.read_text(),re.I)
        if match is None:raise ValueError('Missing source constant '+key)
        constants[key]=float(match[1].lower().replace('d','e'))
    table=[x.split()for x in (run/'native/atoms/O.txt').read_text().splitlines()[7:48]]
    if [int(x[0])for x in table]!=list(range(41)):raise ValueError('Invalid full level table')
    E=np.array([float(x[1])for x in table]);g=np.array([float(x[2])for x in table]);stage=np.array([int(x[3])for x in table])
    atom_file=pin(require_codex_path('experiments.bride_oxygen_atom'))
    if replay['inputs'].get(str(atom_file))!=fingerprint(atom_file):raise ValueError('Wrong replay source atom')
    atom=read_multi_atom(atom_file);model=to_lightweaver(atom)
    np.testing.assert_array_equal(g,[l.g for l in atom.levels]);np.testing.assert_array_equal(stage,[l.stage for l in atom.levels])
    with np.load(native/'arrays.npz')as d:f=d['consumed_state_cgs'][::-1].astype(float);star=d['native_nstar_cm3'][:,::-1].copy()
    with np.load(a.replay/'rates.npz')as d:
        np.testing.assert_array_equal(d['nstar'],star*1e6);C=d['replay_C_s_inv'].copy();pop=d['populations'].copy()
    T=f[:,2];ne=f[:,0];total=star.sum(axis=0)
    native_log,aux=native_lte_weights(E,g,stage,T,ne,constants)
    full_log,_=native_lte_weights(E,g,stage,T,ne,constants,round_native=False)
    no_shield_log,_=native_lte_weights(E,g,stage,T,ne,constants,shield=False)
    lwc=dict(h=c.HPlanck,k=c.KBoltzmann,m=c.MElectron,q=c.QElectron,epsilon=c.Epsilon0)
    lw_log=held_lw_weights([l.E_SI for l in model.levels],g,stage,T,ne*1e6,lwc)
    candidates={name:normalized_weights(log,total)for name,log in [('native_FP32',native_log),('native_FP64',full_log),('native_no_shield',no_shield_log),('held_LW',lw_log)]}
    actual_lw=atomic_set.lte_pops(model,T,ne*1e6,total,debye=True)
    lw_error=float(np.max(abs(candidates['held_LW']/actual_lw-1)))
    native_error=float(np.max(abs(candidates['native_FP32']/star-1)))
    # Arithmetic screens inherited from the prior same-stage LTE audit, declared
    # before this run; they are not physical uncertainty or admission thresholds.
    if native_error>1e-12 or lw_error>1e-12:raise ValueError('Full LTE formula reconstruction failed')
    ratio_error=float(np.max(abs(np.exp(native_log-native_log[0])/(star/star[0])-1)))
    row_data=[];statistics={}
    for name,values in candidates.items():
        stats=dict(max_relative_population_change=float(np.max(abs(values/star-1))),max_relative_ground_ratio_change=float(np.max(abs((values/values[0])/(star/star[0])-1))),collision_with_unchanged_C=collision_summary(C,values))
        # DIAGNOSTIC ONLY: preserve every upward coefficient; reproject reverse
        # ratios algebraically to exhibit the required ownership coupling.
        rebound=np.zeros_like(C)
        for j in range(1,41):
            for i in range(j):
                rebound[j,i]=C[j,i];rebound[i,j]=C[j,i]*values[i]/values[j]
        stats['algebraic_reverse_projection_only']=collision_summary(rebound,values)
        statistics[name]=stats
        for ion in [1,2]:
            fraction=values[stage==ion].sum(axis=0)/total;original=star[stage==ion].sum(axis=0)/total
            for z in range(len(T)):row_data.append(dict(profile=name,stage=ion,depth_top_to_bottom=z,T_K=T[z],native_fraction=original[z],candidate_fraction=fraction[z],relative_fraction_change=fraction[z]/original[z]-1))
    from pipeline.thermodynamic_canary import bb_thermal_defect,occupation
    active_inputs={}
    for name,values in candidates.items():
        centre=[];inverse_changes={'BB':[],'BF':[]}
        for kind,transitions in [('BB',atom.lines),('BF',atom.continua)]:
            for t in transitions:
                i=t.i-1;j=t.j-1
                change=(values[i]/values[j])/(star[i]/star[j])-1
                inverse_changes[kind].append(change)
                if kind=='BB':
                    beta=values[j]/values[i]*g[i]/g[j]
                    nu=(model.levels[j].E_SI-model.levels[i].E_SI)/c.HPlanck
                    centre.append(bb_thermal_defect(beta,occupation(c.HPlanck*nu/(c.KBoltzmann*T))))
        active_inputs[name]=dict(BB_transitions=len(inverse_changes['BB']),BF_transitions=len(inverse_changes['BF']),max_relative_inverse_LTE_factor_change={k:float(np.max(abs(v)))for k,v in inverse_changes.items()},max_BB_held_LW_line_centre_thermal_defect=float(np.max(abs(centre))),scope='Actual candidate populations and source graph; BB line-centre canary only, not finite-profile or actual RT. BF LTE-factor response requires joint inverse-opacity/emissivity/rate updates, not a standalone BF correctness claim.')
    pin(__file__)
    for path in [ROOT/'pipeline/oxygen_lte_contract.py',ROOT/'pipeline/thermal_profile_contract.py',ROOT/'pipeline/thermodynamic_canary.py',ROOT/'pipeline/multi_atom.py',Path(c.__file__),Path(atomic_set.__file__)]:pin(path)
    report=dict(generator=fingerprint(__file__),inputs=checks,computational_completed=True,scientific_admission='HOLD',levels=41,depths=len(T),stages=[1,2],all_stage_native_reconstruction_max_relative=native_error,all_ground_ratio_reconstruction_max_relative=ratio_error,held_LW_independent_formula_max_relative=lw_error,source_auxiliary_ranges={k:[float(v.min()),float(v.max())]for k,v in aux.items()},candidate_statistics=statistics,active_inverse_inputs=active_inputs,source_constants_cgs=constants,held_LW_constants_SI=lwc,normalization='Each diagnostic keeps the measured native total at each depth; complete finite-level shape reconstructed, absolute EOS/element normalization/partition function NOT independently validated. Actual run LEVEL_PARTF=F uses external partition functions; their common normalization cancels in level ratios.',scope='Read-only all-level Saha/Boltzmann/Debye/precision and collision inverse-input contract. Both native and held-LW formulas audited; no preferred profile, RT, SE, solver population/collision/opacity change or physical closure/admission.',seconds=time.perf_counter()-start)
    for path,expected in checks.items():
        if fingerprint(path)!=expected:raise ValueError('Changed input during full LTE audit')
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'ion_fractions.csv').open('w',newline='')as out:
        writer=csv.DictWriter(out,fieldnames=list(row_data[0]),lineterminator='\n');writer.writeheader();writer.writerows(row_data)
    np.savez_compressed(a.output/'diagnostic_populations.npz',native=star,populations_unchanged=pop,**candidates,**aux)
    # Preserve producing source snapshots, not source atomic input tables.
    for path in [Path(__file__),ROOT/'pipeline/oxygen_lte_contract.py']:(a.output/(path.stem+'_source.py')).write_bytes(path.read_bytes())
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k]for k in ['all_stage_native_reconstruction_max_relative','all_ground_ratio_reconstruction_max_relative','held_LW_independent_formula_max_relative','candidate_statistics','active_inverse_inputs']},indent=2))


if __name__=='__main__':main()

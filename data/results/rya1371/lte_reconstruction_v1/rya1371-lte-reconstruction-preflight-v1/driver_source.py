"""Reconstruct the held native LTE convention and audit collision reciprocity.

No RT, SE, coefficient, population, backend or physical constant changes.
"""
import argparse,csv,json,re,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint,load_plan_engine
    from config.constants import require_codex_path
    from pipeline.multi_atom import read_multi_atom,to_lightweaver
    from pipeline.thermodynamic_canary import collision_balance
    import lightweaver.constants as c
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['native-analysis','native-run','replay','output']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    _,app,_,engine=load_plan_engine(dict(engine='openmulti3d',profile_variant='accurate-kernel-control'))
    source=app/'HOP/const_mod.f90';atom_source=app/'SRC/atom_mod.f90'
    constants={}
    for name in ['hh','bk','cc','ee']:
        match=re.search(r'\b'+name+r'\s*=\s*([0-9.+\-dDeE]+)',source.read_text(),re.I)
        if match is None:raise ValueError('Missing source constant')
        constants[name]=float(match[1].lower().replace('d','e'))
    nr=json.loads((a.native_analysis/'report.json').read_text())
    native=json.loads((a.native_run/'report.json').read_text())
    replay=json.loads((a.replay/'report.json').read_text())
    if nr['engine_receipt']!=engine or nr['native_report']!=fingerprint(a.native_run/'report.json')or not replay['computational_completed']:
        raise ValueError('Unmatched native/replay provenance')
    checks={}
    for path in [source,atom_source,Path(c.__file__),Path(__file__),ROOT/'pipeline/thermodynamic_canary.py',ROOT/'pipeline/multi_atom.py']:
        checks[str(path)]=fingerprint(path)
    for directory,names in [(a.native_analysis,['report.json','arrays.npz']),(a.native_run,['report.json']),(a.replay,['report.json','rates.npz'])]:
        for name in names:
            path=directory/name;checks[str(path)]=fingerprint(path)
    for name in ['arrays.npz','report.json']:
        path=a.native_analysis/name
        if replay['inputs'].get(str(path))!=fingerprint(path):raise ValueError('Wrong replay state')
    native_atom=a.native_run/'native/atoms/O.txt'
    entry=next(x for x in native['output_files']if x['path']=='native/atoms/O.txt')
    if fingerprint(native_atom)!=entry['receipt']:raise ValueError('Changed native atom table')
    checks[str(native_atom)]=entry['receipt']
    lines=native_atom.read_text().splitlines()
    if list(map(int,lines[5].split()))[:2]!=[41,256]:raise ValueError('Unsupported native atom layout')
    table=[x.split()for x in lines[7:48]]
    if [int(x[0])for x in table]!=list(range(41)):raise ValueError('Invalid native level indices')
    ev=np.array([float(x[1])for x in table]);g=np.array([float(x[2])for x in table]);stage=np.array([int(x[3])for x in table])
    atom_path=require_codex_path('experiments.bride_oxygen_atom');checks[str(atom_path)]=fingerprint(atom_path)
    if replay['inputs'].get(str(atom_path))!=fingerprint(atom_path):raise ValueError('Wrong source atom')
    atom=read_multi_atom(atom_path);model=to_lightweaver(atom)
    np.testing.assert_array_equal(g,[x.g for x in atom.levels]);np.testing.assert_array_equal(stage,[x.stage for x in atom.levels])
    with np.load(a.native_analysis/'arrays.npz')as d:
        T=d['consumed_state_cgs'][::-1,2].astype(float);star=d['native_nstar_cm3'][:,::-1].copy()*1e6
    with np.load(a.replay/'rates.npz')as d:
        np.testing.assert_array_equal(d['nstar'],star)
        C=d['replay_C_s_inv'].copy();native_C=d['native_C_s_inv'].copy()
    np.testing.assert_array_equal(C,native_C)
    # atom_mod takes REAL(4) temperature and e_kT; LTE sums and ratios are
    # REAL(8). Same-stage ratios cancel common normalization and ionization.
    # This explicitly excludes a Saha/ionization-lowering reconstruction.
    e_k=constants['ee']/constants['bk']
    factor32=np.asarray(e_k/T.astype(np.float32).astype(float),dtype=np.float32).astype(float)
    rows=[];errors={name:[]for name in ['native_fp32','native_fp64','unrounded_HOP','source_LW']}
    for i in range(41):
        for j in range(i+1,41):
            if stage[i]!=stage[j]:continue
            actual=star[j]/star[i]*g[i]/g[j]
            gap_cm=atom.levels[j].E_cm-atom.levels[i].E_cm
            predictors=dict(native_fp32=np.exp(-(ev[j]-ev[i])*factor32),
                native_fp64=np.exp(-(ev[j]-ev[i])*e_k/T),
                unrounded_HOP=np.exp(-gap_cm*constants['hh']*constants['cc']/(constants['bk']*T)),
                source_LW=np.exp(-gap_cm*100*c.HPlanck*c.CLight/(c.KBoltzmann*T)))
            for name,predicted in predictors.items():errors[name].append(actual/predicted-1)
            for z in range(len(T)):
                row=dict(lower_1based=i+1,upper_1based=j+1,depth_top_to_bottom=z,stage=int(stage[i]),temperature_K=T[z],native_ratio=actual[z])
                row.update({name+'_relative_error':float(actual[z]/v[z]-1)for name,v in predictors.items()});rows.append(row)
    maxima={name:float(np.max(abs(np.asarray(values))))for name,values in errors.items()}
    if maxima['native_fp32']>1e-12:raise ValueError('Source LTE reconstruction does not close')
    defect,present=collision_balance(C,star)
    upper=np.triu(np.ones((41,41),dtype=bool),1)[:,:,None];selected=present&upper
    if not np.any(selected):raise ValueError('No collision reciprocity evidence')
    # A diagnostic population only: preserve each stage's total independently.
    # It cannot qualify Saha balance and is NEVER supplied to a solver.
    candidate=np.zeros_like(star)
    for ion in np.unique(stage):
        indices=np.flatnonzero(stage==ion);ground=indices[0]
        energy=np.array([model.levels[k].E_SI-model.levels[ground].E_SI for k in indices])
        weights=g[indices,None]*np.exp(-energy[:,None]/(c.KBoltzmann*T))
        candidate[indices]=weights/weights.sum(axis=0)*star[indices].sum(axis=0)
    counter,mask=collision_balance(C,candidate);np.testing.assert_array_equal(mask,present)
    report=dict(generator=fingerprint(__file__),inputs=checks,native_constants_cgs=constants,
        same_stage_level_pairs=len(rows)//len(T),depths=len(T),max_relative_LTE_ratio_errors=maxima,
        collision_pairs_at_any_depth=int(np.any(selected,axis=2).sum()),collision_pair_depths=int(selected.sum()),
        max_native_star_collision_pair_defect=float(np.max(abs(defect[selected]))),
        max_LW_within_stage_counterfactual_collision_pair_defect=float(np.max(abs(counter[selected]))),
        source_precision='Native e_k/T rounded to REAL(4), then REAL(8) level exponent and ratio; source table energies in eV at six decimals.',
        scientific_admission='HOLD',computational_completed=True,
        scope='Same-stage Boltzmann ratio reconstruction and exact fixed-C reciprocity. Counterfactual preserves each ion-stage total; no Saha, EOS, ionization-lowering or integrated radiative/SE validation. No solver inputs, constants, populations or coefficients changed. Collision inverse rates are constructed from native LTE; reciprocity is not independent validation of microscopic coefficients.')
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed audit input')
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'same_stage_ratios.csv').open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    np.savez_compressed(a.output/'collision_canaries.npz',native_defect=defect,present=present,counterfactual_defect=counter)
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items()if k.startswith('max_')or k.startswith('collision_')}))


if __name__=='__main__':main()

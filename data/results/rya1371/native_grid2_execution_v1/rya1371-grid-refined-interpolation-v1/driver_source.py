"""Frozen-field continuum interpolation counterfactuals using compiled native HOP."""
import argparse,csv,json,os,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint
    from pipeline.multi_atom import read_multi_atom,to_lightweaver
    from pipeline.adaptive_continuum_grid import nested_parabolic_measure
    import lightweaver.constants as c
    import lightweaver.atomic_model as am
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['build','output']:p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--refined-native-state',action='store_true')
    a=p.parse_args()
    if os.environ.get('OMP_NUM_THREADS')!='1'or os.environ.get('OPENBLAS_NUM_THREADS')!='1':raise ValueError('Require single-thread diagnostic')
    start=time.monotonic();checks={}
    def pin(path):
        path=Path(path);checks[str(path)]=fingerprint(path);return path
    pin(__file__);pin(c.__file__);pin(am.__file__);pin(ROOT/'pipeline/multi_atom.py');pin(ROOT/'pipeline/multi_continuum.py');pin(ROOT/'pipeline/adaptive_continuum_grid.py')
    build=json.loads(pin(a.build/'receipt.json').read_text());exe=pin(a.build/'continuum_probe')
    if not build['completed']or fingerprint(exe)!=build['executable']:raise ValueError('Unqualified oracle')
    for path,fp in build['inputs'].items():
        if fingerprint(path)!=fp:raise ValueError('Changed oracle source/object')
    checks.update(build['inputs'])
    if len(build.get('analytic_controls',[]))!=2:raise ValueError('Require analytic interpolation controls')
    from pipeline.reference_engines import load_engine
    _,source_app,_,_=load_engine('openmulti3d')
    for name in ['trans_mod.f90','interpol_mod.f90','const_mod.f90']:
        original=pin(source_app/'HOP'/name)
        linked=[Path(path)for path in build['inputs']if path.endswith('/HOP/'+name)]
        if len(linked)!=1 or fingerprint(linked[0])!=fingerprint(original):raise ValueError('Native interpolation source differs from linked oracle')
    pin(source_app/'SRC/atom_mod.f90');pin(source_app/'SRC/rt_patch_mod.f90')
    replay=require_codex_path('experiments.grid_refined_adaptive_replay'if a.refined_native_state else'experiments.bride_background_adaptive_v1')
    doc=json.loads(pin(replay/'report.json').read_text());pin(replay/'rates.npz')
    if not doc['computational_completed']:raise ValueError('Incomplete fixed-population producer')
    native=require_codex_path('experiments.openmulti3d_oxygen_grid_refined_analysis'if a.refined_native_state else'experiments.openmulti3d_oxygen_shared16_analysis')
    for name in ['arrays.npz','report.json']:
        path=pin(native/name)
        if doc['inputs'].get(str(path))!=fingerprint(path):raise ValueError('Mismatched native state')
    with np.load(native/'arrays.npz')as x:T=x['consumed_state_cgs'][::-1,2].astype(float)
    with np.load(replay/'rates.npz')as x:d={k:x[k].copy()for k in x.files}
    base=require_codex_path('experiments.grid_refined_background_grid'if a.refined_native_state else'experiments.bride_background_knot_grid')
    # Grid packet uses inherited source nodes; source-knot union is the dyadic base.
    with np.load(pin(base/'background_grid.npz'))as x:base_wave=x['wave_vac_nm'].copy()
    audit=require_codex_path('experiments.grid_refined_final_audit'if a.refined_native_state else'experiments.bride_background_final_audit_v1')
    ad=json.loads(pin(audit/'report.json').read_text())
    if not ad['computational_completed']or ad['inputs'].get(str(replay/'rates.npz'))!=fingerprint(replay/'rates.npz'):raise ValueError('Wrong parabolic audit')
    with np.load(pin(audit/'parabolic_rates.npz'))as x:baseline=x['R_s_inv'].copy();saved_bf=x['BF_s_inv'].copy()
    n=d['populations'];C=d['replay_C_s_inv'];gross=np.sum((baseline+C)*n[None,:,:],axis=1)+(baseline+C).sum(axis=0)*n
    atom_path=pin(require_codex_path('experiments.bride_oxygen_atom'));atom=read_multi_atom(atom_path);model=to_lightweaver(atom)
    consumed=pin(require_codex_path('experiments.openmulti3d_oxygen_grid_refined_run'if a.refined_native_state else'experiments.openmulti3d_oxygen_shared16_run')/'native/atoms/O.txt')
    lines=consumed.read_text().splitlines();begin=next(i for i,s in enumerate(lines)if s.split()==['i','lam_air','alpha','j'])+1
    tables={}
    for line in lines[begin:]:
        fields=line.split()
        if len(fields)!=4:break
        i,w,al,j=int(fields[0])+1,float(fields[1]),float(fields[2]),int(fields[3])+1
        tables.setdefault((i,j),[]).append((w,al))
    if set(tables)!={(t.i,t.j)for t in atom.continua}:raise ValueError('Consumed continuum graph mismatch')
    a.output.mkdir(parents=True,exist_ok=False)
    angular=float(np.sum(doc['context_audit']['angles']['weight']))
    matrices={k:baseline-saved_bf for k in ['native_interpolation_raw_table','native_interpolation_consumed_table']}
    rows=[];cross=[]
    for t,m in zip(atom.continua,model.continua):
        i,j=t.i-1,t.j-1
        wave=d['wave_vac_nm'];local=wave[(wave>=m.minLambda)&(wave<=m.lambdaEdge)]
        b=base_wave[(base_wave>=m.minLambda)&(base_wave<=m.lambdaEdge)]
        weight=nested_parabolic_measure(b,local)[:,None]
        idx=np.searchsorted(wave,local);J=d['Jnu'][idx];nu=c.CLight/(local[:,None]*1e-9)
        beta=np.exp(-c.HPlanck*nu/(c.KBoltzmann*T));bare=2*c.HPlanck*nu**3/c.CLight**2
        def integrate(alpha):
            pref=4*np.pi/c.HPlanck*alpha[:,None]/local[:,None]
            return np.sum(weight*pref*J,axis=0),d['nstar'][i]/d['nstar'][j]*np.sum(weight*pref*beta*(bare*angular+J),axis=0)
        own=m.alpha(local);up,down=integrate(own)
        np.testing.assert_allclose(up,saved_bf[j,i],rtol=5e-14,atol=0)
        np.testing.assert_allclose(down,saved_bf[i,j],rtol=5e-14,atol=0)
        raw=np.column_stack([t.wavelength[::-1],t.alpha[::-1]])
        native_table=np.array(tables[(t.i,t.j)])[::-1]
        for key,table in [('native_interpolation_raw_table',raw),('native_interpolation_consumed_table',native_table)]:
            request=f'{len(table)} {len(local)}\n'+'\n'.join(f'{w:.17e} {al:.17e}'for w,al in table)+'\n'+'\n'.join(f'{w*10:.17e}'for w in local)+'\n'
            name=f'{key}_{t.i}_{t.j}'
            (a.output/(name+'.in')).write_text(request)
            result=subprocess.run([str(exe.resolve())],input=request,text=True,capture_output=True,timeout=max(.1,600-(time.monotonic()-start)),check=False)
            (a.output/(name+'.out')).write_text(result.stdout);(a.output/(name+'.err')).write_text(result.stderr)
            if result.returncode:raise RuntimeError('Native interpolation oracle failed; output retained')
            alpha=np.fromstring(result.stdout,sep=' ')*1e-4
            if alpha.shape!=local.shape or not np.isfinite(alpha).all()or np.any(alpha<0):raise ValueError('Invalid native cross sections')
            u,v=integrate(alpha);matrices[key][j,i]+=u;matrices[key][i,j]+=v
            cross.append(dict(case=key,lower_1based=t.i,upper_1based=t.j,points=len(local),outside_consumed_support=int(np.count_nonzero((local*10<table[0,0])|(local*10>table[-1,0]))),max_alpha_relative_difference=float(np.max(abs(alpha/own-1)))))
            for depth in range(len(T)):
                delta=(u[depth]-up[depth])*n[i,depth]-(v[depth]-down[depth])*n[j,depth]
                rows.append(dict(case=key,lower_1based=t.i,upper_1based=t.j,depth_top_to_bottom=depth,up_relative_change=float(u[depth]/up[depth]-1),down_relative_change=float(v[depth]/down[depth]-1),signed_upper_balance_change=float(delta/gross[j,depth]),signed_lower_balance_change=float(-delta/gross[i,depth])))
    def net(R):return np.sum(R*n[None,:,:],axis=1)-R.sum(axis=0)*n
    case_rows=[]
    focus=(40,0)if a.refined_native_state else(39,3)
    for key,R in [('baseline',baseline),*matrices.items()]:
        r=net(R+C)/gross;peak=np.unravel_index(np.argmax(abs(r)),r.shape)
        case_rows.append(dict(case=key,max_absolute_balance_fixed_normalization=float(abs(r).max()),worst_level_1based=int(peak[0]+1),worst_depth=int(peak[1]),signed_balance_at_prior_peak_level40_depth3=float(r[39,3]),signed_balance_at_focus=float(r[focus]),max_absolute_balance_change=float(np.max(abs(net(R-baseline)/gross)))))
    for name,values in [('rates.csv',rows),('cross_sections.csv',cross)]:
        with(a.output/name).open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(values[0]),lineterminator='\n');w.writeheader();w.writerows(values)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed continuum audit input')
    peak_contributions=sorted([r for r in rows if r['depth_top_to_bottom']==focus[1] and r['upper_1based']==focus[0]+1],key=lambda r:-abs(r['signed_upper_balance_change']))
    for key,R in matrices.items():
        total=sum(r['signed_upper_balance_change']for r in peak_contributions if r['case']==key)
        np.testing.assert_allclose(total,net(R-baseline)[focus]/gross[focus],rtol=1e-9,atol=1e-14)
    report=dict(generator=fingerprint(__file__),inputs=checks,cases=case_rows,prior_peak_pair_contributions=peak_contributions,computational_completed=True,scientific_admission='HOLD',seconds=time.monotonic()-start,scope='Compiled native Hermite-frequency interpolation versus held LW WENO-wavelength alpha at same19963-node fixed J/n/nStar/C and positive parabolic measure. Raw-table case isolates interpolation including native FP32 alpha output; consumed-table case additionally includes native printed table rounding/support. Both evaluated only on original bridge continuum support: no claim of complete native-domain reconstruction. No fresh opacity/RT/SE feedback, no native J reconstruction, no physical adoption, no convergence/error-bound claim. All22continua/56depths; unchanged BB coefficients. Fixed baseline gross normalization retained for every case.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(case_rows))


if __name__=='__main__':main()

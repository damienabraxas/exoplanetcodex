"""DRY source-indicator Gauss4 radiation probe; no RT/SE or fitted node choice."""
import argparse,json,os,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path,codex_path
    from pipeline.reference_engines import fingerprint
    from pipeline.multi_atom import read_multi_atom
    from pipeline.native_continuum_control import load_native_continuum_control
    from pipeline.joint_oxygen_lte_control import validate_joint_plan
    from pipeline.radiation_probe_grid import gaussian_probe_grid,physics_plan,METHOD,registered_probe_destination
    from pipeline.saved_field_quadrature import interpolate_saved_field
    import lightweaver.constants as c
    import lightweaver.LwCompiled as compiled
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['fine-plan','coarse-plan','fine','coarse','output','resource']:p.add_argument('--'+key,required=True,type=str if key=='resource'else Path)
    a=p.parse_args();start=time.perf_counter();checks={};plans={};data={};docs={}
    def pin(path,expected=None):
        path=Path(path).resolve();fp=fingerprint(path)
        if expected is not None and fp!=expected:raise ValueError('Changed source '+str(path))
        checks[str(path)]=fp
        return path
    if os.environ.get('OMP_NUM_THREADS')!='1'or os.environ.get('OPENBLAS_NUM_THREADS')!='1':raise ValueError('Require one thread')
    registered_probe_destination(a.output,a.resource,resolver=codex_path)
    for name,plan_path,directory in [('fine',a.fine_plan,a.fine),('coarse',a.coarse_plan,a.coarse)]:
        plan=json.loads(pin(plan_path).read_text());validate_joint_plan(plan)
        if plan['purpose']!='full-source-O41-HOP-fixed-population-rates'or plan['seconds']!=600 or plan['backend']!=fingerprint(compiled.__file__):raise ValueError('Wrong source diagnostic/resources')
        doc=json.loads(pin(directory/'report.json').read_text());pin(directory/'rates.npz')
        if not doc['computational_completed']or doc['inputs'].get(str(plan_path.resolve()))!=fingerprint(plan_path):raise ValueError('Wrong source replay plan')
        pin(c.__file__,doc['inputs'][str(Path(c.__file__))]);plans[name]=plan;docs[name]=doc
        with np.load(directory/'rates.npz')as d:data[name]={k:d[k].copy()for k in d.files}
    if physics_plan(plans['fine'])!=physics_plan(plans['coarse']):raise ValueError('Physics differs between source indicator fields')
    fine=data['fine'];coarse=data['coarse']
    for key in ['populations','nstar','replay_C_s_inv','native_R_s_inv','native_C_s_inv']:np.testing.assert_array_equal(fine[key],coarse[key])
    idx=np.searchsorted(fine['wave_vac_nm'],coarse['wave_vac_nm']);np.testing.assert_array_equal(fine['wave_vac_nm'][idx],coarse['wave_vac_nm']);np.testing.assert_array_equal(fine['Jnu'][idx],coarse['Jnu'])
    native=require_codex_path(plans['fine']['native_analysis_resource']);pin(native/'arrays.npz',plans['fine']['native_arrays']);pin(native/'report.json',plans['fine']['native_report'])
    for doc in docs.values():
        for name in ['arrays.npz','report.json']:
            if doc['inputs'].get(str(native/name))!=fingerprint(native/name):raise ValueError('Wrong native state owner')
    with np.load(native/'arrays.npz')as d:T=d['consumed_state_cgs'][::-1,2].astype(float);native_star=d['native_nstar_cm3'][:,::-1].copy()*1e6
    atom_file=pin(require_codex_path('experiments.bride_oxygen_atom'),plans['fine']['continuum_interpolation_control']['atom']);atom=read_multi_atom(atom_file)
    for path in [Path(__file__),ROOT/'pipeline/radiation_probe_grid.py',ROOT/'pipeline/saved_field_quadrature.py',ROOT/'pipeline/native_continuum_control.py',ROOT/'pipeline/multi_atom.py',ROOT/'pipeline/joint_oxygen_lte_control.py',Path(compiled.__file__)]:pin(path)
    a.output.mkdir(parents=True,exist_ok=False)
    for path in [Path(__file__),ROOT/'pipeline/radiation_probe_grid.py']:(a.output/(path.stem+'_source.py')).write_bytes(path.read_bytes())
    control,more=load_native_continuum_control(plans['fine']['continuum_interpolation_control'],atom_file,atom,resolver=require_codex_path,output=a.output/'native_alpha_queries',deadline=start+600);checks.update(more)
    wave=fine['wave_vac_nm'];scores=np.zeros(len(wave)-1);rows=[]
    for t in atom.continua:
        pair=(t.i-1,t.j-1);table=control.tables[pair][:,0]*.1
        locations=np.flatnonzero((wave>=table[0])&(wave<=table[-1]));x=wave[locations]
        if len(x)<2 or x[0]!=table[0]or x[-1]!=table[-1]or np.any(np.diff(locations)!=1):raise ValueError('Wrong continuum support')
        field_difference=abs(fine['Jnu'][locations]-interpolate_saved_field(coarse['wave_vac_nm'],coarse['Jnu'],x))
        alpha=control.evaluate(pair,x);pref=(4*np.pi/c.HPlanck*alpha/x)[:,None]
        ratio=fine['nstar'][pair[0]]/fine['nstar'][pair[1]];beta=np.exp(-c.HPlanck*c.CLight/(x[:,None]*1e-9*c.KBoltzmann*T[None,:]))
        up=pref*field_difference;down=ratio*pref*beta*field_difference
        up_score=np.max(.5*np.diff(x)[:,None]*(up[1:]+up[:-1])/fine['replay_R_s_inv'][pair[1],pair[0]][None,:],axis=1)
        down_score=np.max(.5*np.diff(x)[:,None]*(down[1:]+down[:-1])/fine['replay_R_s_inv'][pair[0],pair[1]][None,:],axis=1)
        local=np.maximum(up_score,down_score);scores[locations[:-1]]=np.maximum(scores[locations[:-1]],local)
        rows.append(dict(lower_1based=t.i,upper_1based=t.j,max_relative_panel_indicator=float(local.max())))
    expanded,probe,weights,selected=gaussian_probe_grid(wave,scores,32)
    control.verify_inputs()
    for path,expected in checks.items():
        if fingerprint(path)!=expected:raise ValueError('Source changed during probe preparation')
    np.savez_compressed(a.output/'probe_grid.npz',wave_vac_nm=expanded,inherited_wave_vac_nm=wave,panel_scores=scores,probe_wave_vac_nm=probe,probe_weights_nm=weights,selected_panels=selected,populations=fine['populations'],native_nstar=native_star,declared_nstar=fine['nstar'])
    report=dict(generator=fingerprint(__file__),inputs=checks,computational_completed=True,scientific_admission='HOLD',method=METHOD,budget=32,output_resource=a.resource,physics_plan=physics_plan(plans['fine']),fine_rates=str((a.fine/'rates.npz').resolve()),coarse_rates=str((a.coarse/'rates.npz').resolve()),inherited_nodes=len(wave),expanded_nodes=len(expanded),panel_count=len(selected),selected_panels=selected.tolist(),selected_panel_indicators=scores[selected].tolist(),all_BF_indicators=rows,seconds=time.perf_counter()-start,scope='DRY all22BF/all56depth/up-down rate-normalized represented-field-curvature indicator. Select8highest positive fine intervals, stable wavelength ties,4interior Gaussian nodes each. Every ancestor retained, no native residual/observed fit or physical/RT/SE result.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    plan=dict(plans['fine']);plan['wavelength_grid_control']=dict(method=METHOD,resource=a.resource,packet=fingerprint(a.output/'probe_grid.npz'),producer=fingerprint(a.output/'report.json'));plan['scope']+=' Explicit32new-interior-radiation probe; all19963ancestor nodes retained,19995total under20000cap. Source-indicator selection only, no new adopted quadrature or coupled SE.'
    (a.output/'probe_plan.json').write_text(json.dumps(plan,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['method','inherited_nodes','expanded_nodes','panel_count','selected_panels','selected_panel_indicators','seconds']}))


if __name__=='__main__':main()

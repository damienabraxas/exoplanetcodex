"""Matched saved-J interpolant BF quadrature budget; no new RT/SE or operator."""
import argparse,csv,json,os,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint
    from pipeline.multi_atom import read_multi_atom
    from pipeline.native_continuum_control import load_native_continuum_control
    from pipeline.joint_oxygen_lte_control import validate_joint_plan
    from pipeline.bride_hop_transfer import native_lobatto5
    from pipeline.absolute_rate_integrals import bf_rates
    from pipeline.saved_field_quadrature import gaussian_panel_batches,interpolate_saved_field
    import lightweaver.constants as c
    import lightweaver.LwCompiled as compiled
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['plan','coarse-plan','fine','coarse','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();start=time.perf_counter();deadline=start+600;checks={};plans={};data={};docs={}
    if os.environ.get('OMP_NUM_THREADS')!='1'or os.environ.get('OPENBLAS_NUM_THREADS')!='1':raise ValueError('Require original single-thread diagnostic')
    def pin(path,expected=None):
        path=Path(path).resolve();fp=fingerprint(path)
        if expected is not None and fp!=expected:raise ValueError('Changed source '+str(path))
        checks[str(path)]=fp
        return path
    for name,plan_path,directory in [('fine',a.plan,a.fine),('coarse',a.coarse_plan,a.coarse)]:
        plan=json.loads(pin(plan_path).read_text());validate_joint_plan(plan)
        if plan['purpose']!='full-source-O41-HOP-fixed-population-rates'or plan['seconds']!=600 or plan['backend']!=fingerprint(compiled.__file__):raise ValueError('Wrong fixed-state backend/resources')
        doc=json.loads(pin(directory/'report.json').read_text());pin(directory/'rates.npz')
        if not doc['computational_completed']or doc['inputs'].get(str(plan_path.resolve()))!=fingerprint(plan_path)or not doc.get('joint_lte_control'):raise ValueError('Wrong/incomplete saved joint replay')
        plans[name]=plan;docs[name]=doc
        with np.load(directory/'rates.npz')as d:data[name]={key:d[key].copy()for key in d.files}
    if {k:v for k,v in plans['fine'].items()if k not in ['scope','wavelength_grid_control']}!={k:v for k,v in plans['coarse'].items()if k not in ['scope','wavelength_grid_control']}:raise ValueError('Unmatched physics plans')
    for key in ['populations','nstar','replay_C_s_inv','native_R_s_inv','native_C_s_inv']:np.testing.assert_array_equal(data['fine'][key],data['coarse'][key])
    fine=data['fine'];coarse=data['coarse'];idx=np.searchsorted(fine['wave_vac_nm'],coarse['wave_vac_nm']);np.testing.assert_array_equal(fine['wave_vac_nm'][idx],coarse['wave_vac_nm']);np.testing.assert_array_equal(fine['Jnu'][idx],coarse['Jnu'])
    native=require_codex_path(plans['fine']['native_analysis_resource']);pin(native/'arrays.npz',plans['fine']['native_arrays']);pin(native/'report.json',plans['fine']['native_report'])
    with np.load(native/'arrays.npz')as d:T=d['consumed_state_cgs'][::-1,2].astype(float)
    for doc in docs.values():
        if doc['inputs'].get(str(native/'arrays.npz'))!=fingerprint(native/'arrays.npz'):raise ValueError('Wrong consumed atmosphere')
    atom_path=pin(require_codex_path('experiments.bride_oxygen_atom'),plans['fine']['continuum_interpolation_control']['atom']);atom=read_multi_atom(atom_path)
    _,mu_weights,angles=native_lobatto5();W=float(mu_weights.sum())
    if any(doc['context_audit']['angles']!=angles for doc in docs.values()):raise ValueError('Changed actual angular convention')
    for path in [Path(__file__),ROOT/'pipeline/saved_field_quadrature.py',ROOT/'pipeline/native_continuum_control.py',ROOT/'pipeline/absolute_rate_integrals.py',ROOT/'pipeline/multi_atom.py',ROOT/'pipeline/bride_hop_transfer.py',ROOT/'pipeline/joint_oxygen_lte_control.py',Path(c.__file__),Path(compiled.__file__)]:pin(path)
    a.output.mkdir(parents=True,exist_ok=False)
    control,more=load_native_continuum_control(plans['fine']['continuum_interpolation_control'],atom_path,atom,resolver=require_codex_path,output=a.output/'native_alpha_queries',deadline=deadline);checks.update(more)
    manifest=dict(generator=fingerprint(__file__),inputs=checks,computational_completed=False,scientific_admission='HOLD',scope='Pure integration of declared saved physical J_nu interpolated linearly in vacuum nm; all22BF, orders4/8/16, no RT/SE/model change. Native raw-source frequency-Hermite FP32 alpha, held-LW thermal/photon factors and actual angular sum retained. Split all saved-J and source-alpha knots. No true-between-node radiation or physical-error certificate.')
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for path in [Path(__file__),ROOT/'pipeline/saved_field_quadrature.py']:(a.output/(path.stem+'_source.py')).write_bytes(path.read_bytes())
    matrices={name:{order:np.zeros_like(data[name]['replay_R_s_inv'])for order in [0,4,8,16]}for name in ['fine','coarse']};rows=[];max_nodes={name:0 for name in data};bf_mask=np.zeros_like(fine['replay_R_s_inv'],dtype=bool)
    for name in ['fine','coarse']:
        d=data[name];wave=d['wave_vac_nm'];star=d['nstar']
        for t in atom.continua:
            pair=(t.i-1,t.j-1);knots=control.tables[pair][:,0]*.1;select=(wave>=knots[0])&(wave<=knots[-1]);x=wave[select];J=d['Jnu'][select]
            if len(x)<2 or x[0]!=knots[0]or x[-1]!=knots[-1]:raise ValueError('Missing exact continuum support endpoints')
            max_nodes[name]=max(max_nodes[name],len(x));ratio=star[pair[0]]/star[pair[1]]
            alpha=control.evaluate(pair,x)
            up,down=bf_rates(x,alpha,J,nstar_ratio=ratio,T_K=T,h=c.HPlanck,k=c.KBoltzmann,c=c.CLight,angular_sum=W)
            matrices[name][0][pair[1],pair[0]]=up;matrices[name][0][pair[0],pair[1]]=down;bf_mask[pair[1],pair[0]]=True;bf_mask[pair[0],pair[1]]=True
            panels=np.unique(np.r_[x,knots])
            for order in [4,8,16]:
                upward=np.zeros(len(T));downward=upward.copy()
                for query,measure in gaussian_panel_batches(panels,order):
                    if time.perf_counter()>deadline:raise TimeoutError('Original600s saved-field quadrature cap exceeded')
                    al=control.evaluate(pair,query);field=interpolate_saved_field(x,J,query);nu=c.CLight/(query[:,None]*1e-9);bare=2*c.HPlanck*nu**3/c.CLight**2;beta=np.exp(-c.HPlanck*nu/(c.KBoltzmann*T[None,:]));weight=(4*np.pi/c.HPlanck*measure/query*al)[:,None]
                    upward+=np.sum(weight*field,axis=0);downward+=ratio*np.sum(weight*beta*(bare*W+field),axis=0)
                matrices[name][order][pair[1],pair[0]]=upward;matrices[name][order][pair[0],pair[1]]=downward
            for z in range(len(T)):
                row=dict(grid=name,lower_1based=t.i,upper_1based=t.j,depth_top_to_bottom=z,saved_nodes=len(x),integration_panels=len(panels)-1)
                for order in [4,8,16]:
                    row[f'up_G{order}_over_trapezoid_minus_one']=float(matrices[name][order][pair[1],pair[0],z]/up[z]-1);row[f'down_G{order}_over_trapezoid_minus_one']=float(matrices[name][order][pair[0],pair[1],z]/down[z]-1)
                rows.append(row)
    cases={}
    def residual(R,C,n):
        inc=np.sum((R+C)*n[None,:,:],axis=1);out=np.sum(R+C,axis=0)*n;local=abs(inc-out)/(inc+out);worst=np.unravel_index(np.argmax(local),local.shape)
        return dict(max_relative_local=float(local[worst]),level_1based=int(worst[0]+1),depth_top_to_bottom=int(worst[1]))
    for name in ['fine','coarse']:
        d=data[name];m=matrices[name];mask=bf_mask
        if any(not np.isfinite(v).all()or np.any(v[mask]<=0)for v in m.values()):raise ValueError('Incomplete/nonpositive quadrature rates')
        error=float(np.max(abs(m[0][mask]/d['replay_R_s_inv'][mask]-1)));ops=16*max_nodes[name]*len(mu_weights)*2;bound=ops*np.finfo(float).eps/(1-ops*np.finfo(float).eps)
        if error>bound:raise ValueError('Original BF trapezoid does not reconstruct compiled source rates')
        shadow={order:d['replay_R_s_inv']-m[0]+m[order]for order in [4,8,16]}
        if any(np.any(v<0)or not np.array_equal(v>0,d['replay_R_s_inv']>0)for v in shadow.values()):raise ValueError('Invalid whole-graph shadow matrix')
        cases[name]=dict(saved_nodes=len(d['wave_vac_nm']),BF_trapezoid_compiled_max_relative=error,arithmetic_screen=bound,max_relative_G8_vs_G4=float(np.max(abs(m[8][mask]/m[4][mask]-1))),max_relative_G16_vs_G8=float(np.max(abs(m[16][mask]/m[8][mask]-1))),orders={str(order):dict(max_relative_vs_trapezoid=float(np.max(abs(m[order][mask]/m[0][mask]-1))),whole_graph_shadow_residual=residual(shadow[order],d['replay_C_s_inv'],d['populations']))for order in [4,8,16]})
    # Exact telescoping budget of compiled coarse-to-fine BF differences.
    tc=matrices['coarse'][0];tf=matrices['fine'][0];gc=matrices['coarse'][16];gf=matrices['fine'][16]
    represented=gf-gc;defect_change=(tf-gf)-(tc-gc);observed=tf-tc
    scale=abs(tf)+abs(tc)+abs(gf)+abs(gc);mask=bf_mask
    closure=float(np.max(abs((represented+defect_change-observed)[mask])/scale[mask]))
    if closure>128*np.finfo(float).eps:raise ValueError('Telescoping budget arithmetic failed')
    budget=dict(max_relative_compiled_BF_refinement=float(np.max(abs(tf[mask]/tc[mask]-1))),max_relative_G16_represented_field_change=float(np.max(abs(gf[mask]/gc[mask]-1))),max_integration_defect_change_over_coarse_rate=float(np.max(abs(defect_change[mask]/tc[mask]))),telescoping_roundoff_normalized=closure,scope='Tfine-Tcoarse=(Gfine-Gcoarse)+[(Tfine-Gfine)-(Tcoarse-Gcoarse)]. Same source alpha and matched n/nStar/C; G integrates each declared saved-J interpolant. This separates represented-field and trapezoid-defect response, not unknown true radiation error. Maxima are not additive and not fit criteria.')
    control.verify_inputs()
    for path,expected in checks.items():
        if fingerprint(path)!=expected:raise ValueError('Changed quadrature audit input')
    report=dict(manifest,computational_completed=True,cases=cases,budget=budget,continua=len(atom.continua),depths=len(T),angular_weight_sum=W,max_oracle_query=max(row['inside_points']for row in control.calls),oracle_calls=len(control.calls),seconds=time.perf_counter()-start)
    with(a.output/'BF_quadrature.csv').open('w',newline='')as out:writer=csv.DictWriter(out,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    np.savez_compressed(a.output/'BF_quadrature.npz',**{f'{name}_G{order}':m for name,orders in matrices.items()for order,m in orders.items()},represented_field_change=represented,integration_defect_change=defect_change,observed_trapezoid_change=observed)
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='inputs'},indent=2))


if __name__=='__main__':main()

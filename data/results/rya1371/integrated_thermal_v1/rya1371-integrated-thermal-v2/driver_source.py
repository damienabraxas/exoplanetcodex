"""Quadrature-integrated thermal canaries using actual held compiled O profiles.

Independent integral audit only: no formal solve, radiation/SE iteration or
change to any context population, coefficient, constant or boundary.
"""
import argparse,csv,json,os,re,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint,load_plan_engine
    from pipeline.multi_atom import read_multi_atom
    from pipeline.multi_width_law import load_width_law
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.bride_hop_transfer import make_context
    from pipeline.thermodynamic_canary import occupation,integrated_defect,native_split_planck_defect
    import lightweaver.constants as c
    import lightweaver.LwCompiled as compiled
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['plan','replay','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();plan=json.loads(a.plan.read_text());replay=json.loads((a.replay/'report.json').read_text())
    if os.environ.get('OMP_NUM_THREADS')!='1'or os.environ.get('OPENBLAS_NUM_THREADS')!='1'or plan['seconds']!=600 or plan['backend']!=fingerprint(compiled.__file__)or not replay['computational_completed']:
        raise ValueError('Require unchanged bounded single-thread held backend')
    if plan['purpose']!='full-source-O41-HOP-fixed-population-rates'or plan['formal_solver']!='piecewise_linear_1d':raise ValueError('Require qualified source profile plan')
    checks={}
    def pin(path,expected=None):
        path=Path(path);fp=fingerprint(path)
        if expected is not None and fp!=expected:raise ValueError('Changed pinned input '+str(path))
        checks[str(path)]=fp
    pin(a.plan);pin(a.replay/'report.json');pin(a.replay/'rates.npz');pin(compiled.__file__);pin(c.__file__)
    controls=ROOT/plan['controls_source'];pin(controls,plan['controls_receipt'])
    if json.loads(controls.read_text())['LW_constants']!=fingerprint(c.__file__):raise ValueError('Changed inherited LW constants')
    native=require_codex_path(plan['native_analysis_resource'])
    for name,key in [('arrays.npz','native_arrays'),('report.json','native_report')]:
        path=native/name;pin(path,plan[key])
        if replay['inputs'].get(str(path))!=fingerprint(path):raise ValueError('Wrong matched replay state')
    doc=json.loads((native/'report.json').read_text());_,app,_,engine=load_plan_engine(dict(engine='openmulti3d',profile_variant=doc['profile_variant']))
    if doc['engine_receipt']!=engine:raise ValueError('Changed reference engine')
    source=app/'HOP/const_mod.f90';pin(source);pin(app/'SRC/atom_mod.f90');pin(app/'SRC/rt_patch_mod.f90');native_constants={}
    for name,scale in [('hh',1e-7),('bk',1e-7),('cc',.01)]:
        match=re.search(r'\b'+name+r'\s*=\s*([0-9.+\-dDeE]+)',source.read_text(),re.I)
        if match is None:raise ValueError('Missing native constant')
        native_constants[name]=float(match[1].lower().replace('d','e'))*scale
    if native_constants['cc']!=c.CLight:raise ValueError('Frequency coordinate mismatch')
    with np.load(native/'arrays.npz')as d:
        f=d['consumed_state_cgs'][::-1].astype(float);star=d['native_nstar_cm3'][:,::-1].copy()*1e6
        pop=d['native_populations_cm3'][:,::-1].copy()*1e6
        state=dict(height_m=d['height_m'][::-1].copy(),temperature_K=f[:,2],ne_m3=f[:,0]*1e6,rho_kg_m3=f[:,1]*1e3,vturb_m_s=f[:,7]*1e3)
    with np.load(a.replay/'rates.npz')as d:
        wave=d['wave_vac_nm'].copy();np.testing.assert_array_equal(star,d['nstar'])
    atom_path=require_codex_path('experiments.bride_oxygen_atom');ref_path=require_codex_path('experiments.bride_oxygen_reference');ref=json.loads(ref_path.read_text());model_path=require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name']
    pin(atom_path,ref['atom']);pin(ref_path);pin(model_path,ref['atmosphere_file'])
    atom=read_multi_atom(atom_path);law,more=load_width_law(plan['original_width_law'],atom_path,atom);checks.update(more)
    for name in ['bride_hop_transfer.py','hop_packet.py','hop_abi.py','multi_atom.py','multi_width_law.py','thermodynamic_canary.py']:
        pin(ROOT/'pipeline'/name)
    pin(__file__)
    start=time.perf_counter()
    # Detailed-static mode supplies exactly the source profiles without requiring
    # an active SE operator. No formal solver is called and no field is injected.
    ctx,audit=make_context(state,read_marcs_structure(model_path),atom,pop,star,extra_vac_nm=wave,formal_solver=plan['formal_solver'],helium_broadening=False,vdw_broadening=True,original_width_law=law)
    oxygen=ctx.detailedAtoms[0];np.testing.assert_array_equal(ctx.spect.wavelength,wave)
    rows=[];T=state['temperature_K'];mu_weights=np.asarray(audit['angles']['weight']);counts={'BB':0,'BF':0}
    for transition in oxygen.trans:
        if time.perf_counter()-start>plan['seconds']:raise TimeoutError('Original600s diagnostic cap exceeded')
        local=np.asarray(transition.wavelength);measure=np.empty_like(local)
        measure[0]=.5*(local[1]-local[0]);measure[-1]=.5*(local[-1]-local[-2]);measure[1:-1]=.5*(local[2:]-local[:-2])
        nu=c.CLight/(local[:,None]*1e-9);x=c.HPlanck*nu/(c.KBoltzmann*T[None,:]);n=occupation(x)
        hop_n=occupation(native_constants['hh']*nu/(native_constants['bk']*T[None,:]));i=transition.i;j=transition.j
        if transition.type=='Line':
            counts['BB']+=1;model=transition.transModel;nu0=c.CLight/model.lambda0_m
            phi=np.asarray(transition.phi)
            if phi.shape!=(len(local),len(mu_weights),2,len(T)):raise ValueError('Unexpected compiled profile shape')
            profile=np.einsum('lrsk,r->lk',phi,mu_weights)*.5
            weights=measure[:,None]*profile*(model.lambda0/local[:,None])
            beta=star[j]/star[i]*model.atom.levels[i].g/model.atom.levels[j].g
            x0=c.HPlanck*nu0/(c.KBoltzmann*T);consistent=np.exp(-x0)
            occupation_grid=(nu/nu0)**3*n
            hop_grid=(native_constants['hh']/c.HPlanck)*(nu/nu0)**3*hop_n
            centre=np.broadcast_to(occupation(x0),n.shape)
            cases=dict(centre_identity=integrated_defect(centre,consistent[None,:]*(1+centre),weights),
                LW_Planck_consistent_star=integrated_defect(occupation_grid,consistent[None,:]*(1+occupation_grid),weights),
                LW_Planck_native_star=integrated_defect(occupation_grid,beta[None,:]*(1+occupation_grid),weights),
                HOP_Planck_native_star=integrated_defect(hop_grid,beta[None,:]*(1+hop_grid),weights))
            # Source native rt_patch_mod/atom_mod integrates a thermal part and
            # a deviation from B. At J=B, dev=0 and its inverse relation is tied
            # to the imported LTE ratio. Common profile/Einstein factors cancel.
            thermal=np.sum(weights*hop_grid,axis=0)/(beta/(1-beta))
            cases['native_split_Planck_identity']=native_split_planck_defect(beta,thermal)
            cases['native_split_FP32_thermal_identity']=native_split_planck_defect(beta,thermal.astype(np.float32).astype(float))
            kind='BB'
        elif transition.type=='Continuum':
            counts['BF']+=1;bare=2*c.HPlanck*nu**3/c.CLight**2;J=bare*n;hop_J=2*native_constants['hh']*nu**3/c.CLight**2*hop_n
            weights=np.broadcast_to(measure[:,None]/local[:,None]*np.asarray(transition.alpha)[:,None],J.shape)
            beta=np.exp(-x)
            cases=dict(LW_Planck_consistent_star=integrated_defect(J,beta*(bare+J),weights),HOP_Planck_native_star=integrated_defect(hop_J,beta*(bare+hop_J),weights))
            kind='BF'
        else:raise ValueError('Unsupported thermal transition')
        for case,values in cases.items():
            for z,value in enumerate(values):rows.append(dict(kind=kind,lower_1based=i+1,upper_1based=j+1,case=case,depth_top_to_bottom=z,signed_integrated_defect=float(value)))
    if counts!={'BB':256,'BF':22}:raise ValueError('Incomplete source transition inventory')
    np.testing.assert_array_equal(oxygen.n,pop);np.testing.assert_array_equal(oxygen.nStar,star)
    maxima={}
    for row in rows:
        key=row['kind']+'/'+row['case'];maxima[key]=max(maxima.get(key,0),abs(row['signed_integrated_defect']))
    for key in ['BB/centre_identity','BF/LW_Planck_consistent_star','BB/native_split_Planck_identity','BB/native_split_FP32_thermal_identity']:
        if maxima[key]>1e-12:raise ValueError('Integrated thermal identity failed')
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed thermal audit input')
    report=dict(generator=fingerprint(__file__),inputs=checks,counts=counts,max_absolute_integrated_defects=maxima,context_audit=audit,seconds=time.perf_counter()-start,computational_completed=True,scientific_admission='HOLD',scope='Independent rate-integral ratios using actual compiled profiles/continuum cross sections and source trapezoid/angular measures. Common depth/transition factors cancel; not absolute rates or a direct compiled rate accumulator comparison. BB uses constant line-centre Einstein ratio (CRD) over its finite profile: a frequency-varying Planck field need not exactly balance it. BF uses per-frequency exponential; imported LTE factors cancel at n=nStar. No RT, SE, field/population/coefficient/constant change or scientific admission.')
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'integrated_defects.csv').open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(maxima))


if __name__=='__main__':main()

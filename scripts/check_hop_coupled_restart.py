"""Independent saved-transition R+C reconstruction and fresh-process restart check."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','restart','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();checks={};cases=[];populations=[]
    for root in [a.source,a.restart]:
        report=root/'report.json';r=json.loads(report.read_text())
        if not r['computational_completed']or r['dry_run']:raise ValueError('Unqualified completed coupled solve')
        for path,fp in r['inputs'].items():
            if fingerprint(path)!=fp:raise ValueError('Changed solver input '+path)
            checks[path]=fp
        for path in [report,root/'populations.npz',root/'rates.npz']:checks[str(path)]=fingerprint(path)
        with np.load(root/'populations.npz')as d:n=d['populations'].copy();star=d['nstar'].copy()
        with np.load(root/'rates.npz')as d:
            C=d['collision_offdiag'].copy();i=d['transition_i'];j=d['transition_j'];up=d['radiative_up_s_inv'];down=d['radiative_down_s_inv'];stored=d['physical_R_plus_C'].copy();stored_net=d['net'].copy()
        if C.shape!= (41,41,56)or len(i)!=278 or up.shape!=(278,56)or np.any(C[np.arange(41),np.arange(41)]!=0):raise ValueError('Invalid saved full-source R/C dimensions')
        off=C.copy()
        for lo,hi,u,dn in zip(i,j,up,down):off[hi,lo]+=u;off[lo,hi]+=dn
        incoming=np.sum(off*n[None,:,:],axis=1);outgoing=np.sum(off,axis=0)*n;net=incoming-outgoing;gross=incoming+outgoing
        physical=off.copy();physical[np.arange(41),np.arange(41)]=-np.sum(off,axis=0)
        np.testing.assert_allclose(physical,stored,rtol=256*np.finfo(float).eps,atol=0)
        disagreement=float(np.max(abs(net-stored_net)/np.maximum(gross,np.finfo(float).tiny)));res=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)));conservation=float(np.max(abs(n.sum(axis=0)/star.sum(axis=0)-1)))
        if disagreement>256*np.finfo(float).eps or res>=r['controls']['rtol']or conservation>1e-10:raise ValueError('Independent physical R+C/conservation check failed')
        cases.append(dict(source=str(root),independent_physical_residual=res,saved_net_disagreement_normalized=disagreement,conservation=conservation));populations.append((n,star,r))
    n,s,r=populations[0];rn,rs,rr=populations[1];np.testing.assert_array_equal(s,rs);np.testing.assert_array_equal(r['source_atom'],rr['source_atom']);np.testing.assert_array_equal(r['controls'],rr['controls'])
    restart_response=float(np.max(abs(rn/n-1)))
    if restart_response>=r['controls']['rtol']:raise ValueError('Fresh-process restart moves beyond original tolerance')
    from config.constants import require_codex_path
    from pipeline.multi_atom import read_multi_atom,to_lightweaver
    native_path=Path(next(p for p in r['inputs']if p.endswith('/arrays.npz')))
    atom_file=require_codex_path('experiments.bride_oxygen_atom');matrix_report=require_codex_path('experiments.native_oxygen_collision_matrix')/'report.json'
    if str(matrix_report)not in r['inputs']:raise ValueError('Wrong native source-constant report')
    table=Path(next(p for p in r['inputs']if p.endswith('/native/atoms/O.txt')))
    text=table.read_text().splitlines();h=next(i for i,x in enumerate(text)if x.split()==['ev','g','ion','label']);energy=np.array([float(x.split()[1])for x in text[h+1:h+42]]);g=np.array([float(x.split()[2])for x in text[h+1:h+42]]);ek=json.loads(matrix_report.read_text())['native_e_k'];model=to_lightweaver(read_multi_atom(atom_file))
    with np.load(native_path)as data:tau=data['source_logtau500_on_native_height'][::-1].copy();T=data['consumed_state_cgs'][::-1,2].astype(float);initial=data['native_populations_cm3'][:,::-1].copy()*1e6
    legacy=[]
    for line in model.lines:
        expected=g[line.j]/g[line.i]*np.exp(-(energy[line.j]-energy[line.i])*ek/T);legacy.append(float(np.max(abs((s[line.j]/s[line.i])/expected-1))))
    response=n/initial;where=np.unravel_index(np.argmin(response),response.shape);band=(tau>=-3)&(tau<=0)
    population_comparison=dict(min_native_population_ratio=float(response.min()),min_ratio_level_1based=int(where[0]+1),min_ratio_logtau500=float(tau[where[1]]),max_native_population_ratio=float(response.max()),triplet_max_relative_in_predeclared_logtau_minus3_to0=float(np.max(abs(response[[5,7,8,9]][:,band]-1))),native_table_Boltzmann_max_relative=float(max(legacy)),scope='Inherited logtau interval is a diagnostic band, not a measured contribution function. Large outer-layer differences remain unresolved; no layer exclusion from SE.')
    checks[str(atom_file)]=fingerprint(atom_file);checks[str(matrix_report)]=fingerprint(matrix_report);checks[str(table)]=fingerprint(table);checks[str(native_path)]=fingerprint(native_path);checks[str(ROOT/'pipeline/multi_atom.py')]=fingerprint(ROOT/'pipeline/multi_atom.py')
    a.output.mkdir(parents=True,exist_ok=False)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed independent analysis input')
    report=dict(generator=fingerprint(__file__),inputs=checks,cases=cases,population_comparison=population_comparison,restart_max_relative_population_change=restart_response,scientific_admission='HOLD',scope='Saved per-transition R plus C independently summed with incoming/outgoing flows. Numerical restart/rate closure only; shared source data, fixed historical-constant EOS/LTE and strong outer-layer differences remain unvalidated physics.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(cases=cases,restart_max_relative_population_change=restart_response),indent=2))

if __name__=='__main__':main()

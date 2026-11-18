"""Quantify installed constant/LTE convention mixing without changing any model."""
import argparse,csv,json,re,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint,load_plan_engine
    from pipeline.multi_atom import read_multi_atom,to_lightweaver
    from pipeline.thermodynamic_canary import occupation,bb_thermal_defect,bf_thermal_defect
    from config.constants import require_codex_path
    import lightweaver.constants as c
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--native-analysis',type=Path,required=True);p.add_argument('--replay',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    _,app,_,engine=load_plan_engine(dict(engine='openmulti3d',profile_variant='accurate-kernel-control'));source=app/'HOP/const_mod.f90';text=source.read_text();native={}
    for name,scale in [('hh',1e-7),('bk',1e-7),('cc',.01)]:
        m=re.search(r'\b'+name+r'\s*=\s*([0-9.+\-dDeE]+)',text,re.I)
        if m is None:raise ValueError('Missing source thermal constant')
        native[name]=float(m[1].lower().replace('d','e'))*scale
    if native['cc']!=c.CLight:raise ValueError('Frequency coordinate constants differ; separate control needed')
    doc=json.loads((a.native_analysis/'report.json').read_text());replay=json.loads((a.replay/'report.json').read_text())
    if doc['engine_receipt']!=engine or not replay['computational_completed']:raise ValueError('Unqualified matched state')
    checks={str(q):fingerprint(q)for q in [source,Path(c.__file__),a.native_analysis/'report.json',a.native_analysis/'arrays.npz',a.replay/'report.json',a.replay/'rates.npz',ROOT/'pipeline/thermodynamic_canary.py']}
    for name in ['arrays.npz','report.json']:
        q=a.native_analysis/name
        if replay['inputs'].get(str(q))!=fingerprint(q):raise ValueError('Thermal canary does not use the exact replay state')
    header=require_codex_path('engines.bride_eval')/'lw_src/Source/Constants.hpp'
    build=require_codex_path('experiments.bride_native_ghost_plugin')/'receipt.json';receipt=json.loads(build.read_text())
    if receipt['source_headers'].get('Constants.hpp')!=fingerprint(header):raise ValueError('Changed held LW constant header')
    for name,value in [('HPlanck',c.HPlanck),('KBoltzmann',c.KBoltzmann),('CLight',c.CLight)]:
        match=re.search(r'constexpr double '+name+r'\s*=\s*([0-9.+\-eE]+)',header.read_text())
        if match is None or float(match[1])!=value:raise ValueError('LW header/Python thermal convention mismatch')
    checks[str(header)]=fingerprint(header);checks[str(build)]=fingerprint(build)
    with np.load(a.native_analysis/'arrays.npz')as d:T=d['consumed_state_cgs'][::-1,2].astype(float);star=d['native_nstar_cm3'][:,::-1].copy()*1e6
    with np.load(a.replay/'rates.npz')as d:wave=d['wave_vac_nm'].copy();np.testing.assert_array_equal(d['nstar'],star)
    atom_file=require_codex_path('experiments.bride_oxygen_atom');atom=read_multi_atom(atom_file);model=to_lightweaver(atom);checks[str(atom_file)]=fingerprint(atom_file);rows=[];bb_consistent=[];bb_foreign=[];bb_native_field=[]
    for index,line in enumerate(model.lines):
        nu=c.CLight/line.lambda0_m;x=c.HPlanck*nu/(c.KBoltzmann*T);nb=occupation(x);native_nb=occupation(native['hh']*nu/(native['bk']*T));beta=star[line.j]*model.levels[line.i].g/(star[line.i]*model.levels[line.j].g)
        consistent=bb_thermal_defect(np.exp(-x),nb);foreign=bb_thermal_defect(beta,nb);mixed_field=bb_thermal_defect(beta,(native['hh']/c.HPlanck)*native_nb)
        bb_consistent.append(consistent);bb_foreign.append(foreign);bb_native_field.append(mixed_field)
        for z in range(len(T)):rows.append(dict(line_index=index,lower_1based=line.i+1,upper_1based=line.j+1,vacuum_nm=line.lambda0,depth_top_to_bottom=z,LW_blackbody_foreign_star_signed_defect=float(foreign[z]),HOP_blackbody_foreign_star_signed_defect=float(mixed_field[z])))
    nu=c.CLight/(wave[:,None]*1e-9);shape=(len(wave),len(T));x=c.HPlanck*nu/(c.KBoltzmann*T[None,:]);bare=np.broadcast_to(2*c.HPlanck*nu**3/c.CLight**2,shape);LW=bare*occupation(x);HOP=2*native['hh']*nu**3/native['cc']**2*occupation(native['hh']*nu/(native['bk']*T[None,:]));beta=np.exp(-x)
    bf_consistent=bf_thermal_defect(LW,bare,beta);bf_mixed=bf_thermal_defect(HOP,bare,beta)
    bb_consistent=np.asarray(bb_consistent);bb_foreign=np.asarray(bb_foreign);bb_native_field=np.asarray(bb_native_field)
    if np.max(abs(bb_consistent))>64*np.finfo(float).eps or np.max(abs(bf_consistent))>64*np.finfo(float).eps:raise ValueError('Analytic consistent-constant identity fails')
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'bb_line_depth.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    np.savez_compressed(a.output/'canaries.npz',wave_vac_nm=wave,temperature_K=T,BB_foreign_star=np.array(bb_foreign),BB_native_blackbody=np.array(bb_native_field),BF_mixed=bf_mixed)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed thermodynamic audit input')
    report=dict(generator=fingerprint(__file__),inputs=checks,native_constants_SI=native,LW_constants_SI=dict(h=c.HPlanck,k=c.KBoltzmann,c=c.CLight),relative_h_over_k_difference=float((native['hh']/native['bk'])/(c.HPlanck/c.KBoltzmann)-1),max_BB_consistent_defect=float(np.max(abs(bb_consistent))),max_BB_foreign_star_LW_blackbody_defect=float(np.max(abs(bb_foreign))),max_BB_foreign_star_HOP_blackbody_defect=float(np.max(abs(bb_native_field))),max_BF_consistent_defect=float(np.max(abs(bf_consistent))),max_BF_mixed_HOP_blackbody_LW_coefficients_defect=float(np.max(abs(bf_mixed))),max_relative_blackbody_difference=float(np.max(abs(HOP/LW-1))),scientific_admission='HOLD',computational_completed=True,scope='Analytic BB line-centre/BF per-frequency LTE identities at existing state. No RT, SE, field/constant/star/rate modification, no modern-constant substitution or fitted parameter. Separate installed constant sets; not full profile/frequency-integrated rates, microscopic collision validation, or unique attribution of the NLTE residual. Any coherent convention change requires joint LTE/background/collision/atomic requalification.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in report.items()if k.startswith('max_')or k=='relative_h_over_k_difference'}))
if __name__=='__main__':main()

"""Sample actual LW damping formulas on the matched native consumed atmospheric state."""
import argparse,json,os,sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.multi_atom import read_multi_atom,to_lightweaver
    from pipeline.hop_packet import compute_packet
    from pipeline.bride_hop_transfer import Perturbers,NativePerturberBroadener,verify_helium_reads
    from pipeline.bride_native_table_control import apply_table_control
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());checks={str(a.plan.resolve()):fingerprint(a.plan)}
    import lightweaver.broadening as broadening
    import lightweaver.atomic_model as atomic_model
    import lightweaver.atomic_table as atomic_table
    import lightweaver.constants as constants
    import lightweaver.LwCompiled as compiled
    if fingerprint(Path(compiled.__file__))!=plan['backend']:raise ValueError('Changed qualified LW backend')
    for module in [broadening,atomic_model,atomic_table,constants,compiled]:
        path=Path(module.__file__);checks[str(path)]=fingerprint(path)
    if os.environ.get('OMP_NUM_THREADS')!='1'or plan.get('helium_broadening')is not False or plan.get('table_control',{}).get('mode')!='centres_and_strengths':raise ValueError('Require declared matched optical table/He-off plan')
    state_file=require_codex_path(plan['native_analysis_resource'])/'arrays.npz'
    if fingerprint(state_file)!=plan['native_arrays']:raise ValueError('Changed native state')
    checks[str(state_file)]=fingerprint(state_file)
    with np.load(state_file)as d:state=d['consumed_state_cgs'][::-1].astype(float)
    atom_file=require_codex_path('experiments.bride_oxygen_atom');ref_file=require_codex_path('experiments.bride_oxygen_reference');ref=json.loads(ref_file.read_text());model_file=require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name'];table=require_codex_path(plan['table_control']['resource'])/'native/atoms/O.txt'
    if fingerprint(atom_file)!=ref['atom']or fingerprint(model_file)!=ref['atmosphere_file']or fingerprint(table)!=plan['table_control']['receipt']:raise ValueError('Changed physical/table input')
    for path in [atom_file,ref_file,model_file,table]:checks[str(path)]=fingerprint(path)
    for name in ['pipeline/bride_native_table_control.py','pipeline/bride_hop_transfer.py','pipeline/hop_packet.py','pipeline/multi_atom.py']:
        path=ROOT/name;checks[str(path)]=fingerprint(path)
    structure=read_marcs_structure(model_file);atom=read_multi_atom(atom_file);packet,packet_meta=compute_packet(state[:,2],state[:,0]*1e6,state[:,1]*1e3,np.array([500.,845.]),structure['abundances_dex'],atom,vmic_m_s=float(state[0,7]*1e3))
    atmos=SimpleNamespace(temperature=state[:,2],ne=state[:,0]*1e6,vturb=state[:,7]*1e3);model=to_lightweaver(atom);audit=apply_table_control(model,table,'centres_and_strengths');perturbers=Perturbers(packet,helium=False);arrays=[]
    for record in audit['changes']:
        line=model.lines[record['line_index']];line.broadening.elastic=[NativePerturberBroadener(x,perturbers)for x in line.broadening.elastic];line.setup(model);damping,elastic=line.damping(atmos,perturbers)
        if damping.shape!=(56,)or not np.isfinite(damping).all()or np.any(damping<=0):raise ValueError('Invalid actual damping support')
        arrays.append(damping)
    verify_helium_reads(perturbers.reads,False);samples=np.asarray(arrays);a.output.mkdir(parents=True,exist_ok=False);np.savez_compressed(a.output/'damping.npz',a_damp=samples,v_broad_m_s=model.vBroad(atmos),line_indices=np.array([x['line_index']for x in audit['changes']]),temperature_K=atmos.temperature)
    for path,expected in checks.items():
        if fingerprint(path)!=expected:raise ValueError('Changed damping input')
    report=dict(generator=fingerprint(__file__),inputs=checks,packet=packet_meta,table_control=audit,shape=list(samples.shape),damping_range=[float(samples.min()),float(samples.max())],output=fingerprint(a.output/'damping.npz'),scientific_admission='HOLD',limitations=['Actual LW broadening formulas at matched56 native consumed states with declared table and He-off control; not native stored damping values.','No RT or SE solve, no populations or source changed; native versus LW broadening coefficients remain to compare.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(shape=report['shape'],damping_range=report['damping_range'])))
if __name__=='__main__':main()

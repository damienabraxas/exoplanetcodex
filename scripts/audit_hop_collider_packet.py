"""Compare a read-only native H-ground export against the original HOP packet."""
import argparse,json,os,sys
from pathlib import Path
import ctypes as ct
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.hop_packet import compute_packet,checked_runtime
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());native=require_codex_path(plan['native_analysis_resource'])/'arrays.npz'
    if fingerprint(native)!=plan['native_arrays']or os.environ.get('OMP_NUM_THREADS')!='1':raise ValueError('Wrong native input/thread state')
    atom_file=require_codex_path('experiments.bride_oxygen_atom');reference=require_codex_path('experiments.bride_oxygen_reference');ref=json.loads(reference.read_text());model_file=require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name']
    if fingerprint(atom_file)!=ref['atom']or fingerprint(model_file)!=ref['atmosphere_file']:raise ValueError('Changed source atom/atmosphere')
    checks={str(q):fingerprint(q)for q in [a.plan.resolve(),native,atom_file,reference,model_file,ROOT/'pipeline/hop_packet.py',ROOT/'pipeline/hop_abi.py']}
    with np.load(native)as d:state=d['consumed_state_cgs'][::-1].astype(float)
    atom=read_multi_atom(atom_file);structure=read_marcs_structure(model_file)
    args=(state[:,2],state[:,0]*1e6,state[:,1]*1e3,np.array([500.,777.,845.]),structure['abundances_dex'],atom)
    original,om=compute_packet(*args,vmic_m_s=state[0,7]*1e3)
    collider,cm=compute_packet(*args,vmic_m_s=state[0,7]*1e3,provider_variant='collider')
    for key,value in original.items():np.testing.assert_array_equal(collider[key],value)
    ground=collider['H_ground_m3'];total=collider['H_I_m3']
    if np.any(ground<=0)or np.any(collider['H_ground_occupation_probability']*ground>total)or np.any(collider['H_II_m3']<=0)or np.any(collider['H_-_m3']<=0):raise ValueError('Invalid positive physical collider fields')
    # Fresh library is closed by compute_packet. Guard checks must leave sentinels intact.
    api,_,runtime,_,_,_=checked_runtime('collider');checks.update(runtime)
    refusal=[]
    for n in [-1,0,len(ground)]:
        out=np.full(len(ground),-123.,dtype=float);occupation=out.copy();status=api['hop_get_h_ground'](n,out.ctypes.data_as(ct.POINTER(ct.c_double)),occupation.ctypes.data_as(ct.POINTER(ct.c_double)))
        if status==0:raise ValueError('Unloaded/invalid ground export accepted')
        np.testing.assert_array_equal(out,-123.);np.testing.assert_array_equal(occupation,-123.);refusal.append(dict(n=n,status=status))
    # A second collider packet proves rejected calls did not contaminate the next EOS.
    replay,rm=compute_packet(*args,vmic_m_s=state[0,7]*1e3,provider_variant='collider')
    for key,value in collider.items():np.testing.assert_array_equal(replay[key],value)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed collider audit source')
    a.output.mkdir(parents=True,exist_ok=False);np.savez_compressed(a.output/'packet.npz',**collider)
    report=dict(generator=fingerprint(__file__),inputs=checks,original_packet=om,collider_packet=cm,replay_packet=rm,
        common_fields_bit_identical=list(original),refusals=refusal,ground_over_neutral_total_range=[float((ground/total).min()),float((ground/total).max())],
        occupation_weighted_ground_over_HI_range=[float((ground*collider['H_ground_occupation_probability']/total).min()),float((ground*collider['H_ground_occupation_probability']/total).max())],scientific_admission='HOLD',
        scope='Read-only native ground-H/proton/H-minus packet; no fake H-level populations, rate update, coupled SE, closed molecular budget, or solar admission.')
    report['packet']=fingerprint(a.output/'packet.npz');(a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k]for k in ['ground_over_neutral_total_range','occupation_weighted_ground_over_HI_range','refusals']}))

if __name__=='__main__':main()

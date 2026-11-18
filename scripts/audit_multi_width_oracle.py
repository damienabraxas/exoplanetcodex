"""Compare exact original MULTI widths with actual LW prescriptions at matched native states."""
import argparse,json,os,subprocess,sys,time
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
    from pipeline.bride_hop_transfer import Perturbers,NativePerturberBroadener
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--oracle',type=Path,required=True);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());receipt=a.oracle/'receipt.json';build=json.loads(receipt.read_text());exe=a.oracle/'multi_width_probe'
    if os.environ.get('OMP_NUM_THREADS')!='1'or not build['completed']or fingerprint(exe)!=build['executable']:raise ValueError('Require verified bounded native oracle')
    checks={str(receipt):fingerprint(receipt),str(exe):fingerprint(exe),str(a.plan.resolve()):fingerprint(a.plan)}
    for name,expected in build['inputs'].items():
        if fingerprint(name)!=expected:raise ValueError('Changed oracle source')
        checks[name]=expected
    native=require_codex_path(plan['native_analysis_resource'])/'arrays.npz'
    if fingerprint(native)!=plan['native_arrays']:raise ValueError('Changed matched atmospheric state')
    checks[str(native)]=fingerprint(native)
    with np.load(native)as d:state=d['consumed_state_cgs'][::-1].astype(float)
    atom_file=require_codex_path('experiments.bride_oxygen_atom');ref_file=require_codex_path('experiments.bride_oxygen_reference');reference=json.loads(ref_file.read_text());model_file=require_codex_path('grids.marcs_standard')/reference['atmosphere_file']['name']
    if fingerprint(atom_file)!=reference['atom']or fingerprint(model_file)!=reference['atmosphere_file']:raise ValueError('Changed held atom/atmosphere')
    for f in [atom_file,ref_file,model_file]:checks[str(f)]=fingerprint(f)
    import lightweaver.LwCompiled as compiled
    import lightweaver.broadening as broadening
    import lightweaver.atomic_model as atomic_model
    import lightweaver.constants as constants
    if fingerprint(Path(compiled.__file__))!=plan['backend']:raise ValueError('Changed LW backend')
    for m in [compiled,broadening,atomic_model,constants]:checks[str(Path(m.__file__))]=fingerprint(m.__file__)
    for rel in ['pipeline/multi_atom.py','pipeline/hop_packet.py','pipeline/bride_hop_transfer.py']:
        f=ROOT/rel;checks[str(f)]=fingerprint(f)
    atom=read_multi_atom(atom_file)
    if atom.symbol!='O'or len(atom.levels)!=41 or len(atom.lines)!=256 or len(state)!=56:raise ValueError('Unqualified atom/state scope')
    structure=read_marcs_structure(model_file);packet,packet_receipt=compute_packet(state[:,2],state[:,0]*1e6,state[:,1]*1e3,np.array([500.,845.]),structure['abundances_dex'],atom,vmic_m_s=float(state[0,7]*1e3));a.output.mkdir(parents=True,exist_ok=False);start=time.perf_counter()
    rows=[f'41 256 56 {atom.weight:.17e}']
    rows += [f'{l.E_cm:.17e} {l.g:.17e} {l.stage:d}'for l in atom.levels]
    rows += [f'{x.i:d} {x.j:d} {x.gRad:.17e} {x.gVdw:.17e} {x.gStark:.17e} {1e8/(atom.levels[x.j-1].E_cm-atom.levels[x.i-1].E_cm):.17e}'for x in atom.lines]
    rows += [f'{state[k,2]:.17e} {state[k,0]:.17e} {packet["H_I_m3"][k]*1e-6:.17e} {state[k,7]*1e3:.17e}'for k in range(56)]
    (a.output/'input.txt').write_text('\n'.join(rows)+'\n');manifest=dict(generator=fingerprint(__file__),inputs=checks,packet=packet_receipt,scientific_admission='HOLD',scope='Original exact DAMP with aggregate neutral-H storage only; actual source O41,56 matched states, all256 lines; no RT/SE or physical width default adoption.')
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    with(a.output/'input.txt').open()as fi,(a.output/'output.txt').open('w')as fo,(a.output/'run.log').open('w')as log:subprocess.run([str(exe.resolve())],cwd=a.output,stdin=fi,stdout=fo,stderr=log,timeout=600,check=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
    data=np.loadtxt(a.output/'output.txt',skiprows=1)
    if data.shape!=(256*56,7)or not np.isfinite(data).all():raise ValueError('Incomplete/nonfinite native widths')
    np.testing.assert_array_equal(data[:,0],np.repeat(np.arange(1,257),56));np.testing.assert_array_equal(data[:,1],np.tile(np.arange(1,57),256));native_widths=data[:,2:].reshape(256,56,5)
    model=to_lightweaver(atom);atmos=SimpleNamespace(temperature=state[:,2],ne=state[:,0]*1e6,vturb=state[:,7]*1e3);comparisons={};vbroads=[]
    for helium in [False,True]:
        current=to_lightweaver(atom);perturbers=Perturbers(packet,helium=helium);gamma=[];damping=[]
        for line in current.lines:
            line.broadening.elastic=[NativePerturberBroadener(x,perturbers)for x in line.broadening.elastic];line.setup(current);width=line.broadening.broaden(atmos,perturbers);gamma.append(width.Qelast);damping.append(line.damping(atmos,perturbers)[0])
        comparisons['LW_He_on'if helium else'LW_He_off']=dict(gamma=np.asarray(gamma),damping=np.asarray(damping));vbroads.append(current.vBroad(atmos))
    np.testing.assert_array_equal(native_widths[:,:,1],np.broadcast_to([x.gRad for x in atom.lines],(56,256)).T);np.testing.assert_array_equal(native_widths[:,:,3],0);np.testing.assert_allclose(native_widths[:,:,0],native_widths[:,:,1]+native_widths[:,:,2],rtol=1e-14)
    stats={}
    for name,d in comparisons.items():
        groups={}
        for family,mask in [('Unsold',np.array([x.gVdw==1 for x in atom.lines])),('ABO',np.array([x.gVdw>=20 for x in atom.lines]))]:
            ratio=d['gamma'][mask]/native_widths[mask,:,2];groups[family]=dict(n_lines=int(mask.sum()),LW_over_original_min=float(ratio.min()),LW_over_original_max=float(ratio.max()))
        stats[name]=groups
    constline=(a.output/'output.txt').read_text().splitlines()[0];consts=dict(zip(['BK_cgs','UU_g','HH_cgs','CC_cgs','EE_erg_per_eV','PI'],map(float,constline.split()[1:])))
    np.savez_compressed(a.output/'widths.npz',native_widths=native_widths,LW_He_off_gamma=comparisons['LW_He_off']['gamma'],LW_He_on_gamma=comparisons['LW_He_on']['gamma'],LW_He_off_adamp=comparisons['LW_He_off']['damping'],LW_He_on_adamp=comparisons['LW_He_on']['damping'],temperature_K=state[:,2],native_HI_m3=packet['H_I_m3'],native_HeI_m3=packet['He_I_m3'],source_gvdw=np.array([x.gVdw for x in atom.lines]))
    for path,f in checks.items():
        if fingerprint(path)!=f:raise ValueError('Changed width probe input')
    report=dict(**manifest,seconds=time.perf_counter()-start,native_constants=consts,source_atomic_mass_amu=atom.weight,LW_atomic_mass_amu=float(model.element.mass),statistics=stats,output=fingerprint(a.output/'widths.npz'),limitations=['Original source comments name fixed He/H0.1 in legacy Unsold recipe; runtime total-H input is not proof of helium-free physical coefficients.','Modern LW and legacy prescriptions are compared without rescaling or selecting by agreement; no physical width uncertainty/adoption inferred solely from equality.','Native H aggregate slot is solely for O DAMP reduction, not a collider-level H population packet. Fresh balance/EOS/collider/solar/3D/all-band gates stay open.']);(a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(stats))
if __name__=='__main__':main()

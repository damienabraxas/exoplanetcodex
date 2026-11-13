"""Original DAMP component-law controls with frozen EOS densities, not new stellar models."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--base',type=Path,required=True);p.add_argument('--oracle',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if os.environ.get('OMP_NUM_THREADS')!='1':raise ValueError('Require one thread')
    build_file=a.oracle/'receipt.json';build=json.loads(build_file.read_text());exe=a.oracle/'multi_width_probe';base_report=a.base/'report.json';base=json.loads(base_report.read_text());base_arrays=a.base/'widths.npz';base_input=a.base/'input.txt'
    if not build['completed']or fingerprint(exe)!=build['executable']or fingerprint(base_arrays)!=base['output']:raise ValueError('Unqualified base/oracle')
    for source,f in base['inputs'].items():
        if fingerprint(source)!=f:raise ValueError('Changed baseline source')
    inputs={str(x):fingerprint(x)for x in [exe,build_file,base_report,base_arrays,base_input]};lines=base_input.read_text().splitlines();nk,nl,nd=map(int,lines[0].split()[:3]);mass=float(lines[0].split()[3]);offset=1+nk+nl;state=np.array([[float(v)for v in line.split()]for line in lines[offset:]])
    if (nk,nl,nd)!=(41,256,56)or state.shape!=(56,4):raise ValueError('Unexpected baseline support')
    with np.load(base_arrays)as d:original=d['native_widths'].copy();gvdw=d['source_gvdw'].copy()
    cases={'half_neutral_H':dict(h_scale=.5),'plus_10K_frozen_EOS':dict(delta_T=10.),'xi_0p25kms':dict(vmic_m_s=250.),'xi_2kms':dict(vmic_m_s=2000.),'xi_4kms':dict(vmic_m_s=4000.)};a.output.mkdir(parents=True,exist_ok=False);start=time.perf_counter();results=[];arrays={};c=base['native_constants'];lam=np.array([float(x.split()[-1])for x in lines[1+nk:offset]])
    def speed(s):return np.sqrt(2*(c['BK_cgs']*1e-7)*s[:,0]/(mass*c['UU_g']*1e-3)+s[:,3]**2)
    unit_a=original[:,:,0]*lam[:,None]*1e-10/(4*c['PI']*speed(state)[None,:]);unit_error=float(np.max(abs(original[:,:,4]/unit_a-1)));bound=float(8*np.finfo(np.float32).eps)
    if unit_error>bound:raise ValueError('Native DAMP normalization/unit mapping fails source-literal precision bound')
    manifest=dict(generator=fingerprint(__file__),inputs=inputs,cases=cases,seconds_limit=600,threads=1,scientific_admission='HOLD',interpretation='Frozen-EOS width-component laws only: no new physical atmosphere, RT or SE solve.');(a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for name,control in cases.items():
        target=a.output/name;target.mkdir();s=state.copy();s[:,0]+=control.get('delta_T',0);s[:,2]*=control.get('h_scale',1)
        if 'vmic_m_s'in control:s[:,3]=control['vmic_m_s']
        payload='\n'.join(lines[:offset])+'\n'+'\n'.join(' '.join(f'{v:.17e}'for v in row)for row in s)+'\n';(target/'input.txt').write_text(payload)
        with(target/'input.txt').open()as fi,(target/'output.txt').open('w')as fo,(target/'run.log').open('w')as log:subprocess.run([str(exe.resolve())],stdin=fi,stdout=fo,stderr=log,cwd=target,timeout=600-(time.perf_counter()-start),check=True)
        raw=np.loadtxt(target/'output.txt',skiprows=1)
        if raw.shape!=(256*56,7)or not np.isfinite(raw).all():raise ValueError('Invalid law output')
        current=raw[:,2:].reshape(256,56,5);expected=original[:,:,2].copy()*control.get('h_scale',1)
        if control.get('delta_T',0):
            power=np.where(gvdw>=20,(1-(gvdw-np.floor(gvdw)))/2,.3);expected*=((s[:,0]/state[:,0])[None,:])**power[:,None]
        relative=float(np.max(abs(current[:,:,2]/expected-1)))
        if relative>bound:raise ValueError('Native width functional law failed')
        np.testing.assert_array_equal(current[:,:,1],original[:,:,1]);np.testing.assert_array_equal(current[:,:,3],original[:,:,3])
        if 'vmic_m_s'in control:np.testing.assert_array_equal(current[:,:,:4],original[:,:,:4])
        expected_a=current[:,:,0]*lam[:,None]*1e-10/(4*c['PI']*speed(s)[None,:]);a_error=float(np.max(abs(current[:,:,4]/expected_a-1)))
        if a_error>bound:raise ValueError('Native Doppler/microturbulence units failed')
        arrays[name]=current;results.append(dict(case=name,max_vdw_law_relative_error=relative,max_adamp_SI_relative_error=a_error,vdw_ratio_min=float(np.min(current[:,:,2]/original[:,:,2])),vdw_ratio_max=float(np.max(current[:,:,2]/original[:,:,2]))))
    np.savez_compressed(a.output/'laws.npz',**arrays)
    for path,f in inputs.items():
        if fingerprint(path)!=f:raise ValueError('Changed source during laws')
    report=dict(**manifest,seconds=time.perf_counter()-start,baseline_adamp_SI_relative_error=unit_error,source_literal_roundoff_bound=bound,results=results,output=fingerprint(a.output/'laws.npz'),limits='Component-law controls do not qualify atmospheric T/xi parameter inference or collisional H ground-state semantics; no physical/default adoption.');(a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(baseline_adamp_SI_error=unit_error,cases=results)))
if __name__=='__main__':main()

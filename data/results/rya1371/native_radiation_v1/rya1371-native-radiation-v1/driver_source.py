"""Read saved native J and reconstruct continuum rates before field comparisons."""
import argparse,csv,json,os,re,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint,load_engine
    from pipeline.openmulti3d_output import read_patch_set
    from pipeline.multi_atom import read_multi_atom
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if os.environ.get('OMP_NUM_THREADS')!='1'or os.environ.get('OPENBLAS_NUM_THREADS')!='1':raise ValueError('Require one-thread audit')
    checks={};start=time.monotonic()
    def pin(path):
        path=Path(path);checks[str(path)]=fingerprint(path);return path
    for path in [Path(__file__),ROOT/'pipeline/openmulti3d_output.py',ROOT/'pipeline/multi_atom.py']:pin(path)
    run=require_codex_path('experiments.openmulti3d_oxygen_shared16_run');doc=json.loads(pin(run/'report.json').read_text())
    if not doc['computational_completed']:raise ValueError('Incomplete native producer')
    for row in doc['output_files']:
        if fingerprint(pin(run/row['path']))!=row['receipt']:raise ValueError('Changed native producer output')
    _,app,_,engine=load_engine('openmulti3d')
    # Source arithmetic is pinned explicitly; this is not a new formal solve.
    for path in [app/'HOP/const_mod.f90',app/'HOP/spectrum_mod.f90',app/'SRC/rt_patch_mod.f90',app/'SRC/atom_mod.f90',app.parents[1]/'tasks/patch_mod.f90']:pin(path)
    text=(app/'HOP/const_mod.f90').read_text()
    constants={key:float(re.search(r'\b'+key+r'\s*=\s*([0-9.+\-dDeE]+)',text,re.I)[1].lower().replace('d','e'))for key in ['hh','cc','bk']}
    h,cc,k=constants['hh'],constants['cc'],constants['bk']
    build=json.loads(pin(a.build/'receipt.json').read_text());exe=pin(a.build/'continuum_probe')
    if not build['completed']or len(build.get('analytic_controls',[]))!=2 or fingerprint(exe)!=build['executable']:raise ValueError('Unqualified compiled interpolation oracle')
    for path,fp in build['inputs'].items():
        if fingerprint(path)!=fp:raise ValueError('Changed oracle input')
    checks.update(build['inputs'])
    save=run/'native/save';wave=np.fromfile(pin(save/'all_lam.bin'),dtype='<f8');pin(save/'Jnu_meta.txt')
    patches=read_patch_set(save,'Jnu');sizes=[patches[(1,i)].shape[0]for i in [1,2]]
    if sizes!=[19921,19920]or len(wave)!=sum(sizes)or np.any(np.diff(wave)<=0):raise ValueError('Unexpected native window/grid layout')
    if set(patches)!={(z,w)for z in range(1,9)for w in [1,2]}:raise ValueError('Incomplete native radiation mesh')
    windows=[]
    for w in [1,2]:
        parts=[patches[(z,w)]for z in range(1,9)]
        if any(x.shape!=(sizes[w-1],1,1,7,1)for x in parts):raise ValueError('Unexpected native depth/set shape')
        windows.append(np.concatenate([x[:,0,0,:,0]for x in parts],axis=1)[:,::-1])
    J=np.concatenate(windows,axis=0)
    # The strict reader checks both saved twin hemispheres are byte/value equal.
    analysis=require_codex_path('experiments.openmulti3d_oxygen_shared16_analysis');pin(analysis/'arrays.npz')
    with np.load(analysis/'arrays.npz')as x:T=x['consumed_state_cgs'][::-1,2].astype(np.float32)
    replay=require_codex_path('experiments.bride_background_adaptive_v1');pin(replay/'report.json')
    with np.load(pin(replay/'rates.npz'))as x:d={key:x[key].copy()for key in x.files}
    table=pin(run/'native/atoms/O.txt');lines=table.read_text().splitlines();begin=next(i for i,s in enumerate(lines)if s.split()==['i','lam_air','alpha','j'])+1;tables={}
    for line in lines[begin:]:
        f=line.split()
        if len(f)!=4:break
        tables.setdefault((int(f[0])+1,int(f[3])+1),[]).append((float(f[1]),float(f[2])))
    atom=read_multi_atom(pin(require_codex_path('experiments.bride_oxygen_atom')))
    if set(tables)!={(t.i,t.j)for t in atom.continua}:raise ValueError('Wrong continuum table graph')
    a.output.mkdir(parents=True,exist_ok=False);rows=[]
    for t in atom.continua:
        tab=np.array(tables[(t.i,t.j)])[::-1];idx=np.flatnonzero((wave>=tab[0,0])&(wave<=tab[-1,0]));local=wave[idx];nu=cc/(local*1e-8)
        if len(local)<2 or np.any(np.diff(idx)!=1):raise ValueError('Incomplete native continuum support')
        alpha=[]
        for number,chunk in enumerate(np.array_split(local,2)):
            if len(chunk)<2 or len(chunk)>20000:raise ValueError('Oracle query outside existing packet size scope')
            request=f'{len(tab)} {len(chunk)}\n'+'\n'.join(f'{w:.17e} {v:.17e}'for w,v in tab)+'\n'+'\n'.join(f'{w:.17e}'for w in chunk)+'\n'
            name=f'continuum_{t.i}_{t.j}_{number}';(a.output/(name+'.in')).write_text(request)
            remaining=600-(time.monotonic()-start)
            if remaining<=0:raise TimeoutError('Original diagnostic budget exceeded')
            result=subprocess.run([str(exe.resolve())],input=request,text=True,capture_output=True,timeout=remaining,check=False)
            (a.output/(name+'.out')).write_text(result.stdout);(a.output/(name+'.err')).write_text(result.stderr)
            if result.returncode:raise RuntimeError('Native interpolation caller refused')
            value=np.fromstring(result.stdout,sep=' ',dtype=np.float32)
            if value.shape!=chunk.shape or not np.isfinite(value).all()or np.any(value<0):raise ValueError('Invalid native alpha')
            alpha.append(value)
        alpha=np.concatenate(alpha);weight=np.empty_like(nu);weight[0]=abs(nu[1]-nu[0])*.5;weight[-1]=abs(nu[-1]-nu[-2])*.5;weight[1:-1]=abs(nu[2:]-nu[:-2])*.5
        # Match explicit native storage casts, half-hemisphere accumulation and
        # window-local FP32 sums; no tolerance or rate correction is applied.
        bare=(2*h/cc**2*nu**3).astype(np.float32);wmunu=(.5*weight/(h/(4*np.pi)*nu)).astype(np.float32)
        wa=(bare*alpha*wmunu).astype(np.float32)
        hkt=np.float32((h/k)/T);hv=(nu[:,None]*hkt).astype(np.float32);beta=np.maximum(np.exp(-hv),np.finfo(np.float32).tiny)
        native_J=J[idx];sums=[]
        for field in [native_J]:
            up=np.zeros(56);down=np.zeros(56)
            for win in [idx<sizes[0],idx>=sizes[0]]:
                if not np.any(win):continue
                absorption=wa[win,None]*field[win]
                emission=(wa[win,None]*beta[win])*(field[win]+np.float32(1))
                up+=2*np.cumsum(absorption,axis=0,dtype=np.float32)[-1].astype(float)
                down+=2*np.cumsum(emission,axis=0,dtype=np.float32)[-1].astype(float)
            down*=d['nstar'][t.i-1]/d['nstar'][t.j-1];sums.append((up,down))
        up,down=sums[0];native_up=d['native_R_s_inv'][t.j-1,t.i-1];native_down=d['native_R_s_inv'][t.i-1,t.j-1]
        for depth in range(56):rows.append(dict(lower_1based=t.i,upper_1based=t.j,depth_top_to_bottom=depth,up_relative_reconstruction=float(up[depth]/native_up[depth]-1),down_relative_reconstruction=float(down[depth]/native_down[depth]-1)))
    with(a.output/'reconstruction.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    maximum=max(abs(r[key])for r in rows for key in ['up_relative_reconstruction','down_relative_reconstruction'])
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed native radiation audit input')
    report=dict(generator=fingerprint(__file__),inputs=checks,computational_completed=True,scientific_admission='HOLD',max_relative_BF_rate_reconstruction=maximum,continua=len(atom.continua),depths=56,native_window_points=sizes,seconds=time.monotonic()-start,scope='Read-only saved native J reconstruction gate only. Native J is in dimensionless photon-occupation work units; thermal/spontaneous terms and two combined twin hemispheres included. No new RT/SE, interpolation onto bridge grid, field adoption or numerical/physical PASS. Stage consistency must be assessed from full rate reconstruction; a successful file read is insufficient.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='inputs'}))


if __name__=='__main__':main()

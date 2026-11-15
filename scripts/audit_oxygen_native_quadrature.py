"""Stream Gaussian BF integrals of a declared saved-J interpolant; no RT/SE."""
import argparse,csv,json,os,re,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint,load_engine
    from pipeline.openmulti3d_output import read_patch_set
    from pipeline.multi_atom import read_multi_atom
    from pipeline.hop_packet import MAX_PACKET_WAVELENGTHS
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['build','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();start=time.monotonic();checks={}
    if os.environ.get('OMP_NUM_THREADS')!='1'or os.environ.get('OPENBLAS_NUM_THREADS')!='1':raise ValueError('Require one-thread diagnostic')
    def pin(path):
        path=Path(path);checks[str(path)]=fingerprint(path);return path
    pin(__file__)
    for name in ['openmulti3d_output.py','multi_atom.py','hop_packet.py']:pin(ROOT/'pipeline'/name)
    run=require_codex_path('experiments.openmulti3d_oxygen_shared16_run');doc=json.loads(pin(run/'report.json').read_text())
    if not doc['computational_completed']:raise ValueError('Incomplete native producer')
    for r in doc['output_files']:
        if fingerprint(pin(run/r['path']))!=r['receipt']:raise ValueError('Changed native output')
    build=json.loads(pin(a.build/'receipt.json').read_text());exe=pin(a.build/'continuum_probe')
    if not build['completed']or len(build.get('analytic_controls',[]))!=2 or fingerprint(exe)!=build['executable']:raise ValueError('Unqualified oracle')
    for path,fp in build['inputs'].items():
        if fingerprint(path)!=fp:raise ValueError('Changed oracle source/object')
    checks.update(build['inputs'])
    _,app,_,_=load_engine('openmulti3d');text=pin(app/'HOP/const_mod.f90').read_text()
    constants={key:float(re.search(r'\b'+key+r'\s*=\s*([0-9.+\-dDeE]+)',text,re.I)[1].lower().replace('d','e'))for key in ['hh','cc','bk']};h,cc,k=[constants[key]for key in ['hh','cc','bk']]
    save=run/'native/save';wave=np.fromfile(pin(save/'all_lam.bin'),dtype='<f8');pin(save/'Jnu_meta.txt');parts=read_patch_set(save,'Jnu')
    if set(parts)!={(z,w)for z in range(1,9)for w in [1,2]}:raise ValueError('Wrong native layout')
    J=np.concatenate([np.concatenate([parts[(z,w)][:,0,0,:,0]for z in range(1,9)],axis=1)[:,::-1]for w in [1,2]],axis=0).astype(float)
    if J.shape!=(len(wave),56):raise ValueError('Wrong native J shape')
    nu=cc/(wave*1e-8);physical=J*(2*h/cc**2*nu**3).astype(np.float32)[:,None]
    native=require_codex_path('experiments.openmulti3d_oxygen_shared16_analysis')
    with np.load(pin(native/'arrays.npz'))as x:T=x['consumed_state_cgs'][::-1,2].astype(float)
    own=require_codex_path('experiments.bride_background_adaptive_v1')
    with np.load(pin(own/'rates.npz'))as x:d={key:x[key].copy()for key in x.files}
    lines=pin(run/'native/atoms/O.txt').read_text().splitlines();begin=next(i for i,s in enumerate(lines)if s.split()==['i','lam_air','alpha','j'])+1;tables={}
    for line in lines[begin:]:
        f=line.split()
        if len(f)!=4:break
        tables.setdefault((int(f[0])+1,int(f[3])+1),[]).append((float(f[1]),float(f[2])))
    atom=read_multi_atom(pin(require_codex_path('experiments.bride_oxygen_atom')))
    if set(tables)!={(t.i,t.j)for t in atom.continua}:raise ValueError('Wrong BF graph')
    a.output.mkdir(parents=True,exist_ok=False);counter=0;max_query=0
    def alpha_at(table,query):
        nonlocal counter,max_query
        if not 2<=len(query)<=MAX_PACKET_WAVELENGTHS:raise ValueError('Oracle request exceeds unchanged packet scope')
        max_query=max(max_query,len(query));remaining=600-(time.monotonic()-start)
        if remaining<=0:raise TimeoutError('Original600s diagnostic cap exceeded')
        request=f'{len(table)} {len(query)}\n'+'\n'.join(f'{x:.17e} {y:.17e}'for x,y in table)+'\n'+'\n'.join(f'{x:.17e}'for x in query)+'\n'
        name=f'oracle_{counter:05d}';counter+=1;(a.output/(name+'.in')).write_text(request)
        r=subprocess.run([str(exe.resolve())],input=request,text=True,capture_output=True,timeout=remaining,check=False)
        (a.output/(name+'.out')).write_text(r.stdout);(a.output/(name+'.err')).write_text(r.stderr)
        if r.returncode:raise RuntimeError('Oracle failed; raw evidence retained')
        out=np.fromstring(r.stdout,sep=' ')
        if out.shape!=query.shape or not np.isfinite(out).all()or np.any(out<0):raise ValueError('Invalid native alpha')
        return out
    matrices={order:np.zeros_like(d['native_R_s_inv'])for order in [0,4,8,16]};rows=[]
    for t in atom.continua:
        table=np.array(tables[(t.i,t.j)])[::-1];lo,hi=table[[0,-1],0]
        sel=(wave>=lo)&(wave<=hi);x=wave[sel];field=physical[sel]
        if len(x)<2 or x[0]!=lo or x[-1]!=hi:raise ValueError('Require exact native/source support endpoints')
        ratio=d['nstar'][t.i-1]/d['nstar'][t.j-1]
        def integrands(q,al,jn):
            freq=cc/(q*1e-8);pref=4*np.pi/h*al[:,None]/freq[:,None];beta=np.exp(-h*freq[:,None]/(k*T));bare=2*h/cc**2*freq[:,None]**3
            return pref*jn,ratio*pref*beta*(bare+jn)
        old_alpha=np.concatenate([alpha_at(table,z)for z in np.array_split(x,int(np.ceil(len(x)/MAX_PACKET_WAVELENGTHS)))])
        u,v=integrands(x,old_alpha,field);freq=cc/(x*1e-8);weight=np.empty_like(freq);weight[0]=abs(freq[1]-freq[0])*.5;weight[-1]=abs(freq[-1]-freq[-2])*.5;weight[1:-1]=abs(freq[2:]-freq[:-2])*.5
        matrices[0][t.j-1,t.i-1]=np.sum(weight[:,None]*u,axis=0);matrices[0][t.i-1,t.j-1]=np.sum(weight[:,None]*v,axis=0)
        # Split at both saved-J and source-alpha knots. Source Hermite alpha is
        # evaluated by the compiled oracle, not reimplemented or clipped.
        panels=np.unique(np.r_[x,table[:,0]]);f=cc/(panels[::-1]*1e-8)
        for order in [4,8,16]:
            nodes,weights=np.polynomial.legendre.leggauss(order);up=np.zeros(56);down=up.copy()
            for first in range(0,len(f)-1,256):
                lower=f[first:first+256];upper=f[first+1:first+257];lower=lower[:len(upper)];span=upper-lower
                qfreq=(lower[:,None]+span[:,None]*(nodes[None,:]+1)*.5).ravel();qw=cc/qfreq/1e-8
                al=alpha_at(table,qw);jn=np.stack([np.interp(qw,x,field[:,depth])for depth in range(56)],axis=1)
                uu,vv=integrands(qw,al,jn);measure=(span[:,None]*weights[None,:]*.5).ravel()[:,None];up+=np.sum(measure*uu,axis=0);down+=np.sum(measure*vv,axis=0)
            matrices[order][t.j-1,t.i-1]=up;matrices[order][t.i-1,t.j-1]=down
        for depth in range(56):
            row=dict(lower_1based=t.i,upper_1based=t.j,depth_top_to_bottom=depth,native_nodes=len(x),integration_panels=len(panels)-1)
            for order in [4,8,16]:
                row[f'up_relative_G{order}_vs_trapezoid']=float(matrices[order][t.j-1,t.i-1,depth]/matrices[0][t.j-1,t.i-1,depth]-1)
                row[f'down_relative_G{order}_vs_trapezoid']=float(matrices[order][t.i-1,t.j-1,depth]/matrices[0][t.i-1,t.j-1,depth]-1)
            rows.append(row)
    Q=d['replay_R_s_inv']+d['replay_C_s_inv'];n=d['populations'];gross=np.sum(Q*n[None,:,:],axis=1)+Q.sum(axis=0)*n
    metrics=[]
    mask=matrices[0]>0
    for order in [4,8,16]:
        delta=matrices[order]-matrices[0];balance=(np.sum(delta*n[None,:,:],axis=1)-delta.sum(axis=0)*n)/gross
        metrics.append(dict(order=order,max_rate_relative_vs_trapezoid=float(np.max(abs(matrices[order][mask]/matrices[0][mask]-1))),balance_shift_at_level40_depth3=float(balance[39,3]),max_absolute_balance_shift=float(abs(balance).max())))
    agreement=float(np.max(abs(matrices[16][mask]/matrices[8][mask]-1)))
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed quadrature audit input')
    with(a.output/'rates.csv').open('w',newline='')as out:w=csv.DictWriter(out,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    np.savez_compressed(a.output/'rates.npz',**{f'BF_G{order}_s_inv':m for order,m in matrices.items()})
    report=dict(generator=fingerprint(__file__),inputs=checks,computational_completed=True,scientific_admission='HOLD',cases=metrics,max_relative_G16_vs_G8=agreement,continua=22,depths=56,max_oracle_query=max_query,oracle_calls=counter,seconds=time.monotonic()-start,scope='Pure integral-of-interpolant diagnostic, not new radiation samples: linear physical J_nu in vacuum wavelength between saved native nodes. Actual compiled frequency-Hermite alpha and source h/k, double exponential, fixed source support/nStar. Compare native-node frequency trapezoid with4/8/16-point Gauss-Legendre panels split at all saved J and alpha knots. Streamed queries remain under existing20000limit, no HOP/RT/SE packet or larger model created. Gaussian agreement qualifies only the declared interpolant integral, not true-between-node radiation accuracy or a self-consistent model. Original native saved-rate reconstruction/FP32/stage discrepancy remains separate. No observed/abundance/physical admission.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({key:value for key,value in report.items()if key!='inputs'}))


if __name__=='__main__':main()

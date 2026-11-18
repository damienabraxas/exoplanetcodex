"""Export pinned native get_lam offsets with held source NQFAC2/QNORM10."""
import argparse,json,subprocess,sys,os
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint,load_plan_engine
    from pipeline.multi_atom import read_multi_atom
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['build','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();_,app,_,engine=load_plan_engine(dict(engine='openmulti3d',profile_variant='accurate-kernel-control'))
    atom_path=require_codex_path('experiments.bride_oxygen_atom');atom=read_multi_atom(atom_path);receipt_file=a.build/'receipt.json';receipt=json.loads(receipt_file.read_text());exe=a.build/'line_grid_probe'
    if not receipt['completed']or fingerprint(exe)!=receipt['executable']or os.environ.get('OMP_NUM_THREADS')!='1':raise ValueError('Unqualified oracle/thread')
    checks={str(receipt_file):fingerprint(receipt_file),str(exe):fingerprint(exe),str(atom_path):fingerprint(atom_path)}
    for path,fp in receipt['inputs'].items():
        if fingerprint(path)!=fp:raise ValueError('Changed oracle input')
        checks[path]=fp
    native_trans=app/'HOP/trans_mod.f90';linked=[Path(x)for x in receipt['inputs']if x.endswith('/trans_mod.f90')]
    if len(linked)!=1 or native_trans.read_bytes()!=linked[0].read_bytes():raise ValueError('Oracle is not the original pinned native grid code')
    params=require_codex_path('experiments.openmulti3d_oxygen_accurate_run')/'native/params.nml';text=params.read_text();import re
    expected={'QNORM':10.,'NQ_FAC':2.,'NQ_MIN':0.}
    for key,value in expected.items():
        match=re.search(r'\b'+key+r'\s*=\s*([^,\s]+)',text)
        if match is None or float(match[1].replace('D','e'))!=value:raise ValueError('Unqualified native '+key)
    if re.search(r'\bLEGACY_FREQL\s*=F\s*,',text)is None:raise ValueError('Unsupported legacy native frequency rule')
    checks[str(params)]=fingerprint(params);checks[str(native_trans)]=fingerprint(native_trans)
    nq=[2*(line.nq-1)+1 for line in atom.lines]
    if any(line.iw!=0 for line in atom.lines):raise ValueError('Unsupported source IW')
    a.output.mkdir(parents=True,exist_ok=False);rows=['256 10 F']+[f'{n} {line.qmax:.17e} {line.q0:.17e}'for n,line in zip(nq,atom.lines)];(a.output/'input.txt').write_text('\n'.join(rows)+'\n')
    with(a.output/'input.txt').open()as fi,(a.output/'output.txt').open('w')as fo,(a.output/'oracle.log').open('w')as log:subprocess.run([str(exe)],stdin=fi,stdout=fo,stderr=log,timeout=600,check=True)
    lines=iter((a.output/'output.txt').read_text().splitlines());grids=[]
    for k,(line,count)in enumerate(zip(atom.lines,nq)):
        i,n=map(int,next(lines).split())
        if i!=k+1 or n!=count:raise ValueError('Wrong native grid record')
        offsets=np.array([float(next(lines))for _ in range(n)])
        if not np.isfinite(offsets).all()or np.any(np.diff(offsets)<=0)or offsets[n//2]!=0:raise ValueError('Invalid native fractional grid')
        grids.append(dict(index=k,i=line.i-1,j=line.j-1,f=line.f,offsets_fraction=offsets.tolist(),nq=count))
    if next(lines,None)is not None:raise ValueError('Trailing native grid output')
    packet=dict(schema='codex.native_line_grid_control/1',atom=fingerprint(atom_path),engine=engine,grids=grids,qnorm_km_s=10.,nq_fac=2.,legacy=False,scientific_admission='HOLD',scope='Native get_lam fractional offsets applied around ORIGINAL source photon centres; no centre, f, level, continuum, collision or population override. Numerical grid control, not physical model admission.')
    (a.output/'grids.json').write_text(json.dumps(packet,indent=2)+'\n')
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed grid oracle input')
    report=dict(computational_completed=True,inputs=checks,generator=fingerprint(__file__),packet=fingerprint(a.output/'grids.json'),source_lines=len(grids),points=sum(nq),scientific_admission='HOLD');(a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['source_lines','points']}))
if __name__=='__main__':main()

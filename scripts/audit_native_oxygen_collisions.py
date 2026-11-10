"""Reconstruct native O collision matrices with native interpolation and explicit H fields."""
import argparse,json,os,re,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint,load_plan_engine
    from pipeline.openmulti3d_output import read_rate_pair
    from config.constants import require_codex_path
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['native-run','plan','collider-packet','oracle','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();plan=json.loads(a.plan.read_text());run=json.loads((a.native_run/'report.json').read_text());_,_,_,engine=load_plan_engine(dict(engine='openmulti3d',profile_variant='accurate-kernel-control'))
    if not run['computational_completed']or run['engine_receipt']!=engine or os.environ.get('OMP_NUM_THREADS')!='1':raise ValueError('Wrong qualified native execution/thread')
    checks={str(a.plan.resolve()):fingerprint(a.plan),str(a.native_run/'report.json'):fingerprint(a.native_run/'report.json')}
    for entry in run['output_files']:
        path=a.native_run/entry['path']
        if fingerprint(path)!=entry['receipt']:raise ValueError('Changed native output')
        checks[str(path)]=entry['receipt']
    table=a.native_run/'native/atoms/O.txt';lines=table.read_text().splitlines();header=next(i for i,x in enumerate(lines)if x.split()[:3]==['key','i','j']);grid=np.array(list(map(float,lines[header].split()[3:])));records=[x.split()for x in lines[header+1:]]
    if len(records)!=1106 or len(grid)!=14:raise ValueError('Unexpected native collision table')
    level_header=next(i for i,x in enumerate(lines)if x.split()==['ev','g','ion','label']);levels=[list(map(float,x.split()[:4]))for x in lines[level_header+1:level_header+42]];energy=np.array([x[1]for x in levels]);g=np.array([x[2]for x in levels])
    params=a.native_run/'native/params.nml';text=params.read_text()
    for key in ['E_EXCITATION','E_IONISATION','H_EXCITATION','H_IONISATION','H_FORBIDDEN']:
        m=re.search(r'\b'+key+r"\s*=\s*'([^']*)'",text)
        if m is None or m[1].strip():raise ValueError('Unqualified classical fill '+key)
    for key in ['CBBESC','CBFESC','CBBHSC','CBFHSC','COLSCALE']:
        m=re.search(r'\b'+key+r'\s*=\s*([^,\s]+)',text)
        if m is None or float(m[1].replace('D','e'))!=1:raise ValueError('Unqualified collision scale')
    if re.search(r'\bCH0_HMIP\s*=\s*F\s*,',text)is None:raise ValueError('Unqualified CH0 energy correction')
    native_arrays=require_codex_path(plan['native_analysis_resource'])/'arrays.npz'
    if fingerprint(native_arrays)!=plan['native_arrays']:raise ValueError('Changed native analysis')
    with np.load(native_arrays)as d:state=d['consumed_state_cgs'][::-1].astype(float);expected_star=d['native_nstar_cm3'][:,::-1].copy()
    packet_file=a.collider_packet/'packet.npz';packet_report=a.collider_packet/'report.json';pr=json.loads(packet_report.read_text())
    if fingerprint(packet_file)!=pr['packet']or pr['inputs'].get(str(native_arrays))!=plan['native_arrays']:raise ValueError('Unmatched collider packet')
    with np.load(packet_file)as d:Hground=d['H_ground_m3'].copy()*1e-6;Htotal=d['H_I_m3'].copy()*1e-6;np.testing.assert_array_equal(d['temperature_K'],state[:,2]);np.testing.assert_array_equal(d['ne_m3'],state[:,0]*1e6)
    receipt=a.oracle/'receipt.json';build=json.loads(receipt.read_text());exe=a.oracle/'collision_interpolation_probe'
    if not build['completed']or fingerprint(exe)!=build['executable']:raise ValueError('Unqualified interpolation oracle')
    for path,fp in build['inputs'].items():
        if fingerprint(path)!=fp:raise ValueError('Changed native interpolation input')
        checks[path]=fp
    for path in [table,params,native_arrays,packet_file,packet_report,receipt,exe,ROOT/'pipeline/openmulti3d_output.py']:checks[str(path)]=fingerprint(path)
    parts=[]
    for path in sorted((a.native_run/'native/save').glob('rates_O_*.bin')):parts.append(read_rate_pair(path,path.with_name(path.name.replace('rates_','colrates_'))))
    if len(parts)!=8:raise ValueError('Missing native depth slabs')
    star=np.concatenate([p['nstar_cm3']for p in parts],axis=1)[:,::-1];native=np.concatenate([p['collision_rates_s_inv']for p in parts],axis=2)[:,:,::-1];np.testing.assert_array_equal(star,expected_star)
    if native.shape!=(41,41,56):raise ValueError('Unexpected collision matrix axes')
    a.output.mkdir(parents=True,exist_ok=False)
    values=np.array([list(map(float,r[4:]))for r in records]);query=state[:,2]
    rows=[f'1106 14 56',' '.join(f'{v:.17e}'for v in grid),' '.join(f'{v:.17e}'for v in query)]+[' '.join(f'{v:.17e}'for v in row)for row in values]
    (a.output/'input.txt').write_text('\n'.join(rows)+'\n');manifest=dict(generator=fingerprint(__file__),inputs=checks,scientific_admission='HOLD');(a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    with(a.output/'input.txt').open()as fi,(a.output/'output.txt').open('w')as fo,(a.output/'oracle.log').open('w')as log:subprocess.run([str(exe.resolve())],stdin=fi,stdout=fo,stderr=log,timeout=600,check=True)
    ek=float((a.output/'output.txt').read_text().splitlines()[0].split()[1]);data=np.loadtxt(a.output/'output.txt',skiprows=1);np.testing.assert_array_equal(data[:,:2],np.array([(i,j)for i in range(1,1107)for j in range(1,57)]));q=data[:,2].reshape(1106,56)
    from pipeline.native_oxygen_collisions import collision_matrix
    checks[str(ROOT/'pipeline/native_oxygen_collisions.py')]=fingerprint(ROOT/'pipeline/native_oxygen_collisions.py')
    controls=['native_contract','total_H_substitution','directed_CH0_only','LW_CH_direction']
    matrices={name:collision_matrix(records,q,query,state[:,0],Hground,Htotal,star,energy,g,ek,control=name)for name in controls}
    counts={key:sum(r[1]==key for r in records)for key in ['CE','CI','CH','CH0']}
    mask=native>0;stats={}
    for name,matrix in matrices.items():
        stats[name]=dict(max_relative_nonzero=float(np.max(abs(matrix[mask]/native[mask]-1))),max_absolute_s_inv=float(np.max(abs(matrix-native))),extra_nonzero_cells=int(np.count_nonzero((native==0)&(matrix!=0))))
    ii,jj=np.tril_indices(41,-1)
    upward_exact=bool(np.array_equal(matrices['native_contract'][ii,jj].astype(np.float32).astype(float),native[ii,jj]))
    bound=float(16*np.finfo(np.float32).eps);qualified=upward_exact and stats['native_contract']['max_relative_nonzero']<=bound and stats['native_contract']['extra_nonzero_cells']==0
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed collision audit input')
    np.savez_compressed(a.output/'matrices.npz',native=native,nstar_cm3=star,coefficients_cm3=q,**matrices)
    report=dict(**manifest,counts=counts,native_e_k=ek,numerical_bound=bound,numerical_contract_pass=qualified,stored_upward_float32_bit_identical=upward_exact,stored_upward_cells=len(ii)*56,statistics=stats,scope='Native stored upward collision reconstruction and labelled negative controls. Downward stored-reader rates are inferred through nStar, not independently saved rates. Not fresh RT/SE or physical collision-data adoption. Directed-only CH0 and alternate CH use the same native interpolation to isolate semantics, not reproduce a different full engine.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(pass_contract=qualified,statistics=stats)))

if __name__=='__main__':main()

"""Fresh bounded build of actual pinned MULTI DAMP; no physical formula substitution."""
import argparse,hashlib,json,os,platform,re,shutil,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));PIN='4d6aa3c39d634b38da772e469a5692435eba1354'

def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if platform.system()!='Darwin'or os.environ.get('OMP_NUM_THREADS')!='1':raise ValueError('Require one-thread Mac oracle')
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.source_root,text=True).strip()!=PIN:raise ValueError('Changed original MULTI pin')
    rel=Path('pymulti/files/multi/mul23/source');source=a.source_root/rel;inputs={}
    def read(name):
        path=source/name;raw=subprocess.check_output(['git','show',f'{PIN}:{rel/name}'],cwd=a.source_root)
        if raw!=path.read_bytes():raise ValueError('Changed source '+name)
        inputs[str(path)]=fingerprint(path);return raw.decode('latin1')
    text=read('mul23_sub.f');match=re.search(r'^      SUBROUTINE DAMP\s*$',text,re.M);end=re.search(r'^      END\s*$',text[match.start():],re.M);routine=text[match.start():match.start()+end.end()]+'\n';bconst=read('mul23_bconst.f');read('make_mul23_dbl.csh');driver=ROOT/'native/multi_width_oracle/driver.f';inputs[str(driver)]=fingerprint(driver);compiler=Path(shutil.which('gfortran'));inputs[str(compiler)]=fingerprint(compiler);a.output.mkdir(parents=True,exist_ok=False);start=time.perf_counter()
    for name in ['CATOM','CATMOS','CATMO2','CCONST','CINPUT','CLU','CTERM','CTRAN','CSLINE','CBCNST','C3D']:(a.output/name).write_text(re.sub(r'\bREAL ', 'REAL(8) ',read(name)))
    (a.output/'PREC').write_text(read('PREC.dp'));param=read('PARAM')
    if 'MK=30'not in param or'MLINE=150'not in param:raise ValueError('Changed source capacities')
    (a.output/'PARAM').write_text(param.replace('MK=30','MK=64').replace('MLINE=150','MLINE=512'))
    generated=a.output/'multi_damp_dp.f';generated.write_text(routine+bconst+driver.read_text(),encoding='latin1');exe=a.output.resolve()/'multi_width_probe';cmd=[str(compiler),'-std=legacy','-O2','-fcheck=bounds','-ffixed-line-length-132',str(generated.resolve()),'-o',str(exe)]
    manifest=dict(generator=fingerprint(__file__),inputs=inputs,source_commit=PIN,routine_sha256=hashlib.sha256(routine.encode('latin1')).hexdigest(),generated=fingerprint(generated),command=cmd,seconds_limit=600,threads=1,scientific_admission='HOLD',completed=False);(a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    with(a.output/'compiler.log').open('w')as f:result=subprocess.run(cmd,cwd=a.output,stdout=f,stderr=subprocess.STDOUT,timeout=600,check=False)
    manifest.update(seconds=time.perf_counter()-start,returncode=result.returncode,completed=result.returncode==0)
    if result.returncode==0:manifest['executable']=fingerprint(exe)
    for path,f in inputs.items():
        if fingerprint(path)!=f:raise ValueError('Changed build source')
    (a.output/'receipt.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(dict(completed=manifest['completed'],seconds=manifest['seconds'])))
    if result.returncode:raise RuntimeError('Native DAMP oracle build refused, inspect preserved compiler.log')
if __name__=='__main__':main()

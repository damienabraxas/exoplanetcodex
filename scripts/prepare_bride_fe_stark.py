"""Check supplied and Cowley Fe I Stark equations through unchanged TS fragment."""
import argparse,json,subprocess,shutil,platform,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def source_rows():
    from scripts.prepare_bride_fe_lines import read_packet
    packet=read_packet(ROOT/'data/results/rya1281/trace_fe_v1/lines.json')
    pressure=json.loads((ROOT/'data/results/rya1281/fe_pressure_v1/source/report.json').read_text())
    return [dict(line_id=l['line_id'],wavelength_air_A=float(e['source_row']['wave_A']),excitation_eV=l['excitation_eV'],
        ionization_eV=pressure['ionization_eV'],log_gamma_stark_10000=float(e['source_row']['stark']))for l,e in zip(packet['lines'],packet['evidence'])]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Mac only')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_fe_stark import stark_gamma
    import lightweaver.constants as c
    root=require_codex_path('engines.ispec');pin='21269948750e2177fe32c5df38e1d9bb81363374'
    sourcefiles=['synthesizer/turbospectrum/source/depth.f','synthesizer/turbospectrum/source/bsyn.f','ispec/lines.py']
    sources={}
    for rel in sourcefiles:
        if subprocess.check_output(['git','show',pin+':'+rel],cwd=root)!=(root/rel).read_bytes():raise ValueError('Changed source '+rel)
        sources[rel]=fingerprint(root/rel)
    text=(root/sourcefiles[0]).read_text();start=text.index('          if (ijon.eq.1) then\n            eion=chi')
    end=text.index('          gamma=gamma+gamstark',start)+len('          gamma=gamma+gamstark')
    fragment=text[start:end]
    driver='''      program bride_stark
      implicit none
      integer n,j,m,i,k,l,ijon
      real temp,pe(1),boltz,gamst,gamstark,gamma
      real chiu,chie,chi,chi2,chi3,xl,eion,eudiff,tempstark,xne
      real air,logst,nelist(3),tlist(3)
      data tlist/4000.,5772.,10000./
      data nelist/0.,1.e12,1.e14/
      read(*,*) n,boltz
      ijon=1
      chi2=0.
      chi3=0.
      i=1
      do j=1,n
      read(*,*) air,chie,chi,logst
      xl=air*1.e-8
      do m=1,2
      do k=1,3
      temp=tlist(k)
      do l=1,3
      xne=nelist(l)
      pe(1)=xne*boltz*temp
      gamst=0.
      if(m.eq.1)gamst=logst
      gamma=0.
'''+fragment+'''
      write(*,*) 'BRIDE_STARK',j,m,k,l,gamma
      enddo
      enddo
      enddo
      enddo
      end
'''
    rows=source_rows();args.output.mkdir(parents=True,exist_ok=False);(args.output/'driver.f').write_text(driver)
    compiler=shutil.which('gfortran')
    if not compiler:raise RuntimeError('Existing gfortran required')
    command=[compiler,'-O2','-fno-automatic','-std=legacy','-ffixed-line-length-none','-fcheck=bounds','-o','stark_oracle','driver.f']
    run=subprocess.run(command,cwd=args.output,capture_output=True,text=True,timeout=60)
    (args.output/'build.log').write_text(run.stdout+run.stderr)
    if run.returncode:raise RuntimeError('Stark driver compile failed')
    inputs=f'{len(rows)} {c.KBoltzmann*1e7:.14e}\n'+''.join(f'{r["wavelength_air_A"]} {r["excitation_eV"]} {r["ionization_eV"]} {r["log_gamma_stark_10000"]}\n'for r in rows)
    (args.output/'input.txt').write_text(inputs)
    run=subprocess.run([str(args.output/'stark_oracle')],cwd=args.output,input=inputs,capture_output=True,text=True,timeout=30)
    (args.output/'stdout.txt').write_text(run.stdout);(args.output/'stderr.txt').write_text(run.stderr)
    if run.returncode:raise RuntimeError('Stark oracle failed')
    errors=[];cases=[]
    for line in run.stdout.splitlines():
        parts=line.split()
        if not parts or parts[0]!='BRIDE_STARK':continue
        j,m,k,i=map(int,parts[1:5]);native=float(parts[5]);T=[4000.,5772.,10000.][k-1];ne=[0.,1e18,1e20][i-1];mode=['provided','cowley'][m-1]
        py=float(stark_gamma(np.array([T]),np.array([ne]),rows[j-1],mode)[0])
        error=abs(py/native-1)if native else abs(py)
        errors.append(error);cases.append(dict(line_id=rows[j-1]['line_id'],mode=mode,temperature_K=T,ne_m3=ne,native_gamma_s_inv=native,python_gamma_s_inv=py))
    if len(cases)!=41*2*3*3 or max(errors)>128*np.finfo(np.float32).eps:raise ValueError('Incomplete/failed native Stark arithmetic check')
    result=dict(schema='codex.bride_fe_stark_source/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        kernel=fingerprint(ROOT/'pipeline/bride_fe_stark.py'),source_commit=pin,sources=sources,rows=rows,cases=cases,
        driver=fingerprint(args.output/'driver.f'),binary=fingerprint(args.output/'stark_oracle'),input=fingerprint(args.output/'input.txt'),native_stdout=fingerprint(args.output/'stdout.txt'),
        compiler=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0],command=command,max_relative_error=max(errors),
        interpretation='Neutral Fe only. Supplied negative-log width and native Cowley fallback are distinct interventions; no full TS EOS/export parity or physical accuracy claimed.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(len(cases),'native controls; max relative error',max(errors))


if __name__=='__main__':main()

"""Execute pinned native FeII pressure/Stark branches in an isolated driver."""
import argparse,csv,json,re,shutil,subprocess,sys,platform
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Mac-only source audit')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_fe_pressure import pressure_gamma
    from pipeline.bride_feii_profile import feii_stark_gamma
    from scripts.run_bride_feii_transfer import selected_line
    import lightweaver.constants as c
    from lightweaver.atomic_table import PeriodicTable
    plan=json.loads((ROOT/'config/experiments/rya1281_feii_transfer.json').read_text());line,evidence=selected_line(plan);src=evidence['source']
    root=require_codex_path('engines.ispec');pin='21269948750e2177fe32c5df38e1d9bb81363374'
    paths=['synthesizer/turbospectrum/source/anstee.f','synthesizer/turbospectrum/source/depth.f','synthesizer/turbospectrum/source/bsyn.f',
           'synthesizer/turbospectrum/DATA/Anstee-1802-2014.dat','synthesizer/turbospectrum/DATA/IRWIN_atoms_v07.3.dat','ispec/lines.py']
    sources={}
    for rel in paths:
        if subprocess.check_output(['git','show',pin+':'+rel],cwd=root)!=(root/rel).read_bytes():raise ValueError('Changed pinned source '+rel)
        sources[rel]=fingerprint(root/rel)
    energies=re.findall(r"^\s*\d+\s+26\.02\s+'[^']+'\s+(\S+)",(root/paths[4]).read_text(),re.M)
    if len(energies)!=1:raise ValueError('Ambiguous FeII ionization energy')
    ionization=float(energies[0]);depth=(root/paths[1]).read_text()
    wblock=depth[depth.index('          gam6= 10.**gamvdw10000'):depth.index("        else if (recipe.eq.'A'")]
    ablock=depth[depth.index('          sigma=sigmacross*2.8002e-17'):depth.index("        else\n          print*,'DEPTH. WARNING!!'")]
    starkblock=depth[depth.index('          if (ijon.eq.1) then\n            eion=chi'):depth.index('          gamma=gamma+gamstark')+len('          gamma=gamma+gamstark')]
    tail=depth[depth.index('      FUNCTION gammln'):]
    driver='''      program feii_width
      implicit none
      integer i,j,k,m,l,ijon
      real xlb,elo,eion,sigmacross,velexp,fdamp,logst
      real boltz,mh,ma,fyrapi,redmas,constb,sigma,relvel
      real gam6,gam6anstee,gamma,gamrad,gamvdw10000,gammln
      real t(1),ph(1),phe(1),ph2(1),pe(1)
      real temp,gamst,gamstark,chiu,chie,chi,chi2,chi3,xl
      real eudiff,tempstark,xne,basegamma
      real tlist(3),hlist(3),elist(3)
      character lo*1,up*1,recipe*1
      data tlist/4000.,5772.,10000./
      data hlist/1.e13,1.e15,1.e17/
      data elist/0.,1.e12,1.e14/
      read(*,*) boltz,mh,ma
      read(*,*) xlb,elo,chi2,fdamp,gamrad,logst,lo,up
      call anstee(26,2,xlb,elo,chi2,sigmacross,velexp,
     & lo,up,recipe)
      if(recipe.eq.'U'.and.fdamp.lt.0.)recipe='W'
      write(*,*) 'FEII_RECIPE',recipe,sigmacross,velexp
      if(recipe.eq.'U')stop 'unsupported recipe'
      i=1
      ijon=2
      chi=0.
      chi3=0.
      chie=elo
      xl=xlb*1.e-8
      fyrapi=16.*atan(1.)
      redmas=1./((1./mh)+(1./ma))
      constb=8.*boltz*4.0/fyrapi/redmas
      do m=1,2
      do j=1,3
      do k=1,3
      t(i)=tlist(j)
      temp=t(i)
      ph(i)=hlist(k)*boltz*temp
      phe(i)=ph(i)*.1
      ph2(i)=0.
      gamvdw10000=fdamp
      if(recipe.eq.'W')then
'''+wblock+'''      else
'''+ablock+'''      endif
      basegamma=gamma
      do ijon=2,2
      do l=1,3
      xne=elist(l)
      pe(1)=xne*boltz*temp
      gamst=0.
      if(m.eq.1)gamst=logst
      gamma=basegamma
'''+starkblock+'''
      write(*,*) 'FEII_WIDTH',m,j,k,l,basegamma,gamma
      enddo
      enddo
      i=1
      enddo
      enddo
      enddo
      end
'''+tail
    # Exact source blocks; i remains the single atmosphere row, l indexes tests.
    args.output.mkdir(parents=True,exist_ok=False);(args.output/'DATA').mkdir()
    (args.output/'driver.f').write_text(driver)
    shutil.copyfile(root/paths[0],args.output/'anstee.f');shutil.copyfile(root/paths[3],args.output/'DATA/Anstee-1802-2014.dat')
    compiler=shutil.which('gfortran')
    if compiler is None:raise RuntimeError('Existing gfortran required')
    command=[compiler,'-O2','-fno-automatic','-std=legacy','-ffixed-line-length-none','-fcheck=bounds','-o','feii_oracle','driver.f','anstee.f']
    build=subprocess.run(command,cwd=args.output,capture_output=True,text=True,timeout=60);(args.output/'build.log').write_text(build.stdout+build.stderr)
    if build.returncode:raise RuntimeError('FeII driver compile failed')
    text=f'{c.KBoltzmann*1e7:.14e} {PeriodicTable["H"].mass*c.Amu*1000:.14e} {PeriodicTable["Fe"].mass*c.Amu*1000:.14e}\n'
    text+=f'{src["wave_A"]} {src["lower_state_eV"]} {ionization} {src["turbospectrum_fdamp"]} {src["turbospectrum_rad"]} {src["stark"]} {src["lower_orbital_type"]} {src["upper_orbital_type"]}\n'
    (args.output/'input.txt').write_text(text)
    run=subprocess.run([str(args.output/'feii_oracle')],cwd=args.output,input=text,capture_output=True,text=True,timeout=30)
    (args.output/'stdout.txt').write_text(run.stdout);(args.output/'stderr.txt').write_text(run.stderr)
    if run.returncode:raise RuntimeError('FeII native execution failed')
    native=[s.split() for s in run.stdout.splitlines() if s.split() and s.split()[0]=='FEII_RECIPE']
    if len(native)!=1:raise ValueError('Missing native recipe')
    recipe=dict(recipe=native[0][1],sigma_au=float(native[0][2]),alpha=float(native[0][3]),log_gamma6_10000=float(src['turbospectrum_fdamp']))
    starkrow=dict(ionization_eV=ionization,excitation_eV=float(src['lower_state_eV']),wavelength_air_A=float(src['wave_A']),log_gamma_stark_10000=float(src['stark']))
    checks=[]
    for s in run.stdout.splitlines():
        fields=s.split()
        if not fields or fields[0]!='FEII_WIDTH':continue
        mode,ti,hi,ei=map(int,fields[1:5]);T=[4000.,5772.,10000.][ti-1];H=[1e19,1e21,1e23][hi-1];ne=[0.,1e18,1e20][ei-1]
        base=float(fields[5]);native=float(fields[6]);method=['provided','cowley'][mode-1]
        pressure=float(pressure_gamma(T,H,.1*H,0.,recipe));st=float(feii_stark_gamma(T,ne,starkrow,method))
        py=line['gamma_rad_s_inv']+pressure+st
        checks.append(dict(mode=method,temperature_K=T,neutral_H_m3=H,ne_m3=ne,native_base_gamma=base,native_gamma=native,python_gamma=py,relative_error=abs(py/native-1)))
    if len(checks)!=54 or max(r['relative_error'] for r in checks)>128*np.finfo(np.float32).eps:raise ValueError('Native width arithmetic parity failed')
    widths=[];frozen=json.loads((ROOT/'config/experiments/rya1281_continuum_anchor.json').read_text());origins=[]
    for case in frozen['cases']:
        directory=require_codex_path(case['resource'])/case['relative'];r=json.loads((directory/'report.json').read_text())
        if fingerprint(directory/'neutral_colliders.csv')!=r['neutral_colliders']:raise ValueError('Changed colliders')
        fields=r['local_binding']['local_binding'];T=np.array(fields['fields']['temperature_K']);ne=np.array(fields['ne_m3'])
        coll=np.genfromtxt(directory/'neutral_colliders.csv',delimiter=',',names=True)
        pg=pressure_gamma(T,coll['neutral_H_m3'],coll['neutral_He_m3'],np.zeros_like(T),recipe)
        supplied=feii_stark_gamma(T,ne,starkrow,'provided');cowley=feii_stark_gamma(T,ne,starkrow,'cowley')
        for i in range(len(T)):
            widths.append(dict(case=case['name'],depth_index=i,logtau500=float(coll['logtau500'][i]),temperature_K=float(T[i]),
                gamma_natural=line['gamma_rad_s_inv'],gamma_pressure=float(pg[i]),gamma_stark_provided=float(supplied[i]),gamma_stark_cowley=float(cowley[i]),
                total_provided=float(line['gamma_rad_s_inv']+pg[i]+supplied[i]),total_export_cowley=float(line['gamma_rad_s_inv']+pg[i]+cowley[i])))
        origins.append(dict(case=case['name'],report=fingerprint(directory/'report.json'),colliders=fingerprint(directory/'neutral_colliders.csv')))
    with (args.output/'widths.csv').open('x') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(widths[0]),lineterminator='\n');writer.writeheader();writer.writerows(widths)
    result=dict(schema='codex.bride_feii_broadening/1',scientific_admission='HOLD',generator=fingerprint(__file__),line=line,evidence=evidence,
        source_commit=pin,sources=sources,recipe=recipe,stark_row=starkrow,native_checks=checks,max_relative_error=max(r['relative_error'] for r in checks),
        binary=fingerprint(args.output/'feii_oracle'),driver=fingerprint(args.output/'driver.f'),input=fingerprint(args.output/'input.txt'),stdout=fingerprint(args.output/'stdout.txt'),
        compiler=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0],command=command,origins=origins,widths=fingerprint(args.output/'widths.csv'),
        scope='Native branch arithmetic and actual-column widths only; no RT or physical broadening accuracy/adoption. H2 omitted. Provided Stark and export-Cowley retained separately.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(recipe=recipe,native_checks=len(checks),max_relative_error=result['max_relative_error'])))


if __name__=='__main__':main()

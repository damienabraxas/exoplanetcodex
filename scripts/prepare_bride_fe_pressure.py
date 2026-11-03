"""Execute isolated pinned TS Anstee selection and audit native pressure widths."""
import argparse,json,re,shutil,subprocess,sys,platform
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Mac-only diagnostic')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from scripts.prepare_bride_fe_lines import read_packet
    from pipeline.bride_fe_pressure import pressure_gamma
    from lightweaver.atomic_table import PeriodicTable
    import lightweaver.constants as c
    packet_path=ROOT/'data/results/rya1281/trace_fe_v1/lines.json';packet=read_packet(packet_path)
    root=require_codex_path('engines.ispec');commit='21269948750e2177fe32c5df38e1d9bb81363374'
    files=['synthesizer/turbospectrum/source/anstee.f','synthesizer/turbospectrum/source/depth.f','synthesizer/turbospectrum/source/bsyn.f','synthesizer/turbospectrum/DATA/Anstee-1802-2014.dat','synthesizer/turbospectrum/DATA/IRWIN_atoms_v07.3.dat']
    sources={}
    for rel in files:
        path=root/rel
        if subprocess.check_output(['git','show',commit+':'+rel],cwd=root)!=path.read_bytes():raise ValueError('Changed source: '+rel)
        sources[rel]=fingerprint(path)
    irwin=(root/files[-1]).read_text();matches=re.findall(r"^\s*\d+\s+26\.01\s+'[^']+'\s+(\S+)",irwin,re.M)
    if len(matches)!=1:raise ValueError('Ambiguous original Fe ionization energy')
    ionization=float(matches[0]);depth=(root/files[1]).read_text()
    # Exact source fragments; driver sets cgs collider pressures from explicit densities.
    wblock=depth[depth.index('          gam6= 10.**gamvdw10000'):depth.index("        else if (recipe.eq.'A'")]
    ablock=depth[depth.index('          sigma=sigmacross*2.8002e-17'):depth.index("        else\n          print*,'DEPTH. WARNING!!'")]
    gammln=depth[depth.index('      FUNCTION gammln'):]
    mh=PeriodicTable['H'].mass*c.Amu*1000;ma=PeriodicTable['Fe'].mass*c.Amu*1000;k=c.KBoltzmann*1e7
    driver='''      program bride_pressure
      implicit none
      integer i,j,n
      real xlb,elo,eion,sigmacross,velexp,fdamp
      real boltz,mh,ma,fyrapi,redmas,constb,sigma,relvel
      real gam6,gam6anstee,gamma,gamrad,gamvdw10000,gammln
      real t(1),ph(1),phe(1),ph2(1),nh,nhe,nh2
      character lo*1,up*1,recipe*1
      read(*,*) n,boltz,mh,ma
      fyrapi=16.*atan(1.)
      redmas=1./((1./mh)+(1./ma))
      constb=8.*boltz*4.0/fyrapi/redmas
      gamrad=0.
      do j=1,n
      read(*,*) xlb,elo,eion,lo,up,fdamp
      call anstee(26,1,xlb,elo,eion,sigmacross,velexp,
     & lo,up,recipe)
      if(recipe.eq.'U'.and.fdamp.lt.0.) recipe='W'
      write(*,*) 'BRIDE_RECIPE',j,recipe,sigmacross,velexp
      do i=1,1
      t(i)=5772.
      nh=1.e17
      nhe=1.e16
      nh2=1.e14
      ph(i)=nh*boltz*t(i)
      phe(i)=nhe*boltz*t(i)
      ph2(i)=nh2*boltz*t(i)
      gamvdw10000=fdamp
      if(recipe.eq.'W')then
'''+wblock+'''      else
'''+ablock+'''      endif
      write(*,*) 'BRIDE_WIDTH',j,gamma
      enddo
      enddo
      end
'''+gammln
    args.output.mkdir(parents=True,exist_ok=False);(args.output/'DATA').mkdir()
    shutil.copyfile(root/files[0],args.output/'anstee.f');shutil.copyfile(root/files[3],args.output/'DATA/Anstee-1802-2014.dat')
    (args.output/'driver.f').write_text(driver)
    compiler=shutil.which('gfortran')
    if compiler is None:raise RuntimeError('Existing gfortran required')
    command=[compiler,'-O2','-fno-automatic','-std=legacy','-ffixed-line-length-none','-fcheck=bounds','-o','pressure_oracle','driver.f','anstee.f']
    build=subprocess.run(command,cwd=args.output,capture_output=True,text=True,timeout=60)
    (args.output/'build.log').write_text(build.stdout+build.stderr)
    if build.returncode:raise RuntimeError('Oracle compile failed')
    inputs=f'{len(packet["lines"])} {k:.14e} {mh:.14e} {ma:.14e}\n'
    for e in packet['evidence']:
        r=e['source_row'];inputs+=f'{r["wave_A"]} {r["lower_state_eV"]} {ionization} {r["lower_orbital_type"]} {r["upper_orbital_type"]} {r["turbospectrum_fdamp"]}\n'
    (args.output/'input.txt').write_text(inputs)
    run=subprocess.run([str(args.output/'pressure_oracle')],input=inputs,cwd=args.output,capture_output=True,text=True,timeout=30)
    (args.output/'stdout.txt').write_text(run.stdout);(args.output/'stderr.txt').write_text(run.stderr)
    if run.returncode:raise RuntimeError('Oracle execution failed')
    recipes={};widths={}
    for line in run.stdout.splitlines():
        parts=line.split()
        if parts and parts[0]=='BRIDE_RECIPE':recipes[int(parts[1])]=(parts[2],float(parts[3]),float(parts[4]))
        if parts and parts[0]=='BRIDE_WIDTH':widths[int(parts[1])]=float(parts[2])
    if len(recipes)!=len(packet['lines']) or len(widths)!=len(recipes):raise ValueError('Incomplete native output')
    selected=[];errors=[]
    for i,(line,e) in enumerate(zip(packet['lines'],packet['evidence']),1):
        recipe,sigma,alpha=recipes[i]
        if recipe not in ('A','B','C','S','W'):raise ValueError('Unsupported actual branch')
        row=dict(line_id=line['line_id'],recipe=recipe,sigma_au=sigma,alpha=alpha,log_gamma6_10000=float(e['source_row']['turbospectrum_fdamp']))
        py=float(pressure_gamma(np.array([5772.]),np.array([1e23]),np.array([1e22]),np.array([1e20]),row)[0])
        error=abs(py/widths[i]-1);errors.append(error);row['native_reference_gamma_s_inv']=widths[i];selected.append(row)
    # Native default REAL: engineering arithmetic check, not a physical accuracy gate.
    guard=128*np.finfo(np.float32).eps
    if max(errors)>guard:raise ValueError('Pressure port fails native REAL arithmetic guard')
    result=dict(schema='codex.bride_fe_pressure_source/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        pressure_kernel=fingerprint(ROOT/'pipeline/bride_fe_pressure.py'),lines=fingerprint(packet_path),source_commit=commit,sources=sources,
        driver=fingerprint(args.output/'driver.f'),binary=fingerprint(args.output/'pressure_oracle'),compiler=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0],command=command,
        input=fingerprint(args.output/'input.txt'),native_stdout=fingerprint(args.output/'stdout.txt'),ionization_eV=ionization,
        constants_cgs=dict(k=k,mH=mh,mFe=ma),recipes=selected,native_width_max_relative_error=max(errors),native_real_guard=float(guard),
        caveats='Native anstee selection/pressure fragments only, not full TS spectrum/EOS. Supplied negative fdamp is secondary to anstee branch. Pressure recipe includes H/He/H2; actual HHe-only control explicitly omits H2/Stark.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    from collections import Counter
    print('Actual recipes',Counter(r['recipe']for r in selected),'max width arithmetic error',max(errors))


if __name__=='__main__':main()

"""Fixed-state CN coordinate factorial; no spectrum or chemistry/provider run."""
import argparse,csv,json,os,shutil,subprocess,sys,time
from pathlib import Path
import numpy as np
from scipy.special import wofz
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
WRAPPER='''program probe
 implicit none
 integer count,i
 double precision wave,centre
 real c,xl,width,a,vt,v,newvoigt
 read(*,*) count
 if(count<1.or.count>20000)stop 2
 do i=1,count
  read(*,*)wave,centre,c,xl,width,a
  vt=(wave-centre)*1.d-8
  vt=c*vt/xl**2
  v=vt/width
  write(*,'(I8,1P,2E25.16)')i,v,newvoigt(a,v)
 enddo
end program probe
'''


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 from config.constants import require_codex_path
 from pipeline.bride_reference import fingerprint
 from pipeline.bride_cn_diagnostic import frequency_offsets,external_state_extinction
 from pipeline.wavelength_util import air_to_vac,vac_to_air
 from pipeline import wavelength_util
 plan_path=ROOT/'config/experiments/rya1281_cn_coordinates.json';plan=json.loads(plan_path.read_text());source_dir=require_codex_path('experiments.bride_cn_coefficients');archive=ROOT/'data/results/rya1281/cn_coefficients_v1'
 for name in ['report.json','CN/states.csv','all/states.csv']:
  if(source_dir/name).read_bytes()!=(archive/name).read_bytes():raise ValueError('Changed accepted native state/report')
 states=np.genfromtxt(source_dir/'CN/states.csv',delimiter=',',names=True);n=len(states)
 if n!=168:raise ValueError('Unexpected CN state shape')
 engine=require_codex_path('engines.ispec');native=engine/'synthesizer/turbospectrum/source/newvoigt.f';native_bytes=native.read_bytes()
 if native_bytes!=subprocess.check_output(['git','show',plan['source_commit']+':synthesizer/turbospectrum/source/newvoigt.f'],cwd=engine):raise ValueError('Native Voigt changed')
 a.output.mkdir(parents=True,exist_ok=False);(a.output/'probe.f90').write_text(WRAPPER);(a.output/'newvoigt.f').write_bytes(native_bytes)
 compiler=shutil.which('gfortran')
 if compiler is None:raise RuntimeError('Existing compiler required')
 program=a.output/'probe';command=[compiler,'-O3','-fno-automatic','-std=legacy','-o',str(program),'probe.f90','newvoigt.f'];before=time.perf_counter();r=subprocess.run(command,cwd=a.output,capture_output=True,text=True,timeout=plan['compile_seconds']);(a.output/'compiler.log').write_text(r.stdout+r.stderr)
 if r.returncode:raise RuntimeError('Coordinate probe compile failed')
 compile_seconds=time.perf_counter()-before
 probe=np.asarray(plan['native_probe_air_A']);pairs=[(wave,row)for row in states for wave in probe]
 text=str(len(pairs))+'\n'+''.join(' '.join(f'{float(x):.17g}'for x in [wave,row['wavelength_air_A'],row['c_cm_s'],row['wavelength_cm'],row['doppler_width_Hz'],row['voigt_a']])+'\n'for wave,row in pairs);(a.output/'input.txt').write_text(text)
 before=time.perf_counter();r=subprocess.run([str(program)],input=text,cwd=a.output,capture_output=True,text=True,timeout=plan['probe_seconds'],env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'));elapsed=time.perf_counter()-before;(a.output/'stdout.txt').write_text(r.stdout);(a.output/'stderr.txt').write_text(r.stderr)
 raw=np.array([[float(x)for x in line.split()]for line in r.stdout.splitlines()])
 if r.returncode or raw.shape!=(len(pairs),3)or not np.isfinite(raw).all():raise ValueError('Invalid coordinate probe')
 np.testing.assert_array_equal(raw[:,0],np.arange(1,len(pairs)+1));v_native=raw[:,1].reshape(n,len(probe)).T
 args=(states['wavelength_air_A'],states['doppler_width_Hz'],states['c_cm_s'],states['wavelength_cm'])
 v_python=frequency_offsets(probe,*args,input_medium='air',mode='native_linear_air');mask=v_native!=0
 if np.any(v_python[~mask]!=0):raise ValueError('Native/Python line centres disagree')
 mapping_error=float(np.max(np.abs(v_python[mask]/v_native[mask]-1)));bound=float(plan['native_arithmetic_eps_multiple']*np.finfo(np.float32).eps)
 lo,hi=plan['range_air_A'];wave=np.linspace(lo,hi,int(round((hi-lo)/plan['step_air_A']))+1);vac=air_to_vac(wave);centrevac=air_to_vac(states['wavelength_air_A']);roundtrip=vac_to_air(vac)-wave
 strength=np.pi*states['e_esu']**2/(states['electron_mass_g']*states['c_cm_s'])*1e-4;pop=states['lower_population_per_g']*states['rho_g_cm3']*1e6
 outputs={};offsets={}
 for case in plan['cases']:
  width=states['doppler_width_Hz'].copy();damping=states['voigt_a'].copy();stim=states['stim'].copy()
  if case=='native_linear_air':grid=wave;medium='air';mode=case
  elif case=='exact_air_frequency':grid=wave;medium='air';mode=case
  else:
   grid=vac;medium='vacuum';mode='exact_vacuum_frequency'
   if case=='exact_vacuum_rescaled_width_a_stim':
    ratio=states['wavelength_cm']/(centrevac*1e-8);width*=ratio;damping/=ratio
    stim=-np.expm1(-states['h_erg_s']*states['c_cm_s']/(states['boltz_erg_K']*centrevac*1e-8*states['temperature_K']))
  v=frequency_offsets(grid,states['wavelength_air_A'],width,states['c_cm_s'],states['wavelength_cm'],input_medium=medium,mode=mode);offsets[case]=v
  outputs[case]=external_state_extinction(v,damping,width,pop,states['gf'],stim,strength)
 summary=[];baseline=outputs['native_linear_air']
 for case,value in outputs.items():
  for j,row in enumerate(states):
   line=row['wavelength_air_A'];select=np.abs(wave-line)<=.5;area=float(np.trapezoid(value[select,j],wave[select]));basearea=float(np.trapezoid(baseline[select,j],wave[select]))
   summary.append(dict(case=case,state_index=j,C=int(row['C']),N=int(row['N']),depth_index=int(row['depth_index']),line_air_A=line,max_absolute_change_over_baseline_peak=float(np.max(np.abs(value[:,j]-baseline[:,j]))/np.max(baseline[:,j])),max_pointwise_relative_change=float(np.max(np.abs(value[:,j]/baseline[:,j]-1))),finite_air_window_integral_m_inv_A=area,relative_integral_change=area/basearea-1))
 with(a.output/'summary.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(summary[0]),lineterminator='\n');w.writeheader();w.writerows(summary)
 with(a.output/'probe_comparison.csv').open('w',newline='')as f:
  w=csv.writer(f,lineterminator='\n');w.writerow(['state_index','wave_air_A','native_v','python_v','native_H','faddeeva_H'])
  w.writerows((j,float(probe[i]),v_native[i,j],v_python[i,j],raw[j*len(probe)+i,2],wofz(v_python[i,j]+1j*states['voigt_a'][j]).real)for j in range(n)for i in range(len(probe)))
 maxima={case:{key:max(abs(row[key])for row in summary if row['case']==case)for key in ['max_absolute_change_over_baseline_peak','max_pointwise_relative_change','relative_integral_change']}for case in plan['cases']}
 if native.read_bytes()!=native_bytes or(source_dir/'CN/states.csv').read_bytes()!=(archive/'CN/states.csv').read_bytes():raise ValueError('Protected source/state changed')
 report=dict(schema='codex.bride_cn_coordinates/1',scientific_admission='HOLD',generator=fingerprint(__file__),plan=fingerprint(plan_path),kernel=fingerprint(ROOT/'pipeline/bride_cn_diagnostic.py'),converter=fingerprint(wavelength_util.__file__),native_source=fingerprint(native),states=fingerprint(source_dir/'CN/states.csv'),state_report=fingerprint(source_dir/'report.json'),wrapper=fingerprint(a.output/'probe.f90'),program=fingerprint(program),compiler_command=command,compile_seconds=compile_seconds,probe_seconds=elapsed,native_calls=len(pairs),mapping_max_relative_difference=mapping_error,prospective_arithmetic_bound=bound,mapping_passed=mapping_error<=bound,grid_rows=len(wave),roundtrip_max_abs_A=float(np.max(np.abs(roundtrip))),line_air_to_vac_shift_A=(centrevac-states['wavelength_air_A']).tolist(),summary_maxima=maxima,
 limitations=['Profiles are extinction m^-1 on fixed external native populations, not emergent flux/EWs or abundance.','AIR->vacuum through canonical owner once; approximate inverse residual measured, no iterative replacement or second conversion.','Width/a/stim rescaling is a separate coupled diagnostic case, not a correction/adopted model. Native constants retained.','Legacy consumer AIR coordinates do not qualify primary line-list frame. No velocity/RV or observed-spectrum/telluric manipulation.','No provider/chemistry/EOS/independent physical reference admission or NLTE solve.'])
 (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k]for k in ['native_calls','mapping_max_relative_difference','mapping_passed','roundtrip_max_abs_A','summary_maxima']}));print(f'Read immutable native states/source; wrote private probe, logs,672summary rows and comparison/report to {a.output}')
 if not report['mapping_passed']:raise RuntimeError('Native coordinate arithmetic failed; no bound adjustment')

if __name__=='__main__':main()

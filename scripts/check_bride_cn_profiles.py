"""Native Voigt calls and external-state CN extinction unit diagnostic."""
import argparse,csv,json,os,shutil,subprocess,sys,time
from pathlib import Path
import numpy as np
from scipy.special import wofz
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
WRAPPER='''program probe
  implicit none
  integer count,i
  real a,v,newvoigt
  read(*,*) count
  if (count<1.or.count>20000) stop 2
  do i=1,count
    read(*,*) a,v
    write(*,'(I8,1P,3E25.16)') i,a,v,newvoigt(a,v)
  enddo
end program probe
'''


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_cn_diagnostic import external_state_extinction
    plan_path=ROOT/'config/experiments/rya1281_cn_profiles.json';plan=json.loads(plan_path.read_text());state_dir=require_codex_path(plan['source_state_key']);state_report=state_dir/'report.json';previous=json.loads(state_report.read_text())
    if previous['source_pin']!=plan['source_commit']or not previous['all_equation_checks_pass']:raise ValueError('Unqualified input state report')
    archived=ROOT/'data/results/rya1281/cn_coefficients_v1/CN/states.csv';states_path=state_dir/'CN/states.csv'
    if states_path.read_bytes()!=archived.read_bytes()or states_path.read_bytes()!=(state_dir/'all/states.csv').read_bytes():raise ValueError('Changed native state')
    states=np.genfromtxt(states_path,delimiter=',',names=True)
    if len(states)!=168:raise ValueError('Unexpected sampled CN state count')
    engine=require_codex_path('engines.ispec');source=engine/'synthesizer/turbospectrum/source/newvoigt.f';source_bytes=source.read_bytes()
    if source_bytes!=subprocess.check_output(['git','show',plan['source_commit']+':synthesizer/turbospectrum/source/newvoigt.f'],cwd=engine):raise ValueError('Changed native Voigt source')
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'newvoigt.f').write_bytes(source_bytes);(a.output/'probe.f90').write_text(WRAPPER)
    compiler=shutil.which('gfortran')
    if compiler is None:raise RuntimeError('Existing compiler required')
    program=a.output/'probe';command=[compiler,'-O3','-fno-automatic','-std=legacy','-o',str(program),'probe.f90','newvoigt.f']
    before=time.perf_counter();r=subprocess.run(command,cwd=a.output,text=True,capture_output=True,timeout=plan['compile_seconds']);(a.output/'compiler.log').write_text(r.stdout+r.stderr)
    if r.returncode:raise RuntimeError('Native Voigt probe compile failed')
    compile_seconds=time.perf_counter()-before
    offsets=np.asarray(plan['signed_offsets'],dtype=np.float32);pairs=[(row['voigt_a'],v)for row in states for v in offsets];sample_count=len(pairs)
    chosen=np.array([np.min(states['voigt_a']),np.median(states['voigt_a']),np.max(states['voigt_a'])],dtype=np.float32);integration=[]
    for order in plan['quadrature_orders']:
        nodes,weights=np.polynomial.legendre.leggauss(order);v=np.tan(np.pi*nodes/2).astype(np.float32);jac=np.pi/2/np.cos(np.pi*nodes/2)**2
        for damping in chosen:
            start=len(pairs);pairs.extend((damping,x)for x in v);integration.append((order,float(damping),start,len(pairs),weights*jac))
    text=str(len(pairs))+'\n'+''.join(f'{float(np.float32(d)):.17g} {float(v):.17g}\n'for d,v in pairs);(a.output/'input.txt').write_text(text)
    before=time.perf_counter();r=subprocess.run([str(program)],input=text,cwd=a.output,text=True,capture_output=True,timeout=plan['probe_seconds'],env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'));elapsed=time.perf_counter()-before
    (a.output/'stdout.txt').write_text(r.stdout);(a.output/'stderr.txt').write_text(r.stderr)
    if r.returncode:raise RuntimeError('Native Voigt probe failed')
    raw=np.array([[float(x)for x in line.split()]for line in r.stdout.splitlines()])
    if raw.shape!=(len(pairs),4)or not np.isfinite(raw).all()or np.any(raw[:,3]<0):raise ValueError('Invalid native Voigt output')
    np.testing.assert_array_equal(raw[:,0],np.arange(1,len(pairs)+1));expected=np.asarray(pairs,dtype=np.float32).astype(float);np.testing.assert_array_equal(raw[:,1:3],expected)
    native=raw[:sample_count,3].reshape(len(states),len(offsets)).T;reference=wofz(offsets[:,None].astype(float)+1j*states['voigt_a'][None,:]).real
    if not np.array_equal(native,native[::-1]):raise ValueError('Native signed-offset symmetry failed')
    difference=np.abs(native/reference-1);max_index=np.unravel_index(np.argmax(difference),difference.shape)
    # Direct native mass-opacity conversion (cm2/g * g/cm3 -> cm^-1 -> m^-1).
    direct=native*states['pre_voigt_amplitude_over_standard'][None,:]*states['standard_opacity_cm2_g'][None,:]*states['rho_g_cm3'][None,:]*100.
    # Separate SI kernel from explicit population per m3 and integrated strength.
    strength=np.pi*states['e_esu']**2/(states['electron_mass_g']*states['c_cm_s'])*1e-4
    diagnostic=external_state_extinction(offsets,states['voigt_a'],states['doppler_width_Hz'],states['lower_population_per_g']*states['rho_g_cm3']*1e6,states['gf'],states['stim'],strength)
    # Remove approximation difference only for the unit-accounting check; retain
    # raw native-versus-Faddeeva comparison independently, with no tuned correction.
    native_shape_si=diagnostic*(native/reference);unit_difference=float(np.max(np.abs(native_shape_si/direct-1)))
    bound=float(plan['unit_check_float32_eps_multiple']*np.finfo(np.float32).eps)
    area=[]
    for order,damping,start,end,weights in integration:
        H=raw[start:end,3];v=raw[start:end,2];a_native=raw[start:end,1]
        area.append(dict(order=order,damping_a=damping,native_H_area_over_sqrt_pi=float(np.sum(weights*H)/np.sqrt(np.pi)),faddeeva_H_area_over_sqrt_pi=float(np.sum(weights*wofz(v+1j*a_native).real)/np.sqrt(np.pi))))
    with (a.output/'profile_samples.csv').open('w',newline='')as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['state_index','offset_doppler','a','native_H','faddeeva_H','native_chi_m_inv','diagnostic_chi_m_inv'])
        w.writerows((j,float(offsets[i]),float(states['voigt_a'][j]),native[i,j],reference[i,j],direct[i,j],diagnostic[i,j])for j in range(len(states))for i in range(len(offsets)))
    if source.read_bytes()!=source_bytes or fingerprint(states_path)!=fingerprint(archived):raise ValueError('Protected source/state changed during probe')
    report=dict(schema='codex.bride_cn_profiles/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),kernel=fingerprint(ROOT/'pipeline/bride_cn_diagnostic.py'),plan=fingerprint(plan_path),source=fingerprint(source),state_report=fingerprint(state_report),states=fingerprint(states_path),wrapper=fingerprint(a.output/'probe.f90'),program=fingerprint(program),compiler_command=command,compile_seconds=compile_seconds,probe_seconds=elapsed,native_calls=len(pairs),sampled_states=len(states),offsets=offsets.tolist(),max_profile_relative_difference=float(difference[max_index]),max_profile_absolute_difference=float(np.max(np.abs(native-reference))),worst_state_index=int(max_index[1]),worst_offset=float(offsets[max_index[0]]),native_signed_symmetry_exact=True,profile_accuracy_admission_threshold=None,unit_accounting_max_relative_difference=unit_difference,unit_accounting_prospective_bound=bound,unit_accounting_passed=unit_difference<=bound,area_checks=area,
      limitations=['Explicit native damping and dimensionless offsets; no AIR/vacuum wavelength conversion or linear-versus-exact frequency mapping is validated.','External native lower populations, isotope weights and native constants retained; no chemical/primary provenance or physical opacity admission.','Unit check separates native profile approximation from SI bookkeeping; diagnostic uses uncorrected SciPy Faddeeva. No fitted profile correction applied.','Finite quadrature orders demonstrate only sampled numerical area behavior, not a physical error budget or complete far-wing validation.','No new RT spectrum, NLTE solve, background provider, per-program peak-memory measurement or abundance.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k]for k in ['native_calls','max_profile_relative_difference','max_profile_absolute_difference','worst_offset','unit_accounting_max_relative_difference','unit_accounting_passed','area_checks']}));print(f'Read unchanged native source and168states; wrote private probe, raw logs, samples and report to {a.output}')
    if not report['unit_accounting_passed']:raise RuntimeError('Unit-accounting check failed; no bound adjustment')

if __name__=='__main__':main()

"""Pinned external comparison engines; this registry confers no product admission."""
import hashlib,json,platform,subprocess
from pathlib import Path
import yaml
ROOT=Path(__file__).resolve().parents[1]


def fingerprint(path):
    path=Path(path);return dict(name=path.name,bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def inside(root,relative):
    p=Path(relative)
    if p.is_absolute()or'..'in p.parts:raise ValueError('Require an owned relative engine path')
    result=(root/p).resolve()
    if not result.is_relative_to(root.resolve()):raise ValueError('Engine path escapes checkout')
    return result


def load_engine(name,*,registry=ROOT/'config/reference_engines.yaml',resolver=None):
    if resolver is None:
        from config.constants import require_codex_path
        resolver=require_codex_path
    data=yaml.safe_load(Path(registry).read_text())
    if data['schema']!='codex.reference_engines/1'or name not in data['engines']:raise ValueError('Unknown reference engine/schema')
    spec=data['engines'][name]
    if spec['status']!='experimental-reference'or spec['production_admission']!='HOLD':raise ValueError('Unexpected reference-engine admission state')
    checkout=resolver(spec['checkout_resource']);app=inside(checkout,spec['application_rel']);exe=inside(app,spec['executable_rel']);build=resolver(spec['build_resource'])/'receipt.json';receipt=json.loads(build.read_text())
    actual=subprocess.check_output(['git','rev-parse','HEAD'],cwd=checkout,text=True).strip()
    dirty=subprocess.check_output(['git','status','--porcelain','--untracked-files=no'],cwd=checkout,text=True)
    if actual!=spec['source_commit']or dirty:raise ValueError('Changed/dirty reference-engine source')
    if receipt['source_commit']!=actual or receipt['tracked_changes']or fingerprint(exe)!=receipt['executable']:raise ValueError('Changed reference-engine executable/build')
    if name=='openmulti3d':
        lib=app/'HOP/standalone/build/libhop.so'
        if fingerprint(lib)!=receipt['libraries']['hop']:raise ValueError('Changed HOP library')
    else:raise ValueError('No executable adapter for this reference engine')
    return spec,app,exe,dict(registry=fingerprint(registry),build_receipt=fingerprint(build),executable=fingerprint(exe),source_commit=actual)


def load_plan_engine(plan,*,resolver=None):
    """Original engine by default; opt-in isolated profile derivative has its own spec."""
    spec,app,exe,original=load_engine(plan['engine'],resolver=resolver)
    variant=plan.get('profile_variant','original')
    if variant=='original':return spec,app,exe,original
    if variant!='accurate-kernel-control':raise ValueError('Unsupported reference profile variant')
    if resolver is None:
        from config.constants import require_codex_path
        resolver=require_codex_path
    registry=ROOT/'config/experiments/openmulti3d_profile_variant.yaml';derived=yaml.safe_load(registry.read_text())
    if derived['schema']!='codex.reference_profile_variant/1'or derived['variant']!=variant or derived['status']!='experimental-reference'or derived['production_admission']!='HOLD'or derived['source_commit']!=spec['source_commit']:raise ValueError('Invalid profile variant spec')
    checkout=resolver(derived['checkout_resource']);app=inside(checkout,spec['application_rel']);exe=inside(app,spec['executable_rel']);receipt_file=resolver(derived['build_resource'])/'receipt.json';receipt=json.loads(receipt_file.read_text());patch=inside(ROOT,derived['patch']);diff=subprocess.check_output(['git','diff','--binary','HEAD'],cwd=checkout)
    actual=subprocess.check_output(['git','rev-parse','HEAD'],cwd=checkout,text=True).strip()
    if actual!=spec['source_commit']or diff!=patch.read_bytes():raise ValueError('Changed profile derivative source')
    if not receipt['completed']or receipt['source_commit']!=actual or receipt['original_reference']!=original or receipt['patch']!=fingerprint(patch)or receipt['tracked_diff_sha256']!=hashlib.sha256(diff).hexdigest():raise ValueError('Unqualified profile derivative build')
    if fingerprint(exe)!=receipt['executable']or fingerprint(app/'HOP/standalone/build/libhop.so')!=receipt['libraries']['hop']:raise ValueError('Changed profile derivative executable/library')
    return spec,app,exe,dict(variant=variant,variant_spec=fingerprint(registry),build_receipt=fingerprint(receipt_file),executable=fingerprint(exe),source_commit=actual,patch=fingerprint(patch),original_reference=original)


MAX_REFERENCE_OUTPUT_POINTS=20000


def output_grid_vacuum_A(grid):
    """Build explicit fine AIR-A or vacuum-nm windows; return native vacuum A.

    This is a coordinate/point-budget contract, not a physical coverage claim.
    """
    import numpy as np
    frame=grid.get('frame')
    key='step_A'if frame=='air_A'else'step_nm'if frame=='vacuum_nm'else None
    if key is None or set(grid)!={'frame',key,'windows'}:raise ValueError('Require unambiguous AIR-A or vacuum-nm output grid')
    step=grid[key];maximum=.01 if frame=='air_A'else .001
    if isinstance(step,bool)or not isinstance(step,(int,float))or not np.isfinite(step)or not 0<step<=maximum:raise ValueError('Require explicit fine positive output sampling')
    windows=grid['windows']
    if not isinstance(windows,list)or not windows:raise ValueError('Require nonempty output windows')
    previous=0.;count_bound=0.
    for pair in windows:
        if not isinstance(pair,(list,tuple))or len(pair)!=2:raise ValueError('Invalid output window')
        lo,hi=pair
        if any(isinstance(v,bool)or not isinstance(v,(int,float))or not np.isfinite(v)for v in pair)or not(2000. if frame=='air_A'else 0.)<lo<hi or lo<=previous:raise ValueError('Invalid/overlapping output windows or AIR below200nm')
        previous=hi;count_bound+=(hi-lo)/step+2
    if count_bound>MAX_REFERENCE_OUTPUT_POINTS:raise ValueError('Output grid exceeds bounded reference scope')
    wave=np.concatenate([np.linspace(lo,hi,int(np.ceil((hi-lo)/step))+1)for lo,hi in windows])
    if frame=='air_A':
        from pipeline.wavelength_util import air_to_vac
        wave=air_to_vac(wave)
    else:wave=wave*10.
    if not np.isfinite(wave).all()or np.any(np.diff(wave)<=0)or len(wave)>MAX_REFERENCE_OUTPUT_POINTS:raise ValueError('Invalid native vacuum output axis')
    return wave,dict(requested_frame=frame,native_frame='vacuum_A',points=len(wave),windows=windows,conversion='canonical AIR->vacuum once'if frame=='air_A'else'nm->A scale10 only; NO air conversion',scientific_coverage='Numerical requested samples only; no full-band/model/instrument admission')


def validate_plan(plan):
    if plan['schema']!='codex.reference_engine_run/1'or plan['engine']!='openmulti3d':raise ValueError('Unsupported run schema/engine')
    if plan['mode']not in ['1D-LTE','1D-NLTE']:raise ValueError('Only explicit1D reference runs are implemented')
    if plan.get('profile_variant','original')not in ['original','accurate-kernel-control']:raise ValueError('Unsupported reference profile variant')
    c=plan['controls']
    if c['threads']!=1 or not 0<c['seconds']<=600 or not 0<c['maxiter']<=1500:raise ValueError('Reference run exceeds authorized resource bounds')
    if not 0<c['rtol']<=1e-6:raise ValueError('Invalid reference convergence threshold')
    if set(plan['atmosphere'])!={'resource','filename','receipt'}:raise ValueError('Incomplete atmosphere identity')
    if set(plan['atom'])not in [{'resource','receipt'},{'supplied_file','receipt'}]:raise ValueError('Require exactly one explicit atom identity')
    for axis in ['ntheta','nphi']:
        if not isinstance(c[axis],int)or c[axis]<=0:raise ValueError('Invalid ray quadrature')
    if c['angle_quad']not in ['lobatto','radau','gauss']:raise ValueError('Unrecognized quadrature declaration')
    if c['ntheta']>(10 if c['angle_quad']=='radau' else 6):raise ValueError('Quadrature exceeds pinned native angular support')
    if not isinstance(c.get('bound_free_in_output',False),bool):raise ValueError('Bound-free output control must be explicit boolean')
    if c.get('molecular_equilibrium','native')not in ['native','suppressed-above-1K-control']:raise ValueError('Unsupported molecular EOS control (tmolim0 crashes pinned native initialization)')
    if c.get('save_species','')not in ['', 'O I,O II,CO,OH,H2O,H I,H II']:raise ValueError('Unsupported diagnostic species export')
    if c.get('output_angle_quad','none')not in ['none','disk_center']:raise ValueError('Unsupported reference output quadrature')
    if c.get('output_solver','quadratic')not in ['linear','quadratic']:raise ValueError('Unsupported reference output formal solver')
    if not isinstance(c.get('disable_stark_control',False),bool):raise ValueError('Require explicit Stark diagnostic boolean')
    if not isinstance(c.get('disable_vdw_control',False),bool):raise ValueError('Require explicit VDW diagnostic boolean')
    if 'decouple_continuum'in c and not isinstance(c['decouple_continuum'],bool):raise ValueError('Require explicit radiation-decoupling boolean')
    if 'output_grid' in plan:output_grid_vacuum_A(plan['output_grid'])
    if not isinstance(plan['expected_element'],str)or not plan['expected_element'].isalpha():raise ValueError('Explicit element required')


def quote_nml(value):
    s=str(value)
    if any(c in s for c in ["'",'\n','\r']):raise ValueError('Unsafe Fortran string')
    return "'"+s+"'"


def prepare_run(plan,output,app,*,resolver):
    """Resolve and receipt inputs before allocating an exclusively owned directory."""
    validate_plan(plan);at=plan['atmosphere'];atmos=inside(resolver(at['resource']),at['filename']);atom_spec=plan['atom']
    atom=resolver(atom_spec['resource'])if 'resource'in atom_spec else inside(app/'input/atoms',atom_spec['supplied_file'])
    for path,expected in [(atmos,at['receipt']),(atom,atom_spec['receipt'])]:
        if fingerprint(path)!=expected:raise ValueError('Changed declared physical input '+str(path))
    entries=[l.strip()for l in atom.read_text().splitlines()if l.strip()and not l.lstrip().startswith('*')]
    if entries[0]!=plan['expected_element']:raise ValueError('Requested element differs from loaded atom')
    # Classical MULTI header identity only; the native backend owns conversion/rates.
    abundance,mass=map(float,entries[1].split());levels,lines,continua,fixed=map(int,entries[2].split())
    if abundance!=plan['abundance']:raise ValueError('Declared/source atom abundance mismatch')
    c=plan['controls'];rtol=c['rtol']if plan['mode']=='1D-NLTE'else 0.;maxiter=c['maxiter']if rtol else 0
    local_limit=rtol/10 if c.get('freeze_depths',True)else 0.
    radiation_control=(', decouple_continuum='+('T'if c['decouple_continuum']else'F'))if 'decouple_continuum'in c else''
    nml=f'''&io_params datadir='native', gb_step=6., do_trace=F /
&timer_params sec_per_report=1d8, verbose=1 /
&atmos_params atmos_format='Marcs', atmos_file={quote_nml(atmos)}, use_abnd=T, use_ne=T, use_rho=T, prune_top=F, save_atmos=T, save_eos=T /
&atom_params atom_file={quote_nml(atom)}, use_atom_abnd=T, qnorm={c['qnorm']}, nq_fac={c['nq_fac']} /
&nlte_params convlim={rtol:.12e}, convmax={local_limit:.12e}, premature_finish=F{radiation_control} /
&rt_iteration_params angle_quad={quote_nml(c['angle_quad'])}, ntheta={c['ntheta']}, nphi={c['nphi']}, maxiter={maxiter}, include_absmet=F, include_linelist=F /
&rt_output_params solver={quote_nml(c.get('output_solver','quadratic'))}, angle_quad={quote_nml(c.get('output_angle_quad','none'))}, include_absmet=F, include_linelist=F /
&output_params save_populations=T, save_rates=1, save_rate_matrix=T, save_net_rates=T, save_Jnu=T, bound_free_in_output={'T'if c.get('bound_free_in_output',False)else'F'} /
&run_params mem_limit_gb=2. /
&m3d_debug_params verbose=1, check_NaN=T /
'''
    composition=[]
    if c.get('disable_stark_control',False):composition.append('disable_stk=T')
    if c.get('disable_vdw_control',False):composition.append('disable_vdw=T')
    if c.get('save_species')or c.get('molecular_equilibrium','native')!='native':
        composition.append('save_species='+quote_nml(c.get('save_species','')))
        if c.get('molecular_equilibrium','native')=='suppressed-above-1K-control':composition.append('tmolim=1')
    if composition:nml+='&composition_params '+', '.join(composition)+' /\n'
    output=Path(output)
    if output.resolve().is_relative_to(app.resolve()):raise ValueError('Run output must be outside engine source')
    output.mkdir(parents=True,exist_ok=False);(output/'input').symlink_to(app/'input',target_is_directory=True);(output/'input.nml').write_text(nml)
    output_grid_audit=None
    if 'output_grid' in plan:
        import numpy as np
        wave,output_grid_audit=output_grid_vacuum_A(plan['output_grid'])
        np.savetxt(output/'output_wavelengths_vac_A.txt',wave,fmt='%.12f')
        nml+="&spectrum_params lam_file='output_wavelengths_vac_A.txt', in_air=F /\n"
        (output/'input.nml').write_text(nml)
    return dict(output_grid=output_grid_audit,atmosphere=fingerprint(atmos),atom=fingerprint(atom),element=entries[0],atom_header=dict(abundance=abundance,mass=mass,levels=levels,lines=lines,continua=continua,fixed=fixed),namelist=fingerprint(output/'input.nml'),controls=c,physical_limitations=['Native HOP EOS/continuous opacity and collision interpolation differ from Bride; explicit matched-input audit is required.','Line blanketing and background bound-bound opacity explicitly off in this first reference plan.','Output quadrature is explicit; none skips synthesis and disk_center requests vertical intensity, not flux or an abundance fit.','Source atom is not silently replaced by the supplied variant; source paths and hashes are checked.'])

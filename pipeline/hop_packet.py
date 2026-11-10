"""Explicit immutable HOP EOS/continuum packets for experimental fixed transfer.

No source-data fallback, active-atom rates, population solve, or product admission.
The HOP instance is closed after each packet: no live global EOS is shared by
Lightweaver contexts. Wavelengths enter in vacuum nm; exported opacity is SI.
"""
import ctypes as ct
import hashlib,json,os,re,subprocess,threading
from pathlib import Path
import numpy as np
from pipeline.hop_abi import bind_header
from pipeline.reference_engines import load_engine,fingerprint
ROOT=Path(__file__).resolve().parents[1]
_LOCK=threading.Lock()


def checked_runtime(provider_variant='trace'):
    from config.constants import require_codex_path
    spec,reference_app,_,reference=load_engine('openmulti3d')
    if provider_variant not in ('trace','collider'):raise ValueError('Unknown explicit HOP provider variant')
    prefix='engines.hop_'+provider_variant+'_provider_'
    source=require_codex_path(prefix+'source');receipt_file=require_codex_path(prefix+'build')/'receipt.json';r=json.loads(receipt_file.read_text());patch=ROOT/('native/hop_trace_provider/trace_continua.patch'if provider_variant=='trace'else'native/hop_collider_provider/trace_and_ground.patch')
    if not r['completed']or r['upstream_commit']!=spec['source_commit']or r['source']!=str(source.resolve())or r['patch']!=fingerprint(patch):raise ValueError('Changed/invalid derived HOP recipe')
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip();diff=subprocess.check_output(['git','diff','--binary'],cwd=source)
    if head!=r['upstream_commit']or hashlib.sha256(diff).hexdigest()!=r['tracked_diff_sha256']:raise ValueError('Changed derived HOP source')
    app=source/spec['application_rel'];library=app/'HOP/standalone/build/libhop.so';header=app/'HOP/hopacity/libhop.h'
    if fingerprint(library)!=r['library']or fingerprint(header)!=r['header']:raise ValueError('Changed derived HOP library/header')
    paths=[reference_app/'input'/name for name in ['absdat/TS_absdat.txt','abund/abund_magg22.txt','isotopes/isotopes_blp.txt','molec/molec_data_BC.txt']]
    auxiliary=dict(vcs_data_dir=reference_app/'input/vcs',anstee_file=reference_app/'input/Anstee-1802-2014.txt')
    extras=[auxiliary['anstee_file'],*(p for p in auxiliary['vcs_data_dir'].rglob('*')if p.is_file())]
    options_file=ROOT/'config/experiments/rya1371_hop_native_options.json';recipe=json.loads(options_file.read_text())
    if recipe['schema']!='codex.hop_fixed_reference_options/1':raise ValueError('Unknown native options recipe')
    checks={str(path):fingerprint(path)for path in [library,header,receipt_file,patch,app/'HOP/const_mod.f90',ROOT/'pipeline/hop_abi.py',options_file,*paths,*extras]}
    return bind_header(ct.CDLL(str(library)),header.read_text()),paths,checks,reference,auxiliary,recipe['options']


def native_planck(wave_nm,temperature,source):
    """Native historical HOP constants and Wien branch, converted to SI I_nu."""
    text=Path(source).read_text();values={}
    for key in ['hh','bk','cc']:
        match=re.search(r'\b'+key+r'\s*=\s*([0-9.+\-dDeE]+)',text,re.I)
        if match is None:raise ValueError('Missing pinned native physical constant')
        values[key]=float(match[1].lower().replace('d','e'))
    nu=values['cc']/(np.asarray(wave_nm)[:,None]*1e-7);T=np.asarray(temperature)[None,:];x=values['hh']*nu/(values['bk']*T);occupation=np.empty_like(x);small=x<80
    occupation[small]=1/np.expm1(x[small]);occupation[~small]=np.exp(-x[~small])
    return 2*values['hh']*nu**3/values['cc']**2*occupation*1e-3,values


def compute_packet(T_K,ne_m3,rho_kg_m3,wave_vac_nm,abundances_dex,atom,*,vmic_m_s,provider_variant='trace'):
    """Fresh prescribed-T/ne/rho EOS and selected background, fixed atmosphere.

    The returned state and coefficients are in the caller's explicit depth order.
    Molecule populations are selected diagnostics of HOP's complete internal EOS.
    Active O level populations/nstar are owned separately by a bound model packet.
    """
    if provider_variant not in ('trace','collider'):raise ValueError('Unknown explicit HOP provider variant')
    if os.environ.get('OMP_NUM_THREADS')!='1':raise ValueError('Require explicit one-thread HOP environment')
    T=np.array(T_K,dtype=float,copy=True,order='C');ne=np.array(ne_m3,dtype=float,copy=True,order='C')*1e-6;rho=np.array(rho_kg_m3,dtype=float,copy=True,order='C')*1e-3;wave=np.array(wave_vac_nm,dtype=float,copy=True)
    if T.ndim!=1 or not 1<=len(T)<=256 or ne.shape!=T.shape or rho.shape!=T.shape or any(not np.isfinite(v).all()or np.any(v<=0)for v in [T,ne,rho]):raise ValueError('Invalid complete thermodynamic binding')
    if wave.ndim!=1 or not 1<=len(wave)<=20000 or not np.isfinite(wave).all()or np.any(wave<=0)or np.any(np.diff(wave)<=0):raise ValueError('Invalid bounded vacuum-nm packet grid')
    abundance=np.ascontiguousarray(np.asarray(abundances_dex,np.float32).astype(float))
    if abundance.shape!=(92,)or abundance[0]!=12 or not np.isfinite(abundance).all():raise ValueError('Require explicit full MARCS abundance set')
    if atom.symbol!='O'or not np.isfinite(vmic_m_s)or vmic_m_s<0:raise ValueError('Only explicitly qualified Oxygen packet implemented')
    from pipeline.multi_atom import CM1_TO_EV
    stages=np.ascontiguousarray([atom.level(c.i).stage for c in atom.continua],dtype=np.int32)
    exc=np.ascontiguousarray([(atom.level(c.i).E_cm-min(l.E_cm for l in atom.levels if l.stage==atom.level(c.i).stage))*CM1_TO_EV for c in atom.continua]);vac=np.ascontiguousarray(wave*10)
    with _LOCK:
        api,files,checks,reference,auxiliary,options=checked_runtime(provider_variant);ptr=lambda x:x.ctypes.data_as(ct.POINTER(ct.c_double));fptr=lambda x:x.ctypes.data_as(ct.POINTER(ct.c_float))
        def call(name,*args):
            status=api[name](*args)
            if status:
                buf=ct.create_string_buffer(2048);api['hop_last_error'](buf,2048);raise RuntimeError(name+': '+buf.value.decode())
        initialized=False
        try:
            if api['hop_get_n']()!=0:raise ValueError('Another live native EOS state is present')
            for key,value in options.items():
                encoded=('true'if value else'false')if isinstance(value,bool)else str(value)
                call('hop_set_option',key.encode(),encoded.encode())
            call('hop_set_external_abundances',len(abundance),ptr(abundance));call('hop_set_option',b'abundances',f'{atom.symbol}={atom.abundance:.17g}'.encode())
            for key,path in auxiliary.items():call('hop_set_option',key.encode(),str(path).encode())
            call('hop_init',*[str(v).encode()for v in files],0.,-99.,0.,0.);initialized=True
            count=ct.c_int(-1);call('hop_exclude_model_continua',atom.symbol.encode(),atom.weight,len(exc),stages.ctypes.data_as(ct.POINTER(ct.c_int)),ptr(exc),ct.byref(count))
            if count.value!=1:raise ValueError('Unexpected qualified O background ownership count')
            names=['O I','O II','H I','H II','H -','He I','CO','OH','H2O']
            for name in names:call('hop_keep_species',name.encode())
            call('hop_load_state',len(T),ptr(T),ptr(ne),ptr(rho),float(vmic_m_s/1000),1,1)
            state={}
            for name in ['temp','ne','rho']:
                v=np.empty_like(T);call('hop_get_field',name.encode(),len(T),ptr(v));state[name]=v
            for expected,name in [(T,'temp'),(np.asarray(ne_m3)*1e-6,'ne'),(rho,'rho')]:np.testing.assert_array_equal(state[name],expected)
            grph=api['hop_get_grph']()
            if not np.isfinite(grph)or grph<=0:raise ValueError('Invalid native mass per H nucleus')
            for name in names:
                v=np.empty_like(T);call('hop_get_species',name.encode(),len(T),ptr(v));state[name]=v
            if provider_variant=='collider':
                ground=np.empty_like(T);occupation=np.empty_like(T);call('hop_get_h_ground',len(T),ptr(ground),ptr(occupation))
                # Native init_H_levels keeps w outside level%n. Collision nh1 is raw;
                # only occupation*raw is bounded by total HI. Never renormalize raw.
                if not np.isfinite(ground).all()or np.any(ground<0)or not np.isfinite(occupation).all()or np.any(occupation<0)or np.any(occupation>1)or np.any(occupation*ground>state['H I']):raise ValueError('Invalid native occupation-aware H ground field')
                state['H_ground']=ground;names.append('H_ground')
            chi=np.zeros((len(T),len(wave)),np.float32);scat=np.zeros_like(chi);independent=np.zeros_like(chi)
            call('hop_opacity',len(wave),ptr(vac),fptr(chi),fptr(scat),0,0,100.,0);call('hop_scattering',len(wave),ptr(vac),fptr(independent));np.testing.assert_array_equal(scat,independent)
            if any(not np.isfinite(v).all()or np.any(v<0)for v in [chi,scat,*state.values()])or np.any(scat>chi):raise ValueError('Invalid native EOS/opacity packet')
        finally:
            if initialized:call('hop_free')
    const_file=next(k for k in checks if k.endswith('/HOP/const_mod.f90'));B,constants=native_planck(wave,T,const_file);chi=chi.T.astype(float)*100;scat=scat.T.astype(float)*100;eta=(chi-scat)*B
    if not np.isfinite(eta).all():raise ValueError('Invalid native LTE emissivity')
    for path,receipt in checks.items():
        if fingerprint(path)!=receipt:raise ValueError('HOP packet source changed')
    arrays=dict(temperature_K=T.copy(),ne_m3=np.asarray(ne_m3).copy(),rho_kg_m3=np.asarray(rho_kg_m3).copy(),nHTot_m3=rho/grph*1e6,wave_vac_nm=wave.copy(),chi_m_inv=chi,scat_m_inv=scat,eta_SI=eta,Bnu_SI=B,**{name.replace(' ','_')+'_m3':state[name]*1e6 for name in names})
    if provider_variant=='collider':arrays['H_ground_occupation_probability']=occupation.copy()
    for value in arrays.values():value.flags.writeable=False
    meta=dict(reference=reference,inputs=checks,module_options_reset=options,source_constants_cgs=constants,atom_trace_selection=dict(symbol=atom.symbol,abundance=atom.abundance,weight=atom.weight,stages=stages.tolist(),excitation_eV=exc.tolist()),excluded_O_records=count.value,frame='vacuum_nm',opacity='total extinction m-1 including scattering exactly once',source='native HOP historical-constant LTE I_nu; scattering source added by RT backend',molecular_equilibrium='native HOP, LTE; no molecular line-list opacity',scientific_admission='HOLD')
    return arrays,meta

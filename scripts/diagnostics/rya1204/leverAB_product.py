"""RYA-1204 Lever A, measured alone: the RYA-759 near-UV Fe product run with jonabs
component 17 REPLACED by the resonance-averaged Bautista/NORAD cross-section on a 20 A
grid (jonabs_LEVER_A.dat). Molecules stay OFF, so this is Lever A by itself.

Nothing is scaled to a target: the values are Bautista's own, in jonabs's measured
convention. Whatever it does to A(Fe) is the answer.
"""
import os, sys, glob, shutil, subprocess, tempfile
import numpy as np
SP=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,'/Users/ryanschmitt/codex/ispec'); sys.path.insert(0,'/Users/ryanschmitt/codex/rya1195')
sys.path.insert(0,'/Users/ryanschmitt/codex/rya1195/scripts')
TS='/Users/ryanschmitt/codex/ispec/synthesizer/turbospectrum'
TABLE=os.path.join(SP, os.environ.get('TABLE','jonabs_LEVER_A.dat'))
import ispec, ispec.synth.turbospectrum as tsmod
# A+B: also swap in the near-UV molecular lists and turn use_molecules on.
MOL=os.path.join(SP,'molecules')
_osym=os.symlink
def _sym(src,dst,*a,**k):
    if str(dst).endswith('/molecules'): src=MOL
    return _osym(src,dst,*a,**k)
os.symlink=_sym
_gen=ispec.generate_spectrum
def _g(*a,**k):
    k['use_molecules']=True
    return _gen(*a,**k)
ispec.generate_spectrum=_g
import pipeline.abundances_derive as _ad
_ad.ispec.generate_spectrum=_g
def patched(atmosphere_layers_file, abundances, MH, microturbulence_vel,
            wave_base, wave_top, wave_step, verbose=0, opacities_filename=None,
            tmp_dir=None, is_marcs_model=True):
    if opacities_filename is None:
        f=tempfile.NamedTemporaryFile(mode='wt',delete=False,dir=tmp_dir); f.close()
        opacities_filename=f.name
    run=tempfile.mkdtemp(dir=tmp_dir); os.makedirs(os.path.join(run,'DATA'))
    for p in glob.glob(os.path.join(TS,'DATA','*')):
        os.symlink(p, os.path.join(run,'DATA',os.path.basename(p)))
    os.remove(os.path.join(run,'DATA','jonabs_vac_v19.2.dat'))
    shutil.copy(TABLE, os.path.join(run,'DATA','jonabs_vac_v19.2.dat'))
    inp =f"'LAMBDA_MIN:'  '{wave_base*10.}'\n'LAMBDA_MAX:'  '{wave_top*10.}'\n"
    inp+=f"'LAMBDA_STEP:' '{wave_step*10.}'\n'MODELINPUT:' '{atmosphere_layers_file}'\n"
    inp+=f"'MARCS-FILE    :' '{'.true.' if is_marcs_model else '.false.'}'\n"
    inp+=f"'MODELOPAC:' '{opacities_filename}'\n'METALLICITY:'    '0.00'\n"
    inp+="'ALPHA/Fe   :'    '0.00'\n'HELIUM     :'    '0.00'\n"
    inp+="'R-PROCESS  :'    '0.00'\n'S-PROCESS  :'    '0.00'\n"
    at=abundances[abundances['code']<=92]
    inp+=f"'INDIVIDUAL ABUNDANCES:'   '{len(at)}'\n"
    for a in at: inp+="%i  %.2f\n"%(a['code'],12.036+a['Abund'])
    inp+=f"'XIFIX:' 'T'\n{microturbulence_vel}\n"
    subprocess.run([os.path.join(TS,'bin','babsma_lu')],input=inp.encode(),cwd=run,
                   capture_output=True,timeout=1800)
    shutil.rmtree(run,ignore_errors=True)
    return opacities_filename
tsmod.calculate_opacities=patched
sys.argv=['rya759_nearuv_fe_product.py','--limit',os.environ.get('LIMIT','40'),
          '--tag','rya1204_AB','--out',os.path.join(SP,'nearuv_AB.json')]
import runpy
runpy.run_path('/Users/ryanschmitt/codex/rya1195/scripts/rya759_nearuv_fe_product.py',
               run_name='__main__')

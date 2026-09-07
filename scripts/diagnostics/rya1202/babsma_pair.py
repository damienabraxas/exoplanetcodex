"""RYA-1202 — measure what Fe I bound-free photoionization contributes to the near-UV
continuous opacity, by running babsma twice through the REAL code path with only the
Fe I bf cross-sections changed.

The control is the point: `stock` uses a byte-identical copy of the shipped table, so a
difference between `stock` and the production default would mean the harness itself moved
something, and a difference between `stock` and `fezero` is attributable to Fe I bf alone.
"""
import os, subprocess, sys, shutil, glob
sys.path.insert(0,'/Users/ryanschmitt/codex/ispec')
sys.path.insert(0,'/Users/ryanschmitt/codex/rya1195')
import ispec
from pipeline.abundances_derive import _ISPEC_SOLAR_ABUND_FILE

SP  = os.path.dirname(os.path.abspath(__file__))
TS  = '/Users/ryanschmitt/codex/ispec/synthesizer/turbospectrum'
LO, HI, STEP = 3000.0, 3800.0, 2.0          # the near-UV band policy window
TEFF, LOGG, MH, VT = 5772.0, 4.44, 0.0, 1.0

ab = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
atoms = ab[ab['code'] <= 92]
assert len(atoms) == 92, len(atoms)

def build_input(opac):
    s  = f"'LAMBDA_MIN:'  '{LO}'\n'LAMBDA_MAX:'  '{HI}'\n'LAMBDA_STEP:' '{STEP}'\n"
    s += f"'MODELINPUT:' '{SP}/solar_marcs.mod'\n'MARCS-FILE    :' '.true.'\n"
    s += f"'MODELOPAC:' '{opac}'\n"
    s += "'METALLICITY:'    '0.00'\n'ALPHA/Fe   :'    '0.00'\n'HELIUM     :'    '0.00'\n"
    s += "'R-PROCESS  :'    '0.00'\n'S-PROCESS  :'    '0.00'\n"
    s += f"'INDIVIDUAL ABUNDANCES:'   '{len(atoms)}'\n"
    for a in atoms:                              # iSpec's SPECTRUM->TS offset
        s += "%i  %.2f\n" % (a['code'], 12.036 + a['Abund'])
    s += f"'XIFIX:' 'T'\n{VT}\n"
    return s

def run(tag, table):
    run_dir = os.path.join(SP, f'run_{tag}')
    shutil.rmtree(run_dir, ignore_errors=True)
    os.makedirs(os.path.join(run_dir, 'DATA'))
    for f in glob.glob(os.path.join(TS, 'DATA', '*')):
        os.symlink(f, os.path.join(run_dir, 'DATA', os.path.basename(f)))
    os.remove(os.path.join(run_dir, 'DATA', 'jonabs_vac_v19.2.dat'))
    shutil.copy(table, os.path.join(run_dir, 'DATA', 'jonabs_vac_v19.2.dat'))
    opac = os.path.join(run_dir, 'opac.out')
    p = subprocess.run([os.path.join(TS,'bin','babsma_lu')],
                       input=build_input(opac).encode(), cwd=run_dir,
                       capture_output=True, timeout=900)
    open(os.path.join(run_dir,'babsma.log'),'wb').write(p.stdout + b'\n--STDERR--\n' + p.stderr)
    size = os.path.getsize(opac) if os.path.exists(opac) else 0
    print(f'{tag:8s} rc={p.returncode} opac={size}B')
    return opac if size else None

a = run('stock',  os.path.join(SP,'jonabs_STOCK.dat'))
b = run('fezero', os.path.join(SP,'jonabs_FeI_ZEROED.dat'))
print('both produced output:', bool(a and b))

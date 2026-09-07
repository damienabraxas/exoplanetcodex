"""RYA-1202 item 3 — the PAYOFF, measured rather than argued.

Synthesise a near-UV Fe I window through the production path twice: once with the shipped
continuous-opacity table, once with the Fe I bound-free cross-sections zeroed. The flux
difference is what the whole Fe I photoionization component is worth in this band, which
BOUNDS what refining it (uncut resonances, a finer grid) could ever recover.

The swap is a monkeypatch on iSpec's own `calculate_opacities`, so everything else -- model,
abundances, line list, NLTE, wave grid -- is the production call, unchanged.
"""
import os, sys, glob, shutil, subprocess, tempfile
import numpy as np
sys.path.insert(0,'/Users/ryanschmitt/codex/ispec')
sys.path.insert(0,'/Users/ryanschmitt/codex/rya1195')
import ispec
import ispec.synth.turbospectrum as tsmod

SP  = os.path.dirname(os.path.abspath(__file__))
TS  = '/Users/ryanschmitt/codex/ispec/synthesizer/turbospectrum'
TABLE = {'name': None}

_orig = tsmod.calculate_opacities
def patched(atmosphere_layers_file, abundances, MH, microturbulence_vel,
            wave_base, wave_top, wave_step, verbose=0, opacities_filename=None,
            tmp_dir=None, is_marcs_model=True):
    """Same babsma call iSpec makes, but with our chosen jonabs table in DATA."""
    if opacities_filename is None:
        f = tempfile.NamedTemporaryFile(mode='wt', delete=False, dir=tmp_dir); f.close()
        opacities_filename = f.name
    run = tempfile.mkdtemp(dir=tmp_dir)
    os.makedirs(os.path.join(run,'DATA'))
    for p in glob.glob(os.path.join(TS,'DATA','*')):
        os.symlink(p, os.path.join(run,'DATA',os.path.basename(p)))
    os.remove(os.path.join(run,'DATA','jonabs_vac_v19.2.dat'))
    shutil.copy(TABLE['name'], os.path.join(run,'DATA','jonabs_vac_v19.2.dat'))
    inp  = f"'LAMBDA_MIN:'  '{wave_base*10.}'\n'LAMBDA_MAX:'  '{wave_top*10.}'\n"
    inp += f"'LAMBDA_STEP:' '{wave_step*10.}'\n'MODELINPUT:' '{atmosphere_layers_file}'\n"
    inp += f"'MARCS-FILE    :' '{'.true.' if is_marcs_model else '.false.'}'\n"
    inp += f"'MODELOPAC:' '{opacities_filename}'\n'METALLICITY:'    '0.00'\n"
    inp += "'ALPHA/Fe   :'    '0.00'\n'HELIUM     :'    '0.00'\n"
    inp += "'R-PROCESS  :'    '0.00'\n'S-PROCESS  :'    '0.00'\n"
    at = abundances[abundances['code'] <= 92]
    inp += f"'INDIVIDUAL ABUNDANCES:'   '{len(at)}'\n"
    for a in at: inp += "%i  %.2f\n" % (a['code'], 12.036 + a['Abund'])
    inp += f"'XIFIX:' 'T'\n{microturbulence_vel}\n"
    subprocess.run([os.path.join(TS,'bin','babsma_lu')], input=inp.encode(),
                   cwd=run, capture_output=True, timeout=1800)
    shutil.rmtree(run, ignore_errors=True)
    return opacities_filename
tsmod.calculate_opacities = patched

from pipeline.abundances_derive import _ISPEC_SOLAR_ABUND_FILE, _SYNTH_ISOTOPE_FILE
TEFF, LOGG, MH, VT = 5772.0, 4.44, 0.0, 1.0
LO_nm, HI_nm = 335.0, 337.0                    # a 20 A near-UV window, mid-band
pack = ispec.load_modeled_layers_pack('/Users/ryanschmitt/codex/ispec/input/atmospheres/MARCS.GES/')
atm  = ispec.interpolate_atmosphere_layers(pack, {'teff':TEFF,'logg':LOGG,'MH':MH,'alpha':0.0}, code='turbospectrum')
ab   = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
iso  = ispec.read_isotope_data(_SYNTH_ISOTOPE_FILE)
ll   = ispec.read_atomic_linelist(os.environ['NEARUV_LL'])
ll   = ll[(ll['wave_nm'] >= LO_nm-1) & (ll['wave_nm'] <= HI_nm+1)]
print(f'window {LO_nm*10:.0f}-{HI_nm*10:.0f} A, {len(ll)} lines in list')
wave = np.arange(LO_nm, HI_nm, 0.0002)

out = {}
for tag, table in (('stock', f'{SP}/jonabs_STOCK.dat'), ('fezero', f'{SP}/jonabs_FeI_ZEROED.dat')):
    TABLE['name'] = table
    f = ispec.generate_spectrum(wave, atm, TEFF, LOGG, MH, 0.0, ll, iso, ab, [],
                                microturbulence_vel=VT, macroturbulence=0.0, vsini=0.0,
                                R=0, verbose=0, code='turbospectrum', tmp_dir=SP)
    out[tag] = np.asarray(f['flux'] if hasattr(f,'dtype') and f.dtype.names and 'flux' in f.dtype.names else f)
    print(f'  {tag:7s} mean flux {out[tag].mean():.5f}  min {out[tag].min():.5f}')

np.save(f'{SP}/flux_stock.npy', out['stock']); np.save(f'{SP}/flux_fezero.npy', out['fezero'])
np.save(f'{SP}/wave.npy', wave)
d = out['fezero'] - out['stock']
print(f'\nflux change from REMOVING Fe I bf: mean {d.mean():+.5f}  max {np.abs(d).max():.5f}')

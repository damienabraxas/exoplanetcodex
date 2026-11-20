"""Run MULTI 2.3 (our build) on a Gerber zladder record's MARCS model. argv: record-glob abund_mode(solar|model) outdir [IOSMET]
Writes triplet b vs depth and compares with the Gerber record."""
import sys, os, json, glob, subprocess, shutil, re
sys.path.insert(0, '.')
import numpy as np
from pathlib import Path
from config.constants import require_codex_path
from pipeline.bride_reference import read_marcs_structure, marcs_number_densities

M = Path('~/codex/bride_inputs/rya1282/mul23_build').expanduser()
PM = Path('~/codex/bride_inputs/rya1281/pymulti-source-v2/pymulti').expanduser()
sys.path.insert(0, str(PM.parent));
import importlib.util
spec = importlib.util.spec_from_file_location('pmm', PM / 'multi.py')
src = (PM / 'multi.py').read_text()
ns = {}
exec('import numpy as np\nfrom scipy.io import FortranFile\n' + src[src.index('def mulrd'):src.index('def initiate_dscale')], ns)
mulrd = ns['mulrd']

glob_, mode, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
iosmet = sys.argv[4] if len(sys.argv) > 4 else '0'
osscale = sys.argv[5] if len(sys.argv) > 5 else None   # 'feh' -> 10**[Fe/H]
rec_path = sorted(glob.glob(os.path.expanduser('~/codex/bride_inputs/rya1282/bench_gerber2023_O_MARCS_zladder/' + glob_)))[0]
rec = json.loads(Path(rec_path).read_text()); ref = rec['reference']
S = read_marcs_structure(require_codex_path('grids.marcs_standard') / rec['atmosphere_file']['name'])
sol = read_marcs_structure(require_codex_path('grids.marcs_standard') / 'p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod')
AO = float(ref['abundance'])
run = out / (f"{Path(rec_path).stem}_{mode}_os{iosmet}" + ("" if osscale is None else "_sc")); shutil.rmtree(run, ignore_errors=True); run.mkdir(parents=True)

lt = np.asarray(S['logtau_500'], float); T = np.asarray(S['temperature_K'], float)
ne_m3, _ = marcs_number_densities(S); ne = ne_m3 * 1e-6
ndep = len(lt)
with open(run / 'ATMOS', 'w') as f:
    f.write(f"  {S.get('model_id','MARCS')[:20]}\n  TAU(5000) SCALE\n* LG G\n  {S['logg_cgs']:.2f}\n* NDEP\n  {ndep}\n*\n")
    for a, b_, c in zip(lt, T, ne):
        f.write(f"  {a:14.6E}  {b_:14.6E}  {c:14.6E}  0.000000E+00  {S['vturb_m_s']/1e3:14.6E}\n")
with open(run / 'DSCALE', 'w') as f:
    f.write(f"  {S.get('model_id','MARCS')[:20]}\n  TAU(5000) SCALE\n  {ndep} {lt[0]-0.1:.1f}\n")
    for a in lt: f.write(f"  {a:13.6E}\n")
# ABUND: MULTI's 21 ABSDAT elements; solar = Asplund09 file (pymulti), model = the MARCS model's mix. O = node in both.
els = 'H HE C N O NE NA MG AL SI S K CA SC TI V CR MN FE CO NI'.split()
Z = dict(H=1, HE=2, C=6, N=7, O=8, NE=10, NA=11, MG=12, AL=13, SI=14, S=16, K=19, CA=20, SC=21, TI=22, V=23, CR=24, MN=25, FE=26, CO=27, NI=28)
a09 = {l.split()[0].upper(): float(l.split()[1]) for l in (PM / 'files/multi/mul23/run/abund').read_text().splitlines() if l.strip()}
with open(run / 'ABUND', 'w') as f:
    for e in els:
        v = a09.get(e, sol['abundances_dex'][Z[e]-1]) if mode == 'solar' else S['abundances_dex'][Z[e]-1]
        if e == 'O': v = AO
        f.write(f"{e:3s} {v:6.2f}\n")
atom = Path(require_codex_path('experiments.bride_oxygen_atom')).read_text().splitlines()
assert atom[1].split()[0] == '8.70', atom[1]
atom[1] = f"     {AO:5.2f}   15.999"
(run / 'ATOM').write_text('\n'.join(atom) + '\n')
(run / 'INPUT').write_text(
    "DIFF=5.0,ELIM1=1.0E-4,ELIM2=1.0E-5,QNORM=10.00,THIN=0.1,\n"
    "IATOM2=0,ICONV=0,IHSE=0,ILAMBD=1,IOPAC=1,ISTART=0,ISUM=0,\n"
    "ITMAX=200,ITRAN=0,NMU=5,\nIWABND=0,IWATMS=0,IWATOM=0,IWCHAN=0,IWDAMP=0,IWEMAX=1,IWEQW=0,\n"
    "IWEVEC=0,IWHEAD=0,IWHSE=0,IWLGMX=1,IWLINE=0,IWLTE=0,IWN=0,IWNIIT=0,\n"
    "IWOPAC=0,IWRAD=0,IWRATE=0,IWSTRT=0,IWTAUQ=0,IWTEST=0,IWWMAT=0,\n"
    f"IWJFIX=0,IWARN=2,IOPACL=0,ISCAT=1,INCRAD=0,INGACC=1,ICRSW=0,\nIOSMET={iosmet},EOSMET=0.5,\n"
    "IDL1=1,IDLNY=0,IDLCNT=0,IDLOPC=0\n")
if osscale is not None:
    feh = float(re.search(r'_([-+][0-9.]+)_a', Path(rec_path).stem).group(1))
    (run / 'OSSCALE').write_text(f"{10**feh:.6e}\n")
for n in ('absdat', 'absmet'):
    os.symlink(PM / f'files/multi/mul23/input/{n}', run / n.upper())
p = subprocess.run([str(M / 'source_dp/mul23lus.x')], cwd=run, capture_output=True, text=True, timeout=3600)
(run / 'stdout.txt').write_text(p.stdout + p.stderr)
if not (run / 'IDL1').exists():
    print('FAILED', rec_path, (p.stdout + p.stderr)[-1500:]); sys.exit(1)
r = mulrd(str(run / 'IDL1'))
b = np.asarray(r['n'], float) / np.asarray(r['nstar'], float).reshape(np.shape(r['n']))
b = b if b.shape[0] == 41 else b.T
bref = np.asarray(ref['departures'], float); bref = bref if bref.shape[0] == 41 else bref.T
ltr = np.log10(np.asarray(ref['tau'], float)); phot = (ltr >= -3) & (ltr <= 0)
sys.path.insert(0, 'scripts'); import bride_qa_benchmark as Q
TL = list(Q.TRIPLET_LEVELS)
assert len(ltr) == b.shape[1], (len(ltr), b.shape)
d = np.log10(b[TL]) - np.log10(bref[TL])
np.save(run / 'b.npy', b)
conv = re.findall(r'EMAX.*', (run / 'JOBLOG').read_text() if (run / 'JOBLOG').exists() else '')[-1:]
print(f"{Path(rec_path).stem[:22]} {mode:5s} os{iosmet}{'' if osscale is None else ' scaled'}  MULTI-vs-Gerber triplet phot max|d| {np.abs(d[:, phot]).max():.4f}  median|d| {np.median(np.abs(d[:, phot])):.4f}")
for tau in (-2, -1, 0):
    j = int(np.argmin(abs(ltr - tau)))
    print(f"   logtau {ltr[j]:+.2f} b_multi {np.round(b[TL, j], 4)} b_ref {np.round(bref[TL, j], 4)}")

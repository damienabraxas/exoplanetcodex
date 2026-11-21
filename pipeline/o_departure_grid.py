"""O I/II departure-coefficient grids on the Gerber+2023 MARCS nodes (RYA-1282, 2026-11-20).

Two engines on the same nodes, the same O41f atom and the same output format, so each is a
drop-in replacement for the Gerber+2023 deck (NLTEgrid4TS_O_MARCS + auxData):

  'bride'   the Bride (Lightweaver, MARCS/TS background with the UV hydrogen opacities).
            The frontier product; runs on the Mac only.
  'multi23' MULTI 2.3, the Gerber grid's own code, with its two metallicity artefacts
            removed (docs/audit/rya1282/lightweaver_differences.md): ABUND = the model's own
            mix (Gerber/pymulti: always solar -> n(H) and the H collisions off by up to x12),
            and the solar-mix ABSMET background-line table scaled by 10^[Fe/H] (Gerber:
            unscaled). The independent second engine; may run on Sirius (support job).

In both, A(O) of the node sets the atom abundance AND the O entry of the background mix.

Per MARCS model, one .npz with every A(O) of that model (resumable: existing files are
skipped). `assemble` writes the TS binary + aux in the Gerber layout that
pipeline.gerber_nlte.read_deck_node and TS's interpol_modeles_nlte read: a 1000-byte text
header, then per record a 500-byte model id, int32 ndep, int32 nk, float64 linear tau
(ndep), float64 departures stored (nk, ndep); aux pointers are 1-based byte offsets.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

HEADER_BYTES = 1000
ID_BYTES = 500
MULTI_ELEMENTS = 'H HE C N O NE NA MG AL SI S K CA SC TI V CR MN FE CO NI'.split()
_Z = dict(H=1, HE=2, C=6, N=7, O=8, NE=10, NA=11, MG=12, AL=13, SI=14, S=16, K=19, CA=20, SC=21,
          TI=22, V=23, CR=24, MN=25, FE=26, CO=27, NI=28)
GERBER_O_AUX = 'auxData_O_MARCS_May-21-2021.txt'


# ----------------------------------------------------------------------------- nodes
def gerber_nodes(aux_path) -> dict:
    """{model_id: dict(teff, logg, feh, alpha, vturb, abundances=[A(O)...])} from a Gerber aux."""
    out = {}
    for line in Path(aux_path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        f = line.split()
        mid = f[0].strip("'")
        d = out.setdefault(mid, dict(teff=float(f[1]), logg=float(f[2]), feh=float(f[3]),
                                     alpha=float(f[4]), vturb=float(f[6]), abundances=[]))
        d['abundances'].append(round(float(f[7]), 4))
    for d in out.values():
        d['abundances'] = sorted(set(d['abundances']))
    return out


def select_nodes(nodes: dict, *, teff=(0, 1e9), logg_min=-9., a_step=None) -> dict:
    """Subset by Teff range / minimum log g; a_step thins the A(O) axis to every a_step dex
    measured from the lowest offset A(O) - [Fe/H] (offsets are exact 0.1-dex steps)."""
    sel = {}
    for mid, d in nodes.items():
        if not (teff[0] <= d['teff'] <= teff[1]) or d['logg'] < logg_min:
            continue
        ab = d['abundances']
        if a_step:
            off0 = round(ab[0] - d['feh'], 4)
            ab = [a for a in ab if abs(((a - d['feh'] - off0) / a_step) - round((a - d['feh'] - off0) / a_step)) < 1e-6]
        sel[mid] = dict(d, abundances=ab)
    return sel


def priority(d: dict) -> tuple:
    """Run order: FGK dwarfs/subgiants first, nearest the Sun first."""
    return (d['logg'] < 3.0, abs(d['teff'] - 5750) / 250 + abs(d['logg'] - 4.5) + abs(d['feh']) / 0.5)


# ----------------------------------------------------------------------------- engines
def _load_model(model_id):
    from config.constants import require_codex_path
    from pipeline.bride_reference import read_marcs_structure
    return read_marcs_structure(require_codex_path('grids.marcs_standard') / f'{model_id}.mod')


def _atom():
    from config.constants import require_codex_path
    from pipeline.multi_atom import read_multi_atom
    return read_multi_atom(require_codex_path('experiments.bride_oxygen_atom'))


def solve_bride(S, A_O, atom=None):
    """Bride departures (nk, ndep) at one node; standard column settings (bride_qa_benchmark)."""
    from pipeline.bride_columns import column_context, marcs_column, solve_column
    atom = replace(atom or _atom(), abundance=float(A_O))
    ab = np.array(S['abundances_dex'], float); ab[7] = A_O
    ctx = column_context(marcs_column(S), atom, ng=True, abundances_dex=ab, logg_cgs=S['logg_cgs'],
                         vturb_m_s=S['vturb_m_s'], rays=3)
    n, hist, _, _ = solve_column(ctx, maxiter=6000, seconds=3000.)
    return n / ctx.activeAtoms[0].nStar, len(hist)


def _read_idl1(idl1):
    """pymulti's mulrd (the reference IDL1 reader), loaded from the held pymulti source."""
    from config.constants import require_codex_path
    src_path = require_codex_path('engines.mul23') / 'pymulti_multi.py'
    src = src_path.read_text()
    ns = {}
    exec('import numpy as np\nfrom scipy.io import FortranFile\n'
         + src[src.index('def mulrd'):src.index('def initiate_dscale')], ns)
    r = ns['mulrd'](str(idl1))
    n = np.asarray(r['n'], float)
    nstar = np.asarray(r['nstar'], float).reshape(n.shape)
    b = n / nstar
    return b if b.shape[0] != len(r['temp']) else b.T


def solve_multi23(S, A_O, *, workdir=None, scale_lines=True, solar_abund=False, iosmet=1):
    """MULTI 2.3 departures (nk, ndep) at one node, Gerber's code with the artefacts removed
    (solar_abund=True, scale_lines=False reproduces the Gerber records to 0.001-0.003 dex)."""
    from config.constants import require_codex_path
    from pipeline.bride_reference import marcs_number_densities
    home = require_codex_path('engines.mul23')
    exe = home / 'source_dp' / 'mul23lus.x'
    tmp = Path(tempfile.mkdtemp(prefix='mul23_', dir=workdir))
    try:
        lt = np.asarray(S['logtau_500'], float); T = np.asarray(S['temperature_K'], float)
        ne = marcs_number_densities(S)[0] * 1e-6
        mid = str(S.get('model_id', 'MARCS'))[:20]
        rows = ''.join(f"  {a:14.6E}  {t:14.6E}  {e:14.6E}  0.000000E+00  {S['vturb_m_s'] / 1e3:14.6E}\n"
                       for a, t, e in zip(lt, T, ne))
        (tmp / 'ATMOS').write_text(f"  {mid}\n  TAU(5000) SCALE\n* LG G\n  {S['logg_cgs']:.2f}\n"
                                   f"* NDEP\n  {len(lt)}\n*\n" + rows)
        (tmp / 'DSCALE').write_text(f"  {mid}\n  TAU(5000) SCALE\n  {len(lt)} {lt[0] - 0.1:.1f}\n"
                                    + ''.join(f"  {a:13.6E}\n" for a in lt))
        a09 = {l.split()[0].upper(): float(l.split()[1])
               for l in (home / 'input' / 'abund_asplund09').read_text().splitlines() if l.strip()}
        ab = S['abundances_dex']
        lines = []
        for e in MULTI_ELEMENTS:
            v = A_O if e == 'O' else (a09[e] if solar_abund and e in a09 else ab[_Z[e] - 1])
            lines.append(f"{e:3s} {v:6.2f}\n")
        (tmp / 'ABUND').write_text(''.join(lines))
        atom = Path(require_codex_path('experiments.bride_oxygen_atom')).read_text().splitlines()
        if atom[1].split()[0] != '8.70':
            raise ValueError(f'unexpected O41f abundance line: {atom[1]!r}')
        atom[1] = f"     {A_O:5.2f}   15.999"
        (tmp / 'ATOM').write_text('\n'.join(atom) + '\n')
        (tmp / 'INPUT').write_text(
            "DIFF=5.0,ELIM1=1.0E-4,ELIM2=1.0E-5,QNORM=10.00,THIN=0.1,\n"
            "IATOM2=0,ICONV=0,IHSE=0,ILAMBD=1,IOPAC=1,ISTART=0,ISUM=0,\n"
            "ITMAX=200,ITRAN=0,NMU=5,\nIWABND=0,IWATMS=0,IWATOM=0,IWCHAN=0,IWDAMP=0,IWEMAX=1,IWEQW=0,\n"
            "IWEVEC=0,IWHEAD=0,IWHSE=0,IWLGMX=1,IWLINE=0,IWLTE=0,IWN=0,IWNIIT=0,\n"
            "IWOPAC=0,IWRAD=0,IWRATE=0,IWSTRT=0,IWTAUQ=0,IWTEST=0,IWWMAT=0,\n"
            f"IWJFIX=0,IWARN=2,IOPACL=0,ISCAT=1,INCRAD=0,INGACC=1,ICRSW=0,\nIOSMET={iosmet},EOSMET=0.5,\n"
            "IDL1=1,IDLNY=0,IDLCNT=0,IDLOPC=0\n")
        if scale_lines:
            (tmp / 'OSSCALE').write_text(f"{10 ** float(S['feh']):.6e}\n")   # run_model sets S['feh']
        for n in ('absdat', 'absmet'):
            os.symlink(home / 'input' / n, tmp / n.upper())
        p = subprocess.run([str(exe)], cwd=tmp, capture_output=True, text=True, timeout=1800)
        if not (tmp / 'IDL1').exists():
            raise RuntimeError(f'MULTI failed: {(p.stdout + p.stderr)[-400:]}')
        out = (tmp / 'OUT').read_text(errors='replace')
        emax = [float(x) for x in re.findall(r'LOG\(EMAX \(.*?\)\)=\s*([-0-9.Ee+]+)', out)]
        if not emax or emax[-1] > -4.0:
            raise RuntimeError(f'MULTI not converged: last log EMAX {emax[-1] if emax else None}')
        return _read_idl1(tmp / 'IDL1'), len(emax)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


ENGINES = {'bride': solve_bride, 'multi23': solve_multi23}


def run_model(model_id: str, node: dict, engine: str, out_dir, workdir=None) -> Path:
    """Solve every A(O) of one MARCS model; write <out_dir>/<model_id>.npz (skip if present)."""
    out = Path(out_dir) / f'{model_id}.npz'
    if out.exists():
        return out
    S = _load_model(model_id)
    S['feh'] = node['feh']
    lt = np.asarray(S['logtau_500'], float)
    atom = _atom() if engine == 'bride' else None
    B, ok, steps, secs, err = [], [], [], [], []
    for A in node['abundances']:
        t0 = time.perf_counter()
        try:
            b, k = (solve_bride(S, A, atom) if engine == 'bride'
                    else solve_multi23(S, A, workdir=workdir))
            b = np.asarray(b, float)
            if b.shape != (41, len(lt)) or not np.all(np.isfinite(b)) or np.any(b <= 0):
                raise RuntimeError(f'bad departures shape {b.shape} / non-finite / non-positive')
            B.append(b); ok.append(True); steps.append(k); err.append('')
        except Exception as e:                                   # noqa: BLE001 - recorded per node
            B.append(np.full((41, len(lt)), np.nan)); ok.append(False); steps.append(-1); err.append(str(e)[:300])
        secs.append(time.perf_counter() - t0)
    tmp = out.with_suffix('.tmp.npz')
    np.savez_compressed(tmp, model_id=model_id, engine=engine, teff=node['teff'], logg=node['logg'],
                        feh=node['feh'], alpha=node['alpha'], vturb=node['vturb'],
                        abundances=np.array(node['abundances']), logtau=lt, b=np.array(B),
                        ok=np.array(ok), steps=np.array(steps), seconds=np.array(secs), error=np.array(err))
    os.replace(tmp, out)
    return out


# ----------------------------------------------------------------------------- TS writer
def assemble(npz_dir, bin_path, aux_path, *, header_text: str) -> dict:
    """Write every converged node of npz_dir as a TS NLTE grid (Gerber layout). Returns counts."""
    files = sorted(Path(npz_dir).glob('*.npz'))
    n_rec = n_bad = 0
    head = header_text.encode('ascii')
    if len(head) > HEADER_BYTES:
        raise ValueError('header too long')
    with open(bin_path, 'wb') as fb, open(aux_path, 'w') as fa:
        fa.write('#atmos ID, Teff [K], log(g) [cgs], [Fe/H], [alpha/Fe], mass, Vturb [km/s], A(X), pointer\n')
        fb.write(head.rjust(HEADER_BYTES))
        for f in files:
            z = np.load(f)
            mid = str(z['model_id'])
            tau = 10. ** np.asarray(z['logtau'], float)
            mass = float(re.search(r'_m([0-9.]+)_', mid).group(1)) if re.search(r'_m([0-9.]+)_', mid) else 0.
            for A, b, ok in zip(z['abundances'], z['b'], z['ok']):
                if not ok:
                    n_bad += 1
                    continue
                ptr = fb.tell() + 1
                nk, ndep = b.shape
                fb.write(mid.encode('ascii').ljust(ID_BYTES))
                fb.write(np.array([ndep, nk], '<i4').tobytes())
                fb.write(tau.astype('<f8').tobytes())
                fb.write(np.ascontiguousarray(b, dtype='<f8').tobytes())      # (nk, ndep) as Gerber
                fa.write(f"'{mid}'\t{float(z['teff']):.4f}\t{float(z['logg']):.4f}\t{float(z['feh']):.4f}\t"
                         f"{float(z['alpha']):.4f}\t{mass:.2f}\t{float(z['vturb']):.2f}\t{float(A):.4f}"
                         f"{ptr:>22d} \n")
                n_rec += 1
    return dict(records=n_rec, failed_nodes=n_bad, models=len(files))


def read_record(bin_path, pointer: int):
    """(model_id, log tau, departures (nk, ndep)) at a 1-based aux pointer - the layout check."""
    with open(bin_path, 'rb') as fh:
        fh.seek(pointer - 1)
        mid = fh.read(ID_BYTES).decode('ascii').strip()
        ndep, nk = np.frombuffer(fh.read(8), '<i4')
        tau = np.frombuffer(fh.read(8 * ndep), '<f8')
        dep = np.frombuffer(fh.read(8 * ndep * nk), '<f8').reshape(nk, ndep)
    return mid, np.log10(tau), dep

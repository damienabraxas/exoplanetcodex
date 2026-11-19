"""RYA-1282 Bride G5 engine check: OpenMULTI3D 1D NLTE on the same STAGGER columns.

Each column (as the Bride reads it, `bride_columns.stagger_column`) is written as an
OpenMULTI3D MULTI-1D model (`atmos.<name>`: z [cm, outward], T [K], ne [cm^-3],
vz [km/s, outward], vmic [km/s] = 0) and solved with the same held O atom (o41f_ref,
sha d058b75b) at the same A(O). Disk-centre (mu = 1) line and continuum intensities are
measured on the Bride's EW windows. Engine differences that remain by construction:
OM derives n_H from its own EOS given (T, ne); the Bride uses the snapshot's n_H;
background opacities and collision-table extrapolation above 10000 K are each engine's own.

Heights MUST be in cm: OM reads a 1D model whose z spans < 1e4 as log tau.

--format multi1d (above) lets OM rebuild the density from (T, ne) with its own EOS: that
is NOT the Bride's input. --format text writes `<name>.txt` rows (z cm, T, Pe = ne k T cgs,
rho g/cm^3, vmic, vz km/s) and runs with use_rho=T use_ne=T, so OM gets the snapshot's
density and electron density exactly as the Bride does. rho = n_Htot * MASS_PER_H_AMU * amu,
MASS_PER_H_AMU = 1.365315 from the Bride's MARCS 5750/4.5 abundances (lightweaver PeriodicTable).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ATOM_SHA = 'd058b75bdec80c5317523a6d7f185a64ca59144d3d494a17eeabee5f10d04b57'
ENGINE_COMMIT = '32085d26c46e9919a986397b5514a9728abe2074'
MASS_PER_H_AMU = 1.365315
AMU_G = 1.66053906660e-24
K_B_CGS = 1.380649e-16


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_multi1d(col, path, name):
    rows = ['* RYA-1282: STAGGER column written for OpenMULTI3D (multi1d); z cm outward, vz km/s outward',
            name, 'Height', '4.44', str(col['ndep'])]
    for z, t, ne, vz in zip(col['height_m'] * 100., col['temperature_K'], col['ne_m3'] * 1e-6,
                            col['vlos_m_s'] * 1e-3):
        rows.append(f'{z:.6e} {t:.3f} {ne:.6e} {vz:.6f} 0.0')
    Path(path).write_text('\n'.join(rows) + '\n')


def write_text1d(col, path, name):
    rows = ['* RYA-1282: STAGGER column for OpenMULTI3D (text): z cm, T K, Pe cgs, rho g/cm3, vmic km/s, vz km/s',
            name, str(col['ndep'])]
    ne = col['ne_m3'] * 1e-6
    rho = col['nHtot_m3'] * 1e-6 * MASS_PER_H_AMU * AMU_G
    for z, t, n, r, vz in zip(col['height_m'] * 100., col['temperature_K'], ne, rho, col['vlos_m_s'] * 1e-3):
        rows.append(f'{z:.6e} {t:.3f} {n * K_B_CGS * t:.6e} {r:.6e} 0.0 {vz:.6f}')
    Path(path).write_text('\n'.join(rows) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--columns', nargs='+', required=True, help='ix,iy')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--abundance', type=float, default=8.70)
    p.add_argument('--convlim', type=float, default=1e-4)
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--format', choices=('text', 'multi1d'), default='text')
    a = p.parse_args()
    work = Path(os.environ.get('CODEX_WORK_ROOT', Path.home() / 'codex'))
    app = work / 'bride_inputs/rya1281/openmulti3d-reference-v1/experiments/OpenMULTI3D'
    py = work / 'venv_openmulti3d/bin/python'
    if Path(sys.prefix).resolve() != py.parent.parent.resolve() and not os.environ.get('_OM3D_CHILD'):
        env = dict(os.environ, _OM3D_CHILD='1',
                   PYTHONPATH=f"{app}:{app / 'HOP'}:{ROOT}:{Path.home() / 'bride_rya1118'}")
        sys.exit(subprocess.call([str(py), __file__, *sys.argv[1:]], env=env))
    sys.path[:0] = [str(ROOT), str(Path.home() / 'bride_rya1118')]
    import tumult as tm
    from multi3d_read import read_atmos_multi3d
    from pipeline.bride_columns import EW_WINDOWS_AIR_A, stagger_column
    from pipeline.wavelength_util import air_to_vac
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=app, text=True).strip()
    if head != ENGINE_COMMIT:
        raise ValueError(f'engine at {head}, pinned {ENGINE_COMMIT}')
    atom = work / 'bride_inputs/rya1281/atom.o41f'
    if sha(atom) != ATOM_SHA:
        raise ValueError('held O atom changed')
    stem = 't5777g44m0005_00091_80x80x240'
    base = Path.home() / 'bride_rya1118'
    cube = read_atmos_multi3d(str(base / f'mesh.{stem}'), str(base / f'atm3d.{stem}'))
    a.output.mkdir(parents=True, exist_ok=False)
    wave = np.arange(7770.5, 7777.0 + 0.0025, 0.005)
    lam = a.output / 'lam_vac_A.txt'
    np.savetxt(lam, air_to_vac(wave), fmt='%.6f')
    summary = dict(schema='codex.openmulti3d_columns/1', ticket='RYA-1282', engine_commit=head,
                   atom_sha256=ATOM_SHA, abundance_O=a.abundance, atmos_format=a.format,
                   mass_per_h_amu=MASS_PER_H_AMU, convlim=a.convlim, snapshot=stem,
                   ew_windows_air_A=EW_WINDOWS_AIR_A, columns=[])
    for c in a.columns:
        ix, iy = (int(v) for v in c.split(','))
        name = f'x{ix:02d}_y{iy:02d}'
        d = a.output / name
        d.mkdir()
        (d / 'input').symlink_to(app / 'input', target_is_directory=True)
        col = stagger_column(cube, ix, iy)
        if a.format == 'text':
            atm = f'{name}.txt'
            write_text1d(col, d / atm, name)
            extra = " atmos_format='text' use_rho=T use_ne=T"
        else:
            atm = f'atmos.{name}'
            write_multi1d(col, d / atm, name)
            extra = ''
        (d / 'run.nml').write_text('\n'.join([
            "&io_params datadir='data' gb_step=6.0 do_trace=F /",
            "&timer_params sec_per_report=1d8 verbose=1 /",
            f"&atmos_params atmos_file='{atm}'{extra} /",
            f"&atom_params atom_file='{atom}' abundance={a.abundance:.4f} output_from_line_list=F /",
            f"&nlte_params convlim={a.convlim:.1e} /",
            f"&spectrum_params lam_file='{lam.resolve()}' in_air=F /",
            "&run_params mem_limit_gb=4. /"]) + '\n')
        t0 = time.perf_counter()
        with open(d / 'native.log', 'w') as log:
            rc = subprocess.run([str(app / 'dispatch.x'), 'run.nml'], cwd=d, stdout=log, stderr=subprocess.STDOUT,
                                env=dict(os.environ, OMP_NUM_THREADS=str(a.threads), OMP_STACKSIZE='256M')).returncode
        rec = dict(column=[ix, iy], returncode=rc, seconds=time.perf_counter() - t0)
        if rc == 0:
            run = tm.read(str(d / 'data/run'))
            if abs(float(run.mus[0]) - 1.) > 1e-9:
                raise ValueError('first OM ray is not mu = 1')
            w = np.asarray(run.lam, float)          # tumult: air A
            line = np.asarray(run.ie, float)[:, np.asarray(run.imus) == 0].mean(axis=1)
            cont = np.asarray(run.ie_cnt, float)[:, np.asarray(run.imus) == 0].mean(axis=1)
            r = 1 - line / cont
            rec['ew_mA'] = [float(np.trapezoid(r[(w >= lo) & (w < hi)], w[(w >= lo) & (w < hi)]) * 1e3)
                            for lo, hi in EW_WINDOWS_AIR_A]
            np.savez_compressed(a.output / f'{name}.npz', wavelength_air_A=w, line=line, continuum=cont)
        summary['columns'].append(rec)
        (a.output / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
        print(name, rc, round(rec['seconds'], 1), np.round(rec.get('ew_mA', []), 3), flush=True)


if __name__ == '__main__':
    main()

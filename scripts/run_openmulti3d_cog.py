"""OpenMULTI3D 1D LTE+NLTE curve of growth for any element: atom x abundance grid -> per-line EWs.

One native run per (atom, abundance) gives NLTE, LTE and continuum flux (tumult reader).
Per-line EW = integral of (1 - F/Fc) over the line's declared integration window, so LTE and
NLTE share one definition and the NLTE correction at fixed EW is differential.

    python scripts/run_openmulti3d_cog.py --spec config/experiments/<spec>.json --output <dir>

Requires the venv with tumult's dependencies (CODEX_OPENMULTI3D_PY, default
$CODEX_WORK_ROOT/venv_openmulti3d/bin/python); this driver re-executes itself under it.
"""
import argparse, csv, hashlib, json, os, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def namelist(spec, atom, abundance, lam_file, atmos):
    ll = spec.get('linelist')
    lines = [
        "&io_params datadir='data' gb_step=6.0 do_trace=F /",
        "&timer_params sec_per_report=1d8 verbose=1 /",
        f"&atmos_params atmos_file='{atmos}' /",
        f"&atom_params atom_file='{atom}' abundance={abundance:.4f} output_from_line_list=F /",
        f"&nlte_params convlim={spec['convlim']:.1e} /",
        f"&spectrum_params lam_file='{lam_file}' in_air=F /",
        "&run_params mem_limit_gb=4. /",
    ]
    if ll:
        lines.insert(5, f"&linelist_params linelist_file='{ll}' match_nlte_levels=T /")
    return '\n'.join(lines) + '\n'


def output_axis_vac(windows_air, step):
    import numpy as np
    sys.path.insert(0, str(ROOT))
    from pipeline.wavelength_util import air_to_vac
    air = np.concatenate([np.arange(lo, hi + step / 2, step) for lo, hi in windows_air])
    return air_to_vac(air)


def measure(run_dir, spec):
    import numpy as np
    import tumult as tm
    run = tm.read(str(run_dir))
    lam = np.asarray(run.lam, float)  # tumult: standard air, Angstrom
    out = {}
    for kind, f in [('nlte', run.flux), ('lte', run.flux_lte)]:
        r = np.asarray(f, float) / np.asarray(run.flux_cnt, float)
        for line in spec['lines']:
            lo, hi = line['integrate_air_A']
            m = (lam >= lo) & (lam <= hi)
            if m.sum() < 20:
                raise ValueError(f"too few samples in {line['name']} window")
            out[(line['name'], kind)] = 1e3 * float(np.trapezoid(1 - r[m], lam[m]))
    return out, list(map(float, run.mus))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spec', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    work = Path(os.environ.get('CODEX_WORK_ROOT', Path.home() / 'codex'))
    py = os.environ.get('CODEX_OPENMULTI3D_PY', str(work / 'venv_openmulti3d/bin/python'))
    spec = json.loads(a.spec.read_text())
    app = Path(spec['engine_app'].replace('$CODEX_WORK_ROOT', str(work)))
    if Path(sys.prefix).resolve() != Path(py).parent.parent.resolve() and not os.environ.get('_OM3D_CHILD'):
        env = dict(os.environ, _OM3D_CHILD='1', PYTHONPATH=f"{app}:{app/'HOP'}:{ROOT}")
        sys.exit(subprocess.call([py, __file__, *sys.argv[1:]], env=env))
    sys.path.insert(0, str(ROOT))
    import numpy as np
    exe = app / 'dispatch.x'
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=app, text=True).strip()
    if head != spec['engine_commit']:
        raise ValueError(f'engine at {head}, spec pins {spec["engine_commit"]}')
    atmos = app / spec['atmosphere_rel']
    if sha(atmos) != spec['atmosphere_sha256']:
        raise ValueError('atmosphere hash changed')
    a.output.mkdir(parents=True, exist_ok=False)
    np.savetxt(a.output / 'lam_vac_A.txt', output_axis_vac(spec['windows_air_A'], spec['step_air_A']), fmt='%.6f')
    rows, receipts = [], {}
    for atom_name, atom_path in spec['atoms'].items():
        atom_path = Path(atom_path.replace('$CODEX_WORK_ROOT', str(work)))
        h = sha(atom_path)
        if h != spec['atom_sha256'][atom_name]:
            raise ValueError(f'{atom_name} hash changed')
        receipts[atom_name] = h
        for A in spec['abundances']:
            tag = f'{atom_name}_A{A:.3f}'
            d = a.output / 'runs' / tag
            d.mkdir(parents=True)
            (d / 'input').symlink_to(app / 'input', target_is_directory=True)
            (d / 'run.nml').write_text(namelist(spec, atom_path, A, (a.output / 'lam_vac_A.txt').resolve(), atmos))
            t0 = time.perf_counter()
            with (d / 'native.log').open('w') as log:
                rc = subprocess.run([str(exe), 'run.nml'], cwd=d, stdout=log, stderr=subprocess.STDOUT,
                                    env=dict(os.environ, OMP_STACKSIZE='256M', OMP_NUM_THREADS=str(spec.get('threads', 4))),
                                    timeout=600).returncode
            logtxt = (d / 'native.log').read_text()
            if rc or 'terminated successfully' not in logtxt:
                raise RuntimeError(f'{tag} failed rc={rc}')
            iters = [l for l in logtxt.splitlines() if '/ 99' in l or '/ 9' in l]
            ews, mus = measure(d / 'data' / 'run', spec)
            for line in spec['lines']:
                rows.append(dict(atom=atom_name, A=A, line=line['name'], wavelength_air_A=line['wavelength_air_A'],
                                 ew_lte_mA=round(ews[(line['name'], 'lte')], 4), ew_nlte_mA=round(ews[(line['name'], 'nlte')], 4),
                                 iterations=len(iters) - 1, seconds=round(time.perf_counter() - t0, 2)))
            print(tag, {k[0] + ':' + k[1]: round(v, 2) for k, v in ews.items()}, flush=True)
    with (a.output / 'cog.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    (a.output / 'report.json').write_text(json.dumps(dict(
        schema='codex.openmulti3d_cog/1', spec=str(a.spec), spec_sha256=sha(a.spec), generator_sha256=sha(__file__),
        engine_commit=head, executable_sha256=sha(exe), atmosphere=str(atmos), atmosphere_sha256=sha(atmos),
        atom_sha256=receipts, flux_mus=mus, ew_definition='1e3*trapz(1-F/Fc) over integrate_air_A, flux, tumult air axis',
        triplet_gf_source='atom transitions (output_from_line_list=F); background lines from linelist if declared'), indent=2) + '\n')


if __name__ == '__main__':
    main()

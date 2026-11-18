"""Build an OpenMULTI3D 1D NLTE departure deck on a MARCS model, for the ORCH Engine-B route.

    <openmulti3d venv python> scripts/build_openmulti3d_deck.py --key O --out <deck dir>

For each abundance on the key's axis: one native NLTE solve (no output spectrum) on the
MARCS model, then b(level, depth) mapped onto the MARCS file's own 56 layers by COLUMN MASS.
Mapping by the engine's recomputed log tau500 would misplace the departures: its tau scale
differs from MARCS's lgTau5 (T(tau500) off by up to 74 K near tau=1 at the solar node),
while column mass reproduces the MARCS temperatures to a few K. The temperature residual of
the mapping is recorded and gated. Layers outside the engine's mesh hold the edge value.

Run by `pipeline.openmulti3d_nlte.ensure_deck`; needs the tumult venv (engines.openmulti3d_python)
and PYTHONPATH=<engine>:<engine>/HOP:<repo>.
"""
import argparse, hashlib, json, os, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
MAX_MAPPING_DT_K = 25.0


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def marcs_layers(path):
    txt = Path(path).read_text().splitlines()

    def block(key):
        i = next(k for k, l in enumerate(txt) if key in l)
        rows = []
        for l in txt[i + 1:]:
            p = l.split()
            if not p or not p[0].lstrip('-').isdigit():
                break
            rows.append([float(x) for x in p])
        return np.array(rows)
    b1, b2 = block('lgTau5'), block('RHOX')
    if len(b1) != len(b2):
        raise ValueError('MARCS blocks differ in depth count')
    return dict(lgTauR=b1[:, 1], lgTau5=b1[:, 2], T=b1[:, 4], rhox=b2[:, 7])


def departures_on_marcs(run_dir, marcs):
    import tumult as tm
    run = tm.read(str(run_dir))
    _, dep = run.get_depart()
    dep = np.asarray(dep, float)                      # (ndepth, nlevel), engine order
    at = run.atmos
    z = np.asarray(at.zz, float).squeeze() * 1e8      # Mm -> cm
    rho = np.asarray(at.rho, float).squeeze()
    T = np.asarray(at.temp, float).squeeze()
    top_first = rho[0] < rho[-1]
    if not top_first:
        z, rho, T, dep = z[::-1], rho[::-1], T[::-1], dep[::-1]
    m = marcs['rhox'][0] + np.concatenate([[0.0], np.cumsum(0.5 * (rho[1:] + rho[:-1]) * np.abs(np.diff(z)))])
    lm, lmm = np.log10(m), np.log10(marcs['rhox'])
    inside = (lmm >= lm[0]) & (lmm <= lm[-1])
    dT = np.interp(lmm, lm, T) - marcs['T']
    b = np.column_stack([np.interp(lmm, lm, dep[:, k]) for k in range(dep.shape[1])])   # edge-held outside
    return b, dict(max_abs_dT_inside_K=float(np.max(np.abs(dT[inside]))), rms_dT_inside_K=float(np.sqrt(np.mean(dT[inside] ** 2))),
                   n_layers_outside_engine_mesh=int((~inside).sum()), engine_ndepth=int(len(m)),
                   engine_column_mass_range=[float(m[0]), float(m[-1])])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--key', required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    from pipeline import openmulti3d_nlte as om
    spec = om.deck_spec(a.key)
    app, exe, engine = om.engine()
    atom, atmos = Path(spec['atom_path']), Path(spec['atmosphere'])
    if sha(atom) != spec['atom_sha256']:
        raise ValueError(f'{atom} hash differs from the registered {spec["atom_sha256"]}')
    marcs = marcs_layers(atmos)
    work = a.out.with_name(a.out.name + '.runs')
    a.out.mkdir(parents=True, exist_ok=False)
    work.mkdir(parents=True, exist_ok=True)
    bs, audits = [], []
    for A in spec['abundance_axis']:
        d = work / f'A{A:.3f}'
        if d.exists():
            raise FileExistsError(d)
        d.mkdir()
        (d / 'input').symlink_to(app / 'input', target_is_directory=True)
        (d / 'run.nml').write_text(
            "&io_params datadir='data' gb_step=6.0 do_trace=F /\n"
            f"&atmos_params atmos_file='{atmos}' /\n"
            f"&atom_params atom_file='{atom}' abundance={A:.4f} /\n"
            f"&nlte_params convlim={spec['convlim']:.1e} /\n"
            "&rt_output_params solver='none' /\n&run_params mem_limit_gb=8. /\n")
        t0 = time.perf_counter()
        with (d / 'native.log').open('w') as log:
            rc = subprocess.run([str(exe), 'run.nml'], cwd=d, stdout=log, stderr=subprocess.STDOUT, timeout=spec['seconds'],
                                env=dict(os.environ, OMP_STACKSIZE='256M', OMP_NUM_THREADS=str(spec.get('threads', 4)))).returncode
        logtxt = (d / 'native.log').read_text()
        if rc or 'terminated successfully' not in logtxt:
            raise RuntimeError(f'{a.key} A={A}: native solve failed (rc={rc}); see {d}/native.log')
        b, audit = departures_on_marcs(d / 'data' / 'run', marcs)
        if audit['max_abs_dT_inside_K'] > MAX_MAPPING_DT_K:
            raise ValueError(f'column-mass mapping misplaces depths: max |dT| {audit["max_abs_dT_inside_K"]:.1f} K')
        if b.shape != (len(marcs['T']), spec['nk']) or not np.all(np.isfinite(b)) or np.any(b <= 0):
            raise ValueError(f'invalid departures {b.shape}')
        audit.update(A=A, seconds=round(time.perf_counter() - t0, 2),
                     iterations=sum(1 for l in logtxt.splitlines() if '/ 99' in l) - 1)
        bs.append(b); audits.append(audit)
        print(json.dumps(audit), flush=True)
    np.savez_compressed(a.out / 'deck.npz', abundance_axis=np.array(spec['abundance_axis']), departures=np.array(bs),
                        lgTau5=marcs['lgTau5'], lgTauR=marcs['lgTauR'])
    (a.out / 'deck.json').write_text(json.dumps(dict(
        schema='codex.openmulti3d_deck/1', key=a.key, spec=spec, engine=engine, atmosphere_sha256=sha(atmos),
        deck_npz_sha256=sha(a.out / 'deck.npz'), generator_sha256=sha(__file__), audits=audits,
        mapping='b(level, MARCS layer) by log column mass; engine log tau500 NOT used (its scale differs from MARCS lgTau5)',
        built_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())), indent=2) + '\n')


if __name__ == '__main__':
    main()

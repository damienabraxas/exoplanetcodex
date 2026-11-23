#!/usr/bin/env python
"""RYA-1282: build an O departure grid on the Gerber+2023 MARCS nodes (pipeline.o_departure_grid).

  run      --engine bride|multi23 --out DIR [--teff 4000:7000] [--logg-min 3] [--a-step 0.2]
           [--feh -1.5:1.0] [--workers N] [--shard i/n] [--limit K]
  run      --engine bride --out RETRY_DIR --retry-from BASE_DIR [--aux FILE] [--workers N]
  assemble --npz DIR --bin FILE --aux-out FILE --header TEXT [--aux FILE]
  audit    --out DIR [--aux FILE]
  status   --out DIR

The Bride engine runs on the Mac only (Ryan 11-20); multi23 is a support job and may run on
Sirius. The node list is the Gerber O aux file (grids.gerber_ts).
"""
import argparse
import json
import os
import sys
import time
import subprocess
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np                                                          # noqa: E402
import pipeline.o_departure_grid as G                                        # noqa: E402


def _nodes(a):
    from config.constants import require_codex_path
    aux = Path(a.aux) if a.aux else require_codex_path('grids.gerber_ts') / G.GERBER_O_AUX
    nodes = G.gerber_nodes(aux)
    t0, t1 = (float(x) for x in a.teff.split(':'))
    sel = G.select_nodes(nodes, teff=(t0, t1), logg_min=a.logg_min, a_step=a.a_step)
    f0, f1 = (float(x) for x in a.feh.split(':'))
    # Ryan 2026-11-21: grid bounds [Fe/H] -1.5..+1.0. Gerber's -5..+1 node list was built for Galactic
    # surveys; no known planet host is below -1.0 (NASA Exoplanet Archive, 4197 hosts) and the
    # Codex targets span -0.49 (tau Cet) .. +0.31 (55 Cnc A).
    sel = {m: v for m, v in sel.items() if f0 - 1e-6 <= v['feh'] <= f1 + 1e-6}
    order = sorted(sel, key=lambda m: G.priority(sel[m]))
    if a.shard:
        i, n = (int(x) for x in a.shard.split('/'))
        order = order[i::n]
    return sel, order[:a.limit] if a.limit else order


def _audit(out, sel, order, engine):
    """Check an exact selection and retain the abundance and error for each failed node."""
    missing, failed, files = [], {}, []
    for model_id in order:
        path = out / f'{model_id}.npz'
        if not path.exists():
            missing.append(model_id)
            continue
        expected = sel[model_id]['abundances']
        with np.load(path) as z:
            if (str(z['model_id']) != model_id or str(z['engine']) != engine
                    or not np.array_equal(z['abundances'], expected)
                    or z['ok'].shape != (len(expected),)
                    or z['error'].shape != (len(expected),)):
                raise ValueError(f'{path}: result does not match the selected grid')
            failures = [
                dict(abundance=abundance, error=str(error))
                for abundance, ok, error in zip(expected, z['ok'], z['error']) if not ok
            ]
        if failures:
            failed[model_id] = failures
        files.append(path)
    return files, missing, failed


def _manifest_selection(out, aux):
    manifest = json.loads((out / 'manifest.json').read_text())
    if 'retry_from' in manifest:
        engine, sel, _ = _manifest_selection(Path(manifest['retry_from']), aux)
        retry = {m: dict(sel[m], abundances=ab) for m, ab in manifest['selection'].items()}
        return engine, retry, list(retry)
    selection = argparse.Namespace(
        aux=aux, teff=manifest['teff'], feh=manifest['feh'], logg_min=manifest['logg_min'],
        a_step=manifest['a_step'], shard=manifest.get('shard'), limit=manifest.get('limit'),
    )
    sel, order = _nodes(selection)
    if (len(order) != manifest['models']
            or sum(len(sel[m]['abundances']) for m in order) != manifest['nodes']):
        raise ValueError('aux selection does not match the manifest totals')
    return manifest['engine'], sel, order


def _audit_manifest(out, aux):
    engine, sel, order = _manifest_selection(out, aux)
    return _audit(out, sel, order, engine)


def _retry_selection(source, out, aux, engine):
    if source.resolve() == out.resolve():
        raise ValueError('retry output must be separate from the original grid')
    source_engine, sel, order = _manifest_selection(source, aux)
    if engine != source_engine:
        raise ValueError('retry engine must match the original grid')
    _, missing, failed = _audit(source, sel, order, engine)
    if missing:
        raise ValueError(f'original grid still has {len(missing)} missing models')
    retry = {
        model_id: dict(sel[model_id], abundances=[node['abundance'] for node in failures])
        for model_id, failures in failed.items()
    }
    return retry, list(retry)



def cmd_audit(a):
    files, missing, failed = _audit_manifest(Path(a.out), a.aux)
    print(json.dumps(dict(present_models=len(files), missing_models=missing, failed=failed), indent=2))


def _one(args):
    mid, node, engine, out, workdir = args
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    t0 = time.perf_counter()
    f = G.run_model(mid, node, engine, out, workdir)
    z = np.load(f)
    return mid, int(z['ok'].sum()), len(z['ok']), time.perf_counter() - t0


def _free_gb():
    """Free + inactive + speculative memory (GB) on the Mac; inf where vm_stat is absent."""
    try:
        out = subprocess.check_output(['vm_stat'], text=True)
    except (OSError, subprocess.CalledProcessError):
        return float('inf')
    page = int(out.split('page size of ')[1].split()[0])
    pages = {l.split(':')[0]: int(l.split(':')[1].strip().rstrip('.')) for l in out.splitlines()[1:] if ':' in l}
    return sum(pages.get(k, 0) for k in ('Pages free', 'Pages inactive', 'Pages speculative')) * page / 1e9


def cmd_run(a):
    if a.engine == 'bride' and os.uname().sysname != 'Darwin':
        sys.exit('the Bride grid runs on the Mac only (Ryan 2026-11-20)')
    out = Path(a.out)
    if a.retry_from:
        sel, order = _retry_selection(Path(a.retry_from), out, a.aux, a.engine)
        manifest = dict(engine=a.engine, retry_from=str(Path(a.retry_from).resolve()),
                        selection={m: sel[m]['abundances'] for m in order})
    else:
        sel, order = _nodes(a)
        manifest = dict(engine=a.engine, teff=a.teff, feh=a.feh, logg_min=a.logg_min,
                        a_step=a.a_step, shard=a.shard, limit=a.limit)
    manifest.update(models=len(order), nodes=sum(len(sel[m]['abundances']) for m in order))
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = out / 'manifest.json'
    if a.retry_from and manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        raise ValueError('retry directory already holds a different selection')
    if a.retry_from:
        _audit(out, sel, order, a.engine)
    todo = [m for m in order if not (out / f'{m}.npz').exists()]
    manifest_path.write_text(json.dumps(manifest, indent=1))
    print(f'{a.engine}: {len(order)} models, {len(todo)} to do, {a.workers} workers', flush=True)
    log = open(out / 'progress.log', 'a')
    # A Lightweaver context leaks ~17 MB per Bride solve and a pooled worker never returns it
    # (11-21: six workers reached 13-19 GB each and killed the 16 GB Mac). Each model therefore
    # gets a fresh process (max_tasks_per_child=1), at most `workers` models are in flight, and
    # no model starts while free memory is below --min-free-gb.
    pending, k, queue = set(), 0, list(todo)
    with ProcessPoolExecutor(a.workers, max_tasks_per_child=1) as ex:
        while queue or pending:
            while queue and len(pending) < a.workers:
                if (pending or a.retry_from) and _free_gb() < a.min_free_gb:
                    break
                m = queue.pop(0)
                pending.add(ex.submit(_one, (m, sel[m], a.engine, str(out), a.workdir)))
            if not pending:
                time.sleep(30)
                continue
            done, pending = wait(pending, timeout=60, return_when=FIRST_COMPLETED)
            for f in done:
                k += 1
                try:
                    mid, nok, n, dt = f.result()
                    log.write(f'{time.strftime("%H:%M:%S")} {k}/{len(todo)} {mid} ok {nok}/{n} {dt:.0f}s\n')
                except Exception as e:                                 # noqa: BLE001
                    log.write(f'{time.strftime("%H:%M:%S")} {k}/{len(todo)} ERROR {e}\n')
                log.flush()


def cmd_status(a):
    files = sorted(Path(a.out).glob('*.npz'))
    ok = n = 0; secs = []
    for f in files:
        z = np.load(f); ok += int(z['ok'].sum()); n += len(z['ok']); secs += list(z['seconds'])
    m = json.loads((Path(a.out) / 'manifest.json').read_text()) if (Path(a.out) / 'manifest.json').exists() else {}
    print(f'models {len(files)}/{m.get("models")}  nodes ok {ok}/{n} (target {m.get("nodes")})  '
          f'median {np.median(secs) if secs else float("nan"):.1f} s/node')


def cmd_assemble(a):
    files, missing, failed = _audit_manifest(Path(a.npz), a.aux)
    if missing or failed:
        sys.exit(f'Grid incomplete: {len(missing)} missing models, '
                 f'{sum(map(len, failed.values()))} failed abundance nodes. Run audit for details.')
    print(G.assemble(a.npz, a.bin, a.aux_out, header_text=a.header, files=files))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    s = p.add_subparsers(dest='cmd', required=True)
    r = s.add_parser('run'); r.add_argument('--engine', choices=sorted(G.ENGINES), required=True)
    r.add_argument('--out', required=True); r.add_argument('--aux'); r.add_argument('--teff', default='0:99999')
    r.add_argument('--feh', default='-99:99', help='[Fe/H] range lo:hi')
    r.add_argument('--logg-min', type=float, default=-9.); r.add_argument('--a-step', type=float)
    r.add_argument('--workers', type=int, default=4); r.add_argument('--shard'); r.add_argument('--limit', type=int)
    r.add_argument('--retry-from'); r.add_argument('--workdir'); r.add_argument('--min-free-gb', type=float, default=4.0)
    t = s.add_parser('status'); t.add_argument('--out', required=True)
    d = s.add_parser('audit'); d.add_argument('--out', required=True); d.add_argument('--aux')
    m = s.add_parser('assemble'); m.add_argument('--aux'); m.add_argument('--npz', required=True); m.add_argument('--bin', required=True)
    m.add_argument('--aux-out', required=True); m.add_argument('--header', required=True)
    a = p.parse_args()
    dict(run=cmd_run, status=cmd_status, audit=cmd_audit, assemble=cmd_assemble)[a.cmd](a)


if __name__ == '__main__':
    main()

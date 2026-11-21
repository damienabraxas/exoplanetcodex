#!/usr/bin/env python
"""RYA-1282: build an O departure grid on the Gerber+2023 MARCS nodes (pipeline.o_departure_grid).

  run      --engine bride|multi23 --out DIR [--teff 4000:7000] [--logg-min 3] [--a-step 0.2]
           [--workers N] [--shard i/n] [--limit K]
  assemble --npz DIR --bin FILE --aux FILE --header TEXT
  status   --out DIR

The Bride engine runs on the Mac only (Ryan 11-20); multi23 is a support job and may run on
Sirius. The node list is the Gerber O aux file (grids.gerber_ts).
"""
import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
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
    order = sorted(sel, key=lambda m: G.priority(sel[m]))
    if a.shard:
        i, n = (int(x) for x in a.shard.split('/'))
        order = order[i::n]
    return sel, order[:a.limit] if a.limit else order


def _one(args):
    mid, node, engine, out, workdir = args
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    t0 = time.perf_counter()
    f = G.run_model(mid, node, engine, out, workdir)
    z = np.load(f)
    return mid, int(z['ok'].sum()), len(z['ok']), time.perf_counter() - t0


def cmd_run(a):
    if a.engine == 'bride' and os.uname().sysname != 'Darwin':
        sys.exit('the Bride grid runs on the Mac only (Ryan 2026-11-20)')
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    sel, order = _nodes(a)
    todo = [m for m in order if not (out / f'{m}.npz').exists()]
    (out / 'manifest.json').write_text(json.dumps(dict(engine=a.engine, teff=a.teff, logg_min=a.logg_min,
                                                       a_step=a.a_step, models=len(order),
                                                       nodes=sum(len(sel[m]['abundances']) for m in order)), indent=1))
    print(f'{a.engine}: {len(order)} models, {len(todo)} to do, {a.workers} workers', flush=True)
    log = open(out / 'progress.log', 'a')
    with ProcessPoolExecutor(a.workers) as ex:
        futs = [ex.submit(_one, (m, sel[m], a.engine, str(out), a.workdir)) for m in todo]
        for k, f in enumerate(as_completed(futs), 1):
            try:
                mid, nok, n, dt = f.result()
                log.write(f'{time.strftime("%H:%M:%S")} {k}/{len(todo)} {mid} ok {nok}/{n} {dt:.0f}s\n')
            except Exception as e:                                     # noqa: BLE001
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
    print(G.assemble(a.npz, a.bin, a.aux_out, header_text=a.header))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    s = p.add_subparsers(dest='cmd', required=True)
    r = s.add_parser('run'); r.add_argument('--engine', choices=sorted(G.ENGINES), required=True)
    r.add_argument('--out', required=True); r.add_argument('--aux'); r.add_argument('--teff', default='0:99999')
    r.add_argument('--logg-min', type=float, default=-9.); r.add_argument('--a-step', type=float)
    r.add_argument('--workers', type=int, default=4); r.add_argument('--shard'); r.add_argument('--limit', type=int)
    r.add_argument('--workdir')
    t = s.add_parser('status'); t.add_argument('--out', required=True)
    m = s.add_parser('assemble'); m.add_argument('--npz', required=True); m.add_argument('--bin', required=True)
    m.add_argument('--aux-out', required=True); m.add_argument('--header', required=True)
    a = p.parse_args()
    dict(run=cmd_run, status=cmd_status, assemble=cmd_assemble)[a.cmd](a)


if __name__ == '__main__':
    main()

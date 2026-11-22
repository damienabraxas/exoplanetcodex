"""RYA-1395 v76: correlation-layer sizing for the O I MCDHF-EOL stage (Sirius; native rcsfgenerate ONLY).

Declared on RYA-1395 before build. MR = our s0 final lists (rya1354_mcdhf_s0_sequential, final 'ci' step);
SD substitutions from 2p and the outer electron (1s, 2s closed; Li et al. 2023 Sect. 2); J ranges of the s0 lists.
Layer k = 1..6: s,p <= 6+k; d,f <= 5+k; g <= 3+k (k >= 2); h <= 1+k (k >= 5) -> k=6 is Li's {12s,12p,11d,11f,9g,7h}.
Per-block cost is EXTRAPOLATED from the v52 measurement (odd J=1, 33 575 CSFs: rangular 282 s, rci 495 s; ~N^2.2).
CSF lists over KEEP_BYTES are deleted after the streamed inventory (sha256 + stdin kept: regenerable).
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1354_mcdhf_s0_sequential as S  # noqa: E402
from rya1354_li_scheme_sizing import inventory_stream  # noqa: E402
from rya1354_solver_smoke import run, sha  # noqa: E402

ANCHOR_N, ANCHOR_RANG_S, ANCHOR_RCI_S, POWER = 33575, 282.0, 495.0, 2.2
KEEP_BYTES = 200 * 2 ** 20


def active_set(k):
    top = {'s': 6 + k, 'p': 6 + k, 'd': 5 + k, 'f': 5 + k}
    if k >= 2:
        top['g'] = 3 + k
    if k >= 5:
        top['h'] = 1 + k
    return ','.join(f'{top[l]}{l}' for l in 'spdfgh' if l in top)


def mr(parity):
    k = len(S.STEPS[parity]) - 1
    fams = S.included(parity, k)
    return S.config_lines(parity, fams), S.twoj_range(fams)


def stdin_for(parity, k):
    lines, twoj = mr(parity)
    return '*\n0\n' + '\n'.join(lines) + '\n\n' + active_set(k) + '\n' + twoj + '\n2\nn\n'


def cost(n):
    f = (n / ANCHOR_N) ** POWER
    return dict(rangular_s=round(ANCHOR_RANG_S * f), rci_s=round(ANCHOR_RCI_S * f))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--layers', type=int, nargs='+', default=[1, 2, 3, 4, 5, 6])
    a = ap.parse_args()
    if platform.node().lower().split('.')[0] != 'sirius':
        raise RuntimeError('Atomic tools belong on Sirius')
    a.out = a.out.resolve()
    a.out.mkdir(parents=True, exist_ok=False)
    receipt = dict(schema='rya1395.corr-layer-sizing.v76', status='RUNNING', scientific_admission='HOLD',
                   started_utc=datetime.now(timezone.utc).isoformat(), cost_anchor=dict(
                       N=ANCHOR_N, rangular_s=ANCHOR_RANG_S, rci_s=ANCHOR_RCI_S, power=POWER, source='v52 odd J=1 n5'),
                   cases={})
    try:
        for k in a.layers:
            for parity in ('even', 'odd'):
                name = f'{parity}_k{k}'
                out = a.out / name
                out.mkdir()
                stdin = stdin_for(parity, k)
                step = run(a.grasp_bin / 'rcsfgenerate', stdin, out, 7200)
                inv = inventory_stream(out / 'rcsf.out', parity, 'mcdhf')
                size = (out / 'rcsf.out').stat().st_size
                case = dict(active_set=active_set(k), twoJ=mr(parity)[1], wall_s=step['wall_s'],
                            rcsf_bytes=size, rcsf_sha256=sha(out / 'rcsf.out'), **inv,
                            block_cost_extrapolated={b: cost(n) for b, n in inv['blocks'].items()})
                if size > KEEP_BYTES:
                    (out / 'rcsf.out').unlink()
                    case['rcsf_deleted_regenerable'] = True
                receipt['cases'][name] = case
                (a.out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
                big = max(inv['blocks'].values())
                print(f'{name:9s} {active_set(k):28s} CSFs {inv["CSFs"]:>9d}  largest block {big:>8d}  '
                      f'~rangular {cost(big)["rangular_s"]/3600:.1f} h  ~rci {cost(big)["rci_s"]/3600:.1f} h', flush=True)
        receipt['status'] = 'SIZING_COMPLETED_NO_ATOMIC_PREDICTION'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        (a.out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()

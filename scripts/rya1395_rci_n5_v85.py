"""RYA-1395 v85: weight-selected RCI on the converged n5 model (Sirius). HOLD, report only.

Declared on RYA-1395 before launch, with the v84 sizes (runs/rci_sizing_v84, comment cc64ba31).

Per parity, on the v77 n5 orbitals (rwfn.out, read only; byte-identical to state.w):
    v84 DC-reduced RCI list (1s closed, 2s opened, SD from n=2 + outer into the n5 active set; MR = the
    weight > 0.05 selection, which equals the v77 MCDHF MR; even MR + 2s2p4 ns/nd) -> rangular -> rci
    (Dirac-Coulomb only, as the v77 MCDHF; the v77 targeted levels per J) -> jj2lsj.
Then rbiotransform + rtransition E1 from the CI mixing (.cm) and the Li 2023 Table A2 report.

Gates (the first failure stops the run; it is kept, FAILED_RETAINED):
  (a) each RCI list is the v84 DC-reduced list, byte-identical (sha256 against a fresh hash of the v84 file);
  (b) disk: free space before rangular >= 1.5 x the v84 angular estimate for that parity + 10 GB;
  (c) identity: the leading physical families per J equal the targeted Li states (the v77 gate, E.identity);
  (d) rbiotransform reports normal orbital ordering;
  (e) E1 complete and finite: every v77 n5 transition present, no extra, all log gf finite.
A native command over its budget is a MEASURED stop (kept, reviewed), not a failure.
Reported, not gated: RCI vs v77 n5 MCDHF energies, labels and weights; level energies vs NIST
(Li Table A1); E1 vs v77 n5 and vs Li (by class, both gauges, energy-scaled gauge gap, 777 / key lines).
One niced process; parities in sequence.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import sys
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1354_mcdhf_eol as E  # noqa: E402
import rya1395_order_control_v83 as V  # noqa: E402

SCHEMA = 'rya1395.rci-n5-v85'
BUDGET = {'rangular': 43200, 'rci': 172800, 'jj2lsj': 3600}
V.BUDGET = 172800        # rbiotransform / rtransition inside V.e1 (estimate ~16 h for rtransition)
DISK_MARGIN_GB, DISK_FACTOR = 10, 1.5


def build(a, receipt, parity, sizing):
    out = a.out / parity
    out.mkdir()
    case = receipt['stages'][parity] = dict(steps={})
    src = a.v77 / f'n5_{parity}'
    v84_list = a.v84 / parity / 'DC' / 'rcsf.out'
    shutil.copyfile(v84_list, out / 'rcsf.inp')
    case['list_sha256'] = E.sha(out / 'rcsf.inp')
    if case['list_sha256'] != E.sha(v84_list):
        raise ValueError(f'gate (a) {parity}: RCI list differs from the v84 DC-reduced list')
    case['blocks'] = {k: b['dc_reduced'] for k, b in sizing[parity]['blocks'].items()}
    shutil.copyfile(src / 'isodata', out / 'isodata')
    if E.sha(src / 'rwfn.out') != E.sha(src / 'state.w'):
        raise ValueError(f'{src}: rwfn.out and state.w differ')
    shutil.copyfile(src / 'rwfn.out', out / 'state.w')
    case['orbitals_sha256'] = E.sha(out / 'state.w')
    need = DISK_FACTOR * sum(b['estimate']['angular_gb'] for b in sizing[parity]['blocks'].values()) + DISK_MARGIN_GB
    free = shutil.disk_usage(out).free / 1024 ** 3
    case['disk'] = dict(free_gb=round(free, 1), required_gb=round(need, 1))
    if free < need:
        raise ValueError(f'gate (b) {parity}: {free:.1f} GB free < {need:.1f} GB required')

    def step(exe, stdin):
        case['steps'][exe] = E.timed(a.grasp_bin / exe, stdin, out, BUDGET[exe], E.MEM_CAP_BYTES)
        V.save(a, receipt)

    step('rangular', 'y\n')
    os.link(out / 'rcsf.inp', out / 'state.c')
    lv = E.levels_stdin(E.block_order(out / 'rcsf.inp'), parity)
    step('rci', f'y\nstate\n{V.RCI_PHYSICS}{lv}')
    step('jj2lsj', V.JJ2LSJ_CI)
    case['leading_families'] = E.identity((out / 'state.lsj.lbl').read_text(), parity)
    case['vs_v77_mcdhf'] = V.compare_levels(src / 'state.lsj.lbl', out / 'state.lsj.lbl', gate=False)
    case['energy_shift_vs_v77'] = shifts(src / 'state.lsj.lbl', out / 'state.lsj.lbl')
    angular = sum(p.stat().st_size for p in out.glob('mcp.*'))
    case['angular_gb_measured'] = round(angular / 1024 ** 3, 2)
    V.save(a, receipt)
    return out


def shifts(ref, got):
    """Per level E(RCI) - E(MCDHF) in cm-1; the spread is what moves excitation energies."""
    ha = 219474.6313632
    r, g = V.levels(ref), V.levels(got)
    d = {k: (g[k][0] - r[k][0]) * ha for k in r}
    return dict(min_cm=min(d.values()), max_cm=max(d.values()), median_cm=statistics.median(d.values()),
                spread_cm=max(d.values()) - min(d.values()))


def summary(path):
    rows = list(csv.DictReader(open(path, newline='')))
    out = {}
    for c in 'ABCDE':
        rs = [r for r in rows if r['li_acc'] == c]
        if not rs:
            continue
        out[c] = dict(n=len(rs),
                      median_dB=statistics.median(float(r['d_B']) for r in rs),
                      median_abs_dB=statistics.median(abs(float(r['d_B'])) for r in rs),
                      median_dC=statistics.median(float(r['loggf_C']) + 2 * math.log10(float(r['sigma_model']) / float(r['sigma_nist']))
                                                  - float(r['li_loggf_B']) for r in rs),
                      gauge_gap_scaled=statistics.median(V.gap(r) for r in rs))
    key = [dict(nm=round(float(r['lam_vac_nm']), 3), line=f"{r['lower']} -> {r['upper']}", d_B=round(float(r['d_B']), 4),
                loggf_B=round(float(r['loggf_B']), 4), li=float(r['li_loggf_B']))
           for r in rows if any(abs(float(r['lam_vac_nm']) - w) < 0.6 for w in V.KEY_NM)]
    return dict(by_class=out, key_lines=sorted(key, key=lambda x: x['nm']))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin', type=Path, required=True)
    ap.add_argument('--v77', type=Path, required=True)
    ap.add_argument('--v84', type=Path, required=True)
    ap.add_argument('--nist-table', type=Path, required=True)
    ap.add_argument('--li-table', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if platform.node().lower().split('.')[0] != 'sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    for k, v in vars(a).items():
        if isinstance(v, Path):
            setattr(a, k, v.resolve())
    sizing = json.loads((a.v84 / 'receipt.json').read_text())
    if sizing['status'] != 'COMPLETE':
        raise ValueError('v84 sizing is not COMPLETE')
    a.out.mkdir(parents=True, exist_ok=False)
    for f in (__file__, E.__file__, V.__file__, V.REP.__file__):
        shutil.copyfile(f, a.out / Path(f).name)
    receipt = dict(schema=SCHEMA, status='RUNNING', scientific_admission='HOLD', adoption_allowed=False,
                   declared='RYA-1395 (v85 declaration)', started_utc=datetime.now(timezone.utc).isoformat(),
                   budgets_s=dict(BUDGET, rbiotransform=V.BUDGET, rtransition=V.BUDGET),
                   v84_estimate_totals=sizing['estimate_totals'], stages={}, transitions={}, comparisons={}, stopped=None)
    try:
        built = {p: build(a, receipt, p, sizing) for p in ('even', 'odd')}
        rci_csv = V.e1(a, receipt, 'n5', 'native', built['even'], built['odd'])
        ref_csv = a.v77 / 'n5_E1' / 'loggf_vs_li2023.csv'
        cmp = receipt['comparisons'] = dict(rci_vs_v77_mcdhf=V.paired(ref_csv, rci_csv),
                                            v77_mcdhf_vs_li=summary(ref_csv), rci_vs_li=summary(rci_csv))
        V.save(a, receipt)
        r = cmp['rci_vs_v77_mcdhf']
        if r['only_ref'] or r['only_got'] or not all(map(math.isfinite, (r['max_abs_dB'], r['max_abs_dC']))):
            raise ValueError(f'gate (e): RCI E1 incomplete or non-finite: {r["shared"]} shared, {r["only_ref"]} v77-only, '
                             f'{r["only_got"]} RCI-only, max|dB| {r["max_abs_dB"]}')
        receipt['status'] = 'RCI_COMPLETED_REVIEW_REQUIRED'
    except E.BudgetExceeded as exc:
        receipt.update(status='MEASURED_STOP_REVIEW_REQUIRED', stopped=str(exc))
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        receipt['files'] = {str(p.relative_to(a.out)): dict(sha256=E.sha(p), bytes=p.stat().st_size)
                            for p in sorted(a.out.rglob('*'))
                            if p.is_file() and p.name != 'receipt.json' and not p.name.startswith('mcp.')}
        V.save(a, receipt)


if __name__ == '__main__':
    main()

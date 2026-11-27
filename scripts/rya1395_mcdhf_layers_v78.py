"""RYA-1395 v78: resume the O I n6 correlation layer after the v77 MEASURED_STOP (Sirius).

Declared on RYA-1395 before build. v77 n6_even hit the declared 24 h rmcdhf budget at iteration 15 (residuals
2e-8..1.4e-6). Here:
  1. n6 even: CSF list, isodata and angular files hard-linked from v77 n6_even (identical list; rangular not redone),
     v77 iteration-15 rwfn.out as rwfn.inp, rmcdhf with the v77 stdin (levels, varied 6s/6d/6f/6g/6h, method 3,
     strict ACCY), then jj2lsj + identity gate;
  2. n6 odd: a normal v77 stage from v77 n5_odd;
  3. n6 E1 + Li 2023 report.
rmcdhf budget 72 h per parity; v77 floor admission unchanged; first failure stops and is kept.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1395_mcdhf_layers_v77 as V77  # noqa: E402  (patches the v53 driver: grid, strict, floor rule)

E, S, REP = V77.E, V77.S, V77.REP
SCHEMA = 'rya1395.mcdhf-layers-v78'
E.BUDGET_S = dict(E.BUDGET_S, rmcdhf=72 * 3600)


def resume_even(a, receipt, v77_dir):
    name = 'n6_even'
    src, out = v77_dir / name, a.out / name
    out.mkdir()
    case = dict(parity='even', layer=6, resumed_from=str(src / 'rwfn.out'), steps={})
    receipt['stages'][name] = case
    case['disk_free_before'] = E.disk_guard(a.out)
    links = ['rcsf.out', 'rcsf.inp', 'isodata'] + sorted(p.name for p in src.glob('mcp.*'))
    for f in links:
        os.link(src / f, out / f)
    shutil.copyfile(src / 'rwfn.out', out / 'rwfn.inp')
    case['linked'] = links
    case['rwfn_inp_sha256'] = E.sha(out / 'rwfn.inp')
    case.update(active_set=E.active_set('even', 6), **E.inventory_stream(out / 'rcsf.out', 'even', 'mcdhf'))
    levels = E.levels_stdin(E.block_order(out / 'rcsf.out'), 'even')
    varied = ','.join(E.new_orbitals('even', 6))
    stdin = V77.rmcdhf_stdin(levels, varied, '')
    if stdin != (src / 'rmcdhf.stdin').read_text():
        raise ValueError('rmcdhf stdin differs from v77 n6_even')
    case['varied'] = varied
    case['steps']['rmcdhf'] = E.timed(a.strict_bin / 'rmcdhf', stdin, out, E.BUDGET_S['rmcdhf'], E.MEM_CAP_BYTES)
    diag = E.scf_diagnostics((out / 'rmcdhf.stdout').read_text(), (out / 'rmcdhf.stderr').read_text())
    case['scf'] = diag
    if diag['iteration_limit_hit'] or not diag['all_reported_functions_below_ACC']:
        raise ValueError('Native orbital criterion unmet in n6_even (resume)')
    for s_, ext in (('rcsf.inp', 'c'), ('rwfn.out', 'w'), ('rmix.out', 'm')):
        os.link(out / s_, out / f'state.{ext}')
    case['steps']['jj2lsj'] = E.timed(a.grasp_bin / 'jj2lsj', E.JJ2LSJ_STDIN, out, E.BUDGET_S['jj2lsj'], E.MEM_CAP_BYTES)
    case['leading_families'] = E.identity((out / 'state.lsj.lbl').read_text(), 'even')
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin', type=Path, required=True)
    ap.add_argument('--strict-bin', type=Path, required=True)
    ap.add_argument('--v77', type=Path, required=True)
    ap.add_argument('--nist-table', type=Path, required=True)
    ap.add_argument('--li-table', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if platform.node().lower().split('.')[0] != 'sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    for k, v in vars(a).items():
        if isinstance(v, Path):
            setattr(a, k, v.resolve())
    a.v47 = None
    a.out.mkdir(parents=True, exist_ok=False)
    for f in (__file__, V77.__file__, E.__file__, S.__file__, REP.__file__):
        shutil.copyfile(f, a.out / Path(f).name)
    receipt = dict(schema=SCHEMA, status='RUNNING', scientific_admission='HOLD', adoption_allowed=False,
                   started_utc=datetime.now(timezone.utc).isoformat(), v77=str(a.v77), budgets_s=E.BUDGET_S,
                   stages={}, transitions={}, stopped=None)
    try:
        even = resume_even(a, receipt, a.v77)
        (a.out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        odd = E.stage(a, receipt, 'odd', 6, a.v77 / 'n5_odd')
        (a.out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        key = 'n6_E1'
        receipt['transitions'][key] = E.transitions(a.grasp_bin, even, odd, a.out / key, 'E1', False)
        REP.main(a.out / key, a.nist_table, a.li_table, a.out / key / 'loggf_vs_li2023.csv')
        receipt['last_completed_stage'] = key
        receipt['status'] = 'N6_COMPLETED_REVIEW_REQUIRED'
    except E.BudgetExceeded as exc:
        receipt.update(status='MEASURED_STOP_REVIEW_REQUIRED', stopped=str(exc))
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        (a.out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()

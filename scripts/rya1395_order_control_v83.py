"""RYA-1395 v83: native vs REVERSED orbital order through GRASP's own path (Sirius only). HOLD, report only.

Declared on RYA-1395 (comment 0f981784, 11-26) before build.

GRASP rbiotransform (KAPDATA) takes the orbital order from the shell order of the CSF lists and accepts only
all-normal or all-reversed order, identical in both parities. The native way to obtain reversed order is
rcsfgenerate's first prompt ('*' -> 'r'). The mixing vectors must then come from a CI in that basis, so for
every layer and parity this driver runs, on the CONVERGED native orbitals (read only):

    rcsfgenerate (v77/v78 stdin, first answer '*' or 'r') -> rangular -> rci (Dirac-Coulomb only,
    the rmcdhf targeted levels) -> jj2lsj

in two arms that differ ONLY in that first answer: 'native' ('*') and 'reverse' ('r'). Then rbiotransform +
rtransition E1 (CI mixing) per arm, and the Li 2023 Table A2 report.

Gates, reverse arm against the native arm of the same layer (the first failure stops the run; it is kept):
  (a) CSF count per J block and the multiset of configurations equal;
  (b) rbiotransform reports 'Reverse orbital ordering' (native arm: 'Normal orbital ordering');
  (c) v83b: level identity: each reversed level within 10% of its nearest same-J gap (a swap needs 50%) of the native level
      (an absolute 1e-8/1e-7 Ha bound failed on 5e-9-relative round-off; max |dE| is reported);
  (d) v83b: REPORTED only: canonical leading labels and weights (a reversed coupling tree is another LS basis);
  (e) n3 and n5 only: Li class A/B and key lines |d log gf| <= 0.005 dex in both gauges.
At n6, E1 is the measurement (reverse vs native vs v79), not a gate.
Reported, not gated: RCI vs the native rmcdhf energies, labels and E1.

One niced process; layers and arms in sequence; 24 h budget per native command.
"""
import argparse
import csv
import json
import math
import os
import platform
import re
import shutil
import statistics
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1354_mcdhf_eol as E  # noqa: E402
import rya1395_s0_e1_report as REP  # noqa: E402

SCHEMA = 'rya1395.order-control-v83b'
BUDGET = 86400
RCI_PHYSICS = 'n\nn\nn\nn\nn\nn\n'          # transverse, photon freq, VP, NMS, SMS, self-energy: Dirac-Coulomb only
JJ2LSJ_CI = 'state\ny\nn\ny\n'               # as E.JJ2LSJ_STDIN, but mixing from the CI (.cm)
GAP_FRACTION, GF_TOL = 0.1, 0.005   # v83b: gate (c) is level identity: |dE| < 10% of the nearest same-J gap (a swap needs >= 50%)
KEY_NM = (130.4, 135.6, 615.8, 645.6, 777.4, 844.8, 926.6, 1128.9, 1316.4)
ORBITAL = re.compile(r'(\d+[a-z]-?)\s*\(\s*(\d+)\)')
ORDER_MSG = {'native': 'Normal orbital ordering', 'reverse': 'Reverse orbital ordering'}


def blocks(path):
    """{J-block index: Counter(sorted configuration)} from a native CSF list (3-line records, '*' separators)."""
    out, block, buf = {}, 0, []
    with open(path) as f:
        for _ in range(5):
            next(f)
        for line in f:
            if line.lstrip().startswith('*'):
                block, buf = block + 1, []
                continue
            buf.append(line)
            if len(buf) == 3:
                out.setdefault(block, Counter())[tuple(sorted(ORBITAL.findall(buf[0])))] += 1
                buf = []
    return out


def gate_csfs(native, reverse):
    a, b = blocks(native), blocks(reverse)
    counts = {k: sum(v.values()) for k, v in a.items()}
    if sorted(map(sorted, (c.items() for c in a.values()))) != sorted(map(sorted, (c.items() for c in b.values()))):
        raise ValueError(f'gate (a): configuration content differs: {counts} vs '
                         f'{ {k: sum(v.values()) for k, v in b.items()} }')
    return dict(blocks=len(a), csfs=sum(counts.values()), per_block=counts)


def levels(path):
    """{(J, pos): (energy, leading label, leading weight)} from a jj2lsj .lsj.lbl."""
    out, lines = {}, Path(path).read_text().splitlines()
    for i, line in enumerate(lines):
        m = E.LEVEL.match(line)
        if m:
            comp = lines[i + 1].split()
            out[int(m[2]), int(m[1])] = (float(m[4]), comp[-1], float(comp[1]))
    return out


def canonical(label):
    """Coupling-order-free label: sorted subshell occupations with their own terms, plus the final term.
    Native '1s(2).2s(2).2p(3)4S.3p_3P' and reversed '3p.2p(3)4S_3P.2s(2).1s(2)_3P' both give
    (('1s(2)', '2p(3)', '2s(2)', '3p'), '3P'); subshell and intermediate terms depend on the order and are dropped."""
    body, final = label.rsplit('_', 1)
    return tuple(sorted(re.match(r'\d+[a-z](\(\d+\))?', t)[0] for t in body.split('.'))), final


def compare_levels(ref, got, gate):
    """Levels are identified by (J, position). Gate (c), energies, is the identity proof. Labels and
    weights are reported only, because a reversed coupling tree is a different LS basis (v83b)."""
    a, b = levels(ref), levels(got)
    if a.keys() != b.keys():
        raise ValueError(f'level sets differ: {sorted(a.keys() ^ b.keys())}')
    de = max(abs(a[k][0] - b[k][0]) for k in a)
    dw = max(abs(a[k][2] - b[k][2]) for k in a)
    labels = [k for k in a if canonical(a[k][1]) != canonical(b[k][1])]
    ratio = {}
    for k in a:
        gaps = [abs(a[k][0] - a[o][0]) for o in a if o[0] == k[0] and o != k]
        ratio[k] = abs(a[k][0] - b[k][0]) / min(gaps) if gaps else 0.0
    worst = max(ratio, key=ratio.get)
    res = dict(levels=len(a), max_dE_Ha=de, median_dE_Ha=statistics.median(abs(a[k][0] - b[k][0]) for k in a),
               max_dweight=dw, max_dE_over_gap=ratio[worst], worst_level=str(worst),
               canonical_label_mismatches=[str(k) for k in labels])
    if gate and ratio[worst] >= GAP_FRACTION:
        raise ValueError(f'gate (c): level {worst} |dE| is {ratio[worst]:.3g} of its nearest same-J gap (>= {GAP_FRACTION})')
    return res


def build(a, receipt, layer, arm, parity, source, gen_src):
    name = f'{layer}_{arm}_{parity}'
    out = a.out / name
    out.mkdir()
    case = receipt['stages'][name] = dict(source=str(source), generator_source=str(gen_src), steps={})
    if (gen_src / 'rcsf.out').read_bytes() != (source / 'rcsf.out').read_bytes():
        raise ValueError(f'{source}: CSF list differs from its generator source {gen_src}')

    def step(exe, stdin):
        case['steps'][exe] = E.timed(a.grasp_bin / exe, stdin, out, BUDGET, E.MEM_CAP_BYTES)

    gen = (gen_src / 'rcsfgenerate.stdin').read_text()
    if not gen.startswith('*\n'):
        raise ValueError(f'{source}: rcsfgenerate stdin does not start with the default-order answer')
    step('rcsfgenerate', ('r' if arm == 'reverse' else '*') + gen[1:])
    case['csfs'] = gate_csfs(source / 'rcsf.out', out / 'rcsf.out')
    case['header'] = (out / 'rcsf.out').read_text().splitlines()[3].split()
    shutil.copyfile(source / 'isodata', out / 'isodata')
    os.link(out / 'rcsf.out', out / 'rcsf.inp')
    step('rangular', 'y\n')
    os.link(out / 'rcsf.inp', out / 'state.c')
    shutil.copyfile(source / 'rwfn.out', out / 'state.w')
    case['orbitals_sha256'] = E.sha(out / 'state.w')
    lv = E.levels_stdin(E.block_order(out / 'rcsf.out'), parity)
    step('rci', f'y\nstate\n{RCI_PHYSICS}{lv}')
    step('jj2lsj', JJ2LSJ_CI)
    if arm == 'native':
        case['leading_families'] = E.identity((out / 'state.lsj.lbl').read_text(), parity)
        case['vs_rmcdhf'] = compare_levels(source / 'state.lsj.lbl', out / 'state.lsj.lbl', gate=False)
    else:
        # v83b: jj2lsj writes reversed-order labels in reversed coupling order, which the v77 identity parser
        # cannot read. The reversed levels are matched to the native arm by (J, position) under gate (c);
        # the native arm's labels (identity-gated) then label them for the E1 report.
        native = a.out / f'{layer}_native_{parity}' / 'state.lsj.lbl'
        os.rename(out / 'state.lsj.lbl', out / 'state.reverse.lsj.lbl')
        case['gate_c'] = compare_levels(native, out / 'state.reverse.lsj.lbl', gate=True)
        shutil.copyfile(native, out / 'state.lsj.lbl')
    save(a, receipt)
    return out


def e1(a, receipt, layer, arm, even, odd):
    out = a.out / f'{layer}_{arm}_E1'
    out.mkdir()
    for tag, src in (('lower', even), ('upper', odd)):
        for ext in ('c', 'cm', 'w'):
            shutil.copyfile(src / f'state.{ext}', out / f'{tag}.{ext}')
        shutil.copyfile(src / 'state.lsj.lbl', out / f'{tag}.lsj.lbl')
    shutil.copyfile(even / 'isodata', out / 'isodata')
    steps = {}
    for exe, stdin in (('rbiotransform', 'y\ny\nlower\nupper\ny\n'), ('rtransition', 'y\ny\nlower\nupper\nE1\n')):
        steps[exe] = E.timed(a.grasp_bin / exe, stdin, out, BUDGET, E.MEM_CAP_BYTES)
    log = (out / 'rbiotransform.stdout').read_text()
    if ORDER_MSG[arm] not in log:
        raise ValueError(f'gate (b): {out.name} rbiotransform did not report "{ORDER_MSG[arm]}"')
    table = out / 'lower.upper.ct'
    if not table.is_file() or table.stat().st_size == 0:
        raise ValueError(f'{out.name}: empty or missing CI transition table')
    os.link(table, out / 'lower.upper.t')      # the report reads lower.upper.t; same bytes
    REP.main(out, a.nist_table, a.li_table, out / 'loggf_vs_li2023.csv')
    receipt['transitions'][out.name] = dict(steps=steps, order=ORDER_MSG[arm])
    save(a, receipt)
    return out / 'loggf_vs_li2023.csv'


def load(path):
    with open(path, newline='') as f:
        return {(r['lower'], r['upper']): r for r in csv.DictReader(f)}


def gap(r):
    return abs(float(r['loggf_B']) - float(r['loggf_C']) - 2 * math.log10(float(r['sigma_model']) / float(r['sigma_nist'])))


def paired(ref_csv, got_csv):
    ref, got = load(ref_csv), load(got_csv)
    keys = ref.keys() & got.keys()
    dB = {k: float(got[k]['loggf_B']) - float(ref[k]['loggf_B']) for k in keys}
    dC = {k: float(got[k]['loggf_C']) - float(ref[k]['loggf_C']) for k in keys}
    graded = [k for k in keys if ref[k]['li_acc'] in ('A', 'B')
              or any(abs(float(ref[k]['lam_vac_nm']) - w) < 0.6 for w in KEY_NM)]
    gaps = {c: statistics.median([gap(got[k]) for k in keys if got[k]['li_acc'] == c] or [float('nan')])
            for c in 'ABCDE'}
    key777 = sorted((k for k in keys if abs(float(ref[k]['lam_vac_nm']) - 777.4) < 0.6),
                    key=lambda k: float(ref[k]['lam_vac_nm']))
    return dict(shared=len(keys), only_ref=len(ref.keys() - keys), only_got=len(got.keys() - keys),
                max_abs_dB=max(map(abs, dB.values())), max_abs_dC=max(map(abs, dC.values())),
                median_abs_dB=statistics.median(map(abs, dB.values())),
                graded_lines=len(graded),
                graded_max_abs_dB=max(abs(dB[k]) for k in graded), graded_max_abs_dC=max(abs(dC[k]) for k in graded),
                energy_scaled_gauge_gap_median_by_class=gaps,
                line_777=[dict(line=f'{k[0]} -> {k[1]}', ref_B=float(ref[k]['loggf_B']), got_B=float(got[k]['loggf_B']))
                          for k in key777])


def generator(a, layer, parity):
    """v78 n6 even resumed v77's n6 even and hard-linked its CSF list; its rcsfgenerate stdin is in v77."""
    d = (a.v78 if layer == 'n6' else a.v77) / f'{layer}_{parity}'
    return d if (d / 'rcsfgenerate.stdin').is_file() else a.v77 / f'{layer}_{parity}'


def save(a, receipt):
    (a.out / 'receipt.json').write_text(json.dumps(receipt, indent=2, default=str) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin', type=Path, required=True)
    ap.add_argument('--v77', type=Path, required=True, help='runs/mcdhf_layers_v77 (n3, n5 sources)')
    ap.add_argument('--v78', type=Path, required=True, help='runs/mcdhf_layers_v78 (n6 sources)')
    ap.add_argument('--v79', type=Path, required=True, help='runs/mcdhf_e1_v79 (native n6 E1 report)')
    ap.add_argument('--nist-table', type=Path, required=True)
    ap.add_argument('--li-table', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--layers', nargs='+', default=['n3', 'n5', 'n6'])
    a = ap.parse_args()
    if platform.node().lower().split('.')[0] != 'sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    for k, v in vars(a).items():
        if isinstance(v, Path):
            setattr(a, k, v.resolve())
    a.out.mkdir(parents=True, exist_ok=False)
    for f in (__file__, E.__file__, REP.__file__):
        shutil.copyfile(f, a.out / Path(f).name)
    sources = dict(n3=a.v77, n5=a.v77, n6=a.v78)
    receipt = dict(schema=SCHEMA, status='RUNNING', scientific_admission='HOLD', adoption_allowed=False,
                   declared='RYA-1395 comment 0f981784', started_utc=datetime.now(timezone.utc).isoformat(),
                   tolerances=dict(dE_over_nearest_gap=GAP_FRACTION, gf_dex=GF_TOL), stages={}, transitions={},
                   comparisons={}, stopped=None)
    try:
        for layer in a.layers:
            arms = {}
            for arm in ('native', 'reverse'):
                built = {p: build(a, receipt, layer, arm, p, sources[layer] / f'{layer}_{p}',
                                  generator(a, layer, p)) for p in ('even', 'odd')}
                arms[arm] = e1(a, receipt, layer, arm, built['even'], built['odd'])
            cmp = receipt['comparisons'][layer] = dict(
                reverse_vs_native_rci=paired(arms['native'], arms['reverse']),
                native_rci_vs_rmcdhf=paired((a.v79 if layer == 'n6' else a.v77 / f'{layer}_E1') / 'loggf_vs_li2023.csv',
                                            arms['native']))
            save(a, receipt)
            r = cmp['reverse_vs_native_rci']
            if r['only_ref'] or r['only_got'] or not all(map(math.isfinite, (r['max_abs_dB'], r['max_abs_dC']))):
                raise ValueError(f'gate (e) {layer}: reversed E1 not comparable: {r["shared"]} shared, '
                                 f'{r["only_ref"]} native-only, max|dB| {r["max_abs_dB"]} (NaN never fails a > test)')
            if layer != 'n6':
                if not max(r['graded_max_abs_dB'], r['graded_max_abs_dC']) <= GF_TOL:
                    raise ValueError(f'gate (e) {layer}: graded |d log gf| B {r["graded_max_abs_dB"]:.3g} '
                                     f'C {r["graded_max_abs_dC"]:.3g} > {GF_TOL}')
            receipt['last_completed_layer'] = layer
        receipt['status'] = 'ORDER_CONTROL_COMPLETED_REVIEW_REQUIRED'
    except E.BudgetExceeded as exc:
        receipt.update(status='MEASURED_STOP_REVIEW_REQUIRED', stopped=str(exc))
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        files = {}
        for p in sorted(a.out.rglob('*')):
            if p.is_file() and p.name != 'receipt.json' and not p.name.startswith('mcp.'):
                files[str(p.relative_to(a.out))] = dict(sha256=E.sha(p), bytes=p.stat().st_size)
        receipt['files'] = files
        save(a, receipt)


if __name__ == '__main__':
    main()

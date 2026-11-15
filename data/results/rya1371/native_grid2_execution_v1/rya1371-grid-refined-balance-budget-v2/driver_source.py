"""Localize saved continuum mesh changes in physical level-balance units; no RT/SE."""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def net(matrix, populations):
    return np.sum(matrix * populations[None, :, :], axis=1) - matrix.sum(axis=0) * populations


def main():
    from config.constants import require_codex_path
    from pipeline.multi_atom import read_multi_atom
    from pipeline.reference_engines import fingerprint

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    choice = p.add_mutually_exclusive_group()
    choice.add_argument('--refined-native-state', action='store_true')
    choice.add_argument('--second-refined-native-state', action='store_true')
    a = p.parse_args()
    checks = {str(Path(__file__)): fingerprint(__file__)}

    def load(group, filename):
        path = a.raw_root / ('rya1371-' + group) / filename
        checks[str(path)] = fingerprint(path)
        if path.suffix == '.json':
            doc = json.loads(path.read_text())
            if not doc['computational_completed']:
                raise ValueError('Incomplete producer')
            return doc
        with np.load(path) as packet:
            return {k: packet[k].copy() for k in packet.files}

    groups = ['background-midpoint-v1', 'background-adaptive-v1']
    audits = ['background-adaptive-plan-v1', 'background-final-audit-v1']
    if a.refined_native_state:
        groups = ['hop-grid-refined-midpoint-v1', 'hop-grid-refined-adaptive-v1']
        audits = ['grid-refined-midpoint-audit-v1', 'grid-refined-final-audit-v1']
    if a.second_refined_native_state:
        groups = ['hop-grid-refined2-midpoint-v1', 'hop-grid-refined2-adaptive-v1']
        audits = ['grid-refined2-midpoint-audit-v1', 'grid-refined2-final-audit-v1']
    data, high = [], []
    for group, audit in zip(groups, audits):
        doc = load(group, 'report.json')
        packet = load(group, 'rates.npz')
        audit_doc = load(audit, 'report.json')
        if audit_doc['inputs'].get(str(a.raw_root / ('rya1371-' + group) / 'rates.npz')) != checks[str(a.raw_root / ('rya1371-' + group) / 'rates.npz')]:
            raise ValueError('Audit belongs to a different rate packet')
        data.append(packet)
        high.append(load(audit, 'parabolic_rates.npz')['R_s_inv'])
    old, new = data
    for key in ['populations', 'nstar', 'replay_C_s_inv', 'native_R_s_inv']:
        np.testing.assert_array_equal(old[key], new[key])
    idx = np.searchsorted(new['wave_vac_nm'], old['wave_vac_nm'])
    np.testing.assert_array_equal(new['wave_vac_nm'][idx], old['wave_vac_nm'])
    np.testing.assert_array_equal(new['Jnu'][idx], old['Jnu'])
    atom_path = require_codex_path('experiments.bride_oxygen_atom')
    checks[str(atom_path)] = fingerprint(atom_path)
    checks[str(ROOT / 'pipeline/multi_atom.py')] = fingerprint(ROOT / 'pipeline/multi_atom.py')
    atom = read_multi_atom(atom_path)
    n = new['populations']
    C = new['replay_C_s_inv']
    R = high[-1]
    native = new['native_R_s_inv']
    Q = R + C
    gross = np.sum(Q * n[None, :, :], axis=1) + Q.sum(axis=0) * n
    if not np.isfinite(gross).all() or np.any(gross <= 0):
        raise ValueError('Invalid gross physical normalization')
    residual = net(Q, n) / gross
    peak = np.unravel_index(np.argmax(abs(residual)), residual.shape)
    mesh_delta = high[1] - high[0]
    native_delta = R - native
    bf_mesh = np.zeros_like(R)
    bf_native = np.zeros_like(R)
    pair_rows = []
    for t in atom.continua:
        i, j = t.i - 1, t.j - 1
        for dst, src in [(bf_mesh, mesh_delta), (bf_native, native_delta)]:
            dst[j, i] += src[j, i]
            dst[i, j] += src[i, j]
        # Net transfer into upper level; opposite transfer enters lower level.
        dm = mesh_delta[j, i] * n[i] - mesh_delta[i, j] * n[j]
        dn = native_delta[j, i] * n[i] - native_delta[i, j] * n[j]
        for depth in range(n.shape[1]):
            pair_rows.append(dict(lower_1based=i + 1, upper_1based=j + 1,
                depth_top_to_bottom=depth,
                up_relative_mesh_step=float(mesh_delta[j, i, depth] / high[0][j, i, depth]),
                down_relative_mesh_step=float(mesh_delta[i, j, depth] / high[0][i, j, depth]),
                signed_upper_mesh_balance_step=float(dm[depth] / gross[j, depth]),
                signed_lower_mesh_balance_step=float(-dm[depth] / gross[i, depth]),
                signed_upper_native_gap=float(dn[depth] / gross[j, depth]),
                signed_lower_native_gap=float(-dn[depth] / gross[i, depth])))
    # These are exact finite differences, not error bars or Richardson bounds.
    native_gap = net(native_delta, n) / gross
    bf_gap = net(bf_native, n) / gross
    bb_gap = native_gap - bf_gap
    mesh = net(mesh_delta, n) / gross
    bf_mesh_net = net(bf_mesh, n) / gross
    l1 = np.zeros_like(n)
    for row in pair_rows:
        d = row['depth_top_to_bottom']
        l1[row['lower_1based'] - 1, d] += abs(row['signed_lower_mesh_balance_step'])
        l1[row['upper_1based'] - 1, d] += abs(row['signed_upper_mesh_balance_step'])
    np.testing.assert_allclose(net(mesh_delta, n), net(high[1], n) - net(high[0], n), atol=float(gross.max()) * 5e-15, rtol=0)
    bb_mesh_net = mesh - bf_mesh_net
    np.testing.assert_allclose(np.sum(net(bf_native, n), axis=0), 0, atol=float(gross.max()) * 5e-15, rtol=0)
    level_rows = []
    for i in range(n.shape[0]):
        for d in range(n.shape[1]):
            level_rows.append(dict(level_1based=i + 1, depth_top_to_bottom=d,
                signed_residual=float(residual[i, d]),
                signed_native_rate_gap=float(native_gap[i, d]),
                signed_BF_native_gap=float(bf_gap[i, d]),
                signed_BB_native_gap=float(bb_gap[i, d]),
                signed_mesh_step=float(mesh[i, d]),
                signed_BF_mesh_step=float(bf_mesh_net[i, d]),
                signed_BB_mesh_step=float(bb_mesh_net[i, d]),
                BF_pair_L1_mesh_step=float(l1[i, d])))
    selected = level_rows[int(peak[0]) * n.shape[1] + int(peak[1])]
    relevant = [dict(r, signed_peak_native_gap=r['signed_upper_native_gap'] if r['upper_1based'] == peak[0] + 1 else r['signed_lower_native_gap'],
                     signed_peak_mesh_step=r['signed_upper_mesh_balance_step'] if r['upper_1based'] == peak[0] + 1 else r['signed_lower_mesh_balance_step'])
                for r in pair_rows if r['depth_top_to_bottom'] == peak[1] and peak[0] + 1 in [r['upper_1based'], r['lower_1based']]]
    relevant.sort(key=lambda r: -abs(r['signed_peak_native_gap']))
    report = dict(generator=fingerprint(__file__), inputs=checks,
        computational_completed=True, scientific_admission='HOLD',
        nodes=[len(x['wave_vac_nm']) for x in data],
        peak_balance=selected, peak_BF_pairs=relevant,
        max_L1_BF_mesh_balance_step=float(l1.max()),
        max_absolute_mesh_balance_step=float(abs(mesh).max()),
        max_absolute_BB_mesh_balance_step=float(abs(bb_mesh_net).max()),
        max_absolute_BF_native_balance_gap=float(abs(bf_gap).max()),
        max_absolute_BB_native_balance_gap=float(abs(bb_gap).max()),
        scope='Fixed exact n/nStar/C/native rates; inherited J bit-identical. Same-node positive parabolic BF with compiled BB from each corresponding grid. Both BF and BB mesh changes are accounted separately; adding continuum nodes also changes line sampling. All 22 continua/56 depths/41 levels accounted. One final-step difference is not a numerical error bound, convergence certificate or isolated physics cause. Normalization uses final own gross R+C throughout; native own residual uses a different denominator. No RT, SE, abundance or production change.')
    for path, fp in checks.items():
        if fingerprint(path) != fp:
            raise ValueError('Input changed during balance audit')
    a.output.mkdir(parents=True, exist_ok=False)
    for name, rows in [('pairs.csv', pair_rows), ('levels.csv', level_rows)]:
        with (a.output / name).open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
            w.writeheader()
            w.writerows(rows)
    (a.output / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ['inputs', 'peak_BF_pairs']}))


if __name__ == '__main__':
    main()

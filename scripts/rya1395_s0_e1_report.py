"""RYA-1395 v73 report: s0-stage E1 line strengths from our GRASP model vs Li 2023 (Table A2) and NIST.

Reads lower.upper.t + lower/upper .lsj.lbl from a GRASP transitions directory (file f1 = even, f2 = odd; rows run
upper -> lower in either direction),
labels every level by its leading LSJ component, and converts each line strength S (Babushkin B and Coulomb C)
to gf with NIST energies (Li Table A1 ENIST): gf = 2 sigma S / (3 * 219474.63), sigma in cm-1.
Report only: no grade, no adoption.

Usage: rya1395_s0_e1_report.py E1_DIR TABLE_A1 TABLE_A2 OUT_CSV
"""
import csv
import math
import re
import sys
from pathlib import Path

HC = 219474.6313632
LVL = re.compile(r'^\s*(\d+)\s+(\d+(?:/\d)?)\s+([+-])\s+(-?\d+\.\d+)\s+([\d.]+)%')
TROW = re.compile(r'^\s*(f[12])\s+(\d+)\s+(\d+)\s+([+-])\s+(f[12])\s+(\d+)\s+(\d+)\s+([+-])\s+([\d.]+)\s+([CB])\s+(\S+)\s+(\S+)\s+(\S+)')
BROW = re.compile(r'^\s+([CB])\s+(\S+)\s+(\S+)\s+(\S+)\s*$')


def fnum(s):
    return float(s.replace('D', 'E'))


def levels(path):
    out, lines = {}, Path(path).read_text().splitlines()
    for i, l in enumerate(lines):
        m = LVL.match(l)
        if m:
            out[(int(m[1]), int(m[2]), m[3])] = key_from_ours(lines[i + 1].split()[-1], int(m[2]))
    return out


def key_from_ours(label, j):
    """'1s(2).2s(2).2p(3)4S.3s_5S' -> ('3s', '5S', j); '1s(2).2s(2).2p(4)3P_3P' -> ('2p(4)', '3P', j)."""
    conf, term = label.rsplit('_', 1)
    last = conf.split('.')[-1]
    if last.startswith('2p(4)'):
        last = '2p(4)'
    return (re.sub(r'^\d?[A-Z]', '', last) if not last.startswith('2p') else last, term, j)


def key_from_li(conf, term, j):
    return (conf.strip().split('.')[-1], term.strip(), int(j))


def table_rows(path):
    rows = [r for r in csv.reader((l for l in open(path) if not l.startswith('#')), delimiter='\t')]
    return [r for r in rows if len(r) > 5 and r[0].strip().isdigit()]   # repeated header blocks occur mid-table


def main(e1_dir, a1, a2, out_csv):
    e1 = Path(e1_dir)
    lo, up = levels(e1 / 'lower.lsj.lbl'), levels(e1 / 'upper.lsj.lbl')
    nist = {key_from_li(r[2], r[3], r[4]): float(r[7]) for r in table_rows(a1)}
    li = {}
    for r in table_rows(a2):
        li[(key_from_li(r[1], r[2], r[3]), key_from_li(r[5], r[6], r[7]))] = dict(
            lam_nm=float(r[9]), loggf_B=float(r[12]), loggf_C=float(r[13]), dT=float(r[16]), acc=r[19].strip())
    text = (e1 / 'lower.upper.t').read_text().splitlines()
    res, cur = [], None
    for l in text:
        m = TROW.match(l)
        if m:
            files = {'f1': lo, 'f2': up}   # f1 = even (lower.*), f2 = odd (upper.*); either can hold the upper level
            ku = files[m[1]].get((int(m[2]), int(m[3]), m[4]))
            kl = files[m[5]].get((int(m[6]), int(m[7]), m[8]))
            cur = dict(upper=ku, lower=kl, sigma_model=float(m[9]))
            cur[f'S_{m[10]}'] = fnum(m[13])
            res.append(cur)
            continue
        b = BROW.match(l)
        if b and cur is not None:
            cur[f'S_{b[1]}'] = fnum(b[4])
    rows = []
    for t in res:
        if t['upper'] is None or t['lower'] is None:
            continue
        eu, el = nist.get(t['upper']), nist.get(t['lower'])
        if eu is None or el is None:
            continue
        sig = eu - el
        flipped = sig < 0   # separately optimised parities: the model can order a pair opposite to NIST (s0: 3d < 3p)
        if flipped:
            t['upper'], t['lower'], sig = t['lower'], t['upper'], -sig
        rec = dict(upper=' '.join(map(str, t['upper'])), lower=' '.join(map(str, t['lower'])), sigma_nist=sig,
                   lam_vac_nm=1e7 / sig, sigma_model=t['sigma_model'], model_order_flipped=flipped)
        for g in 'BC':
            gf = 2 * sig * t[f'S_{g}'] / (3 * HC)
            rec[f'loggf_{g}'] = math.log10(gf) if gf > 0 else float('nan')
        L = li.get((t['upper'], t['lower']))
        if L:
            rec.update(li_loggf_B=L['loggf_B'], li_acc=L['acc'], d_B=rec['loggf_B'] - L['loggf_B'],
                       gauge_dex=rec['loggf_B'] - rec['loggf_C'])
        rows.append(rec)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ('upper', 'lower', 'lam_vac_nm'), k))
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    matched = [r for r in rows if 'd_B' in r]
    print(f'{len(rows)} lines with NIST energies; {len(matched)} matched to Li Table A2')
    for acc in sorted({r['li_acc'] for r in matched}):
        d = sorted(r['d_B'] for r in matched if r['li_acc'] == acc)
        print(f"  Li class {acc:2s} n={len(d):3d}  median d {d[len(d)//2]:+.3f}  "
              f"abs68 {sorted(abs(x) for x in d)[int(0.68*(len(d)-1))]:.3f}")
    shown = set()
    for lam in (130.4, 135.6, 615.8, 645.6, 777.4, 844.8, 926.6, 1128.9, 1316.4):
        near = [r for r in matched if abs(r['lam_vac_nm'] - lam) < 1.0]
        for r in sorted(near, key=lambda r: r['lam_vac_nm'])[:4]:
            if id(r) in shown:
                continue
            shown.add(id(r))
            print(f"  {r['lam_vac_nm']:9.3f} nm  {r['upper']:12s} -> {r['lower']:12s}  ours B {r['loggf_B']:+.3f} "
                  f"C {r['loggf_C']:+.3f}  Li {r['li_loggf_B']:+.3f} ({r['li_acc']})  d {r['d_B']:+.3f}")


if __name__ == '__main__':
    main(*sys.argv[1:5])

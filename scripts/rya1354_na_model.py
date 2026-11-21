"""RYA-1354 v61: our own Na I atomic model (FAC), scored against experiment and independent theory,
then used to fill every Na I line in the solar line list with a log gf and a graded sigma.

Declared on RYA-1354 before the native run (2026-11-19). Status: CANDIDATE / HOLD (RYA-1352).

  run    (Sirius)  write FAC inputs for each variant x gauge, run the pinned sfac, keep lev/tr tables
  score  (any)     parse FAC + Froese Fischer (NIST MCHF collection [9] B-spline CI, [16] MCHF-CV),
                   correct to experimental energies, score, grade, fill the solar line list

Referees
  experiment : 3p 2P lifetimes; KP08 (Kelleher & Podobedova 2008, JPCRD 37, 267, Table 4) adopted D-line
               A-values 6.16e7 / 6.14e7 s-1 (class AA), from Volz+1996, Jones+1996, Oates+1996.
  theory     : Froese Fischer 2002 B-spline CI (what NIST adopts for every solar Na I line) and MCHF-CV.

Energy correction: computed line strengths are rescaled to NIST experimental level energies (standard practice;
length form gf ~ dE, velocity form gf ~ 1/dE).
"""
import argparse
import csv
import gzip
import json
import math
import re
import statistics as st
import subprocess
from collections import defaultdict
from pathlib import Path

import rya1354_theory_grade as G

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / 'data/reference/na_mchf_ff2002/raw'
OUT = ROOT / 'data/audit/rya1354/na_model_v61'
SFAC = '/mnt/codex-data/codex/rya1354/toolchains/fac-9ab6def/sfac/sfac'
NMAX = 9
L = 'spdf'
HC_EV_CM = 8065.543937  # cm-1 per eV
A_COEF = 6.6702e15      # A [s-1] = A_COEF * gf / (g_up * lambda_vac[A]^2)
EXP_D = {('3s', 1, '3p', 3): (6.16e7, 0.005), ('3s', 1, '3p', 1): (6.14e7, 0.005)}  # KP08 Table 4, AA
STATUS = 'CANDIDATE_HOLD_NOT_ADOPTED'

VARIANTS = {
    # v61a: valence only, Dirac-Fock-Slater potential optimised on the n=3 average
    'v61a': dict(closed='1s 2s 2p', configs=[(f'n{n}', f'{n}*1') for n in range(3, NMAX + 1)], opt=['n3']),
    # v61b: v61a + 2p5 core-excited perturbers (first attempt at core-valence correlation)
    'v61b': dict(closed='1s 2s', configs=[(f'n{n}', f'2p6 {n}*1') for n in range(3, NMAX + 1)]
                 + [('cv33', '2p5 3*2'), ('cv34', '2p5 3*1 4*1')], opt=['n3']),
}


def fac_input(v, gauge):
    s = VARIANTS[v]
    g = [c[0] for c in s['configs']]
    lines = ['SetAtom("Na")', f'Closed("{s["closed"]}")']
    lines += [f'Config("{n}","{c}")' for n, c in s['configs']]
    lines += [f'OptimizeRadial({json.dumps(s["opt"])})', f'Structure("t.lev.b",{json.dumps(g)})', 'MemENTable("t.lev.b")',
              'PrintTable("t.lev.b","t.lev",1)', f'SetTransitionOptions({gauge},0)',
              f'TransitionTable("t.tr.b",{json.dumps(g)},{json.dumps(g)})', 'PrintTable("t.tr.b","t.tr",1)']
    return '\n'.join(lines) + '\n'


def run(out):
    out = Path(out)
    for v in VARIANTS:
        for gname, gauge in (('B', 2), ('C', 1)):
            d = out / f'{v}_{gname}'
            d.mkdir(parents=True, exist_ok=True)
            (d / 'in.sf').write_text(fac_input(v, gauge))
            with open(d / 'sfac.stdout', 'w') as f:
                rc = subprocess.run(['nice', '-n', '19', SFAC, 'in.sf'], cwd=d, stdout=f, stderr=subprocess.STDOUT, timeout=3600).returncode
            print(v, gname, 'rc', rc, flush=True)
            if rc:
                raise SystemExit(f'{v}_{gname} failed rc={rc}')


# ---------- parsing ----------
def nist_levels():
    lev = {}
    for r in csv.DictReader(open(REF / 'nist_asd_na1_levels.csv')):
        clean = {k: (v or '').strip('="[]') for k, v in r.items() if k}
        m = re.fullmatch(r'2p6\.(\d+)([spdf])', clean['Configuration'])
        if not m or not re.fullmatch(r'\d+/2', clean['J']) or not clean['Level (cm-1)']:
            continue
        lev[(f'{m[1]}{m[2]}', int(clean['J'].split('/')[0]))] = float(clean['Level (cm-1)'])
    return lev


def parse_ff(name):
    """NIST MCHF .lsj.db: blocks of (2J E label) x2, wavenumber line, length row, velocity row."""
    txt = gzip.open(REF / name, 'rt').read().split('\n')
    out = {}
    i = 0
    while i < len(txt):
        a = txt[i].split()
        if len(a) == 3 and a[2].count('_') == 1 and i + 4 < len(txt) and 'E1' in txt[i + 3]:
            b = txt[i + 1].split()
            nl_lo, nl_up = a[2].split('_')[0], b[2].split('_')[0]
            ln = re.findall(r'[-+]?\d\.\d+e[-+]\d+', txt[i + 3])
            ve = re.findall(r'[-+]?\d\.\d+e[-+]\d+', txt[i + 4])
            lo, up = (nl_lo, int(a[0])), (nl_up, int(b[0]))
            e_lo, e_up = float(a[1]), float(b[1])
            if e_up < e_lo:
                lo, up = up, lo
            sig = float(txt[i + 2].split()[0])
            out[(lo[0], lo[1], up[0], up[1])] = dict(sigma_th=sig, gf_L=float(ln[1]), gf_V=float(ve[1]))
            i += 5
        else:
            i += 1
    return out


def fac_name(tok):
    """FAC level name like '3p1' (closed core) or '2p6 3p1' -> '3p'; core-excited levels -> None."""
    m = re.fullmatch(r'(?:2p6\s*)?(\d+[spdf])1', tok.strip())
    return m[1] if m else None


def _open(p):
    return gzip.open(str(p) + '.gz', 'rt') if not p.exists() and Path(str(p) + '.gz').exists() else open(p)


def parse_fac(d):
    lev = {}
    for line in _open(d / 't.lev'):
        # ILEV IBASE ENERGY P VNL 2J, then (2+ space separated) group/short name, nonrel config, rel config
        m = re.match(r'\s*(\d+)\s+(-?\d+)\s+(\S+)\s+(\d)\s+(\d+)\s+(\d+)\s+(.*)$', line)
        if not m:
            continue
        parts = re.split(r'\s{2,}', m[7].strip())
        if len(parts) < 3:
            continue
        lev[int(m[1])] = dict(E_eV=float(m[3]), twoJ=int(m[6]), nl=fac_name(parts[-2]))
    tr = []
    for line in _open(d / 't.tr'):
        f = line.split()
        if len(f) == 8 and f[0].isdigit():
            up, lo = lev[int(f[0])], lev[int(f[2])]
            if up['nl'] and lo['nl']:
                if up['E_eV'] < lo['E_eV']:
                    up, lo = lo, up
                tr.append(dict(lo=(lo['nl'], lo['twoJ']), up=(up['nl'], up['twoJ']), dE_eV=float(f[4]), gf=float(f[5])))
    return lev, tr


def corrected(gf, de_th_cm, de_exp_cm, gauge):
    return gf * (de_exp_cm / de_th_cm if gauge == 'B' else de_th_cm / de_exp_cm)


def model(v, levels):
    """Both gauges, energy-corrected, keyed (lo_nl, 2J, up_nl, 2J)."""
    res = {}
    for gname in 'BC':
        _, tr = parse_fac(OUT / 'native' / f'{v}_{gname}')
        for t in tr:
            k = (t['lo'][0], t['lo'][1], t['up'][0], t['up'][1])
            if (t['lo'] not in levels) or (t['up'] not in levels):
                continue
            de_exp = levels[t['up']] - levels[t['lo']]
            res.setdefault(k, {})[f'gf_{gname}'] = corrected(t['gf'], t['dE_eV'] * HC_EV_CM, de_exp, gname)
    return {k: v for k, v in res.items() if 'gf_B' in v and 'gf_C' in v}


def rate(gf, k, levels):
    lam = 1e8 / (levels[(k[2], k[3])] - levels[(k[0], k[1])])
    return A_COEF * gf / ((k[3] + 1) * lam ** 2)


def dex(a, b):
    return round(math.log10(a / b), 4)


def robust(d):
    if len(d) < 2:
        return None
    m = st.median(d)
    return round(1.4826 * st.median([abs(x - m) for x in d]), 4)


def score(v, levels, ff9, ff16):
    m = model(v, levels)
    # 1. experiment: D-line rates
    exp = {f'{k[2]}_{k[3]}/2': dict(A_exp=a, frac_sigma=s, A_B=rate(m[k]['gf_B'], k, levels), A_C=rate(m[k]['gf_C'], k, levels),
                                    dlogA_B=dex(rate(m[k]['gf_B'], k, levels), a), dlogA_C=dex(rate(m[k]['gf_C'], k, levels), a),
                                    ff9_dlogA=dex(rate(ff9[k]['gf_L'] * (levels[(k[2], k[3])] - levels[(k[0], k[1])]) / ff9[k]['sigma_th'], k, levels), a)
                                    if k in ff9 else None)
           for k, (a, s) in EXP_D.items() if k in m}
    # 2/3. line by line vs FF, grouped by our class
    rows = []
    for k, g in m.items():
        dT = G.dT(g['gf_B'], g['gf_C'])
        a_b = rate(g['gf_B'], k, levels)
        cls = G.accuracy_class(dT, 1.0)  # single-CSF valence levels: CF not defined -> 1 (see doc)
        row = dict(lower=f'{k[0]}_{k[1]}/2', upper=f'{k[2]}_{k[3]}/2', wl_vac_A=round(1e8 / (levels[(k[2], k[3])] - levels[(k[0], k[1])]), 3),
                   loggf_B=round(math.log10(g['gf_B']), 4), loggf_C=round(math.log10(g['gf_C']), 4), dT=round(dT, 4), A_B=a_b, cls=cls)
        for tag, ff in (('ff9', ff9), ('ff16', ff16)):
            if k in ff:
                e = levels[(k[2], k[3])] - levels[(k[0], k[1])]
                ffgf = ff[k]['gf_L'] * e / ff[k]['sigma_th']  # FF length form on the experimental energy
                row[f'{tag}_loggf'] = round(math.log10(ffgf), 4)
                row[f'd_{tag}'] = round(row['loggf_B'] - row[f'{tag}_loggf'], 4)
        rows.append(row)
    # dT~ by A-magnitude group (Li 2023)
    tl = G.dT_tilde([dict(A_B=r['A_B'], A_C=rate(10 ** r['loggf_C'], (r['lower'].split('_')[0], int(r['lower'].split('_')[1][:-2]),
                                                                       r['upper'].split('_')[0], int(r['upper'].split('_')[1][:-2])), levels))
                     for r in rows])
    for r, t in zip(rows, tl):
        r['dT_tilde'] = round(t['dT_tilde'], 4)
        r['cls'] = G.accuracy_class(t['dT_tilde'], 1.0)
    by = defaultdict(list)
    for r in rows:
        if 'd_ff9' in r:
            by[r['cls']].append(r['d_ff9'])
    vs_ff = {c: dict(n=len(d), median=round(st.median(d), 4), robust_sigma=robust(d), rms=round(math.sqrt(sum(x * x for x in d) / len(d)), 4))
             for c, d in sorted(by.items())}
    d_all = [r['d_ff9'] for r in rows if 'd_ff9' in r]
    return m, rows, dict(experiment=exp, vs_ff9_by_class=vs_ff,
                         vs_ff9_all=dict(n=len(d_all), median=round(st.median(d_all), 4), robust_sigma=robust(d_all)) if d_all else {},
                         n_lines=len(rows), class_counts={c: sum(r['cls'] == c for r in rows) for c in 'ABCDE'})


def ff_spread(levels, ff9, ff16):
    """Independent-theory spread: FF B-spline CI [9] vs FF MCHF-CV [16], length forms on experimental energies."""
    d = []
    for k in set(ff9) & set(ff16):
        if (k[0], k[1]) in levels and (k[2], k[3]) in levels:
            d.append(math.log10(ff9[k]['gf_L'] / ff9[k]['sigma_th']) - math.log10(ff16[k]['gf_L'] / ff16[k]['sigma_th']))
    gauge9 = [abs(math.log10(v['gf_L'] / v['gf_V'])) for v in ff9.values() if v['gf_V'] > 0]
    return dict(n_common=len(d), median=round(st.median(d), 4), robust_sigma=robust(d),
                ff9_gauge_median_abs_dex=round(st.median(gauge9), 4))


def solar_lines(levels):
    """Every Na I row in the solar line list, collapsed over HFS components onto fine-structure lines."""
    up_by = {}
    for k, e in levels.items():
        up_by[k] = e
    agg = {}
    for r in csv.DictReader(open(ROOT / 'data/linelists/linelist_solar.csv')):
        if r['element'] != 'Na' or r['ion'] != 'I':
            continue
        air = float(r['wavelength_air_A'])
        vac = air * (1 + 8.336624212083e-5 + 2.408926869968e-2 / (1.301065924522e2 - (1e4 / air) ** 2)
                     + 1.599740894897e-4 / (3.892568793293e1 - (1e4 / air) ** 2))
        elo = float(r['excitation_potential_eV']) * HC_EV_CM
        lo = [k for k, e in levels.items() if abs(e - elo) < 2.0]
        if not lo:
            continue
        eup = elo + 1e8 / vac
        hit = [(lk, uk) for lk in lo for uk, e in levels.items() if abs(e - eup) < 1.5 and L.index(uk[0][-1]) - L.index(lk[0][-1]) in (-1, 1)]
        if len(hit) != 1:
            continue
        key = (hit[0][0][0], hit[0][0][1], hit[0][1][0], hit[0][1][1])
        a = agg.setdefault(key, dict(air_A=[], gf=0.0, depth=0.0, sources=set()))
        a['air_A'].append(air)
        a['gf'] += 10 ** float(r['log_gf'])
        a['depth'] = max(a['depth'], float(r['central_depth'] or 0))
        a['sources'].add(r['loggf_source'])
    return agg


def main_score():
    levels = nist_levels()
    ff9, ff16 = parse_ff('Na_11.9.U.lsj.db.gz'), parse_ff('Na_11.16.mchf_cv.lsj.db.gz')
    summary = dict(schema='rya1354.na-model.v61', status=STATUS, n_nist_levels=len(levels), ff_spread=ff_spread(levels, ff9, ff16), variants={})
    best = None
    for v in VARIANTS:
        if not (OUT / 'native' / f'{v}_B' / 't.tr.gz').exists() and not (OUT / 'native' / f'{v}_B' / 't.tr').exists():
            continue
        m, rows, s = score(v, levels, ff9, ff16)
        summary['variants'][v] = s
        with open(OUT / f'{v}_lines.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=sorted({k for r in rows for k in r}), lineterminator='\n')
            w.writeheader()
            w.writerows(sorted(rows, key=lambda r: r['wl_vac_A']))
        best = (v, rows) if best is None or abs(s['vs_ff9_all'].get('median', 9)) < abs(summary['variants'][best[0]]['vs_ff9_all'].get('median', 9)) else best
    # fill: every solar-line-list Na I fine-structure line
    if best:
        v, rows = best
        byk = {(r['lower'], r['upper']): r for r in rows}
        fill = []
        for k, a in sorted(solar_lines(levels).items(), key=lambda kv: min(kv[1]['air_A'])):
            r = byk.get((f'{k[0]}_{k[1]}/2', f'{k[2]}_{k[3]}/2'))
            fill.append(dict(lower=f'{k[0]}_{k[1]}/2', upper=f'{k[2]}_{k[3]}/2', air_A=round(st.median(a['air_A']), 3), n_hfs_rows=len(a['air_A']),
                             max_depth=a['depth'], linelist_loggf=round(math.log10(a['gf']), 3), linelist_sources='|'.join(sorted(a['sources'])),
                             model=v, model_loggf=r['loggf_B'] if r else None, model_class=r['cls'] if r else None,
                             model_dT_tilde=r['dT_tilde'] if r else None, ff9_loggf=r.get('ff9_loggf') if r else None,
                             d_model_minus_linelist=round(r['loggf_B'] - math.log10(a['gf']), 3) if r else None))
        with open(OUT / 'solar_fill.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(fill[0]), lineterminator='\n')
            w.writeheader()
            w.writerows(fill)
        summary['fill'] = dict(model=v, n_lines=len(fill), n_filled=sum(x['model_loggf'] is not None for x in fill),
                               n_measurable=sum(x['max_depth'] >= 0.01 for x in fill))
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['run', 'score'])
    ap.add_argument('--out')
    a = ap.parse_args()
    if a.cmd == 'run':
        run(a.out)
    else:
        print(json.dumps(main_score(), indent=1))

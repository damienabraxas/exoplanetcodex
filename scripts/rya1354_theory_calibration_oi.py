"""Calibrate THEORETICAL O I gf uncertainty against every experimental anchor that exists, and rank
O I line candidates for a Theoretical Codex grade (RYA-1354 / RYA-1352). Mac-side, no solver.

Inputs (all in the repo):
  data/reference/li2023_oi/raw/tablea{1,2}.tsv   Li, Joensson, Amarsi, Li & Grumer 2023 (A&A 674, A54),
                                                 CDS J/A+A/674/A54 via the VizieR CfA mirror
  data/reference/li2023_oi/table2_experimental_lifetimes.csv, table3_experimental_gf.csv
                                                 experimental values AS QUOTED in Li 2023 Tables 2-3
                                                 (secondary transcription; primaries not re-verified)
  data/reference/nist/oi_lines_asd_rya1354.csv   NIST ASD O I lines (Wiese, Fuhr & Deters 1996 evaluation;
                                                 its O I A-values are themselves CALCULATED, mostly
                                                 Hibbert et al. 1991 = TP ref T5803)
  data/linelists/linelist_solar.csv              our solar line list (VALD3), O I rows

Three tests:
  1. Direct gf: Li vs measured gf (12 UV lines, 98.9-130.6 nm; Table 3).
  2. Lifetimes: Li vs measured lifetimes (13 levels; Table 2). A lifetime tests the SUM of A out of a
     level, i.e. the dominant decay lines; those lines (incl. VIS/IR 777, 844.6, 615.8, 1130, 1316 nm)
     are listed per level with branching fractions.
  3. Independent theory: Li (MCDHF/RCI) vs NIST (Hibbert-based CIV3/LS evaluation) on every common
     line, by Li accuracy class, A-magnitude and band. Spread / sqrt(2) = per-method scatter under the
     equal-quality assumption (labelled as such).

Status: CANDIDATE, HOLD. Nothing here writes a grade or sigma into a product.
"""
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import statistics as st

ROOT = Path(__file__).resolve().parents[1]
LI = ROOT / 'data/reference/li2023_oi'
OUT = ROOT / 'data/audit/rya1354/theory_calibration_oi'
CM_PER_EV = 8065.54429
CLASS_FRAC = {'A': 0.03, 'B': 0.10, 'C': 0.25, 'D': 0.50, 'E': 1.00}   # Li Table 5 upper limits (gf, A)


def dex(frac):
    """Fractional gf uncertainty -> dex (symmetric, log10(1+f))."""
    return math.log10(1 + frac)


def robust_sigma(x):
    m = st.median(x)
    return 1.4826 * st.median(abs(v - m) for v in x)


def rms(x):
    return math.sqrt(sum(v * v for v in x) / len(x))


def read_vizier(path):
    lines = [l for l in path.read_text().splitlines() if l and not l.startswith('#')]
    head = lines[0].split('\t')
    rows = []
    for l in lines[3:]:
        v = [s.strip() for s in l.split('\t')]
        if len(v) == len(head):
            rows.append(dict(zip(head, v)))
    return rows


def conf(c):
    return c.replace('1s(2).', '').strip()


def levels():
    out = {}
    for r in read_vizier(LI / 'raw/tablea1.tsv'):
        key = (conf(r['State']), r['Term'], int(r['J']), r['Parity'])
        out[key] = dict(no=int(r['No']), E_nist=float(r['ENIST']),
                        tau_B=float(r['taul']) if r['taul'] else None, tau_C=float(r['tauv']) if r['tauv'] else None)
    if len(out) != 81:
        raise ValueError(f'Li Table A1 should hold 81 states, got {len(out)}')
    return out


def transitions(lev):
    out = []
    for r in read_vizier(LI / 'raw/tablea2.tsv'):
        u = (conf(r['Upper']), r['Termu'], int(r['Ju']), r['Pu'])
        l = (conf(r['Lower']), r['Terml'], int(r['Jl']), r['Pl'])
        if u not in lev or l not in lev:
            raise ValueError(f'Li transition references an unknown level: {u} / {l}')
        out.append(dict(upper=u, lower=l, E_up=lev[u]['E_nist'], E_lo=lev[l]['E_nist'],
                        lambda_vac_nm=float(r['lambda']), loggf_B=float(r['loggf-B']), loggf_C=float(r['loggf-C']),
                        A_B=float(r['A-B']), A_C=float(r['A-C']), dT=float(r['dT']), CF_B=float(r['CF-B']),
                        CF_C=float(r['CF-C']), acc=r['Acc']))
    return out


def nist_lines():
    clean = lambda s: s.replace('=', '').strip('"').strip()
    rows = []
    for r in csv.DictReader(open(ROOT / 'data/reference/nist/oi_lines_asd_rya1354.csv')):
        r = {k: clean(v) for k, v in r.items() if k}
        if not r.get('fik') or not r.get('Ei(cm-1)') or not r.get('Ek(cm-1)'):
            continue
        try:
            Ei, Ek = float(r['Ei(cm-1)'].strip('[]()')), float(r['Ek(cm-1)'].strip('[]()'))
            Ji, Jk = r['J_i'], r['J_k']
            gi = int(r['g_i'])
        except ValueError:
            continue
        rows.append(dict(Ei=Ei, Ek=Ek, J_i=Ji, J_k=Jk, gf=gi * float(r['fik']), acc=r['Acc'], src=r['tp_ref'],
                         lambda_vac_nm=r['ritz_wl_vac(nm)'] or r['obs_wl_vac(nm)'],
                         lower=f"{r['conf_i']} {r['term_i']}{r['J_i']}", upper=f"{r['conf_k']} {r['term_k']}{r['J_k']}"))
    return rows


def match_nist(tr, nist, tol=0.06):
    """Same lower and upper level: NIST energies within tol cm-1 of Li's NIST-energy column AND same J."""
    idx = defaultdict(list)
    for n in nist:
        idx[(round(n['Ei']), round(n['Ek']))].append(n)
    for t in tr:
        hits = []
        for di in (-1, 0, 1):
            for dk in (-1, 0, 1):
                for n in idx.get((round(t['E_lo']) + di, round(t['E_up']) + dk), []):
                    if (abs(n['Ei'] - t['E_lo']) <= tol and abs(n['Ek'] - t['E_up']) <= tol
                            and n['J_i'] == str(t['lower'][2]) and n['J_k'] == str(t['upper'][2])):
                        hits.append(n)
        t['nist'] = hits[0] if len(hits) == 1 else None
        t['nist_ambiguous'] = len(hits) > 1


def band(nm_vac):
    return 'UV' if nm_vac < 300 else 'VIS' if nm_vac < 700 else 'IR'


def group_a(a):
    return 'A<1e3' if a < 1e3 else '1e3-1e5' if a < 1e5 else '1e5-1e6' if a < 1e6 else 'A>=1e6'


def fit_sigma_th(d, sig):
    """Extra (theory) scatter s >= 0 such that sum d^2/(sig^2+s^2) = N (zero-bias test). 0 if already consistent."""
    chi = lambda s: sum(x * x / (e * e + s * s) for x, e in zip(d, sig))
    if chi(0) <= len(d):
        return 0.0
    lo, hi = 0.0, 5.0
    for _ in range(80):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if chi(mid) > len(d) else (lo, mid)
    return round((lo + hi) / 2, 4)


def nf_mixed(t):
    """Upper 2p3(4S)nf: 5F/3F nearly degenerate and strongly mixed (near jK coupling), so per-component LS
    labels differ between calculations; compare and grade these as blended multiplets."""
    return t['upper'][0].endswith('f') or t['lower'][0].endswith('f')


def summarise(d):
    return dict(n=len(d), median=round(st.median(d), 4), robust_sigma=round(robust_sigma(d), 4),
                rms=round(rms(d), 4), per_method_robust_sigma=round(robust_sigma(d) / math.sqrt(2), 4)) if len(d) >= 3 else dict(n=len(d))


def test_direct_gf(tr, lev):
    rows = list(csv.DictReader(open(LI / 'table3_experimental_gf.csv')))
    by = defaultdict(list)
    for r in rows:
        by[(r['upper_conf'], r['upper_term'], r['upper_J'], r['lower_term'], r['lower_J'], r['lambda_vac_nm'])].append(r)
    out = []
    for k, rs in by.items():
        w = [1 / float(r['exp_sigma']) ** 2 for r in rs]
        mean = sum(wi * float(r['exp_gf']) for wi, r in zip(w, rs)) / sum(w)
        sig = 1 / math.sqrt(sum(w))
        chi = sum(((float(r['exp_gf']) - mean) / float(r['exp_sigma'])) ** 2 for r in rs)
        birge = math.sqrt(chi / (len(rs) - 1)) if len(rs) > 1 else 1.0
        sig *= max(1.0, birge)                  # inflate when the experiments disagree with each other
        li_B = float(rs[0]['li_gf_B'])
        # accuracy class from Li Table A2 (the J=1+2 blend takes the worse class of the two)
        Js = [int(j) for j in k[2].split('+')]
        accs = [t['acc'] for t in tr if t['upper'][0] == k[0] and t['upper'][1] == k[1] and t['upper'][2] in Js
                and t['lower'][0] == '2s(2).2p(4)' and t['lower'][1] == k[3] and t['lower'][2] == int(k[4])]
        acc = max(accs) if accs else None
        d = math.log10(li_B / mean)
        sig_exp_dex = sig / mean / math.log(10)
        out.append(dict(line=f'{k[0]} {k[1]}{k[2]} - 2p4 {k[3]}{k[4]}', lambda_vac_nm=float(k[5]), band='UV',
                        li_gf_B=li_B, exp_gf=round(mean, 5), exp_sigma=round(sig, 5), n_exp=len(rs), birge=round(birge, 2),
                        li_class=acc, dlog_gf_li_minus_exp=round(d, 4), exp_sigma_dex=round(sig_exp_dex, 4),
                        z_exp_only=round(d / sig_exp_dex, 2)))
    return out


def test_lifetimes(tr, lev):
    rows = list(csv.DictReader(open(LI / 'table2_experimental_lifetimes.csv')))
    by = defaultdict(list)
    for r in rows:
        by[(r['conf'], r['term'])].append(r)
    out = []
    for (c, term), rs in by.items():
        keys = [k for k in lev if k[0] == c and k[1] == term]
        if not keys:
            raise ValueError(f'No Li level {c} {term}')
        unit = 1e-6 if rs[0]['unit'] == 'us' else 1e-9
        tau_li = st.mean(lev[k]['tau_B'] for k in keys) / unit
        if abs(tau_li / float(rs[0]['li_tau_B']) - 1) > 0.01:     # transcription check vs Table A1
            raise ValueError(f'Table 2 vs A1 lifetime mismatch for {c} {term}: {tau_li} vs {rs[0]["li_tau_B"]}')
        # decay channels out of these upper levels (Li Babushkin A), with branching fractions
        ch = defaultdict(float)
        for t in tr:
            if t['upper'][0] == c and t['upper'][1] == term:
                ch[(t['lower'][0], t['lower'][1], round(t['lambda_vac_nm'], 1))] += t['A_B']
        tot = sum(ch.values())
        main = sorted(((a / tot, k) for k, a in ch.items()), reverse=True)
        meas = []
        for r in rs:
            e, s = float(r['exp_tau']), float(r['exp_sigma'])
            meas.append(dict(ref=r['exp_ref'], tau=e, sigma=s, dlogA_li_minus_exp=round(math.log10(e / tau_li), 4),
                             z=round((e - tau_li) / s, 2)))
        w = [1 / m['sigma'] ** 2 for m in meas]
        mean = sum(wi * m['tau'] for wi, m in zip(w, meas)) / sum(w)
        sig = 1 / math.sqrt(sum(w))
        chi = sum(((m['tau'] - mean) / m['sigma']) ** 2 for m in meas)
        birge = math.sqrt(chi / (len(meas) - 1)) if len(meas) > 1 else 1.0
        out.append(dict(level=f'{c} {term}', unit=rs[0]['unit'], li_tau_B=round(tau_li, 4), exp_weighted_mean=round(mean, 4),
                        exp_sigma=round(sig * max(1, birge), 4), birge=round(birge, 2), n_exp=len(meas),
                        dlogA_li_minus_exp=round(math.log10(mean / tau_li), 4),
                        z=round((mean - tau_li) / (sig * max(1, birge)), 2), measurements=meas,
                        decay_channels=[dict(lower=f'{k[0]} {k[1]}', lambda_vac_nm=k[2], band=band(k[2]), branching=round(b, 3))
                                        for b, k in main if b >= 0.05]))
    return out


def test_theory(tr, nist):
    m = [t for t in tr if t['nist']]
    rows = []
    for t in m:
        n = t['nist']
        rows.append(dict(lambda_vac_nm=t['lambda_vac_nm'], band=band(t['lambda_vac_nm']), li_class=t['acc'], nf_label_mixed=nf_mixed(t), dT=t['dT'],
                         CF_B=t['CF_B'], A_B=t['A_B'], A_group=group_a(t['A_B']), li_loggf_B=t['loggf_B'], li_loggf_C=t['loggf_C'],
                         nist_loggf=round(math.log10(n['gf']), 4), nist_acc=n['acc'], nist_src=n['src'],
                         d_li_minus_nist=round(t['loggf_B'] - math.log10(n['gf']), 4),
                         lower=f"{t['lower'][0]} {t['lower'][1]}{t['lower'][2]}", upper=f"{t['upper'][0]} {t['upper'][1]}{t['upper'][2]}"))
    res = dict(n_li=len(tr), n_matched=len(m), n_ambiguous=sum(t['nist_ambiguous'] for t in tr),
               n_nf_label_mixed=sum(r['nf_label_mixed'] for r in rows))
    clean = [r for r in rows if not r['nf_label_mixed']]
    for key in ('li_class', 'A_group', 'band', 'nist_acc'):
        res[f'by_{key}'] = {g: summarise([r['d_li_minus_nist'] for r in clean if r[key] == g])
                            for g in sorted({r[key] for r in clean})}
    res['by_li_class_strong_A>=1e6'] = {g: summarise([r['d_li_minus_nist'] for r in clean if r['li_class'] == g and r['A_B'] >= 1e6])
                                         for g in 'ABCDE'}
    res['by_li_class_incl_nf_components'] = {g: summarise([r['d_li_minus_nist'] for r in rows if r['li_class'] == g]) for g in 'ABCDE'}
    res['nf_multiplet_sums'] = {}
    for lo in ('3d 5D', '3d 3D'):
        # ALL components on both sides (label swaps between 5F and 3F move gf between matched pairs)
        sel = [t for t in tr if f"{t['lower'][0][-2:]} {t['lower'][1]}" == lo and t['upper'][0].endswith('4f')]
        e_lo = {round(t['E_lo'], 2) for t in sel}
        nsel = [n for n in nist if any(abs(n['Ei'] - e) <= 0.06 for e in e_lo) and '4f' in n['upper']]
        a = sum(10 ** t['loggf_B'] for t in sel)
        b = sum(n['gf'] for n in nsel)
        res['nf_multiplet_sums'][f'{lo} - 4f (1.80/1.82 um)'] = dict(n_li=len(sel), n_nist=len(nsel), li_log_sum_gf=round(math.log10(a), 4),
                                                                   nist_log_sum_gf=round(math.log10(b), 4), d=round(math.log10(a / b), 4))
    res['gauge_B_minus_C_by_band'] = {b: summarise([r['li_loggf_B'] - r['li_loggf_C'] for r in rows if r['band'] == b]) for b in ('UV', 'VIS', 'IR')}
    res['claimed_class_sigma_dex'] = {c: round(dex(f), 4) for c, f in CLASS_FRAC.items()}
    return rows, res


def candidates(tr, lifetimes, direct):
    """Our solar line list O I rows matched to Li by lower-level energy and vacuum wavenumber."""
    ll = [r for r in csv.DictReader(open(ROOT / 'data/linelists/linelist_solar.csv')) if r['element'] == 'O' and r['ion'] in ('1', 'I')]
    anchored_levels = {l['level']: l for l in lifetimes}
    out = []
    for r in ll:
        lam_air = float(r['wavelength_air_A']) / 10
        sig2 = (1e3 / lam_air) ** 2
        n_air = 1 + 0.0000834254 + 0.02406147 / (130 - sig2) + 0.00015998 / (38.9 - sig2) if lam_air > 200 else 1.0
        lam_vac = lam_air * n_air
        E_lo = float(r['excitation_potential_eV']) * CM_PER_EV
        hit = [t for t in tr if abs(t['E_lo'] - E_lo) < 3 and abs(t['lambda_vac_nm'] - lam_vac) < 0.02]
        t = min(hit, key=lambda t: abs(t['lambda_vac_nm'] - lam_vac)) if hit else None
        up = f"{t['upper'][0]} {t['upper'][1]}" if t else None
        anchor = []
        if t and any(abs(d['lambda_vac_nm'] - t['lambda_vac_nm']) < 0.01 for d in direct):
            anchor.append('direct_gf')
        if up in anchored_levels and any(abs(c['lambda_vac_nm'] - round(t['lambda_vac_nm'], 1)) < 0.06 and c['branching'] >= 0.10
                                          for c in anchored_levels[up]['decay_channels']):
            anchor.append('upper_lifetime')
        blend = None
        if t and nf_mixed(t):
            # The nf components of one lower TERM span <=0.12 nm (5F/3F upper levels 0.094 cm-1 apart = 0.03 nm at
            # 1.82 um), i.e. one blended solar feature (line width ~0.03-0.04 nm there). LS-coupled VALD/NIST and
            # intermediate-coupled Li distribute gf differently among components (up to 0.9 dex per feature) but
            # agree on the TERM total, so the graded quantity is the term-to-configuration multiplet sum.
            li_c = [x for x in tr if x['lower'][:2] == t['lower'][:2] and x['upper'][0] == t['upper'][0]]
            va_c = [x for x in ll if abs(float(x['excitation_potential_eV']) - float(r['excitation_potential_eV'])) < 2e-3
                    and abs(float(x['wavelength_air_A']) / 10 - lam_air) < 0.2]
            blend = dict(n_li=len(li_c), n_vald=len(va_c), li=round(math.log10(sum(10 ** x['loggf_B'] for x in li_c)), 3),
                         vald=round(math.log10(sum(10 ** float(x['log_gf']) for x in va_c)), 3))
        depth = float(r['central_depth'] or 0)
        cls = t['acc'] if t else None
        cf = t['CF_B'] if t else None
        if t is None:
            tier = 'NOT_IN_LI_SET'
        elif anchor:
            tier = 'T1_EXPERIMENTAL_ANCHOR'
        elif (cls in ('A', 'B') or (nf_mixed(t) and blend)) and cf >= 0.1 and depth >= 0.01:   # nf: graded as blend
            tier = 'T2_THEORETICAL_CODEX_CANDIDATE'
        elif depth >= 0.01:
            tier = 'T3_MEASURABLE_LOW_THEORY_CLASS'
        else:
            tier = 'T4_TOO_WEAK_IN_SUN'
        out.append(dict(lambda_air_nm=round(lam_air, 4), lambda_vac_nm=round(lam_vac, 4), band=band(lam_vac),
                        chi_eV=float(r['excitation_potential_eV']), vald_loggf=float(r['log_gf']), solar_depth=depth,
                        li_loggf_B=t['loggf_B'] if t else None, li_loggf_C=t['loggf_C'] if t else None,
                        li_class=cls, li_dT=t['dT'] if t else None, li_CF_B=cf,
                        nist_loggf=round(math.log10(t['nist']['gf']), 3) if t and t['nist'] else None,
                        nist_acc=t['nist']['acc'] if t and t['nist'] else None,
                        d_vald_minus_li=round(float(r['log_gf']) - t['loggf_B'], 3) if t else None,
                        transition=f"{t['lower'][0]} {t['lower'][1]}{t['lower'][2]} - {t['upper'][0]} {t['upper'][1]}{t['upper'][2]}" if t else None,
                        anchor='+'.join(anchor) or None, tier=tier, grade_as_multiplet=bool(t and nf_mixed(t)),
                        blend_li_loggf=blend['li'] if blend else None, blend_vald_loggf=blend['vald'] if blend else None,
                        blend_n_li=blend['n_li'] if blend else None, blend_n_vald=blend['n_vald'] if blend else None,
                        ground_observable=lam_air >= 300))
    return sorted(out, key=lambda c: (c['tier'], -c['solar_depth']))


def own_model_score(life):
    """Score OUR sealed FAC v41 models (RYA-1354 pilot) with the same yardsticks: Li's dT class from the two
    gauges, and the experimental 3s 3S lifetime (130 nm upper level), the one quantity they compute end to end."""
    v41 = json.loads((ROOT / 'data/audit/rya1354/source_curation_v2/fac_radial_rates_v41.json').read_text())
    exp = next(l for l in life if l['level'] == '2s(2).2p(3).3s 3S')
    out = {}
    for case, g in v41['paired_n3_radial_comparison'].items():
        tb, tc = g['B']['radial_lifetime_ns'], g['C']['radial_lifetime_ns']
        dT = 1 - min(tb, tc) / max(tb, tc)                # A ~ 1/tau
        p = 100 * dT
        cls = 'ABCDE'[0 if p <= 3 else 1 if p <= 10 else 2 if p <= 25 else 3 if p <= 50 else 4]
        out[case] = dict(tau_B_ns=round(tb, 3), tau_C_ns=round(tc, 3), dT=round(dT, 3), li_scale_class=cls,
                         dlogA_B_minus_exp=round(math.log10(exp['exp_weighted_mean'] / tb), 3))
    return dict(source='fac_radial_rates_v41.json (sealed exploratory FAC models)', level='3s 3S (130.2/130.5/130.6 nm)',
                exp_tau_ns=exp['exp_weighted_mean'], li_tau_ns=exp['li_tau_B'], cases=out)


def main():
    lev = levels()
    tr = transitions(lev)
    nist = nist_lines()
    match_nist(tr, nist)
    direct = test_direct_gf(tr, lev)
    life = test_lifetimes(tr, lev)
    theory_rows, theory = test_theory(tr, nist)
    cand = candidates(tr, life, direct)
    OUT.mkdir(parents=True, exist_ok=True)
    for name, rows in (('direct_gf', direct), ('li_vs_nist', theory_rows), ('candidates', cand)):
        with open(OUT / f'{name}.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
            w.writeheader()
            w.writerows(rows)
    d = [r['dlog_gf_li_minus_exp'] for r in direct]
    dl = [l['dlogA_li_minus_exp'] for l in life]
    summary = dict(
        schema='rya1354.theory-calibration-oi.v1', status='CANDIDATE_HOLD_NOT_ADOPTED',
        sources=dict(li2023='10.1051/0004-6361/202245645; CDS J/A+A/674/A54', nist='NIST ASD O I (Wiese, Fuhr & Deters 1996)',
                     experimental='as quoted in Li 2023 Tables 2-3; primaries not re-verified'),
        direct_gf=dict(n=len(d), mean=round(st.mean(d), 4), rms=round(rms(d), 4), robust_sigma=round(robust_sigma(d), 4),
                       chi2_exp_only=round(sum(r['z_exp_only'] ** 2 for r in direct), 2), band='UV only (98.9-130.6 nm)',
                       extra_theory_sigma_dex=fit_sigma_th(d, [r['exp_sigma_dex'] for r in direct]),
                       by_li_class={c: dict(n=len(x), rms=round(rms(x), 4)) for c in 'ABC'
                                    for x in [[r['dlog_gf_li_minus_exp'] for r in direct if r['li_class'] == c]] if x}),
        lifetimes=dict(n_levels=len(life), median_dlogA=round(st.median(dl), 4), rms_dlogA=round(rms(dl), 4),
                       extra_theory_sigma_dex=fit_sigma_th(dl, [l['exp_sigma'] / l['exp_weighted_mean'] / math.log(10) for l in life]),
                       extra_theory_sigma_dex_birge_ok=fit_sigma_th([l['dlogA_li_minus_exp'] for l in life if l['birge'] <= 2],
                                                                    [l['exp_sigma'] / l['exp_weighted_mean'] / math.log(10) for l in life if l['birge'] <= 2]),
                       per_level=[{k: l[k] for k in ('level', 'li_tau_B', 'exp_weighted_mean', 'exp_sigma', 'birge', 'dlogA_li_minus_exp', 'z', 'decay_channels')} for l in life]),
        theory_vs_theory=theory,
        candidates={t: sum(c['tier'] == t for c in cand) for t in sorted({c['tier'] for c in cand})})
    summary['own_model_v41'] = own_model_score(life)
    floor = summary['lifetimes']['extra_theory_sigma_dex_birge_ok']
    spread = theory['by_li_class']
    summary['proposed_sigma_dex_by_class'] = dict(
        rule='sigma = sqrt(floor^2 + spread^2): floor = extra scatter the experimental LIFETIMES need beyond their own '
             'errors (levels whose experiments agree, Birge<=2); spread = per-method Li-vs-NIST robust scatter for the '
             'class (strong and weak lines, nf components excluded). CANDIDATE: needs Ryan + RYA-1352 ruling.',
        floor_dex=floor,
        classes={c: dict(li_claimed_max_dex=round(dex(CLASS_FRAC[c]), 3), method_spread_dex=spread[c]['per_method_robust_sigma'],
                         proposed_dex=round(math.hypot(floor, spread[c]['per_method_robust_sigma']), 3),
                         proposed_percent=round(100 * (10 ** math.hypot(floor, spread[c]['per_method_robust_sigma']) - 1), 1))
                 for c in 'ABCDE'})
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


if __name__ == '__main__':
    s = main()
    print(json.dumps({k: v for k, v in s.items() if k not in ('lifetimes', 'theory_vs_theory')}, indent=1))
    for l in s['lifetimes']['per_level']:
        print(f"{l['level']:24s} Li {l['li_tau_B']:9.3f} exp {l['exp_weighted_mean']:9.3f}±{l['exp_sigma']:<7} birge {l['birge']:5} dlogA {l['dlogA_li_minus_exp']:+.3f} z {l['z']:+.1f}  "
              + ', '.join(f"{c['lambda_vac_nm']}nm({c['branching']})" for c in l['decay_channels']))
    t = s['theory_vs_theory']
    print({k: t[k] for k in ('n_li', 'n_matched', 'n_ambiguous')})
    print({k: s['lifetimes'][k] for k in ('median_dlogA', 'rms_dlogA', 'extra_theory_sigma_dex', 'extra_theory_sigma_dex_birge_ok')})
    for k in ('by_li_class', 'by_li_class_strong_A>=1e6', 'by_li_class_incl_nf_components', 'nf_multiplet_sums', 'by_A_group', 'by_band', 'by_nist_acc', 'gauge_B_minus_C_by_band'):
        print(k, json.dumps(t[k]))

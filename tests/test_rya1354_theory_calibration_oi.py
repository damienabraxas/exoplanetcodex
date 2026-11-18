"""RYA-1354 O I theory calibration: sources pinned, transcriptions validated against the machine-readable
CDS tables, and the committed summary must equal a fresh recomputation (no artifact drift)."""
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import rya1354_theory_calibration_oi as T  # noqa: E402

RAW_SHA = {'data/reference/li2023_oi/raw/tablea1.tsv': '127c0ff567f61c2bef450d4c37bb463c97f3aa17f1cfeee1ba67890554e7fe9b',
           'data/reference/li2023_oi/raw/tablea2.tsv': '943672d7e96a688da56e2383acdc6ecab0cb7fbad49d06cad533f6d8be16326f',
           'data/reference/nist/oi_lines_asd_rya1354.csv': 'f4e652ebcc16c107aa486cdb60cef24910d9801751c15af3bee917601717d9ab'}


def test_sources_pinned():
    for p, h in RAW_SHA.items():
        assert hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == h, p


def test_li_tables_complete_and_every_transition_resolves_to_a_level():
    lev = T.levels()
    tr = T.transitions(lev)
    assert len(lev) == 81 and len(tr) == 989     # Li 2023: 81 states, 989 E1 transitions (Sect. 3.2)


def test_table3_transcription_matches_cds_table_a2():
    """Li's printed Table 3 gf (Babushkin) must equal 10**loggf-B of the same transition in Table A2."""
    tr = T.transitions(T.levels())
    checked = 0
    for r in csv.DictReader(open(ROOT / 'data/reference/li2023_oi/table3_experimental_gf.csv')):
        Js = [int(j) for j in r['upper_J'].split('+')]
        sel = [t for t in tr if t['upper'][0] == r['upper_conf'] and t['upper'][1] == r['upper_term'] and t['upper'][2] in Js
               and t['lower'][0] == r['lower_conf'] and t['lower'][1] == r['lower_term'] and t['lower'][2] == int(r['lower_J'])]
        assert len(sel) == len(Js), r
        gf = sum(10 ** t['loggf_B'] for t in sel)
        assert abs(gf / float(r['li_gf_B']) - 1) < 0.03, (r['lambda_vac_nm'], gf, r['li_gf_B'])
        checked += 1
    assert checked == 27


def test_nist_match_is_by_level_identity_and_unambiguous():
    lev = T.levels()
    tr = T.transitions(lev)
    T.match_nist(tr, T.nist_lines())
    m = [t for t in tr if t['nist']]
    assert len(m) == 347 and not any(t['nist_ambiguous'] for t in tr)
    for t in m:
        assert abs(t['nist']['Ei'] - t['E_lo']) <= 0.06 and abs(t['nist']['Ek'] - t['E_up']) <= 0.06


def test_fit_sigma_th_recovers_known_scatter():
    assert T.fit_sigma_th([0.01, -0.01], [0.05, 0.05]) == 0.0
    s = T.fit_sigma_th([0.3, -0.3, 0.3, -0.3], [0.1, 0.1, 0.1, 0.1])
    assert abs(math.hypot(0.1, s) - 0.3) < 1e-3


def test_committed_outputs_equal_a_fresh_run(tmp_path, monkeypatch):
    committed = json.loads((T.OUT / 'summary.json').read_text())
    monkeypatch.setattr(T, 'OUT', tmp_path)
    fresh = json.loads(json.dumps(T.main()))
    assert fresh == committed
    for name in ('direct_gf', 'li_vs_nist', 'candidates'):
        assert (tmp_path / f'{name}.csv').read_text() == (ROOT / 'data/audit/rya1354/theory_calibration_oi' / f'{name}.csv').read_text()


def test_key_findings():
    s = json.loads((T.OUT / 'summary.json').read_text())
    life = {l['level']: l for l in s['lifetimes']['per_level']}
    # 844.6 nm upper level: the Bischel two-photon LIF lifetime says theory A is ~0.08 dex too strong
    assert 0.07 < life['2s(2).2p(3).3p 3P']['dlogA_li_minus_exp'] < 0.10 and life['2s(2).2p(3).3p 3P']['z'] > 4
    # 615.8 nm upper level
    assert life['2s(2).2p(3).4d 5D']['dlogA_li_minus_exp'] > 0.1
    # direct UV gf: consistent with zero extra theory scatter at the experiments' precision
    assert s['direct_gf']['extra_theory_sigma_dex'] == 0.0
    # 1.8 um 3d-4f: components differ, multiplet totals agree within 0.02 dex
    for v in s['theory_vs_theory']['nf_multiplet_sums'].values():
        assert abs(v['d']) < 0.02
    # independent-theory scatter orders with Li's class
    sp = [s['theory_vs_theory']['by_li_class'][c]['per_method_robust_sigma'] for c in 'ACDE']
    assert sp == sorted(sp)
    assert s['status'] == 'CANDIDATE_HOLD_NOT_ADOPTED'

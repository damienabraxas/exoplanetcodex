"""RYA-1354 Fe I theory calibration (Kurucz K17 vs primary lab): source pinned, matching by level identity,
committed outputs reproduce."""
import gzip
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import rya1354_theory_calibration_fe as F  # noqa: E402


def test_kurucz_source_pinned_and_pure_theory():
    raw = gzip.open(F.KURUCZ).read()
    assert hashlib.sha256(raw).hexdigest() == 'f06b83d4f3afb68864734fe46a3515da0e878a28a4d85ef832d7a95f1ebbfc4f'
    rows, _ = F.kurucz()          # raises on any non-K17 (substituted) row
    assert len(rows) == 126288


def test_every_lab_line_matches_one_k17_line_by_both_levels():
    _, idx = F.kurucz()
    matched, cal = F.calibrate(idx)
    assert cal['n_lab'] == 493 and cal['n_matched'] == 493
    for r in matched:
        assert abs(r['d'] - (r['k17_loggf'] - r['lab_loggf'])) < 1e-9


def test_k_class_rules():
    assert F.k_class('(4F)5s e5F', '(4F)5p 5G', 0.5) == 'dS0|strong'
    assert F.k_class('(4F)4p y5F', '4s6D5s e7D', -8.3) == 'dS|weak'
    assert F.k_class('4s6D4d e7G', 's4f 3+[5+]', -2.4) == 'jK|mid'


def test_committed_outputs_equal_a_fresh_run(tmp_path, monkeypatch):
    committed = json.loads((F.OUT / 'summary.json').read_text())
    monkeypatch.setattr(F, 'OUT', tmp_path)
    assert json.loads(json.dumps(F.main())) == committed
    for name in ('k17_vs_lab', 'candidates'):
        assert (tmp_path / f'{name}.csv').read_text() == (ROOT / 'data/audit/rya1354/theory_calibration_fe' / f'{name}.csv').read_text()


def test_key_findings():
    s = json.loads((F.OUT / 'summary.json').read_text())
    cs = s['calibration']['class_sigma']
    # quality ladder: strong spin-conserving lines best; weak predicted lines unusable
    assert cs['dS0|strong']['sigma68'] < cs['dS0|mid']['sigma68'] < cs['dS|mid']['sigma68'] < cs['dS|weak']['sigma68']
    assert 0.08 < cs['dS0|strong']['bias'] < 0.2           # K17 too strong for its best class
    ir = s['calibration']['class_sigma_ir_gt_1um']['dS0|strong']
    assert ir['extra_theory_sigma_about_median'] < 0.03    # H band: after the bias, consistent with lab errors
    assert s['candidates']['by_tier_and_band']['T2_THEORETICAL_CODEX_CANDIDATE']['IR>1um'] > 200

"""RYA-1354 observed-spectrum screening of Theoretical Codex candidates. Needs the Mac atlas mirror (KP 1984,
IAG Baker 2020) and the repo CRIRES+ products; skips with the reason when the atlases are absent."""
import csv
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import rya1354_candidate_spectrum_check as S  # noqa: E402

HAVE = S.IAG_BAKER.exists() and S.KP_DIR.exists()
need = pytest.mark.skipif(not HAVE, reason='solar atlas mirror (KP 1984 / IAG Baker 2020) not on this host')


@need
def test_positive_control_asplund_fe_lines_pass_the_geometric_tests():
    A = S.load_atlases()
    ws = [float(r['wavelength_air_A']) for r in csv.DictReader(open(ROOT / 'data/reference/line_sets/fe_agss21_FeI.csv'))]
    ws += [7771.944, 7774.166, 7775.388]
    fails = []
    for w in ws:
        m = S.measure(S.pick(A, w), w, 0.3)
        fails += [f for f in m['fails'].split('+') if f]
    assert len(ws) == 43
    assert not {'offset', 'asymmetric', 'telluric', 'shallow'} & set(fails)
    assert fails.count('close_neighbour') <= 4


def test_multiplet_components_merge_into_one_feature():
    c = [dict(species='O I', air_A=6155.971, pred_depth=0.002), dict(species='O I', air_A=6155.989, pred_depth=0.002),
         dict(species='O I', air_A=6156.778, pred_depth=0.004)]
    f = S.features(c)
    assert [x['n_components'] for x in f] == [2, 1]
    assert abs(f[0]['pred_depth'] - (1 - 0.998 ** 2)) < 1e-4


@need
def test_committed_outputs_equal_a_fresh_run(tmp_path, monkeypatch):
    committed = json.loads((S.OUT / 'summary.json').read_text())
    monkeypatch.setattr(S, 'OUT', tmp_path)
    assert json.loads(json.dumps(S.main())) == committed
    assert (tmp_path / 'checked.csv').read_text() == (ROOT / 'data/audit/rya1354/candidate_spectrum_check/checked.csv').read_text()

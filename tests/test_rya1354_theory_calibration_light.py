"""RYA-1354 light-element theory-vs-lab (Mg MCHF, Al MCDHF/Kurucz/TOPbase)."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import rya1354_theory_calibration_light as L  # noqa: E402


def test_committed_summary_equals_fresh_run(tmp_path, monkeypatch):
    committed = json.loads((L.OUT / 'summary.json').read_text())
    monkeypatch.setattr(L, 'OUT', tmp_path)
    assert json.loads(json.dumps(L.main())) == committed


def test_transcription_matches_paper_row_count_and_split():
    mg = L.rows(ROOT / 'data/reference/mg_gf_lab/mg1_lab_loggf_pehlivan2017.csv')
    assert len(mg) == 26                      # 34 lab gf - 2 blended (authors) - 6 5f rows without calc
    assert sum(r['tau_theory'] == 'False' for r in mg) == 6


def test_ab_initio_beats_kurucz_bias_on_al():
    al = L.main()['results']['Al I (vs Burheim 2023 lab)']
    assert abs(al['Papoulia 2019 MCDHF']['median']) < 0.02
    assert al['Kurucz 1995']['median'] > 0.05

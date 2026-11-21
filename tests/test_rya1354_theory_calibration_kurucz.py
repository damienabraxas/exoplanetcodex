"""RYA-1354 multi-species Kurucz theory-vs-lab calibration: sources pinned, outputs reproduce, findings hold."""
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import rya1354_theory_calibration_kurucz as K  # noqa: E402

SHA = {
    'gf1400.pos.gz': 'b430ebee37a0378a64e6690ce3d8817daf1348351fca162801acaf6e38d84d0d',
    'gf2000.pos.gz': '2fa943aee4cf13dce819a09c1eec240d22158bc077e51314ce56b000e7e5d9b0',
    'gf2200.pos.gz': 'd7e30f7a055ee66852cf7764ab967ebf8e4a4b9581346a659889168f187d937c',
    'gf2201.pos.gz': 'c3bdf6800ec62580ff912a40f85a74285ea018c0e8dcb4dd9ecb384dd8621299',
    'gf2300.pos.gz': '11e3d136ec0fb15bcf5a94c21378fa26bbe4908bf321ad30b7dd8bf8655894f2',
    'gf2301.pos.gz': 'eb8d34ab3c402f99e5829945e42ca2ccaec6e4b49ec254205e8b6abc21607683',
    'gf2400.pos.gz': '940f8c849aa1e6405a468461bb56b0ca8fb4a8baa1028c1060faf14feb581e9f',
    'gf2401.pos.gz': '82578c40694ff1c57b428e670cb26052ce1576a868b70c8c61898d96e60a824c',
    'gf2500.pos.gz': '83be62d95fa5ab9487c641b621332ca8b2c2cc9d0cfb6f1edcf8bbfc4f6faf88',
    'gf2501.pos.gz': 'c8e15bf4e9d5dcbc3f0aa7ff883752b252c4bab0d83e6d51a1468a8c28bb094a',
    'gf2601.pos.gz': '484b0e6d0a8f999d4dc1bc52d278a0729383e1608efa1119085b56fb51e6776f',
    'gf2700.pos.gz': '31141da448c48c6b64bbd06881efca860987bb3de6168545bfece851a5f84a3a',
    'gf2701.pos.gz': 'a35b1cdcce82801ac83d4f8075c9dcfb0e3eb50ac5c5a368372d3f1ff1169eb2',
    'gf2800.pos.gz': 'ce893f14a0f8f9b59b9afce59490d7e978e82ba3c1968c24b43c1d1fd00ca77f',
}


def test_kurucz_files_pinned():
    for name, h in SHA.items():
        assert hashlib.sha256((ROOT / 'data/reference/kurucz' / name).read_bytes()).hexdigest() == h, name
    assert set(f'gf{c}.pos.gz' for c in K.SPECIES.values() if c != '2600') == set(SHA)


def test_committed_outputs_equal_a_fresh_run(tmp_path, monkeypatch):
    committed = json.loads((K.OUT / 'summary.json').read_text())
    monkeypatch.setattr(K, 'OUT', tmp_path)
    assert json.loads(json.dumps(K.main())) == committed
    for name in ('kurucz_vs_lab', 'candidates'):
        assert (tmp_path / f'{name}.csv').read_text() == (ROOT / 'data/audit/rya1354/theory_calibration_kurucz' / f'{name}.csv').read_text()


def test_fe_reproduces_the_fe_only_calibration():
    s = json.loads((K.OUT / 'summary.json').read_text())
    fe = json.loads((ROOT / 'data/audit/rya1354/theory_calibration_fe/summary.json').read_text())['calibration']['class_sigma']
    assert s['per_species']['Fe I']['n_matched'] == 493
    assert s['per_species']['Fe I']['classes'] == {k: v for k, v in fe.items()}


def test_key_findings():
    s = json.loads((K.OUT / 'summary.json').read_text())
    p = s['pooled_classes']
    assert p['dS0|strong']['sigma68'] < p['dS0|mid']['sigma68'] < p['dS|mid']['sigma68'] < p['dS|weak']['sigma68']
    for sp in ('Ti II', 'V II', 'Fe II'):       # ions: Kurucz best class good to <= 0.07 dex
        assert s['per_species'][sp]['classes']['dS0|strong']['sigma68'] <= 0.07
    # Cr II: main's Cr lab intake (2026-11-20 rebase) doubled the lab set (183 -> 358 lines; strong class n 62 -> 115)
    # and Kurucz's best class widened 0.039 -> 0.074 dex. Recorded as measured, not tuned away.
    cr = s['per_species']['Cr II']
    assert cr['n_lab'] >= 358 and cr['classes']['dS0|strong']['sigma68'] <= 0.08
    for sp, v in s['per_species'].items():
        assert v['n_unmatched'] <= 0.1 * v['n_lab'] + 1, sp
    assert s['candidates']['by_tier']['T2_THEORETICAL_CODEX_CANDIDATE'] > 1000

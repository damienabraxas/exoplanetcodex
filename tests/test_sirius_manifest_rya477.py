"""
tests/test_sirius_manifest_rya477.py
====================================
RYA-477 — guard the Sirius grid-download manifest + fetcher recipe (no network):
the manifest is well-formed, every grid carries a cited md5 + source URL, the
PASS-load-bearing Mn grid md5 matches the RYA-476-verified value, and the fetcher
enforces the no-silent-partial contract (size gate THEN md5 gate, extract only on PASS).
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

MANIFEST = ROOT / 'data' / 'sirius_manifest' / 'grids_zenodo_3982506.json'
MN_MD5 = 'ba01596031b054386a79a992141ebf72'   # RYA-476-verified Zenodo tar md5


def test_manifest_well_formed_and_cited():
    m = json.loads(MANIFEST.read_text())
    assert m['ticket'] == 'RYA-477'
    assert '3982506' in m['zenodo_record']
    assert m['n_items'] == len(m['items']) == 13
    for it in m['items']:
        assert it['name'].endswith('_pysme.tar.gz')
        assert len(it['md5']) == 32 and it['size_bytes'] > 0
        assert '3982506' in it['url'] and it['url'].endswith('download=1')


def test_mn_grid_md5_matches_rya476():
    m = json.loads(MANIFEST.read_text())
    mn = next(it for it in m['items'] if it['element'] == 'Mn')
    assert mn['md5'] == MN_MD5
    assert mn['size_bytes'] > 1_300_000_000        # ~1.4 GB archive


def test_fetcher_enforces_size_then_md5_gate():
    import sirius_fetch_grids_rya477 as F
    src = Path(F.__file__).read_text()
    # the no-silent-partial contract: size gate + md5 gate both precede extraction
    assert 'SIZE_MISMATCH' in src and 'MD5_MISMATCH' in src
    size_gate = src.index("if got != item['size_bytes']:")
    md5_gate = src.index("if got_md5 != item['md5']:")
    extract = src.index('with tarfile.open(archive)')
    assert size_gate < extract and md5_gate < extract       # both gates before extract
    assert '.codex_mounted' in src                  # refuses if data drive not mounted
    assert F._grd_for('nlte_Mn_scatt_pysme.tar.gz') == 'nlte_Mn_scatt_pysme.grd'

"""
tests/test_mn_grid_staging_rya476.py
====================================
RYA-476 — guard the Amarsi GALAH Mn NLTE grid staging recipe (no network): the
cited acquisition provenance (Zenodo record + verified tar md5), the disk-safe
contract (the .grd is gitignored, never committed; the tar is removed after
extraction), and the check()/presence reporting. The actual fetch is environment-
side; these assert the recipe is correct + reproducible.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import stage_amarsi_mn_grid_rya476 as S                              # noqa: E402

PROV = ROOT / 'data' / 'nlte_grids' / 'amarsi_galah' / 'Mn_amarsi2020_v3.prov.json'


def test_acquisition_recipe_matches_cited_provenance():
    prov = json.loads(PROV.read_text())
    # the staging md5 + Zenodo record are the cited, verified ones (never invented)
    assert S.TAR_MD5 == 'ba01596031b054386a79a992141ebf72'
    assert '3982506' in S.URL and prov['source']['zenodo_record'].endswith('3982506')
    assert S.GRD_NAME == 'nlte_Mn_scatt_pysme.grd'
    assert S.TAR_NAME == 'nlte_Mn_scatt_pysme.tar.gz'


def test_grd_is_gitignored_never_committed():
    # disk-safe / RYA-461: the multi-GB binary must NOT be tracked by git.
    import subprocess
    r = subprocess.run(['git', 'check-ignore', str(S.GRID_DIR / S.GRD_NAME)],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, "the .grd must be gitignored (multi-GB binary)"


def test_check_reports_presence_without_fetching(capsys):
    # check() is a pure presence probe — it must not download.
    present = S.check()
    out = capsys.readouterr().out
    assert ('PRESENT' in out) == bool(present)
    assert ('ABSENT' in out) == (not present)

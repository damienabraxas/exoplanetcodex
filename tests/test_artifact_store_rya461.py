"""
tests/test_artifact_store_rya461.py
===================================
RYA-461 — the artifact-preservation helper contract.

Pins the save-before-clean behaviour: save_artifact copies a gitignored artifact
into the canonical store, de-duplicates by md5, records provenance, never mutates
or deletes the source, and never clobbers a same-name different-content file.
"""
import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import artifact_store as A  # noqa: E402


@pytest.fixture
def store(tmp_path, monkeypatch):
    monkeypatch.setenv('CODEX_ARTIFACT_STORE', str(tmp_path))
    A.ensure_store()
    return tmp_path


def _write(p: Path, text: str) -> Path:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    return p


def test_kinds_and_store_layout(store):
    for k in A.KINDS:
        assert (store / k).is_dir()
    assert set(A.KINDS) == {'plots', 'data', 'diagnostics', 'grids', 'atlases'}


def test_save_copies_and_records_provenance(store, tmp_path):
    src = _write(tmp_path / 'wt' / 'results' / 'plots' / 'diag.png', 'PNGDATA')
    dest = A.save_artifact(src, 'plots', provenance={'worktree': 'wt-x', 'branch': 'b'})
    assert dest.exists() and dest.read_text() == 'PNGDATA'
    assert dest.parent == store / 'plots'
    assert src.exists()                                   # source untouched
    rows = list(csv.DictReader(open(store / A._MANIFEST_NAME)))
    assert len(rows) == 1
    assert rows[0]['md5'] == A.md5sum(src)
    assert rows[0]['source_worktree'] == 'wt-x' and rows[0]['source_branch'] == 'b'


def test_dedup_same_md5_not_copied_twice(store, tmp_path):
    a = _write(tmp_path / 'a' / 'x.png', 'SAME')
    b = _write(tmp_path / 'b' / 'x.png', 'SAME')          # identical content, different tree
    d1 = A.save_artifact(a, 'plots', provenance={'worktree': 'a', 'branch': 'b'})
    d2 = A.save_artifact(b, 'plots', provenance={'worktree': 'b', 'branch': 'b'})
    assert d1 == d2                                       # deduped to the same stored file
    assert len(list((store / 'plots').glob('*.png'))) == 1
    assert len(list(csv.DictReader(open(store / A._MANIFEST_NAME)))) == 1


def test_same_name_different_content_does_not_clobber(store, tmp_path):
    a = _write(tmp_path / 'a' / 'x.png', 'ONE')
    b = _write(tmp_path / 'b' / 'x.png', 'TWO')           # same name, different content
    d1 = A.save_artifact(a, 'plots')
    d2 = A.save_artifact(b, 'plots')
    assert d1 != d2                                       # disambiguated, not overwritten
    assert d1.read_text() == 'ONE' and d2.read_text() == 'TWO'
    assert len(list((store / 'plots').glob('*.png'))) == 2


def test_bad_kind_rejected(store, tmp_path):
    src = _write(tmp_path / 'x.csv', 'data')
    with pytest.raises(ValueError):
        A.save_artifact(src, 'not_a_kind')


def test_missing_source_raises(store):
    with pytest.raises(FileNotFoundError):
        A.save_artifact(store / 'nope.png', 'plots')

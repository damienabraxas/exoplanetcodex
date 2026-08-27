"""
tests/test_artifact_store_single_root_rya1090.py
================================================
RYA-1090 — the artifact store must have exactly ONE root, and must say so when it
does not.

The defect these pin: `store_root()` was `Path(__file__).resolve().parents[2]`, so
the "canonical" store was a function of which worktree the code ran from. Worktrees
moved from ~/Documents/Exoplanet Codex to ~/codex and the store forked in two — 85
artifacts unreachable, 27 stored twice, and a save-before-clean rescue that reported
success while the verifier read the other manifest (RYA-1087).

`test_store_root_is_identical_from_two_worktree_locations` is the direct regression
test: it FAILS against the old derivation and passes against the named constant.
"""
import csv
import importlib.util
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import artifact_store as A  # noqa: E402

SRC = ROOT / 'pipeline' / 'artifact_store.py'


def _load_from(worktree: Path, modname: str):
    """Import artifact_store as if it lived in `worktree`/pipeline/."""
    pkg = worktree / 'pipeline'
    pkg.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SRC, pkg / 'artifact_store.py')
    spec = importlib.util.spec_from_file_location(modname, pkg / 'artifact_store.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- the invariant

def test_store_root_is_identical_from_two_worktree_locations(tmp_path, monkeypatch):
    """THE regression test. Two worktrees under DIFFERENT parents must agree."""
    monkeypatch.delenv('CODEX_ARTIFACT_STORE', raising=False)
    a = _load_from(tmp_path / 'parentA' / 'wt1', 'as_a')
    b = _load_from(tmp_path / 'parentB' / 'wt2', 'as_b')
    assert a.store_root() == b.store_root() == A.CANONICAL_STORE_ROOT
    # and each still SEES its own parent as the legacy-derived root, which is what
    # keeps the guard honest after a future relocation
    assert a._legacy_derived_root() != b._legacy_derived_root()


def test_store_root_ignores_cwd(tmp_path, monkeypatch):
    monkeypatch.delenv('CODEX_ARTIFACT_STORE', raising=False)
    before = A.store_root()
    monkeypatch.chdir(tmp_path)
    assert A.store_root() == before


def test_env_override_still_wins(tmp_path, monkeypatch):
    monkeypatch.setenv('CODEX_ARTIFACT_STORE', str(tmp_path))
    assert A.store_root() == tmp_path


def test_canonical_root_is_not_inside_a_worktree_parent(monkeypatch):
    """~/codex is a worktree parent — a durable store must not live in one."""
    monkeypatch.delenv('CODEX_ARTIFACT_STORE', raising=False)
    assert A.CANONICAL_STORE_ROOT not in [p.expanduser() for p in A.LEGACY_STORE_ROOTS]


# ------------------------------------------------------------ divergence guard

def _stage_store(root: Path, n_rows: int, kind: str = 'plots') -> Path:
    for k in A.KINDS:
        (root / k).mkdir(parents=True, exist_ok=True)
    rows = []
    for i in range(n_rows):
        f = root / kind / f'a{i}.png'
        f.write_text(f'DATA{i}')
        rows.append({'date': '2026-01-01', 'kind': kind, 'stored_name': f.name,
                     'md5': A.md5sum(f), 'bytes': f.stat().st_size,
                     'source_worktree': 'wt', 'source_branch': 'br',
                     'source_relpath': f'/x/a{i}.png'})
    with open(root / A._MANIFEST_NAME, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=A._MANIFEST_HEADER)
        w.writeheader()
        w.writerows(rows)
    return root


def test_guard_fires_on_a_second_populated_store(tmp_path, monkeypatch):
    """POSITIVE CONTROL: stage a real fork, assert the guard actually complains."""
    monkeypatch.delenv('CODEX_ARTIFACT_STORE', raising=False)
    canon, legacy = tmp_path / 'canon', tmp_path / 'legacy'
    _stage_store(canon, 1)
    _stage_store(legacy, 3)
    monkeypatch.setattr(A, 'CANONICAL_STORE_ROOT', canon)
    monkeypatch.setattr(A, 'LEGACY_STORE_ROOTS', (legacy,))
    monkeypatch.setattr(A, '_legacy_derived_root', lambda: canon)

    assert A.divergent_stores() == [(legacy, 3)]
    with pytest.raises(A.StoreDivergenceError) as ei:
        A.assert_single_store()
    assert str(legacy) in str(ei.value) and '3 manifest rows' in str(ei.value)


def test_guard_quiet_when_only_one_store(tmp_path, monkeypatch):
    """NEGATIVE control — the guard must not cry wolf, or it will be ignored."""
    monkeypatch.delenv('CODEX_ARTIFACT_STORE', raising=False)
    canon, legacy = tmp_path / 'canon', tmp_path / 'legacy'
    _stage_store(canon, 2)
    (legacy / 'plots').mkdir(parents=True)              # exists but holds no manifest
    monkeypatch.setattr(A, 'CANONICAL_STORE_ROOT', canon)
    monkeypatch.setattr(A, 'LEGACY_STORE_ROOTS', (legacy,))
    monkeypatch.setattr(A, '_legacy_derived_root', lambda: canon)
    assert A.divergent_stores() == []
    A.assert_single_store()


def test_empty_legacy_manifest_is_not_divergence(tmp_path, monkeypatch):
    """A header-only manifest is a retired store, not a live fork."""
    monkeypatch.delenv('CODEX_ARTIFACT_STORE', raising=False)
    canon, legacy = tmp_path / 'canon', tmp_path / 'legacy'
    _stage_store(canon, 1)
    legacy.mkdir(parents=True)
    with open(legacy / A._MANIFEST_NAME, 'w', newline='') as fh:
        csv.DictWriter(fh, fieldnames=A._MANIFEST_HEADER).writeheader()
    monkeypatch.setattr(A, 'CANONICAL_STORE_ROOT', canon)
    monkeypatch.setattr(A, 'LEGACY_STORE_ROOTS', (legacy,))
    monkeypatch.setattr(A, '_legacy_derived_root', lambda: canon)
    assert A.divergent_stores() == []


def test_guard_skipped_under_explicit_override(tmp_path, monkeypatch):
    """An explicit $CODEX_ARTIFACT_STORE is deliberate; the guard targets accidents."""
    legacy = _stage_store(tmp_path / 'legacy', 3)
    monkeypatch.setenv('CODEX_ARTIFACT_STORE', str(tmp_path / 'canon'))
    monkeypatch.setattr(A, 'LEGACY_STORE_ROOTS', (legacy,))
    A.assert_single_store()                              # must not raise


def test_ensure_store_raises_while_forked(tmp_path, monkeypatch):
    monkeypatch.delenv('CODEX_ARTIFACT_STORE', raising=False)
    canon, legacy = tmp_path / 'canon', tmp_path / 'legacy'
    _stage_store(canon, 1)
    _stage_store(legacy, 2)
    monkeypatch.setattr(A, 'CANONICAL_STORE_ROOT', canon)
    monkeypatch.setattr(A, 'LEGACY_STORE_ROOTS', (legacy,))
    monkeypatch.setattr(A, '_legacy_derived_root', lambda: canon)
    with pytest.raises(A.StoreDivergenceError):
        A.ensure_store()
    A.ensure_store(check_divergence=False)               # the migration's escape hatch


def test_save_artifact_refuses_while_forked(tmp_path, monkeypatch):
    """The fork must stop preservation, not silently write to the wrong root."""
    monkeypatch.delenv('CODEX_ARTIFACT_STORE', raising=False)
    canon, legacy = tmp_path / 'canon', tmp_path / 'legacy'
    _stage_store(canon, 1)
    _stage_store(legacy, 2)
    monkeypatch.setattr(A, 'CANONICAL_STORE_ROOT', canon)
    monkeypatch.setattr(A, 'LEGACY_STORE_ROOTS', (legacy,))
    monkeypatch.setattr(A, '_legacy_derived_root', lambda: canon)
    src = tmp_path / 'new.png'
    src.write_text('X')
    with pytest.raises(A.StoreDivergenceError):
        A.save_artifact(src, 'plots')


def test_bridged_legacy_root_is_not_divergence(tmp_path, monkeypatch):
    """A legacy root SYMLINKED to the canonical store is the same store, not a fork.

    Worktrees that predate this fix still resolve store_root() to the old parent. The
    stopgap is to bridge that parent into the canonical store with symlinks; the guard
    must not then report the bridge as a fork, or it fires on the very mitigation.
    """
    monkeypatch.delenv('CODEX_ARTIFACT_STORE', raising=False)
    canon, legacy = tmp_path / 'canon', tmp_path / 'legacy'
    _stage_store(canon, 4)
    legacy.mkdir(parents=True)
    (legacy / A._MANIFEST_NAME).symlink_to(canon / A._MANIFEST_NAME)
    for k in A.KINDS:
        (legacy / k).symlink_to(canon / k, target_is_directory=True)
    monkeypatch.setattr(A, 'CANONICAL_STORE_ROOT', canon)
    monkeypatch.setattr(A, 'LEGACY_STORE_ROOTS', (legacy,))
    monkeypatch.setattr(A, '_legacy_derived_root', lambda: canon)

    assert A.divergent_stores() == []
    A.assert_single_store()
    # and a save through the BRIDGE lands in the canonical store, deduped against it
    src = tmp_path / 'x.png'
    src.write_text('DATA0')                       # same content as a staged row
    monkeypatch.setattr(A, 'CANONICAL_STORE_ROOT', legacy)   # old code's view
    dest = A.save_artifact(src, 'plots')
    assert dest.resolve().parent == (canon / 'plots').resolve()
    assert len(list(csv.DictReader(open(canon / A._MANIFEST_NAME)))) == 4   # deduped

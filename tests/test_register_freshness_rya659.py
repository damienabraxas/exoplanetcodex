"""
tests/test_register_freshness_rya659.py
=======================================
RYA-659 — register-freshness guard + state-surface registry.

The register drifted eleven state-changing tickets (RYA-556..652) because the
obligation "bump the register" lived only in prose: RYA-643 changed element state
(data/audit/element_status_tracker.csv) and did not touch CODEX_STATE_REGISTER.md.
These tests pin the guard that makes that drift loud.

The teeth are proved on a DISPOSABLE git repo built in a tmp dir, not on this
checkout — the guard reads git history, so proving it fails requires a real commit
touching a state surface without the register. Note that `touch`-ing a file does
NOT trip the guard by design: an uncommitted mtime is not a state change.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from pipeline import state_surfaces as ss

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_register_freshness.py"


# --------------------------------------------------------------------------
# The registry itself (single source of truth)
# --------------------------------------------------------------------------

def test_surface_registry_is_non_empty_and_unique():
    paths = [s.path for s in ss.STATE_SURFACES]
    assert paths, "the state-surface registry must not be empty"
    assert len(paths) == len(set(paths)), "duplicate surface path in the registry"


def test_every_ledger_is_also_a_state_surface():
    # LEDGERS.md members are a subset of the surfaces the guard watches.
    assert set(ss.LEDGER_PATHS) <= set(ss.SURFACE_PATHS)
    assert ss.LEDGER_PATHS, "at least the element tracker + catalogs are ledgers"


def test_registered_surfaces_exist_in_this_checkout():
    # A surface that has been renamed away would silently stop being watched.
    missing = [s.path for s in ss.STATE_SURFACES if not (ss.REPO_ROOT / s.path).exists()]
    assert not missing, f"registered surfaces missing from the repo: {missing}"


def test_register_and_ledgers_index_exist():
    assert (ss.REPO_ROOT / ss.REGISTER).exists()
    assert (ss.REPO_ROOT / ss.LEDGERS_INDEX).exists()


def test_ledgers_index_names_every_ledger_path():
    # LEDGERS.md must actually point at the ledger files, not drift from them.
    text = (ss.REPO_ROOT / ss.LEDGERS_INDEX).read_text()
    for path in ss.LEDGER_PATHS:
        assert path in text, f"{ss.LEDGERS_INDEX} does not name {path}"


# --------------------------------------------------------------------------
# The guard's teeth, on a disposable repo
# --------------------------------------------------------------------------

def _git(repo: Path, *args: str, when: str | None = None) -> None:
    """Run git in `repo`, optionally pinning the commit timestamp.

    Timestamps must be pinned: git commit times have 1-second resolution, so two
    commits made in the same second are indistinguishable to the guard (which
    treats "same second" as fresh — the bias that lets a single commit touching
    both a surface and the register pass).
    """
    env = None
    if when is not None:
        env = {**os.environ, "GIT_AUTHOR_DATE": when, "GIT_COMMITTER_DATE": when}
    subprocess.run(
        ["git", *args], cwd=repo, check=True, capture_output=True, text=True, env=env
    )


@pytest.fixture
def fake_repo(tmp_path: Path) -> Path:
    """A minimal git repo carrying the register + one state surface."""
    repo = tmp_path / "codex"
    surface = ss.LEDGER_PATHS[0]          # the element tracker
    (repo / Path(surface).parent).mkdir(parents=True)
    (repo / surface).write_text("element,status\nFe,PASS\n")
    (repo / ss.REGISTER).write_text("# Codex State Register\n\n**Version: v1**\n")

    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "config", "user.email", "test@example.com")
    _git(repo, "config", "user.name", "test")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "seed: register + surface together",
         when="2026-07-01T12:00:00")
    return repo


def _guard_at(root: Path, *argv: str) -> subprocess.CompletedProcess:
    """Run the guard's full main() with its ROOT pointed at `root`.

    Goes through main() (not history_mode directly) so the arg parsing and the
    GitError -> exit 2 contract are exercised as shipped.
    """
    code = (
        "import sys\n"
        f"sys.path.insert(0, {str(SCRIPT.parents[1])!r})\n"
        "import scripts.check_register_freshness as g\n"
        "from pathlib import Path\n"
        f"g.ROOT = Path({str(root)!r})\n"
        f"sys.argv = ['check_register_freshness.py', *{list(argv)!r}]\n"
        "sys.exit(g.main())\n"
    )
    return subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True
    )


def test_passes_when_register_committed_with_the_surface(fake_repo):
    r = _guard_at(fake_repo)
    assert r.returncode == 0, r.stderr
    assert "Register freshness OK" in r.stdout


def test_fails_when_a_surface_moves_without_the_register(fake_repo):
    # The RYA-643 hole, reproduced: element state changes, register does not.
    surface = ss.LEDGER_PATHS[0]
    (fake_repo / surface).write_text("element,status\nFe,PASS\nCo,PASS\n")
    _git(fake_repo, "commit", "-q", "-a", "-m", "change element state only",
         when="2026-07-10T12:00:00")

    r = _guard_at(fake_repo)
    assert r.returncode == 1, f"guard must FAIL on drift; stdout={r.stdout}"
    assert "REGISTER FRESHNESS GUARD FAILED" in r.stderr
    assert surface in r.stderr


def test_recovers_once_the_register_is_bumped(fake_repo):
    surface = ss.LEDGER_PATHS[0]
    (fake_repo / surface).write_text("element,status\nFe,PASS\nCo,PASS\n")
    _git(fake_repo, "commit", "-q", "-a", "-m", "change element state only",
         when="2026-07-10T12:00:00")

    (fake_repo / ss.REGISTER).write_text("# Codex State Register\n\n**Version: v2**\n")
    _git(fake_repo, "commit", "-q", "-a", "-m", "bump register to v2",
         when="2026-07-11T12:00:00")

    r = _guard_at(fake_repo)
    assert r.returncode == 0, r.stderr


def test_touch_alone_does_not_trip_the_guard(fake_repo):
    # Documented behaviour: the guard reads git history, never mtimes.
    (fake_repo / ss.LEDGER_PATHS[0]).touch()
    r = _guard_at(fake_repo)
    assert r.returncode == 0, "an uncommitted mtime is not a state change"


# --------------------------------------------------------------------------
# PR mode (--since-main)
# --------------------------------------------------------------------------

def test_pr_mode_fails_when_state_changes_without_the_register(fake_repo):
    _git(fake_repo, "branch", "-q", "base")
    (fake_repo / ss.LEDGER_PATHS[0]).write_text("element,status\nFe,PASS\nCo,PASS\n")
    _git(fake_repo, "commit", "-q", "-a", "-m", "PR: element state only",
         when="2026-07-10T12:00:00")

    r = _guard_at(fake_repo, "--since-main", "base")
    assert r.returncode == 1, f"PR mode must FAIL; stdout={r.stdout}"
    assert "PR changes state but not the register" in r.stderr
    assert ss.LEDGER_PATHS[0] in r.stderr


def test_pr_mode_passes_when_the_register_moves_too(fake_repo):
    _git(fake_repo, "branch", "-q", "base")
    (fake_repo / ss.LEDGER_PATHS[0]).write_text("element,status\nFe,PASS\nCo,PASS\n")
    (fake_repo / ss.REGISTER).write_text("# Codex State Register\n\n**Version: v2**\n")
    _git(fake_repo, "commit", "-q", "-a", "-m", "PR: element state + register bump",
         when="2026-07-10T12:00:00")

    r = _guard_at(fake_repo, "--since-main", "base")
    assert r.returncode == 0, r.stderr
    assert "register updated" in r.stdout


def test_pr_mode_passes_when_no_state_surface_is_touched(fake_repo):
    _git(fake_repo, "branch", "-q", "base")
    (fake_repo / "notes.txt").write_text("not a state surface\n")
    _git(fake_repo, "add", "-A")
    _git(fake_repo, "commit", "-q", "-m", "PR: docs only",
         when="2026-07-10T12:00:00")

    r = _guard_at(fake_repo, "--since-main", "base")
    assert r.returncode == 0, r.stderr
    assert "0 state surface(s) changed" in r.stdout


def test_git_failure_is_loud_not_a_pass(tmp_path):
    # A non-repo directory must ERROR (exit 2), never silently report OK.
    r = _guard_at(tmp_path)
    assert r.returncode == 2, f"expected loud error, got {r.returncode}: {r.stdout}"
    assert "ERRORED (not a pass)" in r.stderr

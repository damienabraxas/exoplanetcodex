"""
RYA-686 -- a result artifact must not land without its generating harness.

Two halves:
  * BOTH DIRECTIONS, hermetically. A throwaway git repo is built with a result JSON and
    no harness (guard must FAIL), then the harness and its manifest entry are added
    (guard must PASS). The guard reads `git ls-files`, so these exercise the real code
    path rather than a mocked one.
  * THE LIVE MANIFEST. The real repo must be clean, and the set of artifacts exempted
    from having a generator is FROZEN here. Adding an exemption means editing this test
    -- which is the point: RYA-559 happened silently, and the next one should not be able
    to.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.check_result_generators import (  # noqa: E402
    MANIFEST,
    RESULTS_DIR,
    check,
    declared_long_flags,
    invocation_long_flags,
    load_manifest,
    repo_root,
)

ROOT = repo_root()

# The COMPLETE set of landed artifacts that have no committed generator, frozen as of
# RYA-686. Every entry is a defect or a deliberate hand-authored record, never a
# convenience. To add one you must edit this test AND say why in the manifest note.
FROZEN_UNREPRODUCIBLE = {"solar_ba_synthesis_rya559.json"}
FROZEN_HAND_AUTHORED = {
    "sr2_line_selection_rya430.json",
    "rya342_corrected_solar_fe.txt",
    # RYA-939, added deliberately when the guard was extended to data/processed/.
    # Gaia-ESO Survey pre-stored solar EWs, transcribed under RYA-196. It has no
    # generator because no code in this repo produced the numbers -- inventing one
    # would be the fabrication this manifest exists to prevent. It had been tracked
    # and entirely unguarded since RYA-196; extending the scan is what surfaced it.
    "data/processed/solar_ew_ges_reference.csv",
    # RYA-847, added deliberately. Neither is a program's output and neither can be:
    #   README.md          reading instructions for the directory, including the warning
    #                      that these cells must not be diffed against the rya845-era
    #                      published matrix.
    #   pregate_control    the sweep's PRE-GATE numbers. Its generator existed and its
    #                      output was destroyed by a tree re-sync, so the numbers are
    #                      carried as data with their source named rather than invented a
    #                      harness for after the fact. They are not merely asserted:
    #                      scripts/rya847_gate_diff.py checks every row against the
    #                      artifact it describes (n_post + n_caught == n_pre) and fails if
    #                      a control row finds no product.
    "rya847/README.md",
    "rya847/rya847_pregate_control.csv",
    # RYA-925: narrative/reporting surfaces, not machine-generated measurements.
    # Their machine sources remain registered separately in GENERATORS.yaml.
    "rya925/REPORT.md",
    "rya925/tracker.html",
}


# ── Hermetic both-directions demonstration ────────────────────────────────────────────

def _git(repo: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(repo), *args], check=True,
                   capture_output=True, text=True)


@pytest.fixture()
def fake_repo(tmp_path: Path) -> Path:
    """A minimal git repo carrying one result artifact and an empty manifest."""
    repo = tmp_path / "repo"
    (repo / RESULTS_DIR).mkdir(parents=True)
    (repo / "scripts").mkdir()
    (repo / RESULTS_DIR / "widget_synthesis_rya999.json").write_text(
        json.dumps({"A_LTE": 1.234})
    )
    (repo / MANIFEST).write_text("artifacts: []\n")
    _git(repo, "init", "-q")
    _git(repo, "add", "-A")
    return repo


def test_guard_fails_when_a_result_lands_without_its_harness(fake_repo: Path):
    """DIRECTION 1 -- this is the RYA-559 landing, replayed."""
    problems = check(fake_repo)
    assert len(problems) == 1
    assert "widget_synthesis_rya999.json" in problems[0]
    assert "NO entry" in problems[0]


def test_guard_still_fails_when_the_entry_names_a_harness_that_is_not_there(fake_repo: Path):
    """An entry is not a substitute for the file. Pointing at a missing script fails."""
    (fake_repo / MANIFEST).write_text(
        "artifacts:\n"
        "  - artifact: widget_synthesis_rya999.json\n"
        "    generator: scripts/rya999_widget_synth_sirius.py\n"
        "    ticket: RYA-999\n"
    )
    _git(fake_repo, "add", "-A")
    problems = check(fake_repo)
    assert len(problems) == 1
    assert "does not exist in the repo" in problems[0]


def test_guard_passes_once_the_harness_is_committed_with_an_entry(fake_repo: Path):
    """DIRECTION 2 -- same artifact, harness present and linked. Clean."""
    (fake_repo / "scripts" / "rya999_widget_synth_sirius.py").write_text(
        "import argparse\n"
        "def main():\n"
        "    ap = argparse.ArgumentParser()\n"
        "    ap.add_argument('--deblend', action='store_true')\n"
        "    ap.parse_args()\n"
    )
    (fake_repo / MANIFEST).write_text(
        "artifacts:\n"
        "  - artifact: widget_synthesis_rya999.json\n"
        "    generator: scripts/rya999_widget_synth_sirius.py\n"
        "    invocation: python scripts/rya999_widget_synth_sirius.py --deblend\n"
        "    ticket: RYA-999\n"
    )
    _git(fake_repo, "add", "-A")
    assert check(fake_repo, check_invocations=True) == []


def test_guard_catches_the_rya560_case_flag_documented_but_not_declared(fake_repo: Path):
    """A harness that is committed but not invocable AS DOCUMENTED.

    scripts/rya560_zr2_synth_sirius.py shipped with no argparse at all, so RYA-585's
    `--deblend` entrypoint did not exist. This is the static form of that check: the
    recorded invocation names a flag the source does not declare.
    """
    (fake_repo / "scripts" / "rya999_widget_synth_sirius.py").write_text(
        "print('no argparse here')\n"
    )
    (fake_repo / MANIFEST).write_text(
        "artifacts:\n"
        "  - artifact: widget_synthesis_rya999.json\n"
        "    generator: scripts/rya999_widget_synth_sirius.py\n"
        "    invocation: python scripts/rya999_widget_synth_sirius.py --deblend\n"
        "    ticket: RYA-999\n"
    )
    _git(fake_repo, "add", "-A")

    # Without --check-invocations the harness is merely present: clean.
    assert check(fake_repo) == []

    problems = check(fake_repo, check_invocations=True)
    assert len(problems) == 1
    assert "declares no argparse arguments at all" in problems[0]


def test_guard_rejects_an_unexplained_exemption(fake_repo: Path):
    """UNREPRODUCIBLE without a note is a hiding place, not a record."""
    (fake_repo / MANIFEST).write_text(
        "artifacts:\n"
        "  - artifact: widget_synthesis_rya999.json\n"
        "    status: UNREPRODUCIBLE\n"
        "    generator: null\n"
        "    ticket: RYA-999\n"
    )
    _git(fake_repo, "add", "-A")
    problems = check(fake_repo)
    assert len(problems) == 1
    assert "requires a `note:`" in problems[0]


def test_guard_rejects_a_stale_entry(fake_repo: Path):
    """The reverse direction: an entry for an artifact that is not tracked."""
    (fake_repo / MANIFEST).write_text(
        "artifacts:\n"
        "  - artifact: widget_synthesis_rya999.json\n"
        "    status: UNREPRODUCIBLE\n"
        "    generator: null\n"
        "    ticket: RYA-999\n"
        "    note: harness never committed\n"
        "  - artifact: deleted_thing.json\n"
        "    generator: scripts/whatever.py\n"
        "    ticket: RYA-999\n"
    )
    _git(fake_repo, "add", "-A")
    problems = check(fake_repo)
    assert any("deleted_thing.json" in p and "not tracked" in p for p in problems)


# ── Static helpers ────────────────────────────────────────────────────────────────────

def test_declared_long_flags_reads_argparse_statically(tmp_path: Path):
    src = tmp_path / "h.py"
    src.write_text(
        "import argparse\n"
        "ap = argparse.ArgumentParser()\n"
        "ap.add_argument('-o', '--out', default=None)\n"
        "ap.add_argument('--only')\n"
        "ap.add_argument('positional')\n"
    )
    flags, has_argparse = declared_long_flags(src)
    assert has_argparse is True
    assert flags == {"--out", "--only"}

    bare = tmp_path / "b.py"
    bare.write_text("print('hi')\n")
    assert declared_long_flags(bare) == (set(), False)


def test_invocation_long_flags_strips_values():
    got = invocation_long_flags("python x.py --only 4077 --out=a.json pos -v")
    assert got == {"--only", "--out"}


# ── The live repo ─────────────────────────────────────────────────────────────────────

def test_live_repo_is_clean():
    problems = check(ROOT, check_invocations=True)
    assert problems == [], "\n".join(problems)


def test_exempted_artifacts_are_exactly_the_frozen_set():
    """The audit, pinned. Widening it must be a deliberate, reviewed edit."""
    entries = load_manifest(ROOT)
    unreproducible = {e["artifact"] for e in entries
                      if e.get("status") == "UNREPRODUCIBLE"}
    hand_authored = {e["artifact"] for e in entries
                     if e.get("status") == "HAND_AUTHORED"}
    assert unreproducible == FROZEN_UNREPRODUCIBLE
    assert hand_authored == FROZEN_HAND_AUTHORED


def test_rya559_ba_harness_is_still_absent_and_recorded_as_such():
    """Guards the honesty of the RYA-559 record itself.

    If someone ever commits a `scripts/rya559_*` file, this fails -- and the right
    response is to determine whether it is genuinely the original harness (it is not,
    unless recovered from outside this repo) rather than to quietly relabel the entry
    COMMITTED. RYA-581's rebuild generates ITS OWN artifact, not this one.
    """
    assert not list((ROOT / "scripts").glob("rya559*")), (
        "a scripts/rya559* file now exists -- re-adjudicate the UNREPRODUCIBLE entry "
        "for solar_ba_synthesis_rya559.json rather than assuming it is the original"
    )
    entry = next(e for e in load_manifest(ROOT)
                 if e["artifact"] == "solar_ba_synthesis_rya559.json")
    assert entry["status"] == "UNREPRODUCIBLE"
    assert entry["generator"] is None

"""RYA-1195 — the CNO ingest must not write a clock reading into a tracked file.

`scripts/ingest_amarsi2021_cno_rya1136.py` stamped `retrieved_at = str(date.today())` into
`data/reference/amarsi2021_cno/manifest.json`, which git tracks. Every full-suite run
dirtied the tree, because `tests/test_cno_intake_rya1136.py` shells out to the ingest.

🔴 IT ALREADY GOT COMMITTED, WHICH IS THE WHOLE ARGUMENT. `c314879` (2026-08-30, RYA-1136)
vendored the raw files and recorded the true retrieval date; `df17f8b` -- RYA-1178, about
`Fe.json`'s publication schema and nothing to do with CNO -- carried a stray 2026-09-03
into the record and it was merged. A non-deterministic write to a tracked file is not a
hygiene annoyance; it is a provenance corruption with a ride-along vector.
"""
from __future__ import annotations

import ast
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
INGEST = ROOT / "scripts/ingest_amarsi2021_cno_rya1136.py"
MANIFEST = ROOT / "data/reference/amarsi2021_cno/manifest.json"
HOLDING = "data/reference/amarsi2021_cno"
#: Every tracked path this ingest writes. The audit dir is included deliberately: the
#: ingest also writes `summary.json` there, and the clean-tree claim is worth nothing if
#: it only covers the directory the bug happened to be in.
INGEST_WRITES = (HOLDING, "data/audit/rya1136_cno_intake")

#: The date the vendoring commit recorded. Not "when the test ran".
TRUE_RETRIEVAL = "2026-08-30"


def _git(*args) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True,
                          text=True).stdout.strip()


def test_the_ingest_leaves_the_tracked_tree_clean_across_two_runs():
    """The spec's smoke test. Run it twice; `git status` clean both times."""
    assert _git("status", "--porcelain", *INGEST_WRITES) == "", (
        "a path this ingest writes was already dirty before the test ran")
    for run in (1, 2):
        r = subprocess.run([sys.executable, str(INGEST)], cwd=ROOT,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
        dirty = _git("status", "--porcelain", *INGEST_WRITES)
        assert dirty == "", f"run {run} dirtied a tracked file:\n{dirty}"


def test_the_ingest_writes_no_clock_reading():
    """AST, not grep: the prose above says `date.today` and must not trip this, and a
    docstring promise must not satisfy it either."""
    tree = ast.parse(INGEST.read_text())
    assert not [n for n in ast.walk(tree)
                if isinstance(n, ast.Import)
                for a in n.names if a.name == "datetime"], "datetime imported again"
    bad = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Call):
            fn = getattr(n.func, "attr", None) or getattr(n.func, "id", None)
            if fn in {"today", "now", "utcnow", "time"}:
                bad.append(n.lineno)
    assert not bad, (
        f"the ingest reads the clock at line(s) {bad}. `raw/` is vendored and this script "
        f"makes no network call, so there is nothing to timestamp -- a re-run would "
        f"overwrite a real acquisition date with a test-run date, in a tracked file.")


def test_the_recorded_date_is_the_one_the_vendoring_commit_wrote():
    """A static date is only better than a clock if it is the TRUE date. Anchored to the
    commit that actually added the raw files, so the constant cannot drift into fiction."""
    got = json.loads(MANIFEST.read_text())["retrieved_at"]
    assert got == TRUE_RETRIEVAL

    added = _git("log", "--diff-filter=A", "--format=%ad", "--date=short",
                 "--", f"{HOLDING}/raw/table2.dat")
    if added:                      # skip gracefully in a shallow or export-only checkout
        assert added.splitlines()[-1] == TRUE_RETRIEVAL, (
            f"the raw table was vendored on {added}, but the manifest claims "
            f"{TRUE_RETRIEVAL}. One of the two is wrong.")

    src = INGEST.read_text()
    assert f'RETRIEVED_AT = "{TRUE_RETRIEVAL}"' in src, (
        "the constant and the artifact must not drift apart (RYA-1101)")


def test_the_manifest_stays_tracked():
    """The alternative fixes were to gitignore the manifest or move the field to a
    sidecar. Both were wrong: this file carries the DOI, the bibcode and the sha256 of
    the vendored raw files, so it is a provenance record and belongs in version control.
    Only the clock reading was the problem."""
    tracked = _git("ls-files", str(MANIFEST.relative_to(ROOT)))
    assert tracked, "the manifest must stay tracked — it is provenance, not scratch"
    doc = json.loads(MANIFEST.read_text())
    for key in ("doi", "ads_bibcode", "source_files", "citation"):
        assert doc.get(key), f"{key} missing — this is why the file is tracked"
    assert set(doc["source_files"]) == {"raw/ReadMe", "raw/table2.dat"}


def test_the_guard_would_fire_on_the_defect_it_was_written_for(tmp_path):
    """Mutation control: the AST scan must object to the exact shape that was fixed."""
    bad = tmp_path / "bad.py"
    bad.write_text("from datetime import date\n"
                   "m = {'retrieved_at': str(date.today())}\n")
    tree = ast.parse(bad.read_text())
    hits = [n.lineno for n in ast.walk(tree)
            if isinstance(n, ast.Call)
            and (getattr(n.func, "attr", None) or getattr(n.func, "id", None))
            in {"today", "now", "utcnow", "time"}]
    assert hits, "the scanner does not detect date.today()"

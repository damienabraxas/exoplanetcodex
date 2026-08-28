"""
tests/test_line_key_guard_rya1037.py — RYA-1037
===============================================
Keying a line match on wavelength alone has been fixed one call site at a time across seven
tickets and keeps coming back. The enforcement is the point of this ticket, so the controls
matter more than the check: **a guard that has never been shown to fail proves nothing.**
"""
from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.line_match import LineMatchError, match  # noqa: E402

AUDIT = ROOT / "scripts" / "audit_line_keys_rya1037.py"


def _run(root: Path):
    return subprocess.run([sys.executable, str(AUDIT), "--check", "--root", str(root)],
                          capture_output=True, text=True)


def _tree(tmp_path: Path, body: str, name: str = "offender.py") -> Path:
    """A minimal repo-shaped tree with one staged offender."""
    for d in ("pipeline", "scripts", "config"):
        (tmp_path / d).mkdir(parents=True, exist_ok=True)
    (tmp_path / "pipeline" / name).write_text(textwrap.dedent(body))
    return tmp_path


# ── the guard, in both directions ────────────────────────────────────────────

def test_the_real_tree_passes():
    """Every finding on main is either fixed or waived WITH a ticket. If this fails, a new
    wavelength-only key was just introduced — which is the whole point."""
    r = _run(ROOT)
    assert r.returncode == 0, r.stdout[-3000:]


def test_control_a_wavelength_only_merge_FAILS(tmp_path):
    """POSITIVE CONTROL 1, required by the ticket."""
    t = _tree(tmp_path, '''
        import pandas as pd
        def join(a, b):
            return a.merge(b, on="wavelength_air_A", how="left")
    ''')
    r = _run(t)
    assert r.returncode == 1, r.stdout
    assert "WAVE_ONLY_MERGE" in r.stdout


def test_control_a_rounded_wavelength_key_FAILS(tmp_path):
    """POSITIVE CONTROL 2, required by the ticket. RYA-1033: round() and np.round()
    disagree on ...x5, so a rounded wavelength is not stable even against itself."""
    t = _tree(tmp_path, '''
        def build(rows):
            out = {}
            for r in rows:
                out[round(float(r["wavelength_air_A"]), 3)] = r
            return out
    ''')
    r = _run(t)
    assert r.returncode == 1, r.stdout
    assert "ROUNDED_WAVE_KEY" in r.stdout


def test_control_nearest_wins_FAILS(tmp_path):
    """argmin on a wavelength distance is first-hit-in-file-order wearing another hat."""
    t = _tree(tmp_path, '''
        import numpy as np
        def pick(cand, wave_A):
            return cand.iloc[np.abs(cand["wavelength_air_A"] - wave_A).argmin()]
    ''')
    r = _run(t)
    assert r.returncode == 1 and "ARGMIN_ON_WAVE" in r.stdout


def test_control_rounding_for_DISPLAY_is_not_flagged(tmp_path):
    """The counter-control. The rule is about a rounded value USED AS A KEY — flagging every
    round(wave, 1) written for a report gave 100 findings, almost all display formatting,
    and a check that cries wolf gets waived into uselessness."""
    t = _tree(tmp_path, '''
        def report(w):
            return {"coverage_A": [round(float(w.min()), 1), round(float(w.max()), 1)]}
    ''')
    assert _run(t).returncode == 0


def test_control_a_band_filter_is_not_flagged(tmp_path):
    """Selecting a WINDOW is not identifying a LINE."""
    t = _tree(tmp_path, '''
        def trim(sp, lo, hi):
            return sp[(sp.wave_A >= lo) & (sp.wave_A <= hi)]
    ''')
    assert _run(t).returncode == 0


def test_control_a_new_occurrence_in_a_WAIVED_file_still_fails(tmp_path):
    """Waivers are keyed (kind, file) WITH A COUNT, so a waived file cannot absorb new
    defects. Line numbers would be brittle; a count survives edits and still refuses."""
    for d in ("pipeline", "config"):
        (tmp_path / d).mkdir(parents=True, exist_ok=True)
    (tmp_path / "pipeline" / "x.py").write_text(textwrap.dedent('''
        import pandas as pd
        def a(a_, b): return a_.merge(b, on="wavelength_air_A")
        def b(a_, b): return a_.merge(b, on="wavelength_air_A")
    '''))
    (tmp_path / "config" / "line_key_waivers.yaml").write_text(textwrap.dedent('''
        waivers:
          - kind: WAVE_ONLY_MERGE
            file: pipeline/x.py
            count: 1
            ticket: RYA-1037
            reason: one known site
    '''))
    r = _run(tmp_path)
    assert r.returncode == 1
    assert "MORE findings than the waiver allows" in r.stdout


def test_a_waiver_without_a_ticket_is_rejected(tmp_path):
    """A waiver without an owner is a suppression, and suppression is how this defect class
    survived seven tickets."""
    for d in ("pipeline", "config"):
        (tmp_path / d).mkdir(parents=True, exist_ok=True)
    (tmp_path / "config" / "line_key_waivers.yaml").write_text(
        "waivers:\n  - kind: WAVE_ONLY_MERGE\n    file: pipeline/x.py\n")
    r = _run(tmp_path)
    assert r.returncode != 0
    assert "suppression" in r.stdout + r.stderr


# ── the canonical matcher — RYA-1033's, extended, NOT a second one ───────────
#
# 🔴 `pipeline/line_match.py` ALREADY EXISTED. RYA-1037 asks for "one canonical matcher" and
# the honest answer to that ask is that RYA-1033 already built it — the work here is the
# AUDIT and the ENFORCEMENT, plus one narrow gap closed below. Writing a second matcher
# would have been the "two homes" defect (RYA-954/368) committed while fixing another one.

def test_the_canonical_matcher_is_rya1033s_and_there_is_only_one():
    """A second module named like a matcher is the defect, not the fix."""
    found = sorted(p.name for p in (ROOT / "pipeline").glob("*line_match*.py"))
    assert found == ["line_match.py"], found
    src = (ROOT / "pipeline" / "line_match.py").read_text()
    assert "RYA-1033" in src


def test_strict_mode_refuses_a_LONE_candidate_without_EP():
    """The gap RYA-1037 closes. `match()` records AMBIGUOUS for >1 candidate, but a SINGLE
    candidate with no EP resolved on wavelength alone — and one candidate is not the same as
    the right transition. RYA-853: crosscheck_nist stamped a grade on Fe I 6065.490 from a
    lone in-window row whose EP was 2.608 eV against our line's 4.956."""
    with pytest.raises(LineMatchError):
        match([6149.246], [6149.2465], require_ep=True)


def test_strict_mode_is_OPT_IN_so_existing_callers_are_unchanged():
    """Additive by design: RYA-1033's ten tests still pass, and no existing behaviour moves.
    A keying fix that silently changed which lines match would change abundances."""
    assert list(match([6149.246], [6149.2465]).index) == [0]
    assert list(match([6149.246], [6149.2465], want_ep=[3.889], src_ep=[3.889],
                      require_ep=True).index) == [0]


def test_more_than_one_candidate_is_still_recorded_not_resolved():
    """RYA-1033's own contract, pinned here because RYA-1037 depends on it."""
    assert match([6149.246], [6149.2445, 6149.2475]).ambiguous


def test_the_matcher_never_rounds_a_wavelength():
    """AST, not a substring search — the docstring quotes `round(6136.615, 2)` to explain
    the defect, and a naive text check flags its own explanation."""
    import ast
    tree = ast.parse((ROOT / "pipeline" / "line_match.py").read_text())
    calls = [n for n in ast.walk(tree)
             if isinstance(n, ast.Call)
             and (getattr(n.func, "id", "") == "round"
                  or getattr(n.func, "attr", "") in ("round", "around"))]
    assert not calls, f"line_match.py rounds at line(s) {[c.lineno for c in calls]}"

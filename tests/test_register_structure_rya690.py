"""
tests/test_register_structure_rya690.py
=======================================
RYA-690 — the register's `**Version:` header must be single-valued.

WHY THIS GUARD EXISTS (ratified by Ryan 2026-08-08: "add the guard, three times is
enough"). Three consecutive tickets produced the identical defect — RYA-673, then
RYA-684, then RYA-682 — each bumping the header without adding a changelog row. Two
of the three happened *after* the ticket describing the defect was already filed. By
the time RYA-690 ran, `origin/main` carried **eight** header lines, `v36` was claimed
by two different tickets, and three landings had no changelog row of their own.

It recurred because the correct merge resolution is not the intuitive one. Every
ticket edits the header at the same anchor, so concurrent branches always conflict,
and "accept both" is CORRECT for an append-only log (SEQUENCE.md) but WRONG for a
single-valued header. Nothing signalled which was which — so the safe-looking local
choice propagated the bug. This guard is the signal.

The four assertions, each traced to an observed failure:

  1. exactly one header line                    — eight were on main
  2. no changelog version appears twice         — `v29` identified two things
  3. the header's version has a changelog row   — the orphan class
  4. the header names the NEWEST row            — closes 3's blind spot: RYA-673's
     header said v29 and RYA-675's said v36; both versions existed, but owned by
     RYA-581 and RYA-676, so "is there a row?" passed while the landing was orphaned

Deliberately NOT asserted, and both were checked against the real file first:
contiguity (the changelog legitimately skips v25) and ordering (it already reads
v36, v35, v34, v31, v32 — non-monotonic). A guard that failed on legitimate history
would invite someone to fabricate a row to silence it.
"""
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import check_register_freshness as guard  # noqa: E402

REGISTER = ROOT / 'CODEX_STATE_REGISTER.md'

ONE_HEADER = "**Version: v9** · _Last updated: 2026-08-08 · By: Mr. Code — pointer._"
LOG_HEAD = "## Changelog"


def _register(headers, versions):
    """Build a minimal register: given header versions and changelog versions."""
    out = ["# Codex State Register", ""]
    for v in headers:
        out.append(f"**Version: v{v}** · _Last updated: 2026-08-08 · By: Mr. Code "
                   f"(**RYA-{600 + v} — narrative.**)_")
    out += ["", LOG_HEAD]
    for v in versions:
        out.append(f"- **v{v}** (2026-08-08) — **RYA-{600 + v} — narrative.**")
    return "\n".join(out) + "\n"


def _run(tmp_path, text, monkeypatch):
    p = tmp_path / 'CODEX_STATE_REGISTER.md'
    p.write_text(text)
    monkeypatch.setattr(guard, 'ROOT', tmp_path)
    return guard.structure_mode()


# ── the real file must satisfy the invariant ─────────────────────────────────

def test_committed_register_passes_the_structural_guard():
    assert guard.structure_mode() == 0, (
        "CODEX_STATE_REGISTER.md violates the RYA-690 structural invariant")


def test_committed_register_has_exactly_one_header_line():
    n = sum(1 for l in REGISTER.read_text().splitlines()
            if l.startswith('**Version:'))
    assert n == 1, f"expected exactly one '**Version:' header line, found {n}"


# ── 1. header multiplicity ───────────────────────────────────────────────────

def test_two_header_lines_fail(tmp_path, monkeypatch, capsys):
    """The exact shape main was in: a merge that kept both sides."""
    rc = _run(tmp_path, _register([9, 8], [9, 8]), monkeypatch)
    assert rc == 1
    err = capsys.readouterr().err
    assert 'exactly ONE' in err and 'found 2' in err


def test_eight_header_lines_fail(tmp_path, monkeypatch):
    """Main really did reach eight."""
    assert _run(tmp_path, _register([9, 8, 7, 6, 5, 4, 3, 2],
                                    [9, 8, 7, 6, 5, 4, 3, 2]), monkeypatch) == 1


def test_zero_header_lines_fail(tmp_path, monkeypatch):
    assert _run(tmp_path, _register([], [9]), monkeypatch) == 1


# ── 2. changelog uniqueness ──────────────────────────────────────────────────

def test_duplicate_changelog_version_fails(tmp_path, monkeypatch, capsys):
    text = _register([9], [9, 8]) + "- **v8** (2026-08-08) — **RYA-999 — other.**\n"
    assert _run(tmp_path, text, monkeypatch) == 1
    assert 'more than' in capsys.readouterr().err


# ── 3. the orphan class ──────────────────────────────────────────────────────

def test_header_version_without_a_changelog_row_fails(tmp_path, monkeypatch, capsys):
    """RYA-682's exact failure: header claims v33, no v33 row anywhere."""
    rc = _run(tmp_path, _register([9], [8]), monkeypatch)
    assert rc == 1
    err = capsys.readouterr().err
    assert 'NO changelog row' in err or 'newest changelog row' in err


# ── 4. header must name the newest row ───────────────────────────────────────

def test_header_pointing_at_an_older_row_fails(tmp_path, monkeypatch, capsys):
    """RYA-673/RYA-675's failure: the version exists, but belongs to someone else.

    Assertion 3 alone passes here — the row exists — which is precisely how two of
    the three orphans slipped through. Assertion 4 is what catches it.
    """
    rc = _run(tmp_path, _register([8], [9, 8]), monkeypatch)
    assert rc == 1
    assert 'newest changelog row' in capsys.readouterr().err


def test_header_naming_the_newest_row_passes(tmp_path, monkeypatch):
    assert _run(tmp_path, _register([9], [9, 8, 7]), monkeypatch) == 0


# ── what must NOT be asserted ────────────────────────────────────────────────

def test_gaps_in_the_version_sequence_are_allowed(tmp_path, monkeypatch):
    """The changelog legitimately skips v25; demanding contiguity would fail on
    real history and tempt someone to fabricate a row to silence the guard."""
    assert _run(tmp_path, _register([9], [9, 7, 6, 3]), monkeypatch) == 0


def test_non_monotonic_changelog_order_is_allowed(tmp_path, monkeypatch):
    """The real file reads v36, v35, v34, v31, v32 — already out of order."""
    assert _run(tmp_path, _register([9], [9, 6, 7, 8]), monkeypatch) == 0


def test_real_register_is_indeed_non_contiguous_and_unordered():
    """Pin the two facts the omissions above rest on, so a future tightening of
    the guard trips here first rather than in CI on legitimate state."""
    import re
    vs = [int(m.group(1)) for l in REGISTER.read_text().splitlines()
          if (m := re.match(r'^-\s*\*\*v(\d+)\*\*\s*\(\d{4}-\d{2}-\d{2}\)', l))]
    assert len(vs) == len(set(vs)), "changelog versions must be unique"
    assert sorted(vs) != list(range(min(vs), max(vs) + 1)), (
        "the changelog is contiguous now — re-check whether a contiguity rule is safe")
    assert vs != sorted(vs, reverse=True), (
        "the changelog is fully ordered now — re-check whether an ordering rule is safe")


# ── the guard is wired into the CI entry point ───────────────────────────────

def test_structure_check_runs_in_the_default_mode():
    """It must not be reachable only via an opt-in flag."""
    src = (ROOT / 'scripts' / 'check_register_freshness.py').read_text()
    i_def = src.index('def main(')
    assert 'structure_mode()' in src[i_def:], (
        "structure_mode must be called from main(), not only from a sub-mode")


def test_structure_only_flag_works_standalone():
    r = subprocess.run([sys.executable, str(ROOT / 'scripts' / 'check_register_freshness.py'),
                        '--structure-only'], capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr


# ── §3D: the two files get deliberately different merge strategies ───────────

def test_sequence_md_has_union_merge_and_the_register_does_not():
    """'Accept both' is right for the log and wrong for the header — union on the
    register would have automated the bug rather than prevented it."""
    ga = (ROOT / '.gitattributes').read_text()
    seq = [l for l in ga.splitlines()
           if l.strip().startswith('SEQUENCE.md') and not l.strip().startswith('#')]
    assert seq and 'merge=union' in seq[0], "SEQUENCE.md must declare merge=union"
    reg = [l for l in ga.splitlines()
           if 'CODEX_STATE_REGISTER.md' in l and not l.strip().startswith('#')]
    assert not reg, (
        "CODEX_STATE_REGISTER.md must NOT get a merge strategy — union on a "
        "single-valued header automates the RYA-690 defect")


def test_header_carries_no_narrative():
    """§3D: prose lives in the changelog, so the header is one short token line.

    This is what shrinks the conflict surface: two branches still touch the line,
    but there is nothing to merge except the version itself.
    """
    hdr = next(l for l in REGISTER.read_text().splitlines()
               if l.startswith('**Version:'))
    assert len(hdr) < 400, (
        f"the header line is {len(hdr)} chars — narrative belongs in the changelog "
        f"row, not the header (RYA-690 §3D)")

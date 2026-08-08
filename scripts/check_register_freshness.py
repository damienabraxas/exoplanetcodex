#!/usr/bin/env python3
"""
Register-freshness guard (RYA-659).

Loud-fails if a STATE-CHANGING surface has been modified more recently than the
State Register -- i.e. the register drifted. This is the RYA-643 hole: element
state changed (data/audit/element_status_tracker.csv, 2026-07-27, 97485e9) while
CODEX_STATE_REGISTER.md sat at 2026-07-17 (08265f9), eleven tickets behind. The
register exists so "a gate cannot sign off while its rows are stale"; nothing was
checking that.

Two modes:

  --since-main : PR mode. Compare `git diff --name-only <base>...HEAD`. If any
                 state surface changed but CODEX_STATE_REGISTER.md did NOT change
                 in the same range -> FAIL.
  (default)    : history mode. If the newest commit touching any state surface is
                 newer than the newest commit touching the register -> FAIL.

Loud-fail, no silent pass: every git invocation raises on failure rather than
degrading to an empty result, because an empty result would look like "clean".

NOTE ON TIME: both modes read GIT history, never filesystem mtimes. `touch`-ing a
state surface therefore does NOT trip this guard, by design -- an uncommitted
mtime is not a state change. To exercise the teeth, commit a change to a state
surface without touching the register (see tests/test_register_freshness_rya659.py).

The surface registry is single-sourced in pipeline/state_surfaces.py -- do not
add a second copy here or in RYA-632's ledger_consistency_guard.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.state_surfaces import (  # noqa: E402
    REGISTER,
    STATE_SURFACES,
    existing_surfaces,
)


class GitError(RuntimeError):
    """A git invocation failed. Never swallowed -- a failed check is not a pass."""


def _run(*args: str) -> str:
    """Run git, raising loudly on failure. No silent fallback to ''."""
    try:
        out = subprocess.run(
            ["git", *args],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=True,
        )
    except FileNotFoundError as exc:  # git not installed
        raise GitError(f"git executable not found: {exc}") from exc
    except subprocess.CalledProcessError as exc:
        raise GitError(
            f"git {' '.join(args)} failed (exit {exc.returncode}): "
            f"{(exc.stderr or '').strip()}"
        ) from exc
    return out.stdout.strip()


def _last_commit_epoch(path: str) -> int:
    """Committer epoch of the newest commit touching `path`.

    Returns 0 only when git reports NO commit for the path (never committed).
    That is a real answer, not a swallowed error -- _run raises on git failure.
    """
    out = _run("log", "-1", "--format=%ct", "--", path)
    return int(out) if out else 0


def _fmt_gap(seconds: int) -> str:
    return f"{seconds // 86400}d {(seconds % 86400) // 3600}h"


def structure_mode() -> int:
    """RYA-690 structural guard — the register's header must be single-valued.

    RATIFIED by Ryan 2026-08-08 ("add the guard, three times is enough") after the
    same defect landed three times in a row: RYA-673, RYA-684 and RYA-682 each
    bumped the `**Version:` header without adding a changelog row. By the time
    RYA-690 ran, main carried EIGHT header lines, `v36` was claimed twice, and
    three landings had no changelog row of their own.

    It recurred because the natural merge resolution is wrong here. Every ticket
    edits the header at the same anchor, so concurrent branches always conflict;
    "accept both" is CORRECT for an append-only log like SEQUENCE.md and WRONG for
    a single-valued header. Same dialog, same gesture, opposite correctness. This
    guard is what tells the two apart.

    Three assertions, each derived from an observed failure:

      1. exactly one `**Version:` line          (eight were on main)
      2. its version has a changelog row        (kills the orphan class outright)
      3. no version appears twice in the log    (`v29` identified two things)

    Deliberately NOT asserted:

      * CONTIGUITY. The changelog legitimately skips v25, and a guard demanding an
        unbroken run would fail on real history and tempt someone to fabricate a
        row to silence it. Uniqueness is the invariant; completeness is not.
      * ORDERING. Verified against the file before deciding: the changelog reads
        v36, v35, v34, v31, v32, v30 ... — already non-monotonic. An ordering rule
        would fail on existing, legitimate state.

    Loud-fail per RYA-518: this raises the exit code, it never warns and continues.
    """
    path = ROOT / REGISTER
    if not path.exists():
        print(f"REGISTER STRUCTURE GUARD FAILED: {REGISTER} not found.", file=sys.stderr)
        return 1
    text = path.read_text(encoding="utf-8")
    lines = text.split("\n")

    headers = [
        (i, int(m.group(1)))
        for i, l in enumerate(lines, 1)
        if (m := re.match(r"^\*\*Version:\s*v(\d+)\*\*", l))
    ]
    log = [
        (i, int(m.group(1)))
        for i, l in enumerate(lines, 1)
        if (m := re.match(r"^-\s*\*\*v(\d+)\*\*\s*\(\d{4}-\d{2}-\d{2}\)", l))
    ]

    problems: list[str] = []

    if len(headers) != 1:
        where = ", ".join(f"line {i} (v{v})" for i, v in headers) or "none found"
        problems.append(
            f"expected exactly ONE '**Version:' header line, found {len(headers)}: {where}.\n"
            f"      A merge that kept both sides of this line is the usual cause. The header "
            f"is single-valued: keep the NEWER version and make sure the older one has its "
            f"own changelog row before dropping its line."
        )

    counts: dict[int, int] = {}
    for _, v in log:
        counts[v] = counts.get(v, 0) + 1
    dupes = sorted(v for v, c in counts.items() if c > 1)
    if dupes:
        problems.append(
            f"changelog version(s) {', '.join('v%d' % v for v in dupes)} appear more than "
            f"once. A version number identifies exactly one landing."
        )

    # 4. The header is DERIVED, not allocated: it must name the newest changelog row.
    #
    # This is the RYA-690 §3D structural fix made enforceable, and it closes the hole
    # the other three assertions leave. RYA-673's header said v29 and RYA-675's said
    # v36 — both versions EXIST in the changelog, but owned by RYA-581 and RYA-676
    # respectively, so an "is there a row?" check passes while the landing is still
    # orphaned. Requiring header == max(changelog) means the only way to bump the
    # header is to add a row, which is the behaviour the orphan class violated three
    # times. It also removes the independently-allocated counter that made two
    # concurrent branches collide on a number every single time.
    if log and len(headers) == 1:
        newest = max(v for _, v in log)
        hv = headers[0][1]
        if hv != newest:
            problems.append(
                f"header declares v{hv} but the newest changelog row is v{newest}. The "
                f"header is a POINTER derived from the changelog, not a number you "
                f"allocate: add your changelog row first, then point the header at it. "
                f"(If v{hv}'s row belongs to a different ticket, yours is orphaned.)"
            )

    log_versions = {v for _, v in log}
    for i, v in headers:
        if v not in log_versions:
            problems.append(
                f"header line {i} declares v{v}, which has NO changelog row. This is the "
                f"orphan class: the header is that landing's only record, so the line "
                f"cannot be dropped without losing it. Add a `- **v{v}** (YYYY-MM-DD) — ...` "
                f"row carrying the narrative, then collapse the header."
            )

    if problems:
        print("REGISTER STRUCTURE GUARD FAILED (RYA-690):", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1

    print(
        f"Register structure OK: 1 header line (v{headers[0][1]}), "
        f"{len(log)} changelog rows, all version numbers unique."
    )
    return 0


def history_mode() -> int:
    reg_epoch = _last_commit_epoch(REGISTER)
    if reg_epoch == 0:
        print(
            f"REGISTER FRESHNESS GUARD FAILED: no commit found for {REGISTER}. "
            "The register must be tracked in git -- its history IS its version record.",
            file=sys.stderr,
        )
        return 1

    surfaces = existing_surfaces(ROOT)
    stale = []
    for surf in surfaces:
        epoch = _last_commit_epoch(surf.path)
        if epoch > reg_epoch:
            stale.append((surf, epoch - reg_epoch))

    if stale:
        print(
            "REGISTER FRESHNESS GUARD FAILED (register is behind state surfaces):",
            file=sys.stderr,
        )
        for surf, gap in sorted(stale, key=lambda t: -t[1]):
            print(
                f"  - {surf.path} is newer than {REGISTER} by {_fmt_gap(gap)}"
                f"  [owns: {surf.owns}]",
                file=sys.stderr,
            )
        print(
            f"  Bump {REGISTER} (v-next + changelog) to reflect these, "
            "or cite why no state changed.",
            file=sys.stderr,
        )
        return 1

    skipped = len(STATE_SURFACES) - len(surfaces)
    note = f" ({skipped} not present on disk, skipped)" if skipped else ""
    print(
        f"Register freshness OK: {REGISTER} is at or ahead of all "
        f"{len(surfaces)} state surfaces{note}."
    )
    return 0


def since_main_mode(base: str) -> int:
    changed = set(_run("diff", "--name-only", f"{base}...HEAD").splitlines())
    touched = sorted(s for s in STATE_SURFACES if s.path in changed)

    if touched and REGISTER not in changed:
        print(
            "REGISTER FRESHNESS GUARD FAILED (PR changes state but not the register):",
            file=sys.stderr,
        )
        for surf in touched:
            print(
                f"  - {surf.path} changed in this PR  [owns: {surf.owns}]",
                file=sys.stderr,
            )
        print(
            f"  This PR must also update {REGISTER} "
            "(or document why state is unchanged).",
            file=sys.stderr,
        )
        return 1

    state = "updated" if REGISTER in changed else "not required"
    print(
        f"Register freshness OK (PR mode vs {base}): {len(touched)} state "
        f"surface(s) changed; register {state}."
    )
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument(
        "--since-main",
        metavar="BASE",
        nargs="?",
        const="origin/main",
        default=None,
        help="PR mode: check the diff BASE...HEAD (default base: origin/main).",
    )
    ap.add_argument(
        "--structure-only",
        action="store_true",
        help="run only the RYA-690 structural check (no git history needed).",
    )
    a = ap.parse_args()

    # RYA-690: the structural check runs in EVERY mode and FIRST. It reads the file,
    # not git history, so it costs nothing and cannot be skipped by choosing a mode.
    # Both failures are reported rather than short-circuiting, so one CI run tells
    # you everything that is wrong.
    structure = structure_mode()
    if a.structure_only:
        return structure

    try:
        freshness = since_main_mode(a.since_main) if a.since_main else history_mode()
    except GitError as exc:
        print(f"REGISTER FRESHNESS GUARD ERRORED (not a pass): {exc}", file=sys.stderr)
        return 2
    return structure or freshness


if __name__ == "__main__":
    sys.exit(main())

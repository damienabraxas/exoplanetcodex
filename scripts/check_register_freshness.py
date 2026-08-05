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
    a = ap.parse_args()
    try:
        return since_main_mode(a.since_main) if a.since_main else history_mode()
    except GitError as exc:
        print(f"REGISTER FRESHNESS GUARD ERRORED (not a pass): {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())

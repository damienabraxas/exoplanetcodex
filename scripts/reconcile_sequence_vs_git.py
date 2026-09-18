#!/usr/bin/env python3
"""
SEQUENCE-vs-git reconciler (RYA-1184, spec 2.1 / 2.4).

Git is merge truth; Linear pills are not. This walks `git log --merges
--first-parent` on a ref and builds the authoritative set of RYA tickets that
LANDED, then subtracts the set that SEQUENCE.md records as landings. The
remainder is drift, and the job is to drive it to EMPTY.

Two modes:

  (default)    : audit mode. Report the missing set for a range and exit 1 if it
                 is non-empty.
  --since-main : PR/hook mode. Same check against origin/main, used by the
                 pre-push hook so post-merge drift is caught by a machine
                 instead of by an audit two weeks later (RYA-1221 [K]).

NOTE ON "NAMED": a landing is recorded only by a LANDING LINE -- a list item
whose bullet opens with `- **RYA-XXX**`. An incidental mention inside another
ticket's summary does NOT count. RYA-1217 was the case that forced this: its
only occurrence in SEQUENCE.md was inside RYA-1134's sentence ("RYA-1217 handoff
remains held"), which records nothing about RYA-1217 having landed. A substring
search would have called that drift-free.

LOUD-FAIL, NEVER WARN-AND-PASS (RYA-1184 permanent rule): every git invocation
raises on failure rather than degrading to an empty result, because an empty
result is indistinguishable from "clean" and would silently pass the gate.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEQUENCE = ROOT / "SEQUENCE.md"

TICKET_RE = re.compile(r"RYA-(\d+)")
# A landing line: "- **RYA-1185** — ...", "- **RYA-587 / RYA-1213** — ...",
# "- **RYA-1143/1144/1150** — ...", "- **RYA-1185** (Part A) — ..."
LANDING_LINE_RE = re.compile(r"^\s*-\s+\*\*(?P<head>[^*]+)\*\*")
# SEQUENCE.md declares its own high-water mark: "> **Current as of `main` f3688d5**".
# That stamp IS the reconcile baseline -- advancing it is part of doing the backfill,
# so the gate can actually go green instead of re-litigating all history every push.
STAMP_RE = re.compile(r"Current as of\s+`main`\s+(?P<sha>[0-9a-f]{7,40})", re.I)


def _git(*args: str) -> str:
    """Run git, raising on any failure -- an empty result must never look clean."""
    proc = subprocess.run(
        ["git", "-C", str(ROOT), *args],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"git {' '.join(args)} failed ({proc.returncode}): {proc.stderr.strip()}"
        )
    return proc.stdout


def declared_baseline() -> str | None:
    """The SHA SEQUENCE.md claims to be current as of, or None if unstamped."""
    if not SEQUENCE.exists():
        return None
    m = STAMP_RE.search(SEQUENCE.read_text(encoding="utf-8"))
    return m.group("sha") if m else None


def landed_tickets(ref: str, since: str | None) -> dict[str, list[tuple[str, str, str]]]:
    """{ticket: [(sha, date, subject)]} for every first-parent merge on `ref`.

    A ticket is attributed from the merge subject, which carries the source
    branch ("Merge pull request #552 from owner/codex/rya-1222-orchestrator-mvp").
    Branch names are the reliable carrier; commit bodies cite unrelated tickets.
    """
    fmt = "%H\x1f%ad\x1f%s"
    args = ["log", "--merges", "--first-parent", f"--format={fmt}", "--date=short"]
    if since:
        args.append(f"{since}..{ref}")
    else:
        args.append(ref)

    out: dict[str, list[tuple[str, str, str]]] = {}
    for line in _git(*args).splitlines():
        if not line.strip():
            continue
        sha, date, subject = line.split("\x1f", 2)
        # attribute from the branch portion only (after "from "), not the whole
        # subject -- a PR title can name tickets it did not land.
        branch = subject.split(" from ", 1)[1] if " from " in subject else subject
        for num in TICKET_RE.findall(branch.upper().replace("RYA-", "RYA-")):
            out.setdefault(f"RYA-{num}", []).append((sha[:8], date, subject))
        for num in re.findall(r"rya[-_]?(\d+)", branch, flags=re.I):
            out.setdefault(f"RYA-{num}", []).append((sha[:8], date, subject))
    # de-duplicate the two passes
    for k in out:
        seen, uniq = set(), []
        for rec in out[k]:
            if rec[0] not in seen:
                seen.add(rec[0])
                uniq.append(rec)
        out[k] = uniq
    return out


def recorded_tickets() -> set[str]:
    """Tickets SEQUENCE.md records as LANDINGS (landing lines only)."""
    if not SEQUENCE.exists():
        raise RuntimeError(f"{SEQUENCE} is missing -- cannot reconcile against nothing")
    found: set[str] = set()
    for raw in SEQUENCE.read_text(encoding="utf-8").splitlines():
        m = LANDING_LINE_RE.match(raw)
        if not m:
            continue
        for num in TICKET_RE.findall(m.group("head")):
            found.add(f"RYA-{num}")
        # "- **RYA-1143/1144/1150**" -- bare numbers after the first ticket
        head = m.group("head")
        if "RYA-" in head:
            tail = head.split("RYA-", 1)[1]
            for num in re.findall(r"(?<![\d-])(\d{2,5})(?![\d])", tail):
                found.add(f"RYA-{num}")
    return found


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--since", default=None,
                    help="SHA/ref/date to reconcile from (exclusive). Default: all history.")
    ap.add_argument("--ref", default="origin/main", help="Ref to treat as merge truth.")
    ap.add_argument("--since-main", action="store_true",
                    help="Hook/PR mode: reconcile origin/main, same teeth.")
    args = ap.parse_args()

    ref = "origin/main" if args.since_main else args.ref

    since = args.since
    baseline_src = "--since"
    if since is None:
        since = declared_baseline()
        baseline_src = "SEQUENCE.md 'Current as of' stamp"
    if since is None:
        baseline_src = "(none -- SEQUENCE.md carries no stamp; reconciling ALL history)"

    try:
        landed = landed_tickets(ref, since)
        recorded = recorded_tickets()
    except RuntimeError as exc:
        print(f"RECONCILER FAILED (not a clean result): {exc}", file=sys.stderr)
        return 2

    missing = {t: recs for t, recs in landed.items() if t not in recorded}

    print(f"SEQUENCE-vs-git reconciler (RYA-1184)")
    print(f"  ref             : {ref} @ {_git('rev-parse', '--short', ref).strip()}")
    print(f"  since           : {since or '(all history)'}  [{baseline_src}]")
    print(f"  landed tickets  : {len(landed)}")
    print(f"  recorded in SEQUENCE.md (landing lines): {len(recorded)}")
    print(f"  MISSING SET     : {len(missing)}")

    if not missing:
        print("\nOK: every ticket landed on this ref is named by a landing line in SEQUENCE.md.")
        return 0

    print("\nDRIFT -- these landed on the ref and SEQUENCE.md does not record them:")
    for ticket in sorted(missing, key=lambda t: int(t.split("-")[1])):
        for sha, date, subject in missing[ticket]:
            print(f"  {ticket:<10} {sha}  {date}  {subject}")
    print("\nFAIL: bump SEQUENCE.md (and CODEX_STATE_REGISTER.md) in the same PR (RYA-659),")
    print("      then advance SEQUENCE.md's 'Current as of `main` <sha>' stamp to the ref above.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

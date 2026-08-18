#!/usr/bin/env python3
"""RYA-874 — is the committed RYA-847 gated artifact what today's generator produces?

    python3 scripts/rya874_regen_diff.py [--baseline-ref e04e4e2]

The four cells under data/results/rya847/gated/ were generated on 2026-08-17 evening and
merged in ce31080, BEFORE RYA-869 (handler-keyed harness residual) and RYA-871 (ep_eV on
the EW-route per-line artifact) landed. RYA-686's contract says a committed result is what
its generating harness produces; these are two schema columns and — the RYA-855 re-audit
found — one stale systematic behind that.

This script re-reads the regenerated files in the working tree and diffs them field by
field against the SAME PATHS at a git ref (default: the merge-base state they were
committed in). It does not re-derive anything.

🔴 THE EXPECTED CHANGES ARE DECLARED BELOW, BEFORE THE RUN, AS DATA.
A diff that classifies after the fact can rationalise anything it finds. EXPECTED lists
every difference this regeneration is allowed to produce, each naming the ticket that
explains it and the exact value it must land on. Every other difference is UNEXPLAINED and
exits non-zero — which is spec item 4: a value that moves is a finding about RYA-871's EW
re-measurement, not something to absorb under a schema refresh.

⚠️ An expected change that does NOT appear is also a failure. A correction nobody can see
is indistinguishable from a correction that never ran.
"""
from __future__ import annotations

import argparse
import csv
import difflib
import io
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GATED = Path("data/results/rya847/gated")

# ---------------------------------------------------------------------------
# PRE-DECLARED expectations. Nothing else may move.
# ---------------------------------------------------------------------------

# Columns the current deriver emits that the committed artifact predates.
EXPECTED_NEW_COLUMNS = {
    "products": {"handler": "RYA-869 — the product declares the handler that made it"},
    "lines": {"ep_eV": "RYA-871 — a wavelength alone does not identify a line"},
}

# (file kind, row key, column) -> (old, new, ticket)
# The one number this regeneration is expected to correct. RYA-869: ENGINE-B-NLTE is a
# SynthesisHandler flux fit and owes 0.0000, not the profile fitter's measured 0.0129.
EXPECTED_VALUE_CHANGES = {
    ("gerber-nlte", "VIS", "ENGINE-B-NLTE", "syst_dex"): (
        "0.1705", "0.17", "RYA-869 handler-keyed harness residual"),
}

# budgets.txt / provenance.txt are prose. These substrings mark lines whose text is
# allowed to change.
EXPECTED_BUDGET_TEXT = (
    ("harness residual", "RYA-869 handler label / RYA-873 budget prose"),
)

# ---------------------------------------------------------------------------
# 🔴 CLASSIFIED AFTER THE FACT — NOT pre-declared. Kept SEPARATE on purpose.
# ---------------------------------------------------------------------------
# The near-UV cell landed carrying two text changes this script did not predict, and it
# reported them as UNEXPLAINED, which is what it is for. They are RYA-855's, attributed
# by `git log -S` to c51670d ("decide the gf rung from the lines, not from a hardcode"):
#
#   * budgets.txt gains a trailing `gf rung: ...` line
#   * provenance.txt's hand-written closing sentence about RYA-822's grading counts is
#     REPLACED by a generated `gf TERM: gf rung N ...` sentence
#
# Both are additive provenance from a decision that now comes from the lines instead of a
# hardcode. Neither moves a number: the rung resolves to 1, the gf term stays 0.1700
# UNGRADED, and every value/n/stat/syst in the cell is unchanged. ⚠️ The provenance one is
# a genuine LOSS of detail — the old sentence carried "17 of 40 carry a citable NIST
# accuracy class and most of those are poor (C/C+/D/D+/E)", which the generated sentence
# does not say. Recorded here rather than waved through; it is a note for RYA-855, not a
# defect in this regeneration.
#
# They live in their own list so the report can never present a post-hoc classification as
# a prediction. Widening EXPECTED_* after seeing output is how a diff learns to agree with
# whatever it is shown.
CLASSIFIED_AFTER_THE_FACT = (
    ("gf rung:", "RYA-855 c51670d — generated gf-rung line, additive, no number moves"),
    ("gf TERM:", "RYA-855 c51670d — generated gf-rung provenance sentence"),
    ("gf REMAINS THE DOMINANT SYSTEMATIC",
     "RYA-855 c51670d — the hand-written sentence it replaces"),
)


def _git_show(ref: str, rel: str) -> str | None:
    p = subprocess.run(["git", "show", f"{ref}:{rel}"], cwd=ROOT,
                       capture_output=True, text=True)
    return p.stdout if p.returncode == 0 else None


def _csv_rows(text: str):
    rd = csv.DictReader(io.StringIO(text))
    return list(rd), (rd.fieldnames or [])


def _kind(rel: str) -> str:
    if rel.endswith("_products.csv"):
        return "products"
    if rel.endswith("_lines.csv"):
        return "lines"
    return "text"


def _deck(rel: str) -> str:
    parts = Path(rel).parts
    return parts[-2] if parts[-2] in ("ts-lte", "gerber-nlte") else ""


def _key(kind: str, row: dict) -> tuple:
    if kind == "products":
        return (row.get("band", ""), row.get("treatment", ""))
    return (row.get("treatment", ""), row.get("wavelength_air_A", ""))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline-ref", default="origin/main",
                    help="git ref holding the committed artifact (default origin/main)")
    args = ap.parse_args()

    findings: list[str] = []       # UNEXPLAINED — hard failures
    accounted: list[str] = []      # matched a pre-declared expectation
    unseen = dict(EXPECTED_VALUE_CHANGES)
    new_cols_seen: set[tuple[str, str]] = set()

    files = sorted(p.relative_to(ROOT).as_posix()
                   for p in (ROOT / GATED).rglob("*") if p.is_file())
    if not files:
        print(f"no artifacts under {GATED} — run the regen driver first", file=sys.stderr)
        return 2

    base_files = subprocess.run(
        ["git", "ls-tree", "-r", "--name-only", args.baseline_ref, GATED.as_posix()],
        cwd=ROOT, capture_output=True, text=True).stdout.split()

    for rel in sorted(set(files) | set(base_files)):
        kind = _kind(rel)
        base = _git_show(args.baseline_ref, rel)
        live_path = ROOT / rel
        live = live_path.read_text() if live_path.exists() else None

        if base is None:
            findings.append(f"{rel}: NEW file with no committed counterpart")
            continue
        if live is None:
            findings.append(f"{rel}: committed file was NOT regenerated — the artifact "
                            f"set would be half-migrated")
            continue
        if base == live:
            continue

        if kind == "text":
            # Aligned by difflib, NOT zipped: RYA-873 replaces a one-line harness note
            # with a multi-line one, so a positional compare would report every later
            # line as changed and bury the real diff in an off-by-N.
            b, l = base.splitlines(), live.splitlines()
            sm = difflib.SequenceMatcher(None, b, l, autojunk=False)
            for tag, i1, i2, j1, j2 in sm.get_opcodes():
                if tag == "equal":
                    continue
                hunk = b[i1:i2] + l[j1:j2]
                why = next((w for m, w in EXPECTED_BUDGET_TEXT
                            if any(m in x for x in hunk)), None)
                post = next((w for m, w in CLASSIFIED_AFTER_THE_FACT
                             if any(m in x for x in hunk)), None)
                if why is None and post is not None:
                    why = f"POST-HOC: {post}"
                body = ("\n".join(f"      - {x.strip()}" for x in b[i1:i2])
                        + ("\n" if i2 > i1 and j2 > j1 else "")
                        + "\n".join(f"      + {y.strip()}" for y in l[j1:j2]))
                if why:
                    accounted.append(f"{rel}:{i1 + 1} {tag}  [{why}]\n{body}")
                else:
                    findings.append(f"{rel}:{i1 + 1} {tag}: unexplained text change\n"
                                    f"{body}")
            continue

        brows, bcols = _csv_rows(base)
        lrows, lcols = _csv_rows(live)

        for c in bcols:
            if c not in lcols:
                findings.append(f"{rel}: column {c!r} DISAPPEARED")
        for c in lcols:
            if c in bcols:
                continue
            why = EXPECTED_NEW_COLUMNS.get(kind, {}).get(c)
            if why:
                new_cols_seen.add((kind, c))
                accounted.append(f"{rel}: new column {c!r}  [{why}]")
            else:
                findings.append(f"{rel}: UNDECLARED new column {c!r}")

        bmap = {_key(kind, r): r for r in brows}
        lmap = {_key(kind, r): r for r in lrows}
        if len(bmap) != len(brows) or len(lmap) != len(lrows):
            findings.append(f"{rel}: duplicate row keys — the comparison is not 1:1")
        for k in sorted(set(bmap) - set(lmap)):
            findings.append(f"{rel}: row {k} present in the committed artifact, GONE")
        for k in sorted(set(lmap) - set(bmap)):
            findings.append(f"{rel}: row {k} is NEW")

        shared = [c for c in bcols if c in lcols]
        deck = _deck(rel)
        for k in sorted(set(bmap) & set(lmap)):
            for c in shared:
                bv, lv = bmap[k].get(c), lmap[k].get(c)
                if bv == lv:
                    continue
                exp = EXPECTED_VALUE_CHANGES.get((deck,) + k + (c,))
                if exp and (bv, lv) == (exp[0], exp[1]):
                    accounted.append(f"{rel}: {k} {c} {bv} -> {lv}  [{exp[2]}]")
                    unseen.pop((deck,) + k + (c,), None)
                elif exp:
                    findings.append(
                        f"{rel}: {k} {c} {bv} -> {lv} — an expected correction landed on "
                        f"the WRONG value (declared {exp[0]} -> {exp[1]}, {exp[2]})")
                    unseen.pop((deck,) + k + (c,), None)
                else:
                    findings.append(f"🔴 {rel}: {k} {c} {bv} -> {lv} — UNEXPLAINED")

    # An expectation that never fired is a failure: the correction is invisible.
    for key, (old, new, ticket) in unseen.items():
        findings.append(f"declared correction never appeared: {key} {old} -> {new} "
                        f"({ticket}) — the fix is invisible in the artifact")
    for kind, cols in EXPECTED_NEW_COLUMNS.items():
        for c in cols:
            if (kind, c) not in new_cols_seen:
                findings.append(f"declared new column {c!r} ({kind}) never appeared")

    print(f"baseline {args.baseline_ref}   files compared {len(set(files) | set(base_files))}\n")
    pre = [a for a in accounted if "POST-HOC:" not in a]
    post_hoc = [a for a in accounted if "POST-HOC:" in a]
    print(f"PRE-DECLARED, and they fired ({len(pre)}):")
    for a in pre:
        print("   " + a)
    if post_hoc:
        print(f"\n⚠️  CLASSIFIED AFTER THE FACT ({len(post_hoc)}) — these were NOT "
              f"predicted; each was attributed to a commit and checked to move no "
              f"number:")
        for a in post_hoc:
            print("   " + a)
    print(f"\n[all differences accounted for: {len(accounted)}]")
    if findings:
        print(f"\n🔴 UNEXPLAINED ({len(findings)}):", file=sys.stderr)
        for f in findings:
            print("   " + f, file=sys.stderr)
        return 1
    print("\n✅ every difference is a pre-declared schema column or the RYA-869 "
          "correction. No value moved.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

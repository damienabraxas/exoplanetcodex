#!/usr/bin/env python
"""
Guard: a result artifact must not land without its generating harness (RYA-686).

THE FAILURE THIS PREVENTS
-------------------------
Sirius computes; only the result JSON comes back; the harness that produced it stays
in a scratch directory and is never committed. Two verified instances:

  * RYA-559 shipped `data/results/solar_ba_synthesis_rya559.json` (merge 564824a) plus a
    phase_c hook and a test -- but NOT the synthesis harness. `scripts/rya559_ba2_synth_sirius.py`
    has never existed in this repo on any branch, and is not on Sirius either. The
    measurement that put Ba on the board at 2.410 is NOT REPRODUCIBLE, and RYA-581 had
    to rebuild the harness from scratch off the RYA-551 pattern plus a commit message.
  * `scripts/rya560_zr2_synth_sirius.py` IS committed but declares no argparse at all,
    so the `--deblend` entrypoint RYA-585's brief called for did not exist.

Per RYA-567 the compute-on-Sirius rule is correct; what was missing is the preservation
half. This guard is that half.

THE LINKAGE MECHANISM: a manifest (`data/results/GENERATORS.yaml`)
-----------------------------------------------------------------
Three mechanisms were considered. The manifest wins on four counts that the other two
cannot cover -- the reasoning is written up in docs/CONVENTIONS.md and summarised here
so it travels with the code:

  1. FORMAT COVERAGE. 14 of the 32 tracked artifacts are `.csv`/`.txt`. They cannot carry
     a `_meta.generator` JSON field, so that route needs a second mechanism for them.
  2. NO ARTIFACT MUTATION. A `_meta` route means editing historical result files --
     including bare wavelength-keyed dicts like sr2_synthesis_rya551.json, where a new
     top-level key changes what every `for w, d in obj.items()` consumer iterates over.
     A manifest touches no artifact.
  3. HONEST NEGATIVES. The manifest can say `status: UNREPRODUCIBLE, generator: null`
     with a ticket and a reason. A naming convention has no way to express that, so it
     needs an exemption file anyway -- i.e. a manifest, arrived at by a longer road.
  4. THE INVOCATION. Only the manifest carries HOW to run the harness, which is what
     makes the RYA-560 class partially catchable (see `--check-invocations`).

The standing objection to a manifest is "another thing to keep in sync". That is answered
by making the guard BIDIRECTIONAL: an artifact with no entry fails, and an entry naming a
missing artifact or a missing generator fails. Drift is a build break, not a silent
divergence -- the same discipline RYA-654 applies to the element status tracker.

A naming convention (`*_ryaNNN.*` -> `scripts/*ryaNNN*.py`) was rejected as ALSO PROVING
THE WRONG THING: RYA-485 shipped two artifacts and two scripts, so a ticket-token match
cannot tell which script generates which artifact, and would report green if either one
went missing.

USAGE
-----
    python scripts/check_result_generators.py                     # full guard
    python scripts/check_result_generators.py --check-invocations # + argparse flag check
    python scripts/check_result_generators.py --root /path/to/repo

Exit 0 = clean, exit 1 = one or more violations (each printed as one line).
"""

from __future__ import annotations

import argparse
import ast
import shlex
import subprocess
import sys
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover - yaml is a hard requirement of the repo
    print("::error::PyYAML is required by scripts/check_result_generators.py")
    raise

#: Directory the convention governs, repo-relative and POSIX-style.
RESULTS_DIR = "data/results"

#: Every directory whose TRACKED contents this guard covers.
#:
#: RYA-939 added data/processed/. The convention was scoped to data/results/ while the
#: single most consequential product in the repo -- the normalised HARPS solar spectrum
#: that every committed EW derives from -- sat in data/processed/, gitignored, tracked
#: nowhere and guarded by nothing. Deliberately committing it (RYA-939) without
#: extending the guard would have recreated exactly the RYA-559 hole this file exists to
#: close: an artifact on the board that nobody can re-run.
SCANNED_DIRS = (RESULTS_DIR, "data/processed")

#: The manifest itself, which lives beside the artifacts it describes (the same shape as
#: ARTIFACT_MANIFEST.csv at the RYA-461 store root). Excluded from its own coverage check.
MANIFEST = f"{RESULTS_DIR}/GENERATORS.yaml"

#: Files under data/results/ that are not result artifacts and need no generator.
NON_ARTIFACTS = frozenset({MANIFEST, "data/processed/.gitkeep"})

#: Allowed values of an entry's `status`. Three, because the RYA-686 audit found three
#: genuinely different things sitting in data/results/ and collapsing them would lie:
#:
#:   COMMITTED      machine output, harness in the repo. The normal case and the default.
#:   HAND_AUTHORED  a human wrote this file -- a decision record or an analysis summary,
#:                  not a program's stdout. There is no harness to commit, so it must
#:                  instead name the `sources:` its numbers came from.
#:   UNREPRODUCIBLE machine output whose harness was NEVER committed. The honest record
#:                  of a defect, not an exemption to be reached for. Adding one is a
#:                  visible diff and the frozen-set test (RYA-686) makes it deliberate.
STATUS_COMMITTED = "COMMITTED"
STATUS_HAND_AUTHORED = "HAND_AUTHORED"
STATUS_UNREPRODUCIBLE = "UNREPRODUCIBLE"
VALID_STATUSES = (STATUS_COMMITTED, STATUS_HAND_AUTHORED, STATUS_UNREPRODUCIBLE)

#: Statuses that must name an existing generator file.
STATUS_REQUIRES_GENERATOR = (STATUS_COMMITTED,)

#: Statuses that must NOT name a generator -- claiming one would be the fabrication the
#: convention exists to prevent.
STATUS_FORBIDS_GENERATOR = (STATUS_HAND_AUTHORED, STATUS_UNREPRODUCIBLE)


def repo_root(explicit: str | None = None) -> Path:
    if explicit:
        return Path(explicit).resolve()
    return Path(__file__).resolve().parents[1]


def tracked_artifacts(root: Path) -> list[str]:
    """Repo-relative POSIX paths of every TRACKED file under data/results/.

    Tracked, not on-disk: the convention is about what LANDS. A local Sirius-run output
    sitting untracked in a working tree is not a violation of anything -- committing it
    without its harness is.
    """
    out = subprocess.run(
        ["git", "-C", str(root), "ls-files", "--", *SCANNED_DIRS],
        capture_output=True, text=True, check=True,
    ).stdout
    return sorted(p for p in out.splitlines() if p and p not in NON_ARTIFACTS)


def load_manifest(root: Path) -> list[dict]:
    path = root / MANIFEST
    if not path.exists():
        raise FileNotFoundError(f"result-generator manifest missing: {MANIFEST}")
    doc = yaml.safe_load(path.read_text()) or {}
    entries = doc.get("artifacts")
    if entries is None:
        raise ValueError(f"{MANIFEST} has no top-level `artifacts:` list")
    return list(entries)


def declared_long_flags(script: Path) -> tuple[set[str], bool]:
    """(long option strings declared via argparse, whether ANY add_argument was found).

    Static AST read of string literals passed to `*.add_argument(...)`. Deliberately
    conservative: it reports what the source literally declares and never imports or
    executes the harness.
    """
    flags: set[str] = set()
    found = False
    tree = ast.parse(script.read_text())
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "add_argument"):
            continue
        found = True
        for arg in node.args:
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str) \
                    and arg.value.startswith("--"):
                flags.add(arg.value)
    return flags, found


def invocation_long_flags(invocation: str) -> set[str]:
    """Long flags a recorded invocation string actually uses (`--x=1` -> `--x`)."""
    return {tok.split("=", 1)[0] for tok in shlex.split(invocation)
            if tok.startswith("--") and tok != "--"}


def check(root: Path, check_invocations: bool = False) -> list[str]:
    """Return one human-readable violation string per problem. Empty list = clean."""
    problems: list[str] = []
    entries = load_manifest(root)

    by_artifact: dict[str, dict] = {}
    for i, e in enumerate(entries):
        if not isinstance(e, dict):
            problems.append(f"{MANIFEST}: entry #{i} is not a mapping")
            continue
        name = e.get("artifact")
        if not name:
            problems.append(f"{MANIFEST}: entry #{i} has no `artifact:` key")
            continue
        if name in by_artifact:
            problems.append(f"{MANIFEST}: duplicate entry for artifact '{name}'")
            continue
        by_artifact[name] = e

    tracked = tracked_artifacts(root)
    # An entry's `artifact:` is relative to data/results/ (the original convention) OR
    # repo-relative when it names one of the other scanned dirs. Each tracked file gets
    # ONE canonical key -- the form its entry is expected to use -- plus the alternates
    # accepted on lookup. Registering both forms as separate keys would make every
    # data/results artifact look unaccounted-for under its full path.
    def _keys(path: str) -> list[str]:
        if path.startswith(RESULTS_DIR + "/"):
            return [path[len(RESULTS_DIR) + 1:], path]
        return [path]

    tracked_rel = {k: p for p in tracked for k in _keys(p)}

    # (1) THE CORE GATE -- every landed artifact is accounted for.
    for full in tracked:
        rel = _keys(full)[0]
        if not any(k in by_artifact for k in _keys(full)):
            problems.append(
                f"{full}: result artifact has NO entry in {MANIFEST}. A result must not "
                f"land without its generating harness (RYA-686). Commit the harness and "
                f"add an entry; if the measurement genuinely cannot be reproduced, say so "
                f"with status: {STATUS_UNREPRODUCIBLE} rather than inventing a generator."
            )

    # (2) The reverse direction -- no entry may describe an artifact that is not there.
    for rel in sorted(by_artifact):
        if rel not in tracked_rel:
            problems.append(
                f"{MANIFEST}: entry '{rel}' names an artifact that is not tracked under "
                f"{'/, '.join(SCANNED_DIRS)}/ (renamed or deleted?)"
            )

    # (3) Per-entry integrity.
    for rel in sorted(by_artifact):
        e = by_artifact[rel]
        where = f"{MANIFEST}: entry '{rel}'"
        status = e.get("status", STATUS_COMMITTED)
        if status not in VALID_STATUSES:
            problems.append(f"{where}: status '{status}' not one of {VALID_STATUSES}")
            continue
        if not e.get("ticket"):
            problems.append(f"{where}: `ticket:` is mandatory (no entry without a ticket)")

        generator = e.get("generator")

        if status in STATUS_FORBIDS_GENERATOR:
            if generator:
                problems.append(
                    f"{where}: status {status} must carry `generator: null`. "
                    f"If a generator exists the status is {STATUS_COMMITTED}."
                )
            if not e.get("note"):
                problems.append(
                    f"{where}: status {status} requires a `note:` recording WHY there is no "
                    f"generator. An unexplained exemption is a hiding place for the next "
                    f"missing harness."
                )
            if status == STATUS_HAND_AUTHORED and not e.get("sources"):
                problems.append(
                    f"{where}: status {STATUS_HAND_AUTHORED} requires a non-empty `sources:` "
                    f"list. A human-written results file still has to say where its numbers "
                    f"came from, or it is an unsourced claim."
                )
            continue

        if status in STATUS_REQUIRES_GENERATOR:
            if not generator:
                problems.append(f"{where}: `generator:` is missing (status {status})")
                continue
            if not (root / generator).exists():
                problems.append(
                    f"{where}: generator '{generator}' does not exist in the repo. This is "
                    f"the RYA-559 failure mode -- the artifact landed, the harness did not."
                )
                continue

        # (4) The RYA-560 half -- is the harness invocable AS DOCUMENTED?
        invocation = e.get("invocation")
        if check_invocations and invocation and e.get("invocation_check") != "skip":
            wanted = invocation_long_flags(invocation)
            if not wanted:
                continue
            gen_path = root / generator
            if gen_path.suffix != ".py":
                continue
            declared, has_argparse = declared_long_flags(gen_path)
            if not has_argparse:
                problems.append(
                    f"{where}: recorded invocation uses {sorted(wanted)} but "
                    f"'{generator}' declares no argparse arguments at all -- it is committed "
                    f"but not invocable as documented (the RYA-560 case)."
                )
                continue
            missing = sorted(wanted - declared)
            if missing:
                problems.append(
                    f"{where}: recorded invocation uses {missing}, which '{generator}' "
                    f"does not declare via add_argument()."
                )

    return problems


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--root", default=None,
                    help="repo root to check (default: the repo this script lives in)")
    ap.add_argument("--check-invocations", action="store_true",
                    help="also verify each recorded invocation's long flags are declared "
                         "by the generator's argparse (catches the RYA-560 case)")
    args = ap.parse_args(argv)

    root = repo_root(args.root)
    problems = check(root, check_invocations=args.check_invocations)

    if problems:
        for p in problems:
            print(f"::error::{p}")
        print(f"\nFAIL: {len(problems)} result-generator violation(s) "
              f"(RYA-686; see docs/CONVENTIONS.md).")
        return 1

    n = len(tracked_artifacts(root))
    print(f"OK: all {n} tracked artifacts under {', '.join(d + '/' for d in SCANNED_DIRS)} "
          f"are accounted for in {MANIFEST}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

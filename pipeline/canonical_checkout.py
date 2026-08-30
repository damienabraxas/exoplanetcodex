"""RYA-1138 — a checkout must be able to PROVE it can produce a valid baseline.

RYA-1090 fixed one worktree-relative resolver: `store_root()` derived the artifact
store from `Path(__file__).resolve().parents[2]`, the parent of whatever worktree the
module happened to be imported from. Move the worktrees and the store forked in two,
silently, and 85 artifacts became unreachable.

That fix was store-specific. The MECHANISM was not, and it is still live in every
other resolver that reads through a worktree's parent. Three instances fired in one
session (2026-08-30), each of which INVERTS or NULLIFIES a verdict rather than
producing a visible error:

  1. FALSE ABSENCE   — a search of the committed tree plus Sirius reported per-line
     artifacts "do not exist". They were in `~/codex/rya1015`, the uncommitted
     worktree THE FEED'S OWN PROVENANCE NAMES.
  2. FALSE MUTATION  — `git checkout --` used to undo a mutation reverted to HEAD and
     wiped the uncommitted fix under test. The failures were real but were not
     measuring the intended mutation.
  3. FALSE BASELINE  — a baseline worktree sited outside `~/codex` produced 25 phantom
     collection errors, because `ISPEC_DIR` is `ROOT.parent / 'ispec'`. Diffed against
     that baseline, ELEVEN REAL FAILURES WOULD HAVE READ AS AN IMPROVEMENT.

🔴 WHY A LOCATION CHECK IS NOT ENOUGH, AND WHAT THIS MODULE CHECKS INSTEAD. The
obvious guard -- "is the baseline under `~/codex`?" -- encodes today's accident. It
passes a checkout that sits in the right place but cannot import `ispec`, and it fails
a perfectly good checkout the day the canonical parent moves, which is exactly the
event RYA-1090 was created by. Location is a proxy.

What actually went wrong is that the two sides of the diff HAD DIFFERENT CAPABILITIES:
one could import `ispec` and collect the full suite, the other could not. So the check
here is a CAPABILITY FINGERPRINT, and the rule is not "be canonical" but:

    THE TWO SIDES OF A SET-DIFF MUST HAVE IDENTICAL CAPABILITY FINGERPRINTS.

That is falsifiable, it degrades correctly when the layout changes, and it catches the
whole family rather than the one member we happened to trip over -- including a
baseline that silently SKIPS tests the branch runs, which no location check can see
(a skipped test is in neither the pass set nor the fail set).

Location is still reported, because a non-canonical checkout is a good early WARNING,
but it is advisory. A fingerprint mismatch is what refuses.
"""
from __future__ import annotations

import json
import re
import os
import subprocess
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path

#: RYA-1138 — the canonical worktree parent, NAMED rather than derived, for the same
#: reason RYA-1090 named `CANONICAL_STORE_ROOT`: a location that other things resolve
#: THROUGH cannot itself be inferred from ambient context without forking the moment
#: that context moves. Advisory here -- it drives a warning, never a refusal.
CANONICAL_WORKTREE_PARENT = Path(
    os.environ.get("CODEX_WORKTREE_PARENT") or (Path.home() / "codex")
).expanduser()

#: The external dependencies whose ABSENCE is silent -- each one turns tests into
#: collection errors or skips rather than failures, which is what makes a broken
#: baseline read as a clean one. Extend this when a new silent dependency appears.
#:
#: ⚠️ `ispec` is the one that caused instance 3. It is resolved by
#: `config.constants.ISPEC_DIR`, which defaults to `ROOT.parent / 'ispec'` -- above
#: the repo root, and therefore worktree-relative in exactly RYA-1090's sense.
SILENT_DEPENDENCIES = ("ispec", "kp_atlas")


class NonComparableCheckoutsError(RuntimeError):
    """The two sides of a set-diff cannot be compared (RYA-1138).

    Raised INSTEAD of emitting a verdict. A verdict from mismatched checkouts is worse
    than no verdict, because it is trusted: RYA-908's eleven regressions would have
    been reported as three fewer failures than main.
    """


@dataclass
class Capabilities:
    """What a checkout can actually DO, as opposed to where it sits."""

    root: str
    is_canonical_location: bool
    ispec_importable: bool
    ispec_dir: str
    kp_atlas_present: bool
    kp_atlas_dir: str
    python_version: str
    collected_tests: int | None = None
    collection_errors: int | None = None
    notes: list[str] = field(default_factory=list)

    def fingerprint(self) -> dict:
        """The subset that MUST match between the two sides of a diff.

        🔴 `is_canonical_location` is deliberately NOT in here. Two checkouts that are
        both non-canonical in the same way are still comparable to each other, and a
        canonical checkout that cannot import `ispec` is not comparable to one that
        can. Comparability is about capability, not address.

        `collected_tests` is also excluded: a branch that ADDS tests legitimately
        collects more, and folding that in would make every test-adding branch
        unmergeable. Collection ERRORS are included, because those are never
        legitimate -- they mean a file failed to import and its tests are in neither
        the pass set nor the fail set.
        """
        return {
            "ispec_importable": self.ispec_importable,
            "kp_atlas_present": self.kp_atlas_present,
            "python_version": self.python_version,
            "collection_errors": self.collection_errors,
        }


def _canonical_location(root: Path) -> bool:
    try:
        return root.resolve().parent == CANONICAL_WORKTREE_PARENT.resolve()
    except OSError:
        return False


def probe(root: Path, *, python: str | None = None) -> Capabilities:
    """Measure a checkout's capabilities by RUNNING it, never by reading its config.

    ⚠️ The probe runs in a SUBPROCESS rooted at `root`. Importing the checkout's
    `config.constants` into this process would measure THIS checkout and report it as
    that one -- the same shape as the RYA-853 referee that ended up scoring its own
    values, and as the vendor echo in RYA-1035 where our own stdin came back as
    provenance. A probe that can return the prober's state is not a probe.
    """
    root = Path(root).resolve()
    python = python or sys.executable
    probe_src = (
        "import json, sys, os\n"
        "from pathlib import Path\n"
        "out = {'ispec': False, 'ispec_dir': '', 'kp': False, 'kp_dir': '',\n"
        "       'py': '%d.%d.%d' % sys.version_info[:3]}\n"
        "sys.path.insert(0, os.getcwd())\n"
        "try:\n"
        "    from config.constants import ISPEC_DIR\n"
        "    out['ispec_dir'] = str(ISPEC_DIR)\n"
        "    sys.path.insert(0, str(ISPEC_DIR))\n"
        "    import ispec\n"
        "    out['ispec'] = True\n"
        "except Exception:\n"
        "    pass\n"
        "try:\n"
        "    from scripts.measure_band_ew import KP_DIR\n"
        "    out['kp_dir'] = str(KP_DIR); out['kp'] = Path(KP_DIR).is_dir()\n"
        "except Exception:\n"
        "    env = os.environ.get('CODEX_KP_ATLAS')\n"
        "    if env:\n"
        "        out['kp_dir'] = env; out['kp'] = Path(env).is_dir()\n"
        "print(json.dumps(out))\n"
    )
    r = subprocess.run([python, "-c", probe_src], cwd=root,
                       capture_output=True, text=True)
    notes: list[str] = []
    try:
        d = json.loads(r.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        d = {"ispec": False, "ispec_dir": "", "kp": False, "kp_dir": "", "py": "?"}
        notes.append(f"probe produced no parseable result: {r.stderr.strip()[:200]}")

    caps = Capabilities(
        root=str(root),
        is_canonical_location=_canonical_location(root),
        ispec_importable=bool(d["ispec"]),
        ispec_dir=d["ispec_dir"],
        kp_atlas_present=bool(d["kp"]),
        kp_atlas_dir=d["kp_dir"],
        python_version=d["py"],
        notes=notes,
    )
    if not caps.is_canonical_location:
        caps.notes.append(
            f"checkout is not under the canonical worktree parent "
            f"({CANONICAL_WORKTREE_PARENT}); advisory only, but note that ISPEC_DIR "
            f"resolves through the parent directory, so this is the condition under "
            f"which RYA-1138's false baseline occurred")
    if not caps.ispec_importable:
        caps.notes.append(
            f"ispec is NOT importable from here (ISPEC_DIR={caps.ispec_dir!r}) — every "
            f"test that imports it becomes a COLLECTION ERROR, which is in neither the "
            f"pass set nor the fail set")
    return caps


def collect(root: Path, target: str = "tests/", *,
            python: str | None = None) -> tuple[int | None, int | None]:
    """Collected-test and collection-error counts, via `pytest --collect-only`.

    Cheap next to a full run, and it is the measurement that would have refused
    RYA-908's phantom baseline outright: 25 collection errors on one side and 0 on the
    other, before a single test executed.
    """
    r = subprocess.run(
        [python or sys.executable, "-m", "pytest", target, "-q", "--collect-only"],
        cwd=Path(root).resolve(), capture_output=True, text=True)
    return parse_collect_output((r.stdout or "") + (r.stderr or ""))


def parse_collect_output(out: str) -> tuple[int | None, int | None]:
    """Split out as a pure function so the parsing can be tested against real pytest
    output shapes without spending a full collection run to produce each one."""
    # pytest reports collection three different ways depending on what happened:
    #   "3527 tests collected in 12.34s"
    #   "3500 tests collected, 25 errors in 9.30s"
    #   "!!!! Interrupted: 25 errors during collection !!!!"   (no count at all)
    # ⚠️ The third form is the dangerous one and the form RYA-908 actually hit: there
    # is NO test count to compare, so a parser that only looks for "N tests collected"
    # returns None and a caller that treats None as "unknown, carry on" reproduces the
    # original defect. None is therefore propagated and `assert_comparable` treats a
    # None-vs-0 mismatch as a refusal like any other.
    m_tests = re.search(r"(\d+)\s+tests?\s+collected", out)
    n_tests = int(m_tests.group(1)) if m_tests else None

    m_err = (re.search(r"collected,\s*(\d+)\s+errors?", out)
             or re.search(r"(\d+)\s+errors?\s+during\s+collection", out)
             or re.search(r"^(\d+)\s+errors?\s+in\s", out, re.M))
    if m_err:
        n_errors = int(m_err.group(1))
    elif n_tests is not None:
        n_errors = 0
    else:
        n_errors = None
    return n_tests, n_errors


def assert_comparable(branch: Capabilities, baseline: Capabilities) -> None:
    """🔴 REFUSE to proceed unless the two checkouts have identical fingerprints.

    This is the loud-fail RYA-869 asks for, applied to the merge gate itself. It is
    the whole point of the ticket: a set-diff whose two sides had different
    capabilities is not a weak result, it is an INVERTED one.
    """
    fb, fl = branch.fingerprint(), baseline.fingerprint()
    if fb == fl:
        return
    diffs = [f"    {k}: branch={fb[k]!r}  baseline={fl[k]!r}"
             for k in sorted(fb) if fb[k] != fl[k]]
    raise NonComparableCheckoutsError(
        "REFUSING to emit a set-diff verdict — the two checkouts do not have the same "
        "capabilities, so a difference in their failure sets cannot be attributed to "
        "the branch.\n"
        f"  branch   {branch.root}\n"
        f"  baseline {baseline.root}\n"
        + "\n".join(diffs)
        + "\n\n  This is RYA-1138: on 2026-08-30 a baseline outside ~/codex could not "
          "import ispec, produced 25 phantom collection errors, and would have made "
          "eleven real regressions read as an IMPROVEMENT.\n"
          "  Fix the baseline (site it under the canonical worktree parent, or set "
          "ISPEC_DIR / CODEX_KP_ATLAS so both sides resolve alike) and re-run. Do not "
          "interpret the diff until the fingerprints match."
    )


def describe(caps: Capabilities) -> str:
    lines = [f"  {caps.root}",
             f"    canonical location : {caps.is_canonical_location}",
             f"    ispec importable   : {caps.ispec_importable}  ({caps.ispec_dir})",
             f"    KP atlas present   : {caps.kp_atlas_present}  ({caps.kp_atlas_dir})",
             f"    python             : {caps.python_version}"]
    if caps.collected_tests is not None:
        lines.append(f"    collected          : {caps.collected_tests} tests, "
                     f"{caps.collection_errors} collection error(s)")
    lines += [f"    ⚠️ {n}" for n in caps.notes]
    return "\n".join(lines)


def as_dict(caps: Capabilities) -> dict:
    return asdict(caps)

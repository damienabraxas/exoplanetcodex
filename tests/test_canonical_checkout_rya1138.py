"""RYA-1138 — the set-diff guard, and POSITIVE CONTROLS that it actually fires.

A guard that has never been seen to fail is indistinguishable from a guard that
cannot fail. That is not a hypothetical worry in this repo: RYA-853's referee ended up
comparing DH19 against DH19 and reported sd 0.000 as a verdict, and RYA-1080's own
value guard spent a period diffing `origin/main` against itself. Both read as passing.

So every refusal here is tested twice -- once that it permits the comparable case, and
once that it REFUSES the exact scenario from 2026-08-30, where a baseline that could
not import `ispec` reported 25 collection errors and no failures, and would have made
eleven real regressions read as an improvement.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import canonical_checkout as cc                        # noqa: E402


def _caps(**over):
    """A healthy checkout, with named overrides."""
    base = dict(root="/Users/x/codex/branch", is_canonical_location=True,
                ispec_importable=True, ispec_dir="/Users/x/codex/ispec",
                kp_atlas_present=True, kp_atlas_dir="/atlas",
                python_version="3.12.1", collected_tests=3500,
                collection_errors=0)
    base.update(over)
    return cc.Capabilities(**base)


# ── the comparable case is permitted ──────────────────────────────────────────

def test_identical_checkouts_are_comparable():
    cc.assert_comparable(_caps(), _caps(root="/Users/x/codex/base"))


def test_a_branch_that_ADDS_tests_is_still_comparable():
    """A branch legitimately collects more tests than main. Folding the collected
    COUNT into the fingerprint would make every test-adding branch unmergeable, which
    is how a guard earns the reputation that gets it disabled."""
    cc.assert_comparable(_caps(collected_tests=3502),
                         _caps(root="/base", collected_tests=3500))


def test_two_equally_NON_canonical_checkouts_are_comparable():
    """🔴 THE ANTI-PROXY CONTROL. The tempting guard is "is it under ~/codex?", and it
    is the wrong question: two checkouts that are both non-canonical IN THE SAME WAY
    still measure the same thing. Encoding location would also fail every checkout the
    day the canonical parent moves -- which is the exact event that created RYA-1090.
    """
    cc.assert_comparable(_caps(is_canonical_location=False),
                         _caps(root="/base", is_canonical_location=False))


# ── POSITIVE CONTROLS: the guard fires ────────────────────────────────────────

def test_POSITIVE_CONTROL_the_rya908_false_baseline_is_REFUSED():
    """🔴 THE EXACT 2026-08-30 SCENARIO. Baseline outside ~/codex, so `ISPEC_DIR`
    (`ROOT.parent / 'ispec'`) does not exist, 25 modules fail to import, and the
    baseline reports collection errors instead of failures. The branch had eleven real
    failures; against this baseline that reads as an IMPROVEMENT.

    Without this control the whole ticket is unverified.
    """
    with pytest.raises(cc.NonComparableCheckoutsError) as e:
        cc.assert_comparable(
            _caps(),
            _caps(root="/private/tmp/scratch/mainbase", is_canonical_location=False,
                  ispec_importable=False, collected_tests=None, collection_errors=25))
    msg = str(e.value)
    assert "ispec_importable" in msg
    assert "collection_errors" in msg
    assert "REFUSING" in msg


def test_POSITIVE_CONTROL_a_CANONICAL_but_broken_baseline_is_still_refused():
    """Location is not competence. A baseline sitting in exactly the right place that
    nonetheless cannot import `ispec` produces exactly the same inverted verdict, and
    a location-based guard waves it through."""
    with pytest.raises(cc.NonComparableCheckoutsError):
        cc.assert_comparable(_caps(),
                             _caps(root="/Users/x/codex/base",
                                   is_canonical_location=True,
                                   ispec_importable=False))


def test_POSITIVE_CONTROL_a_missing_atlas_is_refused():
    """`CODEX_KP_ATLAS` absent kills the Mac suite at COLLECTION. Same class."""
    with pytest.raises(cc.NonComparableCheckoutsError):
        cc.assert_comparable(_caps(), _caps(root="/base", kp_atlas_present=False))


def test_POSITIVE_CONTROL_collection_errors_alone_are_refused():
    """Even with every dependency resolving, a side with collection errors is not
    comparable: those modules' tests are in NEITHER the pass set nor the fail set, so
    the diff silently stops covering them."""
    with pytest.raises(cc.NonComparableCheckoutsError):
        cc.assert_comparable(_caps(), _caps(root="/base", collection_errors=3))


def test_an_unknown_collection_count_is_not_silently_accepted():
    """⚠️ pytest's `Interrupted: N errors during collection` form prints NO test count.
    A parser returning None there, and a caller reading None as "unknown, carry on",
    reproduces the original defect exactly."""
    n_tests, n_errors = cc.parse_collect_output(
        "!!!!!! Interrupted: 25 errors during collection !!!!!!")
    assert n_tests is None and n_errors == 25
    with pytest.raises(cc.NonComparableCheckoutsError):
        cc.assert_comparable(_caps(), _caps(root="/base", collection_errors=25))


@pytest.mark.parametrize("text,expected", [
    ("3527 tests collected in 12.34s", (3527, 0)),
    ("3500 tests collected, 25 errors in 9.30s", (3500, 25)),
    ("!!! Interrupted: 25 errors during collection !!!", (None, 25)),
    ("1 test collected in 0.1s", (1, 0)),
])
def test_collect_output_parses_every_pytest_shape(text, expected):
    assert cc.parse_collect_output(text) == expected


def test_the_parser_is_not_fooled_by_ITS_OWN_test_ids():
    """🔴 THE INSTRUMENT ENTERED ITS OWN MEASUREMENT — a real bug, found by this
    harness refusing a healthy branch.

    `--collect-only` prints every test id BEFORE the summary. The parametrized cases
    above carry ids containing "3527 tests collected" and "25 errors during
    collection", so a parser using `re.search` over the whole output -- which returns
    the FIRST match -- read its own fixtures as the summary and reported 3527 tests
    with 25 collection errors for a run that collected 3563 with none. The set-diff
    then refused a perfectly good branch, citing numbers that came from itself.

    This is RYA-1112's lesson in a new place, and the reason the parser now scans from
    the END with whole-line anchoring rather than searching the whole blob.
    """
    decoy = "\n".join([
        "tests/test_canonical_checkout_rya1138.py::test_shapes[3527 tests collected in 12.34s-0]",
        "tests/test_canonical_checkout_rya1138.py::test_shapes[25 errors during collection]",
        "tests/test_canonical_checkout_rya1138.py::test_shapes[3500 tests collected, 25 errors in 9.30s]",
        "",
        "3563 tests collected in 2.45s",
    ])
    assert cc.parse_collect_output(decoy) == (3563, 0)


def test_the_parser_reads_the_LAST_summary_not_the_first():
    """The general form of the above: whatever precedes it, the trailing summary wins."""
    assert cc.parse_collect_output(
        "1 tests collected in 0.1s\n...noise...\n999 tests collected in 9s") == (999, 0)


# ── the probe cannot echo the prober ──────────────────────────────────────────

def test_the_probe_measures_the_TARGET_not_the_caller():
    """🔴 A probe that can return the prober's state is not a probe (RYA-1035's vendor
    echo: our own stdin came back through a binary and was recorded as provenance).
    Probing a directory that is not a checkout at all must report NOT importable, even
    though `ispec` is importable from the process running this test."""
    caps = cc.probe(Path("/"))
    assert caps.ispec_importable is False
    assert caps.root == "/"


# ── the resolver audit (item 5) ───────────────────────────────────────────────

def test_the_resolver_audit_FINDS_the_real_ispec_escape():
    """Non-vacuity: the audit must actually detect the resolver that caused this
    ticket. An audit reporting zero findings is the same output as a broken audit."""
    from scripts.rya1138_worktree_resolver_audit import audit
    hits = {(f["file"], f["line"]) for f in audit()}
    assert ("config/constants.py", 1119) in hits, (
        "the audit no longer finds `ISPEC_DIR = ROOT.parent / 'ispec'`, the resolver "
        "that produced RYA-1138's false baseline")


def test_the_resolver_audit_does_NOT_flag_correct_in_repo_resolution():
    """🔴 THE OFF-BY-ONE CONTROL, and it caught a real bug in this audit.

    `Path(__file__).resolve().parent.parent` from `pipeline/x.py` lands EXACTLY on the
    repo root and is correct. The first version of this audit counted the file as
    depth 1 instead of 2 and reported 204 escapes, nearly all of them correct code --
    a finding list that size does not get read, and the two that mattered would have
    been lost in it.
    """
    from scripts.rya1138_worktree_resolver_audit import audit
    findings = audit()
    files = {f["file"] for f in findings}
    for correct in ("pipeline/anchor_pools.py", "pipeline/gf_resolver.py",
                    "pipeline/audit/gf_store_consistency.py"):
        assert correct not in files, (
            f"{correct} resolves to its OWN repo root and must not be flagged; the "
            f"audit's depth arithmetic is off by one again")
    assert len(findings) < 40, (
        f"{len(findings)} findings — this is meant to be a reviewable list of real "
        f"escapes, not a grep dump")


def test_the_resolver_audit_excludes_itself_BY_NAME():
    """⚠️ The instrument must not enter its own measurement, and the exclusion must be
    by NAME. RYA-1112's auditor matched its own grep the moment it was tracked, and a
    pattern-based exclusion would also silence unrelated files that merely resemble
    it."""
    src = (ROOT / "scripts" / "rya1138_worktree_resolver_audit.py").read_text()
    assert 'SELF = "scripts/rya1138_worktree_resolver_audit.py"' in src
    from scripts.rya1138_worktree_resolver_audit import audit
    assert not [f for f in audit()
                if f["file"] == "scripts/rya1138_worktree_resolver_audit.py"]


# ── the harness must not measure NOTHING and call it agreement ────────────────

def test_POSITIVE_CONTROL_a_vacuous_failure_set_is_REFUSED():
    """🔴 THE BUG THIS HARNESS SHIPPED WITH, AND THE ONE IT EXISTS TO PREVENT.

    `-r` REPLACES pytest's default report characters rather than adding to them. The
    first version passed `-rs` to collect skips and thereby switched OFF the default
    `fE`, so pytest printed no FAILED lines, the regex matched nothing, and BOTH sides
    returned an empty failure set. The harness reported "failure set IDENTICAL --
    nothing moved" while the two runs were 10 failed and 8 failed.

    A guard whose two sides converge, inside the tool built to stop exactly that. The
    cross-check below is the structural fix: pytest's own summary count must agree
    with what was parsed out of its output, or the run is refused.
    """
    from scripts.suite_set_diff import run_suite, VacuousMeasurementError
    import scripts.suite_set_diff as ssd

    real = ssd.subprocess.run

    class _Fake:
        returncode = 1
        # A summary claiming failures, with the FAILED lines suppressed -- exactly
        # what `-rs` produced.
        stdout = "........F...\n\n10 failed, 3484 passed, 51 skipped in 369.62s\n"
        stderr = ""

    ssd.subprocess.run = lambda *a, **k: _Fake()
    try:
        with pytest.raises(VacuousMeasurementError) as e:
            run_suite(Path("/tmp"), "tests/", sys.executable)
        assert "10 failed" in str(e.value)
        assert "0 were parsed" in str(e.value)
    finally:
        ssd.subprocess.run = real


def test_the_suite_runner_asks_pytest_for_failures_AND_skips():
    """The one-character root cause, pinned. `-rs` alone silences FAILED lines."""
    src = (ROOT / "scripts" / "suite_set_diff.py").read_text()
    assert '"-rfEs"' in src, "the report flags must request f, E and s together"
    assert '"-rs"' not in src, (
        "`-rs` REPLACES the default `fE` — it silences the FAILED lines the set-diff "
        "is built from")


# ── the harness refuses end-to-end, not just in unit form ─────────────────────

def test_the_set_diff_harness_REFUSES_a_non_comparable_baseline_end_to_end():
    """The unit controls above test `assert_comparable`. This runs the actual CLI
    against a baseline that is not a checkout at all, and requires a REFUSAL exit code
    (2) rather than a verdict -- because the failure mode being guarded is a harness
    that prints a normal-looking answer."""
    r = subprocess.run(
        [sys.executable, "scripts/suite_set_diff.py",
         "--baseline-root", "/", "--target", "tests/test_canonical_checkout_rya1138.py"],
        cwd=ROOT, capture_output=True, text=True, timeout=600)
    assert r.returncode == 2, (
        f"expected a REFUSAL (exit 2), got {r.returncode}:\n{r.stdout[-3000:]}")
    assert "REFUSING" in r.stdout

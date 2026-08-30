"""
tests/test_feed_repo_reconciliation_rya1080.py — RYA-1080
=========================================================
The feed publishes; `data/results/band_products/` is what people grep. When the two
disagree, a check of the tree answers a question about the feed with an answer about the
disk — which cost two wrong diagnoses in one session, the second ("Fe II has never been
re-run") nearly triggering a from-scratch re-run of finished work.

The positive control is the point of this file: a guard that has never been shown to fail
is not evidence of anything.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import feed_repo_reconciliation as frr  # noqa: E402

FEED = ROOT / "data" / "products" / "solar" / "Fe.json"
BLOCKED = ROOT / "data" / "results" / "rya1080" / "rya1080_blocked.csv"

#: RESOLVED by RYA-1084 — this set is now EMPTY, and that is the assertion.
#: RYA-1080 blocked ten products whose source file's sha256 no longer matched the feed's
#: record. RYA-1084 identified the writer (a deterministic scripts/derive_band_products.py
#: run: the 2026-08-27 13:34 run reproduced RYA-1051's committed bytes exactly), classified
#: the diffs (six a new `deck` column from RYA-1044/1045; one a 1-ULP rounding boundary,
#: now pinned by pipeline.error_budget.round_dex), and re-ingested them. The evidence of
#: what was wrong lives in data/results/rya1080/rya1080_blocked.csv.
KNOWN_BLOCKED: set[tuple] = set()


@pytest.fixture(scope="module")
def findings():
    return frr.check()


# ── the guard's own contract ──────────────────────────────────────────────────

def test_every_live_product_has_a_committed_checksum_matching_artifact(findings):
    """Was a strict xfail while RYA-1080 held ten products; RYA-1084 cleared their writer
    and re-ingested them, so this is a plain assertion again. 75 of 75."""
    blocking = [f for f in findings if f.blocking]
    assert not blocking, frr.report(findings)


def test_the_blocking_set_is_exactly_the_known_ten(findings):
    """The xfail above says 'some drift is known'. THIS says WHICH — so an eleventh
    unreconciled product is a hard failure today, not a slightly larger number nobody
    reads."""
    got = {f.key() for f in findings if f.blocking}
    assert got == KNOWN_BLOCKED, (
        f"the blocking set has changed.\n  new: {sorted(got - KNOWN_BLOCKED)}\n"
        f"  gone: {sorted(KNOWN_BLOCKED - got)}")


def test_the_guard_exits_zero_now_that_nothing_is_blocked():
    """CI reads the exit code, so it has to be right in BOTH directions. That it goes
    non-zero on drift is proved by the positive controls below, not asserted here."""
    r = subprocess.run([sys.executable, "pipeline/feed_repo_reconciliation.py"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout


# ── POSITIVE CONTROLS — required by the ticket ────────────────────────────────

def _sandbox(tmp_path: Path) -> Path:
    """A throwaway repo-shaped tree with one reconciled product, so a control can break
    something without touching the real tree."""
    doc = json.loads(FEED.read_text())
    good = next(p for p in doc["products"]
                if p["provenance"].get("copied_to")
                and (ROOT / p["provenance"]["copied_to"]).exists())
    rel = good["provenance"]["copied_to"]
    (tmp_path / Path(rel).parent).mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / rel, tmp_path / rel)
    feed_dir = tmp_path / "data" / "products" / "solar"
    feed_dir.mkdir(parents=True, exist_ok=True)
    (feed_dir / "Fe.json").write_text(json.dumps({"products": [good]}, indent=2))
    # git-track it, since the guard requires tracking and band_products is gitignored
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "add", "-f", rel], cwd=tmp_path, check=True)
    return tmp_path


def test_control_the_sandbox_passes_before_it_is_broken(tmp_path):
    """The control's own control. If the untouched sandbox already failed, breaking it
    would prove nothing."""
    s = _sandbox(tmp_path)
    assert not [f for f in frr.check(root=s) if f.blocking]


def test_control_a_renamed_artifact_makes_the_guard_fail(tmp_path):
    """POSITIVE CONTROL 1 — the exact drift this ticket exists to catch: the feed says the
    product is there, the file is not."""
    s = _sandbox(tmp_path)
    rel = json.loads((s / "data/products/solar/Fe.json").read_text())["products"][0][
        "provenance"]["copied_to"]
    (s / rel).rename(s / (rel + ".renamed"))
    kinds = {f.kind for f in frr.check(root=s) if f.blocking}
    assert "ARTIFACT_ABSENT" in kinds, kinds


def test_control_a_modified_artifact_makes_the_guard_fail(tmp_path):
    """POSITIVE CONTROL 2 — the file is present and is not the one the feed measured. This
    is the case that A-value equality would have waved through."""
    s = _sandbox(tmp_path)
    rel = json.loads((s / "data/products/solar/Fe.json").read_text())["products"][0][
        "provenance"]["copied_to"]
    (s / rel).write_bytes((s / rel).read_bytes() + b"\n")
    kinds = {f.kind for f in frr.check(root=s) if f.blocking}
    assert "SHA_MISMATCH" in kinds, kinds


def test_control_a_copied_to_outside_the_repo_is_not_accepted(tmp_path):
    """POSITIVE CONTROL 3 — 🔴 the trap that made the exposure look like 67 instead of 75.
    Eight rows carried a non-null copied_to pointing into /private/tmp/g3d/ and
    /private/tmp/sirius_orphans/, which READS as reconciled and is not committed at all."""
    s = _sandbox(tmp_path)
    f = s / "data/products/solar/Fe.json"
    doc = json.loads(f.read_text())
    doc["products"][0]["provenance"]["copied_to"] = "/private/tmp/g3d/somewhere.csv"
    f.write_text(json.dumps(doc, indent=2))
    kinds = {x.kind for x in frr.check(root=s) if x.blocking}
    assert "COPIED_TO_OUTSIDE_REPO" in kinds, kinds


def test_control_an_untracked_artifact_makes_the_guard_fail(tmp_path):
    """POSITIVE CONTROL 4 — band_products is gitignored (.gitignore:87), so a file being
    present on disk is NOT the same as it being committed. Present-but-untracked must
    fail, or a local run would satisfy the guard for everyone else."""
    s = _sandbox(tmp_path)
    rel = json.loads((s / "data/products/solar/Fe.json").read_text())["products"][0][
        "provenance"]["copied_to"]
    subprocess.run(["git", "rm", "--cached", "-q", rel], cwd=s, check=True)
    kinds = {f.kind for f in frr.check(root=s) if f.blocking}
    assert "ARTIFACT_UNTRACKED" in kinds, kinds


# ── the ticket's CRITICAL conditions ──────────────────────────────────────────

#: The state of the feed immediately BEFORE RYA-1080's reconciliation — the last commit
#: on main that the reconcile branch was cut from. PINNED, not `origin/main`.
#:
#: 🔴 It was `origin/main` and that made this test VACUOUS the moment RYA-1080 merged: it
#: compared `git show origin/main:Fe.json` against the working file, which are the same
#: file once the branch is in. A guard whose two sides converge cannot fail — the same
#: shape as the RYA-853 referee that ended up scoring its own values.
BASELINE_SHA = "6c1529fa7561f347b62eefcc01ec017ee78fd5c4"

#: Keys RYA-1080 and RYA-1084 are ALLOWED to write. Everything else must match baseline.
_RECONCILE_KEYS = {"copied_to", "sha256", "reingested_by", "reingest_reason"}

#: RYA-1084 re-ingested ten products after clearing their writer, and exactly ONE published
#: field legitimately moved with them. Named, not tolerated as a class.
SANCTIONED = {("Fe", "I", "VIS", "sigma_stat"): (0.0217, 0.0218)}

#: 🔴 RYA-1100 — `display` IS NOT A PUBLISHED VALUE, AND THIS GUARD NO LONGER PINS IT.
#: This check's own docstring scopes it to "was a published NUMBER edited?" -- reconcile
#: the artifacts, not the numbers. `display` is neither: it is a LABEL COMPUTED from
#: `treatment`, `route` and `gf` by `publish_product.display_name`, all of which remain
#: fully guarded here. Pinning a derived field to a baseline makes it impossible to fix
#: the CODE THAT DERIVES IT, and that is not hypothetical -- it blocked the correction of
#: 22 live ENGINE-A products that were published as "EW" measurements while running on
#: route=SYNTH, plus three ⟨3D⟩ products rendering as raw tokens.
#:
#: This is NOT the guard being widened to let a change through, which is the antipattern
#: RYA-1088 refused to commit. Ownership MOVES, and to a stricter test: RYA-1100's
#: `test_every_published_display_is_DERIVABLE_from_its_own_axes` asserts every live
#: `display` EQUALS what the axis registry derives for that row. A pin is satisfied by
#: never touching the field; a derivability assertion cannot be satisfied by a hand edit
#: at all, so the field is more constrained after this change than before it.
DERIVED_FIELDS = {"display"}


def _index(doc: dict) -> dict:
    out = {}
    for p in doc["products"]:
        prov = {k: v for k, v in p["provenance"].items() if k not in _RECONCILE_KEYS}
        key = (p["element"], p["ion"], p["band"], p["instrument"], p.get("tier"),
               p["treatment"], prov.get("path"))
        out[key] = {**{k: v for k, v in p.items() if k != "provenance"}, "_prov": prov}
    return out


def published_value_edits(baseline: dict, current: dict) -> dict:
    """Fields that EXISTED at baseline and have since changed value.

    ⚠️ Scoped deliberately. This asks "was a published number edited?", NOT "has the
    product schema grown?". A key ADDED later is another ticket's business — RYA-1088's
    `sigma_params` is the case in point, and the unscoped version of this check failed on
    it and forced that ticket to revert its own deliverable. Absent-at-baseline is
    therefore skipped; present-and-changed is the defect.
    """
    before, after = _index(baseline), _index(current)
    edits = {}
    for key, b in before.items():
        a = after.get(key)
        if a is None:          # a product removed later is a different question
            continue
        for fk, bv in b.items():
            if fk == "_prov" or fk in DERIVED_FIELDS or fk not in a:
                continue
            if a[fk] != bv:
                edits[(key[0], key[1], key[2], fk)] = (bv, a[fk])
    return edits


def _baseline_doc():
    r = subprocess.run(["git", "show", f"{BASELINE_SHA}:data/products/solar/Fe.json"],
                       cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        pytest.skip(f"baseline {BASELINE_SHA[:9]} not available (shallow clone?)")
    return json.loads(r.stdout)


def test_no_published_value_was_edited_to_reconcile():
    """🔴 CRITICAL: reconcile the artifacts, not the numbers."""
    edits = published_value_edits(_baseline_doc(), json.loads(FEED.read_text()))
    unexpected = {k: v for k, v in edits.items() if SANCTIONED.get(k) != v}
    assert not unexpected, (
        f"a published value changed outside RYA-1084's re-ingest: {unexpected}")


def test_control_the_baseline_actually_differs_from_today():
    """The control that stops this going vacuous again. If the pinned baseline and the
    live feed were identical the comparison would pass by construction and prove nothing —
    so assert the two differ SOMEWHERE (they must: RYA-1080 wrote copied_to on 75 rows)."""
    assert _baseline_doc() != json.loads(FEED.read_text()), (
        "baseline and live feed are identical — this guard has gone vacuous")


def test_control_an_edited_published_value_is_caught():
    """POSITIVE CONTROL. Move an abundance and the check must see it."""
    live = json.loads(FEED.read_text())
    live["products"][0]["A"] = float(live["products"][0]["A"]) + 0.1
    edits = published_value_edits(_baseline_doc(), live)
    assert any(fk == "A" for (*_, fk) in edits), edits


def test_control_a_newly_added_field_is_NOT_caught():
    """POSITIVE CONTROL for the scoping itself — RYA-1088's exact case. Adding a field to
    a product is a schema addition, not an edited value, and must not fail this guard."""
    live = json.loads(FEED.read_text())
    live["products"][0]["sigma_params"] = 0.001
    live["products"][0]["sigma_params_reason"] = "solar logg & [Fe/H] fixed by definition"
    edits = published_value_edits(_baseline_doc(), live)
    unexpected = {k: v for k, v in edits.items() if SANCTIONED.get(k) != v}
    assert not unexpected, (
        f"adding a field tripped the value guard — it is still over-scoped: {unexpected}")


def test_control_the_DERIVED_field_exemption_is_NARROW():
    """RYA-1100. The exemption must cover the derived label and nothing else -- above all
    not a measurement. A guard that skips a field must prove what it still catches."""
    assert DERIVED_FIELDS == {"display"}
    live = json.loads(FEED.read_text())
    for fk in ("A", "sigma_stat", "sigma_syst", "n_lines", "treatment", "route"):
        assert fk not in DERIVED_FIELDS, f"{fk} is a published value, not a derived label"
    # and the exempted field is not simply unguarded: RYA-1100 asserts it is DERIVABLE.
    from scripts.publish_product import display_name
    for p in live["products"]:
        # RYA-1106: ask for the product's OWN pool, not a hardcoded "kurucz". Forcing
        # kurucz here overrode the pool LEGACY declares for the label, so a product
        # whose axis is gf="lab" (the Amarsi leg, after RYA-1104 measured 67/67 lines
        # primary-laboratory) derived a name without its `· lab-gf` suffix and this
        # assertion compared the feed against a name nothing publishes.
        assert p["display"] == display_name(
            p["treatment"], gf=p.get("gf"), route=p["route"])


def test_regenerability_gaps_are_recorded_not_silent(findings):
    """Committing the artifact is the floor; regenerability is the net (RYA-1011). Where
    the holding is Mac-local the product cannot be rebuilt on the committed runner, and
    that must appear in the guard's own output rather than being discovered later."""
    gaps = [f for f in findings if f.kind == "REGENERABILITY_GAP"]
    assert all("RYA-1011" in f.detail for f in gaps)

    # 🔴 RYA-1100 — WAS `assert len(gaps) >= 50`, A CONSTANT PINNED TO THE LIVE
    # POPULATION. RYA-1100 withdrew 13 duplicate products to `superseded[]` and the count
    # fell to 46, failing a check that nothing had actually broken: the threshold was
    # measuring how many products exist, not whether gaps are recorded.
    #
    # Replaced with the invariant the generator actually implements -- ONE gap per live
    # product produced on the Mac (`pipeline/feed_repo_reconciliation.py`, the
    # `host == "mac"` branch). This is EXACT rather than a floor, so it is strictly
    # stronger than the constant: it now catches a gap going missing for a single
    # product, which `>= 50` could never see, and it cannot be satisfied by the
    # population merely being large. NOT a lowered tolerance (RYA-161) -- a different and
    # tighter assertion, and it would fail if RYA-1100 had dropped a gap rather than a
    # duplicate.
    live = json.loads(FEED.read_text())["products"]
    on_mac = [p for p in live if str(p["provenance"].get("host", "")).lower() == "mac"]
    assert len(gaps) == len(on_mac), (
        f"{len(on_mac)} live products were produced on the Mac but {len(gaps)} "
        f"regenerability gaps were recorded — a gap went silent (RYA-1011)")
    assert gaps, "no regenerability gaps at all — the check has gone vacuous"


def test_the_blocked_ledger_records_both_checksums():
    """A blocked product must say what the feed recorded AND what is on disk, or the next
    person cannot tell a cosmetic rewrite from a changed measurement."""
    if not BLOCKED.exists():
        pytest.skip("blocked ledger absent")
    import pandas as pd
    d = pd.read_csv(BLOCKED)
    assert len(d) == 10
    assert d.feed_sha256.notna().all() and d.actual_sha256.notna().all()
    assert (d.feed_sha256 != d.actual_sha256).all()

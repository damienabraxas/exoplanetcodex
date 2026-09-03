"""RYA-1178 — the re-emitted Fe feed, and the xi pairing guard that feeds its derivatives."""
import json
import math
from pathlib import Path

import pytest

from pipeline.xi_pairing import (LEG_MINUS, LEG_PLUS, XiPairingError, assert_pair,
                                 read_stamp, span_kms, write_stamp)

NOMINAL, STEP = 1.0, 0.10


def _pair(tmp_path, *, lo=None, hi=None, lo_leg=LEG_MINUS, hi_leg=LEG_PLUS):
    a, b = tmp_path / "minus", tmp_path / "plus"
    a.mkdir(); b.mkdir()
    write_stamp(a, xi_kms=NOMINAL - STEP if lo is None else lo, leg=lo_leg,
                xi_nominal=NOMINAL, step_kms=STEP) if lo is None else \
        (a / "xi_run.json").write_text(json.dumps(
            {"xi_kms": lo, "leg": lo_leg, "xi_nominal_kms": NOMINAL, "step_kms": STEP}))
    write_stamp(b, xi_kms=NOMINAL + STEP if hi is None else hi, leg=hi_leg,
                xi_nominal=NOMINAL, step_kms=STEP) if hi is None else \
        (b / "xi_run.json").write_text(json.dumps(
            {"xi_kms": hi, "leg": hi_leg, "xi_nominal_kms": NOMINAL, "step_kms": STEP}))
    return a, b


# ── scope A: the pairing guard ────────────────────────────────────────────────
def test_a_correct_pair_is_accepted_and_yields_the_real_span(tmp_path):
    """The POSITIVE control. A guard that only ever refuses proves nothing."""
    a, b = _pair(tmp_path)
    lo, hi = assert_pair(a, b, xi_nominal=NOMINAL, step_kms=STEP)
    assert lo["xi_kms"] == pytest.approx(0.9) and hi["xi_kms"] == pytest.approx(1.1)
    assert span_kms(lo, hi) == pytest.approx(0.2)


def test_a_mis_paired_run_is_refused(tmp_path):
    """🔴 THE REQUIRED REFUSAL TEST. A leg run at the WRONG xi still produces per-line
    abundances that difference to a plausible number — nothing downstream could tell.
    Before the stamp, only worktree isolation stood between this and a published
    derivative; the near-UV worktree collision is what that is worth."""
    a, b = _pair(tmp_path, hi=1.30)          # plus leg run at +0.30, not +0.10
    with pytest.raises(XiPairingError, match="plus leg was run at xi"):
        assert_pair(a, b, xi_nominal=NOMINAL, step_kms=STEP)


def test_a_missing_stamp_is_refused_not_guessed_from_the_directory(tmp_path):
    a, b = _pair(tmp_path)
    (a / "xi_run.json").unlink()
    with pytest.raises(XiPairingError, match="does not record which xi"):
        assert_pair(a, b, xi_nominal=NOMINAL, step_kms=STEP)


def test_a_duplicated_leg_is_refused(tmp_path):
    """Two copies of one leg difference to ~0 — 'xi does not matter', stated with
    total confidence. The RYA-853 vacuous-guard shape."""
    a, b = _pair(tmp_path, hi=0.9, hi_leg=LEG_MINUS)
    with pytest.raises(XiPairingError, match="both runs are stamped leg"):
        assert_pair(a, b, xi_nominal=NOMINAL, step_kms=STEP)


def test_a_swapped_legs_are_refused_because_they_flip_the_sign(tmp_path):
    a, b = _pair(tmp_path, lo=1.1, lo_leg=LEG_PLUS, hi=0.9, hi_leg=LEG_MINUS)
    with pytest.raises(XiPairingError, match="wrong way round"):
        assert_pair(a, b, xi_nominal=NOMINAL, step_kms=STEP)


def test_a_a_stamp_cannot_contradict_its_own_leg(tmp_path):
    d = tmp_path / "leg"; d.mkdir()
    with pytest.raises(XiPairingError, match="refusing to stamp"):
        write_stamp(d, xi_kms=1.1, leg=LEG_MINUS, xi_nominal=NOMINAL, step_kms=STEP)


def test_a_the_campaign_derivative_refuses_an_unpaired_run(tmp_path):
    """The guard is wired into dA_dxi, not merely available beside it."""
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from rya1120_xi_campaign import dA_dxi
    a, b = _pair(tmp_path, hi=1.30)
    with pytest.raises(XiPairingError):
        dA_dxi([], [], step_kms=STEP, minus_dir=a, plus_dir=b, xi_nominal=NOMINAL)


# ── the re-emitted feed ───────────────────────────────────────────────────────
#: 🔴 ABSOLUTE, NEVER RELATIVE. Other modules in this suite chdir, so a relative read
#: passes when this file runs alone and dies in the full run.
ROOT = Path(__file__).resolve().parents[1]
FEED = ROOT / "data/products/solar/Fe.json"


@pytest.fixture(scope="module")
def feed():
    return json.loads(FEED.read_text())


def test_every_product_carries_the_part1_schema(feed):
    need = ["atlas", "telluric_state", "line_set_resolved", "grade", "wavelength_range_A",
            "generated_at", "code_commit", "xi_value_kms", "delta_xi_kms", "xi_state",
            "sigma_syst_complete", "sigma_reported", "science_provenance", "provenance"]
    for p in feed["products"]:
        for f in need:
            assert p.get(f) not in (None, "", [], {}), f"{f} missing on {p['holding']}/{p['treatment']}"


def test_line_set_stays_derived_and_is_never_stored_on_our_own_products(feed):
    """🔴 RYA-1127's call, not this ticket's. The ticket asked to populate `line_set`;
    `test_line_set_identity_rya1127.py` asserts `"line_set" not in ours` because a stored
    copy is a second source of truth free to drift from the one the identity key resolves.
    The readable value is published as `line_set_resolved`, explicitly labelled derived."""
    from pipeline.model_registry import LINE_SETS
    from pipeline.reference_lineset import line_set_for_product
    for p in feed["products"]:
        assert "line_set" not in p, "our own products must not STORE the axis (RYA-1127)"
        assert p["line_set_resolved"] in LINE_SETS
        assert p["line_set_resolved"] == line_set_for_product(p)
        assert "DERIVED" in p["line_set_basis"]


def test_the_retired_consistent_tier_never_appears(feed):
    """RYA-1105 retired it; model_registry.LINE_SETS omits it deliberately."""
    blob = json.dumps(feed)
    assert '"consistent"' not in blob.lower().replace("self-consistent", "")
    for p in feed["products"]:
        assert "Consistent" not in str(p.get("grade"))


def test_sigma_reported_is_the_quadrature_it_claims(feed):
    for p in feed["products"]:
        st, sy, rep = p.get("sigma_stat"), p.get("sigma_syst_complete"), p.get("sigma_reported")
        want = math.sqrt(sum(t * t for t in (st, sy) if t))
        assert rep == pytest.approx(want, abs=1e-6)


def test_sigma_xi_enters_exactly_once(feed):
    """Through sigma_syst_complete, never also added separately — double counting a term
    is as wrong as omitting it, and harder to see."""
    for p in feed["products"]:
        sx, sy, comp = p.get("sigma_xi"), p.get("sigma_syst"), p.get("sigma_syst_complete")
        if sx and sy:
            assert comp == pytest.approx(math.hypot(sy, sx), abs=1e-6)


def test_a_bar_missing_its_xi_term_says_so(feed):
    for p in feed["products"]:
        if p.get("xi_state") in ("UNMEASURED", "NOT_IN_CAMPAIGN"):
            assert "LOWER BOUND" in p.get("sigma_reported_caveat", "")


def test_aliased_means_a_provably_different_pool(feed):
    """Not merely 'shares a band-less campaign key' — that is ambiguity, and it has its
    own state so the two are not confused."""
    for p in feed["products"]:
        if p.get("xi_state") == "ALIASED":
            assert "was measured on a pool of" in p["xi_note"]
        if p.get("xi_state") == "MEASURED_KEY_AMBIGUOUS":
            assert "cannot be proven" in p["xi_note"]


def test_part3_is_deferred_rather_than_published_around_the_gate(feed):
    """🔴 Part 3 asked to fold the four RYA-1106 Asplund products into the feed. Appending
    them from this script was TRIED and it broke 24 tests — RYA-1092's eligibility gate,
    RYA-1034's 'every product states where it came from', and RYA-1111's own entrypoint
    guard each caught a row that had bypassed them. A product enters through the canonical
    publisher or not at all. The assembled record is still built, so the follow-up has it."""
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    import rya1178_emit_fe_schema as M
    assert not [p for p in feed["products"] if p.get("line_set") == "asplund"]
    hold, inst, models, _ = M.load_sources()
    rows = M.asplund_products(models, hold)
    assert len(rows) == 4
    for r in rows:
        assert r["line_set"] == "asplund" and r["grade"] == "Asplund Grade"
        assert r["science_provenance"]["gf_citation"].startswith("Asplund")


def test_the_schema_bump_is_additive_only(feed):
    """/2 is only a safe superset if nothing was dropped. The one reader that pins the
    identifier accepts both BECAUSE of that, so the claim has to be checked."""
    import subprocess
    old = json.loads(subprocess.run(
        ["git", "show", "HEAD:data/products/solar/Fe.json"],
        cwd=ROOT, capture_output=True, text=True, check=True).stdout)
    assert feed["schema"] == "codex.element_product/2"
    assert old["schema"] == "codex.element_product/1"
    assert set(old.keys()) <= set(feed.keys())
    for a, b in zip(old["products"], feed["products"]):
        assert set(a.keys()) <= set(b.keys()), "a field was dropped — /2 is not a superset"


def test_ir_products_name_their_irreducible(feed):
    """RYA-1178 B / RYA-1164 — a wide bar with a stated cause, not a silent one."""
    ir = [p for p in feed["products"] if p["band"] == "NIR"]
    assert len(ir) == 6
    got = {}
    for p in ir:
        d = p["irreducible_dispersion"]
        assert "laboratory gf" in d["reducible_by"]
        assert d["dispersion_dex"] == pytest.approx(
            p["sigma_stat"] * math.sqrt(p["n_lines"]), abs=5e-4)
        got[(p["holding"], p["treatment"])] = d["dispersion_dex"]
    assert got[("solar_crires_plus_y_wide_rya1054", "1D-LTE")] == pytest.approx(0.138, abs=1e-3)
    assert got[("solar_iag", "1D-LTE")] == pytest.approx(0.794, abs=1e-3)
    assert got[("solar_kpno_molecfit_corrected", "1D-LTE")] == pytest.approx(0.890, abs=1e-3)
    assert got[("solar_kpno_molecfit_corrected", "ENGINE-A")] == pytest.approx(1.315, abs=1e-3)


def test_the_widest_ir_bars_carry_the_most_lines(feed):
    """The RCA, asserted: the bar tracks gf quality, not small-n. If this ever inverts,
    'more lines will not fix it' has stopped being true and the note must be re-derived."""
    ir = {(p["holding"], p["treatment"]): p for p in feed["products"] if p["band"] == "NIR"}
    crires = ir[("solar_crires_plus_y_wide_rya1054", "1D-LTE")]
    kp = ir[("solar_kpno_molecfit_corrected", "1D-LTE")]
    assert crires["n_lines"] < kp["n_lines"]
    assert crires["irreducible_dispersion"]["dispersion_dex"] < \
           kp["irreducible_dispersion"]["dispersion_dex"] / 5


def test_kp_nir_graded_tier_is_confirmed_against_its_pool(feed):
    """RYA-1114 F4. The pool is reproduced, not asserted — and it comes back 100%
    primary-laboratory, so the GRADED tier is earned."""
    kp = [p for p in feed["products"]
          if p["band"] == "NIR" and p["holding"] == "solar_kpno_molecfit_corrected"
          and p["treatment"] == "1D-LTE"][0]
    tp = kp["tier_provenance"]
    assert tp["pool_reproduced"] and tp["fraction_primary_lab"] == 1.0
    assert tp["tier_verdict"].startswith("CONFIRMED")
    assert set(tp["gf_references"]) == {"PRIMARY LAB DenHartog2014", "PRIMARY LAB Ruffoni2014"}


def test_every_withdrawn_record_states_why(feed):
    for sec, keys in (("superseded", ("superseded_reason", "superseded_reason_code")),
                      ("quarantine", ("quarantine_reason",)),
                      ("archive", ("superseded_reason", "quarantine_reason", "quarantine_codes"))):
        for r in feed[sec]:
            assert any(str(r.get(k) or "").strip() for k in keys), f"{sec} row has no reason"


def test_the_mean3d_pair_is_documented_not_silently_identical(feed):
    """Part 0. The molecfit pair publishes the same A; the record must say why."""
    pair = [p for p in feed["products"]
            if p["holding"] == "solar_kpno_molecfit_corrected"
            and "mean3D" in p["treatment"]]
    assert len(pair) == 2
    assert pair[0]["A"] == pair[1]["A"]
    for p in pair:
        assert "median" in str(p.get("aggregate_collision_note", "")).lower()

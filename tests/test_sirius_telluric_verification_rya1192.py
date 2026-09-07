"""RYA-1192 (Sirius leg) — pin what the RAW comparison settled, and what it did not.

VERIFICATION ONLY; the fixes are separate (RYA-1191, RYA-1194).

The Mac leg left 131 lines unjudged and gave CRIRES+ a proxy verdict it explicitly
labelled "not proof". These tests are written around the four ways the Sirius answer can
be misread:

  * treating a CLEARED window as if the whole arm were cleared, or the reverse;
  * reading a quiet correlation as "no telluric" in a window where the test had no power
    to see any (the positive control is what separates those, and it is asserted here);
  * conflating a HOLDING's claim with a LINE's state, which is how kurucz2005 would
    otherwise be recorded as verified on the strength of lines no telluric band reaches;
  * letting the two RYA-1192 artifacts disagree silently about the H arm.
"""
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "data/results/rya1192/rya1192_sirius_verification.json"
PERLINE = ROOT / "data/results/rya1192/rya1192_sirius_per_line.csv"
WINDOWS = ROOT / "data/results/rya1192/rya1192_sirius_crires_windows.csv"
SRC = ROOT / "scripts/rya1192_sirius_verification.py"


@pytest.fixture(scope="module")
def doc():
    if not DOC.exists():
        pytest.skip("RYA-1192 Sirius artifact absent")
    return json.loads(DOC.read_text())


@pytest.fixture(scope="module")
def lines():
    if not PERLINE.exists():
        pytest.skip("RYA-1192 Sirius artifact absent")
    return pd.read_csv(PERLINE)


@pytest.fixture(scope="module")
def wins():
    if not WINDOWS.exists():
        pytest.skip("RYA-1192 Sirius artifact absent")
    return pd.read_csv(WINDOWS)


# ── what the Sirius leg exists to close ──────────────────────────────────────────────
def test_no_measured_line_is_left_unverifiable(lines):
    """🔴 THE POINT OF THE SIRIUS RUN. The Mac left 131 lines (88 kurucz2005 + 43 iag)
    recorded UNVERIFIABLE-HERE — a statement about that machine, not a verdict. Every one
    of them resolves here, and nothing may regress into that state."""
    bad = lines[lines.state.astype(str).str.startswith("UNVERIFIABLE")]
    assert bad.empty, f"{len(bad)} lines still unjudged: {bad.holding.value_counts().to_dict()}"
    assert len(lines) > 400, "the line inventory shrank — re-read before trusting the rates"


def test_iags_claim_is_verified_against_the_raw_sibling_nobody_had_used(lines, doc):
    """🔴 THE SIBLING NOBODY HAD USED. `solar_iag` (Baker+2020, corrected) and
    `solar_iag_reiners2016` (uncorrected) are the SAME IAG FTS atlas. Reiners is
    span-capped at 5001.1 A so SELECTION can never prefer it to the corrected atlas — a
    policy, not the file's extent — and reading through the cap is what makes the claim
    checkable at all. Verified in all five discriminating bands, O2 A-band 2.1% below 0.8
    against the raw sibling's 66.8%."""
    k = lines[lines.holding == "solar_iag"]
    assert len(k) == 43, f"IAG line count moved: {len(k)}"
    assert (k.raw_sibling == "solar_iag_reiners2016").all()
    claim = doc["holding_claim_verdicts"]["solar_iag"]["claim_verdict"]
    assert claim.startswith("VERIFIED-CORRECTED"), claim
    row = next(k for k, v in doc["band_verdicts"].items()
               if k.startswith("solar_iag|O2 A-band"))
    v = doc["band_verdicts"][row]
    assert v["verdict"] == "VERIFIED-CORRECTED" and v["anchor"] == "solar_iag_reiners2016"


def test_the_difference_test_is_not_assumed_to_be_the_stronger_one(doc):
    """🔴 THE MISTAKE THAT PUT IAG'S CLAIM ON THE WRONG SIDE. Baker and Reiners are
    separate reductions sampled 5.6x apart and differ by 0.083 in flux in telluric-free
    windows before any telluric question is asked. Ten times that floor is more absorption
    than the O2 A-band contains — so the difference test could not have returned CORRECTED
    there for ANY data, and duly returned RAW-LIKE on a holding that is plainly corrected.
    Which test applies is decided by comparing the bar with the ceiling, per band."""
    floors = doc["reduction_floors_from_control_windows"]
    assert floors["solar_kpno_molecfit_corrected"] == 0.0, (
        "KP's two products are one array — an exact zero is what makes its zero decisive")
    assert floors["solar_iag"] > 0.05, (
        "if IAG's floor has collapsed, the fallback reasoning must be re-derived")
    iag = {k: v for k, v in doc["band_verdicts"].items()
           if k.startswith("solar_iag|") and v["verdict"].startswith("VERIFIED")}
    assert iag and all(v["test"].startswith("DEPTH") for v in iag.values()), (
        "IAG must fall back to depth, and say why")
    assert all("no power on this pair" in v["fell_back_to_depth_because"]
               for v in iag.values())
    harps = {k: v for k, v in doc["band_verdicts"].items()
             if k.startswith("solar_harps") and "O2 B-band" in k}
    assert all(v["test"].startswith("DIFFERENCE") for v in harps.values()), (
        "HARPS's floor is 3e-5 — the difference test is the strong one there")


def test_kurucz2005_holding_claim_and_line_state_are_NOT_conflated(lines, doc):
    """🔴 TWO DIFFERENT QUESTIONS, AND ONLY ONE OF THEM IS ABOUT THESE LINES.

    Kurucz ships only the corrected residual atlas, so no raw sibling exists and the 1984
    KP atlas is a different observation, not one (RYA-929 'provenance-distinct'). The
    HOLDING's delivered-corrected claim is therefore judged at band level against the
    uncorrected anchor — and it holds. But not one of its 88 measured Fe lines sits in a
    declared telluric band, so that verdict is NOT evidence about them, and they must not
    borrow it."""
    k = lines[lines.holding == "solar_kpno_kurucz2005_corrected"]
    assert len(k) == 88, f"kurucz2005 line count moved: {len(k)}"
    assert (k.state == "NO-TELLURIC-BAND-DECLARED").all(), k.state.value_counts().to_dict()
    assert k.raw_sibling.isna().all(), "kurucz2005 has no raw sibling — none may appear"
    # the holding-level claim is separately recorded, and it is NOT a clean pass
    claim = doc["holding_claim_verdicts"]["solar_kpno_kurucz2005_corrected"]["claim_verdict"]
    assert claim.startswith("⚠️ MIXED"), claim
    assert (k.holding_claim == claim).all(), (
        "each line must carry the holding claim explicitly, so it is never inferred")


def test_delivered_corrected_is_not_uniform_across_bands(doc):
    """🔴 'DELIVERED TELLURIC-CORRECTED' IS NOT ONE STATE. Kurucz 2005 is corrected in the
    O2 B-band, O2 A-band and both far-red H2O complexes — and RAW-LIKE in H2O 7160-7340,
    where it retains 2.95% below 0.8 against the anchor's 12.35% (a ratio of 4.2, under
    the 5x bar). One holding label is describing two different states, which is the same
    shape of defect as the KP one the ticket opened on."""
    v = doc["band_verdicts"]["solar_kpno_kurucz2005_corrected|H2O|7160-7340"]
    assert v["verdict"] == "VERIFIED-RAW-LIKE", v
    assert v["ratio_anchor_over_holding"] < 5.0
    o = doc["holding_claim_verdicts"]["solar_kpno_kurucz2005_corrected"]
    assert len(o["corrected_bands"]) == 4 and len(o["raw_like_bands"]) == 1


def test_an_undeclared_band_is_not_an_absent_one(lines):
    """RYA-833. NO-TELLURIC-BAND-DECLARED says our registry declares nothing here — and
    that registry has already been caught incomplete once (the missing o2gamma, RYA-1193).
    The label and its basis must both say so rather than reading as 'clean'."""
    k = lines[lines.state == "NO-TELLURIC-BAND-DECLARED"]
    assert len(k), "expected the kurucz2005 lines here"
    assert k.basis.str.contains("RYA-833").all()
    assert not k.state.str.contains("CORRECTED").any()


def test_kp_molecfit_is_still_raw_on_every_measured_line(lines):
    """The Mac finding, REPRODUCED on the machine that holds everything: 146 of 146
    near-UV + VIS lines byte-identical to raw under a name that claims a correction."""
    k = lines[lines.holding == "solar_kpno_molecfit_corrected"]
    assert len(k) == 146, f"KP line count moved: {len(k)}"
    assert (k.state == "VERIFIED-RAW").all(), k.state.value_counts().to_dict()
    assert (k.max_abs_diff.fillna(0) == 0).all()


# ── the method, and the control that licenses it ─────────────────────────────────────
def test_the_null_is_measured_by_displacement_not_by_a_handful_of_clean_windows(doc, wins):
    """🔴 THE ESTIMATOR I HAD TO REPLACE. The first cut took the null to be max(corr_r)
    over the arm's clean windows. The H arm has TWO of those and both were slightly
    negative, so the threshold landed at -0.015 and flagged 12 of 18 telluric windows. A
    null from two samples is not a null. Every window now carries its own, from the same
    statistic on the same data with the template displaced beyond any physical residual."""
    for arm in doc["crires_by_arm"]:
        assert "|dv| >= 150 km/s" in arm["null_estimator"]
        assert "per window" in arm["null_estimator"]
    assert "NULL_SHIFTS_KMS" in SRC.read_text()
    m = wins[wins.state == "MEASURED"]
    assert m.null_sd.notna().all() and (m.null_sd > 0).all(), (
        "a window with a zero or missing null sd cannot yield a z-score")
    for col in ("raw_z", "corr_z", "ratio_z"):
        assert m[col].notna().all(), f"{col} missing on some measured window"


def test_a_quiet_product_is_only_evidence_where_the_raw_frame_is_loud(doc, wins):
    """🔴 THE GUARD THAT MAKES 'CLEAN' MEAN ANYTHING. corr_z near zero says 'no telluric'
    only where raw_z is large; where the template is flat there is nothing to detect and
    the same number means 'the test has no power here'. Power is asserted PER WINDOW —
    gating the whole arm on a median once sent 18 well-powered windows to UNDETERMINED
    because the arm median missed a chosen multiple by 0.1 sigma."""
    src = SRC.read_text()
    assert "POWER IS A PROPERTY OF EACH WINDOW" in src
    for arm in doc["crires_by_arm"]:
        if arm["verdict"].startswith(("VERIFIED", "⚠️ LOCALISED")):
            assert arm["n_strong_windows_with_power"] >= 1
            assert arm["verdict_detail"]["raw_z_median_in_strong"] >= 3.0, (
                f"{arm['arm']}: a verdict was reached without a passing positive control")


def test_the_ratio_is_the_positive_form_of_the_same_claim(doc):
    """Two statistics, opposite in sign, from independent parts of the data: the product
    must NOT correlate with the telluric template, and raw/product MUST. A correction that
    is really there shows up both ways, and only one of them depends on the raw frame's
    flattening — so agreement between them is not an assumption of the normaliser."""
    for arm in doc["crires_by_arm"]:
        d = arm["verdict_detail"]
        if not d or d.get("ratio_z_median_in_strong") is None:
            continue
        assert d["ratio_z_median_in_strong"] > d["corr_z_median_in_strong"], (
            f"{arm['arm']}: the removed transmission is not more significant than the "
            f"residual — the correction claim does not hold up both ways")


# ── CRIRES+, the priority ask ────────────────────────────────────────────────────────
def test_the_crires_verdict_is_per_arm_and_the_arms_differ(doc):
    """One CRIRES+ verdict averages an almost telluric-free arm with a contaminated one.
    The Y arm's template barely bites at all; the H arm is bracketed by CO2 and the 1.9um
    H2O wing. They must be reported separately."""
    arms = {a["arm"]: a for a in doc["crires_by_arm"]}
    assert set(arms) == {"Y", "H"}
    assert arms["Y"]["verdict"] != arms["H"]["verdict"]
    assert arms["Y"]["n_strong_windows"] < arms["H"]["n_strong_windows"], (
        "the Y arm is the clean one — if that inverts, re-read before trusting either")


def test_the_raw_comparand_is_the_SAME_OBSERVATIONS(doc):
    """🔴 WHAT MADE CHECK (a) POSSIBLE AT ALL, and it had been recorded as impossible.
    `normalize_vesta_ir.py` still says a sweep found ZERO Vesta FITS on Sirius, so molecfit
    'had nothing to run on'. Eighteen IDPs are staged there now, and Elgueta's own
    table1.dat dates the solar spectrum to 2022-11-22 00:23 — the same night and programme
    as those IDPs. This is a raw/corrected pair, not two reductions of two nights."""
    assert "SAME OBSERVATIONS" in doc["supersedes"]["now_judged_by"] or \
           "SAME OBSERVATIONS" in SRC.read_text()
    assert "60.A-9051" in SRC.read_text()
    assert doc["telluric_template_files"], "no telluric template — the leg cannot run"
    lo, hi = doc["telluric_template_span_A"]
    assert lo < 15007 and hi > 17277.5, (
        f"template {lo}-{hi} does not span the graded H window")


def test_the_H_arm_finding_is_localised_and_names_its_windows(doc):
    """Whatever the H arm's verdict, it must be stated window by window. 'The H arm is
    dirty' and 'the H arm is clean' are both wrong: 18 windows carry telluric and they do
    not behave alike."""
    h = next(a for a in doc["crires_by_arm"] if a["arm"] == "H")
    assert h["n_strong_windows"] >= 10
    if h["flagged_windows"]:
        for f in h["flagged_windows"]:
            assert f["corr_z"] >= 3.0 and f["raw_z"] >= 3.0
        assert len(h["flagged_windows"]) < h["n_strong_windows"], (
            "if every window flags, suspect the method before the data")
        assert "LOCALISED" in h["verdict"]


def test_the_superseded_finding_is_named_not_left_to_be_noticed(doc):
    """🔴 TWO ARTIFACTS, ONE TICKET. `rya1192_verification.json` is committed on main with
    a CRIRES+ H candidate residual reached from the catalogue-excess proxy. This leg
    re-judges it against raw. The stale artifact is the one a reader finds first, so the
    replacement has to say what it replaces."""
    s = doc["supersedes"]
    assert s["artifact"].endswith("rya1192_verification.json")
    for w in ("15700", "16000", "17300"):
        assert w in s["prior_finding"], f"the prior flagged window {w} is not named"
    assert s["H_arm_now"] and s["Y_arm_now"]
    assert "displaced" in s["now_judged_by"] or "MTRANS" in s["now_judged_by"]


def test_the_named_graded_line_separates_measurability_from_telluric_state(doc):
    """RYA-1094's Ruffoni line 16009.610 was flagged by the Mac leg as sitting in a
    candidate-residual window AND as unserved by any arm. Those are different problems and
    the artifact must not let one stand in for the other."""
    n = doc["named_line_checks"]["16009.610"]
    assert n["in_canonical_gf_LAB"] is True
    assert "DIFFERENT questions" in n["note"]
    assert "its_window" in n


# ── the band signature, and what makes it a calibration ──────────────────────────────
def test_the_signature_test_is_anchored_and_the_controls_must_AGREE(doc):
    """⚠️ WHAT SEPARATES 'TELLURIC REMOVED' FROM 'SPECTRUM FLATTENED'. A holding that is
    shallower than the anchor everywhere has lost solar structure, not tellurics. The
    control windows are where a real correction must MATCH the uncorrected anchor, and
    both corrected KP-class holdings do: 15.7% vs 16.9% below 0.8 at 5000-5100 A."""
    ctrl = [r for r in doc["band_signature_table"] if r["kind"] == "CONTROL"]
    assert len(ctrl) >= 2, "one control window cannot show agreement"
    for row in ctrl:
        a = row["holdings"]["solar_kpno"]
        if a.get("state") != "MEASURED":
            continue
        for hold in ("solar_kpno_kurucz2005_corrected", "solar_iag"):
            c = row["holdings"].get(hold, {})
            if c.get("state") != "MEASURED":
                continue
            assert c["pct_below_0.8"] > 0.5 * a["pct_below_0.8"], (
                f"{hold} is shallower than the uncorrected anchor in a CONTROL window "
                f"({c['pct_below_0.8']} vs {a['pct_below_0.8']}) — that is a flattened "
                f"spectrum, not a telluric correction")


def test_a_band_the_anchor_cannot_discriminate_says_so(doc):
    """🔴 THE GATE THAT KILLED FOUR FALSE VERDICTS. The O2 gamma-band came back
    VERIFIED-RAW-LIKE on all four holdings at once — because 6% of pixels below 0.8 at
    6270 A is ORDINARY SOLAR LINE DENSITY, and it only looked like telluric when compared
    with a global control 270 A away. A band now has to be deeper than its OWN SIDE
    WINDOWS before it is allowed to judge anything."""
    vs = doc["band_verdicts"]
    assert vs, "no band verdicts"
    allowed = {"VERIFIED-CORRECTED", "VERIFIED-RAW-LIKE",
               "ANCHOR-SHOWS-NO-LOCAL-TELLURIC", "NO-LOCAL-BASELINE",
               "NO-ANCHOR", "NOT-COVERED"}
    assert {v["verdict"] for v in vs.values()} <= allowed
    gamma = {k: v for k, v in vs.items() if "6270-6300" in k}
    assert gamma, "the o2gamma band must still be probed"
    assert not any(v["verdict"] == "VERIFIED-RAW-LIKE"
                   and v.get("anchor_is_own_raw_sibling") is False
                   for v in gamma.values()), (
        "a depth-based verdict in the O2 gamma-band is the artifact this gate removes")


def test_a_holding_with_a_raw_sibling_is_judged_by_DIFFERENCE_not_depth(doc):
    """🔴 A DEPTH RATIO ACROSS TWO INSTRUMENTS IS NOT A MEASURE OF A CORRECTION. HARPS and
    IAG each have their own uncorrected counterpart, so the question 'did anything change
    in this band' is answerable directly, and the answer needs no view on how much of the
    absorption was telluric. Only kurucz2005, which has no sibling at all, falls back to a
    depth comparison — and its rows say so."""
    vs = doc["band_verdicts"]
    for k, v in vs.items():
        hold = k.split("|")[0]
        if v["verdict"] not in ("VERIFIED-CORRECTED", "VERIFIED-RAW-LIKE"):
            continue
        if hold == "solar_kpno_kurucz2005_corrected":
            assert v["anchor_is_own_raw_sibling"] is False, (
                "kurucz2005 is the one holding with no sibling at all")
            assert v["anchor"] == "solar_kpno"
        else:
            assert v["anchor_is_own_raw_sibling"] is True, f"{k} used a foreign anchor"
        # whichever test answered, the row must name it and name its anchor
        assert v["test"].startswith(("DIFFERENCE", "DEPTH")), k
        if v["test"].startswith("DEPTH") and v["anchor_is_own_raw_sibling"]:
            assert "no power on this pair" in v["fell_back_to_depth_because"], k


def test_each_baseline_is_the_right_KIND_of_baseline(doc):
    """⚠️ TWO BASELINES, AND SWAPPING THEM BREAKS BOTH TESTS.

    The DEPTH test needs a LOCAL one: solar line density runs from 16.9% below 0.8 at
    5000 A to 3.0% at 6000 A, so a band compared with a control 270 A away is compared
    with a different line forest — which is how the O2 gamma-band passed for telluric.

    The DIFFERENCE test needs a GLOBAL one: the reduction floor is a property of the pair,
    and side windows sit in a telluric band's WINGS, which the correction also touches.
    Judged against side windows, HARPS's O2 B-band came back "raw-like" at a ratio of 3.8
    on a pair whose in-band correction is four orders above its floor."""
    src = SRC.read_text()
    assert "SIDE_PAD_FRACTION" in src and "SIDE_GAP_A" in src
    assert "Solar line density varies enormously with wavelength" in src
    assert "THE REDUCTION FLOOR IS A GLOBAL PROPERTY" in src
    for k, v in doc["band_verdicts"].items():
        t = v.get("test", "")
        if t.startswith("DIFFERENCE"):
            assert "reduction_floor_from_controls" in v, k
            assert "anchor_side_pct_below_0.8" not in v, k
        elif t.startswith("DEPTH"):
            assert "anchor_side_pct_below_0.8" in v, k


def test_the_iag_raw_sibling_is_read_through_a_SELECTION_cap_deliberately(doc):
    """The Reiners cap at 5001.1 A keeps measurement from preferring the uncorrected atlas.
    The audit reads past it on purpose and the code has to say why, or the next reader
    will think the cap was a bug and remove it."""
    src = SRC.read_text()
    assert "SELECTION" in src and "policy, not the file's extent" in src
    assert "load_iag_reiners_window" in src
    # and the raw sibling really did speak above the cap
    row = next(r for r in doc["band_signature_table"]
               if r["kind"] == "TELLURIC" and r["lo_A"] == 7594.0)
    r = row["holdings"]["solar_iag_reiners2016"]
    assert r["state"] == "MEASURED" and r["pct_below_0.8"] > 20, (
        "Reiners must show the O2 A-band it never corrected, or it is not the raw sibling")


# ── the discipline this ticket is about ──────────────────────────────────────────────
def test_it_verifies_and_never_corrects():
    """No correction applied, no value moved. Reading a product or a linelist is the whole
    job; WRITING anywhere but this ticket's own results directory is what must not happen,
    so the check is on the write calls and not on the mention of a path."""
    import ast as _ast
    tree = _ast.parse(SRC.read_text())
    dests = []
    for node in _ast.walk(tree):
        if not (isinstance(node, _ast.Call) and isinstance(node.func, _ast.Attribute)):
            continue
        if node.func.attr not in ("write_text", "to_csv", "write_bytes"):
            continue
        # the DESTINATION is the receiver for write_text and the first arg for to_csv
        target = (node.args[0] if node.func.attr == "to_csv" and node.args
                  else node.func.value)
        dests.append(_ast.unparse(target))
    assert dests, "no writes found at all — the artifact would be empty"
    for d in dests:
        assert d.startswith("out /"), f"a write escapes this ticket's results dir: {d}"
    assert "VERIFICATION ONLY" in SRC.read_text()
    # and the destination really is this ticket's own directory
    assert 'default="data/results/rya1192"' in SRC.read_text()


def test_the_registry_evidence_for_the_crires_claim_is_recorded_as_unsound():
    """🔴 THE CLAIM WAS TRUE AND ITS STATED EVIDENCE WAS NOT. The holdings registry backs
    `elgueta2026_vizier` telluric_applied=applied with '0.10% of window points below 0.5 in
    the O2 A-band vs 51.3% for Kitt Peak in the same window'. The O2 A-band is 7594-7685 A;
    the Elgueta Y arm starts at 9796.5 A. That statistic cannot have come from these
    spectra — it is RYA-794's Y science-window number relabelled. Verifying the claim does
    not repair the reasoning, and the script says so."""
    src = SRC.read_text()
    assert "7594-7685" in src and "9796.5" in src
    assert "relabelled" in src


# ── the coverage finding this leg had to correct ─────────────────────────────────────
def test_the_inter_order_gap_finding_was_a_threshold_artifact(doc):
    """🔴 THE CORRECTION IS BIGGER THAN THE FINDING. The artifact on main reports 12 graded
    Ruffoni lines "in CRIRES+ H inter-order gaps" and 16 of 29 unmeasurable. It asked
    `len(flux) >= 20` in a +/-0.5 A window of a product sampled at 16.3 points/A — so a
    NORMALLY SAMPLED line has ~16 points and fails. The threshold was selecting
    order-OVERLAP regions, not served lines."""
    g = doc["graded_ir_line_coverage"]
    assert g["n_graded_beyond_10796"] == 29
    assert g["n_true_gaps_zero_points"] == 0, (
        "a line with zero points would be a real gap — if one appears, the correction "
        "below is wrong and must be re-derived")
    assert g["n_that_the_prior_20_point_test_called_gaps"] == 12, (
        "the twelve lines the prior test rejected are the finding being corrected")
    assert g["n_served_by_a_live_product"] == 25
    assert g["n_outside_every_arm"] == 4
    assert "ARTIFACTS OF A >=20-POINT THRESHOLD" in g["correction"]


def test_the_H_product_has_almost_no_inter_order_gap_structure(doc):
    """The claim "inter-order gap" is checkable directly: a gap is a STEP in wavelength.
    The H product's largest step is 1.18 A and only two steps in the whole arm exceed
    0.5 A, so there is nearly nothing to fall into."""
    h = doc["graded_ir_line_coverage"]["product_sampling"]["H"]
    assert h["n_steps_over_0.5A"] <= 3, h["steps_over_0.5A"]
    assert h["points_per_A"] < 20, (
        "if sampling now exceeds the old 20-point threshold the artifact reasoning changes")


def test_the_corrected_count_matches_the_holdings_own_registered_count(doc):
    """RYA-1094 registered the H product as carrying 25 primary-lab Fe I lines. The audit
    said 13. A holding disagreeing with its own registry about its own contents is the
    signal that one of them is wrong — and here it was the audit."""
    reg = ROOT / "data/catalog/holdings_manifest_registry.csv"
    note = pd.read_csv(reg).set_index("holding_id").loc[
        "solar_crires_plus_h_rya1094", "notes"]
    assert "25 primary-lab Fe I lines" in str(note)
    assert doc["graded_ir_line_coverage"]["n_served_by_a_live_product"] == 25


def test_reduced_sampling_is_reported_as_precision_not_as_coverage(doc):
    """⚠️ Half-density sampling is still worse sampling. Correcting "gap" to "served" must
    not erase that, so every line carries its point count and its density relative to the
    arm mean."""
    for r in doc["graded_ir_line_coverage"]["lines"]:
        if r["arm"]:
            assert "n_points_pm0.5A" in r and "sampling_vs_arm_mean" in r
    thin = [r for r in doc["graded_ir_line_coverage"]["lines"]
            if r["arm"] and r["sampling_vs_arm_mean"] is not None
            and r["sampling_vs_arm_mean"] < 0.8]
    assert thin, "some lines really are thinly sampled — that fact must survive"

def test_not_one_measured_line_sits_where_a_correction_would_have_mattered(doc, lines):
    """🔴 THE SIZE OF THE KP DEFECT, IN SCIENCE TERMS. 146 lines are demonstrably on raw
    data under a name that claims a correction — and not one of the 462 measured Fe lines
    in this inventory sits inside a registered telluric band. That makes it a provenance
    defect rather than an abundance one, and the artifact has to say so in both directions
    rather than leaving a reader to infer either."""
    lim = doc["per_line_inventory_limits"]
    assert lim["n_inside_a_registered_telluric_band"] == 0
    assert not lines.correction_expected_here.any()
    assert "provenance defect, not an abundance one" in lim["finding"]


def test_the_line_inventory_states_what_it_does_not_cover(doc):
    """⚠️ AND THE STATEMENT ABOVE IS BOUNDED. The inventory is built from per-line
    artifacts on disk, and none exists for ANY red-optical, NIR or CRIRES+ Fe product —
    which is precisely where the registered telluric complexes live. Without that limit
    recorded, 'no measured line is in a telluric band' reads as a claim about all our Fe
    products, and it is not one (RYA-833)."""
    lim = doc["per_line_inventory_limits"]
    assert "red-optical" in lim["limit"] and "NIR" in lim["limit"]
    assert lim["wavelength_span_A"][1] < 6910, (
        "the inventory now reaches past VIS — the stated limit must be re-derived")


def test_lines_above_the_noise_floor_are_not_called_a_correction(doc):
    """⚠️ 28 HARPS lines clear their pair's reduction floor. The largest is 0.0017 in flux
    against 0.185 for the real O2 B-band correction — two orders down — and every one of
    them sits outside every telluric band. The artifact must carry that scale, or the
    count reads as 28 telluric corrections."""
    c = doc["per_line_corrected_but_outside_every_telluric_band"]
    assert c["n"] > 0
    assert c["max_abs_diff_max"] < 0.01, c["max_abs_diff_max"]
    assert "never as" in c["note"]


def test_a_window_verdict_is_carried_to_the_lines_that_sit_in_it(doc):
    """🔴 A WINDOW-LEVEL VERDICT ONLY MATTERS WHERE A LINE IS, and the answer moved with
    the windows. The prior artifact made this point about 16009.610; that line's window is
    now clear at −1.0σ, while seven graded Ruffoni lines in the H arm's blue end sit inside
    windows this leg flagged. Reporting flagged windows without saying which graded lines
    they contain leaves the consequence for someone else to work out."""
    g = doc["graded_lines_in_a_flagged_window"]
    h = next(a for a in doc["crires_by_arm"] if a["arm"] == "H")
    if not h["flagged_windows"]:
        assert g["n"] == 0
        return
    assert g["n"] > 0, "flagged windows exist — the lines inside them must be enumerated"
    for r in g["lines"]:
        assert r["window_A"][0] <= r["wavelength_air_A"] < r["window_A"][1]
        assert r["lab_source_tag"], "a graded line must name its lab source"
    assert g["n"] < g["of_n_graded_served"], (
        "if every served line were flagged this would be an arm verdict, not a localised one")
    assert "not a disqualification" in g["note"]


def test_the_named_line_16009_is_now_clear_and_the_artifact_says_which_way(doc):
    """The line the directive names by hand. Its window is the most strongly telluric in
    the RAW frame of the whole arm (raw_z 17.6) and the product sits BELOW its own null
    there — which is the strongest form the clearance can take, not the weakest."""
    n = doc["named_line_checks"]["16009.610"]
    w = n["its_window"]
    assert w["raw_z"] > 10, "the positive control must be loud for the clearance to mean anything"
    assert w["corr_z"] < 3.0, w
    assert not any(r["wavelength_air_A"] == 16009.61
                   for r in doc["graded_lines_in_a_flagged_window"]["lines"])

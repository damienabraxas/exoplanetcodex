"""RYA-1192 — the telluric label lies; these pin what verification actually found.

VERIFICATION ONLY. The fixes are separate (RYA-1191).

The tests are written around the three ways this map can be misread: treating
UNVERIFIABLE-HERE as a verdict, reading "0% corrected" as "wrong number" in every band
equally, and trusting a band registry that two sources disagree about.
"""
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "data/results/rya1192/rya1192_verification.json"
PERLINE = ROOT / "data/results/rya1192/rya1192_per_line.csv"
SRC = ROOT / "scripts/rya1192_telluric_verification.py"


@pytest.fixture(scope="module")
def doc():
    if not DOC.exists():
        pytest.skip("RYA-1192 artifact absent")
    return json.loads(DOC.read_text())


@pytest.fixture(scope="module")
def lines():
    if not PERLINE.exists():
        pytest.skip("RYA-1192 artifact absent")
    return pd.read_csv(PERLINE)


def test_every_kp_molecfit_measured_line_is_on_RAW_data(lines):
    """🔴 THE FINDING, and it is larger than the sample that prompted the ticket.
    RYA-1190 saw 7 of 10 sampled red-optical windows byte-identical to raw. At line level
    across near-UV + VIS it is EVERY ONE: 146 of 146, `max|corrected - raw| == 0`."""
    k = lines[lines.holding == "solar_kpno_molecfit_corrected"]
    assert len(k) > 100, "the KP line inventory shrank — re-read before trusting the rate"
    assert (k.state == "VERIFIED-RAW").all(), (
        f"some KP lines are now corrected: {k.state.value_counts().to_dict()}")
    assert (k.max_abs_diff.fillna(0) == 0).all()


def test_a_nonzero_difference_is_not_automatically_a_correction(lines, doc):
    """🔴 THE OVERSTATEMENT I MADE AND THEN CAUGHT. Judging on `max|diff| > 0` alone put
    184 of 185 HARPS lines in VERIFIED-CORRECTED. But HARPS's corrected and raw products
    are SEPARATELY REDUCED and differ by ~9e-5 EVERYWHERE as a reduction floor, while its
    real correction is FOUR ORDERS larger (mean |diff| 0.264 in the O2 B-band). A 1e-4
    difference at a measured line is the reduction, not a telluric correction.

    Each holding now carries its own baseline from control windows, and KP's comes back
    EXACTLY 0.0 — which is itself the proof that KP's two products are one array with an
    in-band correction, so a zero there is decisive in a way it would not be for HARPS."""
    b = doc["reduction_baselines"]
    assert b["solar_kpno_molecfit_corrected"] == 0.0, (
        "KP's corrected and raw are no longer the same array outside the fitted bands — "
        "the exact-zero verdict loses its force")
    assert b["solar_harps_molecfit_corrected"] > 1e-6, (
        "HARPS no longer shows a reduction floor — re-check before trusting per-line "
        "'corrected' verdicts on it")
    h = lines[lines.holding == "solar_harps_molecfit_corrected"]
    assert (h.state == "NO-MATERIAL-DIFFERENCE").sum() > 100, (
        "the HARPS lines are back to reading as corrected on reduction noise")
    assert (h.state == "VERIFIED-CORRECTED").sum() < 20


def test_the_correction_is_larger_where_telluric_actually_bites(lines):
    """⚠️ The two HARPS populations must stay SEPARATED by the materiality threshold. If a
    "materially corrected" line ever sits below a "reduction floor" one, the threshold is
    no longer distinguishing a correction from the reduction and every HARPS verdict is
    back in question."""
    h = lines[lines.holding == "solar_harps_molecfit_corrected"]
    mat = h[h.state == "VERIFIED-CORRECTED"].max_abs_diff
    noise = h[h.state == "NO-MATERIAL-DIFFERENCE"].max_abs_diff
    assert len(mat) and len(noise), "one of the two populations is gone"
    assert mat.min() > noise.max(), (
        "the materially-corrected and reduction-floor populations now overlap — the "
        "threshold has stopped separating a correction from the reduction")


def test_unverifiable_is_never_recorded_as_corrected_or_raw(lines, doc):
    """⚠️ MISREAD 1. KP2005, IAG and the CRIRES+ raw IDPs are staged on Sirius and cannot
    be read here. That is a limit of THIS MACHINE, not a verdict, and calling it either
    way would fabricate one (RYA-833)."""
    for hold in ("solar_iag", "solar_kpno_kurucz2005_corrected"):
        s = set(lines[lines.holding == hold].state)
        assert s == {"UNVERIFIABLE-HERE"}, f"{hold} acquired a verdict it cannot have: {s}"
    av = doc["availability_on_this_machine"]
    assert av["solar_crires_plus_y_wide_rya1054"]["probe"] == "UNVERIFIABLE-HERE"
    assert "Sirius" in SRC.read_text() and "UNVERIFIABLE-HERE" in SRC.read_text()


def test_the_two_band_registries_disagree_and_it_is_reported(doc):
    """🔴 The holding label says "six registered bands" and telluric_policy lists six;
    RYA-940 shipped SEVEN fit manifests. The extra is o2gamma 6270-6300 A — a band
    molecfit corrected and the policy does not know about. Reported, not resolved."""
    dis = doc["band_registry_disagreement"]
    assert [6270.0, 6300.0] in [list(x) for x in dis["fitted_not_in_policy"]]
    assert dis["policy_not_fitted"] == []
    assert len(doc["bands_molecfit_ACTUALLY_fitted"]) == 7
    assert len(doc["bands_the_policy_LISTS"]) == 6


def test_the_fit_pad_is_not_corrected(doc):
    """⚠️ `fit_window_A` is 25 A wider each side than `band_A`. Reading the fit window as
    the corrected span called four KP lines "inside a fitted band" while the comparison
    found them byte-identical to raw. The comparison was right and the membership test was
    wrong — the pad constrains the fit and is not itself corrected."""
    pad = doc["pad_is_not_corrected"]
    assert len(pad["lines_in_pad_only"]) >= 4
    for w in (6315.811, 6843.655, 6855.161, 6858.148):
        assert any(abs(x - w) < 0.01 for x in pad["lines_in_pad_only"]), w


def test_kp_band_coverage_is_reported_per_band_not_as_one_number(doc):
    """⚠️ MISREAD 2. The four bands are NOT equally affected, and a single headline would
    hide that. near-UV is 0% corrected but is also where the fitted molecules barely
    absorb; red-optical is 25% corrected and is where RYA-1190 measured a real depression
    in the untouched windows. Label defect everywhere, science impact unevenly."""
    cov = {(c["holding"], c["band"], tuple(c["window_A"])): c
           for c in doc["band_coverage_by_product"]}
    kp = {k: v for k, v in cov.items() if k[0] == "solar_kpno_molecfit_corrected"}
    assert kp, "the KP coverage rows are gone"
    fr = {k[1]: v["fraction_of_band_corrected"] for k, v in kp.items()}
    assert fr.get("near-UV") == 0.0, "near-UV is no longer 0% corrected"
    assert 0.0 < fr.get("red-optical", 0) < 0.5
    # and a holding with a global correction must NOT be given a fraction
    harps = [c for c in doc["band_coverage_by_product"]
             if c["holding"] == "solar_harps_molecfit_corrected"]
    assert all(c["fraction_of_band_corrected"] is None for c in harps)
    assert all("not band-limited" in c["basis"] for c in harps)


def test_the_crires_verdict_says_what_it_does_and_does_not_establish(doc):
    """The priority ask. Across the full 9800-10796 span the deep-absorption fraction
    tracks the KP1984 RAW atlas to ~0.003, so no telluric residual REMAINS. ⚠️ That is
    check (b): with the raw IDP sibling unreadable here, it is not proof the correction
    RAN, and in a window this clean the two would look alike either way."""
    v = doc["crires_verdict"]
    assert v.startswith("NO EVIDENCE OF UNCORRECTED TELLURIC")
    assert "NOT PROOF THE CORRECTION RAN" in v
    assert "OWED" in v and "Sirius" in v
    span = doc["crires_full_span"]
    assert len(span) >= 9, "the full span is no longer scanned"
    for w in span:
        a, k = w["crires_y_wide"], w["kp_raw_same_window"]
        if "frac_below_0.80" in a and "frac_below_0.80" in k:
            assert a["frac_below_0.80"] < 0.10, (
                f"CRIRES+ now shows a deep-absorption forest at {w['lo_A']} A")


def test_the_forest_statistic_is_calibrated_not_asserted(doc):
    """A bare 0.80 threshold would be a chosen number. It is calibrated on a pair where
    the answer is known — KP raw vs KP molecfit, same atlas — and the out-of-band CONTROL
    is what makes it a calibration: there the two are the same spectrum."""
    cal = {c["band"]: c for c in doc["forest_calibration"]}
    o2a = cal.get("o2a")
    assert o2a, "the O2 A-band calibration row is gone"
    assert o2a["kp_raw"]["frac_below_0.80"] > 5 * o2a["kp_molecfit"]["frac_below_0.80"], (
        "the statistic no longer separates a known-raw from a known-corrected spectrum")
    ctl = cal.get("CONTROL (no telluric band)")
    assert ctl and ctl["kp_raw"]["frac_below_0.80"] == ctl["kp_molecfit"]["frac_below_0.80"], (
        "the control must be byte-identical between raw and corrected — if it is not, the "
        "statistic is measuring something other than the correction")


def test_it_verifies_and_never_corrects():
    """VERIFICATION ONLY: it must not place a correction, and must write only its own
    results directory. Write sites are parsed, with the target node taken per method."""
    import ast as _ast
    src = SRC.read_text()
    for forbidden in ("molecfit_run", "apply_correction", "correct_exposures",
                      "crires_telluric", "telluric_intake"):
        assert forbidden not in src, f"this audit must not invoke {forbidden}"
    writes = []
    for node in _ast.walk(_ast.parse(src)):
        if not (isinstance(node, _ast.Call) and isinstance(node.func, _ast.Attribute)):
            continue
        if node.func.attr in {"to_csv", "to_json", "savefig"} and node.args:
            writes.append(_ast.unparse(node.args[0]))
        elif node.func.attr in {"write_text", "write_bytes"}:
            writes.append(_ast.unparse(node.func.value))
    assert writes
    for w in writes:
        assert "out_dir" in w, f"writes outside its own results directory: {w}"

"""RYA-1026: the ratified rules must be UNREPRESENTABLE to violate, not remembered."""
import ast
from pathlib import Path

import pytest

from pipeline.prenormalised_guard import (
    LOCALLY_NORMALISED, PRE_NORMALISED_HOLDINGS, RenormalisationError,
    assert_not_renormalising, is_pre_normalised)
from pipeline.telluric_display_policy import (
    REGISTRY_ANOMALIES, TelluricDisplayError, anomaly, assert_displayable, display_state)

KP1984 = "solar_kpno"
KP1984_CORR = "solar_kpno_molecfit_corrected"
KP2005 = "solar_kpno_kurucz2005_corrected"
KP_CLASS = (KP1984, KP1984_CORR, KP2005)
DELBOUILLE = "solar_delbouille_liege"

ROOT = Path(__file__).resolve().parent.parent


def _spec_flags() -> dict[str, bool]:
    """holding_id -> `pre_normalised`, parsed STATICALLY from `measure_band_ew.py`.

    Parsed rather than imported on purpose. That module resolves the Kitt Peak atlas at
    import time and `SystemExit`s when it is not staged, so importing it would make this
    test pass or fail on whether a data drive happens to be mounted -- an environment
    check wearing a policy check's label. The flags are literals in the source; reading
    them with `ast` needs no data at all.
    """
    tree = ast.parse((ROOT / "scripts" / "measure_band_ew.py").read_text())
    flags: dict[str, bool] = {}
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and getattr(node.func, "id", None)
                == "HoldingSpec"):
            continue
        holding_id = node.args[0].value
        for kw in node.keywords:
            if kw.arg == "pre_normalised":
                flags[holding_id] = kw.value.value
                break
        else:                                   # pragma: no cover - defended below
            raise AssertionError(
                f"{holding_id}: HoldingSpec declares no `pre_normalised`. The flag is "
                f"load-bearing (RYA-1026) and must never be left to a default.")
    assert flags, "parsed no HoldingSpecs -- the parser has drifted from the source"
    return flags


# ── the 2005 pre-normalised reversal (corrects RYA-929) ──────────────────────

def test_kp2005_is_registered_pre_normalised():
    """RYA-929 had this False. RYA-1026 reverses it: the product SHIPS NORMALISED."""
    assert is_pre_normalised(KP2005)
    assert KP2005 in PRE_NORMALISED_HOLDINGS


def test_the_WHOLE_kitt_peak_class_is_pre_normalised():
    """Ryan widened deliverable 5 from 'the 2005 file' to the whole KP class: it has now
    bitten twice, once per product (RYA-940 on 1984, RYA-929 on 2005)."""
    for h in KP_CLASS:
        assert is_pre_normalised(h), h


@pytest.mark.parametrize("holding", KP_CLASS)
@pytest.mark.parametrize("op", ["fitting_continuum", "applying_continuum"])
def test_touching_the_continuum_of_a_kp_atlas_RAISES(holding, op):
    """The corruption is SILENT -- the fit succeeds and the EW looks fine -- so the only
    effective guard is one that refuses. Fitting and applying are both refused: the
    ratified rule names the operation, not the mechanism."""
    with pytest.raises(RenormalisationError, match="SHIPS NORMALISED"):
        assert_not_renormalising(holding, pre_normalised=True, where="test",
                                 **{op: True})


def test_reading_without_touching_the_continuum_is_allowed():
    """Reading, or measuring the shipped continuum for a report, is a reading -- the line
    is whether the number we divide by is ours or the product's."""
    assert_not_renormalising(KP2005, pre_normalised=True)


def test_registry_vs_spec_disagreement_is_itself_the_bug():
    """RYA-929's flag said False while the product shipped normalised and NOTHING
    objected for months. Checking one signal would just re-bless whichever was written
    down; checking both makes the drift loud."""
    with pytest.raises(RenormalisationError, match="DRIFT"):
        assert_not_renormalising(KP2005, pre_normalised=False)


def test_registry_matches_every_holding_spec_flag():
    """The two signals only stay independent if they are also kept in step. A spec that
    drifts from the registry is caught at the call site by the DRIFT branch -- but only
    for holdings something actually measures. This catches the rest at import time."""
    specs = _spec_flags()
    for holding_id, flag in specs.items():
        assert (holding_id in PRE_NORMALISED_HOLDINGS) == flag, holding_id
    # Registry entries with no spec yet are allowed, but must be deliberate.
    registry_only = PRE_NORMALISED_HOLDINGS - set(specs)
    assert registry_only == {DELBOUILLE}, (
        f"unexpected registry-only entries {registry_only}: add the HoldingSpec, or "
        f"state here why the product is ratified ahead of being wired")


def test_repinning_a_locally_normalised_arm_to_unity_RAISES():
    """RYA-944: Delbouille is normalised per-window, max 0.9959 over 2.38 M points. A
    tidy-up to exactly 1.0 rescales every point by a number we chose."""
    assert DELBOUILLE in LOCALLY_NORMALISED
    with pytest.raises(RenormalisationError, match="pin the continuum to 1.0"):
        assert_not_renormalising(DELBOUILLE, pre_normalised=True, pinning_unity=True)


def test_unknown_holding_is_not_assumed_normalised():
    """Assuming it would apply unity as a continuum and inflate every EW silently.
    Absence of a claim is not a claim (RYA-833)."""
    assert not is_pre_normalised("some_new_arm")
    assert_not_renormalising("some_new_arm", pre_normalised=False,
                             fitting_continuum=True)      # allowed: nothing claimed


# ── telluric display gate ────────────────────────────────────────────────────

def test_corrected_products_render():
    for h in (KP1984_CORR, "solar_harps_molecfit_corrected", KP2005):
        assert_displayable(h)
        assert display_state(h) == "CLEAN"


def test_display_state_is_DERIVED_from_the_registry_not_hand_written():
    """A hand-written clean-set is a second declaration of one fact (RYA-845) and goes
    stale silently in the permissive direction. The first draft of this module did that
    and called `not-applied` Delbouille clean under an id that does not even exist."""
    from pipeline.telluric_policy import applied_state
    for h in (KP1984_CORR, KP2005, "solar_harps", DELBOUILLE):
        expected = "CLEAN" if applied_state(h) == "applied" else "BLOCKED"
        if h in REGISTRY_ANOMALIES:
            expected = "CLEAN_WITH_ANOMALY"
        assert display_state(h) == expected, h


def test_uncorrected_harps_and_delbouille_are_refused_as_science():
    """Both are registered `not-applied`. Passing `telluric_policy.gate_holding` for
    MEASUREMENT is not permission to ship -- RYA-1026's display rule is stricter."""
    for h in ("solar_harps", DELBOUILLE):
        with pytest.raises(TelluricDisplayError, match="telluric-CORRECTED input only"):
            assert_displayable(h)
        assert display_state(h) == "BLOCKED"


def test_uncorrected_1984_is_refused_as_science():
    with pytest.raises(TelluricDisplayError, match="REFUSING to render as science"):
        assert_displayable(KP1984)


def test_1984_is_allowed_as_the_named_control():
    """KP2005-vs-KP1984 IS the molecfit validation; deleting the uncorrected half would
    destroy the only pair that proves the correction works."""
    assert_displayable(KP1984, as_control=True)
    assert display_state(KP1984) == "CONTROL_ONLY"


def test_control_status_must_be_stated_not_inferred():
    """If a caller could get the waiver by accident the gate would be decorative."""
    with pytest.raises(TelluricDisplayError):
        assert_displayable(KP1984, as_control=False)


def test_a_disputed_applied_state_still_renders_but_carries_its_doubt():
    """RYA-944 found `iag_fts_solar_atlas` catalogued `corrected` while the manifest
    routes to the telluric-RETAINING Reiners file. No value moves under RYA-1026, so IAG
    still renders -- but the doubt must travel with it rather than be resolved silently
    in our favour."""
    assert display_state("solar_iag") == "CLEAN_WITH_ANOMALY"
    assert_displayable("solar_iag")
    assert "Reiners" in (anomaly("solar_iag") or "")
    assert anomaly(KP2005) is None


def test_unregistered_holding_is_distinguished_from_uncorrected():
    """'we hold this uncorrected' and 'we never wrote down what state this is in' are
    different problems with different fixes, so they must not collapse (RYA-833)."""
    with pytest.raises(TelluricDisplayError, match="UNKNOWN"):
        assert_displayable("some_arm_we_never_registered")
    assert display_state("some_arm_we_never_registered") == "UNREGISTERED"


# ── the UNDECLARED half: does the FLUX agree with the flag? ──────────────────
#
# Synthetic fixtures, deliberately. The real controls that FIXED these thresholds --
# KP1984's two columns at identical wavelengths, and KP2005's raw-vs-staged pair -- live
# on data drives that are not mounted everywhere, and a policy test that passes or fails
# on whether a drive is mounted is an environment check wearing a policy check's label.
# The measured numbers those controls produced are recorded beside the constants in
# `prenormalised_guard`; these fixtures reproduce each class's SHAPE.

import numpy as np

from pipeline.prenormalised_guard import (
    NORMALISED, UN_NORMALISED, UNKNOWN, NormalisationStateMismatch,
    UNITY_TOLERANCE, assert_data_matches_declaration, detect_normalisation_state)

_RNG = np.random.default_rng(1026)


def _spectrum(level, slope=0.0, n=20000, lo=5000.0, hi=5100.0):
    """An absorption spectrum whose upper envelope sits at `level` and tilts by `slope`."""
    w = np.linspace(lo, hi, n)
    # Centred on the band, so a tilt does NOT also move the median level. Otherwise the
    # two tests could never be exercised independently: a slope big enough to fail the
    # flatness bound would drag the level out of tolerance too, and the level test would
    # take the credit.
    x = (w - lo) / (hi - lo) - 0.5
    envelope = level * (1.0 + slope * x)
    depth = np.where(_RNG.random(n) < 0.04, _RNG.uniform(0.1, 0.8, n), 0.0)
    return w, envelope * (1.0 - depth)


def test_detects_a_normalised_product():
    """KP1984 col1 measured 0.9800 and KP2005's staged TSV 0.9997."""
    ev = detect_normalisation_state(*_spectrum(1.0))
    assert ev.value == NORMALISED, ev.citation()


def test_detects_absolute_flux():
    """KP1984 col2 measured a median rolling-P95 of 213.2875."""
    ev = detect_normalisation_state(*_spectrum(213.0, slope=0.02))
    assert ev.value == UN_NORMALISED


def test_the_LEVEL_test_catches_what_the_slope_test_cannot():
    """KP2005 `irradthu` as served sits at 2.14 but its slope over 100 A is only -0.044
    -- INSIDE the flatness bound. The level test is what catches it."""
    ev = detect_normalisation_state(*_spectrum(2.14, slope=-0.02))
    assert abs(ev.slope_across_band) <= 0.10          # the slope test would pass it
    assert ev.value == UN_NORMALISED                  # the level test does not


def test_the_SLOPE_test_catches_what_the_level_test_cannot():
    """A blaze or SED riding on a product that happens to sit near unity. Neither test
    alone is sufficient, which is why both are required."""
    ev = detect_normalisation_state(*_spectrum(1.0, slope=0.30))
    assert abs(ev.median_p95 - 1.0) <= UNITY_TOLERANCE   # the level test would pass it
    assert ev.value == UNKNOWN                           # not waved through
    assert "TILTED" in ev.reason


def test_too_few_windows_is_UNKNOWN_never_a_default():
    """`unknown` is a real answer, not a failure to try -- the third value is
    load-bearing exactly as it is in `telluric_intake` (RYA-833)."""
    ev = detect_normalisation_state(*_spectrum(1.0, n=400, lo=5000.0, hi=5004.0))
    assert ev.value == UNKNOWN
    assert "UNKNOWN" in ev.reason


def test_agreement_passes_and_returns_its_evidence():
    ev = assert_data_matches_declaration("x", *_spectrum(1.0), declared=True)
    assert ev.value == NORMALISED
    assert "median rolling-P95" in ev.citation()


def test_flux_disagreeing_with_the_flag_is_a_LOUD_STOP():
    """A declared flag and a MIS-ROUTED file agree with each other perfectly and are both
    wrong -- the `iag_fts_solar_atlas` shape (RYA-944), and what RYA-1026 found on the
    KP2005 route. The detector informs the declared state; it must never auto-fix."""
    with pytest.raises(NormalisationStateMismatch, match="THE FLUX AND THE FLAG DISAGREE"):
        assert_data_matches_declaration("solar_kpno_kurucz2005_corrected",
                                        *_spectrum(2.14), declared=True,
                                        where="test")


def test_unknown_does_not_raise():
    """Turning 'the data could not speak' into a hard failure would make every short
    window a blocker."""
    ev = assert_data_matches_declaration("x", *_spectrum(1.0, n=400, lo=5000.0, hi=5004.0),
                                         declared=True)
    assert ev.value == UNKNOWN


# ── the gate must be CALLED, not merely defined ──────────────────────────────

def test_the_display_gate_is_wired_into_the_render_path():
    """🔴 A GUARD NOBODY INVOKES IS DECORATIVE — this module's own docstring says so, and
    for a while that was exactly what this was: `assert_displayable` existed and was
    called only by these tests.

    The tracker generator is where a product becomes visible, so it is where the gate
    belongs. This pins the wiring rather than the wording, so the tracker cannot quietly
    go back to rendering products with no telluric state attached.
    """
    src = (ROOT / "scripts" / "rya935_live_status.py").read_text()
    assert "telluric_display_policy" in src, "the render path does not import the gate"
    assert "display_state(" in src, "the render path never calls display_state"
    assert "anomaly(" in src, "the render path drops the anomaly text"


def test_the_render_path_CARRIES_the_state_rather_than_raising():
    """Deliberate, and the reasoning must not be lost.

    Most existing Fe products were built on the telluric-UNCORRECTED 1984 atlas -- RYA-1026
    says so plainly -- so calling `assert_displayable` in the tracker would kill the whole
    page over historical products the ticket already knows about. Worse, DROPPING them
    would make a BLOCKED product indistinguishable from an ABSENT one, which is the
    absence-as-conclusion error the same ticket forbids (RYA-833).

    The raise belongs where a NEW product is published; the tracker labels.
    """
    import ast as _ast
    tree = _ast.parse((ROOT / "scripts" / "rya935_live_status.py").read_text())
    called = {getattr(n.func, "id", None) for n in _ast.walk(tree)
              if isinstance(n, _ast.Call)}
    # Parsed, not grepped: the module explains WHY it does not raise, so the name appears
    # in prose. A substring check would fire on the explanation and force the reasoning
    # to be deleted to keep the test green -- which is how a comment that matters gets
    # removed by a test that does not.
    assert "assert_displayable" not in called, (
        "the tracker must LABEL, not raise: one legacy product would take the page down "
        "and a dropped row reads as 'no data'")
    assert "display_state" in called

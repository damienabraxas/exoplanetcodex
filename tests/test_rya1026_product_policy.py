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

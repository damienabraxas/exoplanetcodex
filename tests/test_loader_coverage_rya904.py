"""
tests/test_loader_coverage_rya904.py
====================================
RYA-904 — the CRIRES+ telluric-CORRECTED holdings were UNADDRESSABLE, and the two keys
that made them so.

WHAT THE DEFECT ACTUALLY WAS
----------------------------
Not a telluric decision. `gate_holding()` was correct and is unchanged: it refuses
`solar_vesta_crires_plus_idp` (telluric_applied=not-applied on an instrument registered
telluric_required=yes) and passes both corrected holdings. What was wrong was the
DISPATCH one level up — `_LOADER_HOLDING` mapped one INSTRUMENT to one HOLDING, and for
`crires_plus` it named the raw IDPs. So the harness asked for the holding that should be
refused, was refused, and the two corrected holdings could not be named by anything.

An unreachable holding leaves no trace. A refusal leaves a reason; unreachability looks
exactly like having no data, which is why this survived from RYA-806 to RYA-904 with the
telluric gate green the whole time (RYA-833: an absence is a hypothesis, not a result).

WHAT THESE TESTS PIN
--------------------
  * both keys are HOLDING-level, and they moved TOGETHER. Wiring the corrected Y arm
    while leaving `PRE_NORMALISED` instrument-keyed sets a continuum on an
    already-normalised spectrum — the RYA-713 defect, EWs low by a median 11.7 %.
  * selection prefers a passing holding and NEVER falls back to a refused one.
  * the provenance NAMES the holding, because `crires_plus` alone does not distinguish
    the corrected Y arm from the raw IDPs.
  * the guard is RED on the pre-RYA-904 wiring. A guard nobody has seen fail is not a
    guard (RYA-805/833), so the control is executed here, not asserted in prose.

Everything here runs off the git-tracked RYA-794 product and the registry. No Sirius-only
data is needed, which is the point: the wiring is testable even where the spectra are not.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import measure_band_ew as M  # noqa: E402
from pipeline import loader_coverage as LC  # noqa: E402
from pipeline.telluric_policy import gate_holding  # noqa: E402

Y_HOLDING = "solar_crires_plus_y_rya794"
IDP_HOLDING = "solar_vesta_crires_plus_idp"
#: fully inside the RYA-794 product's 10280-10680 A span
IN_Y = 10400.0


# ── the two keys, and that they are the SAME key ─────────────────────────────

def test_both_maps_are_holding_level_and_agree_on_their_keys():
    """`_LOADER_HOLDING` and `PRE_NORMALISED` were separately instrument-keyed dicts.

    They are now one table. `PRE_NORMALISED` is DERIVED from `_INSTRUMENT_HOLDINGS`
    rather than written a second time, so the pair cannot come apart the way RYA-845's
    two declarations of one constant did.
    """
    assert not hasattr(M, "_LOADER_HOLDING"), (
        "the instrument->holding map is the defect; it must not survive alongside "
        "its replacement, or a caller can still reach the wrong one")
    wired = {h.holding_id for specs in M._INSTRUMENT_HOLDINGS.values() for h in specs}
    assert set(M.PRE_NORMALISED) == wired
    # every wired holding is a real registered holding, not a name someone typed
    import pandas as pd
    reg = set(pd.read_csv(LC.HOLDINGS).holding_id.astype(str))
    assert wired <= reg, f"unregistered holding(s) wired: {sorted(wired - reg)}"


def test_one_instrument_two_holdings_two_normalisation_answers():
    """The reason a per-instrument flag could not be right: CRIRES+ disagrees with itself.

    Raw Vesta IDPs are un-normalised adu (TUNIT2='adu', flux to ~1.9e5); the RYA-794
    corrected Y arm is already continuum-normalised (median 0.9985). A single answer for
    `crires_plus` is wrong for one of them, and being wrong in the
    "renormalise something already normalised" direction is RYA-713.
    """
    specs = {h.holding_id: h for h in M.holdings_for("crires_plus")}
    assert specs[Y_HOLDING].pre_normalised is True
    assert specs[IDP_HOLDING].pre_normalised is False


# ── selection: prefer passing, never fall back to refused ────────────────────

def test_the_gate_itself_is_unchanged_and_still_discriminates():
    """POSITIVE CONTROL on the gate before blaming the dispatch (RYA-896 verified this;
    re-run here so 'the gate is fine' is measured in this file, not quoted)."""
    assert gate_holding(Y_HOLDING, "crires_plus")[0] is True
    assert gate_holding(IDP_HOLDING, "crires_plus")[0] is False
    assert gate_holding("elgueta2026_vizier", "crires_plus")[0] is True


def test_select_prefers_the_corrected_holding_inside_its_span():
    spec = M.select_holding("crires_plus", IN_Y, 1.2)
    assert spec.holding_id == Y_HOLDING
    assert spec.pre_normalised is True


def test_outside_the_corrected_span_it_refuses_rather_than_serving_a_truncated_window():
    """Coverage must be TOTAL. A window half inside a fixed-span product is a truncated
    window, and serving it quietly would be worse than the refusal it replaces."""
    with pytest.raises(M.TelluricNotCorrected):
        M.select_holding("crires_plus", 10679.5, 1.2)     # straddles the red edge
    with pytest.raises(M.TelluricNotCorrected):
        M.select_holding("crires_plus", 10050.0, 1.0)     # blueward of the product


def test_a_refused_holding_is_never_silently_fallen_back_to():
    """Named explicitly, the raw IDP holding still cannot be measured — the gate decides,
    and `holding=` is a way to ASK, not a way to bypass."""
    with pytest.raises(M.TelluricNotCorrected) as e:
        M.select_holding("crires_plus", IN_Y, 1.2, holding=IDP_HOLDING)
    assert IDP_HOLDING in str(e.value)
    assert "not-applied" in str(e.value)


def test_unknown_instrument_is_loud():
    with pytest.raises(LookupError, match="no window loader for instrument"):
        M.select_holding("jwst_nirspec_prism", 10400.0, 1.0)


# ── the window itself ────────────────────────────────────────────────────────

def test_load_window_serves_the_corrected_arm_and_names_the_holding():
    win = M.load_window_ex("crires_plus", IN_Y, 1.2)
    assert win.holding.holding_id == Y_HOLDING
    assert f"holding={Y_HOLDING}" in win.provenance
    assert win.pre_normalised is True
    assert win.wave.size == win.flux.size > 12
    assert np.all(np.diff(win.wave) >= 0)
    # the corrected product is NORMALISED — assert the relationship, not a magic constant
    assert 0.9 < float(np.median(win.flux)) < 1.1
    # ... and the raw IDPs are not, which is the whole reason this is per-holding
    assert M.PRE_NORMALISED[IDP_HOLDING] is False


def test_the_gdsat_caveat_travels_with_the_holding():
    """RYA-904 spec 6. Elgueta's G-dwarf robustness flag certifies ZERO Fe I lines for a
    solar-type star, and RYA-794 quoted no abundance. Reaching the pool does not
    pre-decide that it is publishable, so the caveat is a property of the holding rather
    than a sentence someone remembers to retype in each consumer."""
    spec = M.select_holding("crires_plus", IN_Y, 1.2)
    assert "GDSat" in spec.caveat and "RYA-794" in spec.caveat


def test_the_product_is_air_and_rest_frame_and_that_is_checked_not_assumed():
    """CRIRES+ delivers VACUUM/TOPOCENT natively; `wavelength_air_A` in a DERIVED product
    is a claim. A vacuum copy has the same columns, the same row count and the same flux
    — it just fits every line ~2.9 A away, at the wrong abundance. So the loader runs a
    POSITIVE control on RYA-794's own published Fe I positions.

    THE TEST MUST DISCRIMINATE (RYA-805): shifting the product by the air-vacuum
    difference must make the check fail, or it is not testing anything.
    """
    M._crires_cache.pop("y", None)
    w, f = M.crires_y_spectrum()               # passes on the real product
    with pytest.raises(LookupError, match="AIR wavelengths"):
        M._assert_air_rest_frame(w + 2.9, f)   # a vacuum-scale shift is caught
    M._crires_cache.pop("y", None)


# ── the guard, and the control that proves it fires ──────────────────────────

def test_guard_is_green_on_the_wiring_as_shipped():
    LC.reconcile_loader_coverage()


def test_POSITIVE_CONTROL_guard_is_red_on_the_pre_rya904_wiring():
    """🔴 THE CONTROL. Point `crires_plus` back at ONLY the raw IDP holding — the exact
    `_LOADER_HOLDING` this ticket replaces — and clear the declared gaps, and the guard
    must name BOTH corrected holdings as unreachable.

    This is the state `origin/main` was in, with every telluric test green.
    """
    pre_904 = {"solar_kpno": "kpno_solar_atlas",
               "solar_iag": "iag_fts_solar_atlas",
               IDP_HOLDING: "crires_plus"}
    with pytest.raises(LC.LoaderCoverageError) as e:
        LC.reconcile_loader_coverage(addressable=pre_904, declared_gaps={})
    msg = str(e.value)
    assert Y_HOLDING in msg
    assert "elgueta2026_vizier" in msg
    # and it must NOT accuse the holding that is correctly refused
    assert IDP_HOLDING not in msg


def test_guard_does_not_fire_on_a_correctly_refused_holding():
    """DISCRIMINATION. `solar_vesta_crires_plus_idp` is verified and unaddressable-by-
    measurement, and it must NEVER appear: it does not pass the gate, so its absence from
    the harness is a decision with a reason, not a gap. A guard that also flagged it
    would be pressure to weaken the telluric gate — the opposite of the fix."""
    rows = [b["holding_id"] for b in LC.unreachable_holdings(
        addressable={}, declared_gaps={})]
    assert IDP_HOLDING not in rows
    assert Y_HOLDING in rows                    # ... while the corrected one does appear


def test_every_declared_gap_names_a_reason_and_a_ticket():
    """A gap is a CLAIM, and an unexplained one is indistinguishable from an oversight —
    which is the failure this whole ticket is about. The table can only shrink by
    building something."""
    assert LC.DECLARED_GAPS, "an empty gap table means nothing was declared, not that "\
                             "nothing is missing"
    import re
    for hid, why in LC.DECLARED_GAPS.items():
        assert len(why) > 80, f"{hid}: a gap needs a reason, not a label"
        assert re.search(r"RYA-\d+", why), f"{hid}: name the ticket that owns the gap"


def test_declared_gaps_are_real_holdings_not_typos():
    import pandas as pd
    reg = set(pd.read_csv(LC.HOLDINGS).holding_id.astype(str))
    assert set(LC.DECLARED_GAPS) <= reg


def test_a_gap_that_is_actually_wired_is_itself_a_failure():
    """The table rots BOTH ways. RYA-897 is wiring the HARPS solar arm as this ships;
    the moment it lands, `solar_harps`'s entry becomes an untrue statement about the code
    — and a stale exemption is how a guard quietly stops guarding, because it would wave
    through the next holding to go unreachable."""
    with pytest.raises(LC.LoaderCoverageError, match="are in fact WIRED"):
        LC.reconcile_loader_coverage(
            addressable={**LC.addressable_holdings(), "solar_harps": "harps"})

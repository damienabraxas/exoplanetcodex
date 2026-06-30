"""RYA-494 — tests for the generalized per-star/per-instrument IR telluric module.

Fast tests only: no molecfit run (that is the slow, engine-gated path verified by
scripts/condition_acen_a_ir_rya494.py). These lock the GENERALIZATION contract —
star/instrument-agnostic structure, the branched velocity step, the preserved
permanent IR rules, and the honest α Cen A coverage findings.
"""
import numpy as np
import pytest

from pipeline import ir_telluric as irt
from pipeline.crires_telluric import NM_TO_A, TelluricNotCorrectedError


# ── target registry + velocity-step branch ────────────────────────────────────
def test_velocity_mode_branch_is_per_target():
    assert irt.TARGETS["vesta"].velocity_mode is irt.VelocityMode.REFLECTED_SOLAR
    assert irt.TARGETS["alpha_cen_a"].velocity_mode is irt.VelocityMode.STELLAR
    # α Cen A is a direct star → NOT asteroid ephemeris
    assert irt.TARGETS["alpha_cen_a"].body_key is None


def test_object_normalizer_matches_rya479():
    assert irt.re_norm("HD128620") == "alpha_cen_a"
    assert irt.re_norm("alf Cen A") == "alpha_cen_a"
    assert irt.re_norm("HD128621") == "alpha_cen_b"
    assert irt.re_norm("Star S5") == "other"


# ── permanent IR rules preserved ──────────────────────────────────────────────
def _fake_nirps(telluric_corrected=False):
    f = irt.NirpsFrame(
        path=type("P", (), {"name": "fake.fits"})(), object_hdr="AlphaCenB",
        mjd=60000.0, date_obs="2023", specsys="BARYCENT", berv_kms=15.0,
        wave_A=np.linspace(9700, 9800, 200), flux=np.ones(200),
        err=np.ones(200) * 0.01, atm_transm=np.linspace(0.5, 1.0, 200))
    f.telluric_corrected = telluric_corrected
    return f


def test_rv_condition_refuses_non_telluric_corrected():
    """Permanent rule: no RV shift before telluric correction is verified."""
    t = irt.TARGETS["alpha_cen_a"]
    with pytest.raises(TelluricNotCorrectedError):
        irt._stellar_rv_condition(_fake_nirps(telluric_corrected=False), t)


def test_nirps_telluric_correct_sets_flag_and_needs_real_model():
    """NIRPS 'correction' = real ATM_TRANSM present; an all-ones model is refused."""
    ok = _fake_nirps()
    irt.telluric_correct_nirps(ok)
    assert ok.telluric_corrected is True
    flat = _fake_nirps(); flat.atm_transm = np.ones(200)
    with pytest.raises(Exception):
        irt.telluric_correct_nirps(flat)


def test_stellar_branch_nirps_applies_only_systemic_rv():
    """NIRPS is already BARYCENT → only the systemic RV is removed (no double-BERV)."""
    t = irt.TARGETS["alpha_cen_a"]
    f = _fake_nirps(telluric_corrected=True)
    w0 = f.wave_A.copy()
    irt._stellar_rv_condition(f, t)
    assert f.rest_frame is True
    assert f._rv_applied["berv"] == 0.0                 # not re-applied (already BARYCENT)
    assert f._rv_applied["systemic"] == t.systemic_rv_kms
    # blueshift for negative systemic RV
    assert np.allclose(f.wave_A, w0 / (1 + t.systemic_rv_kms / irt._C_KMS))


# ── honest α Cen A findings (coverage + attribution) ──────────────────────────
def test_oi_844_926_below_ir_blue_edges():
    """The brief's headline O I 844/926 nm fall below CRIRES-Y (949.6) & NIRPS (966)."""
    for wl_A in (8446.5, 9265.9):
        assert wl_A < 9496, "O I line is below the CRIRES-Y blue edge"
        assert wl_A < 9661, "O I line is below the NIRPS blue edge"


def test_nirps_attribution_is_starid_not_object_header():
    """RYA-423 IR RV star-ID (verdict=A) is authoritative over the mislabeled OBJECT."""
    note = irt.ALPHA_CEN_A_NIRPS_NOTE
    assert "verdict=A" in note and "RYA-423" in note
    assert "INVERTS the RYA-479" in note                # the cross-ticket correction


def test_instruments_enumerated():
    assert {i.value for i in irt.Instrument} == {"CRIRES", "NIRPS"}

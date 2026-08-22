"""`SYNTH_BANDS` has ONE home, and the lift changed no number — RYA-967.

The four synthesis windows were declared inside `scripts/derive_band_products.py`, and
`scripts/rya855_rung_audit.py` imported them FROM THAT SCRIPT. A config constant whose
only home is an executable is the second-home shape this project keeps paying for
(RYA-350/353/954), with a sharper cost here: the driver's import chain loads the Kitt Peak
atlas, so reading a band's half-width required an atlas to be on disk.

Two things are pinned:

1. **The lift is value-preserving.** The three pre-existing entries are asserted against
   the literals that were in the script. A move that silently rounded `1.10` to `1.1` and
   re-derived it from the invariant would move the published near-UV 7.487.
2. **The VIS entry is derived, not invented.** `synthesis_route` refused to invent a
   half-width for an uncharacterised regime; the test asserts the new value sits on the
   same invariant the other three do, so it is that refusal being satisfied rather than
   bypassed.
"""
from __future__ import annotations

import pytest

from config.synth_bands import (SYNTH_BANDS, HALF_WIDTH_IN_DOPPLER_SIGMA,
                                STELLAR_SIGMA_KMS, derive_half_width_A, ISPEC_GES_V6)

#: Verbatim from `scripts/derive_band_products.py` before the lift (git 9167398).
PRE_LIFT = {
    "near-UV":     dict(half_width_A=0.40, min_sep_A=4.0, n_lines=40),
    "red-optical": dict(half_width_A=1.10, min_sep_A=4.0, n_lines=40),
    "NIR":         dict(half_width_A=1.40, min_sep_A=4.0, n_lines=40),
}


@pytest.mark.parametrize("band", sorted(PRE_LIFT))
def test_the_lift_changed_no_pre_existing_value(band):
    cfg = SYNTH_BANDS[band]
    for field, want in PRE_LIFT[band].items():
        got = getattr(cfg, field)
        assert got == want, (
            f"{band}.{field} is {got}, was {want} before the RYA-967 lift. Moving "
            f"SYNTH_BANDS to config must not change a number — the near-UV entry keys "
            f"the published 7.487 (RYA-832).")


def test_vis_was_added():
    assert "VIS" in SYNTH_BANDS
    vis = SYNTH_BANDS["VIS"]
    assert (vis.lo_A, vis.hi_A) == (3780.0, 6910.0)
    assert vis.linelist_spec == ISPEC_GES_V6, (
        "Part B compares EW-integration against profile-synthesis with METHOD as the only "
        "variable, so the synth leg must read the same list the EW route inverts against.")


def test_every_band_sits_on_the_shared_invariant():
    """The rule, asserted as a RELATIONSHIP (RYA-870).

    Pinning `VIS.half_width_A == 0.62` would pin the example. What makes 0.62 defensible
    is that it is the same rule the other three already obeyed, so that is what is tested —
    and it stays true if a future band is added or `stellar_sigma_kms` is re-ratified.
    """
    for name, cfg in SYNTH_BANDS.items():
        k = cfg.half_width_in_doppler_sigma
        assert abs(k - HALF_WIDTH_IN_DOPPLER_SIGMA) <= 0.35, (
            f"{name} sits at {k:.2f} Doppler sigma against the shared "
            f"{HALF_WIDTH_IN_DOPPLER_SIGMA:.2f}. Either the entry is wrong or the "
            f"invariant no longer describes the set — say which, in the YAML.")


def test_the_invariant_is_evaluated_at_the_anchor_not_the_band_centre():
    """🔴 THIS DISTINCTION IS WHY THE SET LOOKED BROKEN ON THE FIRST PASS.

    Two anchors are NOT band midpoints — red-optical spans 6910-9199 (centre 8055) and is
    anchored at 9600. Computing sigma_D from the centre reports 24.08 for it and makes a
    consistent set look inconsistent, which is exactly the kind of wrong diagnosis that
    gets a correct value 'fixed'.
    """
    ro = SYNTH_BANDS["red-optical"]
    centre = 0.5 * (ro.lo_A + ro.hi_A)
    assert ro.anchor_A != pytest.approx(centre), "the premise of this test has changed"
    from_centre = ro.half_width_A / (centre * STELLAR_SIGMA_KMS / 299792.458)
    assert abs(from_centre - HALF_WIDTH_IN_DOPPLER_SIGMA) > 1.0
    assert abs(ro.half_width_in_doppler_sigma - HALF_WIDTH_IN_DOPPLER_SIGMA) <= 0.35


def test_vis_half_width_is_what_the_invariant_gives():
    assert derive_half_width_A(SYNTH_BANDS["VIS"].anchor_A) == pytest.approx(0.62, abs=0.005)


def test_reading_a_half_width_does_not_import_an_engine():
    """The lift's whole point. `config.synth_bands` must be importable on a machine with
    no iSpec — the VIS entry's list resolves lazily, on `.linelist` access only."""
    import subprocess, sys, pathlib
    root = pathlib.Path(__file__).resolve().parent.parent
    r = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.modules['ispec']=None\n"
         "from config.synth_bands import SYNTH_BANDS\n"
         "assert SYNTH_BANDS['VIS'].half_width_A == 0.62\n"
         "import sys as s; assert s.modules.get('pipeline.abundances_derive') is None\n"
         "print('ok')"],
        cwd=root, capture_output=True, text=True)
    assert r.returncode == 0 and "ok" in r.stdout, (
        f"importing config.synth_bands pulled in an engine.\n{r.stdout}\n{r.stderr}")


def test_the_driver_no_longer_declares_the_bands():
    """The old home must be gone, not merely unused — an unused duplicate is still a
    duplicate, and the next edit lands in whichever one the editor opened first."""
    src = (__import__("pathlib").Path(__file__).resolve().parent.parent
           / "scripts" / "derive_band_products.py").read_text(encoding="utf-8")
    assert "SYNTH_BANDS: dict[str, SynthBand] = {" not in src
    assert "class SynthBand:" not in src
    assert "from config.synth_bands import" in src


def test_legacy_nearuv_constants_derive_from_the_config():
    """`NEARUV_HALF_WIDTH_A` and friends were literals in the driver AND values in the
    near-UV entry. After the lift the literals would have been a second declaration of the
    same three numbers, so they now read from the config."""
    sys_path_guard = __import__("sys")
    import scripts.derive_band_products as dbp
    assert dbp.NEARUV_HALF_WIDTH_A is SYNTH_BANDS["near-UV"].half_width_A
    assert dbp.NEARUV_MIN_SEP_A is SYNTH_BANDS["near-UV"].min_sep_A
    assert dbp.NEARUV_N_LINES is SYNTH_BANDS["near-UV"].n_lines
    del sys_path_guard

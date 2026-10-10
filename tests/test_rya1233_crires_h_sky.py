"""RYA-1233: the CRIRES+ H telluric term, priced on CRIRES+'s own sky (molecfit of our IDPs)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import crires_h_sky as H  # noqa: E402

pytestmark = pytest.mark.skipif(not list(H.SEG_DIR.glob("*_corrected.npz")),
                                reason="RYA-1233 molecfit segments not present")


def test_every_frame_has_one_velocity():
    vs = {}
    for s in H.segments():
        vs.setdefault(s["frame"], []).append(s["v_kms"])
        assert s["ccf"] > 0.5, s["tag"]
    assert all(max(v) - min(v) <= H.FRAME_TOL_KMS for v in vs.values())


def test_a_clear_line_prices_zero_and_a_sky_saturated_line_cannot():
    clear = H.telluric_line(15376.831, 0.9)
    assert clear["n_telluric_px"] == 0 and clear["residual_flux"] == 0.0
    sat = H.telluric_line(16434.927, 0.9)
    # Si 16434.927 sits on a saturated telluric line: unverifiable pixels are charged full
    # depth, never dropped (which would read as a perfect correction).
    assert sat["n_unverifiable_px"] > 0 and sat["residual_flux"] > 0.3

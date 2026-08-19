"""RYA-897 — the direct-solar HARPS arm is reachable, and stays reachable.

The arm every Codex abundance is differential to. `config` listed it as a solar arm and
the harness had no branch for it, so it was half-wired from RYA-713 until now.

Ported onto RYA-904's holding architecture: the original RYA-897 patch targeted the old
`WINDOW_LOADERS` / `_LOADER_HOLDING` dicts, which no longer exist.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))


def _mbe():
    import measure_band_ew
    return measure_band_ew


def test_harps_is_a_wired_arm():
    m = _mbe()
    holdings = m.holdings_for("harps")
    assert [h.holding_id for h in holdings] == ["solar_harps"]
    assert holdings[0].reader == "harps"


def test_harps_declares_itself_normalised():
    """🔴 THE CONTRACT RYA-911 CHANGED, and the reason is worth keeping attached.

    RYA-897 first served `flux_raw` with pre_normalised=False, reasoning that HARPS
    arrives un-normalised from the archive. RYA-911 measured the result: EWs ~16% low and
    abundances 0.34 dex low on lines the old path also kept. This product is NOT the
    archive S1D -- it is the pipeline's own output and carries the continuum the pipeline
    already fitted. Serving the normalised column uses that continuum instead of placing
    a second one on top of it.
    """
    m = _mbe()
    assert m.holdings_for("harps")[0].pre_normalised is True


def test_harps_serves_normalised_flux_not_raw():
    """Measured, not asserted: residual flux sits near unity; raw counts do not."""
    m = _mbe()
    if not m.HARPS_CSV.exists():
        pytest.skip(f"HARPS spectrum not staged at {m.HARPS_CSV}")
    win = m.load_window_ex("harps", 6147.734, 1.2)
    med = float(np.median(win.flux))
    assert 0.5 < med < 1.5, (
        f"HARPS window median flux is {med:.4g}. The holding declares "
        f"pre_normalised=True, so this must be residual flux near unity. A value in the "
        f"thousands means the RAW column is being served and the harness will place a "
        f"second continuum on it (RYA-713/911).")
    assert win.holding.holding_id == "solar_harps"


def test_harps_refuses_outside_its_span():
    """3782.6-6910 A. Beyond that it must refuse, not serve a truncated window."""
    m = _mbe()
    with pytest.raises(LookupError) as e:
        m.load_window_ex("harps", 9000.0, 1.2)
    assert "harps" in str(e.value)


def test_solar_harps_is_no_longer_a_declared_gap():
    """The exemption must go when the wiring lands, or the guard stops guarding.

    RYA-904 wrote the entry AND the check that catches it going stale -- deleting it
    without wiring would fail, and wiring without deleting would also fail. This asserts
    the pair is now consistent.
    """
    from pipeline import loader_coverage as lc
    assert "solar_harps" not in lc.DECLARED_GAPS
    assert lc.stale_gaps() == []
    lc.reconcile_loader_coverage()

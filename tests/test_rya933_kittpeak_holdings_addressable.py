"""RYA-933/934: all three Kitt Peak holdings must be reachable BY NAME.

Two telluric-corrected Kitt Peak products existed and neither could be named by
the band harness. That is the RYA-904 shape: an unreachable holding leaves no
trace, because it reads to every caller exactly like having no data. The
loader-coverage guard did not catch it either -- both were registered `audited`
rather than `verified`, and the guard only demands the latter.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))


@pytest.fixture(scope="module")
def harness(tmp_path_factory, ):
    kp = tmp_path_factory.mktemp("kp")
    (kp / "lm0296").touch()
    import os
    os.environ.setdefault("CODEX_KP_ATLAS", str(kp))
    import measure_band_ew as M
    return M


def test_all_three_kitt_peak_holdings_are_addressable(harness):
    from pipeline.loader_coverage import addressable_holdings
    a = addressable_holdings()
    for holding in ("solar_kpno",
                    "solar_kpno_molecfit_corrected",
                    "solar_kpno_kurucz2005_corrected"):
        assert a.get(holding) == "kpno_solar_atlas", f"{holding} cannot be named"


def test_the_whole_kitt_peak_class_is_declared_pre_normalised(harness):
    """RATIFIED RYA-1026, REVERSING RYA-929 for the 2005 product.

    Ryan's ratified rule is about the PRODUCT -- DO NOT NORMALISE ANY KITT PEAK ATLAS,
    the whole class, because a second continuum tilts the spectrum and it has bitten
    twice (RYA-940 on 1984, the 2005 double-normalise that forced the VIS re-run).

    RYA-933 fixed the ROUTE at the same time, and that ordering is the lesson: flipping
    the flag alone would have handed ABSOLUTE IRRADIANCE to a fit against a
    unity-normalised synthesis. Flag and reader move together. Independently measured on
    the two files: `irradthu.dat` gives a median rolling-P95 of 2.14 at 5000-5100 A
    (UN-NORMALISED) where the shipped residual `irradrelwl.dat` gives 0.9883 (NORMALISED).

    Which is why the flag is not the only witness -- see the route cross-check below.
    """
    pre = harness.PRE_NORMALISED
    assert pre["solar_kpno"] is True                      # residual flux
    assert pre["solar_kpno_molecfit_corrected"] is True    # same conventions
    # 🔴 REVERSED BY RYA-933/1026, and the reversal is the finding. This asserted False
    # on the premise that Kurucz 2005 "ships NO continuum, so the harness must place one".
    # It does ship one: `irradrelwl.dat`, the RESIDUAL atlas, distributed alongside
    # `irradthu.dat` but ABSENT FROM `0irrad.readme` -- which is why the RYA-929 intake
    # took only the irradiance file and concluded there was no continuum. Placing our own
    # instead tilted the band 4% blue-to-red (1.0238 at 4400 A to 0.9848 at 6800 A against
    # the 1984 atlas) and biased A(Fe I) low by 0.0218 +/- 0.0040 dex, correlated with
    # wavelength at r=+0.373. Reading the shipped residual removes the trend (r=-0.082).
    # The rule this now pins is RYA-911/938's: never re-fit a continuum a product ships.
    assert pre["solar_kpno_kurucz2005_corrected"] is True


def test_the_kurucz_route_is_cross_checked_against_the_flag(harness):
    """The flag alone cannot catch a MIS-ROUTED product: a declared flag and a mis-routed
    file agree with each other perfectly and are both wrong (the `iag_fts_solar_atlas`
    shape, RYA-944). The loader must therefore let the DATA vote."""
    source = (ROOT / "scripts" / "measure_band_ew.py").read_text()
    block = source[source.index("def load_kurucz2005_window"):]
    block = block[: block.index("def load_kp_window")]
    assert "assert_data_matches_declaration" in block
    assert "declared=True" in block


def test_the_corrected_holding_does_not_displace_the_original(harness):
    order = [h.holding_id for h in harness._INSTRUMENT_HOLDINGS["kpno_solar_atlas"]]
    assert order.index("solar_kpno") < order.index("solar_kpno_molecfit_corrected")


def test_kurucz_declares_its_span_so_the_ir_cannot_be_served_from_it(harness):
    """Nothing telluric-free reaches beyond 10000 A; the table must say so."""
    spec = next(h for h in harness._INSTRUMENT_HOLDINGS["kpno_solar_atlas"]
                if h.holding_id == "solar_kpno_kurucz2005_corrected")
    # 2990-10010 A: the RESIDUAL atlas's own declared range (299.000 to 1001.000 nm),
    # not the 3000-10000 rounding the irradiance file was registered with (RYA-933).
    assert spec.span_A == (2990.0, 10010.0)
    assert not spec.covers(11300.0, 2.0), "Kurucz must not claim the IR"
    assert spec.covers(6500.0, 2.0)


def test_the_corrected_1984_holding_refuses_uncorrected_windows(harness):
    """Only six bands were corrected. Elsewhere it must NOT fall back."""
    with pytest.raises(LookupError) as excinfo:
        harness.load_kp1984_corrected_window(5000.0, 2.0)
    message = str(excinfo.value)
    assert "not one of them" in message
    assert "do NOT relabel" in message


def test_kurucz_reader_declares_the_vacuum_grid(harness):
    """RYA-938 measured it; a reader that forgets displaces ~200 sampled pixels."""
    source = (ROOT / "scripts" / "measure_band_ew.py").read_text()
    block = source[source.index("def load_kurucz2005_window"):]
    block = block[: block.index("def load_kp_window")]
    # RYA-933: the reader now serves the RESIDUAL atlas, which is vacuum on the same
    # terms ("vacuum wavelength including gravitational red shift"), so the conversion
    # still has to be explicit -- it just lives in `_read_kurucz2005_residual` rather
    # than in a `True` flag passed to the crosscheck reader. The INVARIANT is unchanged
    # and is what this test protects: a reader that forgets displaces ~200 sampled pixels.
    assert "_read_kurucz2005_residual" in block, "must read the SHIPPED residual atlas"
    conv = source[source.index("def _read_kurucz2005_residual"):]
    conv = conv[: conv.index("def load_kurucz2005_window")]
    assert "vac_to_air" in conv, "the vacuum->air conversion must be explicit"
    assert "VACUUM" in block

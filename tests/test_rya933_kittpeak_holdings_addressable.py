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


def test_each_holding_declares_its_own_continuum_state(harness):
    """`pre_normalised` is holding-level, and these three genuinely disagree."""
    pre = harness.PRE_NORMALISED
    assert pre["solar_kpno"] is True                      # residual flux
    assert pre["solar_kpno_molecfit_corrected"] is True    # same conventions
    # Kurucz 2005 is ABSOLUTE irradiance and ships NO continuum, so the harness
    # must place one. Declaring it True would be the RYA-713 double-normalisation
    # trap in reverse: consuming unity as a continuum that is not there.
    assert pre["solar_kpno_kurucz2005_corrected"] is False


def test_the_corrected_holding_does_not_displace_the_original(harness):
    order = [h.holding_id for h in harness._INSTRUMENT_HOLDINGS["kpno_solar_atlas"]]
    assert order.index("solar_kpno") < order.index("solar_kpno_molecfit_corrected")


def test_kurucz_declares_its_span_so_the_ir_cannot_be_served_from_it(harness):
    """Nothing telluric-free reaches beyond 10000 A; the table must say so."""
    spec = next(h for h in harness._INSTRUMENT_HOLDINGS["kpno_solar_atlas"]
                if h.holding_id == "solar_kpno_kurucz2005_corrected")
    assert spec.span_A == (3000.0, 10000.0)
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
    assert "read_kurucz2005" in block
    assert ", True)" in block, "the vacuum flag must be passed explicitly"
    assert "VACUUM" in block

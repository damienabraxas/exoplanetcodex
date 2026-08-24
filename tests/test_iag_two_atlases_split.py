"""The two IAG atlases are SEPARATE HOLDINGS, and the corrected one wins where it reaches.

🔴 THE DEFECT THIS CLOSES. `iag_fts_solar_atlas` served TWO files with different reach and
different telluric states under ONE instrument_id -- the RYA-904 shape -- which is why
"is IAG telluric-corrected?" had no single answer for months. It has one now, per holding.

The scare was also half wrong, and the half that was wrong is instructive: a note asserted
the INSTRUMENT routes to the raw Reiners file. Only a stale row in
`solar_reference_holdings_rya708.csv` did. `measure_band_ew` has always opened Baker+2020,
and the served flux measures telluric-free. A claim about a catalogue row was carried as a
claim about the data.
"""
import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def _specs() -> dict:
    """holding_id -> {reader, span_A, pre_normalised}, parsed statically.

    Parsed, not imported: `measure_band_ew` resolves the Kitt Peak atlas at import time
    and SystemExits when it is not staged, so importing it would make this pass or fail on
    whether a data drive is mounted.
    """
    tree = ast.parse((ROOT / "scripts" / "measure_band_ew.py").read_text())
    # Module-level float constants, so a span that REFERENCES one (rather than retyping
    # the number) still resolves. The reference is the correct style -- retyping 5001.10
    # in two places is how the two ends of a boundary drift apart -- so the test bends to
    # the code here, not the other way round.
    consts = {t.id: n.value.value
              for n in tree.body if isinstance(n, ast.Assign)
              for t in n.targets
              if isinstance(t, ast.Name) and isinstance(n.value, ast.Constant)
              and isinstance(n.value.value, float)}

    def _num(e):
        if isinstance(e, ast.Constant):
            return e.value
        if isinstance(e, ast.Name):
            return consts.get(e.id)
        return None

    out = {}
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call)
                and getattr(node.func, "id", None) == "HoldingSpec"):
            continue
        kw = {k.arg: k.value for k in node.keywords}
        span = kw.get("span_A")
        out[node.args[0].value] = {
            "reader": getattr(kw.get("reader"), "value", None),
            "span": (tuple(_num(e) for e in span.elts)
                     if isinstance(span, ast.Tuple) else None),
        }
    return out


def test_both_iag_atlases_are_registered_as_separate_holdings():
    s = _specs()
    assert "solar_iag" in s and "solar_iag_reiners2016" in s
    assert s["solar_iag"]["reader"] == "iag"
    assert s["solar_iag_reiners2016"]["reader"] == "iag_reiners"


def test_each_iag_holding_declares_its_own_span():
    """🔴 `solar_iag` declared NO span, so `covers()` answered True for every window ever
    asked -- including ~954 A of blue Baker+2020 cannot reach at all. That is the RYA-767
    defect exactly, and it is what let one instrument_id look like it covered everything."""
    s = _specs()
    for hid in ("solar_iag", "solar_iag_reiners2016"):
        span = s[hid]["span"]
        assert span and all(isinstance(x, float) for x in span), f"{hid} declares no span"


def test_the_corrected_atlas_wins_wherever_it_reaches():
    """Selection order plus the blue arm's capped span. Where Baker+2020 covers a
    wavelength, the raw sibling must be unreachable -- otherwise the split would just be
    a new way to pick the wrong file."""
    s = _specs()
    baker_lo, baker_hi = s["solar_iag"]["span"]
    rei_lo, rei_hi = s["solar_iag_reiners2016"]["span"]
    assert rei_hi <= baker_lo, (
        "the blue arm must stop where the corrected atlas starts, or selection could "
        "serve raw flux for a wavelength the corrected atlas covers")
    assert rei_lo < baker_lo, "the blue arm must reach BELOW the corrected atlas"
    assert baker_hi > rei_hi


def test_the_two_together_leave_no_coverage_gap():
    """This was mistaken for a data gap. It is not one: the arms abut."""
    s = _specs()
    baker_lo, _ = s["solar_iag"]["span"]
    _, rei_hi = s["solar_iag_reiners2016"]["span"]
    assert rei_hi == pytest.approx(baker_lo, abs=1e-6), (
        "a hole between the two IAG arms would be a REAL data gap")


def test_the_spans_are_AIR_not_the_measured_vacuum_numbers():
    """🔴 Both atlases store VACUUM WAVENUMBER; `covers()` and every line list are AIR.
    The file boundary measures 5002.5 A in vacuum and 5001.10 A in air. Declaring the
    vacuum number would misplace the split by ~1.4 A and hand that sliver to the atlas
    that cannot serve it."""
    s = _specs()
    assert s["solar_iag"]["span"][0] == pytest.approx(5001.10, abs=0.01)
    assert s["solar_iag"]["span"][0] < 5002.5, "that is the VACUUM boundary, not air"


def test_the_stale_catalogue_row_no_longer_points_at_the_raw_file():
    """The row was the ONLY thing claiming IAG was uncorrected, and the reader never
    agreed with it."""
    rows = (ROOT / "data" / "catalog"
            / "solar_reference_holdings_rya708.csv").read_text().splitlines()
    iag = [r for r in rows if r.startswith("iag_fts_solar_atlas,")]
    assert len(iag) == 2, "both IAG arms must be catalogued"
    assert any("iag_baker2020" in r for r in iag), "the corrected atlas is not catalogued"
    assert any("iag_reiners2016" in r for r in iag), "the blue arm is not catalogued"


def test_the_retired_iag_anomaly_is_not_reinstated_by_prose():
    """It was retired because it was WRONG about the served flux, not because it was
    inconvenient. The explanation stays in the source; the entry must not come back."""
    from pipeline.telluric_display_policy import REGISTRY_ANOMALIES, anomaly
    assert "solar_iag" not in REGISTRY_ANOMALIES
    assert anomaly("solar_iag") is None

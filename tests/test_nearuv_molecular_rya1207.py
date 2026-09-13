"""RYA-1207 — molecular opacity in the near-UV: wired, scoped, and measured PAIRED.

PRODUCTION CHANGE. These pin the four things that would otherwise go wrong quietly:

  * molecules must reach the near-UV and NOWHERE ELSE. `use_molecules` is a band property;
    if it ever defaults on, every published band moves for an unmeasured reason;
  * the molecular lists must carry the ratified isotopologue codes -- an empty or
    mis-coded list reads as "molecules included" while contributing nothing;
  * the lever must be quoted PAIRED. Differencing against the published values measures
    the molecular change plus every code change since the feed was written;
  * RYA-1204's -0.087 must not be cited as this pool's number: the two pools share 2 lines.
"""
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "data/results/rya1207/nearuv_molecular_lever.json"


@pytest.fixture(scope="module")
def doc():
    if not DOC.exists():
        pytest.skip("RYA-1207 artifact absent")
    return json.loads(DOC.read_text())


def test_molecules_are_on_for_the_near_uv_and_off_for_every_other_band():
    """🔴 THE SCOPE IS THE WHOLE SAFETY ARGUMENT. Turbospectrum picks molecular lists by the
    nm range in the filename, so enabling this in another band would silently pull in
    whichever of the 99 shipped lists overlap there and move that band."""
    from config.synth_bands import SYNTH_BANDS
    assert SYNTH_BANDS["near-UV"].use_molecules is True
    for name, b in SYNTH_BANDS.items():
        if name != "near-UV":
            assert b.use_molecules is False, f"{name} would synthesise with molecules"


def test_the_generator_leaves_the_call_untouched_when_molecules_are_off():
    """The `_atm_file` discipline: an unused option must not appear in the call at all.
    RYA-770 stabilised the LTE path at -0.026 dex against the banked optical answer."""
    src = (ROOT / "pipeline/abundances_derive.py").read_text()
    assert '_mol = {} if not use_molecules else {"use_molecules": True}' in src
    assert src.count("**_atm_file, **_mol") == 2, (
        "both generate_spectrum call sites must splat it, or the NLTE branch silently "
        "synthesises without molecules while the LTE branch uses them")


def test_the_species_codes_name_the_isotopologue(doc):
    """The assignment that was refused for years, and why it is safe: it cannot move a
    line, only the abundance it enters with."""
    from pipeline.nearuv_linelist import MOLECULAR_SPECIES_CODES
    assert MOLECULAR_SPECIES_CODES == {
        "OH": "0108.000016", "NH": "0107.000014",
        "CN": "0607.012014", "CH": "0106.000012"}
    assert doc["molecular_lists"]["species_codes"] == MOLECULAR_SPECIES_CODES


def test_the_filename_carries_the_nm_range_iSpec_filters_on():
    """A correctly-formatted list whose NAME lacks `_<lo>-<hi>` is silently ignored by
    iSpec rather than refused (RYA-1204 flagged this). The name is load-bearing."""
    from pipeline.nearuv_linelist import molecular_bsyn_filename
    import re
    fn = molecular_bsyn_filename("OH", 3000.0, 3780.0)
    assert fn == "OH_300-378.bsyn"
    assert re.match(r"(.*)_(\d+)-(\d+)\.bsyn", fn), "iSpec's own glob pattern must match"


def test_the_lever_is_quoted_paired_not_against_the_published_values(doc):
    """🔴 The published near-UV values predate the feed by weeks. Differencing against them
    gives -0.082..-0.092, of which about -0.035 is code change, not molecules."""
    lev = doc["the_lever_measured_PAIRED"]
    assert "SAME CODE" in lev["method"]
    for v in lev["Fe_I_dex"].values():
        assert -0.06 < v < -0.04, f"Fe I lever out of the measured band: {v}"
    for v in lev["Fe_II_dex"].values():
        assert -0.35 < v < -0.25, f"Fe II lever out of the measured band: {v}"


def test_the_molecules_act_only_through_their_own_lines(doc):
    """The control. Windows with no molecular line must be unchanged to the bit; if they
    move, the patch is a global offset and the interpretation is wrong, not just the size."""
    iv = doc["the_lever_measured_PAIRED"]["independent_validation"]
    assert len(iv["bit_identical_windows"]) >= 3
    assert iv["chi2r_median_on"] < iv["chi2r_median_off"], (
        "fit quality must improve; a lever that only moves the abundance is not validated")


def test_the_shift_is_significant_and_downward(doc):
    p = doc["the_lever_measured_PAIRED"]["per_line_Fe_I_kurucz2005"]
    assert p["paired_t"] < -2.0, "the lever must be significant on the production pool"
    assert p["n_down"] > p["n_up"], "adding opacity must lower the derived abundance"


def test_rya1204s_number_is_not_claimed_as_this_pools_number(doc):
    """🔴 The two pools share 2 lines of 40. Quoting -0.087 here would be citing a
    measurement of a different line set."""
    v = doc["vs_rya1204"]
    assert v["🔴 lines_in_common"] == 2
    assert v["rya1204_reported"] == -0.087
    assert "disjoint" in v["finding"].lower()


def test_the_fe_I_block_is_recorded_and_not_attributed_to_the_molecules(doc):
    """The four Fe I products could not be published. That must stay visible, and must not
    be blamed on this change: the molecules-off run gives the identical refusal."""
    b = doc["🔴 blocked"]
    assert "gf rung 1" in b["why"]
    assert "SAME verdict" in b["not_caused_by_this_ticket"]
    assert set(b["measured_values_if_it_is_resolved"]) == {
        "kurucz2005_1D-LTE", "kurucz2005_ENGINE-A",
        "molecfit_1D-LTE", "molecfit_ENGINE-A"}

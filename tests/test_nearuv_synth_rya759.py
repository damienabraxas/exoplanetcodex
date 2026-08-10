"""RYA-759 — near-UV line-list conversion + the three loud-fail guards.

These run anywhere: nothing here needs iSpec, Turbospectrum or a grid. The parts that
DO need the engine (does Turbospectrum reach 3000 A, does the near-UV list synthesise)
are Sirius smoke steps in scripts/rya759_nearuv_synth.py, because a unit test that
mocked the engine would re-create exactly the failure Move 1 diagnosed: a harness
reporting success while the engine was never asked anything.
"""
import numpy as np
import pytest

from data.linelists.vald_parse import parse_vald_long
from pipeline import nearuv_synth as ns
from pipeline.nearuv_linelist import (EV_TO_CM1, NearUVLinelistError, band_stats,
                                      read_band, to_ispec_array)

# One VALD3 long-format record, verbatim in shape: data line + 2 config lines + ref.
VALD_RECORD = (
    "'Fe 1',       3705.56600, -1.234,  1.4850,  3.0,  4.8300,  4.0,  1.250,  0.700,"
    "  1.800, 8.010,-5.690,-7.770, 0.732,\n"
    "'  LS                                                       3d7.(4F).4s a5F'\n"
    "'  LS                                                3d6.(3D).4s.4p.(3P*) t3F*'\n"
    "'_          Kurucz Fe I 2014   5 wl:K14   5 gf:K14   5 K14    Fe'\n"
)
VALD_RECORD_2 = (
    "'Ti 2',       3759.29600, -0.400,  0.6070,  2.5,  3.9040,  3.5,  1.200,  1.400,"
    "  1.650, 8.260,-6.040,-7.580, 0.910,\n"
)
VALD_HEADER = (
    " 3000.00000, 3780.00000, 2, 2, 1.0 Wavelength region, lines selected\n"
    "Spec Ion       WL_air(A)  log gf* E_low(eV) J lo E_up(eV)  J up\n"
)


@pytest.fixture
def extract(tmp_path):
    p = tmp_path / "vald_nearuv.txt"
    p.write_text(VALD_HEADER + VALD_RECORD + VALD_RECORD_2)
    return p


# ── the parser extension is ADDITIVE ─────────────────────────────────────────

def test_parse_carries_upper_level_and_lande(extract):
    recs, rep = parse_vald_long(str(extract))
    assert rep["n_failures"] == 0
    fe = recs[0]
    # pre-existing keys, unchanged
    assert (fe["element"], fe["ion"]) == ("Fe", "I")
    assert fe["wavelength"] == pytest.approx(3705.566)
    assert fe["e_low_eV"] == pytest.approx(1.4850)
    assert fe["damping_vdW"] == pytest.approx(-7.770)
    # RYA-759 additions
    assert fe["j_low"] == pytest.approx(3.0)
    assert fe["e_up_eV"] == pytest.approx(4.8300)
    assert fe["j_up"] == pytest.approx(4.0)
    assert fe["lande_lower"] == pytest.approx(1.250)
    assert fe["lande_mean"] == pytest.approx(1.800)


def test_read_band_clips_and_reports(extract):
    band, rep = read_band(extract, 3700.0, 3750.0)
    assert [r["element"] for r in band] == ["Fe"]      # the Ti II line is out of band
    assert rep["n_in_band"] == 1 and rep["n_parsed"] == 2


def test_read_band_raises_rather_than_returning_empty(extract):
    with pytest.raises(NearUVLinelistError, match="0 lines"):
        read_band(extract, 3000.0, 3100.0)


# ── VALD -> iSpec conversion ─────────────────────────────────────────────────

CHEM = np.array([("Fe", 26), ("Ti", 22)],
                dtype=[("symbol", "U4"), ("atomic_num", "i4")])


def test_to_ispec_array_physics_columns(extract):
    recs, _ = parse_vald_long(str(extract))
    arr = to_ispec_array(recs, chem_elements=CHEM)
    assert len(arr) == 2
    fe = arr[arr["element"] == "Fe 1"][0]
    assert fe["wave_A"] == pytest.approx(3705.566)
    assert fe["wave_nm"] == pytest.approx(370.5566)
    assert fe["lower_state_cm1"] == pytest.approx(1.4850 * EV_TO_CM1, rel=1e-9)
    # Turbospectrum wants the statistical weight, not J.
    assert fe["upper_g"] == pytest.approx(2 * 4.0 + 1)
    # VALD's Rad. is log10(gamma); TS wants the rate.
    assert fe["turbospectrum_rad"] == pytest.approx(10.0 ** 8.010)
    # vdW is carried through to the field TS actually reads.
    assert fe["turbospectrum_fdamp"] == pytest.approx(-7.770)
    assert fe["ion"] == 1 and fe["molecule"] == "F" and fe["nlte"] == "F"
    assert fe["turbospectrum_species"] == "26.000000"
    assert fe["turbospectrum_support"] == "T"

    ti = arr[arr["element"] == "Ti 2"][0]
    assert ti["ion"] == 2
    assert ti["turbospectrum_species"] == "22.000000"   # stage lives in `ion`, per GES
    assert ti["spectrum_moog_species"] == "22.1"
    assert ti["width_species"] == "22.01"


def test_zero_radiative_damping_is_not_ten_to_the_zero(extract):
    recs, _ = parse_vald_long(str(extract))
    recs[0]["damping_rad"] = 0.0          # VALD supplied none
    arr = to_ispec_array(recs, chem_elements=CHEM)
    fe = arr[arr["element"] == "Fe 1"][0]
    assert fe["turbospectrum_rad"] == 0.0   # NOT 1.0 s^-1


def test_unknown_species_raises_rather_than_thinning_the_blend_forest(extract):
    recs, _ = parse_vald_long(str(extract))
    with pytest.raises(NearUVLinelistError, match="chemical-elements"):
        to_ispec_array(recs, chem_elements=CHEM[:1])   # Ti missing


def test_gf_source_tag_travels_into_reference_code(extract):
    recs, _ = parse_vald_long(str(extract))
    arr = to_ispec_array(recs, chem_elements=CHEM,
                         gf_sources={("Fe", 1, 3705.566): "K14"})
    assert arr[arr["element"] == "Fe 1"][0]["reference_code"] == "K14"
    assert arr[arr["element"] == "Ti 2"][0]["reference_code"] == "VALD3"


def test_band_stats_counts_the_awkward_facts(extract):
    recs, _ = parse_vald_long(str(extract))
    recs[0]["damping_vdW"] = 0.0
    s = band_stats(to_ispec_array(recs, chem_elements=CHEM))
    assert s["n_lines"] == 2 and s["n_species"] == 2
    assert s["n_vdw_zero"] == 1


# ── the three guards: each must RAISE ────────────────────────────────────────

LL_DTYPE = [("wave_A", "<f8"), ("element", "|U4")]


def test_empty_linelist_raises():
    with pytest.raises(ns.NearUVSynthesisError, match="0 lines"):
        ns.assert_linelist_covers(np.zeros(0, dtype=LL_DTYPE), 3000.0, 3780.0)


def test_linelist_that_misses_the_band_raises():
    ll = np.array([(5000.0, "Fe 1")], dtype=LL_DTYPE)
    with pytest.raises(ns.NearUVSynthesisError, match="0 lines"):
        ns.assert_linelist_covers(ll, 3000.0, 3780.0)


def test_band_covered_but_element_absent_raises():
    ll = np.array([(3500.0, "Ti 2")], dtype=LL_DTYPE)
    with pytest.raises(ns.NearUVSynthesisError, match="no Fe line"):
        ns.assert_linelist_covers(ll, 3000.0, 3780.0, element="Fe")


def test_linelist_covering_the_band_returns_the_count():
    ll = np.array([(3500.0, "Fe 1"), (3600.0, "Fe 1"), (5000.0, "Fe 1")],
                  dtype=LL_DTYPE)
    assert ns.assert_linelist_covers(ll, 3000.0, 3780.0, element="Fe") == 2


def test_all_zero_flux_raises():
    with pytest.raises(ns.NearUVSynthesisError, match="ALL-ZERO"):
        ns.assert_usable_flux(np.zeros(100), where="t")


def test_empty_flux_raises():
    with pytest.raises(ns.NearUVSynthesisError, match="EMPTY"):
        ns.assert_usable_flux(np.zeros(0), where="t")


def test_all_nan_flux_raises():
    with pytest.raises(ns.NearUVSynthesisError, match="non-finite"):
        ns.assert_usable_flux(np.full(50, np.nan), where="t")


def test_real_flux_passes():
    f = np.ones(50)
    f[10] = 0.4
    assert ns.assert_usable_flux(f, where="t")[10] == pytest.approx(0.4)


def test_missing_atmosphere_raises():
    with pytest.raises(ns.NearUVSynthesisError, match="returned nothing"):
        ns.assert_atmosphere(None, teff=5772.0, logg=4.438, feh=0.0, vturb=1.0,
                             model_grid="ATLAS9.Castelli")


def test_single_layer_atmosphere_raises():
    with pytest.raises(ns.NearUVSynthesisError, match="not an atmosphere"):
        ns.assert_atmosphere(np.zeros(1), teff=5772.0, logg=4.438, feh=0.0,
                             vturb=1.0, model_grid="ATLAS9.Castelli")


def test_flat_response_raises_even_when_lines_exist():
    """The subtle one: the list has lines, the synthesiser ignored them."""
    with pytest.raises(ns.NearUVSynthesisError, match="not synthesised from it"):
        ns.assert_sensitive(np.ones(100), np.ones(100), where="t")


def test_responsive_band_passes():
    lo, hi = np.ones(100), np.ones(100)
    hi[5] = 0.9
    assert ns.assert_sensitive(lo, hi, where="t") == pytest.approx(0.1)


# ── gf provenance is stated, never inferred ──────────────────────────────────

def test_near_uv_turns_canonical_gf_off_and_says_why():
    p = ns.gf_provenance(3000.0, 3780.0)
    assert p["apply_canonical_gf"] is False
    assert "canonical_gf.csv starts at" in p["detail"]


def test_optical_keeps_canonical_gf():
    p = ns.gf_provenance(4200.0, 6900.0)
    assert p["apply_canonical_gf"] is True
    assert "RYA-353" in p["detail"]

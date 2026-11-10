"""Synthetic boundary controls: detect wrong units/column/total density, not solar fits."""
from types import SimpleNamespace

import numpy as np
import pytest

from pipeline import bride_foundation as b


@pytest.fixture
def cube():
    # Non-square, independently distinguishable columns/depths.
    shape = (4, 2, 3)
    z, y, x = np.indices(shape)
    return dict(nz=4, ny=2, nx=3, z=np.array([8., 5., 2., -1.]),
                y=np.array([0., 2.]), x=np.array([0., 3., 6.]),
                temperature=3000. + 100. * z + 10. * y + x,
                electron_density=np.full(shape, 1.e15),
                hydrogen_density=np.full(shape, 2.e18),
                proton_density=np.full(shape, 3.e17),
                vx=1. + z + x, vy=2. + y + z, vz=3. + x + y + z)


def select(cube):
    return b.column_inputs(cube, star_id="solar", ix=2, iy=1,
                           vturb_m_s=np.zeros(cube["nz"]))


def test_identity_and_all_depths_survive_non_square_column_selection(cube):
    c = select(cube)
    assert (c.star_id, c.ix, c.iy) == ("solar", 2, 1)
    np.testing.assert_array_equal(c.temperature_K, [3012., 3112., 3212., 3312.])
    np.testing.assert_array_equal(c.vz_m_s, [6., 7., 8., 9.])
    np.testing.assert_array_equal(c.vx_m_s, [3., 4., 5., 6.])
    np.testing.assert_array_equal(c.vy_m_s, [3., 4., 5., 6.])


def test_total_hydrogen_includes_protons_and_never_reuses_neutral_only(cube):
    c = select(cube)
    np.testing.assert_array_equal(c.nH_total_m3, np.full(4, 2.3e18))
    assert np.all(c.nH_total_m3 > cube["hydrogen_density"][:, 1, 2])


def test_canonical_cgs_gravity_is_converted_to_si_not_passed_raw(cube, monkeypatch):
    from config import constants
    monkeypatch.setattr(constants, "get_star_params", lambda star: {"logg": 3.25})
    c = select(cube)
    assert c.logg_cgs == 3.25 and c.logg_si == 1.25
    assert 10 ** c.logg_si * 100 == pytest.approx(10 ** c.logg_cgs)
    # Positive wrong-unit control: raw cgs would accelerate by a factor 100.
    assert 10 ** c.logg_cgs / 10 ** c.logg_si == pytest.approx(100.)


def test_boundary_objects_and_gravity_are_explicit_at_actual_constructor_seam(cube):
    captured = {}
    def make(**kwargs):
        captured.update(kwargs)
        return "synthetic atmosphere"
    lw = SimpleNamespace(ScaleType=SimpleNamespace(Geometric="geometry"),
                         Atmosphere=SimpleNamespace(make_1d=make))
    lower, upper = object(), object()
    c = select(cube)
    assert b.make_vertical_atmosphere(c, lightweaver=lw, lower_bc=lower,
                                      upper_bc=upper) == "synthetic atmosphere"
    assert captured["logG"] == c.logg_si
    assert captured["lowerBc"] is lower and captured["upperBc"] is upper
    assert captured["scale"] == "geometry"
    np.testing.assert_array_equal(captured["nHTot"], c.nH_total_m3)
    np.testing.assert_array_equal(captured["vlos"], c.vz_m_s)
    assert "hydrogenPops" not in captured  # a Saha split is not NLTE level populations


def test_backend_mutation_cannot_change_column_or_source_cube(cube):
    c = select(cube)
    original = c.temperature_K.copy()
    def make(**kwargs):
        kwargs["temperature"][:] = 0.
    lw = SimpleNamespace(ScaleType=SimpleNamespace(Geometric=1),
                         Atmosphere=SimpleNamespace(make_1d=make))
    b.make_vertical_atmosphere(c, lightweaver=lw, lower_bc=object(), upper_bc=object())
    np.testing.assert_array_equal(c.temperature_K, original)
    np.testing.assert_array_equal(cube["temperature"][:, 1, 2], original)


@pytest.mark.parametrize("field,value", [("temperature", 0.), ("electron_density", -1.),
    ("hydrogen_density", -1.), ("proton_density", -1.), ("vz", np.nan), ("vx", np.inf)])
def test_invalid_physical_fields_fail_without_masking(cube, field, value):
    cube[field][1, 1, 2] = value
    with pytest.raises(b.FoundationInputError, match=field):
        select(cube)


@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_reversed_coordinate_or_duplicate_depth_is_not_repaired(cube, axis):
    cube[axis] = cube[axis][::-1]
    with pytest.raises(b.FoundationInputError, match=axis):
        select(cube)


def test_shape_mismatch_cannot_be_treated_as_a_different_column(cube):
    cube["temperature"] = cube["temperature"].transpose(0, 2, 1)
    with pytest.raises(b.FoundationInputError, match="temperature"):
        select(cube)


def test_no_default_microturbulence_or_boundaries_are_invented(cube):
    with pytest.raises(TypeError, match="vturb_m_s"):
        b.column_inputs(cube, star_id="solar", ix=0, iy=0)
    with pytest.raises(b.FoundationInputError, match="boundary"):
        b.make_vertical_atmosphere(select(cube), lightweaver=None, lower_bc=None, upper_bc=object())
    with pytest.raises(b.FoundationInputError, match="vturb"):
        b.column_inputs(cube, star_id="solar", ix=0, iy=0, vturb_m_s=-np.ones(4))


@pytest.mark.parametrize("index", [-1, 3, True, 1.5])
def test_invalid_column_index_is_refused(cube, index):
    with pytest.raises(b.FoundationInputError, match="ix"):
        b.column_inputs(cube, star_id="solar", ix=index, iy=0, vturb_m_s=np.zeros(4))


def test_pinned_real_constructor_preserves_supplied_thermodynamics(cube):
    lw = pytest.importorskip("lightweaver")
    from importlib.metadata import version
    if version("lightweaver") != "0.16.1":
        pytest.skip("adapter audited against 0.16.1 only")
    c = select(cube)
    atmosphere = b.make_vertical_atmosphere(c, lightweaver=lw,
        lower_bc=lw.ThermalisedRadiation(), upper_bc=lw.ZeroRadiation())
    np.testing.assert_array_equal(atmosphere.temperature, c.temperature_K)
    np.testing.assert_array_equal(atmosphere.ne, c.ne_m3)
    np.testing.assert_array_equal(atmosphere.nHTot, c.nH_total_m3)
    np.testing.assert_array_equal(atmosphere.height, c.height_m)


def test_requested_star_identity_controls_gravity_without_cross_star_default(cube, monkeypatch):
    from config import constants
    params = {"fixture_a": {"logg": 4.25}, "fixture_b": {"logg": 2.75}}
    requested = []
    def lookup(star):
        requested.append(star)
        return params[star]
    monkeypatch.setattr(constants, "get_star_params", lookup)
    columns = [b.column_inputs(cube, star_id=star, ix=0, iy=0,
                              vturb_m_s=np.zeros(4)) for star in params]
    assert requested == ["fixture_a", "fixture_b"]
    assert [c.star_id for c in columns] == requested
    assert [c.logg_si for c in columns] == [2.25, 0.75]

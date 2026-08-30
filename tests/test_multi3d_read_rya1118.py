"""RYA-1118 — the Multi3D ingest, and the three ways it can silently be wrong.

🔴 WHY A SYNTHETIC FIXTURE HERE. The real snapshot is 36 MB of external STAGGER data
that does not live in this repo (G3), so these tests build a Multi3D file byte-for-byte
to the documented layout and read it back. That is the right fixture for THIS module:
the thing under test is the LAYOUT CONTRACT — six float32 arrays in a fixed order, a
mesh that must agree with them, and a unit convention — not the physics of any
particular cube. A test that could only run beside a 36 MB binary would not run at all.

The layout is taken from Muspel.jl's `read_atmos_multi3d` (Pereira), the reference
implementation, and is pinned independently by byte count: 6 x nx x ny x nz x 4.
"""
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]

import sys  # noqa: E402

if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

import rya1118_multi3d_read as M  # noqa: E402

NX, NY, NZ = 3, 4, 5


def _write(tmp_path, arrays, *, nx=NX, ny=NY, nz=NZ):
    """A Multi3D pair written to the documented layout."""
    mesh = tmp_path / "mesh.test"
    xs = np.arange(nx) * 1.0e7          # cm
    ys = np.arange(ny) * 1.0e7
    zs = np.linspace(5.0e7, -3.0e7, nz)
    with open(mesh, "w") as fh:
        fh.write(f"{nx}\n" + " ".join(f"{v:.7e}" for v in xs) + "\n")
        fh.write(f"{ny}\n" + " ".join(f"{v:.7e}" for v in ys) + "\n")
        fh.write(f"{nz}\n" + " ".join(f"{v:.7e}" for v in zs) + "\n")
    atm = tmp_path / "atm3d.test"
    with open(atm, "wb") as fh:
        for a in arrays:
            # 🔴 order="F", NOT "C". Multi3D files are written by column-major (Julia)
            # code, so **x varies fastest** on disk. The first version of this helper
            # wrote row-major, which is z-fastest — so the FIXTURE encoded the same
            # misunderstanding as the reader under test, and the suite passed on a reader
            # that had x and z transposed. A fixture that shares the code's assumption
            # tests nothing about that assumption.
            fh.write(np.asarray(a, dtype=np.float32).astype("<f4").tobytes(order="F"))
    return mesh, atm


def _uniform(value):
    return np.full((NX, NY, NZ), value, dtype=np.float32)


def _ramped(base, *, dz=100.0, dy=1.0, dx=0.01):
    """A cube whose value encodes its own (x, y, z) index.

    🔴 EVERY EARLIER TEST HERE USED UNIFORM ARRAYS, WHICH CANNOT DETECT A STRIDING ERROR:
    permuting a constant cube returns the same constant cube. The real reader shipped with
    x and z transposed and this suite passed anyway — the bug surfaced only when the
    horizontally-averaged T(z) of a REAL snapshot came back oscillating instead of rising
    monotonically. So the fixture now varies along each axis with a different weight, and
    the read-back is checked ELEMENT BY ELEMENT.
    """
    a = np.empty((NX, NY, NZ), dtype=np.float32)
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                a[i, j, k] = base + dz * k + dy * j + dx * i
    return a


@pytest.fixture
def simple(tmp_path):
    # ne [cm^-3], T [K], vx, vy, vz [km/s], rho [g cm^-3]
    arrays = [_uniform(1.0e13), _uniform(6000.0), _uniform(1.0),
              _uniform(-2.0), _uniform(3.0), _uniform(1.0e-7)]
    return _write(tmp_path, arrays)


def test_units_are_converted_to_SI_not_left_in_cgs(simple):
    """🔴 The unit seam is where this breaks, and it breaks SILENTLY — cgs values are
    plausible-looking numbers, so nothing downstream would raise."""
    a = M.read_atmos_multi3d(*simple)
    # cm^-3 -> m^-3 is a factor 1e6, NOT 1e-6 and not 1
    assert a["electron_density"].max() == pytest.approx(1.0e13 * 1e6, rel=1e-6)
    # km/s -> m/s
    assert a["vx"].max() == pytest.approx(1.0e3, rel=1e-6)
    assert a["vz"].max() == pytest.approx(3.0e3, rel=1e-6)
    # mesh cm -> m
    assert a["z"].max() == pytest.approx(5.0e5, rel=1e-6)
    # temperature is already K and must NOT be scaled
    assert a["temperature"].max() == pytest.approx(6000.0)


def test_the_sixth_array_is_RHO_and_is_split_into_neutral_H_and_protons(simple):
    """The sixth block is mass density, not a hydrogen number density.

    Reading it as nH directly would be wrong by ~24 orders of magnitude and, worse,
    would still produce finite plausible-looking arrays.
    """
    a = M.read_atmos_multi3d(*simple)
    n_h_total = 1.0e-7 / (M.GRPH_G * M.CM_TO_M ** 3)
    assert (a["hydrogen_density"] + a["proton_density"]).max() == pytest.approx(
        n_h_total, rel=1e-6)
    # at 6000 K with this n_e hydrogen is overwhelmingly neutral
    assert a["hydrogen_density"].max() > 100 * a["proton_density"].max()


def test_axis_order_is_nz_ny_nx_not_the_on_disk_order(simple):
    """On disk the arrays are (nx, ny, nz); every consumer wants (nz, ny, nx).

    With NX != NY != NZ a transposition error changes the SHAPE, so this catches it.
    """
    a = M.read_atmos_multi3d(*simple)
    assert a["temperature"].shape == (NZ, NY, NX)
    assert (a["nx"], a["ny"], a["nz"]) == (NX, NY, NZ)


def test_a_variable_ordering_swap_would_be_caught_by_the_values(tmp_path):
    """Order matters: ne, T, vx, vy, vz, rho. Swapping ne and T is not detectable by
    size, only by value — so read a file whose blocks are distinguishable and check
    each landed in the right field."""
    arrays = [_uniform(2.0e13), _uniform(4500.0), _uniform(1.0),
              _uniform(2.0), _uniform(3.0), _uniform(2.0e-7)]
    a = M.read_atmos_multi3d(*_write(tmp_path, arrays))
    assert a["temperature"].max() == pytest.approx(4500.0)
    assert a["electron_density"].max() == pytest.approx(2.0e19, rel=1e-6)
    assert a["vy"].max() == pytest.approx(2.0e3, rel=1e-6)


def test_striding_matches_the_column_major_source_element_by_element(tmp_path):
    """🔴 THE REGRESSION TEST FOR THE BUG THIS SUITE MISSED.

    Multi3D is written by column-major (Julia) code, so on disk **x varies fastest**.
    numpy is row-major, so the natural `reshape(nx, ny, nz)` makes z fastest — the exact
    opposite — and the cube comes back scrambled while every marginal statistic (min, max,
    finite count, z range) stays correct. Only per-element structure exposes it.
    """
    arrays = [_uniform(1.0e13), _ramped(5000.0), _uniform(1.0),
              _uniform(2.0), _uniform(3.0), _uniform(1.0e-7)]
    a = M.read_atmos_multi3d(*_write(tmp_path, arrays))
    T = a["temperature"]
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                assert T[k, j, i] == pytest.approx(5000.0 + 100.0 * k + 1.0 * j + 0.01 * i), (
                    f"value at [z={k}, y={j}, x={i}] is not the one written there — "
                    "the cube is transposed or mis-strided")


def test_the_depth_profile_is_monotonic_not_scrambled(tmp_path):
    """The check that actually caught the real bug, in miniature.

    A scrambled cube keeps its VALUES and loses its STRUCTURE, so assert on the shape of
    T(z): with a pure z-ramp the horizontally-averaged profile must rise at every step.
    On the real snapshot the broken reader gave an oscillation between 7288 K and 7985 K;
    the fixed one rises monotonically from 4024 K to 11761 K across all 239 steps.
    """
    arrays = [_uniform(1.0e13), _ramped(5000.0, dy=0.0, dx=0.0),
              _uniform(1.0), _uniform(2.0), _uniform(3.0), _uniform(1.0e-7)]
    a = M.read_atmos_multi3d(*_write(tmp_path, arrays))
    prof = a["temperature"].mean(axis=(1, 2))

    # ⚠️ ASSERT THE PROFILE, NOT MERELY THAT IT RISES. A monotonic-only check PASSED
    # against the transposed reader on a fixture this small — the scramble happened to
    # leave the z-mean rising, so the test proved nothing about the bug it was written
    # for. (On the real 240-layer snapshot the same diagnostic did catch it: 7288<->7985
    # oscillating vs 4024->11761 rising.) A miniature needs the exact profile.
    expected = 5000.0 + 100.0 * np.arange(NZ)
    assert prof == pytest.approx(expected), (
        f"T(z) is not the ramp that was written: got {prof}, expected {expected}")
    assert np.all(np.diff(prof) > 0)


def test_a_mesh_that_disagrees_with_the_cube_is_REFUSED(tmp_path):
    """🔴 The byte-count guard. A mesh/cube mismatch reinterprets one variable as
    another and yields a full, finite, entirely wrong atmosphere. It must raise."""
    # a VALID mesh against a cube carrying only five of the six blocks
    mesh, atm = _write(tmp_path, [_uniform(1.0e13)] * 5)
    with pytest.raises(ValueError, match="bytes, expected"):
        M.read_atmos_multi3d(mesh, atm)


def test_a_truncated_mesh_names_the_short_read_not_a_negative_leftover(tmp_path):
    """🔴 Found by this suite. A mesh claiming more depth points than it carries used
    to report "-1 numbers left over", which is arithmetic nonsense and describes the
    wrong fault. It must name the short read."""
    mesh, atm = _write(tmp_path, [_uniform(1.0e13)] * 6)
    bad = tmp_path / "mesh.short"
    bad.write_text(mesh.read_text().replace(f"{NZ}\n", f"{NZ + 1}\n", 1))
    with pytest.raises(ValueError, match="claims nz=.* but only"):
        M.read_atmos_multi3d(bad, atm)


def test_a_mesh_with_trailing_numbers_is_REFUSED(tmp_path):
    """Leftover numbers mean the mesh does not describe this cube — never guess."""
    arrays = [_uniform(1.0e13)] * 6
    mesh, atm = _write(tmp_path, arrays)
    bad = tmp_path / "mesh.extra"
    bad.write_text(mesh.read_text() + "\n1.0 2.0 3.0\n")
    with pytest.raises(ValueError, match="left over"):
        M.read_atmos_multi3d(bad, atm)


def test_saha_ionfrac_is_bounded_and_monotonic_in_temperature():
    """The split must stay a fraction. An unbounded Saha ratio is how hydrogen goes
    negative or over-ionises, and both would surface only as a bad opacity."""
    ne = 1.0e19
    f = M.h_ionfrac_saha(np.array([3000.0, 6000.0, 10000.0, 50000.0]), ne)
    assert np.all((f >= 0.0) & (f <= 1.0))
    assert np.all(np.diff(f) > 0), "ionised fraction must rise with temperature"

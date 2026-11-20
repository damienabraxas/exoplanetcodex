"""RYA-1282 G5: the 1.5D column adapter reproduces the validated 1D solve exactly.

Runs under the held Lightweaver interpreter (engines.bride_eval); skips elsewhere.
"""
from __future__ import annotations

import numpy as np
import pytest

lw = pytest.importorskip('lightweaver')

from config.constants import require_codex_path  # noqa: E402
from pipeline.bride_columns import column_context, marcs_column, solve_column, stagger_column  # noqa: E402
from pipeline.bride_reference import make_reference_context, read_marcs_structure  # noqa: E402
from pipeline.bride_parameter_checks import controlled_structure, marcs_filename  # noqa: E402
from pipeline.multi_atom import read_multi_atom  # noqa: E402


@pytest.fixture(scope='module')
def solar():
    path = require_codex_path('grids.marcs_standard') / marcs_filename(5750, 4.5)
    if not path.exists():
        pytest.skip('MARCS 5750/4.5 not held on this host')
    structure = controlled_structure(read_marcs_structure(path), xi=1.0, drop_bottom=0)  # G4 baseline
    atom = read_multi_atom(require_codex_path('experiments.bride_oxygen_atom'))
    return structure, atom


def test_the_column_builder_solves_a_marcs_column_bit_for_bit_like_the_reference(solar):
    structure, atom = solar
    ref = make_reference_context(structure, atom, rays=3, qscale=1, ng=True)
    n_ref, h_ref, *_ = solve_column(ref, maxiter=1500)
    col = column_context(marcs_column(structure), atom, abundances_dex=structure['abundances_dex'],
                         logg_cgs=structure['logg_cgs'], vturb_m_s=structure['vturb_m_s'],
                         rays=3, qscale=1, ng=True, top='vacuum',   # the reference's boundary
                         background='lightweaver')                  # and its background
    n_col, h_col, *_ = solve_column(col, maxiter=1500)
    assert len(h_ref) == len(h_col)
    assert np.array_equal(n_ref, n_col)


def _homogeneous_cube(structure, nx=2, ny=3):
    """A snapshot whose every column is the same MARCS column, in Multi3D conventions."""
    col = marcs_column(structure)
    nz = col['ndep']
    tile = lambda a: np.broadcast_to(a[:, None, None], (nz, ny, nx)).copy()
    return dict(z=col['height_m'].copy(), temperature=tile(col['temperature_K']),
                electron_density=tile(col['ne_m3']),
                hydrogen_density=tile(col['nHtot_m3'] * 0.75), proton_density=tile(col['nHtot_m3'] * 0.25),
                vz=tile(np.zeros(nz)))


def test_every_column_of_a_homogeneous_snapshot_is_the_1d_column(solar):
    structure, _ = solar
    cube = _homogeneous_cube(structure)
    want = marcs_column(structure)
    for ix in range(2):
        for iy in range(3):
            got = stagger_column(cube, ix, iy)
            for key in ('height_m', 'temperature_K', 'ne_m3', 'vlos_m_s'):
                assert np.array_equal(got[key], want[key]), key
            assert np.allclose(got['nHtot_m3'], want['nHtot_m3'], rtol=1e-15, atol=0)


def test_a_snapshot_on_a_height_axis_or_upside_down_is_refused(solar):
    structure, _ = solar
    cube = _homogeneous_cube(structure)
    flipped = dict(cube, z=cube['z'][::-1])
    with pytest.raises(ValueError, match='top first'):
        stagger_column(flipped, 0, 0)
    upside = {k: (v[::-1] if k in ('temperature',) else v) for k, v in cube.items()}
    with pytest.raises(ValueError, match='hotter at depth'):
        stagger_column(upside, 0, 0)


def test_an_upflow_in_the_snapshot_stays_a_positive_upward_velocity(solar):
    structure, _ = solar
    cube = _homogeneous_cube(structure)
    cube['vz'][:] = 1000.  # +z is up in the held snapshot: a 1 km/s upflow
    assert np.all(stagger_column(cube, 0, 0)['vlos_m_s'] == 1000.)


def test_the_held_snapshot_reads_as_height_with_rising_hot_gas():
    import os
    base = os.path.expanduser('~/bride_rya1118/')
    stem = 't5777g44m0005_00091_80x80x240'
    if not os.path.exists(base + 'atm3d.' + stem):
        pytest.skip('STAGGER snapshot not held on this host')
    import sys
    sys.path.insert(0, base)
    from multi3d_read import read_atmos_multi3d
    cube = read_atmos_multi3d(base + 'mesh.' + stem, base + 'atm3d.' + stem)
    k = int(np.argmin(abs(cube['z'])))
    corr = np.corrcoef(cube['temperature'][k].ravel(), cube['vz'][k].ravel())[0, 1]
    assert corr > 0.3            # hot granules rise: +vz is up
    col = stagger_column(cube, 0, 0)
    assert col['height_m'][0] > 0 > col['height_m'][-1]

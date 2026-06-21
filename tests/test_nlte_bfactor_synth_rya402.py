"""
tests/test_nlte_bfactor_synth_rya402.py
=======================================
RYA-402 — the testable core: the PySME DirectAccessFile (.grd) departure-grid
reader, the b-factor->delta extraction math, and the fail-loud/no-fake contract.

The reader is exercised two ways: against a hand-built synthetic .grd fixture
(always runs, pins the binary format) and against the real vendored Na grid if it
is present (skipped in CI where the multi-GB .grd is gitignored). The live
Turbospectrum synthesis + Na/Al/K results are gated on the Gerber-2022 TS
adaptation and validated by the --validate-against Na STOP-gate, not here.
"""
from pathlib import Path

import numpy as np
import pytest

from pipeline import nlte_bfactor_synth as bf

_NUMPY_TO_IDL = {np.dtype('float64'): 5, np.dtype('float32'): 4, np.dtype('int32'): 3}


def _write_directaccess(path: Path, entries: dict):
    """Write a minimal PySME 'DirectAccess file Version 1.10' .grd from
    {key: ndarray}. Mirrors the reader: 64B version, (nblocks u8, dir_length i2,
    ndir u8) header, ndir x (key S256, size[23] i4, pointer i8), then C-order data
    with the IDL-reversed shape recorded."""
    version = b'DirectAccess file Version 1.10 2020-05-28'.ljust(64, b' ')
    ddt = np.dtype([('key', 'S256'), ('size', '<i4', 23), ('pointer', '<i8')])
    n = len(entries)
    dir_bytes = 64 + 8 + 2 + 8 + n * ddt.itemsize    # data starts after the directory
    directory = np.zeros(n, ddt)
    blobs, ptr = [], dir_bytes
    for i, (key, arr) in enumerate(entries.items()):
        arr = np.ascontiguousarray(arr)
        directory['key'][i] = key.encode()
        size = np.zeros(23, '<i4')
        size[0] = arr.ndim
        size[1:1 + arr.ndim] = arr.shape[::-1]         # IDL (column-major) order
        size[1 + arr.ndim] = _NUMPY_TO_IDL[arr.dtype]
        size[2 + arr.ndim] = arr.size
        directory['size'][i] = size
        directory['pointer'][i] = ptr
        blobs.append(arr.tobytes())
        ptr += len(blobs[-1])
    with open(path, 'wb') as f:
        f.write(version)
        np.array([(len(blobs), ddt.itemsize, n)],
                 np.dtype([('a', '<u8'), ('b', '<i2'), ('c', '<u8')])).tofile(f)
        directory.tofile(f)
        for blob in blobs:
            f.write(blob)


def _fixture(path):
    nlevel, ndepth = 4, 6
    # one departure block; deep layers thermalised (b->1)
    b = np.linspace(0.8, 1.0, ndepth)[None, :] * np.ones((nlevel, 1))
    _write_directaccess(path, {
        'teff': np.array([5000.0, 5750.0, 6000.0], np.float64),
        'grav': np.array([4.0, 4.5], np.float64),
        'feh': np.array([-0.5, 0.0], np.float64),
        'abund': np.array([-1.0, 0.0, 1.0], np.float64),
        'energy': np.zeros(nlevel, np.float64),
        'J': np.zeros(nlevel, np.float64),
        't5750_g+4.50_m+0.00_a+0.00': b.astype(np.float64),
    })
    return nlevel, ndepth


# ── PySME .grd reader (synthetic fixture) ────────────────────────────────────

def test_reader_parses_synthetic_grd(tmp_path):
    grd = tmp_path / 'nlte_X_test_pysme.grd'
    nlevel, ndepth = _fixture(grd)
    g = bf.PySMEGrid(grd)
    assert g.version.startswith('DirectAccess file Version 1.10')
    np.testing.assert_allclose(g.teff, [5000, 5750, 6000])
    np.testing.assert_allclose(g.abund, [-1, 0, 1])
    assert g.n_models == 1 and g.n_levels == nlevel
    dep = g.departures(5750, 4.5, 0.0, 0.0)
    assert dep.shape == (nlevel, ndepth)
    assert dep[:, -1].mean() == pytest.approx(1.0)        # deep layer thermalised


def test_reader_nearest_node_snap(tmp_path):
    grd = tmp_path / 'nlte_X_test_pysme.grd'
    _fixture(grd)
    g = bf.PySMEGrid(grd)
    # only one model node exists -> any query snaps to it
    j = g.nearest_node(6000, 4.0, -0.5, 1.0)
    assert g.node_keys[j] == 't5750_g+4.50_m+0.00_a+0.00'


def test_reader_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        bf.PySMEGrid(tmp_path / 'nope.grd')
    with pytest.raises(FileNotFoundError):
        bf.read_amarsi_grid('Na', base_dir=tmp_path)


def test_get_unknown_key_raises(tmp_path):
    grd = tmp_path / 'nlte_X_test_pysme.grd'
    _fixture(grd)
    with pytest.raises(KeyError):
        bf.PySMEGrid(grd).get('not_a_key')


# ── real vendored Na grid (skips if the gitignored .grd is absent) ───────────

def _na_grid_path():
    p = bf._AMARSI_DIR / 'nlte_Na_scatt_pysme.grd'
    return p if p.exists() else None


@pytest.mark.skipif(_na_grid_path() is None, reason="vendored Na .grd not present (gitignored)")
def test_real_na_grid_physics():
    g = bf.read_amarsi_grid('Na')
    assert g.n_models == 43680 and g.n_levels == 140
    assert g.teff.min() == 2500 and g.teff.max() == 8000
    assert g.feh.max() == pytest.approx(1.0)              # covers 55 Cnc
    b = g.departures(5772, 4.44, 0.0, 0.0)
    assert b.shape == (140, 56)
    assert b[:, -1].mean() == pytest.approx(1.0, abs=1e-3)  # deep layers -> LTE (the format check)
    assert 0.1 < b.min() < b.max() < 10                   # physical departures


@pytest.mark.skipif(_na_grid_path() is None, reason="vendored Na .grd not present (gitignored)")
def test_level_and_tau_mapping():
    g = bf.read_amarsi_grid('Na')
    assert g.level_for_energy(2.102) in (1, 2)            # Na I 3p
    assert g.level_for_energy(4.283) in range(9, 15)      # Na I 4d
    with pytest.raises(ValueError):
        g.level_for_energy(99.0)
    tau = g.node_tau(5772, 4.44, 0.0)
    assert tau.shape == (56,) and tau.min() < 0.01 < tau.max()


@pytest.mark.skipif(_na_grid_path() is None, reason="vendored Na .grd not present (gitignored)")
def test_stop_gate_fires_on_shortcut():
    """The guard MUST refuse the departure-only shortcut: it gives ~-0.01, the anchor
    is -0.107, so validate_against raises (no tuning to the anchor, no Al/K from it)."""
    with pytest.raises(RuntimeError, match="STOP-gate"):
        bf.validate_against('Na')


# ── b-factor -> delta extraction (pure math) ─────────────────────────────────

def test_delta_recovers_known_shift_identity_cog():
    assert bf.delta_from_ews(lambda A: A, ew_nlte=5.107, a_lte0=5.0) == pytest.approx(0.107, abs=1e-3)
    assert bf.delta_from_ews(lambda A: A, ew_nlte=4.893, a_lte0=5.0) == pytest.approx(-0.107, abs=1e-3)


def test_delta_with_curved_cog():
    ew = lambda A: 10.0 ** (A - 7.0)
    delta = bf.delta_from_ews(ew, ew_nlte=ew(6.40), a_lte0=6.20, a_lo=5.0, a_hi=7.5)
    assert delta == pytest.approx(0.20, abs=2e-3)


def test_delta_unbracketed_raises():
    with pytest.raises(ValueError, match="not bracketed"):
        bf.delta_from_ews(lambda A: A, ew_nlte=99.0, a_lte0=5.0, a_lo=4.0, a_hi=6.0)


# ── fail-loud / no-fake contract ─────────────────────────────────────────────

def test_live_synth_fails_loud_without_depgrid():
    with pytest.raises((RuntimeError, NotImplementedError)):
        bf.synth_ew_nlte_vs_lte('Al', 3961.52, {'teff': 5772, 'logg': 4.44, 'feh': 0.0})


def test_known_delta_anchor_is_the_inspect_na():
    assert bf._KNOWN_DELTA['Na']['delta'] == pytest.approx(-0.107)
    assert bf.TARGET_ELEMENTS == ['Al', 'K']

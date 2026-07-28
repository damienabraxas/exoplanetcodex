"""
tests/test_grid_availability_rya462.py
======================================
RYA-462 Part A — the NLTE grid-availability matrix.

Pins the audit contract: the matrix is driven by the live registries + the files on
disk (not a hand list), every wired registry grid is present, the C/O synthesis path
is wired while N I atomic NLTE stays owed, the PySME .grd inputs are flagged as
offline-derivation-only, and — the headline gate — there are ZERO present-but-unwired
production grids after the K wiring.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import nlte_grid_availability as M  # noqa: E402


def test_matrix_schema_and_size():
    df = M.build_matrix()
    assert {'element', 'subsystem', 'grid_file', 'present', 'wired',
            'role', 'source_if_absent', 'note'} <= set(df.columns)
    assert len(df) >= 24


def test_no_present_but_unwired_orphans():
    # The K case RYA-462 closes: a CSV grid on disk whose element is in no registry.
    df = M.build_matrix()
    orph = M.orphans(df)
    assert len(orph) == 0, f"PRESENT-but-UNWIRED grids remain: {orph.to_dict('records')}"


def test_K_present_and_wired():
    df = M.build_matrix()
    krow = df[(df.element == 'K') & (df.subsystem == 'registry-nlte')]
    assert len(krow) == 1
    r = krow.iloc[0]
    assert bool(r['present']) and bool(r['wired'])
    assert r['grid_file'] == 'K_Amarsi2020_PySME.csv'
    assert r['role'] == 'production'


def test_fe_and_registry_grids_present_and_wired():
    df = M.build_matrix()
    for el in ('Fe', 'Ca', 'Cr', 'Ti', 'Na', 'Mg', 'Si', 'Al', 'S', 'Ba', 'Sr', 'Mn'):
        r = df[(df.element == el) & (df.subsystem.isin(['fe-nlte', 'registry-nlte']))]
        assert not r.empty and bool(r.iloc[0]['present']) and bool(r.iloc[0]['wired']), el


def test_N_atomic_owed_not_wired_with_source():
    df = M.build_matrix()
    n = df[(df.element == 'N') & (df.subsystem == 'cno-3dnlte')].iloc[0]
    assert not bool(n['wired'])              # N I atomic NLTE owed (RYA-369)
    assert n['source_if_absent']             # the owed source is named
    for el in ('C', 'O'):                    # C/O synthesis NLTE IS wired
        r = df[(df.element == el) & (df.subsystem == 'cno-3dnlte')].iloc[0]
        assert bool(r['wired'])


def test_grd_inputs_are_offline_not_production():
    df = M.build_matrix()
    grd = df[df.subsystem == 'pysme-grd']
    assert not grd.empty
    assert (grd['role'] == 'offline-derivation').all()
    assert (~grd['wired'].astype(bool)).all()   # never counted as a wired production grid


def test_cu_grid_absent_as_csv():
    df = M.build_matrix()
    cu = df[(df.element == 'Cu') & (df.subsystem == 'registry-nlte')].iloc[0]
    assert not bool(cu['present'])           # no production CSV grid derived
    assert cu['source_if_absent']            # documented (line-quality blocker, RYA-395)

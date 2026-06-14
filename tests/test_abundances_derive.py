# tests/test_abundances_derive.py
"""
Layer 1 unit tests for iSpec + Turbospectrum installation and pipeline
integration. Fast where possible; atmosphere interpolation test requires
the full iSpec install (skipped if ISPEC_DIR missing).

RYA-167
"""
import pytest
import sys
import numpy as np
from pathlib import Path
from config.constants import ISPEC_DIR, RADIATIVE_TRANSFER_CODE


# ── Installation checks ───────────────────────────────────────────────────────

def test_ispec_dir_exists():
    assert ISPEC_DIR.exists(), (
        f"ISPEC_DIR not found: {ISPEC_DIR}\n"
        "Run: git clone https://github.com/marblestation/iSpec.git ispec"
    )

def test_ispec_importable():
    sys.path.insert(0, str(ISPEC_DIR))
    import ispec
    assert callable(getattr(ispec, 'interpolate_atmosphere_layers', None)), \
        "ispec.interpolate_atmosphere_layers not found — iSpec install may be incomplete"

def test_turbospectrum_binaries_exist():
    """babsma_lu and bsyn_lu must be compiled before Turbospectrum can run."""
    ts_dir = ISPEC_DIR / 'synthesizer' / 'turbospectrum' / 'exec-gf'
    babsma = ts_dir / 'babsma_lu'
    bsyn   = ts_dir / 'bsyn_lu'
    assert babsma.exists(), f"Turbospectrum babsma_lu not found: {babsma}\nRun: cd {ts_dir} && make"
    assert bsyn.exists(),   f"Turbospectrum bsyn_lu not found: {bsyn}"

def test_atlas9_grid_exists():
    atm = ISPEC_DIR / 'input' / 'atmospheres' / 'ATLAS9.Castelli'
    assert atm.exists(), "ATLAS9.Castelli grid not found — run: tar -xzf input.tar.gz"

def test_marcs_grid_exists():
    atm = ISPEC_DIR / 'input' / 'atmospheres' / 'MARCS.GES'
    assert atm.exists(), "MARCS.GES grid not found — run: tar -xzf input.tar.gz"

def test_solar_abundances_file_exists():
    abund = ISPEC_DIR / 'input' / 'abundances' / 'Asplund.2009' / 'stdatom.dat'
    assert abund.exists(), "Asplund.2009 abundances file missing"

def test_turbospectrum_line_regions_exist():
    regions = ISPEC_DIR / 'input' / 'regions' / '42000_GES' / \
        'limited_but_with_missing_elements_turbospectrum_synth_good_for_abundances_all_extended.txt'
    assert regions.exists(), "Multi-element Turbospectrum line regions file missing"

def test_radiative_transfer_code_is_turbospectrum():
    assert RADIATIVE_TRANSFER_CODE == 'turbospectrum', \
        f"Expected 'turbospectrum', got '{RADIATIVE_TRANSFER_CODE}'"


# ── API integration tests (require iSpec install) ─────────────────────────────

@pytest.fixture(scope='module')
def atlas9_pack():
    sys.path.insert(0, str(ISPEC_DIR))
    import ispec
    grid = ISPEC_DIR / 'input' / 'atmospheres' / 'ATLAS9.Castelli'
    if not grid.exists():
        pytest.skip("ATLAS9.Castelli grid not found")
    return ispec.load_modeled_layers_pack(str(grid))

def test_atmosphere_interpolation_solar(atlas9_pack):
    """ATLAS9 + Turbospectrum should interpolate at solar params."""
    sys.path.insert(0, str(ISPEC_DIR))
    import ispec
    target = {'teff': 5777.0, 'logg': 4.44, 'MH': 0.0, 'alpha': 0.0}
    atm = ispec.interpolate_atmosphere_layers(atlas9_pack, target, code='turbospectrum')
    assert atm is not None
    assert len(atm) > 0, "Atmosphere interpolation returned empty result"
    assert len(atm) == 72, f"Expected 72 atmosphere layers, got {len(atm)}"

def test_atmosphere_valid_target(atlas9_pack):
    sys.path.insert(0, str(ISPEC_DIR))
    import ispec
    solar = {'teff': 5777.0, 'logg': 4.44, 'MH': 0.0, 'alpha': 0.0}
    assert ispec.valid_atmosphere_target(atlas9_pack, solar)

def test_marcs_interpolation_mdwarf():
    """MARCS.GES should interpolate at M dwarf parameters (LHS 3844-like)."""
    sys.path.insert(0, str(ISPEC_DIR))
    import ispec
    grid = ISPEC_DIR / 'input' / 'atmospheres' / 'MARCS.GES'
    if not grid.exists():
        pytest.skip("MARCS.GES grid not found")
    pack = ispec.load_modeled_layers_pack(str(grid))
    target = {'teff': 3080.0, 'logg': 5.06, 'MH': 0.22, 'alpha': 0.0}
    if not ispec.valid_atmosphere_target(pack, target):
        pytest.skip("LHS 3844 params outside MARCS grid — expected for coarse grid")
    atm = ispec.interpolate_atmosphere_layers(pack, target, code='turbospectrum')
    assert atm is not None and len(atm) > 0


# ── EW injection tests ────────────────────────────────────────────────────────

def test_ew_injection_matches_fe_lines():
    """_build_ispec_line_regions should match Fe I lines from solar_ew.csv.

    RYA-297: re-pointed from the renamed _inject_ew_into_regions /
    _LINE_REGIONS_TS_ALL (function + constant were renamed, not removed — the
    EW-injection capability is intact in _build_ispec_line_regions).
    """
    import pandas as pd
    from config.constants import PATHS
    from pipeline.abundances_derive import _build_ispec_line_regions, _LINE_REGIONS_ALL

    ew_path = PATHS['solar_ew']
    if not Path(str(ew_path)).exists():
        pytest.skip("solar_ew.csv not found — run lines_fit.py first")

    ew_df = pd.read_csv(ew_path)
    # RYA-297: use the full Fe set, not .head(50). solar_ew.csv is wavelength-
    # sorted, so head(50) was all blue (<3831 Å) lines that fall outside the
    # curated "good_for_abundances" region file → zero matches. The injection
    # capability is fine: the full set matches ~124 Fe I/II lines.
    fe_ew = ew_df[(ew_df['element'] == 'Fe') & (ew_df['ew_mA'] > 0)]

    if len(fe_ew) == 0:
        pytest.skip("No Fe lines in solar_ew.csv")

    linemasks = _build_ispec_line_regions(fe_ew, _LINE_REGIONS_ALL)
    assert len(linemasks) > 0, "No Fe lines matched iSpec regions"
    assert np.all(linemasks['ew'] > 0), "Some matched lines have zero EW"
    # All should be Fe
    notes = [str(n) for n in linemasks['note']]
    assert all('Fe' in n for n in notes), "Non-Fe lines in Fe-only injection"

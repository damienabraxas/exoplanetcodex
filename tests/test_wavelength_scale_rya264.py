"""
tests/test_wavelength_scale_rya264.py
=====================================
RYA-264 — the third axis of the loader contract: wavelength UNIT (nm/µm→Å) + SCALE
(vacuum→air). Guards the SSOT converter, the cited per-instrument registry, the
unit-sanity band gate (the ×10 / ×10000 catch that would have killed the RYA-263
zero class), the scale sign-check, and the SPIRou/UVES loader wiring.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import frame_object_contract as foc          # noqa: E402
from pipeline import wavelength_util as wu                  # noqa: E402


# ── SSOT: exactly one vac↔air implementation, reused everywhere ───────────────
def test_single_vac_air_implementation_is_shared():
    # uv_conditioning and frame_object_contract both bind the SAME function object
    from pipeline import uv_conditioning as uvc
    assert uvc.vac_to_air is wu.vac_to_air
    assert uvc.air_to_vac is wu.air_to_vac
    assert foc.vac_to_air is wu.vac_to_air
    # uv_line_selection / crires_telluric delegate to the shared util (no 2nd formula)
    from pipeline import uv_line_selection as uvls
    from pipeline import crires_telluric as ct
    assert abs(uvls._air_to_vac(5000.0) - float(wu.air_to_vac(5000.0))) < 1e-9
    assert abs(float(np.asarray(ct._air_to_vac(np.array([23000.0])))[0]
                     - wu.air_to_vac(np.array([23000.0]))[0])) < 1e-9


def test_vac_air_roundtrip_and_boundary():
    w = np.array([3000.0, 6562.79, 15000.0, 23000.0])
    back = wu.air_to_vac(wu.vac_to_air(w))
    assert np.allclose(back, w, atol=1e-3)              # round-trips to ≪1 mÅ
    # below 2000 Å air is undefined → identity (FUV stays vacuum)
    assert wu.vac_to_air(np.array([1550.0]))[0] == 1550.0
    # vacuum > air above the boundary (the physical sign)
    assert wu.vac_to_air(np.array([15000.0]))[0] < 15000.0


# ── unit conversion + scale per the cited registry ────────────────────────────
def test_to_air_angstrom_spirou_vacuum_nm():
    w_nm = np.array([1000.0, 1500.0, 2400.0])
    air = foc.to_air_angstrom(w_nm, 'SPIRou')
    # ×10 then vac→air: offset is ~3–7 Å in the near-IR, air shorter than vacuum
    offset = air - w_nm * 10.0
    assert np.all(offset < 0) and np.all(np.abs(offset) > 2.0) and np.all(np.abs(offset) < 8.0)


def test_to_air_angstrom_uves_air_is_noop():
    w = np.array([4000.0, 6000.0, 9000.0])
    assert np.allclose(foc.to_air_angstrom(w, 'UVES'), w)


def test_to_air_angstrom_um_instrument_scales_by_1e4():
    # iSHELL is µm/vacuum: 2.3 µm → 23000 Å then vac→air (no band gate declared)
    air = foc.to_air_angstrom(np.array([2.3]), 'iSHELL')
    assert 22980.0 < air[0] < 23000.0


# ── the LOUD failures (the whole point) ───────────────────────────────────────
def test_unknown_instrument_raises_no_default():
    with pytest.raises(foc.FrameContractError, match='undeclared'):
        foc.to_air_angstrom(np.array([1500.0]), 'NOSUCH')


def test_unit_sanity_gate_catches_nm_as_angstrom():
    # SPIRou nm array fed as if already Å (≈1000–2400, i.e. ×10 too small) → RAISE
    with pytest.raises(foc.FrameContractError, match='outside the plausible air band'):
        foc.to_air_angstrom(np.array([1000.0, 1500.0, 2400.0]) / 10.0, 'SPIRou')


def test_unit_sanity_gate_catches_um_as_angstrom():
    # APOGEE H-band as µm-magnitude (≈1.5) instead of ~15000 Å → RAISE
    with pytest.raises(foc.FrameContractError):
        foc.to_air_angstrom(np.array([1.5, 1.6, 1.7]), 'APOGEE')


def test_all_finite_required_for_band_gate():
    with pytest.raises(foc.FrameContractError, match='non-finite'):
        foc.to_air_angstrom(np.array([np.nan, np.nan]), 'SPIRou')


# ── registry integrity: every entry is unit/scale-legal and CITED ─────────────
def test_registry_entries_are_legal_and_cited():
    for inst, conv in foc.WAVELENGTH_CONVENTION.items():
        assert conv.native_unit in foc._UNIT_TO_ANGSTROM, inst
        assert conv.native_scale in foc._NATIVE_SCALES, inst
        assert conv.citation and len(conv.citation) > 15, f"{inst} lacks a real citation"
        if conv.band_A is not None:
            lo, hi = conv.band_A
            assert lo < hi and conv.band_ref, inst


def test_verified_conventions_match_the_docs():
    # the ones verified against DRS docs this ticket (incl. the corrections)
    assert foc.WAVELENGTH_CONVENTION['HARPS'].native_scale == 'air'
    assert foc.WAVELENGTH_CONVENTION['UVES'].native_scale == 'air'
    assert foc.WAVELENGTH_CONVENTION['ESPRESSO'].native_scale == 'vacuum'   # corrected from 'air'
    assert foc.WAVELENGTH_CONVENTION['NIRPS'].native_scale == 'vacuum'
    assert foc.WAVELENGTH_CONVENTION['SPIRou'] == foc.WAVELENGTH_CONVENTION['SPIRou']
    assert foc.WAVELENGTH_CONVENTION['SPIRou'].native_unit == 'nm'
    assert foc.WAVELENGTH_CONVENTION['CRIRES+'].native_unit == 'nm'
    assert foc.WAVELENGTH_CONVENTION['iSHELL'].native_unit == 'um'
    assert foc.WAVELENGTH_CONVENTION['APOGEE'].native_scale == 'vacuum'


# ── WavelengthScale dataclass (mirrors VelocityFrame) ─────────────────────────
def test_wavelength_scale_validate_catches_a_contradiction():
    # a loader claiming SPIRou is air (it is vacuum) must fail loud
    with pytest.raises(foc.FrameContractError, match='cited convention'):
        foc.WavelengthScale('SPIRou', native_unit='nm', native_scale='air').validate()


def test_wavelength_scale_validate_passes_and_declares():
    ws = foc.WavelengthScale('SPIRou', native_unit='nm', native_scale='vacuum').validate()
    d = ws.declare()
    assert 'SPIRou' in d and 'vacuum→air' in d and 'nm→Å' in d


def test_wavelength_scale_provided_air_skips_vac_to_air():
    # NIRPS WAVE_AIR path: loader read the DRS air column → unit only, no vac→air
    ws = foc.WavelengthScale('NIRPS', native_unit='A', native_scale='vacuum',
                             provided_air=True, note='read WAVE_AIR column')
    air_col = np.array([10000.0, 15000.0, 19000.0])     # already air from the DRS
    out = ws.to_air_angstrom(air_col)
    assert np.allclose(out, air_col)                    # unchanged (no conversion re-applied)
    assert 'DRS-provided air column' in ws.declare()


def test_wavelength_scale_rejects_illegal_unit_or_scale():
    with pytest.raises(foc.FrameContractError):
        foc.WavelengthScale('UVES', native_unit='parsec', native_scale='air')
    with pytest.raises(foc.FrameContractError):
        foc.WavelengthScale('UVES', native_unit='A', native_scale='plasma')


# ── scale sign-check (the §C analogue of the §B BERV sign-check) ──────────────
def test_verify_vac_to_air_passes_on_known_line():
    # Mg II k 2796.3543 vac → 2795.528 air (Morton/NIST), the RYA-426 anchor
    resid = foc.verify_vac_to_air(2796.3543, 2795.528, tol_A=0.01)
    assert abs(resid) < 0.01


def test_verify_vac_to_air_raises_on_wrong_target():
    with pytest.raises(foc.FrameContractError, match='vac→air check failed'):
        foc.verify_vac_to_air(2796.3543, 2796.3543, tol_A=0.01)   # forgot to convert


# ── loader wiring: SPIRou now converts, UVES declares air ─────────────────────
def test_spirou_loader_declares_vacuum_nm_scale():
    from pipeline.loaders import spirou_loader as sl
    assert sl._SPIROU_SCALE.native_unit == 'nm' and sl._SPIROU_SCALE.native_scale == 'vacuum'
    sl._SPIROU_SCALE.validate()                          # consistent with the registry
    # the loader converts in the observed frame BEFORE BERV (documented order):
    converted = sl._SPIROU_SCALE.to_air_angstrom(np.array([2000.0]))   # 2000 nm
    assert converted[0] < 20000.0 and converted[0] > 19990.0           # ×10 + small vac→air


def test_uves_loader_declares_air_noop_scale():
    from pipeline.loaders import uves_loader as ul
    assert ul._UVES_SCALE.native_scale == 'air'
    ul._UVES_SCALE.validate()
    assert np.allclose(ul._UVES_SCALE.to_air_angstrom(np.array([5000.0])), [5000.0])

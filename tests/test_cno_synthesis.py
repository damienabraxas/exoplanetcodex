"""
tests/test_cno_synthesis.py
===========================
Tests for pipeline/cno_synthesis.py — the region-aware C/N/O synthesis engine
(RYA-237). These exercise the engine's logic WITHOUT running Turbospectrum
(synthesis is covered by the end-to-end solar VIS validate run): the diagnostic
registry, the multi-element fixed-abundance builder, the pluggable NLTE backends
(VIS LTE-by-design vs the loud-failing Amarsi grid), molecular-list coverage, the
preflight no-silent-fallback guards, and the solar-gate validation.
"""
import numpy as np
import pytest

from pipeline import cno_synthesis as cno


# ── Diagnostic registry ───────────────────────────────────────────────────────

class TestRegistry:
    def test_elements_are_cno(self):
        assert {d.element for d in cno.VIS_DIAGNOSTICS} == {'C', 'N', 'O'}

    def test_primary_band_per_element(self):
        primaries = {d.element for d in cno.VIS_DIAGNOSTICS if d.role == 'primary'}
        assert primaries == {'C', 'N', 'O'}

    def test_windows_inside_harps_range(self):
        for d in cno.VIS_DIAGNOSTICS:
            for lo, hi in d.windows_A:
                assert lo < hi
                assert cno.HARPS_VIS.wave_min_A <= lo
                assert hi <= cno.HARPS_VIS.wave_max_A

    def test_cn_depends_on_carbon(self):
        cn = next(d for d in cno.VIS_DIAGNOSTICS if d.key == 'CN_red')
        assert 'C' in cn.depends_on        # equilibrium coupling: N needs A(C)

    def test_oi_depends_on_carbon_and_pins_nickel(self):
        oi = next(d for d in cno.VIS_DIAGNOSTICS if d.key == 'OI_6300')
        assert 'C' in oi.depends_on        # [O I] tied to A(C) via CO
        assert 'Ni' in oi.pinned_blends    # joint [O I] + Ni I 6300.34 synthesis

    def test_no_red_arm_lines_present(self):
        # O I 7771-5 and N I 7468/8216 are red of the 6910 HARPS edge (351-gated).
        for d in cno.VIS_DIAGNOSTICS:
            for lo, hi in d.windows_A:
                assert hi <= 6910.0


# ── NLTE flags (VIS LTE-by-design, justified per physics) ─────────────────────

class TestNLTEFlags:
    def test_carbon_optical_lines_flagged_lte_assumed(self):
        for key in ('CI_5052', 'CI_5380'):
            d = next(x for x in cno.VIS_DIAGNOSTICS if x.key == key)
            assert d.nlte_flag == 'cI_vis_lte_assumed'   # revisit when RYA-359 lands

    def test_forbidden_oxygen_lte_insensitive(self):
        oi = next(d for d in cno.VIS_DIAGNOSTICS if d.key == 'OI_6300')
        assert oi.nlte_flag == 'lte_forbidden_insensitive'

    def test_molecular_bands_lte(self):
        for key in ('CH_Gband', 'CN_red', 'C2_Swan'):
            d = next(x for x in cno.VIS_DIAGNOSTICS if x.key == key)
            assert d.nlte_flag == 'lte_molecular_band'


# ── Pluggable NLTE backends ───────────────────────────────────────────────────

class TestNLTEBackends:
    def test_vis_lte_backend_applies_zero_correction(self):
        d = cno.VIS_DIAGNOSTICS[0]
        a_nlte, delta, flag, ref = cno.vis_lte_backend(d, 8.46, {})
        assert a_nlte == 8.46 and delta == 0.0
        assert flag == d.nlte_flag

    def test_vis_backend_does_not_loud_fail_without_grid(self):
        # The whole point of the run-scope clarification: VIS must NOT require the
        # RYA-359 C/O grid.
        d = next(x for x in cno.VIS_DIAGNOSTICS if x.key == 'CH_Gband')
        cno.vis_lte_backend(d, 8.5, {})        # no raise

    def test_amarsi_backend_loud_fails_until_grid_lands(self):
        d = cno.VIS_DIAGNOSTICS[0]
        with pytest.raises(NotImplementedError):
            cno.amarsi_grid_backend(d, 8.46, {})

    def test_both_backends_registered(self):
        assert set(cno.NLTE_BACKENDS) == {'lte_by_design', 'amarsi_grid'}


# ── Multi-element fixed-abundance builder ─────────────────────────────────────

class TestFixedAbundances:
    CODES = {'C': 6, 'N': 7, 'O': 8, 'Ni': 28}

    def test_all_elements_present(self):
        fa = cno._fixed_ab({'C': 8.46, 'N': 7.83, 'O': 8.69, 'Ni': 6.2}, self.CODES)
        assert set(fa['element']) == {'C', 'N', 'O', 'Ni'}
        assert set(fa['code']) == {6, 7, 8, 28}

    def test_ispec_scale_offset_applied(self):
        fa = cno._fixed_ab({'C': 8.46}, self.CODES)
        # A(X) -> log(N/Ntot) iSpec scale = A_X - 12.036
        assert fa['Abund'][0] == pytest.approx(8.46 - cno._ISPEC_SCALE_OFFSET)

    def test_varying_free_element_changes_only_its_row(self):
        base = cno._fixed_ab({'C': 8.46, 'O': 8.69}, self.CODES)
        bumped = cno._fixed_ab({'C': 8.76, 'O': 8.69}, self.CODES)
        oi = list(base['element']).index('O')
        assert base['Abund'][oi] == bumped['Abund'][oi]      # O unchanged
        ci = list(base['element']).index('C')
        assert bumped['Abund'][ci] - base['Abund'][ci] == pytest.approx(0.30)


# ── Region config + molecular coverage ────────────────────────────────────────

class TestRegion:
    def test_harps_vis_lsf_and_gates(self):
        r = cno.HARPS_VIS
        assert r.R == 115000
        assert r.telluric_correction_required is False     # optical
        assert r.nlte_backend == 'lte_by_design'

    def test_molecules_cover_real_bands(self):
        # CH (~430 nm), CN (~613 nm), C2 (~516 nm) all have .bsyn lists.
        for key in ('CH_Gband', 'CN_red', 'C2_Swan'):
            d = next(x for x in cno.VIS_DIAGNOSTICS if x.key == key)
            assert cno._molecules_cover(d.windows_A), key

    def test_molecules_do_not_cover_out_of_range(self):
        assert cno._molecules_cover(((200.0, 201.0),)) is False


# ── Preflight no-silent-fallback guards ───────────────────────────────────────

class TestPreflight:
    def test_telluric_required_region_loud_fails(self):
        ir = cno.RegionConfig(name='ir', instrument='CRIRES+', R=100000,
                              wave_min_A=9600, wave_max_A=23000,
                              telluric_correction_required=True,
                              nlte_backend='amarsi_grid')
        with pytest.raises(RuntimeError, match='telluric'):
            cno.preflight(ir, 'solar', ())     # empty diagnostics → reach telluric gate

    def test_unknown_star_broadening_loud_fails(self):
        with pytest.raises((KeyError, Exception)):
            cno.preflight(cno.HARPS_VIS, 'no_such_star_xyz', ())

    def test_uncovered_molecular_band_loud_fails(self):
        bogus = cno.Diagnostic(
            key='bogus', element='C', kind='molecular_band',
            windows_A=((200.0, 201.0),), use_molecules=True, role='primary',
            nlte_flag='lte_molecular_band', nlte_ref='x')
        with pytest.raises(FileNotFoundError):
            cno.preflight(cno.HARPS_VIS, 'solar', (bogus,))


# ── Solar-VIS gate validation logic ───────────────────────────────────────────

def _result(aC, aN, aO):
    co = round(10 ** (aC - aO), 3)
    r = cno.CNOResult(star_id='solar', region='vis')
    r.abundances = {'C': aC, 'N': aN, 'O': aO, 'C/O': co}
    r.uncertainty = {e: {'stat': 0.02, 'sys': 0.02, 'tot': 0.03}
                     for e in ('C', 'N', 'O')}
    r.per_band = [
        {'key': 'CH_Gband', 'A_X': aC}, {'key': 'CI_5052', 'A_X': aC + 0.02},
        {'key': 'CI_5380', 'A_X': aC - 0.02},
    ]
    r.iterations, r.converged = 3, True
    return r


class TestValidate:
    def test_on_target_passes(self):
        assert cno.validate_solar(_result(8.46, 7.83, 8.69)) is True

    def test_off_target_fails(self):
        # A(C) a full dex high → C and C/O gates fail.
        assert cno.validate_solar(_result(9.46, 7.83, 8.69)) is False

    def test_gates_cover_cno_and_ratio(self):
        assert set(cno.SOLAR_VIS_GATES) == {'C', 'N', 'O', 'C/O'}


if __name__ == '__main__':
    raise SystemExit(pytest.main([__file__, '-v']))

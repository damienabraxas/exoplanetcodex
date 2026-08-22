"""RYA-963 — the CRIRES+ stellar telluric driver.

These pin the invariants the ticket had to establish the hard way. Each one corresponds
to a failure that reported SUCCESS: the run finished, wrote products, and returned 0
while producing an uncorrected or unmodelled spectrum. That class is why the assertions
are in the code and not only in a comment.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline import crires_stellar_telluric as cst   # noqa: E402
from pipeline.crires_telluric import (                # noqa: E402
    TELLURIC_MOLECULES_BY_BAND, molecules_for_band)
from pipeline.telluric.esorex_runtime import (        # noqa: E402
    SUPPRESS_PREFIX, esorex_env, gdas_dirs)


# ── esorex runtime: one resolver, and the flag that must never be inherited ──
def test_suppress_prefix_is_explicit_and_true():
    """RYA-939: esorex takes --suppress-prefix from ~/.esorex/esorex.rc when not given,
    and the two machines disagree (Mac: no rc; Sirius source kit: FALSE). Under FALSE the
    recipe succeeds and writes out_0000.fits… so BEST_FIT_MODEL.fits is simply absent and
    the run reads as 'failed (rc=0)'."""
    assert SUPPRESS_PREFIX == "--suppress-prefix=TRUE"


def test_gdas_dirs_include_the_registered_eso_prefix():
    """The retriever globbed only Homebrew and /usr/share/esopipes, so on Sirius — where
    the ESO source kit ships the same Paranal tarball — GDASUnavailable would have fired
    on a profile sitting on the disk. Loud-failing for the wrong reason is still wrong."""
    from pipeline.telluric.esorex_runtime import eso_root
    pats = gdas_dirs()
    root = str(eso_root())
    assert any(p.startswith(root) for p in pats), (
        f"the registered eso_pipelines root {root} must be searched")
    assert pats[0].startswith(root), "the registered root must be searched FIRST"
    assert any('homebrew' in p for p in pats), "the Mac install must stay reachable too"


def test_esorex_env_adds_the_kit_lib_dir(tmp_path):
    (tmp_path / 'bin').mkdir()
    (tmp_path / 'lib').mkdir()
    exe = tmp_path / 'bin' / 'esorex'
    exe.write_text('#!/bin/sh\n')
    env = esorex_env(str(exe), base={'PATH': '/usr/bin'})
    assert str(tmp_path / 'bin') in env['PATH']
    assert str(tmp_path / 'lib') in env['LD_LIBRARY_PATH']


# ── molecule sets ────────────────────────────────────────────────────────────
def test_every_crires_band_declares_its_molecules():
    for band in 'YJHK':
        assert molecules_for_band(band), band


def test_unknown_band_is_loud_not_defaulted():
    """The K-band set is not a safe default for a bluer band; fitting a molecule with no
    band in the window is a free parameter with nothing to constrain it."""
    with pytest.raises(ValueError):
        molecules_for_band('Z')


def test_declared_band_molecules_all_have_a_band_table_entry():
    for band, mols in TELLURIC_MOLECULES_BY_BAND.items():
        for m in mols:
            assert m in cst.MOLECULE_BANDS_UM, f"{band}: {m} has no absorption bands listed"


def test_o2_is_fitted_in_J_but_not_in_K():
    """The 1.27 µm O2 band is in J and nowhere near K; the band sets must reflect that
    rather than carrying one molecule list everywhere."""
    assert 'O2' in molecules_for_band('J')
    assert 'O2' not in molecules_for_band('K')


# ── the fit-moved invariant ──────────────────────────────────────────────────
def _frozen_fit():
    """What molecfit actually returned when WLC_CONST put the model 135 Å off the data:
    success, all products written, every column at its prior with uncertainty 0."""
    return {'initial_chi2': (23229447.81, -1.0), 'best_chi2': (23229447.81, -1.0),
            'rel_mol_col_H2O': (1.0, 0.0), 'rel_mol_col_CH4': (1.0, 0.0),
            'rel_mol_col_CO2': (1.0, 0.0), 'rel_mol_col_CO': (1.0, 0.0)}


def test_a_fit_that_never_moved_is_refused():
    with pytest.raises(RuntimeError, match='no molecular column was actually fitted'):
        cst.assert_fit_moved(_frozen_fit(), ['H2O', 'CH4', 'CO2', 'CO'])


def test_chi2_that_did_not_improve_is_refused():
    frozen = _frozen_fit()
    frozen['rel_mol_col_H2O'] = (1.0, 0.02)      # a column moved…
    with pytest.raises(RuntimeError, match='chi2 never improved'):
        cst.assert_fit_moved(frozen, ['H2O'])    # …but the model still did not respond


def test_a_real_fit_passes_and_reports_which_columns_moved():
    good = {'initial_chi2': (1.9e7, -1.0), 'best_chi2': (15101.5, -1.0),
            'rel_mol_col_H2O': (1.1058, 0.0169), 'rel_mol_col_CH4': (0.9013, 0.0162),
            'rel_mol_col_CO2': (1e-05, 0.0)}
    out = cst.assert_fit_moved(good, ['H2O', 'CH4', 'CO2'])
    assert out['fitted_columns'] == ['H2O', 'CH4']   # CO2 pegged at its floor, unc 0


# ── stellar mask + fit planning ──────────────────────────────────────────────
def test_stellar_intervals_are_merged_and_rv_shifted():
    lo, hi = 15000.0, 15100.0
    at_rest = cst.stellar_line_intervals(lo, hi, rv_kms=0.0)
    shifted = cst.stellar_line_intervals(lo, hi, rv_kms=-40.0)
    assert at_rest, "the canonical solar list must reach the H band"
    for a, b in at_rest:
        assert b > a
    for (a, b), (c, d) in zip(at_rest, at_rest[1:]):
        assert c > b, "intervals must be disjoint after merging"
    # -40 km/s at 1.5 µm is -0.20 Å: a real, sub-mask-width shift
    assert not np.allclose([x[0] for x in at_rest][:5], [x[0] for x in shifted][:5])


def test_stellar_mask_selects_only_inside_intervals():
    w = np.array([100.0, 200.0, 300.0])
    assert cst.stellar_mask(w, [(150.0, 250.0)]).tolist() == [False, True, False]


def test_informative_fraction_ignores_saturated_cores():
    """Scoring on 'below 0.97' alone picks the saturated chip, whose cores carry almost
    no information about the column (dI/dN → 0) while dominating chi2."""
    w = np.linspace(10000.0, 10060.0, 400)
    flat = np.ones_like(w)
    saturated = flat.copy(); saturated[100:300] = 1e-3
    informative = flat.copy(); informative[100:300] = 0.5
    stellar = np.zeros_like(w, dtype=bool)
    assert (cst._absorbed_fraction(w, informative, stellar)
            > cst._absorbed_fraction(w, saturated, stellar))


def test_molecule_band_overlap():
    assert cst._overlaps(2.29, 2.35, cst.MOLECULE_BANDS_UM['CO'])
    assert not cst._overlaps(1.95, 2.05, cst.MOLECULE_BANDS_UM['CO'])


# ── the sets registry ────────────────────────────────────────────────────────
def test_unknown_set_is_loud():
    with pytest.raises(ValueError):
        cst.resolve_set('no_such_set')


def test_alpha_cen_set_declares_its_claim_and_its_gate():
    rec = cst.resolve_set('alpha_cen_a_crires')
    assert rec['holding_id'] == 'alpha_cen_a_crires_plus'
    assert rec['claimed_star'] == 'A'        # a CLAIM, tested by the ID gate
    assert rec['id_gate'] == 'acen_ab'
    assert rec['epoch'] == '2022-04-15'


def test_singleton_id_gate_is_not_silently_the_acen_orbit():
    """RYA-965 points this driver at tau Ceti / eps Eri / 55 Cnc, which are singletons.
    An unimplemented gate must say so, not fall through to the AB orbit."""
    with pytest.raises(NotImplementedError):
        cst.identify_star(None, {}, id_gate='rya964_alias')


# ── frame table ──────────────────────────────────────────────────────────────
def test_frame_table_keeps_overlapping_orders_as_separate_chips():
    """CRIRES+ echelle orders OVERLAP toward the blue — Y1029's ord9/det3 and ord8/det1
    share 28 Å. The K band RYA-373 worked in has no such overlap, which is why the first
    cut of this driver asserted that orders tile. Chips are therefore laid end to end in
    starting-wavelength order, NOT globally sorted: a global sort interleaves the two
    chips' pixels and the transmission model would be mapped back across the wrong chip
    boundary in the overlap region."""
    from pipeline.crires_telluric import CriresFrame, CriresSegment

    def seg(lo, hi, o, d, n=100):
        w = np.linspace(lo, hi, n)
        return CriresSegment(order=o, detector=d, wave_A=w, flux=np.ones(n),
                             err=np.ones(n), qual=np.zeros(n, int))

    def frame(segments):
        return CriresFrame(path=Path('x.fits'), wlen_id='Y1029', band='Y', mjd=0.0,
                           date_obs='', ra=0.0, dec=0.0, snr=1.0, specsys='TOPOCENT',
                           fluxcal='UNCALIBRATED', wmin_nm=0.0, wmax_nm=1.0,
                           segments=segments)

    # the real Y1029 overlap, in the order the loader emits them (order desc = wave asc)
    overlapping = frame([seg(9659.52, 9724.42, 8, 1), seg(9629.76, 9687.52, 9, 3)])
    w, f, e, ix = cst._frame_table(overlapping)
    # chips are reordered by starting wavelength, and each stays ONE contiguous block
    assert ix[0] == 1 and ix[-1] == 0, "chips must be laid out blue-first"
    edges = np.flatnonzero(np.diff(ix) != 0)
    assert len(edges) == 1, "each chip must remain a single contiguous row block"
    # and each block is internally monotonic, which is what the FITS extensions need
    for a, b in ((0, edges[0] + 1), (edges[0] + 1, len(w))):
        assert np.all(np.diff(w[a:b]) > 0)
    # the concatenation is deliberately NOT globally monotonic across the overlap
    assert np.any(np.diff(w) < 0)


def test_frame_table_refuses_a_chip_that_is_not_monotonic():
    """A chip must be one monotonic block; that is what makes it one FITS extension and
    what lets the per-chip transmission be mapped back by row slice."""
    from pipeline.crires_telluric import CriresFrame, CriresSegment
    w = np.array([100.0, 101.0, 101.0, 102.0])       # a duplicate wavelength
    bad = CriresFrame(path=Path('x.fits'), wlen_id='K2192', band='K', mjd=0.0,
                      date_obs='', ra=0.0, dec=0.0, snr=1.0, specsys='TOPOCENT',
                      fluxcal='UNCALIBRATED', wmin_nm=0.0, wmax_nm=1.0,
                      segments=[CriresSegment(order=1, detector=1, wave_A=w,
                                              flux=np.ones(4), err=np.ones(4),
                                              qual=np.zeros(4, int))])
    with pytest.raises(RuntimeError, match='non-increasing'):
        cst._frame_table(bad)


# ── CCF railing ──────────────────────────────────────────────────────────────
def test_ccf_peak_reports_a_railed_maximum_rather_than_returning_the_edge():
    """A peak on the first or last velocity sample is not a measurement — the true
    maximum is outside the grid. A frame whose transmission was mostly NaN produced
    exactly this: -80.0 km/s, the grid edge, which the orbit test then read as 'outside
    the alpha Cen bounds' and returned the confident, wrong verdict NOT-ALPHA-CEN."""
    v = np.arange(-60.0, 60.1, 0.25)
    rising = np.linspace(0.0, 1.0, v.size)          # maximum at the last sample
    out = cst._ccf_peak(v, rising)
    assert out['railed'] is True
    assert 'edge' in out['reason']


def test_ccf_peak_interpolates_a_real_maximum():
    v = np.arange(-60.0, 60.1, 0.25)
    ccf = np.exp(-0.5 * ((v + 25.5) / 3.0) ** 2)
    out = cst._ccf_peak(v, ccf)
    assert out['railed'] is False
    assert abs(out['rv_kms'] + 25.5) < 0.2
    assert out['contrast'] > 1.0


def test_ccf_peak_on_a_flat_ccf_has_no_contrast():
    v = np.arange(-10.0, 10.1, 0.25)
    out = cst._ccf_peak(v, np.zeros_like(v))
    assert out['contrast'] == 0.0


def test_the_b_control_set_is_declared_and_matches_a_on_epoch_and_gate():
    """The α Cen B K2192 frame is the POSITIVE CONTROL for the A star-ID: same setting,
    same night, 16 minutes apart. A control only discriminates if it goes through the
    identical path, so the two sets must agree on epoch and on the id gate — a control
    that differs in either proves nothing about the thing it controls."""
    a = cst.resolve_set('alpha_cen_a_crires')
    b = cst.resolve_set('alpha_cen_b_crires')
    assert b['claimed_star'] == 'B' and a['claimed_star'] == 'A'
    assert b['epoch'] == a['epoch'], "same night, or it is not a matched control"
    assert b['id_gate'] == a['id_gate'], "same rule, or it is not the same test"
    assert b['holding_id'] != a['holding_id']


def test_the_two_orbit_branches_are_far_enough_apart_to_decide():
    """The A/B split is only decidable if the branches separate by more than the RV
    tolerance. At the 2022-04-15 epoch they are 6.75 km/s apart against a 2.5 km/s
    tolerance — declare that in advance rather than discovering it from the answer."""
    from astropy.time import Time
    from pipeline.acen_orbit import predicted_rv
    mjd = Time('2022-04-15T04:00:00', format='isot', scale='utc').mjd
    p = predicted_rv(mjd)
    mod = cst._rya423_verdict()
    assert abs(p['delta_AB']) > 2 * mod.RV_TOL, (
        f"A and B are only {p['delta_AB']:.2f} km/s apart at this epoch; the "
        f"{mod.RV_TOL} km/s match tolerance cannot separate them")


# ── shared-box safety ────────────────────────────────────────────────────────
def test_a_killed_recipe_is_not_reported_as_a_failed_fit():
    """rc < 0 is a SIGNAL. Reporting it as 'esorex FAILED' sends the next reader hunting
    a fit problem that does not exist — the process never got to finish. -9 on Sirius is
    the OOM killer; dmesg confirmed it for the Y1029 model at 13.7 GB anon RSS."""
    class P:
        returncode = -9
        stdout = ''
    with pytest.raises(RuntimeError, match='KILLED by SIGKILL'):
        cst._require_product(P(), Path('/nonexistent'), 'X.fits', 'y')


def test_a_killed_recipe_names_the_memory_cap_as_the_likely_cause():
    class P:
        returncode = -9
        stdout = ''
    with pytest.raises(RuntimeError, match='OOM killer'):
        cst._require_product(P(), Path('/nonexistent'), 'X.fits', 'y')


def test_the_memory_cap_is_set_and_leaves_headroom_on_a_shared_box():
    """Sirius has 15 GB and is shared. The cap exists so that a runaway of OURS fails
    with our name on it, instead of the kernel choosing a victim — which on 2026-08-22
    could as easily have been another session's 2.5-hour synthesis."""
    assert 0 < cst._MEM_CAP_GIB <= 12


def test_relative_window_floor_drops_windows_that_carry_little_information():
    """An absolute floor alone kept Y1029 windows at f=0.058 and 0.062 beside a best of
    0.267 — regions where the non-stellar absorption is mostly stellar residual, and
    which cost the most per iteration precisely where they help least."""
    assert cst._RELATIVE_TELLURIC_FLOOR > 0
    best = 0.267
    floor = max(cst._MIN_TELLURIC_FRAC, cst._RELATIVE_TELLURIC_FLOOR * best)
    assert 0.062 < floor and 0.128 >= floor, (
        "the floor must drop Y1029's two noise windows and keep its two real ones")


# ── measured resolving power (the referee molecfit is never told) ────────────
def _fc_with_lsf(fwhm_px, lo_A=11065.9, hi_A=11125.9, step_A=0.03434):
    from pipeline.crires_telluric import CriresFrame
    n = int((hi_A - lo_A) / step_A) + 200
    wave = lo_A - 100 * step_A + step_A * np.arange(n)
    fr = CriresFrame(path=Path('x.fits'), wlen_id='Y1029', band='Y', mjd=0.0,
                     date_obs='', ra=0.0, dec=0.0, snr=1.0, specsys='TOPOCENT',
                     fluxcal='UNCALIBRATED', wmin_nm=0.0, wmax_nm=1.0, segments=[])
    return cst.FrameCorrection(
        frame=fr, wave_A=wave, flux_raw=np.ones(n), err=np.ones(n),
        mtrans=np.ones(n), flux_corr=np.ones(n), seg_index=np.zeros(n, int),
        windows=[{'lo_A': lo_A, 'hi_A': hi_A, 'order': 1, 'detector': 2,
                  'absorbed_frac': 0.3}],
        fit={'gaussfwhm': (fwhm_px, 0.026)})


def test_resolving_power_recovers_an_instrument_property_it_was_never_given():
    """molecfit is told nothing about R; it fits a kernel width in pixels. Converting
    that back must land near CRIRES+'s real resolving power -- nominal 100,000 for the
    0.2in slit, 86,000 as the project's own working value (RYA-373/952)."""
    out = cst.measured_resolving_power(_fc_with_lsf(4.037))
    assert 60_000 < out['R'] < 120_000, out
    assert abs(out['fwhm_A'] - 4.037 * 0.03434) < 1e-6


def test_a_collapsed_lsf_yields_no_resolving_power_rather_than_a_huge_one():
    """FIT_RES_GAUSS=FALSE does not hold the kernel at a width, it DISABLES convolution
    (RYA-931). A zero width would divide to an infinite R, which would read as a
    spectacular instrument rather than as an unfitted kernel."""
    out = cst.measured_resolving_power(_fc_with_lsf(0.0))
    assert not np.isfinite(out['R'])
    assert 'not fitted' in out['reason']


# ── telluric-anchor closure (RYA-373's rule, reused not re-declared) ─────────
def test_the_closure_bound_is_rya373s_constant_not_a_new_one():
    from pipeline.crires_telluric import _TELLURIC_CLOSURE_MAX
    assert _TELLURIC_CLOSURE_MAX == 3.0


def test_a_large_zero_point_is_a_failed_anchor_not_a_zero_point():
    """K2148 returned a -12.78 km/s 'zero-point' at a perfectly ordinary 3.97-sigma CCF
    contrast, while the other five settings closed within 1.86. Contrast does not catch
    it; the physics does -- tellurics are at topocentric rest, so the anchor must close
    to ~0, and subtracting -12.78 fabricates a velocity shift in the stellar RV."""
    from pipeline.crires_telluric import _TELLURIC_CLOSURE_MAX
    assert abs(-12.782) > _TELLURIC_CLOSURE_MAX      # K2148: refused
    for zp in (-1.857, -0.945, -0.031, -0.013, 0.766):
        assert abs(zp) <= _TELLURIC_CLOSURE_MAX      # the five that closed


# ── RYA-973: multi-night sets, and collisions the alpha Cen shape could not expose ──
def test_tau_ceti_set_is_a_multi_directory_singleton_with_no_date_filter():
    """tau Ceti's CRIRES+ frames are split across CRIRES/ and CRIRESPlus/ by how they
    ARRIVED, not by what they are — 25 in the first (all cr2res_obs_nodding, including
    the only Y1029/J1232/H1582/H1559/K2192 frames) and 4 in the second. A set keyed to
    one directory sees a third of the star.

    And no `epochs` filter: identity is verified PER FRAME by the intake gate, which is
    what actually protects us. A date filter never did — it was only ever a proxy for
    "is this the right target", and the gate answers that directly."""
    rec = cst.resolve_set('tau_ceti_crires')
    assert rec['claimed_star'] == 'tau_ceti'
    assert rec['id_gate'] == 'singleton_astrometry', (
        "a singleton has no close pair to split; the alpha Cen orbit rule must not apply")
    assert len(rec['dirs']) == 2
    assert 'epochs' not in rec, "identity is per-frame, not per-date"


def test_the_astrometry_and_catalogue_ids_disagree_as_raw_strings():
    """`audit_crires` identifies this star as `tau_cet` (from the CRIRES astrometry
    reference) while system_catalog / stars.yaml / the holdings registry call it
    `tau_ceti`. Compared raw, the two identity routes disagree on EVERY tau Ceti frame,
    so the check meant to catch a mislabelled star would cry wolf on all of them. Both
    sides must go through the one alias lookup."""
    from pipeline.star_id import resolve_star
    assert 'tau_cet' != 'tau_ceti'                       # the raw-string trap
    assert resolve_star('tau_cet') == resolve_star('tau_ceti') == 'tau_ceti'


def test_every_astrometry_reference_id_resolves_through_the_alias_lookup():
    """Any id the astrometry reference can emit must be resolvable, or the identity
    cross-check silently degrades to 'UNRESOLVED vs something' for that star."""
    import csv
    from pipeline.star_id import resolve_star
    from config.constants import codex_root
    path = Path(codex_root('repo')) / 'data' / 'reference' / 'crires_target_astrometry.csv'
    unresolved = []
    for row in csv.DictReader(open(path)):
        sid = row['star_id']
        if resolve_star(sid) == 'UNRESOLVED':
            unresolved.append(sid)
    # tau_boo is UNRESOLVED on purpose (no star_params_key yet, RYA-957)
    assert set(unresolved) <= {'tau_boo'}, f"unresolvable astrometry ids: {unresolved}"


def test_work_and_product_names_separate_frames_of_the_same_setting():
    """tau Ceti's four frames are ALL K2148, two per night. A work dir or product named
    for setting+date alone collides — two of them fitted against DIFFERENT nights'
    atmospheres, with only the last surviving. alpha Cen could not expose this: it had
    exactly one frame per setting."""
    stems = ['ADP.2025-05-10T13:29:03.307', 'ADP.2025-05-10T13:29:03.310',
             'ADP.2025-05-10T15:25:08.684', 'ADP.2025-05-10T15:25:08.687']
    dates = ['2022-01-06', '2022-01-06', '2022-01-16', '2022-01-16']
    setting_only = {f"K2148_{d}" for d in dates}
    per_frame = {f"tau_ceti_crires_K2148_{d}_{s}" for d, s in zip(dates, stems)}
    assert len(setting_only) == 2, "setting+date collapses four frames onto two names"
    assert len(per_frame) == 4, "set+setting+date+frame keeps all four distinct"


# ── RYA-973: the standing intake identity procedure ─────────────────────────
def test_star_id_namespace_guard_passes_on_the_committed_tree():
    from pipeline.intake_identity import assert_star_id_namespace
    out = assert_star_id_namespace()
    ids = out['sources']['data/reference/crires_target_astrometry.csv']
    assert 'tau_ceti' in ids and 'tau_cet' not in ids, (
        "the astrometry reference must use the canonical star_params_key")


def test_star_id_namespace_guard_actually_fails_on_drift(tmp_path):
    """Prove the tripwire RED before trusting it. A guard that cannot fail is not a
    guard — and this exact drift sat undetected because nothing joined the two
    registries until the first consumer compared them."""
    import csv
    from pipeline.intake_identity import assert_star_id_namespace
    from config.constants import codex_root
    src = Path(codex_root('repo')) / 'data' / 'reference' / 'crires_target_astrometry.csv'
    dst = tmp_path / 'data' / 'reference'
    dst.mkdir(parents=True)
    rows = list(csv.DictReader(open(src, newline='')))
    for r in rows:
        if r['star_id'] == 'tau_ceti':
            r['star_id'] = 'tau_cet'                      # re-introduce the drift
    with open(dst / 'crires_target_astrometry.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    with pytest.raises(AssertionError, match='namespace drift'):
        assert_star_id_namespace(repo_root=tmp_path)


def test_exemptions_carry_a_stated_reason():
    """An exemption list is where a defect hides if entries can be added silently."""
    from pipeline.intake_identity import NAMESPACE_EXEMPT
    assert set(NAMESPACE_EXEMPT) == {'tau_boo', 'vesta'}
    for sid, why in NAMESPACE_EXEMPT.items():
        assert len(why) > 30, f"{sid} exemption needs a stated reason, not a bare entry"


def test_singleton_gate_delegates_rather_than_reimplementing():
    """A second copy of the identity rule is how the two routes drifted apart. The
    driver must call the standing procedure."""
    import inspect
    src = inspect.getsource(cst.identify_singleton)
    assert 'identify_at_intake' in src
    assert 'load_astrometry' not in src, "re-implementing the astrometry route here"


# ── RYA-973: the cited radial velocity, and what it may referee ─────────────
def _target_ref():
    import csv
    from config.constants import codex_root
    path = Path(codex_root('repo')) / 'data' / 'reference' / 'crires_target_astrometry.csv'
    return {r['star_id']: r for r in csv.DictReader(open(path))}


def test_every_reference_rv_carries_its_own_bibcode():
    """An RV quoted without a source is a number from memory. This file exists to be the
    cited referee, so a value without a bibcode is not usable as one."""
    for sid, row in _target_ref().items():
        if row.get('rv_kms'):
            assert row.get('rv_bibcode'), f"{sid} has an RV with no source"


def test_tau_ceti_rv_reference_is_the_gaia_standard():
    """tau Ceti is a Gaia RADIAL-VELOCITY STANDARD (Soubiran+2018, A&A 616 A7), quoted to
    0.0002 km/s. That is what makes it a sharp external check on a K-band mask CCF —
    and what makes the K2148 anchor failure unmistakable: applying the railed -20 km/s
    'zero-point' would land ~20 km/s from a value known to four decimals."""
    row = _target_ref()['tau_ceti']
    assert row['rv_bibcode'] == '2018A&A...616A...7S'
    assert abs(float(row['rv_kms']) - (-16.597)) < 1e-3
    assert float(row['e_rv_kms']) < 0.01
    assert row['rv_quality'] == 'A'


def test_a_catalogue_rv_is_epoch_ambiguous_for_a_binary():
    """alpha Cen A's catalogued RV lies OUTSIDE the Kervella orbit bounds, which is the
    clearest possible demonstration that a single catalogue value cannot be compared to
    an epoch measurement for a binary component without establishing its epoch. tau Ceti,
    a singleton, has no such problem — which is why it referees cleanly and alpha Cen
    does not."""
    from pipeline.acen_orbit import rv_bounds
    lo, hi = rv_bounds()
    rv_a = float(_target_ref()['alpha_cen_a']['rv_kms'])
    assert not (lo <= rv_a <= hi), (
        f"alpha Cen A catalogue RV {rv_a} now sits inside the orbit bounds "
        f"[{lo:.2f},{hi:.2f}] — re-examine whether it can referee the branch question")
    rv_tau = float(_target_ref()['tau_ceti']['rv_kms'])
    assert -20 < rv_tau < -13, "a singleton's catalogue RV is directly comparable"


def test_a_code_fault_aborts_instead_of_becoming_a_per_frame_finding():
    """RYA-973: a module refactored into existence locally but never synced to the
    compute host made every frame run its full molecfit fit and THEN fail on the import.
    80 minutes of compute, and a report reading 'failed 6' as if the DATA were at fault.
    A code fault must abort on the first frame, where it is cheap and unambiguous."""
    assert ImportError in cst._CODE_FAULTS
    assert AttributeError in cst._CODE_FAULTS
    assert NameError in cst._CODE_FAULTS
    # and the genuinely per-frame conditions must NOT be in there
    assert RuntimeError not in cst._CODE_FAULTS
    assert ValueError not in cst._CODE_FAULTS


def test_sigabrt_and_sigkill_are_reported_as_different_things():
    """They mean different things and the distinction is operational. SIGKILL = the
    KERNEL chose this process, and on a shared box that choice could have fallen on
    someone else's job. SIGABRT under our RLIMIT_AS cap = CPL failing an allocation
    ('failed to allocate 16 bytes' — sixteen bytes means address space was exhausted,
    not that the request was large), i.e. our runaway failing with our own name on it.
    Neither is a fit failure."""
    class P:
        stdout = ''
        def __init__(self, rc): self.returncode = rc
    with pytest.raises(RuntimeError, match='kernel OOM killer'):
        cst._require_product(P(-9), Path('/nonexistent'), 'X.fits', 'y')
    with pytest.raises(RuntimeError, match='RLIMIT_AS cap'):
        cst._require_product(P(-6), Path('/nonexistent'), 'X.fits', 'y')
    # and both must steer toward a smaller problem, not a bigger ceiling
    for rc in (-9, -6):
        try:
            cst._require_product(P(rc), Path('/nonexistent'), 'X.fits', 'y')
        except RuntimeError as e:
            assert 'FEWER fit windows' in str(e)


# ── RYA-993: a well-mixed column is not a free parameter ────────────────────
def _best(**cols):
    b = {'initial_chi2': (1e7, -1.0), 'best_chi2': (1e4, -1.0)}
    for k, v in cols.items():
        b[f'rel_mol_col_{k}'] = (v, 0.01)
    return b


def test_a_runaway_high_column_is_flagged():
    """tau Ceti H1559 fitted CH4 = 22.711 where H1582 on the SAME NIGHT got 0.986. CH4 is
    well mixed; 22.7x atmospheric is not a column, it is a runaway parameter, and it
    paints a forest of weak lines the spectrum does not have."""
    out = cst.check_well_mixed_columns(_best(H2O=0.97, CO2=1.10, CH4=22.711),
                                       ['H2O', 'CO2', 'CH4'],
                                       {'H2O': True, 'CO2': True, 'CH4': True})
    assert not out['passed']
    f = [x for x in out['flagged'] if x['molecule'] == 'CH4'][0]
    assert f['direction'] == 'runaway-high'
    assert 'ABOVE continuum' in f['consequence']


def test_a_column_pegged_at_zero_is_flagged_AND_is_the_dangerous_one():
    """🔴 The asymmetry is the point. A HIGH column over-paints and the D1 residual
    catches it. A column at ZERO makes the model OMIT that molecule, leaving its real
    absorption uncorrected — and the residual stays small because H2O dominates the
    scored pixels. Six frames passed the gate in exactly that state."""
    out = cst.check_well_mixed_columns(_best(H2O=0.90, CO2=1.18, CH4=0.0),
                                       ['H2O', 'CO2', 'CH4'],
                                       {'H2O': True, 'CO2': True, 'CH4': True})
    assert not out['passed']
    f = [x for x in out['flagged'] if x['molecule'] == 'CH4'][0]
    assert f['direction'] == 'pegged-at-floor'
    assert 'UNCORRECTED' in f['consequence']


def test_water_is_deliberately_NOT_well_mixed():
    """Precipitable water genuinely varies by an order of magnitude — tau Ceti saw
    1.93 mm and 13-23 mm ten days apart — which is exactly why H2O is the column worth
    fitting freely. Bounding it would reject real weather."""
    assert 'H2O' not in cst.WELL_MIXED
    out = cst.check_well_mixed_columns(_best(H2O=5.0), ['H2O'], {'H2O': True})
    assert out['passed']


def test_a_molecule_held_at_its_prior_is_not_judged():
    """FIT_MOLEC=0 means the column was never fitted; it is the profile value, not a
    measurement, and has nothing to run away with."""
    out = cst.check_well_mixed_columns(_best(CH4=1.0), ['CH4'], {'CH4': False})
    assert out['passed'] and out['checked'] == []


def test_plausible_columns_pass():
    out = cst.check_well_mixed_columns(_best(CO2=1.10, CH4=0.99, CO=1.13),
                                       ['CO2', 'CH4', 'CO'],
                                       {'CO2': True, 'CH4': True, 'CO': True})
    assert out['passed'] and len(out['checked']) == 3


def test_the_bound_is_generous_enough_not_to_reject_real_variation():
    """Seasonal and altitude variation in CO2/CH4/CO is percent-level. The bound is
    half-to-double so that only a runaway trips it, never real air."""
    assert cst.WELL_MIXED_LO <= 0.8 and cst.WELL_MIXED_HI >= 1.3


# ── RYA-993 part 2: a runaway column is HELD at its physical prior, not re-thresholded ──

def test_o2_is_refereed_as_well_mixed():
    """O2 is 20.95% and constant; a fitted O2 column of 8x is an optimiser artefact.

    It is the dominant non-water opacity in the CRIRES+ Y setting, so leaving it out of
    the referee would let exactly the tau Ceti H1559 failure recur one band over.
    """
    assert 'O2' in cst.WELL_MIXED
    out = cst.check_well_mixed_columns(
        {'rel_mol_col_O2': (8.0, 0.4)}, ['O2'], {'O2': True})
    assert not out['passed']
    assert out['flagged'][0]['molecule'] == 'O2'
    assert out['flagged'][0]['direction'] == 'runaway-high'


def test_the_h1559_window_cleared_every_a_priori_threshold():
    """🔴 The a-priori window fix CANNOT work, and this is the frame that proves it.

    tau Ceti H1559's CH4 window scored an absorbed fraction of 0.175 against a relative
    floor of 0.4 x 0.412 = 0.165. It cleared the floor. CH4 still ran to 22.7x. The
    window score is TOTAL absorption -- at 1.5 um that is nearly all H2O -- so it says
    nothing about how much of the window belongs to CH4. Raising the threshold is not a
    fix; catching the runaway after the fit is.
    """
    windows = [0.334, 0.412, 0.339, 0.246, 0.149, 0.175]
    relative_floor = 0.4 * max(windows)
    assert 0.175 >= relative_floor          # the CH4 window PASSED
    assert 0.149 < relative_floor           # only CO2's came in by the exemption
    # ...and CO2, the one the exemption admitted, fitted perfectly well:
    assert cst.check_well_mixed_columns(
        {'rel_mol_col_CO2': (1.100, 0.02)}, ['CO2'], {'CO2': True})['passed']


def test_hold_clears_the_fit_flag_for_the_runaway_molecule_only():
    """Holding CH4 must not silently stop fitting H2O, which genuinely does vary."""
    molecules = ['H2O', 'CO2', 'CH4']
    fit_flags = {m: True for m in molecules}
    wm = cst.check_well_mixed_columns(
        {'rel_mol_col_H2O': (1.4, 0.1), 'rel_mol_col_CO2': (1.1, 0.02),
         'rel_mol_col_CH4': (22.711, 3.0)}, molecules, fit_flags)
    held = [f['molecule'] for f in wm['flagged']]
    assert held == ['CH4']
    refit = {m: (fit_flags[m] and m not in held) for m in molecules}
    assert refit == {'H2O': True, 'CO2': True, 'CH4': False}


def test_a_held_column_is_declared_modelled_not_measured():
    """WMHELD exists so nobody reports a prior as a measurement of the night's air."""
    src = Path(cst.__file__).read_text()
    assert "h['WMHELD']" in src
    assert 'MODELLED, not measured' in src

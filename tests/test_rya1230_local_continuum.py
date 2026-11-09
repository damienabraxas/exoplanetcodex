"""RYA-1230: the standing per-band continuum rule, model-guided estimator."""
import numpy as np
import pytest

from pipeline import local_continuum as lc


def _grid():
    return np.arange(5000.0, 5012.0, 0.01)


def test_recovers_a_pure_normalisation_error():
    w = _grid()
    model = np.ones_like(w)
    obs = 0.993 * model
    rec, cont = lc.fit_model_guided(w, obs, 5006.0, w, model, exclude_half_width_A=0.6)
    assert rec.applied
    assert rec.level_at_centre == pytest.approx(0.993, abs=1e-6)
    assert np.allclose(obs / cont, 1.0)


def test_a_modelled_forest_is_not_divided_out():
    """The defect that moved A(C) VIS by -0.14: an absorber the SYNTHESIS also carries
    (a line wing, a molecular forest) must not be read as a continuum error."""
    w = _grid()
    forest = 0.988 - 0.004 * np.sin(w * 0.7)               # a wing: never reaches 1.0
    obs = forest.copy()                                      # atlas perfectly normalised
    rec, cont = lc.fit_model_guided(w, obs, 5006.0, w, forest, exclude_half_width_A=0.6)
    assert rec.level_at_centre == pytest.approx(1.0, abs=1e-6)
    env, _ = lc.fit(w, obs, 5006.0, apply=False)
    assert env.level_at_centre < 0.999                        # the envelope would have divided it


def test_the_fitted_line_does_not_set_the_level():
    w = _grid()
    model = np.ones_like(w)
    line = 1.0 - 0.3 * np.exp(-0.5 * ((w - 5006.0) / 0.05) ** 2)
    obs = 0.995 * line                                        # line deeper than modelled is irrelevant
    rec, _ = lc.fit_model_guided(w, obs, 5006.0, w, model * line ** 0.5,
                                 exclude_half_width_A=0.6)
    assert rec.level_at_centre == pytest.approx(0.995, abs=2e-4)


def test_a_level_beyond_the_bound_is_refused_not_divided():
    w = _grid()
    model = np.ones_like(w)
    rec, cont = lc.fit_model_guided(w, 0.80 * model, 5006.0, w, model, exclude_half_width_A=0.6)
    assert not rec.applied and cont is None
    assert rec.reason.startswith("CONTINUUM_UNPLACEABLE")


def _wide(level_fn, lo=4900.0, hi=5112.0):
    w = np.arange(lo, hi, 0.01)
    return w, level_fn(w), w, np.ones_like(w)


def test_an_instrumental_sag_is_divided_out():
    # RYA-1315: HARPS' shipped normalisation sits at ~0.92 over 4300-4700 A. The model-guided
    # level IS that error (lines are in the synthesis), so it is applied, not refused --
    # RYA-1350: once the neighbourhood shows the same smooth level.
    w = _grid()
    model = np.ones_like(w)
    ww, wf, rw, rf = _wide(lambda x: 0.92 + 0.0001 * (x - 5006.0))
    check = lc.neighbourhood_check(ww, wf, 5006.0, rw, rf, exclude_half_width_A=0.6)
    rec, cont = lc.fit_model_guided(w, 0.92 * model, 5006.0, w, model, exclude_half_width_A=0.6,
                                    corroborate=check)
    assert rec.applied and rec.level_at_centre == pytest.approx(0.92, abs=1e-6)
    assert np.allclose(0.92 * model / cont, 1.0)
    assert "sag corroborated" in rec.method


def test_a_local_dip_is_not_a_sag():
    # RYA-1350: Si I 9891.725 read 0.92-0.97 while 50 A either side sat at 0.98-1.00;
    # A(Si) 7.56 -> 6.11 when it was divided out. A level the neighbourhood does not share
    # keeps the shipped continuum.
    w = _grid()
    model = np.ones_like(w)
    ww, wf, rw, rf = _wide(lambda x: np.where(np.abs(x - 5006.0) < 6.0, 0.924, 0.995))
    check = lc.neighbourhood_check(ww, wf, 5006.0, rw, rf, exclude_half_width_A=0.6)
    rec, cont = lc.fit_model_guided(w, 0.924 * model, 5006.0, w, model,
                                    exclude_half_width_A=0.6, corroborate=check)
    assert not rec.applied and cont is None
    assert rec.reason.startswith("CONTINUUM_LOCAL") and "0.995" in rec.reason


def test_without_a_neighbourhood_the_three_percent_bound_holds():
    w = _grid()
    model = np.ones_like(w)
    rec, cont = lc.fit_model_guided(w, 0.95 * model, 5006.0, w, model, exclude_half_width_A=0.6)
    assert not rec.applied and cont is None and rec.reason.startswith("CONTINUUM_LOCAL")
    rec, cont = lc.fit_model_guided(w, 0.975 * model, 5006.0, w, model, exclude_half_width_A=0.6)
    assert rec.applied                                    # inside 3%: unchanged behaviour


def test_too_few_neighbour_windows_refuse():
    w = _grid()
    model = np.ones_like(w)
    ww, wf, rw, rf = _wide(lambda x: np.full_like(x, 0.92), lo=4999.0, hi=5023.0)
    check = lc.neighbourhood_check(ww, wf, 5006.0, rw, rf, exclude_half_width_A=0.6)
    rec, _ = lc.fit_model_guided(w, 0.92 * model, 5006.0, w, model, exclude_half_width_A=0.6,
                                 corroborate=check)
    assert not rec.applied and "need 3" in rec.reason


def test_the_envelope_keeps_its_three_percent_bound():
    assert lc.MAX_SHIFT == 0.03 and lc.MAX_SHIFT_MODEL == 0.15
    assert (lc.NEIGH_A, lc.NEIGH_STEP_A, lc.SMOOTH_TOL, lc.MIN_NEIGH) == (50.0, 10.0, 0.02, 3)


def test_record_only_band_measures_but_never_divides():
    w = _grid()
    model = np.ones_like(w)
    rec, cont = lc.fit_model_guided(w, 0.99 * model, 5006.0, w, model,
                                    exclude_half_width_A=0.6, apply=False)
    assert cont is None and not rec.applied
    assert rec.level_at_centre == pytest.approx(0.99, abs=1e-6)


def test_windows_helper_divides_only_its_own_window():
    w = _grid()
    obs = np.full_like(w, 0.99)

    def model_fn(lo, hi):
        mw = np.arange(lo, hi, 0.01)
        return mw, np.ones_like(mw)

    out, recs = lc.apply_to_windows_model_guided(w, obs, [(5005.8, 5006.2)], model_fn)
    assert recs[0]["applied"]
    inside = (w >= 5004.8) & (w <= 5007.2)
    assert np.allclose(out[inside], 1.0) and np.allclose(out[~inside], 0.99)

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
    rec, cont = lc.fit_model_guided(w, 0.92 * model, 5006.0, w, model, exclude_half_width_A=0.6)
    assert not rec.applied and cont is None
    assert rec.reason.startswith("CONTINUUM_UNPLACEABLE")


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

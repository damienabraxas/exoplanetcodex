"""RYA-470 (Ryan 2026-11-18): the IR-HFS special case fits a continuum scale jointly with
A(X), beside -- never instead of -- the standard frozen-continuum fit."""
import inspect
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]

from pipeline import abundances_derive as ad  # noqa: E402
from pipeline.local_continuum import MAX_SHIFT_MODEL  # noqa: E402


def test_the_scale_recovers_a_known_continuum_offset():
    synth = 1.0 - 0.02 * np.exp(-np.linspace(-3, 3, 200) ** 2)
    assert abs(ad._profile_continuum_scale(0.97 * synth, synth) - 0.97) < 1e-12


def test_the_scale_is_bounded_like_the_model_guided_continuum():
    synth = np.ones(50)
    assert ad._profile_continuum_scale(0.5 * synth, synth) == 1.0 - MAX_SHIFT_MODEL
    assert ad._profile_continuum_scale(1.5 * synth, synth) == 1.0 + MAX_SHIFT_MODEL


def test_the_fitter_defaults_to_the_standard_fit():
    p = inspect.signature(ad._fit_synth_flux).parameters["free_continuum"]
    assert p.default is False


def test_the_fit_one_harness_passes_nothing_unless_asked():
    import rya759_nearuv_fe_product as h
    assert inspect.signature(h.fit_one).parameters["free_continuum"].default is None


def test_the_special_case_gets_its_own_stem():
    import derive_band_products as d
    base = dict(element="V", ion="I", lo=14943, hi=17585, instrument="crires_plus",
                holding="solar_crires_plus_h_rya1232", lines_tier="all")
    std = d.synthesis_stem(SimpleNamespace(**base))
    fc = d.synthesis_stem(SimpleNamespace(**base, free_continuum=True))
    assert fc == std + "_FREECONT"


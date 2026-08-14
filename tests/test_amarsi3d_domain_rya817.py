"""
tests/test_amarsi3d_domain_rya817.py — RYA-817
==============================================
The Amarsi 2022 3D-NLTE MLP is reactivated as its own data product. This file guards
the reason it is SAFE to reactivate: a line-parameter domain check that the vendored
code does not have.

Two properties matter and both are asserted on the live artifacts rather than on
fixtures, because a fixture cannot catch the artifact drifting:

  * the recovered training line list really is the training set — proved against the
    moments the vendored StandardScalers carry, which are moments OF THE TRAINING DATA;
  * the domain check DISCRIMINATES — it admits a line the network was trained on and
    refuses a near-IR line, and it refuses that line on the transition-energy axis
    specifically. A check that admitted everything would be worse than none, since it
    would launder an extrapolation as a verified number.

sklearn is not installed in the Sirius CI venv, so anything that touches the pickles is
skipped there rather than failing. The artifact-level checks need no sklearn and always
run — which is deliberate: those are the ones that catch a bad commit.
"""
from __future__ import annotations

import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import amarsi3d  # noqa: E402
from pipeline.band_products import TREATMENTS  # noqa: E402

VENDOR = ROOT / "vendor" / "1L-3NErrors"
TRAINING = ROOT / "data" / "reference" / "amarsi2022_training" / "amarsi2022_training_lines.csv"
SOLAR_CONTROL = (ROOT / "data" / "reference" / "amarsi2022_training"
                 / "amarsi2022_solar_control_lines.csv")

SOLAR = dict(teff=5772.0, logg=4.438, vmic=1.0)

#: A real measured near-IR Fe I line from the RYA-783 band product (8047.618 A).
#: Its E_low and E_up are each individually unremarkable; the PAIR is the problem.
IR_LINE = dict(ion="I", elo_eV=4.988, eup_eV=6.529, loggf=-1.15)
#: A real training line: Fe I 4787.83, the first row of the Jofre golden list and the
#: first row of the vendored test_data.csv.
TRAINING_LINE = dict(ion="I", elo_eV=2.998045, eup_eV=5.586893, loggf=-2.563)


def _sklearn_or_skip():
    pytest.importorskip("sklearn", reason="sklearn absent (Sirius CI venv); "
                                          "MLP-touching checks skipped")


# ── the artifacts ─────────────────────────────────────────────────────────────

def test_training_artifact_matches_the_published_counts_and_range():
    """171 Fe I + 12 Fe II over 4787.83-6810.26 A — the numbers Amarsi 2022 states."""
    df = pd.read_csv(TRAINING)
    assert int((df.species == "Fe1").sum()) == 171
    assert int((df.species == "Fe2").sum()) == 12
    assert df.wavelength_air_A.min() == pytest.approx(4787.83, abs=0.005)
    assert df.wavelength_air_A.max() == pytest.approx(6810.26, abs=0.005)


def test_training_artifact_reproduces_the_vendored_scaler_moments():
    """The scalers carry mean_/scale_ OF THE TRAINING DATA. A wrong list fails this."""
    _sklearn_or_skip()
    df = pd.read_csv(TRAINING)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scalers = {k: pickle.load(open(VENDOR / f, "rb"))[0] for k, f in
                   (("lt02", "fe1_model_lt02.p"), ("gt02", "fe1_model_gt02.p"),
                    ("fe2", "fe2_model.p"))}
    # gt02 (n=154) and fe2 (n=12) land within 0.01; lt02 has 17 lines and the paper's
    # per-model line-overlap cut bites hardest there, so it gets more room.
    tol = {"lt02": 0.13, "gt02": 0.02, "fe2": 0.02}
    for net, sc in scalers.items():
        sub = df[df.network == net]
        assert float(sub.elo_eV.mean()) == pytest.approx(float(sc.mean_[4]),
                                                         abs=tol[net])
        assert float(sub.loggf.mean()) == pytest.approx(float(sc.mean_[6]),
                                                        abs=tol[net])


def test_scaler_check_can_fail_on_a_wrong_line_list():
    """RYA-805: the control must be shown to DISCRIMINATE, not merely to pass."""
    _sklearn_or_skip()
    df = pd.read_csv(TRAINING)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc = pickle.load(open(VENDOR / "fe1_model_gt02.p", "rb"))[0]
    wrong = float(df[df.network == "gt02"].elo_eV.mean()) + 1.0
    assert wrong != pytest.approx(float(sc.mean_[4]), abs=0.02)


def test_solar_control_line_list_is_a_different_list_from_the_training_set():
    """The reactivation control uses Amarsi's SOLAR list, not the training list.

    They are different lists and conflating them costs 0.04 dex on Fe I — the RYA-785
    wrong-referee failure. Asserting they differ keeps that distinction alive.
    """
    tr = pd.read_csv(TRAINING)
    sol = pd.read_csv(SOLAR_CONTROL)
    assert len(sol) != len(tr)
    # the solar list reaches far lower in excitation potential
    assert float(sol[sol.ion == "I"].elo_eV.min()) < 0.5
    assert float(tr[tr.species == "Fe1"].elo_eV.min()) == pytest.approx(0.0, abs=1e-9)
    # and it carries published 1D-LTE abundances, which the training list does not
    assert "a_1d_lte_ap2002" in sol.columns
    assert "a_1d_lte_ap2002" not in tr.columns


# ── the domain check ──────────────────────────────────────────────────────────

def test_treatment_is_registered_in_the_closed_vocabulary():
    """RYA-798 shipped a product whose treatment was not in TREATMENTS and it died at
    build_product AFTER the synthesis had run. Catch that here instead."""
    assert amarsi3d.TREATMENT in TREATMENTS
    assert amarsi3d.TREATMENT != "ENGINE-A"      # never a relabelled ENGINE-A


def test_a_training_line_is_in_domain():
    """The positive half. A check that rejects everything proves nothing."""
    v = amarsi3d.classify_line(**TRAINING_LINE, afe=7.46, **SOLAR)
    assert v.in_domain, v.reason
    assert v.network == "gt02"


def test_a_near_ir_line_is_refused_on_the_transition_energy_axis():
    """The finding, as a test.

    The IR line's E_low and E_up are each inside the trained ranges. It is the PAIR that
    is out of domain, because Eup - Elo IS the wavelength. If this ever starts passing,
    either the training artifact changed or the delta_E axis was dropped — both of which
    would silently turn a refusal into an extrapolated number.
    """
    v = amarsi3d.classify_line(**IR_LINE, afe=7.46, **SOLAR)
    assert not v.in_domain
    assert not v.delta_E_ok
    assert "transition energy" in v.reason

    d = amarsi3d.domains()["gt02"]
    assert d.elo[0] <= IR_LINE["elo_eV"] <= d.elo[1]     # E_low alone is fine
    assert d.eup[0] <= IR_LINE["eup_eV"] <= d.eup[1]     # E_up alone is fine
    assert v.delta_E_eV < d.delta_E[0]                    # the pair is not


def test_no_training_line_reaches_the_measured_ir_band():
    """The domain gap is total, not marginal: the training set's minimum transition
    energy sits above the maximum a 6910-9199 A line can have."""
    df = pd.read_csv(TRAINING)
    hc_ev_A = 12398.419843320026
    # the LARGEST transition energy available anywhere in the IR band (its blue edge)
    s2 = (1e4 / 6910.0) ** 2
    n = 1.0 + (8342.13 + 2406030.0 / (130.0 - s2) + 15997.0 / (38.9 - s2)) * 1e-8
    ir_max_dE = hc_ev_A / (6910.0 * n)
    gap = float(df.delta_E_eV.min()) - ir_max_dE
    assert gap > 0
    # and the edge tolerance cannot close it — the gap is more than an order of
    # magnitude larger, so admitting the IR would take a deliberate, visible change.
    assert gap > 10 * amarsi3d.BOX_TOL


def test_the_domain_gap_WIDENS_with_wavelength_rya762():
    """RYA-762 asked whether the 3D leg can extend past Engine B's 9199.9 A wall.

    It cannot, and the reason is worth pinning: transition energy falls as 1/lambda, so
    every extra Angstrom of reach moves a band FURTHER from the training floor. The
    9199-13000 A extension misses by ~17x what the 6910-9199 A band misses by. If this
    ever inverts, someone has changed the training artifact.
    """
    df = pd.read_csv(TRAINING)
    hc_ev_A = 12398.419843320026

    def max_dE(lambda_blue_edge_A: float) -> float:
        s2 = (1e4 / lambda_blue_edge_A) ** 2
        n = 1.0 + (8342.13 + 2406030.0 / (130.0 - s2) + 15997.0 / (38.9 - s2)) * 1e-8
        return hc_ev_A / (lambda_blue_edge_A * n)

    floor = float(df.delta_E_eV.min())
    gap_ir = floor - max_dE(6910.0)          # RYA-783 band
    gap_far = floor - max_dE(9199.0)         # RYA-762 extension
    assert gap_ir > 0 and gap_far > 0
    assert gap_far > gap_ir                  # further out, not closer
    assert gap_far / gap_ir > 10


def test_out_of_domain_returns_nan_unless_explicitly_allowed():
    """No silent extrapolation, and the escape hatch has to be asked for by name."""
    _sklearn_or_skip()
    refused, v = amarsi3d.aberr_for_line(**IR_LINE, a_1dlte=7.6, **SOLAR)
    assert np.isnan(refused)
    assert not v.in_domain

    forced, _ = amarsi3d.aberr_for_line(**IR_LINE, a_1dlte=7.6, **SOLAR,
                                        allow_out_of_domain=True)
    assert np.isfinite(forced)          # the number exists; it just may not be a product


def test_stellar_box_is_never_stepped_over_even_when_extrapolation_is_allowed():
    """`allow_out_of_domain` relaxes the LINE axes. The authors' own stellar box is a
    hard guard and stays hard — Procyon (6554 K) must be refused."""
    _sklearn_or_skip()
    ab, v = amarsi3d.aberr_for_line(**TRAINING_LINE, a_1dlte=7.6,
                                    teff=6554.0, logg=4.0, vmic=1.8,
                                    allow_out_of_domain=True)
    assert np.isnan(ab)
    assert not v.stellar_ok


# ── the archived leg is unchanged ─────────────────────────────────────────────

def test_pinned_axis_is_opt_in_and_default_behaviour_is_untouched():
    """RYA-817 added `afe3n_axis` to nlte_corrections._apply_aberr_to_line. Every
    existing caller passes nothing, so the default must reproduce the old per-line
    iteration exactly."""
    _sklearn_or_skip()
    from pipeline.nlte_corrections import _apply_aberr_to_line
    args = ("I", TRAINING_LINE["elo_eV"], TRAINING_LINE["eup_eV"],
            TRAINING_LINE["loggf"], 7.60, SOLAR["teff"], SOLAR["logg"], SOLAR["vmic"])
    default = _apply_aberr_to_line(*args)
    explicit_none = _apply_aberr_to_line(*args, afe3n_axis=None)
    assert default == explicit_none
    # and pinning the axis to a different value must actually change something,
    # otherwise the parameter is decorative
    pinned = _apply_aberr_to_line(*args, afe3n_axis=5.0)
    assert pinned != default

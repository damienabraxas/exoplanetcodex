"""
tests/test_acen_orbit_rya423.py
===============================
RYA-423 — the alpha Cen AB orbit ephemeris (pipeline.acen_orbit). Sanity-anchors the derived
RV amplitudes and the predicted A/B RVs against the Kervella 2016 elements + the observed
NIRPS RV of the spectroscopically-confirmed A frames (~-26.7 km/s, 2024-2025).
"""
import numpy as np
from astropy.time import Time
from pipeline import acen_orbit as O


def test_derived_K_amplitudes_reasonable():
    # K_A ~ 4.6, K_B ~ 5.5 km/s from (a, i, P, e, masses); B heavier-pulled? no, B lighter -> larger K_B
    assert 4.0 < O.K_A < 5.2
    assert 5.0 < O.K_B < 6.0
    assert O.K_B > O.K_A                       # lighter star (B) has the larger reflex amplitude


def test_systemic_and_bounds():
    # every predicted RV must lie within gamma +/- max(K) (a bound orbit cannot exceed this)
    for yr in np.linspace(2022, 2026, 40):
        p = O.predicted_rv(Time(yr, format='byear').mjd)
        assert abs(p['rv_A'] - O.GAMMA) <= O.K_A + 1e-6
        assert abs(p['rv_B'] - O.GAMMA) <= O.K_B + 1e-6


def test_A_branch_matches_confirmed_nirps():
    # the omega assignment is pinned so alpha Cen A is the more-negative branch (~-26 in 2024-25),
    # matching the J-depth+flux-confirmed-A NIRPS frames (obs ~-26.7). This is the data anchor.
    p = O.predicted_rv(Time(2024.7, format='byear').mjd)
    assert -27.5 < p['rv_A'] < -25.0          # alpha Cen A near -26.7
    assert p['rv_A'] < p['rv_B']              # A more negative than B now
    assert 7.0 < p['sep_kms'] < 9.5           # ~8 km/s split, matching the observed NIRPS groups


def test_separation_grows_2022_to_2025():
    s22 = O.predicted_rv(Time(2022.0, format='byear').mjd)['sep_kms']
    s25 = O.predicted_rv(Time(2025.0, format='byear').mjd)['sep_kms']
    assert s25 > s22                           # opening since the ~2016 projected minimum

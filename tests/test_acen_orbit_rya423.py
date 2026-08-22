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


def test_A_branch_matches_the_photometric_A_holding():
    """🔴 THIS TEST REPLACES `test_A_branch_matches_confirmed_nirps`, WHICH IS RETIRED.

    THE RETIRED ANCHOR, and why it went. That test pinned alpha Cen A to the more-negative
    branch (~-26.7 in 2024-25) on RYA-384's "J-depth+flux-confirmed-A" NIRPS frames. It was
    not a stale pin — it correctly FAILED when the omega convention was first flipped, and
    that failure was the reason the flip was reverted for a day. It went because the
    DISCRIMINATOR under it was measured and found to have no power:

        NIRPS jdepth, n=57.  Rule: jd < 0.45 => A, jd > 0.55 => B.
          AlphaCenB  n=28  jd 0.131  rv -26.67   reads A
          Star S5    n= 9  jd 0.256  rv -20.81   reads A
          alf Cen A  n=20  jd 0.822  rv -34.55   reads B  <- and OFF-ORBIT (RYA-431)

    37 read A, 20 read B, and every one of the 20 is the off-orbit contaminant. EVERY
    genuine alpha Cen frame lands 0.125-0.333, nowhere near the 0.45/0.55 boundary. The
    statistic separates "alpha Cen-like" from "hot standard"; it does not resolve a
    G2V/K1V split of 560 K. Confirming "A" against a population that contains no B is not
    a confirmation of branch identity, so the anchor had no discriminating power and is
    retired rather than contradicted.

    WHAT REPLACES IT — a discriminator that IS measured to work. alpha Cen A is 1.3 mag
    brighter than B, and the matched CRIRES K2192 pair (same night, 16 min apart, same
    reduction) gives a count-rate ratio of 2.293 against 2.270 expected for correct labels
    and 0.441 for swapped. Photometry separates the components; J-depth does not. So the
    anchor is now: the holding that is PHOTOMETRICALLY A must sit on the orbit's A branch.

    Ruled by Ryan 2026-08-22 (RYA-971). RYA-431's off-orbit NOT-ALPHA-CEN finding is
    unaffected — `rv_bounds`/`consistent_with_orbit` are symmetric about gamma.
    """
    # The CRIRES A-directory frames: 2022-04-15, measured -19.17 km/s (RYA-963), and the
    # photometrically brighter component. This is the acceptance gate of the RYA-971
    # ruling, kept executable so the convention cannot drift back silently.
    p = O.predicted_rv(Time('2022-04-15T04:00:00', format='isot').mjd)
    measured = -19.17
    assert abs(p['rv_A'] - measured) < abs(p['rv_B'] - measured), (
        "the photometrically-brighter alpha Cen A holding must land on the orbit's A "
        "branch; if it lands on B the omega convention is inverted again")
    assert abs(p['rv_A'] - measured) < 1.0     # resid 0.13 at the ruling
    assert p['rv_A'] > p['rv_B']               # under the visual convention A is the LESS negative branch here


def test_the_retired_nirps_anchor_is_not_silently_reinstated():
    """The retired anchor asserted A near -26.7 in 2024.7. Under the ruled convention that
    is the B branch. Asserting the OPPOSITE keeps the retirement visible: if someone flips
    omega back, this fails and points at the ruling rather than at a mystery."""
    p = O.predicted_rv(Time(2024.7, format='byear').mjd)
    assert -27.5 < p['rv_B'] < -25.0, "B, not A, occupies the -26.7 branch after RYA-971"
    assert p['rv_A'] > p['rv_B']


def test_separation_grows_2022_to_2025():
    s22 = O.predicted_rv(Time(2022.0, format='byear').mjd)['sep_kms']
    s25 = O.predicted_rv(Time(2025.0, format='byear').mjd)['sep_kms']
    assert s25 > s22                           # opening since the ~2016 projected minimum


def test_orbit_bounds_gate_rya431():
    # RYA-431: any bound alpha Cen member is confined to gamma +/- max(K). The 20 "NIRPS B"
    # frames sit at -34.6 km/s -- OFF the orbit (a different K star, not alpha Cen B).
    lo, hi = O.rv_bounds()
    assert lo < O.GAMMA - O.K_B and hi > O.GAMMA + O.K_B   # brackets the systemic +/- reflex
    assert O.consistent_with_orbit(-26.7)      # confirmed alpha Cen A
    assert O.consistent_with_orbit(-18.0)      # alpha Cen B near apoapsis branch
    assert not O.consistent_with_orbit(-34.6)  # the off-orbit "NIRPS B" frames
    assert not O.consistent_with_orbit(float('nan'))
    assert not O.consistent_with_orbit(None)

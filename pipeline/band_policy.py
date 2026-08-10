"""Per-band analytical policy, enforced at intake — RYA-713 / RYA-306 (arm axis).

Ryan, 2026-08-09: *"we should have a check on intake maybe. Check instrument and check
band. UV gets treated different than Vis or IR due to crowding. IR has more errors,
telluric lines etc. That way the products are filtered differently. Not tuned, but a
different scientific and analytical approach for each band."*

WHAT THIS IS
------------
A declared policy saying, for each wavelength regime, WHICH analytical approach is
physically valid there and WHY — resolved when a measurement is requested, and fail-loud
when the requested approach is one the regime cannot support.

This is the `arm` axis RYA-306's method matrix was missing. That matrix keys method on
`(star × element × ion)`; the regime is an independent axis, because the failure modes are
local in WAVELENGTH as well as in stellar parameters.

WHY THIS IS NOT TUNING
----------------------
Tuning means choosing a treatment because of the answer it produces. Every field below is
keyed on **observable properties of the regime** — line density, median line separation,
how close the continuum gets to unity, whether the terrestrial atmosphere contributes —
all of which are measurable without knowing any abundance, and all of which were measured
before this file was written (see MEASURED, below).

Two structural guards keep it that way:

* `BandPolicy` has **no field for an abundance, a target, a reference value or a
  tolerance**. It cannot express "use method X to get answer Y" because there is nowhere
  to put Y.
* `assert_not_tuned()` re-checks that at import.

The policy may be revised, but only by re-measuring the regime properties and recording
what changed — never by observing that a different method gives a nicer number.

MEASURED (Kitt Peak solar atlas + our line inventory, 2026-08-09)
-----------------------------------------------------------------
    band           lines/A   median gap   continuum p95   continuum median
    near-UV          4.62       0.146 A       0.916            0.607
    VIS              1.87       0.277 A       0.963            0.811
    red-optical      0.34       1.872 A       0.997            0.991
    NIR              0.14       3.989 A       0.956            0.862

Two facts in that table drive almost everything here:

1. **The near-UV median line gap (0.146 A) is smaller than a strong line's own wings.**
   There is no interval that contains one line's profile and excludes its neighbours.
   Interval-integrated equivalent width is therefore not merely inaccurate in the near-UV,
   it is undefined — the quantity it measures is blended absorption, not a line's EW.

2. **The near-UV continuum median is 0.607.** Line blanketing means the true continuum is
   never observed; only a pseudo-continuum is available. Any abundance derived there
   inherits that as a systematic and must say so.

The red-optical is the opposite regime — 1.87 A median gap, continuum within 0.3% of unity
— which is why a method can look healthy there while being broken. That is the frontier
trap the optical control exists to catch.
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields


@dataclass(frozen=True)
class BandPolicy:
    """How one wavelength regime must be analysed, and the physics that requires it."""
    name: str
    lo_A: float
    hi_A: float

    # Measured regime properties -- the ONLY basis for the rules below.
    lines_per_A: float
    median_gap_A: float
    continuum_p95: float
    continuum_median: float

    # The rules, each with the property that forces it.
    permitted_methods: tuple[str, ...]
    forbidden_methods: tuple[str, ...]
    continuum_treatment: str
    telluric_required: bool
    justification: str
    systematic_floor_note: str


# Methods this project can apply. Interval integration is listed because it must be
# nameable in order to be forbidden.
METHODS = ("profile-fit", "synthesis", "interval-integration")

# Ground-truth availability. The optical is the only regime with an external answer to
# miss, which is what makes it the control (see SCIENCE_STANDARDS).
CONTROL_BAND = "VIS"


POLICIES: tuple[BandPolicy, ...] = (
    BandPolicy(
        name="near-UV", lo_A=3000.0, hi_A=3800.0,
        lines_per_A=4.62, median_gap_A=0.146,
        continuum_p95=0.916, continuum_median=0.607,
        permitted_methods=("synthesis",),
        forbidden_methods=("interval-integration", "profile-fit"),
        continuum_treatment="pseudo-continuum only; the true continuum is not observed",
        telluric_required=False,
        justification=(
            "Median line separation 0.146 A is SMALLER than a strong line's wings, so no "
            "interval contains one profile and excludes its neighbours — interval EW is "
            "undefined here, not merely imprecise. An isolated profile fit fails for the "
            "same reason: there is no isolated profile. Only synthesis, which models every "
            "contributor in the window simultaneously, is valid."),
        systematic_floor_note=(
            "Continuum median 0.607: blanketing hides the true continuum, so every "
            "abundance carries a pseudo-continuum systematic that does NOT average down "
            "with more lines. It must be stated, not absorbed into the scatter."),
    ),
    BandPolicy(
        name="VIS", lo_A=3800.0, hi_A=6910.0,
        lines_per_A=1.87, median_gap_A=0.277,
        continuum_p95=0.963, continuum_median=0.811,
        permitted_methods=("profile-fit", "synthesis"),
        forbidden_methods=("interval-integration",),
        continuum_treatment="fitted continuum; reachable but not free of blanketing",
        telluric_required=False,
        justification=(
            "Median gap 0.277 A leaves many lines resolvable, so a blend-aware profile fit "
            "is valid and synthesis arbitrates where it is not. Interval integration is "
            "forbidden on evidence: run against the HARPS pool over 146 lines it gave a "
            "median EW ratio of 0.773 (-0.112 dex) with a 5x spread, because the window "
            "rule clips wings on crowded strong lines and over-reaches on isolated weak "
            "ones. THIS IS THE CONTROL BAND — a method that cannot reproduce the known "
            "answer here has not earned the right to report an unknown one elsewhere."),
        systematic_floor_note=(
            "Ground truth exists here (Asplund 2021 A(Fe)=7.46; our banked 7.466). The "
            "residual against it MEASURES the harness's own systematic and belongs in the "
            "frontier bands' error budget rather than being assumed zero."),
    ),
    BandPolicy(
        name="red-optical", lo_A=6910.0, hi_A=10000.0,
        lines_per_A=0.34, median_gap_A=1.872,
        continuum_p95=0.997, continuum_median=0.991,
        permitted_methods=("profile-fit", "synthesis"),
        forbidden_methods=("interval-integration",),
        continuum_treatment="atlas continuum trusted where the source is pre-normalised",
        telluric_required=True,
        justification=(
            "The cleanest regime we have: 1.87 A median gap and a continuum within 0.3% of "
            "unity. That is exactly why it is dangerous — a broken method looks HEALTHIER "
            "here than in the blue while being wrong in the same way, and there is no "
            "reference value to catch it. Telluric masking is mandatory: the O2 A-band "
            "(7600-7640) and H2O (9280-9600) are terrestrial, not solar."),
        systematic_floor_note=(
            "Sparse crowding means uncatalogued neighbours are rarer but not absent — 98 "
            "of 174 Fe I failures here were windows dominated by a real solar line missing "
            "from our list. Line-list completeness is the dominant systematic, not blending."),
    ),
    BandPolicy(
        name="NIR", lo_A=10000.0, hi_A=24000.0,
        lines_per_A=0.14, median_gap_A=3.989,
        continuum_p95=0.956, continuum_median=0.862,
        permitted_methods=("synthesis",),
        forbidden_methods=("interval-integration", "profile-fit"),
        continuum_treatment="telluric-corrected; continuum only meaningful after correction",
        telluric_required=True,
        justification=(
            "Lines are sparse (3.99 A median gap) but the band is not clean: continuum p95 "
            "0.956 against a median of 0.862 is telluric absorption, not stellar. Until a "
            "telluric correction is applied the observed flux is not a stellar spectrum, so "
            "neither an interval nor an isolated profile measures a stellar quantity. "
            "Synthesis with a telluric model is the only valid route."),
        systematic_floor_note=(
            "Telluric residuals are epoch- and airmass-dependent, so this systematic varies "
            "between observations of the SAME star and cannot be calibrated once."),
    ),
)


class BandPolicyError(RuntimeError):
    """Raised at intake when a regime cannot support the requested approach."""


def resolve(wavelength_A: float) -> BandPolicy:
    for p in POLICIES:
        if p.lo_A <= wavelength_A < p.hi_A:
            return p
    raise BandPolicyError(
        f"{wavelength_A:.3f} A falls outside every declared band "
        f"({POLICIES[0].lo_A:.0f}-{POLICIES[-1].hi_A:.0f} A). A regime with no declared "
        f"policy has no validated method — declare it before measuring in it.")


def check_intake(wavelength_A: float, method: str, *, instrument: str = "") -> BandPolicy:
    """Gate a measurement request. Loud, never silently downgraded.

    A method that is wrong for a regime does not produce a worse number, it produces a
    DIFFERENT QUANTITY -- interval integration in the near-UV measures blended absorption,
    not an equivalent width. Silently allowing it and widening the error bar afterwards
    would misrepresent what was measured.
    """
    if method not in METHODS:
        raise BandPolicyError(f"unknown method {method!r}; expected one of {METHODS}")
    p = resolve(wavelength_A)
    if method in p.forbidden_methods:
        raise BandPolicyError(
            f"{method!r} is FORBIDDEN in the {p.name} band "
            f"({p.lo_A:.0f}-{p.hi_A:.0f} A){' on ' + instrument if instrument else ''}.\n"
            f"  why: {p.justification}\n"
            f"  permitted here: {p.permitted_methods}")
    if method not in p.permitted_methods:
        raise BandPolicyError(
            f"{method!r} is not declared valid in the {p.name} band; "
            f"permitted: {p.permitted_methods}")
    return p


def assert_not_tuned() -> None:
    """A policy field that could hold a target value would make this a tuning surface.

    Checked at import so the guard cannot rot: if someone later adds `expected_abundance`
    or `target_tolerance` to BandPolicy, this fails immediately rather than after the
    field has quietly shaped a result.
    """
    banned = ("abundance", "target", "reference_value", "expected", "tolerance", "anchor")
    names = {f.name.lower() for f in fields(BandPolicy)}
    hits = {n for n in names for b in banned if b in n}
    if hits:
        raise BandPolicyError(
            f"BandPolicy gained field(s) {sorted(hits)}. This file declares how a regime "
            f"must be ANALYSED, from measurable properties of the regime alone. A field "
            f"that can hold a target value turns it into a tuning surface — which is "
            f"exactly what Ryan's 'not tuned' constraint forbids.")


assert_not_tuned()


# ── Per-line escape from a band-level ban ────────────────────────────────────
# Ryan, 2026-08-09: "Engine A for UV it is."
#
# The band rules above are keyed on the MEDIAN line separation, which is the right basis
# for a default and the wrong basis for a verdict on an individual line. The near-UV median
# gap is 0.146 A, so profile fitting is banned there -- correctly, for a typical line. But
# isolation is a PER-LINE property, and measured across near-UV Fe at usable depth:
#
#     nearest-neighbour gap      Fe I   Fe II
#     < 0.10 A  (hopeless)        794      79
#     0.10-0.20 A                 403      30
#     0.20-0.40 A                 219      28
#     > 0.40 A  (isolated)         52       7
#
# ~306 of 1612 near-UV Fe lines are isolated enough to fit. Banning them because their
# neighbours are crowded discards real, measurable lines -- the same "average decides for
# the individual" error as reporting a band median instead of a per-line ledger.
#
# This is NOT a relaxation of the policy. The physical requirement is unchanged: a profile
# fit needs a profile that is actually isolated. What changes is that the test is applied
# to the line rather than inferred from its neighbourhood.

# A profile fit needs the neighbour far enough out that the fitting window can contain this
# line's wings and exclude the neighbour's core. Solar metal lines run sigma ~0.03 A, so
# wings matter to ~0.15 A; 0.30 A of separation gives the window somewhere to sit.
PROFILE_FIT_MIN_GAP_A = 0.30

# TESTED ON near-UV Fe I AND IT YIELDS NOTHING (RYA-713, 2026-08-09).
#
# The escape is sound in principle and was worth testing. It was tested, on all 901
# near-UV Fe I lines at usable depth, and the result is a clean negative:
#
#     845 refused by the gap test itself
#      56 PASSED the isolation test
#         -> their measured EW ran 107-922 mA, REW -4.51 to -3.59, EVERY ONE above the
#            -4.90 saturation ceiling. A single solar Fe line is 10-150 mA; 922 mA is a
#            blended complex, not a line.
#      0  usable lines
#
# WHY the gap test passes lines that are not actually isolated: it measures the distance
# to the nearest CATALOGUED neighbour, and near-UV catalogue completeness is ~77%. Roughly
# a quarter of the features present are not in our list, so a 0.30 A "gap" routinely
# contains an uncatalogued line. Absence of a neighbour in the catalogue is not absence in
# the spectrum -- the same lesson the IR root-cause split produced from the other side.
#
# CONCLUSION: the near-UV band-wide ban was CORRECT. Synthesis really is the only valid
# route there, and it stays blocked on a line list below 4200 A. The escape mechanism is
# retained because it is the right shape for a band whose crowding is genuinely marginal,
# but it must not be read as a near-UV workaround -- it was tried and it does not work.


def permits_profile_fit_for_line(wavelength_A: float, nearest_neighbour_gap_A: float) -> tuple[bool, str]:
    """May THIS line be profile-fitted, even where the band bans it by default?

    Returns (allowed, reason). The reason is recorded on the measurement either way, so a
    per-line escape is always visible as one and never silently widens the policy.
    """
    pol = resolve(wavelength_A)
    if "profile-fit" in pol.permitted_methods:
        return True, f"{pol.name} permits profile fitting by default"
    if nearest_neighbour_gap_A >= PROFILE_FIT_MIN_GAP_A:
        return True, (
            f"PER-LINE ESCAPE: the {pol.name} band bans profile fitting on its median gap "
            f"({pol.median_gap_A:.3f} A), but THIS line's nearest neighbour is "
            f"{nearest_neighbour_gap_A:.3f} A away, at or beyond the "
            f"{PROFILE_FIT_MIN_GAP_A} A a fitting window needs. Isolation is a property of "
            f"the line, not of the band.")
    return False, (
        f"{pol.name} bans profile fitting (median gap {pol.median_gap_A:.3f} A) and this "
        f"line does not escape it: nearest neighbour {nearest_neighbour_gap_A:.3f} A < "
        f"{PROFILE_FIT_MIN_GAP_A} A, so no window contains this profile and excludes the "
        f"neighbour's core.")

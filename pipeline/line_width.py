"""How wide a line is, and how narrow physics lets it be — RYA-906/911.

🔴 WHY THIS IS A MODULE AND NOT A HELPER IN THE HARNESS THAT NEEDED IT FIRST.

RYA-906/911 (PR #315) established that a Voigt fit may not be judged by its Gaussian
`sigma`: sigma and gamma are degenerate, so a railed sigma records where the optimiser
resolved the degeneracy, not how wide the line is. The corrected test — total Voigt FWHM
against an instrumental (+) thermal (+) microturbulent floor — was written into
`scripts/measure_band_profilefit.py`.

**It landed on ONE of the two profile-fit implementations.** `pipeline/measure/
profile_fit.py` — the registered `ProfileFitHandler` that `resolve_handler()` returns for
VIS and red-optical — kept the refuted `abs(sigma - floor) < 1e-4` test for another day,
because the corrected physics sat in a *script* whose import chain loads the Kitt Peak
atlas and therefore could not be imported from the pipeline at all.

So the physics lives here, depending on nothing but `config`, and both fitters import it.
A guard that two routes must agree on is not a local helper.
"""
from __future__ import annotations

import math

from config.constants import STAR_PARAMS

C_KMS = 299792.458

#: RYA-906/911 — THE IRREDUCIBLE WIDTH OF A SOLAR IRON LINE.
#:
#: Thermal Doppler + microturbulence, from the star's own ratified parameters. This is a
#: LOWER BOUND, not an expectation: MACROturbulence (vmac 3.8 km/s) is deliberately left
#: out, so a line narrower than this is narrower than physics permits regardless of how
#: much macroturbulent broadening the star actually has.
#:
#: Doppler core width Delta_v_D = sqrt(2kT/m + xi^2); the Gaussian sigma is Delta_v_D/sqrt(2).
_K_B, _AMU, _M_FE = 1.380649e-23, 1.66054e-27, 55.845 * 1.66054e-27


class MissingBroadeningTerm(KeyError):
    """A star lacks a parameter one of these widths needs. Never substitute another star's."""


def star_parameter_provenance(star: str) -> dict:
    """Which of this star's parameters are PINNED fundamentals and which are SOLVED.

    🔴 RYA-985 / RYA-957: a harness that presents a SOLVED parameter as a fundamental is
    doing what the source paper forbids. Heiter+2015 Table 10 prints tau Ceti's log g in
    square brackets — "uncertain and should not be used as a reference for calibration or
    validation purposes" — so `stars.yaml` gives it `pin: [teff, xi]`, `solve: [feh, logg]`.
    The value 4.49 is a STARTING POINT, not an adopted fundamental, and any product built on
    it must say so rather than quoting it like alpha Cen A's dynamical log g.

    Returned so a caller can put it in provenance; this module does not decide what to do
    with it, it only refuses to let the distinction go unrecorded.
    """
    from config.constants import get_star_params
    p = get_star_params(star)
    pin = [str(x) for x in (p.get("pin") or [])]
    solve = [str(x) for x in (p.get("solve") or [])]
    return {
        "star": star, "pinned": pin, "solved": solve,
        "logg_is_fundamental": "logg" in pin,
        "note": (f"{star}: pinned {pin or 'none'}; solved {solve or 'none'}. "
                 + ("log g is a PINNED fundamental." if "logg" in pin else
                    "log g is SOLVED — the stars.yaml value is a starting point, NOT an "
                    "adopted fundamental, and must not be quoted as one.")),
        "source": str(p.get("source", ""))[:200],
        "logg_basis": str(p.get("logg_basis", ""))[:200],
    }


def _star(star: str) -> dict:
    """Resolve through `get_star_params` — the single source (RYA-298/355), and it LOUD-FAILS
    on an unknown id rather than falling back to solar."""
    from config.constants import get_star_params
    return get_star_params(star)


def irreducible_sigma_kms(star: str = "solar") -> float:
    """Thermal+microturbulent sigma for Fe, in km/s, from config (never a literal here).

    RYA-985: this now follows the star, which is the debt the previous docstring recorded —
    "a floor derived from the wrong star's temperature is exactly the kind of silent stand-in
    this codebase keeps finding." `star` defaults to solar so every existing caller is
    bit-identical.
    """
    p = _star(star)
    if "xi" not in p:
        raise MissingBroadeningTerm(
            f"{star!r} has no microturbulence `xi` in config/stars.yaml, so the width FLOOR "
            f"cannot be computed. 55 Cnc A carries `xi_init`/`xi_xcheck` instead because it "
            f"SOLVES xi (RYA-957) — a solved parameter is not a config constant, and "
            f"substituting the solar 1.0 would put a different star's turbulence in the floor.")
    return math.sqrt(_K_B * float(p["teff"]) / _M_FE
                     + (float(p["xi"]) * 1000.0) ** 2 / 2.0) / 1000.0


def voigt_fwhm(sigma_A: float, gamma_A: float | None) -> float:
    """Total FWHM of a Voigt profile. Olivero & Longbothum (1977), ~0.02% accurate.

    🔴 THIS IS THE QUANTITY THAT MEANS SOMETHING, AND `sigma` ALONE IS NOT.
    In a Voigt fit the Gaussian sigma and the Lorentzian gamma are DEGENERATE: both
    broaden the line and the optimiser can trade one against the other at fixed total
    width. Judging a Voigt fit by its sigma alone therefore measures where the
    degeneracy happened to land, not how wide the line is.
    """
    f_G = 2.0 * sigma_A * math.sqrt(2.0 * math.log(2.0))
    if not gamma_A:
        return f_G
    f_L = 2.0 * float(gamma_A)
    return 0.5346 * f_L + math.sqrt(0.2166 * f_L * f_L + f_G * f_G)


def physical_floor_fwhm(wavelength_A: float, sigma_inst_A: float, star: str = "solar") -> float:
    """The narrowest TOTAL width physics allows here: instrumental (+) thermal (+) micro."""
    sig_phys = wavelength_A * irreducible_sigma_kms() / C_KMS
    return voigt_fwhm(math.hypot(sigma_inst_A, sig_phys), None)


def gamma_of(popt, ptype: str) -> float | None:
    """The Lorentzian half-width of a fit, or None when the fit has no Lorentzian.

    `popt` is `[x0, depth, sigma, (gamma)]` (`pipeline.lines_fit._fit_profile`). A
    Gaussian fallback has no fourth element, and `None` is the honest answer for it —
    not 0.0, which would read as a measured zero.
    """
    return float(popt[3]) if ptype == "voigt" and len(popt) > 3 else None


def total_width_below_physical_floor(popt, ptype: str, wavelength_A: float,
                                     sigma_inst_A: float) -> bool:
    return voigt_fwhm(float(popt[2]), gamma_of(popt, ptype)) <= physical_floor_fwhm(
        wavelength_A, sigma_inst_A)


def under_physical_width_reason(popt, ptype: str, wavelength_A: float,
                                sigma_inst_A: float) -> str:
    """The exclusion sentence, so both fitters quarantine under ONE coded reason.

    The reason string is read by downstream root-cause attribution and by anyone auditing
    a pool, so two routes emitting two different sentences for one physical verdict is the
    RYA-911 labelling problem in miniature.
    """
    fwhm = voigt_fwhm(float(popt[2]), gamma_of(popt, ptype))
    floor = physical_floor_fwhm(wavelength_A, sigma_inst_A)
    return (f"UNDER-PHYSICAL-WIDTH: the fitted {ptype} profile has total FWHM "
            f"{fwhm:.4f} A, at or below the {floor:.4f} A that thermal + "
            f"microturbulent broadening alone impose at this wavelength "
            f"(instrumental sigma {sigma_inst_A:.4f} A; macroturbulence excluded, so this "
            f"is a floor and not an expectation). A line cannot be this narrow, so "
            f"the integrated EW is not trustworthy.")


# ── THE CEILING (RYA-959) ────────────────────────────────────────────────────
#
# 🔴 THE FLOOR ABOVE HAD NO MIRROR, AND THAT ASYMMETRY IS A MEASURED DEFECT.
#
# `total_width_below_physical_floor` convicts a fit that is too NARROW. Nothing convicted
# a fit that was too WIDE, so a fit that failed for any reason — a blend, a misplaced
# continuum, a misidentification — escaped UPWARD instead: `measure_band_profilefit.
# instrument_sigma` returns `sigma_max = 0.40 A`, which is ~5x the widest Gaussian solar
# physics permits (0.083 A at 5500 A, HARPS), and the optimiser used all of it.
#
# Measured on the first fresh RYA-959 Fe I VIS run: the fitted sigma is BIMODAL — either
# a physical core (0.015-0.04 A) or pinned exactly at 0.400 with gamma also at its bound.
# A pinned pair integrates to an EW of 400-1000 mA on an Fe I line and returns
# `red_chi2` ~ 0.02, i.e. it looks converged. That is how RYA-958's 251-of-444 pool of
# physically impossible line widths was built.
#
# TWO CHECKS, BECAUSE THERE ARE TWO WAYS TO BE TOO WIDE:
#
#  1. `gaussian_sigma_above_physical_ceiling` — the DOPPLER bound. The Gaussian component
#     of a stellar line is thermal (+) microturbulent (+) macroturbulent (+) rotational
#     (+) instrumental, and nothing else. Pressure broadening is LORENTZIAN and belongs
#     to gamma, so gamma is deliberately NOT bounded here: damping wings are real and a
#     ceiling on them would quarantine the strongest lines in the band, which is exactly
#     the population RYA-945's laboratory gf backbone lives in.
#
#  2. `implied_width_exceeds_ceiling` — the INTEGRAL bound, and the only check in this
#     harness that referees the EW itself rather than the profile's parameters. It asks
#     what FWHM a Gaussian of the line's OBSERVED core depth would need in order to hold
#     the EW that was integrated. A runaway gamma inflates the EW without deepening the
#     core, so it is invisible to check 1 and convicted here.
#
# ⚠️ THE OBSERVED DEPTH, NEVER THE PREDICTED ONE. RYA-958 diagnosed the stale pool with
# `central_depth` out of `linelist_solar.csv`, which is a VALD MODEL quantity. Rejecting a
# MEASUREMENT on a model number would let a wrong predicted depth convict a good fit, and
# `verify_feature` already owns the observed-vs-predicted comparison (GF-GHOST). This test
# is measurement against measurement: the integral against the core it came from.
#
# ⚠️ GENEROUS ON PURPOSE, in the same direction the floor is strict on purpose. `vmac` is
# a RADIAL-TANGENTIAL value (3.8 km/s solar) and is used here as though it were a Gaussian
# sigma, which overstates it by ~2.4x; `vsini` is used the same way. The floor excludes
# macroturbulence so it stays a true lower bound; the ceiling overstates it so it stays a
# true upper bound. Neither is an expectation, and neither may be tightened to make a
# number come out better.

#: EW of a Gaussian = depth * FWHM * sqrt(2*pi) / (2*sqrt(2*ln2)). DERIVED, not typed —
#: this is the factor RYA-958's diagnostic wrote as the bare literal 1.06.
EQUIV_RECT_FACTOR = math.sqrt(2.0 * math.pi) / (2.0 * math.sqrt(2.0 * math.log(2.0)))

#: 🔴 NOT A TUNED THRESHOLD, AND THE FIRST VERSION OF IT WAS.
#:
#: This started as `IMPLIED_WIDTH_ALLOWANCE = 3.0`, a multiple of the Doppler FWHM chosen
#: to look generous. Measured on the RYA-959 Fe I VIS run, the implied/Doppler ratio of
#: the sigma-clean fits is a SMOOTH CONTINUUM -- 50/90/95/99th percentiles 1.22 / 2.90 /
#: 3.34 / 4.04, no valley anywhere -- so any multiple is a cut through a populated
#: distribution, i.e. exactly the tuning this project forbids. It was replaced by a bound
#: that is a fact about the harness instead of a taste about the physics.
#:
#: `pipeline.lines_fit._integrate_profile` integrates the model over the FIT WINDOW ONLY
#: (`x = linspace(wav[0], wav[-1], 5000)`), so the EW it returns cannot exceed
#: `model_depth * window_width`. An equivalent rectangle of the OBSERVED depth that is
#: wider than that window is therefore claiming more absorption than the interval can hold
#: at the depth the spectrum actually shows -- the model core is deeper than the real one.
#: That is an arithmetic impossibility, not a judgement call, and it needs no allowance.
def max_stellar_sigma_kms(star: str = "solar") -> float:
    """The widest Gaussian sigma stellar physics permits, in km/s, from config.

    Thermal (+) microturbulent (+) macroturbulent (+) rotational, in quadrature. The
    mirror of `irreducible_sigma_kms`, which omits the last two so as to stay a floor.

    RYA-985: follows the star, same as the floor.

    🔴 A CEILING WITH A MISSING TERM IS NOT A CEILING — IT IS AN UNDERESTIMATE, AND USING ONE
    TO REJECT LINES REJECTS THEM FOR A REASON THAT IS NOT REAL. `55cnc_a` carries no `vmac`
    (RYA-957 adopted only what Heiter+2015 / Jofre+2014 publish, and neither publishes a
    macroturbulence; RYA-988 has since given `tau_ceti` and `eps_eri` cited RT values, and
    left 55 Cnc A to the fuller RYA-974 adoption). Dropping the term silently would shrink
    the ceiling
    and quarantine good lines as "too wide"; substituting the solar 3.8 would put another
    star's atmosphere in the bound. So this REFUSES, and the caller reports the ceiling check
    as unevaluable rather than running it on a number that does not describe this star.
    """
    p = _star(star)
    missing = [k for k in ("vmac", "vsini") if k not in p]
    if missing:
        raise MissingBroadeningTerm(
            f"{star!r} has no {' or '.join(missing)} in config/stars.yaml, so the physical "
            f"width CEILING cannot be computed. Omitting a broadening term makes the ceiling "
            f"too SMALL, which quarantines good lines as over-wide; borrowing the solar value "
            f"puts a different star's atmosphere in the bound. Adopt a cited value for "
            f"{star!r}, or run without the ceiling check and say so.")
    return math.hypot(irreducible_sigma_kms(star),
                      math.hypot(float(p["vmac"]), float(p["vsini"])))


def physical_ceiling_sigma_A(wavelength_A: float, sigma_inst_A: float, star: str = "solar") -> float:
    """The widest GAUSSIAN sigma physics allows here: instrumental (+) all Doppler terms."""
    return math.hypot(sigma_inst_A,
                      wavelength_A * max_stellar_sigma_kms(star) / C_KMS)


def gaussian_sigma_above_physical_ceiling(popt, wavelength_A: float,
                                          sigma_inst_A: float) -> bool:
    return float(popt[2]) >= physical_ceiling_sigma_A(wavelength_A, sigma_inst_A)


def over_physical_width_reason(popt, ptype: str, wavelength_A: float,
                               sigma_inst_A: float, sigma_max_A: float) -> str:
    """The exclusion sentence for a Doppler-impossible Gaussian width.

    Names `sigma_max_A` explicitly because the diagnosis a reader needs is usually not
    "too wide" but "PINNED at the optimiser's bound", and those are different faults.
    """
    sig = float(popt[2])
    ceil = physical_ceiling_sigma_A(wavelength_A, sigma_inst_A)
    pinned = abs(sig - sigma_max_A) < 1e-6
    return (f"OVER-PHYSICAL-WIDTH: the fitted {ptype} profile has Gaussian sigma "
            f"{sig:.4f} A"
            + (f", PINNED at the optimiser's upper bound {sigma_max_A:.4f} A, "
               if pinned else ", ")
            + f"against the {ceil:.4f} A that instrumental + thermal + microturbulent + "
            f"macroturbulent + rotational broadening together permit at this wavelength "
            f"(instrumental sigma {sigma_inst_A:.4f} A; vmac and vsini are counted as "
            f"Gaussian sigmas, which overstates both, so this is a ceiling and not an "
            f"expectation). No Doppler mechanism can make a line this wide, so the "
            f"integrated EW is absorption from something else, not this line's EW.")


def implied_width_A(ew_mA: float, depth: float) -> float:
    """The FWHM a Gaussian of this observed depth would need to hold this EW, in Angstrom.

    This is RYA-958's `EW / depth / 1.06` with the factor derived rather than typed. It is
    a diagnostic of the INTEGRAL: it does not describe the line, it describes how much
    area was attributed to it relative to how deep it actually is.
    """
    if not depth or depth <= 0 or not math.isfinite(depth):
        return float("nan")
    return (float(ew_mA) / 1000.0) / float(depth) / EQUIV_RECT_FACTOR


def implied_width_exceeds_ceiling(ew_mA: float, observed_depth: float,
                                  fit_window_A: float) -> bool:
    w = implied_width_A(ew_mA, observed_depth)
    if not math.isfinite(w):
        return False
    return w > implied_width_ceiling_A(fit_window_A)


def implied_width_ceiling_A(fit_window_A: float) -> float:
    """The widest equivalent rectangle the integration interval can hold — see above.

    `fit_window_A` is the FULL width of the interval `_integrate_profile` was handed
    (2 x the caller's fit half-width), not a half-width: the quantity being bounded is a
    width, so the bound is a width.
    """
    return float(fit_window_A)


def over_implied_width_reason(ew_mA: float, observed_depth: float, wavelength_A: float,
                              sigma_inst_A: float, fit_window_A: float,
                              predicted_depth: float | None) -> str:
    """The exclusion sentence for an EW inconsistent with the core it was measured on.

    Quotes the Doppler FWHM alongside the window bound because the two answer different
    questions: the window bound is what CONVICTS, the Doppler width is what tells a reader
    how far outside physics the number landed.
    """
    w = implied_width_A(ew_mA, observed_depth)
    ceil = implied_width_ceiling_A(fit_window_A)
    gauss = voigt_fwhm(physical_ceiling_sigma_A(wavelength_A, sigma_inst_A), None)
    pred = ("" if predicted_depth is None or not math.isfinite(float(predicted_depth))
            else f" (the line list predicts depth {float(predicted_depth):.3f})")
    return (f"OVER-IMPLIED-WIDTH: {ew_mA:.1f} mA of absorption was integrated onto a core "
            f"only {observed_depth:.4f} deep{pred}, which needs a Gaussian FWHM of "
            f"{w:.4f} A to hold it — wider than the {ceil:.4f} A window the model was "
            f"integrated over, so the fitted core is deeper than the observed one and the "
            f"EW is arithmetically impossible. For scale, Doppler broadening permits "
            f"{gauss:.4f} A here.")

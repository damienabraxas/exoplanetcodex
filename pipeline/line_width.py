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


def irreducible_sigma_kms() -> float:
    """Thermal+microturbulent sigma for Fe, in km/s, from config (never a literal here).

    ⚠️ Reads the SOLAR row because both callers are solar-only today. When either gains a
    `--star`, this must follow it — a floor derived from the wrong star's temperature is
    exactly the kind of silent stand-in this codebase keeps finding.
    """
    sun = STAR_PARAMS["solar"]
    return math.sqrt(_K_B * float(sun["teff"]) / _M_FE
                     + (float(sun["xi"]) * 1000.0) ** 2 / 2.0) / 1000.0


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


def physical_floor_fwhm(wavelength_A: float, sigma_inst_A: float) -> float:
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

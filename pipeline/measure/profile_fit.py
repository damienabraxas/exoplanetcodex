"""Profile-fit handler — VIS and red-optical — RYA-713 / RYA-429.

Wraps the RYA-429 fitter. It does not reimplement it: `pipeline.lines_fit` owns the
continuum renormalisation, the Voigt/Gaussian fit and the model integration, so the
project keeps ONE equivalent-width definition.

CONTROLLED: reproduces the banked HARPS Fe I pool to a median ratio of 0.971
(-0.0129 dex) over 47 lines, MAD 0.060, after three faults were fixed -- interval
integration replaced by fitting, the line width freed from the instrumental value, and
the continuum policy applied to pre-normalised atlases.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from pipeline.band_policy import BandPolicy
from pipeline.band_products import LineMeasurement
from pipeline.lines_fit import _local_renorm, _fit_profile, _integrate_profile, _ew_error
from pipeline.measure.base import MeasurementHandler, register

FWHM_TO_SIGMA = 1.0 / 2.35482
C_KMS = 299792.458
# Thermal + micro + macroturbulent sigma for solar-type metal lines. A starting guess for
# the fit, never an imposed width -- the fit moves freely between the bounds.
STELLAR_SIGMA_KMS = 1.7
FIT_HALF_A = 0.60
CONT_HALF_A = 1.20


class ProfileFitHandler(MeasurementHandler):
    method = "profile-fit"

    def widths(self, resolving_power: float, wavelength_A: float) -> tuple[float, float, float]:
        """(init, min, max) sigma in Angstrom.

        THE LINE WIDTH IS NOT THE INSTRUMENT'S. An observed line is broadened by star AND
        instrument in quadrature, and at high resolving power the star dominates entirely
        (at 5500 A: HARPS instrument 0.0203, Kitt Peak instrument 0.0047, stellar 0.0312).
        Seeding from the instrument alone started the fit 4x below truth on Kitt Peak and
        the optimiser sat on the lower bound -- Fe I 4995 returned 20.6 mA against a pool
        value of 138.

        The MINIMUM is the instrumental sigma and is a hard physical floor: nothing
        observed can be narrower than the instrument's own profile. A fit landing on it is
        reporting a width it never measured, and is quarantined rather than trusted.
        """
        sigma_inst = (wavelength_A / resolving_power) * FWHM_TO_SIGMA
        sigma_star = wavelength_A * STELLAR_SIGMA_KMS / C_KMS
        return float(np.hypot(sigma_inst, sigma_star)), float(sigma_inst), 0.40

    def measure_line(self, wav, flux, *, element, ion, wavelength_A, instrument,
                     policy: BandPolicy, pre_normalised: bool,
                     context: dict[str, Any]) -> LineMeasurement:
        c = float(wavelength_A)
        R = float(context["resolving_power"])

        def quarantine(reason: str) -> LineMeasurement:
            lm = LineMeasurement(element=element, ion=ion, wavelength_air_A=c,
                                 instrument=instrument, ew_mA=float("nan"),
                                 ew_method=f"PROFILE-FIT attempted in {policy.name}")
            lm.in_aggregate = False
            lm.excluded_reason = reason
            return lm

        m = np.abs(wav - c) <= CONT_HALF_A
        wf, ff = wav[m], flux[m]
        if wf.size < 40:
            return quarantine(f"COVERAGE: only {wf.size} points within "
                              f"+/-{CONT_HALF_A} A of the line")

        # Continuum. A pre-normalised atlas already IS residual flux -- unity is the
        # continuum by construction, and re-fitting through percentile-filtered edge
        # strips lands below unity in crowded spectrum because those strips contain
        # lines. That was worth -0.088 dex on its own.
        fn = ff if pre_normalised else _local_renorm(wf, ff, c, window=CONT_HALF_A)[0]

        s_init, s_min, s_max = self.widths(R, c)
        fit_m = np.abs(wf - c) <= FIT_HALF_A
        popt, pcov, ptype, chi2 = _fit_profile(
            wf[fit_m], fn[fit_m], c, sigma_init=s_init, sigma_min=s_min,
            sigma_max=s_max, core_half_A=max(3 * s_init, 0.03))
        if popt is None or ptype == "failed":
            return quarantine("FIT-FAILED: no Voigt or Gaussian solution converged")

        ew = _integrate_profile(wf[fit_m], popt, ptype)
        err = _ew_error(fn[~fit_m], ew, popt, pcov, ptype)
        lm = LineMeasurement(
            element=element, ion=ion, wavelength_air_A=c, instrument=instrument,
            ew_mA=float(ew),
            ew_method=(f"PROFILE-FIT ({ptype}), model integrated; chi2_red={chi2:.4g}; "
                       f"sigma={popt[2]:.4f} A (init {s_init:.4f}, floor {s_min:.4f}); "
                       f"ew_err={err:.2f} mA; continuum="
                       f"{'atlas (pre-normalised)' if pre_normalised else 'local linear'}"))
        if abs(float(popt[2]) - s_min) < 1e-4:
            lm.in_aggregate = False
            lm.excluded_reason = (
                f"FIT-PINNED: fitted sigma {popt[2]:.4f} A sits on the instrumental floor "
                f"{s_min:.4f} A. Nothing observed can be narrower than the instrument "
                f"profile, so the optimiser hit a bound rather than a minimum.")
        return lm


register(ProfileFitHandler())

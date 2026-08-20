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
from pipeline.line_width import (gamma_of, total_width_below_physical_floor,
                                 under_physical_width_reason)
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

        The MINIMUM is the instrumental sigma. It is a legitimate BOUND on the optimiser —
        no observed feature is narrower than the instrument profile — but 🔴 landing on it
        is NOT grounds for quarantine, and this docstring used to say it was. In a Voigt
        fit sigma and gamma are degenerate, so a railed sigma says where the optimiser
        resolved the degeneracy, not how wide the line is. `pipeline.line_width` judges the
        TOTAL width instead (RYA-906/911).
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
        # The fit's own numbers become columns rather than prose. They were computed here,
        # used to decide the verdict below, and thrown away — the defect RYA-911 closed on
        # the band harness and left standing on this site. `profile_sigma_A` is an ANGSTROM
        # and is NOT `sigma_A`, which is one sigma on A in DEX from the chi2 curvature
        # (RYA-847). One character apart, no shared units.
        lm.profile_sigma_A = float(popt[2])
        lm.profile_gamma_A = gamma_of(popt, ptype)
        lm.profile_sigma_floor_A = float(s_min)
        lm.red_chi2 = float(chi2)

        # 🔴 RYA-906/911 — THIS TEST WAS `abs(popt[2] - s_min) < 1e-4` AND IT WAS WRONG.
        #
        # PR #315 refuted it on the band harness and did not reach here, because the
        # corrected physics lived in a script this module cannot import. Measured there:
        # the sigma-only guard fired on 25 HARPS Fe II lines and EVERY ONE was a Voigt fit
        # — never once a Gaussian. Fe II 6084.102 had sigma exactly on the bound and was
        # rejected while gamma 0.0591 put its total FWHM at 0.1395 A, squarely among the
        # lines that were KEPT.
        #
        # sigma and gamma are DEGENERATE: both broaden the line and the optimiser trades
        # one against the other at fixed total width. So ask the question the guard always
        # meant to ask — is the TOTAL width below what physics permits?
        if total_width_below_physical_floor(popt, ptype, c, s_min):
            lm.in_aggregate = False
            lm.excluded_reason = under_physical_width_reason(popt, ptype, c, s_min)
        return lm


register(ProfileFitHandler())

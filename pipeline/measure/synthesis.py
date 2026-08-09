"""Synthesis handler — near-UV and NIR — RYA-713.

Measures an abundance by fitting the OBSERVED FLUX in a window, not by integrating a line.
That distinction is the whole reason this handler exists: in the near-UV the median line
separation is 0.146 Å, smaller than a strong line's own wings, so no interval contains one
profile and excludes its neighbours, and there is no isolated profile to fit. Only a method
that models every contributor in the window at once is valid there.

WHAT IT ACTUALLY DOES
---------------------
For each line: synthesise the window at a grid of trial A(X), compare to the observed
flux, and take the χ² minimum as the measurement. The target element's abundance is the
only free parameter — every other species stays at the converged composition — so the
blends are *modelled* rather than avoided.

It reuses `abundances_derive._synth_flux_at_abund`, the shared synth-v2 spectrum generator
(RYA-287). Building a second synthesiser here would give the project two definitions of a
synthetic spectrum that drift apart, which is the RYA-701 failure.

THE PSEUDO-CONTINUUM, WHICH IS THE NEAR-UV'S REAL PROBLEM
----------------------------------------------------------
Measured on the Kitt Peak atlas: near-UV median flux runs 0.283–0.805 and the 95th
percentile only reaches 0.80–0.98. **The true continuum is never observed** — blanketing
means no wavelength in the window is unabsorbed.

So the observed spectrum has been normalised to a *pseudo*-continuum (an envelope through
the least-absorbed pixels), while a synthetic spectrum comes out of Turbospectrum on the
*true* continuum. Comparing them directly compares two different normalisations, and the
abundance absorbs the difference.

The treatment here is to apply **the same envelope operation to the synthetic spectrum**
before comparing. Both sides then carry the same pseudo-continuum definition and it
divides out of the comparison. What does not divide out is that the envelope depends on
the synthetic itself, so the systematic is recorded on every measurement rather than being
silently absorbed — `PSEUDO_CONTINUUM_SYSTEMATIC_NOTE` travels with the result.

CONTROL
-------
This handler has NOT been controlled. `assert_controlled()` will refuse to let it run in
the near-UV or NIR until it has reproduced the known optical answer in the VIS. The profile
fitter's licence does not transfer: synthesis fails in entirely different ways — incomplete
line lists, wrong blend abundances, broadening, and the pseudo-continuum above.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from pipeline.band_policy import BandPolicy
from pipeline.band_products import LineMeasurement
from pipeline.measure.base import MeasurementHandler, register

# Half-width of the synthesis window. Wider than a profile-fit window on purpose: the
# point is to include the neighbours so they can be modelled, not to exclude them.
SYNTH_HALF_A = 1.00
# Trial abundance grid. Coarse pass then a refined pass around the minimum.
COARSE_STEP_DEX = 0.20
FINE_STEP_DEX = 0.02
DEFAULT_SPAN_DEX = 1.0

PSEUDO_CONTINUUM_SYSTEMATIC_NOTE = (
    "near-UV: the true continuum is not observed (median flux 0.283-0.805). Observed and "
    "synthetic spectra are both reduced by the SAME envelope so the definition divides "
    "out of the comparison, but the envelope depends on the synthetic, so a residual "
    "pseudo-continuum systematic remains and does NOT average down with more lines.")


def envelope(wav: np.ndarray, flux: np.ndarray, *, frac: float = 0.95,
             n_knots: int = 5) -> np.ndarray:
    """Pseudo-continuum: a smooth curve through the least-absorbed pixels.

    Piecewise so it can follow real curvature, percentile-based so a single hot pixel
    cannot lift it. Applied IDENTICALLY to observed and synthetic spectra — that identity
    is the point, not the particular choice of percentile.
    """
    if wav.size < n_knots * 3:
        return np.full_like(flux, float(np.percentile(flux, frac * 100)))
    edges = np.linspace(wav[0], wav[-1], n_knots + 1)
    xs, ys = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (wav >= lo) & (wav <= hi)
        if m.sum() >= 3:
            xs.append(0.5 * (lo + hi))
            ys.append(float(np.percentile(flux[m], frac * 100)))
    if len(xs) < 2:
        return np.full_like(flux, float(np.percentile(flux, frac * 100)))
    return np.interp(wav, np.array(xs), np.array(ys))


class SynthesisHandler(MeasurementHandler):
    method = "synthesis"

    def __init__(self) -> None:
        super().__init__()
        self._synth = None
        self._tmp_dir = '/tmp/ispec_codex_synth'

    def prepare(self, policy: BandPolicy, context: dict[str, Any]) -> None:
        super().prepare(policy, context)
        missing = [k for k in ("atmosphere", "teff", "logg", "feh", "vturb",
                               "linelist", "isotopes", "solar_abund", "atom_code",
                               "resolving_power") if k not in context]
        if missing:
            raise RuntimeError(
                f"SynthesisHandler needs {missing} in the run context. These come from the "
                f"converged stellar composition and the iSpec setup; they are not defaults "
                f"this handler may invent, because inventing them would silently change "
                f"what the fit is measuring.")
        if policy.telluric_required and not context.get("telluric_corrected"):
            raise RuntimeError(
                f"the {policy.name} band requires telluric correction and the context does "
                f"not declare `telluric_corrected`. Before correction the observed flux is "
                f"not a stellar spectrum, so a fit to it does not measure the star.")
        from pipeline.abundances_derive import _synth_flux_at_abund
        # iSpec writes intermediate atmosphere/linelist files here and does NOT create the
        # directory; without it every synthesis dies with a FileNotFoundError naming a
        # random temp file, which reads like a corrupt install rather than a missing dir.
        import os
        tmp = context.get("tmp_dir", "/tmp/ispec_codex_synth")
        os.makedirs(tmp, exist_ok=True)
        self._tmp_dir = tmp
        self._synth = _synth_flux_at_abund

    # ── the fit ───────────────────────────────────────────────────────────────
    def _chi2(self, wav_nm, obs, ctx, element, trial_A, pseudo: bool) -> tuple[float, np.ndarray]:
        syn = self._synth(
            wav_nm, ctx["atmosphere"], ctx["teff"], ctx["logg"], ctx["feh"], ctx["vturb"],
            ctx["linelist"], ctx["isotopes"], ctx["solar_abund"], element,
            ctx["atom_code"], float(trial_A),
            R=float(ctx["resolving_power"]),
            macroturbulence=float(ctx.get("macroturbulence", 0.0)),
            vsini=float(ctx.get("vsini", 0.0)),
            tmp_dir=self._tmp_dir)
        if pseudo:
            # Same operation on both sides. The observed spectrum reaching us has already
            # been divided by an envelope; dividing the synthetic by ITS envelope puts
            # them on one footing instead of comparing two normalisations.
            syn = syn / np.clip(envelope(wav_nm, syn), 1e-6, None)
        resid = obs - syn
        return float(np.sum(resid ** 2) / max(resid.size - 1, 1)), syn

    def measure_line(self, wav, flux, *, element, ion, wavelength_A, instrument,
                     policy: BandPolicy, pre_normalised: bool,
                     context: dict[str, Any]) -> LineMeasurement:
        c = float(wavelength_A)

        def quarantine(reason: str) -> LineMeasurement:
            lm = LineMeasurement(element=element, ion=ion, wavelength_air_A=c,
                                 instrument=instrument, ew_mA=float("nan"),
                                 ew_method=f"SYNTHESIS attempted in {policy.name}")
            lm.in_aggregate = False
            lm.excluded_reason = reason
            return lm

        if self._synth is None:
            return quarantine("SETUP: prepare() was not called, so no synthesiser is bound")

        m = np.abs(wav - c) <= SYNTH_HALF_A
        w, o = wav[m], flux[m]
        if w.size < 30:
            return quarantine(f"COVERAGE: only {w.size} points within "
                              f"+/-{SYNTH_HALF_A} A of the line")

        w_nm = w / 10.0
        pseudo = "pseudo" in policy.continuum_treatment.lower()
        if pseudo:
            o = o / np.clip(envelope(w, o), 1e-6, None)

        a0 = float(context.get("a_start", context["solar_A"]))
        span = float(context.get("span_dex", DEFAULT_SPAN_DEX))

        try:
            coarse = np.arange(a0 - span, a0 + span + 1e-9, COARSE_STEP_DEX)
            cvals = [self._chi2(w_nm, o, context, element, a, pseudo)[0] for a in coarse]
            k = int(np.argmin(cvals))
            if k in (0, len(coarse) - 1):
                return quarantine(
                    f"FIT-RAILED: chi2 minimum sits at the edge of the trial range "
                    f"[{coarse[0]:.2f}, {coarse[-1]:.2f}], so the true minimum is outside "
                    f"it. Widen the range rather than accepting an edge value.")
            fine = np.arange(coarse[k] - COARSE_STEP_DEX, coarse[k] + COARSE_STEP_DEX + 1e-9,
                             FINE_STEP_DEX)
            fvals = [self._chi2(w_nm, o, context, element, a, pseudo)[0] for a in fine]
            j = int(np.argmin(fvals))
            a_best = float(fine[j])
            chi2_min = float(fvals[j])
            # Uncertainty from the chi2 curvature: the abundance step that raises chi2 by
            # its own rms. Not a formal covariance -- it is the width of the minimum, and
            # it is honest about being that.
            arr = np.asarray(fvals)
            thresh = chi2_min + (chi2_min / max(len(o) - 1, 1)) ** 0.5 * chi2_min
            inside = np.where(arr <= max(thresh, chi2_min * 1.01))[0]
            sigma = float(FINE_STEP_DEX * max(len(inside), 1) / 2.0)
            _, syn_best = self._chi2(w_nm, o, context, element, a_best, pseudo)
        except Exception as e:
            return quarantine(f"SYNTH-FAILED: {type(e).__name__}: {str(e)[:140]}")

        # Synthetic EW at the best fit -- not the measurement, but it makes this handler
        # comparable to the profile fitter in the optical control.
        from pipeline._numcompat import trapezoid as _trapz
        ew = float(_trapz(np.clip(1.0 - syn_best, 0.0, None), w)) * 1000.0

        lm = LineMeasurement(
            element=element, ion=ion, wavelength_air_A=c, instrument=instrument,
            ew_mA=ew, abundance=a_best,
            ew_method=(f"SYNTHESIS chi2 flux fit over +/-{SYNTH_HALF_A} A; "
                       f"A={a_best:.3f} +/- {sigma:.3f}; chi2_red={chi2_min:.4g}; "
                       f"continuum={'pseudo (envelope applied to BOTH sides)' if pseudo else 'true'}; "
                       f"R={context['resolving_power']:.0f}; EW is SYNTHETIC at best fit, "
                       f"not the measurement"
                       + (f" | SYSTEMATIC: {PSEUDO_CONTINUUM_SYSTEMATIC_NOTE}" if pseudo else "")))
        return lm


register(SynthesisHandler())

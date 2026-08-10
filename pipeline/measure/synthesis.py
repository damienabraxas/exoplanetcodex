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
        self._weights = None

    def prepare(self, policy: BandPolicy, context: dict[str, Any]) -> None:
        super().prepare(policy, context)
        # macroturbulence and vsini are REQUIRED, not defaulted (RYA-288).
        #
        # `_synth_flux_at_abund` defaults them to 0 because v1 EW matching is
        # broadening-INVARIANT -- convolution preserves equivalent width. A flux fit is
        # not: with macro=0 the synthetic core is too deep and too narrow, so chi2 is
        # minimised by LOWERING the abundance. That is exactly what the first control
        # showed -- A = 7.046 against a banked 7.466, -0.42 dex, with the EW low too
        # (ratio 0.822). Calling the low-level function directly inherited the EW-mode
        # defaults and bypassed the project's own no-silent-broadening guard.
        missing = [k for k in ("atmosphere", "teff", "logg", "feh", "vturb",
                               "linelist", "isotopes", "solar_abund", "atom_code",
                               "resolving_power", "macroturbulence", "vsini")
                   if k not in context]
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
    def _response(self, wav_nm, ctx, element, a_centre: float, pseudo: bool) -> np.ndarray:
        """Per-pixel sensitivity |d flux / d A(X)| -- which pixels carry information.

        WHY THIS EXISTS. A flat chi2 over the whole window is dominated by whatever
        occupies most of it, and for a WEAK line that is continuum and neighbours, not
        the line being measured. The first controlled run showed it plainly:
        corr(EW, abundance offset) = +0.891, with 16-21 mA lines landing 0.46-0.48 dex
        LOW while 44-59 mA lines came out right. The window was answering a different
        question for weak lines than for strong ones.

        Weighting by the response uses the pixels that actually constrain A(X). This is
        not a tuning knob: the weights come from the model's own derivative, computed
        fresh per line, and no threshold is chosen to make an answer come out.
        """
        step = 0.15
        lo = self._synth_at(wav_nm, ctx, element, a_centre - step, pseudo)
        hi = self._synth_at(wav_nm, ctx, element, a_centre + step, pseudo)
        r = np.abs(hi - lo)
        peak = float(r.max())
        if not np.isfinite(peak) or peak <= 0:
            return np.ones_like(r)          # no sensitivity anywhere: fall back to flat
        return r / peak

    def _synth_at(self, wav_nm, ctx, element, trial_A, pseudo: bool) -> np.ndarray:
        syn = self._synth(
            wav_nm, ctx["atmosphere"], ctx["teff"], ctx["logg"], ctx["feh"], ctx["vturb"],
            ctx["linelist"], ctx["isotopes"], ctx["solar_abund"], element,
            ctx["atom_code"], float(trial_A),
            R=float(ctx["resolving_power"]),
            macroturbulence=float(ctx["macroturbulence"]),
            vsini=float(ctx["vsini"]),
            tmp_dir=self._tmp_dir)
        # Same operation as the observed side, unconditionally -- that identity is the
        # point. `pseudo` no longer gates it; it only flags the recorded systematic.
        syn = syn / np.clip(envelope(wav_nm, syn), 1e-6, None)
        return syn

    def _chi2(self, wav_nm, obs, ctx, element, trial_A, pseudo: bool) -> tuple[float, np.ndarray]:
        # Same envelope operation on both sides when the band has no true continuum:
        # the observed spectrum reaching us is already divided by ITS envelope, so
        # dividing the synthetic by its own puts them on one footing rather than
        # comparing two normalisations.
        syn = self._synth_at(wav_nm, ctx, element, trial_A, pseudo)
        resid = obs - syn
        w = self._weights if self._weights is not None else np.ones_like(resid)
        n_eff = max(float(w.sum()), 1e-9)
        return float(np.sum(w * resid ** 2) / n_eff), syn

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
        # THE ENVELOPE IS APPLIED IN EVERY BAND, not only where the continuum is called
        # "pseudo" (RYA-713). The band policy distinction is about whether a TRUE
        # continuum is REACHABLE; the need for both sides to carry the SAME normalisation
        # is universal, and skipping it in the optical is what kept the control failing.
        #
        # Measured over +/-1 A in the VIS: observed p95 ~0.994 against synthetic p95
        # ~0.999. The atlas continuum sits ~0.5% BELOW the synthetic's true continuum
        # while the synthetic lines run deeper (mean syn-obs -0.0050). A flux chi2
        # reconciles both by LOWERING the abundance -- which is exactly the -0.20 dex the
        # control kept reporting, on both angles once they were made comparable.
        #
        # `pseudo` is retained only to decide whether the residual systematic must be
        # RECORDED on the result: in the near-UV the envelope is standing in for a
        # continuum that is never observed, and that does not average down.
        pseudo = "pseudo" in policy.continuum_treatment.lower()
        o = o / np.clip(envelope(w, o), 1e-6, None)

        a0 = float(context.get("a_start", context["solar_A"]))
        span = float(context.get("span_dex", DEFAULT_SPAN_DEX))

        try:
            # Response weights from the model's own derivative at the starting guess.
            self._weights = self._response(w_nm, context, element, a0, pseudo)
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
            # A fit that never leaves its starting guess has not measured anything -- the
            # chi2 surface was flat, or the synthesis did not respond. Reporting the seed
            # as a result is how a control silently passes on lines it never constrained.
            if abs(a_best - a0) < FINE_STEP_DEX * 0.51:
                return quarantine(
                    f"FIT-DID-NOT-MOVE: best A {a_best:.3f} is the starting guess "
                    f"{a0:.3f} to within one grid step. The chi2 surface is flat here, so "
                    f"this line does not constrain the abundance and its value is the "
                    f"seed, not a measurement.")
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

        # Synthetic EW OF THE TARGET LINE ONLY -- not the measurement, but it is what
        # makes this handler comparable to the profile fitter in the optical control.
        #
        # Integrating the whole window would count every neighbour too, and the reference
        # EW it is compared against is a DEBLENDED single-line value from the profile
        # fitter. That mismatch is why the first controlled run showed EW ratios from
        # 0.64 to 1.47 while the abundances looked sane -- the two sides were measuring
        # different things. Differencing against a synthesis with the element suppressed
        # isolates this element's contribution.
        from pipeline._numcompat import trapezoid as _trapz
        try:
            syn_without = self._synth_at(w_nm, context, element, a_best - 4.0, pseudo)
            target_only = np.clip(syn_without - syn_best, 0.0, None)
            ew = float(_trapz(target_only, w)) * 1000.0
        except Exception:
            # Fall back to the whole-window integral, and SAY that it is one.
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

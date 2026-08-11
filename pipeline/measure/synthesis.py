"""Synthesis handler — a THIN ADAPTER over the validated synth-v2 fit — RYA-713.

Ryan, 2026-08-09: *"make it a thin adapter over `_fit_synth_flux`."*

WHY THIS FILE IS SMALL NOW
--------------------------
It used to reimplement the fit: its own χ² loop, its own trial grid, its own window rule,
its own continuum handling. Every one of those was an invention, and each was wrong in a
way that cost a control run to find:

* a fixed ±1.00 Å window where synth-v2 uses wing-wide ±0.237 Å at 45 mÅ — 4× too wide,
  which dilutes the line among blends and continuum;
* response-weighted χ², which breaks the σ-invariance synth-v2 deliberately relies on
  (*"σ is constant, so the χ² minimum location — i.e. the fitted A(X) — is
  σ-scale-invariant"*), and which did not help anyway;
* an envelope applied in every band: harmless at 0.400 Å knot spacing, actively
  self-cancelling at 0.095 Å, because the *synthetic's* envelope deepens with the trial
  abundance and divides out the very sensitivity being fitted (HARPS: A = 7.956, +0.436);
* `macroturbulence = 0`, inherited from EW-mode defaults, bypassing `_resolve_broadening`
  (RYA-288) — correct for EW, which is broadening-invariant, wrong for a flux fit.

`pipeline.abundances_derive._fit_synth_flux` already does this correctly and is banked at
**A(Fe I) = 7.520 on 23 HARPS lines**. So this class no longer computes an abundance. It
resolves policy, checks the synthesiser can see the line, delegates, and translates.

WHAT REMAINS HERE, AND WHY
--------------------------
Only what `_fit_synth_flux` does not know about, because it is the handler's job:

1. **Band policy** — is synthesis legal here, and does the band need telluric correction
   before the observed flux is a stellar spectrum at all?
2. **Coverage** — a line absent from the synthesis line list gives the SAME spectrum at
   every abundance, so χ² is flat and the "fit" returns noise. Measured: 4065.381 Å sits
   below the GES list's 4200 Å blue limit and synthesised to EW 0.00 against a pool value
   of 73.59. That is a coverage gap and is named as one.
3. **Quarantine translation** — `_fit_synth_flux` returns `status ∈ {ok, edge_pinned,
   failed}`; those become retained-with-reason `LineMeasurement`s (RYA-711), never a
   raised exception, because an exception loses the line silently.
4. **The near-UV pseudo-continuum note** — the true continuum is never observed there
   (median flux 0.283–0.805) and that systematic does not average down.

CONTROL
-------
Still gated. `assert_controlled()` refuses the near-UV and NIR until this adapter
reproduces the known optical answer. Delegating to validated code makes that likely; it
is not a substitute for running it.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from pipeline.band_policy import BandPolicy
from pipeline.band_products import LineMeasurement
from pipeline.measure.base import MeasurementHandler, register

# Bracket for the bounded minimiser. Wide enough that a real solar value is never at the
# edge; `edge_pinned` is reported rather than quietly accepted.
A_LO, A_HI = 4.0, 10.0

# Peak |Δflux| across a ±0.5 dex swing below which the line is not usably in the list.
MIN_SENSITIVITY = 0.002

PSEUDO_CONTINUUM_SYSTEMATIC_NOTE = (
    "near-UV: the true continuum is not observed (median flux 0.283-0.805), so the "
    "abundance carries a pseudo-continuum systematic that does NOT average down with "
    "more lines and must be stated separately from the line-to-line scatter.")


class SynthesisHandler(MeasurementHandler):
    method = "synthesis"

    def __init__(self) -> None:
        super().__init__()
        self._fit = None
        self._synth = None
        self._window = None
        self._tmp_dir = "/tmp/ispec_codex_synth"

    def prepare(self, policy: BandPolicy, context: dict[str, Any]) -> None:
        super().prepare(policy, context)
        # Broadening is REQUIRED, never defaulted (RYA-288): a flux fit is not
        # broadening-invariant the way EW matching is, and wrong broadening biases A(X)
        # itself rather than merely its error bar.
        missing = [k for k in ("atmosphere", "teff", "logg", "feh", "vturb",
                               "linelist", "isotopes", "solar_abund", "atom_code",
                               "resolving_power", "macroturbulence", "vsini")
                   if k not in context]
        if missing:
            raise RuntimeError(
                f"SynthesisHandler needs {missing} in the run context. These come from the "
                f"converged composition and the iSpec setup; inventing any of them would "
                f"silently change what the fit measures.")
        if policy.telluric_required and not context.get("telluric_corrected"):
            raise RuntimeError(
                f"the {policy.name} band requires telluric correction and the context does "
                f"not declare `telluric_corrected`. Before correction the observed flux is "
                f"not a stellar spectrum.")

        import os
        from pipeline.abundances_derive import (
            _fit_synth_flux, _synth_flux_at_abund, _wingwide_window_nm)
        self._tmp_dir = context.get("tmp_dir", self._tmp_dir)
        os.makedirs(self._tmp_dir, exist_ok=True)   # iSpec does not create it
        self._fit = _fit_synth_flux
        self._synth = _synth_flux_at_abund
        self._window = _wingwide_window_nm

    def _synth_kw(self, ctx: dict, element: str) -> dict:
        return dict(atmosphere=ctx["atmosphere"], teff=ctx["teff"], logg=ctx["logg"],
                    feh=ctx["feh"], vturb=ctx["vturb"], linelist=ctx["linelist"],
                    isotopes=ctx["isotopes"], solar_abund=ctx["solar_abund"],
                    element=element, atom_code=ctx["atom_code"],
                    R=float(ctx["resolving_power"]),
                    macroturbulence=float(ctx["macroturbulence"]),
                    vsini=float(ctx["vsini"]), tmp_dir=self._tmp_dir)

    def measure_line(self, wav, flux, *, element, ion, wavelength_A, instrument,
                     policy: BandPolicy, pre_normalised: bool,
                     context: dict[str, Any]) -> LineMeasurement:
        c = float(wavelength_A)

        def quarantine(reason: str, ew=float("nan")) -> LineMeasurement:
            lm = LineMeasurement(element=element, ion=ion, wavelength_air_A=c,
                                 instrument=instrument, ew_mA=ew,
                                 ew_method=f"SYNTHESIS ({policy.name}) via _fit_synth_flux")
            lm.in_aggregate = False
            lm.excluded_reason = reason
            return lm

        if self._fit is None:
            return quarantine("SETUP: prepare() was not called")

        # Wing-wide window from the line's own strength -- synth-v2's rule, not a constant.
        ew_hint = float(context.get("ew_hint_mA", 40.0))
        wave_base, wave_top = self._window(c / 10.0, ew_hint)
        hw_A = (wave_top - wave_base) * 10.0 / 2.0

        obs_w_nm = np.asarray(wav, dtype=float) / 10.0
        obs_f = np.asarray(flux, dtype=float)
        m = (obs_w_nm >= wave_base) & (obs_w_nm <= wave_top)
        if m.sum() < 5:
            return quarantine(f"COVERAGE: {int(m.sum())} observed pixels in the "
                              f"+/-{hw_A:.3f} A wing-wide window")

        kw = self._synth_kw(context, element)
        a0 = float(context.get("a_start", context.get("solar_A", 7.5)))
        try:
            lo = self._synth(obs_w_nm[m], trial_A=a0 - 0.5, **kw)
            hi = self._synth(obs_w_nm[m], trial_A=a0 + 0.5, **kw)
            sens = float(np.max(np.abs(hi - lo)))
        except Exception as e:
            return quarantine(f"SYNTH-FAILED probing sensitivity: {type(e).__name__}: {e}")
        if sens < MIN_SENSITIVITY:
            return quarantine(
                f"NOT-IN-SYNTH-LINELIST: a +/-0.5 dex abundance swing moves the synthetic "
                f"flux by at most {sens:.5f}. This line is absent from the synthesis line "
                f"list, or carries a gf too weak to register, so no abundance can be "
                f"fitted from it. Coverage gap, not a measurement failure.")

        # ── delegate to the validated fit ─────────────────────────────────────
        res = self._fit(
            obs_w_nm, obs_f, context["atmosphere"], context["teff"], context["logg"],
            context["feh"], context["vturb"], context["linelist"], context["isotopes"],
            context["solar_abund"], element, context["atom_code"],
            wave_base, wave_top, A_LO, A_HI,
            float(context["resolving_power"]), float(context["macroturbulence"]),
            float(context["vsini"]), tmp_dir=self._tmp_dir)

        if res["status"] == "failed":
            return quarantine(f"FIT-FAILED: {res.get('reason', 'unknown')}")
        a_best = float(res["A_X"])

        # Synthetic EW of the TARGET LINE ONLY, for comparability with the profile fitter.
        # Integrated over a range WIDER than the fit window: the fit window is deliberately
        # tight so the line dominates chi2, but a DIFFERENCE of two syntheses is zero
        # wherever this element does not absorb, so a wide range adds the damping wings and
        # nothing else -- blends cancel. Integrating over the fit window instead truncated
        # the wings and drove the EW angle to a median ratio of 0.250 by construction.
        ew = float("nan")
        try:
            from pipeline._numcompat import trapezoid as _trapz
            ew_hw = max(3.0 * hw_A, 1.0)
            ew_w = np.linspace(c - ew_hw, c + ew_hw, max(int(2 * ew_hw / 0.004), 200))
            with_el = self._synth(ew_w / 10.0, trial_A=a_best, **kw)
            without = self._synth(ew_w / 10.0, trial_A=a_best - 4.0, **kw)
            ew = float(_trapz(np.clip(without - with_el, 0.0, None), ew_w)) * 1000.0
        except Exception:
            pass

        pseudo = "pseudo" in policy.continuum_treatment.lower()
        method = (f"SYNTHESIS via _fit_synth_flux (synth-v2, RYA-287); A={a_best:.3f}; "
                  f"red_chi2={res['red_chi2']}; n_pix={res['n_pix']}; window "
                  f"+/-{hw_A:.3f} A wing-wide from EW {ew_hint:.0f} mA; "
                  f"R={context['resolving_power']:.0f}, vmac={context['macroturbulence']}, "
                  f"vsini={context['vsini']}; EW is the SYNTHETIC target-line EW over "
                  f"+/-{ew_hw:.2f} A, not the measurement"
                  + (f" | SYSTEMATIC: {PSEUDO_CONTINUUM_SYSTEMATIC_NOTE}" if pseudo else ""))

        # FIT-QUALITY GATE (RYA-713). synth-v2's own comment: "the gate is FIT QUALITY
        # (reduced chi2), applied in the aggregation below. Every line is fit; acceptance
        # is on merit." My adapter accepted every status=='ok' line regardless, so badly
        # fitted lines entered the median and inflated line-to-line scatter to 1.181 --
        # against A_X_std = 0.066 for the banked 23-line synth-v2 result.
        #
        # Fitting every line and gating on merit is the right shape: the line is measured
        # and RETAINED with its chi2, it simply does not enter the aggregate (RYA-711).
        from config.constants import SYNTH_CHI2_GATE
        if float(res["red_chi2"]) >= SYNTH_CHI2_GATE:
            return quarantine(
                f"CHI2-GATE: red_chi2 {res['red_chi2']} >= {SYNTH_CHI2_GATE}. The fit "
                f"converged (A={a_best:.3f}) but does not describe the observed flux well "
                f"enough to carry an abundance. Measured and retained; excluded from the "
                f"aggregate only.", ew=ew)

        if res["status"] == "edge_pinned":
            return quarantine(
                f"FIT-EDGE-PINNED: A={a_best:.3f} sits on the [{A_LO}, {A_HI}] bracket "
                f"edge, so the true minimum lies outside it. Reported, never substituted.",
                ew=ew)

        return LineMeasurement(element=element, ion=ion, wavelength_air_A=c,
                               instrument=instrument, ew_mA=ew, abundance=a_best,
                               ew_method=method)


register(SynthesisHandler())

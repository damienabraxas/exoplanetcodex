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
from pipeline import synth_ew                       # RYA-878
from pipeline.band_products import LineMeasurement
from pipeline.fit_constraint import as_float_or_none as _as_f  # RYA-847
from pipeline.constraint_gate import verdict as _constraint_verdict  # RYA-847
from pipeline.measure.base import MeasurementHandler, register
# RYA-770/765: every per-line accept/reject is recorded on the intake timeline. A
# no-op unless a trace is active, so this costs nothing in production (RYA-429
# discipline: a dropped line must never be silent).
from pipeline import intake_debug as dbg

# Bracket for the bounded minimiser. Wide enough that a real solar value is never at the
# edge; `edge_pinned` is reported rather than quietly accepted.
A_LO, A_HI = 4.0, 10.0

# Peak |Δflux| across a ±0.5 dex swing below which the line is not usably in the list.
MIN_SENSITIVITY = 0.002

# ── RYA-770: iSpec's edge-zeroing, and why the adapter must route around it ────
# `ispec.generate_spectrum` returns EXACTLY 0.0 for the first synthesis pixel and the
# last 3-4, on every window (measured: index 0, plus indices n-4..n-1 or n-3..n-1). It
# is a property of the synthesiser's output grid, not of the spectrum.
#
# `_fit_synth_flux` builds a 0.0002 nm grid over the fit window, synthesises on it, and
# interpolates onto the observed pixels. Any observed pixel inside a zeroed span is
# therefore compared against ~0 flux instead of ~1, contributing (1/0.01)^2 = 1e4 to
# chi2 by itself. With sigma = 0.01 and a few dozen pixels in a wing-wide window, ONE
# such pixel dominates the sum.
#
# Measured over the 12 optical control lines: median red_chi2 92.17 as-is, 2.77 with
# only those pixels removed -- against a banked med red_chi2 of 2.976. Individual
# lines moved 301.59 -> 0.84, 162.71 -> 2.72, 16.58 -> 1.80. That is the whole of the
# adapter's CHI2-GATE failure.
#
# The fix lives HERE and not in `_fit_synth_flux` deliberately (RYA-770 §5/§6): the
# core is the validated path banked at A(Fe I) = 7.520 and must not move in this
# ticket. The adapter owns the observed arrays it hands over, so it drops the observed
# pixels that would land on a zeroed synthetic pixel and passes the rest. The core is
# called unchanged and never sees them.
#
# The span is MEASURED per line, never hardcoded: the trailing count varies with the
# window, and a constant would silently rot the day iSpec changes it.
EDGE_ZERO_FLUX = 1e-6      # below this is the ~0.0 sentinel, not a real flux


def clean_span(sw, *flux_arrays, zero_flux: float = EDGE_ZERO_FLUX):
    """(clean_lo, clean_hi, n_zeroed) — the span of `sw` carrying real synthetic flux.

    A pixel counts as clean only if EVERY supplied flux array is above the sentinel
    there, so a span that is real at one trial abundance and zeroed at another is
    treated as zeroed. Returns (None, None, n) when nothing is clean.

    Pure geometry, no iSpec: this is the part of the RYA-770 edge trim that can be
    tested without a synthesiser, which is why it is a module function rather than
    inline in `measure_line`.
    """
    sw = np.asarray(sw, dtype=float)
    good = np.ones(sw.shape, dtype=bool)
    for f in flux_arrays:
        good &= np.asarray(f, dtype=float) > zero_flux
    n_zeroed = int((~good).sum())
    if not good.any():
        return None, None, n_zeroed
    first = int(np.argmax(good))
    last = len(good) - 1 - int(np.argmax(good[::-1]))
    return float(sw[first]), float(sw[last]), n_zeroed


def edge_pixels_to_drop(obs_w, wave_base: float, wave_top: float,
                        clean_lo: float, clean_hi: float):
    """Mask of observed pixels inside the fit window but OUTSIDE the clean span.

    Only in-window pixels can matter: `_fit_synth_flux` masks the rest away itself, so
    dropping an out-of-window pixel would change nothing and would misreport the count.
    """
    obs_w = np.asarray(obs_w, dtype=float)
    in_window = (obs_w >= wave_base) & (obs_w <= wave_top)
    return in_window & ((obs_w < clean_lo) | (obs_w > clean_hi))

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
        self._wstep = None
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
        # TELLURIC (RYA-786). This used to gate on the BAND's telluric_required plus a
        # `telluric_corrected` key in the context, which re-asked a question already
        # answered in three places (RYA-424 routing, the instrument catalog, the RYA-460
        # per-line regime decision) and made every new band run re-collide with it.
        #
        # It also invited the wrong fix: declaring `telluric_corrected` to get past the
        # gate asserts a correction that was never applied. For a reference atlas the
        # honest basis is the instrument flag PLUS per-line clean-line selection, and the
        # per-line part is enforced where it belongs — on the LINE, in measure_line —
        # not by refusing the whole band here.
        # RYA-806 adds the SECOND axis. `telluric_basis` says whether the BAND needs
        # correction; `telluric_applied` says whether THIS PRODUCT already got one, and
        # the second is not derivable from the first. When the context names the holding
        # we ask the question that consumes both; otherwise we fall back to the
        # instrument-only gate, which is all the older callers can answer.
        from pipeline.telluric_policy import gate as _telluric_gate
        from pipeline.telluric_policy import gate_holding as _telluric_gate_holding
        _inst = context.get("instrument") or context.get("instrument_id")
        _holding = context.get("holding") or context.get("holding_id")
        if _holding:
            # TelluricStateUnknown propagates: an unverified telluric state is its own
            # refusal and must not be flattened into the generic RuntimeError below.
            ok, basis = _telluric_gate_holding(_holding, _inst)
            if not ok:
                raise RuntimeError(basis)
            self._telluric_basis = basis
        elif _inst:
            ok, basis = _telluric_gate(_inst, bool(context.get("analysis_ready")))
            if not ok:
                raise RuntimeError(basis)
            self._telluric_basis = basis
        else:
            # No instrument in the context: fall back to the band flag rather than
            # assuming. An unnamed instrument cannot be looked up, and guessing is the
            # failure this ticket closes.
            if policy.telluric_required and not context.get("telluric_corrected"):
                raise RuntimeError(
                    f"the {policy.name} band requires telluric correction and the run "
                    f"context names no instrument, so the RYA-786 policy cannot resolve "
                    f"it. Pass `instrument` in the context.")
            self._telluric_basis = "band flag (no instrument in context)"

        import os
        from pipeline.abundances_derive import (
            _fit_synth_flux, _synth_flux_at_abund, _wingwide_window_nm,
            SYNTH_FIT_WSTEP_NM)
        # Imported, never re-typed: the edge trim below is only correct if it rebuilds
        # the grid the FIT actually uses (RYA-770).
        self._wstep = float(SYNTH_FIT_WSTEP_NM)
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

    def _quarantine_telluric(self, element, ion, c, instrument, why) -> LineMeasurement:
        """A telluric-band line, refused with its reason (RYA-786/711)."""
        dbg.trace_fallback('synth_line_rejected', f'{element} {ion} {c:.3f}: {why}',
                           severity='INFO', species=f'{element} {ion}', wavelength=c,
                           tag='QUARANTINED-TELLURIC')
        # Opts out of the EW-inversion path like every other row this handler builds:
        # it fits flux and inverts no equivalent width. The RYA-770 guard asserts every
        # construction site does so, and it caught this one when it did not.
        # (The guard counts literal occurrences, so do not name the flag in prose here.)
        lm = LineMeasurement(element=element, ion=ion, wavelength_air_A=c,
                             instrument=instrument, ew_mA=float("nan"),
                             ew_method="not measured — telluric band",
                             treatment="ENGINE-B", ew_inversion=False)
        lm.in_aggregate = False
        lm.excluded_reason = why
        return lm

    def measure_line(self, wav, flux, *, element, ion, wavelength_A, instrument,
                     policy: BandPolicy, pre_normalised: bool,
                     context: dict[str, Any]) -> LineMeasurement:
        c = float(wavelength_A)

        # RYA-786: the telluric decision is PER LINE for a reference atlas. A line inside
        # an enumerated O2/H2O band is quarantined as a stated physics exclusion (RYA-777),
        # never measured and never silently corrected. Lines outside them run.
        from pipeline.telluric_policy import exclusion as _telluric_exclusion
        _why = _telluric_exclusion(c, instrument)
        if _why:
            return self._quarantine_telluric(element, ion, c, instrument, _why)

        # RYA-847 — the constraint metrics for THIS line, filled in once the fit has run.
        # A dict rather than locals because `quarantine` is a closure defined before the
        # fit exists and must be able to emit metrics for post-fit rejections while
        # correctly emitting None for the pre-fit ones (SETUP / COVERAGE), where no chi2
        # surface was ever built.
        _fit_metrics: dict = {}

        def _mk_line(**kw) -> LineMeasurement:
            """THE single constructor for every line this handler emits.

            Both exits used to build a `LineMeasurement` independently, so a field added
            to one arrived on accepted lines and not on rejected ones -- and a rejected
            line with no chi2 is indistinguishable from one whose chi2 was never carried.
            One constructor makes that impossible (RYA-845).
            """
            return LineMeasurement(
                element=element, ion=ion, wavelength_air_A=c, instrument=instrument,
                # RYA-770/342 — set HERE, not at the call sites. Every line this handler
                # emits comes from a flux-space chi2 fit and inverts no equivalent width,
                # so the REW saturation ceiling never applies to any of them. Leaving it
                # to each caller meant the invariant held only by everyone remembering;
                # owning it here makes forgetting impossible.
                ew_inversion=False,
                sigma_A=_fit_metrics.get("sigma_A"),
                frac_rise_weaker=_fit_metrics.get("frac_rise_weaker"),
                edge_distance_dex=_fit_metrics.get("edge_distance_dex"),
                red_chi2=_fit_metrics.get("red_chi2"),
                **kw)

        def quarantine(reason: str, ew=float("nan")) -> LineMeasurement:
            # RYA-770: THE single choke point for every rejection this handler makes,
            # so one tracer call covers all of them and none can be added later that
            # escapes the timeline. The tag is the reason's leading token
            # (COVERAGE / CHI2-GATE / NOT-IN-SYNTH-LINELIST / ...), which is what makes
            # the rejections countable per class rather than just readable.
            dbg.trace_fallback(
                'synth_line_rejected', f'{element} {ion} {c:.3f}: {reason}',
                severity='WARN', species=f'{element} {ion}', wavelength=c,
                reject_class=str(reason).split(':', 1)[0].strip(),
                band=policy.name, instrument=str(instrument))
            lm = _mk_line(ew_mA=ew,
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
        # Probe on the CORE's own synthesis grid, not the observed one. Two reasons:
        # the sensitivity being tested is a property of the line, which the fine grid
        # resolves and the observed sampling may straddle; and the same two syntheses
        # then reveal the zeroed edge span for free (see EDGE_ZERO_FLUX above), so the
        # trim costs no extra Turbospectrum call.
        sw = np.arange(wave_base, wave_top + self._wstep * 0.5, self._wstep)
        try:
            lo = self._synth(sw, trial_A=a0 - 0.5, **kw)
            hi = self._synth(sw, trial_A=a0 + 0.5, **kw)
            sens = float(np.max(np.abs(hi - lo)))
        except Exception as e:
            return quarantine(f"SYNTH-FAILED probing sensitivity: {type(e).__name__}: {e}")

        # ── drop observed pixels that would be compared against a zeroed synthetic ──
        clean_lo, clean_hi, n_zero = clean_span(sw, lo, hi)
        if clean_lo is None:
            return quarantine(
                f"SYNTH-ALL-ZERO: every one of the {sw.size} synthesis pixels came back "
                f"at the ~0 sentinel, so there is no spectrum to fit against.")
        drop = edge_pixels_to_drop(obs_w_nm, wave_base, wave_top, clean_lo, clean_hi)
        n_dropped = int(drop.sum())
        if n_dropped:
            dbg.trace_fallback(
                'synth_edge_pixels_trimmed',
                f'{element} {ion} {c:.3f}: dropped {n_dropped} observed pixel(s) that '
                f'fall on iSpec-zeroed synthesis edges ({n_zero} zeroed of {sw.size})',
                severity='INFO', species=f'{element} {ion}', wavelength=c,
                n_obs_dropped=n_dropped, n_synth_zeroed=n_zero, n_synth=int(sw.size))
            obs_w_nm, obs_f = obs_w_nm[~drop], obs_f[~drop]
            m = (obs_w_nm >= wave_base) & (obs_w_nm <= wave_top)
            if m.sum() < 5:
                return quarantine(
                    f"COVERAGE: only {int(m.sum())} observed pixels survive after "
                    f"removing {n_dropped} that sit on zeroed synthesis edges")

        if sens < MIN_SENSITIVITY:
            return quarantine(
                f"NOT-IN-SYNTH-LINELIST: a +/-0.5 dex abundance swing moves the synthetic "
                f"flux by at most {sens:.5f}. This line is absent from the synthesis line "
                f"list, or carries a gf too weak to register, so no abundance can be "
                f"fitted from it. Coverage gap, not a measurement failure.")

        # ── delegate to the validated fit ─────────────────────────────────────
        # RYA-783: the NLTE path REFUSES a line whose levels are unidentified — bsyn
        # would set departure = 1 and return LTE under an NLTE label. That is per-LINE
        # coverage, not a run-level fault, so it is quarantined with its own reason
        # rather than aborting the whole band (RYA-429 no silent drops, RYA-711 never
        # cull). A band where EVERY line is unlabelled then yields an honest n=0 product.
        try:
            res = self._fit(
                obs_w_nm, obs_f, context["atmosphere"], context["teff"], context["logg"],
                context["feh"], context["vturb"], context["linelist"], context["isotopes"],
                context["solar_abund"], element, context["atom_code"],
                wave_base, wave_top, A_LO, A_HI,
                float(context["resolving_power"]), float(context["macroturbulence"]),
                float(context["vsini"]),
            # RYA-798. None keeps this the LTE call RYA-770 stabilised, byte for byte.
            # 'gerber' turns on the TS-native departure deck; `_fit_synth_flux` then
            # asserts the linelist actually carries NLTE labels for this element, because
            # iSpec would otherwise drop it into `nlte_ignored` and synthesise LTE without
            # raising (RYA-764).
                nlte_deck=context.get("nlte_deck"),
            # RYA-1040 — the <3D> route. `nlte_deck_key` NAMES the deck instead of letting
            # `_fit_synth_flux` infer it from the element (which resolves to the 1D deck
            # and can never reach `<El>@mean3D`); `atmosphere_layers_file` hands iSpec the
            # mul23 <3D> model as a FILE, because a five-column <3D> model cannot be
            # written in the MARCS form iSpec's writer produces. Both default to None via
            # `.get`, so every existing caller reaches the fit unchanged.
                nlte_deck_key=context.get("nlte_deck_key"),
                atmosphere_layers_file=context.get("atmosphere_layers_file"),
                tmp_dir=self._tmp_dir)
            # RYA-847 — carry what the fitter measured about whether A(X) was pinned.
            # Placed immediately after the fit and BEFORE any gate, so a rejected line
            # keeps the numbers that justify its rejection: RYA-711 requires a quarantined
            # line to be reported with its evidence, not merely dropped.
            _fit_metrics.update(
                sigma_A=_as_f(res.get("sigma_A")),
                frac_rise_weaker=_as_f(res.get("frac_rise_weaker")),
                edge_distance_dex=_as_f(res.get("edge_distance_dex")),
                red_chi2=_as_f(res.get("red_chi2")))
        except Exception as _e:
            from pipeline.gerber_nlte import GerberDeckError as _GDE
            if isinstance(_e, _GDE):
                return quarantine(f"NLTE-UNLABELLED: {str(_e)[:220]}")
            raise

        if res["status"] == "failed":
            return quarantine(f"FIT-FAILED: {res.get('reason', 'unknown')}")
        a_best = float(res["A_X"])

        # Synthetic EW of the TARGET LINE ONLY, for comparability with the profile fitter.
        # Integrated over a range WIDER than the fit window: the fit window is deliberately
        # tight so the line dominates chi2, but a DIFFERENCE of two syntheses is zero
        # wherever this element does not absorb, so a wide range adds the damping wings and
        # nothing else -- blends cancel. Integrating over the fit window instead truncated
        # the wings and drove the EW angle to a median ratio of 0.250 by construction.
        # RYA-878 — ONE definition, in `pipeline.synth_ew`. It was written out here and
        # the production synth-v2 path computed no synthetic EW at all, which is why the
        # control's ANGLE 1 could never be measured like-for-like: there was nothing on
        # the engine side to compare against. Both sides call this now.
        ew_hw = synth_ew.ew_half_width_A(hw_A)
        ew = synth_ew.synthetic_ew_mA(self._synth, centre_A=c, abundance=a_best,
                                      fit_half_width_A=hw_A, synth_kwargs=kw)

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

        # RYA-847 — THE CONSTRAINT GATE. Applied AFTER the chi2 merit gate and the edge
        # check so the three reasons stay distinguishable in the artifact instead of
        # collapsing into one "rejected": "the model does not describe the flux",
        # "the answer sits on a search bound", and "the data did not pin A(X)" are
        # different findings and a reader needs to know which one fired.
        #
        # ⚠️ THIS IS NOT A NO-OP, and the comment that used to say so was stale by the
        # time the sweep landed. `SYNTH_CONSTRAINT` is None — permanently, because the
        # nine-cell sweep found no transferable threshold — but the gate still applies
        # the NON-MINIMUM check unconditionally: `frac_rise <= 0` is a correctness test
        # with no number to choose, and it excludes lines here today (5333.768 in both
        # VIS decks, 3617.318 in the near-UV, 12638.703 in the NIR).
        _cv = _constraint_verdict(res.get("constraint") or {})
        if not _cv.ok:
            # State the fitted value, as CHI2-GATE and FIT-EDGE-PINNED above both do:
            # `quarantine` blanks `abundance` for every line it rejects, so the reason
            # string is the ONLY place a reader can see what the fit returned. It matters
            # most exactly here — the cut turns on the SHAPE of chi2 and never on the
            # value, so a reader must be able to check that a plausible-looking number
            # was cut on evidence (3617.318 fitted at a perfectly solar 7.477).
            return quarantine(f"{_cv.reason} Fit reported A={a_best:.3f}.", ew=ew)

        # RYA-770: the accepted lines are recorded too. A rejection ledger that only
        # lists rejections cannot tell you 3-of-12 from 3-of-3 — the accept count is
        # half the diagnosis.
        dbg.trace_decision('synth_line_accepted', f'{element} {ion} {c:.3f}',
                           detail=f'A={a_best:.3f}, red_chi2={res["red_chi2"]}, '
                                  f'n_pix={res["n_pix"]}, window +/-{hw_A:.3f} A',
                           species=f'{element} {ion}', wavelength=c, abundance=a_best,
                           red_chi2=res['red_chi2'], n_pix=res['n_pix'],
                           half_window_A=round(hw_A, 4), ew_hint_mA=ew_hint)
        return _mk_line(ew_mA=ew, abundance=a_best, ew_method=method)


register(SynthesisHandler())

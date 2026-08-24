"""Normalisation state determined at INTAKE, from the FLUX (RYA-1030).

THE AXIS THIS OWNS
------------------
`pre_normalised` on a HoldingSpec, and the ratified `PRE_NORMALISED_HOLDINGS` registry
(RYA-1026), record what a product is DECLARED to be. This module answers the orthogonal
question the declaration cannot: **is this product ACTUALLY normalised?** -- read off the
data itself, never asserted from a filename, a label, or a README.

It is the undeclared half of the double-normalise guard. RYA-1026 refuses to renormalise
a product DECLARED pre-normalised; this catches a product that IS normalised when nothing
said so, and the mis-ROUTED case where the declaration is right about the product and
wrong about the file.

WHY A DECLARATION IS NOT ENOUGH
-------------------------------
🔴 A DECLARED FLAG AND A MIS-ROUTED FILE AGREE WITH EACH OTHER PERFECTLY AND ARE BOTH
WRONG. There is no cross-check to be had between them, because they are the same claim
written twice. Two live instances:

  * KP2005 (RYA-929 -> RYA-933/1026). `0irrad.readme` does not list `irradrelwl.dat`, so
    intake took only `irradthu.dat` -- absolute irradiance, W/m2/nm -- and recorded the
    holding as shipping no continuum. The harness then fitted its own, tilting the band
    4% blue-to-red and biasing A(Fe I) low by 0.0218 +/- 0.0040 dex with a wavelength
    trend (r=+0.373). MEASURED here: `irradthu.dat` has a median rolling-P95 of 2.14 at
    5000-5100 A; the shipped residual `irradrelwl.dat` has 0.9883. This module separates
    them on sight.
  * `iag_fts_solar_atlas` (RYA-944, still open). Catalogued `telluric_basis=corrected`
    while the manifest routes to the telluric-RETAINING Reiners+2016 file. Same shape on
    the telluric axis rather than this one.

So the DATA gets a vote. This is the pattern `telluric_intake` established for RYA-806 --
determine the conditioning axis from the product, at intake, and record it per holding --
applied to the third axis. Three now: `telluric_applied` (RYA-806),
`observed_conditioning` (RYA-1006), `normalization_state` (this).

THE DETECTOR
------------
Scan the UPPER ENVELOPE: a rolling high percentile over running windows.

  * normalised / pseudo-residual -> envelope pinned near 1.0 AND FLAT.
  * un-normalised -> envelope far from 1.0, and/or carrying a wavelength trend
    (blaze / SED).

BOTH CONDITIONS ARE REQUIRED AND NEITHER IS SUFFICIENT, which is the whole reason there
are two:

  * KP2005 `irradthu` sits at 2.14 but its slope over 100 A is only -0.044 -- INSIDE the
    flatness bound. The LEVEL test is what catches it.
  * A blazed product can sit near unity over a narrow window. The SLOPE test is what
    catches that.

Both are tested in isolation so neither can quietly stop working while the other covers
for it.

THE RULE
--------
Cross-check DETECTED against DECLARED:

  * agree           -> proceed.
  * disagree        -> LOUD STOP (RYA-833). NEVER an auto-fix, and never an auto-skip.
  * no declaration  -> the product is UNDECLARED and MUST be scanned before any continuum
    stage runs.

The detector INFORMS the declared state; it does not silently override the science. A
genuinely raw spectrum with a coincidentally flat continuum and a normalised one with a
bad blaze both need a human, and auto-correcting either would replace a loud wrong answer
with a quiet one.

THE THIRD VALUE IS LOAD-BEARING
-------------------------------
`unknown` is a real answer, not a failure to try -- exactly as in `telluric_intake`. Too
few windows means the data could not speak. Defaulting it to `normalised` would apply
unity as a continuum and inflate every EW silently; defaulting it to `un-normalised`
would invite a second continuum onto a product that already has one. So it is reported,
never resolved by convention.

A NOTE ON SPELLING, so nobody "fixes" it: the module and the registry column use the
ticket's `normalization`, while the state values and the older API keep the `normalised`
spelling already carried by `pre_normalised` throughout the harness.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

NORMALISED = "normalised"
UN_NORMALISED = "un-normalised"
UNKNOWN = "unknown"
VALUES = (NORMALISED, UN_NORMALISED, UNKNOWN)

#: 🔴 BELOW THIS THERE IS NO TRUE CONTINUUM TO FIND, so this module does not pretend to.
#: Near-UV line blanketing depresses the 95th percentile of EVERY window, whatever the
#: product's normalisation. MEASURED: the KNOWN-NORMALISED KP1984 residual-flux atlas
#: FAILS ITS OWN TEST at every probe below 4200 A (0.86-0.94) and passes uniformly
#: (0.985-0.998) from 4500 A up.
#:
#: This is NOT a new concept and the number is NOT invented. It is the project's existing
#: CONTINUUM_LIMITED class -- "blue-edge no-true-continuum", RYA-451/460 -- already carried
#: per line in `problem_children.csv` for NH 3360, CN 3883, Sc II 4246, Co I 3845 and
#: Sr II 4077/4215 ("Blue HARPS edge; crowded -- large continuum uncertainty"). And 4500 A
#: is already `diagnostics_abundance.BLUE_EDGE`, drawn on the A(Fe I)-vs-wavelength plot
#: as "blue edge -- exclude". Two independent derivations, same number.
#:
#: ⚠️ It is deliberately ABOVE `PIPELINE['blue_edge_warn_A']` (3900 A), which answers a
#: different question: that flags a line as low-SNR, this asks where a CONTINUUM can be
#: located at all. Collapsing them would be the RYA-806 mistake of treating two axes as one.
CONTINUUM_LIMITED_BLUE_EDGE_A = 4500.0

#: DECLARED, not magic -- DERIVED FROM THE MEASURED GAP between the two populations
#: WITHIN THE VALIDATED DOMAIN (>= the blue edge, telluric bands excluded), not chosen as
#: a band around unity. Every number was measured through the harness's own reader.
#:
#: NORMALISED population -- KP1984 residual flux scanned 4500-9900 A in 200 A steps,
#: KP2005's shipped residual `irradrelwl.dat`, HARPS 5346/6885 A, RYA-794 CRIRES+ Y:
#:     0.968 - 1.0063  =>  worst distance from unity 0.032
#: UN-NORMALISED population -- KP2005 `irradthu.dat` (W/m2/nm) and KP1984's absolute
#: column, the latter a PAIRED control: the same spectrum at the same wavelengths in both
#: states, so nothing but the normalisation differs:
#:     1.6963 · 2.1398 · 2.2398 · 213.2875  =>  nearest distance 0.696
#:
#: Nothing lies between 0.032 and 0.696. 0.15 is the geometric midpoint of that gap:
#: 4.7x clear of the worst normalised case and 4.7x clear of the nearest un-normalised
#: one. A cut pressed against either wall is one that misclassifies the next product.
#:
#: 🔴 THE FIRST VERSION OF THIS CONSTANT WAS 0.05 AND IT WAS WRONG -- derived from clean
#: windows at 5000-5100 A only, and it did not survive contact with the full range. The
#: fix was not to widen it until the failures stopped: it was to find out WHY the
#: reference failed (blanketing, tellurics, fill values), exclude those regimes for stated
#: reasons, and only then re-derive the cut inside the domain that remained.
UNITY_TOLERANCE = 0.15
#: An envelope that RISES OR FALLS across the band carries a blaze or an SED whatever its
#: level. Measured: +0.010 to +0.023 for the normalised products, +3.42 for absolute flux.
ENVELOPE_SLOPE_MAX = 0.10
#: Wide enough that most windows contain some true continuum; see the 0.98 note above.
ENVELOPE_WINDOW_A = 2.0
#: Below this the statistic describes noise, so the answer is UNKNOWN -- never a default.
MIN_ENVELOPE_WINDOWS = 5
#: 🔴 TELLURIC BANDS ARE EXCLUDED FROM THE ENVELOPE. Saturated atmospheric absorption
#: pushes the 95th percentile down for a reason that has NOTHING to do with normalisation,
#: and on a telluric-retaining product it fakes the un-normalised signature: scanned across
#: its range, KP1984 -- known-normalised -- failed at 9300 and 9500 A, both inside the
#: registered H2O complex, and nowhere else between 4200 and 9900. The band list is NOT
#: re-enumerated here; it is `telluric_policy.TELLURIC_BANDS`, the authoritative set
#: (RYA-786). A second copy would be the RYA-845 defect shape.
EXCLUDE_TELLURIC_BANDS = True

#: A window must hold at least this many VALID pixels to contribute a percentile.
MIN_PIXELS_PER_WINDOW = 20

#: 🔴 A FILL VALUE IS NOT FLUX. Non-finite and NON-POSITIVE pixels are dropped before the
#: envelope is computed. Found on real data, not imagined: `solar_harps` at 5321-5337 A
#: returns EXACTLY 0.000 across eight consecutive 2 A windows -- a dead region in the
#: product -- and counting those zeros as flux made the envelope ramp 0.00 -> 0.97 across
#: the band, which this module then reported as a BLAZE. The verdict "something is wrong
#: here" was right and the DIAGNOSIS was wrong, which is worse than either being wrong
#: alone: it sends the reader after a continuum defect that does not exist. Zeros are also
#: how a saturated core, a masked telluric band, and an unfilled array all present
#: themselves, so none of them may be averaged in (RYA-844's rule, on pixels).


@dataclass
class NormalisationEvidence:
    """What the FLUX actually said. `value` is the verdict; the rest is the citation."""

    value: str
    median_p95: float | None = None
    slope_across_band: float | None = None
    n_windows: int = 0
    n_dropped: int = 0
    n_telluric_windows: int = 0
    n_blue_windows: int = 0
    max_flux: float | None = None
    span_A: tuple[float, float] | None = None
    reason: str = ""

    @property
    def at_unity(self) -> bool | None:
        if self.median_p95 is None:
            return None
        return abs(self.median_p95 - 1.0) <= UNITY_TOLERANCE

    @property
    def flat(self) -> bool | None:
        if self.slope_across_band is None:
            return None
        return abs(self.slope_across_band) <= ENVELOPE_SLOPE_MAX

    def citation(self) -> str:
        """One line fit to paste into the holdings registry `notes`."""
        bits = [f"normalization_state={self.value}"]
        if self.span_A:
            bits.append(f"span {self.span_A[0]:.2f}-{self.span_A[1]:.2f} A")
        if self.median_p95 is not None:
            bits.append(f"median rolling-P95={self.median_p95:.4f} "
                        f"(tol |x-1|<={UNITY_TOLERANCE})")
        if self.slope_across_band is not None:
            bits.append(f"envelope slope={self.slope_across_band:+.4f} "
                        f"(max {ENVELOPE_SLOPE_MAX})")
        if self.max_flux is not None:
            bits.append(f"max flux={self.max_flux:.4f}")
        bits.append(f"{self.n_windows} x {ENVELOPE_WINDOW_A:.0f} A windows")
        if self.n_blue_windows:
            bits.append(f"{self.n_blue_windows} window(s) skipped below the "
                        f"{CONTINUUM_LIMITED_BLUE_EDGE_A:.0f} A blue edge "
                        f"(no true continuum there, RYA-451/460)")
        if self.n_telluric_windows:
            bits.append(f"{self.n_telluric_windows} window(s) skipped as registered "
                        f"telluric bands (absorption is not a normalisation signal)")
        if self.n_dropped:
            bits.append(f"{self.n_dropped} px dropped (non-finite or non-positive -- "
                        f"a fill value is not flux)")
        if self.reason:
            bits.append(self.reason)
        return "; ".join(bits)


class NormalisationStateMismatch(RuntimeError):
    """The flux says one thing and the declared flag says another."""


def _telluric_bands() -> tuple[tuple[float, float], ...]:
    """(lo, hi) of every registered telluric complex. CALLED, never re-enumerated."""
    from pipeline.telluric_policy import TELLURIC_BANDS
    return tuple((float(b[0]), float(b[1])) for b in TELLURIC_BANDS)


def detect(wave, flux) -> NormalisationEvidence:
    """Classify a spectrum from its UPPER ENVELOPE. See the module docstring."""
    w = np.asarray(wave, float)
    f = np.asarray(flux, float)
    n_raw = w.size
    good = np.isfinite(w) & np.isfinite(f) & (f > 0.0)
    w, f = w[good], f[good]
    dropped = n_raw - w.size
    if w.size < 2:
        return NormalisationEvidence(
            UNKNOWN, reason=(f"fewer than 2 usable pixels ({dropped} of {n_raw} were "
                             f"non-finite or non-positive)"))

    order = np.argsort(w)
    w, f = w[order], f[order]
    if float(w[-1]) < CONTINUUM_LIMITED_BLUE_EDGE_A:
        return NormalisationEvidence(
            UNKNOWN, span_A=(float(w[0]), float(w[-1])), max_flux=float(f.max()),
            n_dropped=int(dropped),
            reason=(f"entirely below the CONTINUUM_LIMITED blue edge "
                    f"({CONTINUUM_LIMITED_BLUE_EDGE_A:.0f} A): near-UV line blanketing "
                    f"leaves no true continuum in any window, so the envelope says "
                    f"nothing about normalisation here -- the known-normalised KP1984 "
                    f"atlas fails this test on itself below 4200 A. This is the existing "
                    f"blue-edge no-true-continuum class (RYA-451/460), not a new limit. "
                    f"UNKNOWN is the honest answer; scan a redder window."))

    skip = _telluric_bands() if EXCLUDE_TELLURIC_BANDS else ()
    p95, n_telluric = [], 0
    edges = np.arange(w[0], w[-1], ENVELOPE_WINDOW_A)
    n_blue = 0
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi < CONTINUUM_LIMITED_BLUE_EDGE_A:
            n_blue += 1
            continue
        if any(not (hi < blo or lo > bhi) for blo, bhi in skip):
            n_telluric += 1
            continue
        m = (w >= lo) & (w < hi)
        if m.sum() > MIN_PIXELS_PER_WINDOW:
            p95.append(float(np.percentile(f[m], 95)))
    p95 = np.asarray(p95)

    if p95.size < MIN_ENVELOPE_WINDOWS:
        return NormalisationEvidence(
            UNKNOWN, n_windows=int(p95.size), span_A=(float(w[0]), float(w[-1])),
            max_flux=float(f.max()), n_dropped=int(dropped),
            n_telluric_windows=int(n_telluric),
            reason=(f"only {p95.size} usable {ENVELOPE_WINDOW_A:.0f} A windows "
                    f"(need {MIN_ENVELOPE_WINDOWS}) -- too few to describe an envelope, "
                    f"so the answer is UNKNOWN, never a default"))

    ev = NormalisationEvidence(
        UNKNOWN,
        n_dropped=int(dropped), n_telluric_windows=int(n_telluric),
        n_blue_windows=int(n_blue),
        median_p95=float(np.median(p95)),
        slope_across_band=float(np.polyfit(np.arange(p95.size), p95, 1)[0] * p95.size),
        n_windows=int(p95.size), max_flux=float(f.max()),
        span_A=(float(w[0]), float(w[-1])))

    if ev.at_unity and ev.flat:
        ev.value = NORMALISED
    elif not ev.at_unity:
        ev.value = UN_NORMALISED
        ev.reason = (f"envelope sits at {ev.median_p95:.4f}, "
                     f"{abs(ev.median_p95 - 1.0):.4f} from unity")
    else:
        ev.reason = (f"envelope is at unity but TILTED "
                     f"({ev.slope_across_band:+.4f} across the band) -- a blaze or SED "
                     f"survives, so this is not cleanly normalised and is not waved "
                     f"through")
    return ev


def cross_check(holding_id: str, wave, flux, *, declared: bool | None,
                where: str = "") -> NormalisationEvidence:
    """Cross-check DETECTED against DECLARED. LOUD on disagreement, never an auto-fix.

    `declared=None` means the product carries NO normalisation declaration. That is not
    an excuse to proceed: an undeclared product MUST be scanned before any continuum
    stage runs, so this raises with the detected state attached, which is the thing
    needed to declare it.

    `unknown` does NOT raise -- the data could not speak, and turning that into a hard
    failure would make every short window a blocker. It is returned to be recorded.
    """
    ev = detect(wave, flux)
    site = f" [{where}]" if where else ""

    if declared is None:
        raise NormalisationStateMismatch(
            f"{holding_id}: UNDECLARED normalisation state{site}. Nothing records whether "
            f"this product ships a continuum, and a continuum stage must not run until "
            f"something does. The scan says {ev.value.upper()} -- {ev.citation()}. "
            f"Record it on the holding (RYA-1030) rather than letting the next stage "
            f"assume; absence of a claim is not a claim (RYA-833).")

    if ev.value == UNKNOWN:
        return ev

    if (ev.value == NORMALISED) == declared:
        return ev

    raise NormalisationStateMismatch(
        f"{holding_id}: THE FLUX AND THE FLAG DISAGREE{site}. Declared "
        f"pre_normalised={declared}, but the data says {ev.value.upper()} -- "
        f"{ev.citation()}. This is a LOUD STOP, not an auto-fix: the detector informs "
        f"the declared state, it must not silently override the science, auto-normalise, "
        f"or auto-skip. A declared flag and a MIS-ROUTED file agree with each other "
        f"perfectly and are both wrong -- CHECK WHICH FILE THE READER OPENS before you "
        f"touch the flag. That is what KP2005 turned out to be (RYA-929 -> RYA-933/1026): "
        f"the flag was right about the product and the route pointed at "
        f"`irradthu.dat`, absolute irradiance.")

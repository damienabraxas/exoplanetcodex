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

#: DECLARED, not magic. DERIVED from four measured controls rather than chosen:
#:
#:   | product                              | median rolling-P95 | slope across band |
#:   | KP1984 col1 (residual flux)          |             0.9800 |            +0.023 |
#:   | KP2005 `irradrelwl.dat` (residual)   |             0.9883 |            +0.017 |
#:   | KP2005 `irradthu.dat`  (irradiance)  |        1.72 - 2.26 |      strong (SED) |
#:   | KP1984 col2 (absolute)               |           213.2875 |            +3.422 |
#:
#: KP1984's two columns are a PAIRED control: the same spectrum at the same wavelengths
#: in both states, so nothing but the normalisation differs between them.
#:
#: The nearest NORMALISED case sits 0.020 from unity and the nearest UN-NORMALISED case
#: 0.72 away, so this carries a ~14x margin on both sides. The margin is the point: a cut
#: pressed against the nearest case is one that misclassifies the next product, and
#: chance cases pile up at the wall. KP1984 col1 lands at 0.98 rather than 1.00 because a
#: narrow window in a line-rich region may contain no true continuum at all -- which is
#: also why the window is wide and the tolerance is not tighter.
UNITY_TOLERANCE = 0.05
#: An envelope that RISES OR FALLS across the band carries a blaze or an SED whatever its
#: level. Measured: +0.010 to +0.023 for the normalised products, +3.42 for absolute flux.
ENVELOPE_SLOPE_MAX = 0.10
#: Wide enough that most windows contain some true continuum; see the 0.98 note above.
ENVELOPE_WINDOW_A = 2.0
#: Below this the statistic describes noise, so the answer is UNKNOWN -- never a default.
MIN_ENVELOPE_WINDOWS = 5
#: A window must hold at least this many pixels to contribute a percentile.
MIN_PIXELS_PER_WINDOW = 20


@dataclass
class NormalisationEvidence:
    """What the FLUX actually said. `value` is the verdict; the rest is the citation."""

    value: str
    median_p95: float | None = None
    slope_across_band: float | None = None
    n_windows: int = 0
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
        if self.reason:
            bits.append(self.reason)
        return "; ".join(bits)


class NormalisationStateMismatch(RuntimeError):
    """The flux says one thing and the declared flag says another."""


def detect(wave, flux) -> NormalisationEvidence:
    """Classify a spectrum from its UPPER ENVELOPE. See the module docstring."""
    w = np.asarray(wave, float)
    f = np.asarray(flux, float)
    good = np.isfinite(w) & np.isfinite(f)
    w, f = w[good], f[good]
    if w.size < 2:
        return NormalisationEvidence(UNKNOWN, reason="fewer than 2 finite pixels")

    order = np.argsort(w)
    w, f = w[order], f[order]
    p95, edges = [], np.arange(w[0], w[-1], ENVELOPE_WINDOW_A)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (w >= lo) & (w < hi)
        if m.sum() > MIN_PIXELS_PER_WINDOW:
            p95.append(float(np.percentile(f[m], 95)))
    p95 = np.asarray(p95)

    if p95.size < MIN_ENVELOPE_WINDOWS:
        return NormalisationEvidence(
            UNKNOWN, n_windows=int(p95.size), span_A=(float(w[0]), float(w[-1])),
            max_flux=float(f.max()),
            reason=(f"only {p95.size} usable {ENVELOPE_WINDOW_A:.0f} A windows "
                    f"(need {MIN_ENVELOPE_WINDOWS}) -- too few to describe an envelope, "
                    f"so the answer is UNKNOWN, never a default"))

    ev = NormalisationEvidence(
        UNKNOWN,
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

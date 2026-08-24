"""RYA-1026: a pre-normalised product may never be re-normalised. Enforced, not recalled.

THE FAILURE THIS EXISTS TO STOP. Fitting or applying a continuum to a product that
already ships one TILTS the spectrum and corrupts every measurement taken from it. It is
a SILENT corruption: the fit succeeds, the EW comes out, and nothing in the output says
the continuum was invented.

IT HAS NOW BITTEN TWICE, WHICH IS WHY THIS IS A GUARD AND NOT A NOTE:
  * RYA-940 -- KP 1984. A free polynomial followed the saturated telluric band down and
    took the continuum with it.
  * RYA-929/1026 -- KP 2005. `pre_normalised` was set False from a plausible-sounding
    reading of the file header, and the harness re-normalised the product for months
    with nothing anywhere objecting. It forced the VIS re-run.

RATIFIED (Ryan, RYA-1026): **DO NOT NORMALIZE ANY KITT PEAK ATLAS.** Both KP products
ship their own continuum -- KP1984 column 2 is pseudo-residual flux where unity IS the
continuum, KP2005 (`irradthu`) is absolute irradiance with its own continuum baked in.
The deliverable was widened from "the 2005 file" to THE WHOLE KP CLASS, and to any
pre-normalised reference arm (IAG, Delbouille are the same shape). The ONLY thing done to
a KP atlas on the way in is TELLURIC CORRECTION (KP1984 -> the RYA-940 sibling). Never a
continuum refit.

WHAT COUNTS AS RE-NORMALISING: fitting a continuum of our own, applying one, or re-pinning
a locally-normalised arm to exactly 1.0. Reading the product, measuring its continuum for
a REPORT, or comparing ours against its shipped one (RYA-911 `reference_continuum`) are
all fine -- those are readings. The line is whether the number we divide by is ours or
the product's.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: Products that arrive with their continuum already established. Keyed by holding_id so
#: it is per-PRODUCT, never per-instrument: one instrument can serve both a normalised
#: and an un-normalised holding (RYA-904 -- crires_plus does exactly that), and keying by
#: instrument re-creates the collapse that defect was about.
#:
#: This must stay in step with the `pre_normalised` flags on `measure_band_ew`'s
#: HoldingSpecs; `assert_not_renormalising` treats disagreement as the bug, and
#: `tests/test_rya1026_product_policy.py` asserts the two sets match exactly.
PRE_NORMALISED_HOLDINGS: frozenset[str] = frozenset({
    # THE KITT PEAK CLASS -- ratified whole (RYA-1026 comment, Ryan 2026-08-24).
    "solar_kpno",                        # KP1984: residual flux, unity IS the continuum
    "solar_kpno_molecfit_corrected",     # RYA-940 corrected sibling, same conventions
    "solar_kpno_kurucz2005_corrected",   # KP2005 irradthu -- REVERSES RYA-929
    # Other pre-normalised reference arms.
    "solar_harps",                       # ships its own fitted continuum (RYA-911)
    "solar_harps_molecfit_corrected",    # RYA-931 corrected sibling, same normaliser
    "solar_iag",                         # Baker+2020 telluric-free, normalised
    "solar_crires_plus_y_rya794",        # RYA-794 science-ready Y arm, normalised
    "solar_delbouille_liege",            # RYA-944 disk-centre intensity, LOCAL normalisation
})

#: Locally-normalised arms: normalised per-window, so there is NO absolute flux and no
#: broad-band continuum shape -- and the maximum over the whole product is BELOW 1.0
#: (Delbouille: 0.9959 over 2.38 M points, RYA-944). Re-pinning these to exactly 1.0 is a
#: re-normalisation wearing a tidy-up's clothes: it rescales every point by a number we
#: chose, which is the same corruption as fitting a continuum, just harder to notice.
LOCALLY_NORMALISED: frozenset[str] = frozenset({
    "solar_delbouille_liege",
})


class RenormalisationError(RuntimeError):
    """A pre-normalised product was about to have a second continuum placed on it."""


def assert_not_renormalising(holding_id: str, *, pre_normalised: bool,
                             fitting_continuum: bool = False,
                             applying_continuum: bool = False,
                             pinning_unity: bool = False,
                             where: str = "") -> None:
    """Raise if `holding_id` ships normalised and the caller is about to touch its continuum.

    Two independent signals are checked against each other, deliberately:

      * the registry (`PRE_NORMALISED_HOLDINGS`) -- what the ratified decision says
      * the holding spec's own `pre_normalised` flag -- what the code is configured with

    A disagreement between them is itself the bug. RYA-929's flag said False while the
    product shipped normalised, and nothing anywhere objected. Checking one signal would
    have re-blessed whichever was written down; checking both makes drift LOUD.

    `fitting_continuum` / `applying_continuum` / `pinning_unity` are separate arguments
    rather than one boolean because they are three different mistakes with three
    different fixes, and the third one does not look like a mistake at all.
    """
    registered = holding_id in PRE_NORMALISED_HOLDINGS
    site = f" [{where}]" if where else ""

    if registered != pre_normalised:
        raise RenormalisationError(
            f"{holding_id}: DRIFT{site}. The RYA-1026 registry says "
            f"pre_normalised={registered} but the holding spec is configured "
            f"{pre_normalised}. These must agree. Do not 'fix' this by editing whichever "
            f"is convenient -- decide which is true for the PRODUCT and change that one, "
            f"because this exact disagreement is how RYA-929 shipped a product that was "
            f"re-normalised for months without anything objecting.")

    if not registered:
        return

    if fitting_continuum or applying_continuum:
        verb = "FIT" if fitting_continuum else "APPLY"
        raise RenormalisationError(
            f"{holding_id}: REFUSING to {verb} a continuum{site}. This product SHIPS "
            f"NORMALISED (RYA-1026 -- for the Kitt Peak class this is ratified WHOLE, and "
            f"for KP2005 it REVERSES RYA-929). A second continuum TILTS the spectrum and "
            f"corrupts every measurement taken from it -- the RYA-940 failure mode, and "
            f"it is SILENT: the fit succeeds and the EW looks fine. Use the continuum the "
            f"product arrives with. The only thing we do to a KP atlas on the way in is "
            f"telluric correction. If you believe this product genuinely ships no "
            f"continuum, that is a decision to re-ratify on the ticket, not a flag to "
            f"flip here.")

    if pinning_unity and holding_id in LOCALLY_NORMALISED:
        raise RenormalisationError(
            f"{holding_id}: REFUSING to pin the continuum to 1.0{site}. This arm is "
            f"LOCALLY normalised (RYA-944): there is no absolute flux and no broad-band "
            f"continuum shape, and the maximum over the whole product is BELOW unity "
            f"(0.9959 over 2.38 M points for Delbouille). Re-pinning it to exactly 1.0 "
            f"rescales every point by a number we chose -- the same corruption as fitting "
            f"a continuum, just harder to notice.")


def is_pre_normalised(holding_id: str) -> bool:
    """Registry lookup. Unknown holdings return False -- an unregistered product is NOT
    assumed normalised, because assuming it would apply unity as a continuum and inflate
    every EW silently. Absence of a claim is not a claim (RYA-833)."""
    return holding_id in PRE_NORMALISED_HOLDINGS


# ── the UNDECLARED half: is this product ACTUALLY normalised? ────────────────
#
# The registry above says what a product is DECLARED to be. It cannot catch a product
# that is mis-ROUTED -- a normalised holding whose reader opens the raw file. That is not
# hypothetical: `iag_fts_solar_atlas` is catalogued `corrected` while the manifest routes
# to the telluric-RETAINING Reiners file (RYA-944), and RYA-1026 found the same shape on
# KP2005 -- `solar_kpno_kurucz2005_corrected` is ratified pre-normalised, and its reader
# opens `irradthu.dat.txt`, ABSOLUTE IRRADIANCE in W/m**2/nm.
#
# A declared flag and a mis-routed file agree with each other perfectly and are both
# wrong. So the DATA gets a vote. RYA-1030 generalises this into `normalization_intake`
# (all spectra, at intake, backfilled onto the registry); this is the narrow version
# wired at the one site RYA-1026 itself re-flagged.

NORMALISED = "normalised"
UN_NORMALISED = "un-normalised"
UNKNOWN = "unknown"

#: DECLARED, not magic (RYA-1030). Derived from four measured controls, not chosen:
#:
#:   | product                              | median rolling-P95 | slope across band |
#:   | KP1984 col1 (residual flux)          |             0.9800 |            +0.023 |
#:   | KP2005 staged normalised TSV         |             0.9997 |            +0.010 |
#:   | KP2005 `irradthu` as served          |        1.72 - 2.26 |     strong (SED)  |
#:   | KP1984 col2 (absolute)               |           213.2875 |            +3.422 |
#:
#: The nearest NORMALISED case sits 0.020 from unity and the nearest UN-NORMALISED case
#: 0.72 away, so 0.05 has a ~14x margin on both sides. That margin is the point: a cut
#: pressed against the nearest case is a cut that will misclassify the next product
#: (RYA-a-borrowed-threshold-is-not-a-control). KP1984 col1 lands at 0.98 rather than
#: 1.00 because a narrow window in a line-rich region may hold no true continuum at all
#: -- which is why the window is wide and the tolerance is not tighter.
UNITY_TOLERANCE = 0.05
#: An envelope that RISES OR FALLS across the band carries a blaze or an SED, whatever
#: its level. Measured: +0.010 / +0.023 for the normalised pair, +3.42 for absolute flux.
ENVELOPE_SLOPE_MAX = 0.10
#: Wide enough that most windows contain some true continuum; see the 0.98 note above.
ENVELOPE_WINDOW_A = 2.0
#: Fewer than this and the statistic is describing noise, so the answer is UNKNOWN --
#: never a default (RYA-833, and the same reason `telluric_intake` keeps a third value).
MIN_ENVELOPE_WINDOWS = 5


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

    def citation(self) -> str:
        """One line fit to paste into a registry `notes` column."""
        bits = [f"normalisation_state={self.value}"]
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


def detect_normalisation_state(wave, flux) -> NormalisationEvidence:
    """Classify a spectrum from its UPPER ENVELOPE: rolling P95 over running windows.

    Normalised products pin the envelope at ~1.0 AND hold it flat. Both conditions are
    required, and the second is what a level test alone would miss: a raw product can
    sit near unity by coincidence of units over a narrow window, but it carries its
    blaze or SED as a trend. Conversely a normalised product with a residual tilt is
    reported UNKNOWN rather than waved through.
    """
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
        if m.sum() > 20:
            p95.append(float(np.percentile(f[m], 95)))
    p95 = np.asarray(p95)

    if p95.size < MIN_ENVELOPE_WINDOWS:
        return NormalisationEvidence(
            UNKNOWN, n_windows=int(p95.size), span_A=(float(w[0]), float(w[-1])),
            max_flux=float(f.max()),
            reason=(f"only {p95.size} usable {ENVELOPE_WINDOW_A:.0f} A windows "
                    f"(need {MIN_ENVELOPE_WINDOWS}) -- too few to describe an envelope, "
                    f"so the answer is UNKNOWN, never a default"))

    median_p95 = float(np.median(p95))
    slope = float(np.polyfit(np.arange(p95.size), p95, 1)[0] * p95.size)
    ev = NormalisationEvidence(
        UNKNOWN, median_p95=median_p95, slope_across_band=slope,
        n_windows=int(p95.size), max_flux=float(f.max()),
        span_A=(float(w[0]), float(w[-1])))

    at_unity = abs(median_p95 - 1.0) <= UNITY_TOLERANCE
    flat = abs(slope) <= ENVELOPE_SLOPE_MAX

    if at_unity and flat:
        ev.value = NORMALISED
    elif not at_unity:
        ev.value = UN_NORMALISED
        ev.reason = (f"envelope sits at {median_p95:.4f}, "
                     f"{abs(median_p95 - 1.0):.4f} from unity")
    else:
        ev.reason = (f"envelope is at unity but TILTED ({slope:+.4f} across the band) -- "
                     f"a blaze or SED survives, so this is not cleanly normalised and is "
                     f"not waved through")
    return ev


def assert_data_matches_declaration(holding_id: str, wave, flux, *,
                                    declared: bool, where: str = "") -> NormalisationEvidence:
    """Cross-check the DETECTED state against the DECLARED flag. LOUD on disagreement.

    A LOUD STOP, never an auto-fix (RYA-1030): a genuinely raw spectrum with a
    coincidentally flat continuum and a normalised one with a bad blaze both need a
    human. The detector INFORMS the declared state; it must not silently override the
    science, auto-normalise, or auto-skip.

    `unknown` does NOT raise -- it is a real answer meaning the data could not speak, and
    treating it as a mismatch would turn "too few windows" into a hard failure. It is
    returned for the caller to record.
    """
    ev = detect_normalisation_state(wave, flux)
    site = f" [{where}]" if where else ""
    if ev.value == UNKNOWN:
        return ev

    detected_normalised = ev.value == NORMALISED
    if detected_normalised == declared:
        return ev

    raise NormalisationStateMismatch(
        f"{holding_id}: THE FLUX AND THE FLAG DISAGREE{site}. Declared "
        f"pre_normalised={declared}, but the data says {ev.value.upper()} -- "
        f"{ev.citation()}. This is a LOUD STOP, not an auto-fix: the detector informs "
        f"the declared state, it must not silently override the science. A declared flag "
        f"and a MIS-ROUTED file agree with each other perfectly and are both wrong, which "
        f"is how `iag_fts_solar_atlas` came to be catalogued `corrected` while pointing at "
        f"the raw Reiners file (RYA-944). Check WHICH FILE the reader opens before you "
        f"touch the flag. Generalised at intake by RYA-1030.")

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
    "solar_iag_reiners2016",             # the blue arm, normalised at source
    "solar_crires_plus_y_rya794",        # RYA-794 science-ready Y arm, normalised
    # RYA-1054: the SAME source and normaliser as the row above, conditioned over the
    # arm's full measured extent (9800-10796 A) instead of Elgueta's adopted window. It
    # is a SIBLING product, so it inherits the normalisation verdict rather than being
    # re-judged -- but it is listed EXPLICITLY, because this set is keyed per PRODUCT and
    # inheriting by instrument is the exact collapse RYA-904 was about.
    "solar_crires_plus_y_wide_rya1054",
    # RYA-1094: the H arm (15007-17494 A), built from sp/Sun_H_rv.dat by the SAME
    # normaliser as the two Y products above (`normalize_vesta_ir.py --arm H`). Only the
    # arm differs. Listed EXPLICITLY for the same reason RYA-1054 is: this set is keyed
    # per PRODUCT, and letting a new arm inherit its sibling's verdict by instrument is
    # the RYA-904 collapse. Continuum after the residual-slope helper: median 0.996,
    # p95 1.013.
    "solar_crires_plus_h_rya1094",
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


# ── the UNDECLARED half lives in `normalization_intake` (RYA-1030) ───────────
#
# The detector that reads normalisation off the FLUX was born here, under RYA-1026, to
# catch the one mis-routed product this module's registry names. RYA-1030 generalises it
# to all spectra at intake, so it MOVED rather than being copied: two implementations of
# one measurement is the RYA-845 defect shape, and the copy that drifts is always the one
# nobody is looking at. Re-exported so RYA-1026's callers keep working.

from pipeline.normalization_intake import (  # noqa: E402,F401
    NORMALISED, UN_NORMALISED, UNKNOWN, ENVELOPE_SLOPE_MAX, ENVELOPE_WINDOW_A,
    MIN_ENVELOPE_WINDOWS, UNITY_TOLERANCE, NormalisationEvidence,
    NormalisationStateMismatch, detect as detect_normalisation_state)


def assert_data_matches_declaration(holding_id: str, wave, flux, *,
                                    declared: bool, where: str = "") -> NormalisationEvidence:
    """RYA-1026's name for `normalization_intake.cross_check`. See that module."""
    from pipeline.normalization_intake import cross_check
    return cross_check(holding_id, wave, flux, declared=declared, where=where)

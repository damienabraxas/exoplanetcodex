"""Instrument-neutral telluric route resolution (RYA-927).

This is deliberately a routing layer, not a correction engine.  It prevents each
loader from inventing its own interpretation of the catalog's two telluric axes.
Correction-required products are sent to the observation-aware molecfit/GDAS path;
reference atlases with terrestrial absorption use explicit clean-line selection;
already-corrected and above-atmosphere products pass through; unresolved modes stop
for audit rather than silently choosing a method.
"""
from __future__ import annotations

from dataclasses import dataclass

from pipeline.telluric_policy import (
    TelluricStateUnknown, applied_state, basis, gate_holding, requires_correction,
)


@dataclass(frozen=True)
class TelluricRoute:
    instrument: str
    basis: str
    engine: str
    correction_required: bool
    analysis_ready_required: bool
    reason: str


def route_for(instrument: str) -> TelluricRoute:
    """Resolve the sole allowed telluric route for a catalog instrument.

    ``audit_required`` is an intentional stop for ``unspecified`` and
    ``mode_dependent`` entries.  It is not a correction result and cannot be used to
    make data science-ready.
    """
    b = basis(instrument)
    required = requires_correction(instrument)
    if instrument == "harps":
        return TelluricRoute(
            instrument, b, "molecfit_gdas+clean_line_selection", True, False,
            "Hybrid optical route: clean visible lines remain eligible by line selection; "
            "windows intersecting a registered telluric band are conditioned with "
            "observation-night molecfit/GDAS.",
        )
    if b == "correction_required" or required:
        return TelluricRoute(
            instrument, b, "molecfit_gdas", True, True,
            "Atmospheric transmission must be modeled from the observation-night profile "
            "before affected-line measurements.",
        )
    if b == "corrected":
        return TelluricRoute(
            instrument, b, "corrected_product", False, False,
            "The held product declares a validated telluric correction; do not divide it twice.",
        )
    if b == "not_applicable":
        return TelluricRoute(
            instrument, b, "none", False, False,
            "The data are above the atmosphere; terrestrial telluric correction is not applicable.",
        )
    if b == "line_selection":
        return TelluricRoute(
            instrument, b, "clean_line_selection", False, False,
            "Tellurics remain in the reference spectrum; exclude affected lines using the "
            "authoritative band list and retain clean lines.",
        )
    return TelluricRoute(
        instrument, b, "audit_required", False, True,
        "Telluric basis is unresolved for this instrument/mode; determine it from the "
        "product and recipe before measurement.",
    )


def route_for_holding(holding_id: str, instrument: str | None = None) -> TelluricRoute:
    """Resolve the route for an actual product, including its applied state.

    The instrument route describes capability.  The holding route is what loaders must
    consume: an already-corrected CRIRES+/NIRPS/other product is a passthrough, while a
    raw product on the same instrument is sent to the correction engine.  Unknown state
    remains an audit stop and is never inferred from the instrument.
    """
    from pipeline.telluric_policy import HOLDINGS
    import pandas as pd

    rows = pd.read_csv(HOLDINGS)
    hit = rows[rows.holding_id.astype(str) == str(holding_id)]
    if not len(hit):
        raise KeyError(f"holding {holding_id!r} is not registered")
    inst = instrument or str(hit.iloc[0].instrument_id)
    state = applied_state(holding_id)
    if state == "applied":
        return TelluricRoute(
            inst, basis(inst), "corrected_product", False, False,
            f"{holding_id}: telluric correction is already present in the held product; "
            "serve its corrected column/product and do not correct twice.",
        )
    if state == "unknown":
        return TelluricRoute(
            inst, basis(inst), "audit_required", False, True,
            f"{holding_id}: telluric_applied is unknown; resolve from product headers and "
            "recipe provenance before measurement.",
        )
    # Exercise the existing holding gate as a consistency check.  A false result is the
    # expected route for an uncorrected correction-required product.
    try:
        gate_holding(holding_id, inst)
    except TelluricStateUnknown:
        return TelluricRoute(inst, basis(inst), "audit_required", False, True,
                             f"{holding_id}: holding gate refuses unresolved telluric state.")
    return route_for(inst)


def correction_needed_at(wave_A: float, instrument: str) -> bool:
    """Whether a route needs correction at a particular wavelength.

    This deliberately returns ``False`` for clean visible lines on a correction-capable
    instrument: the correction framework is still required for the product, while the
    line-level disposition records that no atmospheric feature overlaps this line.
    """
    from pipeline.telluric_policy import TELLURIC_BANDS

    route = route_for(instrument)
    if "molecfit_gdas" not in route.engine:
        return False
    return any(lo <= float(wave_A) <= hi for lo, hi, _ in TELLURIC_BANDS)

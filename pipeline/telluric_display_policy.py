"""RYA-1026: displayed science is telluric-CORRECTED input only. Enforced at RENDER.

RATIFIED (Ryan, 2026-08-23): we do not display non-telluric data any more. Most current
Fe products were built on `kpno_solar_atlas` -- the 1984 telluric-UNCORRECTED atlas --
and that was a mistake. Those products rebuild on the corrected siblings
(`solar_kpno_molecfit_corrected`, RYA-940; `solar_harps_molecfit_corrected`, RYA-931).

WHY A GATE AND NOT A CONVENTION. An uncorrected product does not look wrong. It renders
a number with a bar, and the telluric contamination shows up as a slightly deep line or
a slightly odd scatter -- indistinguishable from physics unless you already know. The
whole class of defect this ticket exists to close shares that shape: WRONG INFORMATION
THAT DOES NOT ANNOUNCE ITSELF. So the rule is enforced where the product becomes visible.

🔴 THIS IS NOT A SECOND COPY OF `telluric_policy.gate_holding` -- IT IS STRICTER, ON
PURPOSE. That gate governs MEASUREMENT and deliberately lets `not-applied` Kitt Peak
through, because per-line clean-line selection (RYA-460/786) is a stated method defined
ON uncorrected data. RYA-1026 governs DISPLAY, and rules that basis out of the shipped
product. Two different questions, so two gates -- but ONE source of the underlying fact:
`telluric_applied` is read through `telluric_policy.applied_state`, never re-parsed here.
A hand-written clean-set would be a second declaration (the RYA-845 defect shape) and
would go stale silently in the permissive direction. A first draft of this module did
exactly that and got two of five entries wrong -- it listed `solar_delbouille`, an id
that does not exist (it is `solar_delbouille_liege`), and called it clean when the
registry records it `not-applied`.

THE ONE WHITELISTED EXCEPTION -- the KP2005 vs KP1984 pair is RETAINED as the telluric
CONTROL. 2005 is telluric-free and 1984 telluric-retaining, so the pair IS the molecfit
validation: removing the uncorrected half would destroy the only thing that proves the
correction works. A control is not a science product, and it is labelled as such.
"""
from __future__ import annotations

from pipeline.telluric_policy import applied_state

#: Uncorrected holdings permitted ONLY as the named control, never as science.
CONTROL_ONLY: dict[str, str] = {
    "solar_kpno": (
        "KP1984 is telluric-RETAINING and is whitelisted ONLY as the 1984-vs-2005 "
        "molecfit control (RYA-1026). Removing it would destroy the only pair that "
        "demonstrates the correction works. It is a control, not science."),
}

#: Holdings the registry calls `applied` where the applied-ness is KNOWN TO BE IN DOUBT.
#: These still render -- no value moves under RYA-1026 -- but the doubt travels with them
#: instead of being resolved silently in our favour.
REGISTRY_ANOMALIES: dict[str, str] = {
    "solar_iag": (
        "🔴 `iag_fts_solar_atlas` is catalogued telluric_basis=corrected, but the "
        "manifest routes it to the telluric-RETAINING Reiners+2016 file (46.25% of the "
        "O2 A-band below 0.5), so `telluric_policy.exclusion()` excludes nothing there. "
        "Two IAG atlases with opposite telluric states sit under one instrument_id. "
        "Found incidentally by RYA-944 and logged, not fixed -- it needs its own ticket. "
        "Until then an IAG product is displayable but its telluric state is UNVERIFIED "
        "in the O2/H2O bands."),
}


class TelluricDisplayError(RuntimeError):
    """A telluric-uncorrected product was about to be rendered as science."""


def assert_displayable(holding_id: str, *, as_control: bool = False,
                       where: str = "") -> None:
    """Gate a product at the point it becomes visible.

    `as_control=True` is the caller stating, explicitly and in code, that this render is
    the telluric control rather than a science claim. It is deliberately not inferable:
    if a caller could get the exception waived by accident, the gate would be decorative.
    """
    site = f" [{where}]" if where else ""

    if holding_id in CONTROL_ONLY:
        if as_control:
            return
        raise TelluricDisplayError(
            f"{holding_id}: REFUSING to render as science{site}. "
            f"{CONTROL_ONLY[holding_id]} Pass as_control=True to render it as the "
            f"control, or rebuild this product on the corrected sibling.")

    try:
        state = applied_state(holding_id)
    except KeyError as exc:
        raise TelluricDisplayError(
            f"{holding_id}: REFUSING to render as science{site}. This holding is not in "
            f"the registry at all, so its telluric state is not merely uncorrected -- it "
            f"is UNKNOWN, and an unknown must never be defaulted either way (RYA-806). "
            f"Register it at intake first. ({exc})") from exc

    if state == "applied":
        return

    raise TelluricDisplayError(
        f"{holding_id}: REFUSING to render as science{site}. Displayed science is "
        f"telluric-CORRECTED input only (RYA-1026); the registry records this holding "
        f"telluric_applied={state!r}. Note this is STRICTER than "
        f"`telluric_policy.gate_holding`, which may still permit this holding for "
        f"MEASUREMENT on the clean-line-selection basis -- passing that gate is not "
        f"permission to ship it. Rebuild on a corrected sibling "
        f"(solar_kpno_molecfit_corrected RYA-940, solar_harps_molecfit_corrected "
        f"RYA-931). Raw or uncorrected data is permitted only inside correction R&D for "
        f"other stars, never on a displayed product.")


def anomaly(holding_id: str) -> str | None:
    """A stated doubt about this holding's registered telluric state, or None.

    Kept separate from the pass/fail decision on purpose: an anomaly is not a refusal,
    and collapsing it into the boolean would lose the only thing that makes it
    actionable.
    """
    return REGISTRY_ANOMALIES.get(holding_id)


def display_state(holding_id: str) -> str:
    """CLEAN / CLEAN_WITH_ANOMALY / CONTROL_ONLY / BLOCKED / UNREGISTERED.

    For the tracker, which must show the REASON rather than silently omitting a blocked
    product -- an omission reads as 'no data', which is the absence-as-conclusion error
    (RYA-833). BLOCKED and UNREGISTERED stay distinct for the same reason: 'we hold this
    uncorrected' and 'we never wrote down what state this is in' are different problems
    with different fixes.
    """
    if holding_id in CONTROL_ONLY:
        return "CONTROL_ONLY"
    try:
        state = applied_state(holding_id)
    except KeyError:
        return "UNREGISTERED"
    if state == "applied":
        return "CLEAN_WITH_ANOMALY" if holding_id in REGISTRY_ANOMALIES else "CLEAN"
    return "BLOCKED"

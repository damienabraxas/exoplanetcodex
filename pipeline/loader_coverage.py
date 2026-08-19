"""Reconcile the HOLDINGS registry against what the band harness can actually reach.

RYA-904 (the guard RYA-897 specified; built here because 897 has not landed).

THE DEFECT CLASS
----------------
The registry says what we hold and what telluric state each holding is in. The band
harness says which holdings it can serve a window from. Nothing compared the two, so a
holding could be `verified`, pass every gate, and still be **unreachable** -- not by any
decision anybody made, but because no code could name it. That is what happened to both
telluric-CORRECTED CRIRES+ solar holdings: `measure_band_ew._LOADER_HOLDING` mapped
`crires_plus` to the ONE holding (`solar_vesta_crires_plus_idp`) that the telluric gate
correctly refuses, so the arm read as "telluric-blocked" when it was really "not wired".

A refusal is a decision and leaves a reason. Unreachability leaves NOTHING -- it looks
exactly like an arm that has no data, which is why it survived from RYA-806 to RYA-904.
This module makes it look like a failure instead.

WHY AN EXPLICIT GAP TABLE RATHER THAN SILENCE
---------------------------------------------
Not every registered holding should be servable by this harness, and a guard that
demanded it would be turned off within a week. But "we deliberately do not serve this
one" is a CLAIM, and RYA-833's rule applies: an absence is a hypothesis to state, not a
silence to leave. So every unreachable holding must appear in `DECLARED_GAPS` with a
reason and a ticket. The guard fails on holdings in neither set -- which is exactly the
state the corrected CRIRES+ holdings were in.

SCOPE. `system_id` defaults to `solar`, because `scripts/measure_band_*` is the SOLAR
band harness -- alpha Cen and Procyon holdings are served by other code and reporting
them here would be noise that gets the guard ignored. It is a parameter, not a constant,
so extending the harness to another system is a keyword rather than a rewrite.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HOLDINGS = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"


class LoaderCoverageError(AssertionError):
    """A holding we hold, may measure, and cannot name."""


#: holding_id -> why the band harness deliberately does not serve it. Every entry is a
#: STATED gap with an owner. Deleting an entry without wiring the holding makes the guard
#: fail, which is the point: this table can only shrink by building something.
#: 🔴 `solar_harps` WAS HERE, AND RYA-911 REMOVED IT because the loader now exists.
#: That removal was not optional and not tidy-up: `reconcile_loader_coverage` fails on a
#: DECLARED gap that is in fact WIRED, and it duly failed the moment the HARPS holding
#: was added -- which is the whole point of making the table rot loudly in both
#: directions (RYA-904). A stale exemption is an untrue statement about the code that
#: would wave through the next holding to go unreachable.
DECLARED_GAPS: dict[str, str] = {
    "elgueta2026_vizier": (
        "BY DESIGN — this is the UPSTREAM VizieR delivery (RYA-789), not a science-ready "
        "spectrum: sp/*.dat across several bands, Sirius-only, gitignored, md5-pinned, "
        "and carrying the unit trap RYA-794 documents (the ReadMe labels the wavelength "
        "column 0.1nm; it is nm). The harness serves its DERIVATIVE for the one window "
        "that has been conditioned -- `solar_crires_plus_y_rya794`, which is wired. "
        "Wiring the raw delivery as well would give one spectrum two doors and let two "
        "products disagree about which one they came from. The other bands (J/H/K) have "
        "no conditioned derivative yet; that is coverage owed, not a dispatch defect."),
}


def addressable_holdings() -> dict[str, str]:
    """holding_id -> instrument, for every holding the band harness can name.

    Read from the harness itself rather than restated here: a second list would be a
    second declaration of the same fact, which is how RYA-845's double-count survived.
    """
    import sys
    if str(ROOT / "scripts") not in sys.path:
        sys.path.insert(0, str(ROOT / "scripts"))
    from measure_band_ew import _INSTRUMENT_HOLDINGS
    return {h.holding_id: inst
            for inst, specs in _INSTRUMENT_HOLDINGS.items() for h in specs}


def unreachable_holdings(*, addressable: dict[str, str] | None = None,
                         holdings: pd.DataFrame | None = None,
                         declared_gaps: dict[str, str] | None = None,
                         system_id: str = "solar") -> list[dict]:
    """Every holding that is verified + gate-passing + unaddressable + undeclared.

    Arguments are injectable so the RYA-904 POSITIVE CONTROL can point an instrument back
    at only its raw holding and prove this fires. A guard that has never been seen red is
    not a guard.
    """
    from pipeline.telluric_policy import gate_holding
    addressable = addressable_holdings() if addressable is None else addressable
    declared_gaps = DECLARED_GAPS if declared_gaps is None else declared_gaps
    df = pd.read_csv(HOLDINGS) if holdings is None else holdings

    out: list[dict] = []
    for r in df.itertuples():
        if system_id is not None and str(r.system_id) != system_id:
            continue
        if str(r.evidence_state) != "verified":
            continue
        if str(r.holding_id) in addressable or str(r.holding_id) in declared_gaps:
            continue
        try:
            ok, why = gate_holding(str(r.holding_id), str(r.instrument_id))
        except Exception as e:
            # A holding whose telluric state nobody has determined is RYA-806's problem,
            # not this guard's. Recorded, not claimed as reachable.
            ok, why = False, f"{type(e).__name__}: {e}"
        if not ok:
            continue
        out.append(dict(holding_id=str(r.holding_id), instrument=str(r.instrument_id),
                        telluric_applied=str(r.telluric_applied), gate_reason=why))
    return out


def stale_gaps(*, addressable: dict[str, str] | None = None,
               declared_gaps: dict[str, str] | None = None) -> list[str]:
    """Declared gaps that are no longer gaps — the table rotting the other way.

    RYA-897 is wiring the HARPS solar arm as this ships. The moment it lands,
    `solar_harps` becomes addressable and its entry here becomes a false statement about
    the codebase that nothing would otherwise notice. A stale exemption is how a guard
    quietly stops guarding: the next holding to go unreachable could be waved through by
    an entry written for a problem that was already fixed.
    """
    addressable = addressable_holdings() if addressable is None else addressable
    declared_gaps = DECLARED_GAPS if declared_gaps is None else declared_gaps
    return sorted(set(declared_gaps) & set(addressable))


def reconcile_loader_coverage(**kw) -> None:
    """Raise `LoaderCoverageError` naming every holding we may measure and cannot reach."""
    stale = stale_gaps(**{k: v for k, v in kw.items()
                          if k in ("addressable", "declared_gaps")})
    if stale:
        raise LoaderCoverageError(
            f"{len(stale)} holding(s) are declared as GAPS in "
            f"pipeline.loader_coverage.DECLARED_GAPS and are in fact WIRED: "
            f"{', '.join(stale)}. Remove the entry — a stale exemption is an untrue "
            f"statement about the code that would wave through the next holding to go "
            f"unreachable.")
    bad = unreachable_holdings(**kw)
    if not bad:
        return
    lines = [
        f"{len(bad)} holding(s) are VERIFIED, PASS the telluric gate, and are "
        f"UNADDRESSABLE by the band harness (RYA-904).",
        "",
        "This is not a telluric refusal and not a coverage gap. Each of these is data we "
        "hold, in a state we are allowed to measure, that no loader can name — so it "
        "reads to every caller as if it did not exist.",
        "",
    ]
    for b in bad:
        lines.append(f"  {b['holding_id']}  ({b['instrument']}, "
                     f"telluric_applied={b['telluric_applied']})")
        lines.append(f"      gate says: {b['gate_reason'][:160]}")
    lines += [
        "",
        "Fix it in ONE of two ways, both of which are a decision someone makes on the "
        "record:",
        "  1. wire it — add a HoldingSpec to measure_band_ew._INSTRUMENT_HOLDINGS, with "
        "its reader and its OWN pre_normalised flag (both keys are holding-level, "
        "RYA-904); or",
        "  2. declare it — add it to pipeline.loader_coverage.DECLARED_GAPS with the "
        "reason and the ticket that owns it.",
        "",
        "What is NOT allowed is leaving it silent, because silence here is "
        "indistinguishable from having no data at all (RYA-833).",
    ]
    raise LoaderCoverageError("\n".join(lines))

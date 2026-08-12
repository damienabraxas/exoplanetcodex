"""Telluric policy — the single source for a question that keeps getting re-asked (RYA-786).

THE RECURRENCE THIS ENDS
------------------------
The telluric decision for a reference atlas was already made and built three times over:

  * **RYA-424** — telluric correction as a standing data-input stage, with
    instrument-aware routing (molecfit / cr2res / APERO) and an analysis-ready flag.
  * **`data/catalog/instrument_catalog.csv`** — `kpno_solar_atlas` is registered
    `telluric_required = no`; **RYA-380** is the molecfit/GDAS recipe for the instruments
    that do need it.
  * **RYA-460 / `config/physics_regime_rya400.yaml`** — per-line KPNO handling, e.g.
    K I "7665 stays in the O2 A-band, 7699 is the clean line".

Despite that, `SynthesisHandler.prepare` refused a KPNO run on a BAND-level
`telluric_required=True`, so every new band run re-collided with a settled question. That
is a single-source-of-truth defect, not a science gap, and this module is the source.

THE DISTINCTION THAT KEEPS GETTING LOST
---------------------------------------
`telluric_required = no` for a reference atlas does **NOT** mean "this atlas has been
telluric-divided". The Kurucz 1984 KPNO atlas HAS telluric absorption in it. It means
**the tellurics are handled by per-line CLEAN-LINE SELECTION rather than by a correction
stage** — the standard method for a reference atlas that has a clean alternative line and
a second arm (IAG) to cross-check against.

So the honest basis for running KPNO is *"instrument flag + per-line selection"*, cited.
It is NOT a `telluric_corrected` declaration, and fabricating one to satisfy a gate is
forbidden (RYA-786): it asserts a correction that was never applied.

WHAT THIS MODULE DECIDES
------------------------
  1. `requires_correction(instrument)` — from the catalog, the single registry of what an
     instrument needs. `yes` routes to the RYA-380/424 molecfit path; the correction
     machinery stays real and nothing is "avoided" as architecture.
  2. `exclusion(wave_A)` — the O2/H2O band set, enumerated ONCE below. A line inside a
     band is QUARANTINED-TELLURIC: a valid physics exclusion (RYA-777), not a cull.
  3. `gate(instrument, analysis_ready)` — what a handler should ask instead of carrying
     its own band flag.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CATALOG = ROOT / "data" / "catalog" / "instrument_catalog.csv"

# ── THE authoritative telluric exclusion set (RYA-786) ───────────────────────
#
# Enumerated once, here, and consumed everywhere. It was previously a three-entry list
# inside a measurement script — O2 A-band 7600-7640 plus two H2O bands — which is both
# incomplete and too narrow: the A-band runs to ~7685, the O2 B-band was absent entirely,
# and two H2O complexes in the red-optical/NIR were missing. A line sitting in an
# unlisted band is measured as if it were clean, which is the silent version of this bug.
#
# Ranges are air wavelengths in Angstrom, inclusive.
TELLURIC_BANDS: tuple[tuple[float, float, str], ...] = (
    (6867.0, 6884.0, "O2 B-band"),
    (7160.0, 7340.0, "H2O"),
    (7594.0, 7685.0, "O2 A-band"),
    (8100.0, 8400.0, "H2O"),
    (9280.0, 9600.0, "H2O"),
    (11120.0, 11560.0, "H2O"),
)

QUARANTINE_TAG = "QUARANTINED-TELLURIC"


def exclusion(wave_A: float) -> str:
    """Reason string if `wave_A` sits inside a telluric band, else ''."""
    for lo, hi, name in TELLURIC_BANDS:
        if lo <= float(wave_A) <= hi:
            return (f"{QUARANTINE_TAG}: inside the {name} ({lo:.0f}-{hi:.0f} A). The "
                    f"observed flux there is not stellar, so the line is excluded by "
                    f"per-line selection (RYA-460/786), not corrected.")
    return ""


def in_telluric_band(wave_A: float) -> bool:
    return bool(exclusion(wave_A))


_catalog_cache: dict = {}


def requires_correction(instrument: str) -> bool:
    """Does this instrument need a telluric CORRECTION STAGE? From the catalog only.

    Loud on an unknown instrument: the telluric state of an instrument the registry does
    not know cannot be asserted, and guessing it is how the fabricated declaration got in.
    """
    if "df" not in _catalog_cache:
        _catalog_cache["df"] = pd.read_csv(CATALOG)
    df = _catalog_cache["df"]
    hit = df[df.instrument_id.astype(str) == str(instrument)]
    if not len(hit):
        raise KeyError(
            f"instrument {instrument!r} is not in {CATALOG.name}; its telluric requirement "
            f"is unknown and must not be assumed. Register it first (RYA-786).")
    v = str(hit.iloc[0].get("telluric_required", "")).strip().lower()
    return v in ("yes", "true", "1", "required")


def gate(instrument: str, analysis_ready: bool = False) -> tuple[bool, str]:
    """(may_run, basis). What a handler asks instead of carrying a band flag.

    A reference atlas registered `telluric_required=no` runs on the per-line selection
    basis. An instrument that DOES require correction runs only once the RYA-424
    analysis-ready flag says the correction was applied and verified.
    """
    if not requires_correction(instrument):
        return True, (f"{instrument} is registered telluric_required=no; tellurics are "
                      f"handled by per-line clean-line selection over "
                      f"{len(TELLURIC_BANDS)} enumerated bands (RYA-460/786), not by a "
                      f"correction stage. No telluric_corrected declaration is made.")
    if analysis_ready:
        return True, (f"{instrument} requires correction and the RYA-424 analysis-ready "
                      f"flag is set: molecfit/GDAS applied and verified (RYA-380).")
    return False, (f"{instrument} is registered telluric_required=yes and the RYA-424 "
                   f"analysis-ready flag is not set. Route it through the RYA-380 "
                   f"molecfit path; before correction the observed flux is not stellar.")

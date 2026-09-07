"""Curated per-line exclusions — lines dropped on their SCIENTIFIC merits, with evidence.

RYA-1191. Ryan's ruling, 2026-09-07: *"Every graded line, ALL bands: telluric-corrected on
a VERIFIED-clean holding + blend/quality-checked. Drop ONLY genuine blends/artifacts."*

🔴 WHY A REGISTRY AND NOT A CUT. Every other exclusion in this pipeline is COMPUTED by the
fitter — a synthesis that did not converge, a COG inversion that did not bisect, a
telluric band the holding does not serve. Those need no registry because the data decides
each time. A blend does not work like that: the evidence is a property of the LINE LIST
and the solar spectrum, not of any one fit, and the judgement is human. So it is written
down once, with what was measured and who ruled on it, and applied everywhere.

⚠️ AND IT IS DELIBERATELY NOT A THRESHOLD. RYA-1191's blend audit ranks all 353 graded
lines by how much of their window's absorption belongs to the target species, and 164 of
them are not the dominant absorber. Auto-dropping on that statistic would remove a third
of the pool on a number nobody ratified — the RYA-981 error, a quota dressed as a cut
(RYA-161: validate, don't tune). The audit RANKS; this registry records the individual
calls, each with the evidence that earned it.

Read it with `excluded(species, wave_A)`; every row must carry a reason, the evidence and
the ticket that ruled, or the loader refuses it.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
REGISTRY = ROOT / "data" / "catalog" / "line_curation_exclusions.csv"

#: Match tolerance, Å. Wide enough to absorb the display rounding that differs between
#: canonical_gf and a synthesis list, far narrower than any fit window (RYA-1033: the key
#: is not a function of the value, so never match on a rounded string).
MATCH_TOL_A = 0.01

_cache: dict = {}


def _table() -> pd.DataFrame:
    if "t" not in _cache:
        if not REGISTRY.exists():
            _cache["t"] = pd.DataFrame(
                columns=["species", "wavelength_air_A", "band", "reason", "evidence",
                         "ticket", "ruled_by"])
        else:
            d = pd.read_csv(REGISTRY, comment="#")
            missing = [c for c in ("species", "wavelength_air_A", "reason", "evidence",
                                   "ticket") if c not in d.columns]
            if missing:
                raise ValueError(f"{REGISTRY.name} is missing {missing}")
            blank = d[["reason", "evidence", "ticket"]].isna().any(axis=1)
            if blank.any():
                raise ValueError(
                    f"{REGISTRY.name}: {int(blank.sum())} row(s) drop a line without a "
                    f"reason, evidence or a ticket. A line removed with no stated why is "
                    f"indistinguishable from a line removed because it was inconvenient "
                    f"(RYA-161).")
            _cache["t"] = d
    return _cache["t"]


def excluded(species: str, wave_A: float) -> str:
    """The reason this line is curated OUT, or '' if it is not.

    `species` is the canonical_gf spelling ('Fe I', 'Fe II')."""
    d = _table()
    if not len(d):
        return ""
    hit = d[(d.species.astype(str).str.strip() == str(species).strip())
            & ((d.wavelength_air_A.astype(float) - float(wave_A)).abs() <= MATCH_TOL_A)]
    if not len(hit):
        return ""
    r = hit.iloc[0]
    return f"CURATED-EXCLUSION ({r.ticket}): {r.reason} — {r.evidence}"


def excluded_wavelengths(species: str) -> list[float]:
    d = _table()
    if not len(d):
        return []
    m = d.species.astype(str).str.strip() == str(species).strip()
    return sorted(float(x) for x in d.loc[m, "wavelength_air_A"])

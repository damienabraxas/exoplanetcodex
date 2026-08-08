"""
State-changing surfaces -- the ONE authoritative registry (RYA-659).

A "state surface" is a file whose modification changes what the Codex claims is
TRUE RIGHT NOW: a verdict, a frozen reference value, an element's status, a
star/instrument/holdings identity, or the constants those are derived from. A
change to any of them obligates a review of CODEX_STATE_REGISTER.md.

WHY THIS MODULE EXISTS (single source of truth):
    The register drifted 11 state-changing tickets (RYA-556..652) because the
    obligation "bump the register" lived only in prose. RYA-643 changed element
    state and did not bump it. The list below is what makes that obligation
    checkable -- so it must exist exactly once. Consumers import from here; they
    do NOT keep their own copy:

      * scripts/check_register_freshness.py  (RYA-659) -- register drift guard
      * pipeline/ledger_consistency_guard.py (RYA-632) -- tracker-vs-verdict
        contradiction guard; it currently defines TRACKER_PATH / PHASE_C_PATH /
        PHYSICS_REGIME_PATH locally. Those three are LEDGER_PATHS entries here.
        When RYA-632 merges, repoint it at this module rather than duplicating.

Naming follows RYA-631: Catalog = enumeration, Register = mutable state ledger,
Tracker = per-item progress, Reference = frozen gold.
"""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

REPO_ROOT = Path(__file__).resolve().parents[1]

#: The mutable current-truth ledger this whole module exists to protect.
REGISTER = "CODEX_STATE_REGISTER.md"

#: The repo-root startup index naming the canonical read-set (RYA-659 Part B).
LEDGERS_INDEX = "LEDGERS.md"

# Named handles for the surfaces a consumer resolves INDIVIDUALLY rather than by
# sweeping the whole registry (RYA-632's guard needs three of them; RYA-654's tracker
# generator needs all four). They are defined here and referenced in STATE_SURFACES
# below, so each path literal exists in the repo exactly once -- renaming a surface is
# a one-line edit that cannot leave a consumer reading a stale path.
TRACKER = "data/audit/element_status_tracker.csv"
TRACKER_EDITORIAL = "data/audit/element_status_tracker_editorial.yaml"
PHASE_C_VERDICT_JSON = "data/audit/cno_synthesis/solar_phase_c_verdict.json"
PHYSICS_REGIME = "config/physics_regime_rya400.yaml"
#: RYA-676. The SSOT for "this row is owed -> this ticket resolves it". Hand-maintained,
#: joined onto the tracker's `refinement_debt` column by pipeline/refinement_debt_join.py.
REFINEMENT_REGISTRY = "data/audit/element_refinement_registry.csv"
#: RYA-663's per-element disposition report, in both its forms. The markdown is the
#: read-set member (a human reads it at session start); the JSON is its machine twin,
#: the same pairing phase_c already has.
DISPOSITION_REPORT_MD = "docs/audit/element_disposition_rya663.md"
DISPOSITION_REPORT_JSON = "data/audit/element_disposition_rya663.json"


class StateSurface(NamedTuple):
    """One state-changing surface. `path` is repo-relative, POSIX-style."""

    path: str
    owns: str          # what state this surface is authoritative for
    is_ledger: bool    # True => it is a LEDGERS.md read-set member


# The registry. Ordered: verdict channel -> frozen reference -> ledgers ->
# derived-constant sources. Every entry is a file a reviewer must notice.
STATE_SURFACES: tuple[StateSurface, ...] = (
    StateSurface(
        "scripts/phase_c_verdict_rya371.py",
        "the verdict GENERATOR -- classification logic, blank-cause tripwire",
        False,
    ),
    StateSurface(
        "docs/audit/solar_phase_c_verdict_rya371.md",
        "the live solar 27-element verdict (counts, per-element values)",
        False,
    ),
    StateSurface(
        PHASE_C_VERDICT_JSON,
        "the machine-readable twin of the verdict above",
        False,
    ),
    StateSurface(
        "scripts/build_solar_reference_v2_rya522.py",
        "the gold-reference BUILDER (freeze/tier logic)",
        False,
    ),
    StateSurface(
        "data/reference/solar/CURRENT",
        "which frozen gold vN is live -- changing it IS a re-freeze",
        False,
    ),
    StateSurface(
        TRACKER,
        "per-element status / tier / verdict (RYA-594; GENERATED since RYA-654)",
        True,
    ),
    StateSurface(
        # The hand-authored half of the tracker (RYA-654). Not a LEDGERS.md read-set
        # member -- you read the generated tracker, not its inputs -- but editing it
        # changes what the Codex claims about an element, so it is a state surface.
        TRACKER_EDITORIAL,
        "analyst half of the tracker: classification / action / vintages (RYA-654)",
        False,
    ),
    StateSurface(
        # A state surface but NOT a read-set member: you read the debt through the
        # tracker column it generates, not the registry itself. Editing it changes what
        # the Codex claims is OWED, which is a state claim like any other.
        REFINEMENT_REGISTRY,
        "which owed rows have a resolving ticket, and which have none (RYA-676)",
        False,
    ),
    StateSurface(
        # READ-SET MEMBER (RYA-676 §2C). RYA-663's report is what answers "can this
        # element flip, and if not what exactly is holding it" — the question RYA-672
        # had to reconstruct by hand because nothing routed a reader to this file.
        DISPOSITION_REPORT_MD,
        "per-element disposition: flip-or-blocked, with the three gates shown (RYA-663)",
        True,
    ),
    StateSurface(
        DISPOSITION_REPORT_JSON,
        "the machine-readable twin of the disposition report above",
        False,
    ),
    StateSurface(
        # A state surface but NOT a LEDGERS.md read-set member: it is consumed by
        # the tracker/verdict rather than read directly at session start.
        PHYSICS_REGIME,
        "per-element physics regime + NLTE-grid wiring claims (RYA-400)",
        False,
    ),
    StateSurface(
        "data/catalog/system_catalog.csv",
        "star identity + pipeline lifecycle stage (RYA-631)",
        True,
    ),
    StateSurface(
        "data/catalog/instrument_catalog.csv",
        "instrument capability + coverage (RYA-652)",
        True,
    ),
    StateSurface(
        "data/catalog/instrument_modes.csv",
        "per-instrument observing modes (RYA-652)",
        True,
    ),
    StateSurface(
        "data/catalog/holdings_manifest_registry.csv",
        "what data we already hold -- the anti-reinvent surface (RYA-652)",
        True,
    ),
    StateSurface(
        "config/stars.yaml",
        "STAR_PARAMS source -- the register's MIRROR block reads from it",
        False,
    ),
    StateSurface(
        "config/constants.py",
        "gates, corrections, grid selection; exposes STAR_PARAMS",
        False,
    ),
)

#: Convenience view: repo-relative paths only, in registry order.
SURFACE_PATHS: tuple[str, ...] = tuple(s.path for s in STATE_SURFACES)

#: The LEDGERS.md read-set members (RYA-632 needs a subset of these).
LEDGER_PATHS: tuple[str, ...] = tuple(s.path for s in STATE_SURFACES if s.is_ledger)


def existing_surfaces(root: Path | None = None) -> tuple[StateSurface, ...]:
    """Surfaces that are actually present on disk.

    A surface may legitimately not exist yet on an older commit or a partial
    checkout. Skipping a missing file is NOT a silent fallback -- the guard
    reports how many it checked, so a shrinking count is visible.
    """
    base = REPO_ROOT if root is None else root
    return tuple(s for s in STATE_SURFACES if (base / s.path).exists())

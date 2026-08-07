"""
Results-ledger consistency guard (RYA-632).

Cross-checks the element status tracker against the verdict artifacts and
loud-fails on any contradiction that is true regardless of which artifact is
canonical. Source-AGNOSTIC: it does not decide who is right; it refuses silent
disagreement. A stale tracker (disagreeing with the verdict artifacts) fails
here -- that is how "tracker updated on every merge" is enforced without a
merge hook.

Generalizes RYA-596's `_assert_blank_cause_is_honest` across all results
artifacts. Documented+ratified disagreements are allowed via an exceptions
file; undocumented disagreement fails.

RYA-436 sibling: that ticket guards the INPUT side (stored constant vs its
cited source); this one guards the OUTPUT side (verdict artifacts vs each
other and vs the human ledger).
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from pipeline import provenance_honesty, state_surfaces

REPO_ROOT = state_surfaces.REPO_ROOT
EXCEPTIONS_PATH = REPO_ROOT / "data" / "audit" / "known_verdict_divergences.yaml"

# Normalized verdict vocabulary. Map each artifact's raw verdict onto these in
# the adapter. counts convention (RYA-596): PASS / NLTE_OWED / CURATION_OWED / DATA_GAP.
#
# RYA-632 audit deviation from the brief: CURATION_OWED is NOT in MEASURED_VERDICTS.
# The brief seeded it there ("element has a real line pool"), but the live ledger
# refutes that: CURATION-OWED is the catch-all not-PASS bucket, and RYA-596 VERIFIED
# that Mg/Y/Zr/Eu are CURATION-OWED *and* genuinely zero-survivor (0 of 5/3/6/1 lines
# kept). Keeping it here made those four rows both claim-no-data and show-data at once,
# i.e. the guard would have manufactured a contradiction inside a single artifact on
# four rows whose state is verified correct. Line-pool evidence is carried by n_lines /
# frozen_value, which is exactly what RYA-596's tripwire keyed on.
MEASURED_VERDICTS = {"PASS", "NLTE_OWED"}     # verdict alone asserts a real line pool
NO_DATA_VERDICTS = {"DATA_GAP", "GET_DATA"}   # element claims zero usable data
# Raw cause strings that assert "no usable line" -- used for the honesty check.
# The canonical graded-cull claim is OWNED by pipeline/provenance_honesty.py
# (RYA-596/653) and imported, never re-spelled here: that module is what both
# EMITTERS assert with, so this guard and the emitters cannot drift apart on the
# wording. Listed below are only the additional spellings that module does not
# own -- the physics_regime ROUTING strings, and the older bare wording.
EXTRA_NO_SURVIVOR_PATTERNS = ("no line survives", "0 measured lines",
                              "no measured lines", "get-data")
OWED_TIERS = {"nlte_owed", "curation_owed", "owed"}


@dataclass(frozen=True)
class ArtifactState:
    artifact: str                    # "tracker" | "phase_c" | "gold" | "physics_regime"
    element: str
    verdict: Optional[str] = None    # one of the normalized verdicts above, or None
    n_lines: Optional[int] = None
    frozen_value: Optional[float] = None
    tier: Optional[str] = None
    raw_cause: Optional[str] = None  # any human cause string, checked for honesty


def _claims_no_data(s: ArtifactState) -> bool:
    if s.verdict in NO_DATA_VERDICTS:
        return True
    if s.raw_cause:
        # The RYA-596/653 claim, in the emitters' own words.
        if provenance_honesty.claims_zero_survivors(s.raw_cause):
            return True
        if any(p in s.raw_cause.lower() for p in EXTRA_NO_SURVIVOR_PATTERNS):
            return True
    return False


def _has_data(s: ArtifactState) -> bool:
    if s.n_lines is not None and s.n_lines > 0:
        return True
    if s.frozen_value is not None:
        return True
    if s.verdict in MEASURED_VERDICTS:
        return True
    return False


def check_element(states: list[ArtifactState], allowed_divergences: dict) -> list[str]:
    """Return a list of contradiction messages for one element (empty == clean)."""
    errs: list[str] = []
    if not states:
        return errs
    element = states[0].element

    # C1 -- honesty tripwire (generalized _assert_blank_cause_is_honest):
    # no artifact may claim zero data while any artifact shows data.
    no_data = [s for s in states if _claims_no_data(s)]
    has_data = [s for s in states if _has_data(s)]
    if no_data and has_data:
        if not _pair_annotated(no_data, has_data, element, allowed_divergences):
            errs.append(
                f"{element}: {[s.artifact for s in no_data]} claim no data while "
                f"{[(s.artifact, s.n_lines, s.frozen_value) for s in has_data]} show data "
                f"-- unrepresentable unless annotated+ratified in {EXCEPTIONS_PATH.name}"
            )

    # C2 -- owed freezes no value (2026-07-05 ratification: suspect -> held, not immortalised).
    for s in states:
        tier = (s.tier or "").strip().lower()
        if tier in OWED_TIERS and s.frozen_value is not None:
            errs.append(f"{element}: tier '{s.tier}' is owed but froze a value "
                        f"({s.frozen_value}) in {s.artifact}")

    # C4 -- cross-artifact verdict agreement (undocumented disagreement fails).
    verdict_by_artifact = {s.artifact: s.verdict for s in states if s.verdict is not None}
    if len(set(verdict_by_artifact.values())) > 1:
        # allow only if EVERY disagreeing pair is annotated+ratified
        if not _all_disagreeing_pairs_annotated(verdict_by_artifact, element,
                                                allowed_divergences):
            errs.append(
                f"{element}: verdict disagreement across artifacts {verdict_by_artifact} "
                f"-- annotate+ratify each pair in {EXCEPTIONS_PATH.name} or reconcile"
            )
    return errs


def _pair_annotated(group_a, group_b, element, allowed) -> bool:
    for a in group_a:
        for b in group_b:
            if a.artifact == b.artifact:
                # A single artifact that both claims no data and shows data is a
                # self-contradiction, not a cross-artifact divergence -- an exceptions
                # entry cannot express it, so it must always fail.
                return False
            if frozenset((a.artifact, b.artifact)) not in allowed.get(element, set()):
                return False
    return True


def _all_disagreeing_pairs_annotated(verdict_by_artifact, element, allowed) -> bool:
    """Only pairs that actually DISAGREE need an exceptions entry.

    (The brief's `_all_pairs_annotated` demanded an entry for every pair including the
    agreeing ones, which contradicted its own comment and would have forced us to
    "ratify" agreements. Intent kept, code matched to it.)
    """
    arts = sorted(verdict_by_artifact)
    for i in range(len(arts)):
        for j in range(i + 1, len(arts)):
            if verdict_by_artifact[arts[i]] == verdict_by_artifact[arts[j]]:
                continue
            if frozenset((arts[i], arts[j])) not in allowed.get(element, set()):
                return False
    return True


def check_counts(tracker_counts: dict, tallied_counts: dict) -> list[str]:
    """C3 -- tracker's declared counts must equal the tally over per-element verdicts.
    A stale tracker trips this. tracker_counts/tallied_counts keyed by verdict class."""
    errs = []
    keys = set(tracker_counts) | set(tallied_counts)
    for k in sorted(keys):
        if tracker_counts.get(k, 0) != tallied_counts.get(k, 0):
            errs.append(f"count mismatch [{k}]: tracker={tracker_counts.get(k, 0)} "
                        f"tallied={tallied_counts.get(k, 0)}")
    return errs


def load_allowed_divergences(path: Path = EXCEPTIONS_PATH) -> dict:
    """Each entry: {element, between:[artifactA,artifactB], reason, ratified_by}.
    An entry with no ratified_by is itself invalid (no annotation without a ticket,
    mirroring RYA-436's documented-divergence rule) -> raise."""
    import yaml
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as fh:
        entries = yaml.safe_load(fh) or []
    allowed: dict = {}
    for e in entries:
        if not e.get("ratified_by"):
            raise ValueError(f"exceptions file entry for {e.get('element')!r} has no ratified_by -- "
                             "no annotation without a ratifying ticket")
        pair = frozenset(e["between"])
        if len(pair) != 2:
            raise ValueError(f"exceptions file entry for {e.get('element')!r} must name exactly "
                             f"two distinct artifacts in `between`, got {e.get('between')!r}")
        allowed.setdefault(e["element"], set()).add(pair)
    return allowed


# ─────────────────────────────────────────────────────────────────────────────
#  ADAPTER -- the ONLY place that knows the live file formats (RYA-632 Part 1).
#  Same containment pattern as RYA-631's STAR_PARAMS accessor: every format
#  assumption lives below this line and nowhere else.
#
#  Confirmed live artifacts (RYA-632 audit, origin/main e4e98bb):
#
#    tracker         data/audit/element_status_tracker.csv        (RYA-594, HAND-EDITED,
#                    no generator exists) -- the human ledger of per-element standing.
#    phase_c         data/audit/cno_synthesis/solar_phase_c_verdict.json  (generated by
#                    scripts/phase_c_verdict_rya371.py) -- the LIVE verdict channel;
#                    authoritative for current verdict / A(X) / n_lines.
#    gold            data/reference/solar/solar_abundances_<CURRENT>.csv (RYA-469
#                    write-once freeze, built by scripts/build_solar_reference_v2_rya522.py)
#                    -- authoritative for what is FROZEN, i.e. tier + staked value.
#    physics_regime  config/physics_regime_rya400.yaml (RYA-400) -- authoritative for
#                    each element's regime + grid/data ROUTING verdict.
#
#  NOT adapted, deliberately: data/audit/element_status_tracker_drift.md (RYA-594) is
#  hand-written narrative markdown with no per-element machine-readable state. Per the
#  brief's "do not force it" rule it is left out rather than prose-parsed.
# ─────────────────────────────────────────────────────────────────────────────

# Resolved from the RYA-659 state-surface registry, NOT re-declared here: all three
# are registered surfaces, and state_surfaces.py is the single source of truth for
# where they live. Importing the names means a rename there reaches this guard.
TRACKER_PATH = REPO_ROOT / state_surfaces.TRACKER
PHASE_C_PATH = REPO_ROOT / state_surfaces.PHASE_C_VERDICT_JSON
PHYSICS_REGIME_PATH = REPO_ROOT / state_surfaces.PHYSICS_REGIME

# RYA-654 retired the tier->verdict inference. The tracker is now GENERATED and carries
# an explicit `verdict` column sourced from phase_c, so the guard reads the verdict it is
# checking instead of deducing it from a `tier` column that mixed two vocabularies
# ("gold PASS" beside "gf_floor"). `tier` now carries only the RYA-522 freeze tiers and is
# still read, for C2 (an owed tier may freeze no value).
#
# The Fe II arbiter row asserts no verdict (blank) and is skipped: phase_c carries a
# single Fe row annotated "62 Fe I + 3 Fe II", so keeping the diagnostic row in the tally
# would compare 27 tracker rows against 26 verdicts.

# physics_regime verdicts are ROUTING states (LOCKED | GET-GRID | GET-3D | GET-DATA |
# HARD-carry-forward | LTE-OK), not measurement verdicts. Only GET-DATA makes a claim
# about data existing; the rest say nothing about whether the element is measured, so
# they normalize to None rather than being forced onto the measurement vocabulary.
PHYSICS_REGIME_VERDICT_MAP = {"GET-DATA": "GET_DATA"}


def _normalize_verdict(raw: Optional[str]) -> Optional[str]:
    """'CURATION-OWED' -> 'CURATION_OWED'. Used for phase_c and gold, which share the
    RYA-371 verdict vocabulary."""
    if raw is None:
        return None
    v = str(raw).strip().upper().replace("-", "_")
    return v or None


def _int_or_none(v) -> Optional[int]:
    import pandas as pd
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return None
    return int(v)


def _float_or_none(v) -> Optional[float]:
    import pandas as pd
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return None
    return float(v)


def _read_tracker(path: Path = TRACKER_PATH):
    """Tracker rows -> (states, counts).

    Fe is carried at a FINER grain than the verdict (separate Fe I and Fe II rows; the
    verdict has one Fe row annotated '62 Fe I + 3 Fe II' -- drift log section B,
    'structural mismatch'). The Fe II arbiter row asserts no verdict, so it is dropped
    and the two tallies compare like with like.
    """
    import pandas as pd
    if not path.exists():
        raise FileNotFoundError(f"element status tracker not found at {path}")
    df = pd.read_csv(path, comment="#")
    if "verdict" not in df.columns:
        raise ValueError(
            f"{path.name} has no `verdict` column. Since RYA-654 the tracker is GENERATED "
            f"and carries the phase_c verdict explicitly -- regenerate it with "
            f"scripts/generate_element_status_tracker_rya654.py rather than reviving the "
            f"tier->verdict inference this replaced.")
    states, counts = [], {}
    for _, r in df.iterrows():
        el = str(r["element"]).strip()
        raw = r["verdict"]
        verdict = _normalize_verdict(_str_or_none(raw))
        if verdict is None:      # Fe II arbiter/diagnostic row
            continue
        counts[verdict] = counts.get(verdict, 0) + 1
        # NOTE: the tracker's `verdict_value` is a RECORDED measurement, not a freeze.
        # It must NOT map onto frozen_value -- C2 ("owed freezes no value") is an
        # invariant of the RYA-522 gold freeze, and mapping it here would fire C2 on
        # the 11 tracker rows that legitimately record a value at an owed tier.
        states.append(ArtifactState(
            artifact="tracker", element=el, verdict=verdict,
            tier=_str_or_none(r.get("tier")),
        ))
    return states, counts


def _read_phase_c(path: Path = PHASE_C_PATH):
    """Live verdict channel -> (states, counts)."""
    import json
    if not path.exists():
        raise FileNotFoundError(
            f"phase_c verdict artifact not found at {path} -- regenerate with "
            f"scripts/phase_c_verdict_rya371.py")
    payload = json.loads(path.read_text(encoding="utf-8"))
    states, counts = [], {}
    for r in payload["verdicts"]:
        verdict = _normalize_verdict(r.get("verdict"))
        if verdict:
            counts[verdict] = counts.get(verdict, 0) + 1
        states.append(ArtifactState(
            artifact="phase_c", element=str(r["element"]), verdict=verdict,
            n_lines=_int_or_none(r.get("n_lines")),
            raw_cause=r.get("channel"),
        ))
    return states, counts


def _str_or_none(v) -> Optional[str]:
    import pandas as pd
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return None
    return str(v).strip() or None


def _read_gold():
    """Frozen gold reference -> states. Read through pipeline.data_namespace so the
    CURRENT pointer stays the single source of truth for which version is live."""
    import pandas as pd

    from pipeline.data_namespace import read_solar_reference
    df, version = read_solar_reference()
    states = []
    for _, r in df.iterrows():
        a = r.get("A_X_nlte")
        a = a if pd.notna(a) else r.get("A_X")
        states.append(ArtifactState(
            artifact="gold", element=str(r["element"]).strip(),
            verdict=_normalize_verdict(_str_or_none(r.get("verdict"))),
            n_lines=_int_or_none(r.get("n_lines")),
            frozen_value=_float_or_none(a),
            tier=_str_or_none(r.get("confidence")),
            raw_cause=_str_or_none(r.get("note")),
        ))
    return states, version


def _read_physics_regime(path: Path = PHYSICS_REGIME_PATH):
    """RYA-400 regime map -> states. `indicators.in_data` carries the data claim in
    prose ('0 measured lines', '2 measured lines (marginal)', True); it is passed
    through as raw_cause so the honesty patterns see it."""
    import yaml
    if not path.exists():
        raise FileNotFoundError(f"physics regime map not found at {path}")
    doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    states = []
    for el, body in (doc.get("elements") or {}).items():
        raw = body.get("verdict")
        in_data = (body.get("indicators") or {}).get("in_data")
        states.append(ArtifactState(
            artifact="physics_regime", element=str(el),
            verdict=PHYSICS_REGIME_VERDICT_MAP.get(str(raw).strip()),
            raw_cause=(in_data if isinstance(in_data, str) else None),
        ))
    return states


def collect_element_states(repo_root: Path = REPO_ROOT):
    """ADAPTER -- read the confirmed live files into the normalized model.

    Returns (states_by_element, tracker_counts, tallied_counts). This is the ONLY
    place that knows the live file formats.
    """
    if Path(repo_root).resolve() != REPO_ROOT:
        # The gold reference is read through pipeline.data_namespace, whose CURRENT
        # pointer is bound to this checkout. Silently reading three artifacts from one
        # root and the fourth from another would be a false clean.
        raise ValueError(
            f"collect_element_states only supports this checkout ({REPO_ROOT}); the gold "
            f"reference is resolved via pipeline.data_namespace's CURRENT pointer, which "
            f"is not re-rootable. Got repo_root={repo_root}.")
    tracker_states, tracker_counts = _read_tracker()
    phase_c_states, tallied_counts = _read_phase_c()
    gold_states, _version = _read_gold()
    regime_states = _read_physics_regime()

    states_by_element: dict[str, list[ArtifactState]] = {}
    for s in (*tracker_states, *phase_c_states, *gold_states, *regime_states):
        states_by_element.setdefault(s.element, []).append(s)
    return states_by_element, tracker_counts, tallied_counts


# ─────────────────────────────────────────────────────────────────────────────
#  RYA-654 -- DOCUMENTED-OWED reds.
#
#  These elements are red ON PURPOSE and the entry below explains why. Read the
#  mechanism carefully: it ANNOTATES the error message and does NOT suppress the
#  error. The exit code is unchanged, the failure still prints, and nothing is
#  moved into known_verdict_divergences.yaml.
#
#  That asymmetry is the whole point. The ratified RYA-654 rule is that an
#  element whose measurement is not yet trusted has NOTHING to diverge from, so
#  annotating it as a ratified divergence would launder an un-done element into
#  a fake reconciliation. But a red with no explanation rots into "the guard is
#  just always red", which is how a real contradiction hides. So: say why, stay
#  red. If you ever find yourself wanting to make this suppress, the honest move
#  is to finish the measurement instead.
# ─────────────────────────────────────────────────────────────────────────────
_OWED_NOT_LAUNDERED = {
    "Sc": ("RYA-654 measurement-not-trusted: the only Sc value anywhere is the RYA-460 "
           "Kitt Peak Sc II 4246 blue-edge HFS single line (3.203), which phase_c itself "
           "holds LOW_CONFIDENCE -- 'HFS-resolved synthesis + a cleaner Sc II line owed "
           "before any PASS'. physics_regime's GET-DATA is CORRECT; there is nothing to "
           "reconcile to. Cleared by measuring Sc, NOT by an exceptions entry."),
}
# Ba's entry is GONE, and its removal is the point. It read "Cleared by the RYA-527 gold
# v3 re-freeze, NOT by this guard" -- RYA-665 performed exactly that freeze, so gold v3
# now carries Ba n_lines=1 + the honest RYA-559 synthesis cause and Ba is no longer an
# offender. Leaving the entry behind would be worse than dead code: it would annotate a
# FUTURE genuine Ba contradiction as "documented-owed" and drop it out of the
# undocumented count, which is the exact hiding place this guard exists to close.


def _annotate_documented(message: str) -> str:
    """Append the documented-owed explanation to an error line, if one applies.

    Matches on the leading 'El: ' the check_* messages already emit.
    """
    element = message.split(":", 1)[0].strip()
    note = _OWED_NOT_LAUNDERED.get(element)
    return f"{message}\n      DOCUMENTED-OWED, still failing: {note}" if note else message


def run_check(repo_root: Path = REPO_ROOT) -> int:
    allowed = load_allowed_divergences()
    states_by_element, tracker_counts, tallied_counts = collect_element_states(repo_root)
    errors = []
    for element in sorted(states_by_element):
        errors.extend(check_element(states_by_element[element], allowed))
    errors.extend(check_counts(tracker_counts, tallied_counts))
    if errors:
        undocumented = [e for e in errors
                        if e.split(":", 1)[0].strip() not in _OWED_NOT_LAUNDERED]
        print("RESULTS-LEDGER CONSISTENCY GUARD FAILED:", file=sys.stderr)
        for e in errors:
            print(f"  - {_annotate_documented(e)}", file=sys.stderr)
        print(f"\n  {len(errors)} failure(s): {len(errors) - len(undocumented)} documented-owed "
              f"(RYA-654/653 -- explained above, deliberately NOT annotated away), "
              f"{len(undocumented)} undocumented.", file=sys.stderr)
        return 1
    n = len(states_by_element)
    print(f"Results-ledger consistency guard OK: {n} elements, counts {tracker_counts}")
    return 0


if __name__ == "__main__":
    sys.exit(run_check())

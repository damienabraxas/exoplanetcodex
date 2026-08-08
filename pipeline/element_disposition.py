"""
pipeline/element_disposition.py  (RYA-663 — the pre-527 straggler sweep)
=======================================================================
ONE QUESTION, ANSWERED PER ELEMENT: *can this element flip to PASS before the
RYA-527 re-run, and if not, what exactly is holding it?*

WHY THIS IS GENERATED AND NOT A HAND AUDIT
------------------------------------------
RYA-524 asked the same question as a hand pass over 27 elements. It is stale:
its verdict is dated 2026-06-29 and the ledger has moved through RYA-553 (Fe
3D), 556 (N), 559 (Ba), 561 (floor promotion), 564 (Co), 592 (Mg), 643 (Sr/Co
re-runs) since. A hand audit of a moving ledger is stale the day after it is
written, which is the whole recurrence class RYA-436/632 exist to end. So this
is a generator: re-run it at the freeze instead of re-auditing.

THIS MODULE DECIDES NOTHING
---------------------------
The promotion rule is ALREADY RATIFIED and ALREADY IMPLEMENTED — Ryan's STRICT
three-gate rule (RYA-561, 2026-07-27), living in
``pipeline.engine_selection.evaluate_floor_promotion``. This module *applies*
that function and shows its working. It does not re-derive, relax or re-tune a
gate, and it must never grow its own copy of one: an element that wants a
different answer needs a ratification, not an edit here.

    gate 1  the Engine-B NLTE atom is RYA-534 anchor-validated
    gate 2  |A(X) - reference| <= 0.10 dex
    gate 3  a REAL cross-engine delta, |dCE| <= 0.10 — STRICT: missing FAILS,
            because a single-engine value has zero independent confirmation,
            and the atom delta may not stand in for it (that is gate 1 under a
            second name — validate-don't-tune, RYA-161).

TWO DISTINCTIONS THIS REPORT REFUSES TO BLUR
--------------------------------------------
1. **"CURATION-OWED" is not one bucket.** Five different states wear that label
   and each implies a different plan — an element whose pool was culled to zero
   (Mg) and an element carrying a perfectly good held value that merely has no
   second engine (Ca) are not the same problem. See ``OwedReason``.
2. **Gate 3 UNEVALUABLE is not gate 3 FAILED.** "No two-engine record exists for
   this element at all" points at the RYA-527 re-run; "a record exists and shows
   no cross-engine delta" points at the element. Reporting both as "failed"
   would hide which ones the re-run is about to fix.

STALENESS IS DETECTED, NOT ASSUMED
----------------------------------
The two-engine artifact is a REVIEW artifact and is known to lag (RYA-592). This
module does not guess: it cross-checks every value it reads against the live
verdict channel and reports the concrete contradictions it finds (see
``detect_stale_inputs``). Every input is stamped with its git provenance so the
report says what it read and how old it was.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from pipeline import data_namespace, state_surfaces
from pipeline.ratified_constraints import (    # RYA-674 emission-time gate
    assert_ratified_constraints_satisfied,
)
from pipeline.engine_selection import (
    FLOOR_PROMOTION,
    TwoEngineError,
    evaluate_floor_promotion,
    is_ratified_excluded_species,
    is_upper_limit_disposition,
    ratified_reported_ion,
)

REPO_ROOT = state_surfaces.REPO_ROOT

# ── inputs ───────────────────────────────────────────────────────────────────
# Each is declared once, with what it is authoritative FOR. phase_c is canonical
# for status (ratified in RYA-654); gold owns frozen/held VALUES; the two-engine
# record owns the cross-engine delta and nothing else.
PHASE_C_PATH = REPO_ROOT / state_surfaces.PHASE_C_VERDICT_JSON
TWO_ENGINE_PATH = REPO_ROOT / "data" / "audit" / "rya527_two_engine" / "solar_two_engine_records.json"
GOLD_HELD_PATH = REPO_ROOT / "data" / "reference" / "solar" / "solar_abundances_v1.csv"


def gold_current_path() -> Path:
    """The gold reference the CURRENT pointer names — resolved, never hardcoded.

    RYA-663 pinned this to ``solar_abundances_v2.csv``. RYA-665 then froze **v3** and
    moved CURRENT, which silently left this report reading a superseded reference: the
    ion column and the frozen values it compares against were one freeze behind, and
    nothing failed. A report whose whole job is to catch stale inputs must not itself
    be pinned to a version that a freeze can retire (RYA-669).
    """
    return data_namespace.reference_path(data_namespace.current_version())

#: Gold v2 BLANKS A_X on every owed row (the ratified RYA-522 tier holds the value
#: rather than freezing it). The held numbers therefore live in v1, and reading v2
#: for them would silently yield NaN and fail gate 2 on absence rather than on
#: physics. v1 is SUPERSEDED as the reference but is still the only record of what
#: the graded cull actually produced — that is why both are read.
HELD_VALUE_COLUMN = "A_X_nlte"


class DispositionError(RuntimeError):
    """Raised when an input is missing or self-contradictory. Never a silent skip."""


# ── the owed-reason taxonomy ─────────────────────────────────────────────────
class OwedReason:
    PASS = "PASS"
    OWED_HELD = "owed-HELD"
    OWED_BLANK = "owed-BLANK"
    MEASURED_UNFROZEN = "measured-awaiting-freeze"
    EW_POOL = "EW-pool"
    UPPER_LIMIT = "upper-limit"


#: What each bucket means, and what would move it. Kept next to the taxonomy so the
#: report never explains a bucket in prose that has drifted from the classifier.
OWED_REASON_PLAN = {
    OwedReason.PASS: "already PASS — nothing owed",
    OwedReason.OWED_HELD: (
        "pool survived and a value exists in gold v1, held unfrozen by the ratified "
        "RYA-522 tier. This is a PROMOTION DECISION: run the three gates."),
    OwedReason.OWED_BLANK: (
        "the graded cull left ZERO survivors — there is no value to promote. Needs "
        "line-pool / gf work (a real second line), not a ratification."),
    OwedReason.MEASURED_UNFROZEN: (
        "measured on a dedicated synthesis / Kitt Peak channel; the value exists but "
        "is not frozen. Clears at the RYA-527 v3 re-freeze."),
    OwedReason.EW_POOL: (
        "a plain EW pool sitting at the gf floor — a TIER question (owed vs gf_floor), "
        "adjudicated against Ti/Cr/Si, not a measurement gap."),
    OwedReason.UPPER_LIMIT: (
        "ratified UPPER_LIMIT disposition (RYA-563) — structurally never a PASS point "
        "value. Excluded from the flip denominator."),
}

#: Gate-3 states. The middle one is the finding the RYA-527 re-run exists to clear.
GATE3_OK = "OK"
GATE3_FAILED = "FAILED"
GATE3_UNEVALUABLE = "UNEVALUABLE"
#: The gates were never run — a ratified veto or an existing PASS short-circuited them.
#: Distinct from UNEVALUABLE, which means "we tried and the input does not exist".
GATE3_NA = "n/a"

#: Markdown table readability only — the JSON always carries the untruncated blocker.
_BLOCKER_MAX = 240
JSON_NAME = "element_disposition_rya663.json"


@dataclass
class ElementDisposition:
    element: str
    verdict: str
    owed_reason: str
    n_lines: Optional[int]
    channel: str
    reference: Optional[float]
    value: Optional[float]
    value_source: str
    cross_engine_delta: Optional[float]
    gate1: Optional[bool] = None
    gate2: Optional[bool] = None
    gate3_state: str = GATE3_UNEVALUABLE
    atom_citation: str = ""
    promoted: bool = False
    blocker: str = ""
    promotion_reason: str = ""

    def as_dict(self) -> dict:
        d = dict(self.__dict__)
        d["plan"] = OWED_REASON_PLAN[self.owed_reason]
        return d


# ── provenance ───────────────────────────────────────────────────────────────
def git_provenance(path: Path) -> dict:
    """Last commit that touched `path`, as {sha, date}.

    Uses the COMMIT time, not the mtime: a checkout or a `touch` rewrites mtimes
    and would make a stale artifact look fresh — the same reasoning RYA-659's
    freshness guard is built on.
    """
    if not path.exists():
        raise DispositionError(f"input not found: {path}")
    try:
        out = subprocess.run(
            ["git", "log", "-1", "--format=%h|%cI", "--", str(path)],
            cwd=REPO_ROOT, capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, OSError) as exc:      # noqa: BLE001
        raise DispositionError(f"cannot read git provenance for {path}: {exc}") from exc
    if not out:
        raise DispositionError(f"{path} is not tracked by git — refusing to report on it")
    sha, _, date = out.partition("|")
    return {"path": str(path.relative_to(REPO_ROOT)), "commit": sha, "committed": date}


# ── loaders ──────────────────────────────────────────────────────────────────
def load_phase_c(path: Path = PHASE_C_PATH) -> tuple[list[dict], dict]:
    """The canonical status channel (RYA-654 ratification)."""
    doc = json.loads(path.read_text(encoding="utf-8"))
    rows = doc.get("verdicts")
    if not rows:
        raise DispositionError(f"{path} carries no 'verdicts' — cannot report")
    return rows, doc.get("summary", {})


_ION_NUMERAL = {1: "I", 2: "II", 3: "III"}


def reported_ion(element: str, gold_ions: dict[str, str]) -> Optional[str]:
    """The ionisation stage the Codex REPORTS for `element`, as a numeral.

    Precedence, and why:

    1. **The NLTE registry lock** (`ratified_reported_ion`, RYA-558/240) — a ratified
       science decision (Cr I, Sr II, Ba II). It wins outright.
    2. **The gold reference's own `ion` column** — for elements the registry does not
       lock (Fe, Si, Ti ...), gold is the artifact that records what was reported.

    Gold is deliberately SECOND, not first: its `ion` column reads ``Sr I`` / ``Ba I``
    while the registry ratifies both as ion 2 — a legacy of the superseded raw-EW leg
    (the same Sr I-vs-Sr II wiring gap RYA-551/643 flagged for RYA-527). Trusting gold
    first would silently evaluate Sr against the wrong species.
    """
    locked = ratified_reported_ion(element)
    if locked is not None:
        return _ION_NUMERAL.get(locked)
    return gold_ions.get(element)


def load_two_engine(path: Path = TWO_ENGINE_PATH,
                    gold_ions: Optional[dict[str, str]] = None) -> dict[str, dict]:
    """Cross-engine deltas, keyed by element on the REPORTED ion.

    The artifact is per-SPECIES (Fe I and Fe II, Ti I and Ti II ...). Evaluating gate 3
    against the wrong species would test a number the Codex does not report — Cr II's
    dCE is -2.677 against Cr I's -1.058 — so the ion is resolved explicitly and a
    genuinely unresolvable element loud-fails rather than defaulting.
    """
    doc = json.loads(path.read_text(encoding="utf-8"))
    records = doc.get("records")
    if not records:
        raise DispositionError(f"{path} carries no 'records'")
    gold_ions = {} if gold_ions is None else gold_ions

    by_element: dict[str, dict] = {}
    seen: dict[str, list] = {}
    for rec in records:
        seen.setdefault(rec["element"], []).append(rec.get("ion"))
    for rec in records:
        el, ion = rec["element"], rec.get("ion")
        want = reported_ion(el, gold_ions)
        if len(seen[el]) > 1:
            if want is None:
                raise DispositionError(
                    f"{el}: two-engine artifact carries ions {seen[el]} and neither the NLTE "
                    f"registry nor the gold reference says which one is reported — refusing "
                    f"to pick one")
            if ion != want:
                continue
        by_element[el] = rec
    return by_element


def load_held_values(path: Path = GOLD_HELD_PATH) -> dict[str, tuple[float, str]]:
    """The owed-HELD values from gold v1, as ``element -> (value, ion)``.

    The ion is carried, not discarded, because gold v1 predates the ratified reported
    ion for at least one element: it holds **Sr I 4.961** (the superseded raw-EW leg)
    while the registry ratifies Sr as **Sr II**. Dropping the ion would let gate 2 test
    an Sr I value against an Sr II reference and report the resulting 2-dex miss as a
    physics failure, when it is a species mismatch. Fe likewise has two rows here.
    """
    df = pd.read_csv(path, comment="#")
    out: dict[str, tuple[float, str]] = {}
    for _, row in df.iterrows():
        val = row.get(HELD_VALUE_COLUMN)
        if pd.notna(val):
            out.setdefault(str(row["element"]), (float(val), str(row.get("ion") or "")))
    return out


def load_gold_ions(path: Optional[Path] = None) -> dict[str, str]:
    """element -> reported ion numeral, as the CURRENT gold reference records it."""
    df = pd.read_csv(path or gold_current_path(), comment="#")
    return {str(r["element"]): str(r["ion"]) for _, r in df.iterrows() if pd.notna(r.get("ion"))}


# ── classification ───────────────────────────────────────────────────────────
def classify_owed_reason(row: dict) -> str:
    """Bucket ONE element from its live verdict row.

    Derived from the row's own fields — never a hardcoded element list, so a
    measurement landing re-buckets the element automatically.
    """
    element = row["element"]
    if str(row.get("verdict", "")).upper() == "PASS":
        return OwedReason.PASS
    if is_upper_limit_disposition(element):
        return OwedReason.UPPER_LIMIT
    channel = str(row.get("channel") or "")
    if not (row.get("n_lines") or 0):
        return OwedReason.OWED_BLANK
    if "HELD at gold tier" in channel:
        return OwedReason.OWED_HELD
    if channel.lower().startswith("ew"):
        return OwedReason.EW_POOL
    return OwedReason.MEASURED_UNFROZEN


def _value_for(row: dict, reason: str, held: dict[str, tuple[float, str]],
               two_engine: dict[str, dict],
               gold_ions: Optional[dict[str, str]] = None) -> tuple[Optional[float], str]:
    """The A(X) gate 2 should test, and where it came from."""
    element = row["element"]
    if reason == OwedReason.OWED_HELD:
        if element in held:
            value, ion = held[element]
            want = reported_ion(element, gold_ions or {})
            if want and ion and ion != want:
                # Refuse to gate across species. This is a finding, not a gate failure.
                return None, (
                    f"ION MISMATCH: the only held value is {element} {ion} ({value:.3f}) "
                    f"but the ratified reported ion is {element} {want} — gate 2 would "
                    f"compare different species. Needs a {element} {want} value in gold.")
            return value, "gold v1 held value (RYA-522 tier holds it unfrozen)"
        return None, "no held value on record"
    if row.get("A_measured") is not None:
        return float(row["A_measured"]), "phase_c A_measured (live verdict channel)"
    rec = two_engine.get(element) or {}
    if rec.get("reported") is not None:
        return float(rec["reported"]), "two-engine reported (REVIEW artifact)"
    return None, "no value on record"


# ── the report ───────────────────────────────────────────────────────────────
def detect_stale_inputs(phase_c_rows: list[dict], two_engine: dict[str, dict]) -> list[str]:
    """Concrete contradictions between the two-engine artifact and the live channel.

    Not a heuristic and not a date comparison: where BOTH carry a value for an
    element and they disagree beyond rounding, the artifact demonstrably predates
    the live measurement. That is evidence, and it is what gets reported.
    """
    live = {r["element"]: r.get("A_measured") for r in phase_c_rows}
    stale = []
    for element, rec in sorted(two_engine.items()):
        reported, measured = rec.get("reported"), live.get(element)
        if reported is None or measured is None:
            continue
        if abs(float(reported) - float(measured)) > 0.001:
            stale.append(
                f"{element}: two-engine reports {float(reported):.3f} but the live verdict "
                f"channel measures {float(measured):.3f} "
                f"(delta {float(reported) - float(measured):+.3f}) — the two-engine artifact "
                f"predates that measurement")
    return stale


def value_disagreements(phase_c_rows: list[dict], two_engine: dict[str, dict],
                        gold_path: Optional[Path] = None,
                        tol: float = 0.001) -> list[dict]:
    """Every element whose A(X) differs across the artifacts that carry one.

    This is the "clean up the stale values" surface, and it exists because the
    RYA-632 ledger guard is **blind to it**: that guard compares verdicts, tiers and
    line counts, never values. Gold's frozen Fe 7.516 against a live 7.466 is invisible
    to it by construction. A disagreement here is not automatically an error — a frozen
    gold legitimately lags the live channel until the next re-freeze — but every one of
    them must be either explained or corrected before the v3 freeze, which is exactly
    what this list is for.
    """
    gold = pd.read_csv(gold_path or gold_current_path(), comment="#")
    gold_values = {str(r["element"]): r.get("A_X") for _, r in gold.iterrows()}
    out = []
    for row in sorted(phase_c_rows, key=lambda r: r["element"]):
        el = row["element"]
        sources = {}
        if row.get("A_measured") is not None:
            sources["phase_c (live)"] = float(row["A_measured"])
        gv = gold_values.get(el)
        if gv is not None and pd.notna(gv):
            sources["gold (frozen)"] = float(gv)
        rec = two_engine.get(el)
        if rec and rec.get("reported") is not None:
            sources["two-engine (review)"] = float(rec["reported"])
        if len(sources) < 2:
            continue
        spread = max(sources.values()) - min(sources.values())
        if spread > tol:
            out.append({"element": el, "spread": round(spread, 4), "values": sources})
    return out


def build_report(phase_c_path: Optional[Path] = None,
                 two_engine_path: Optional[Path] = None,
                 ticket: str = "RYA-663 per-element disposition (pre-527 straggler sweep)",
                 ) -> dict:
    """Assemble the full per-element disposition. Pure read; changes nothing.

    The two inputs are arguments rather than fixed module state because RYA-669 runs
    this report TWICE and the difference between the two runs is the finding: once on
    the 2026-07-18 two-engine record (what RYA-663 saw, gate 3 provisional) and once on
    the fresh Phase 2 record. A second copy of the classifier would let those two
    answers drift apart, which is the failure this module was written to end.
    """
    phase_c_rows, summary = load_phase_c(phase_c_path or PHASE_C_PATH)
    gold_ions = load_gold_ions()
    two_engine = load_two_engine(two_engine_path or TWO_ENGINE_PATH, gold_ions=gold_ions)
    held = load_held_values()

    dispositions: list[ElementDisposition] = []
    for row in sorted(phase_c_rows, key=lambda r: r["element"]):
        element = row["element"]
        reason = classify_owed_reason(row)
        value, value_source = _value_for(row, reason, held, two_engine, gold_ions)
        rec = two_engine.get(element)
        dce = None if rec is None else rec.get("mean_cross_engine_delta")

        d = ElementDisposition(
            element=element,
            verdict=str(row.get("verdict") or ""),
            owed_reason=reason,
            n_lines=row.get("n_lines"),
            channel=str(row.get("channel") or ""),
            reference=row.get("asplund2021"),
            value=value,
            value_source=value_source,
            cross_engine_delta=dce,
        )

        if reason in (OwedReason.PASS, OwedReason.UPPER_LIMIT):
            d.blocker = OWED_REASON_PLAN[reason]
            d.gate3_state = GATE3_NA          # short-circuited, not attempted
            dispositions.append(d)
            continue

        species = f"{element} {rec['ion']}" if rec and rec.get("ion") else None
        try:
            fp = evaluate_floor_promotion(element, value, d.reference, dce, species=species)
        except TwoEngineError as exc:      # loud, never a silent skip
            d.blocker = f"gate evaluation refused: {exc}"
            dispositions.append(d)
            continue

        d.gate1, d.gate2 = fp.gate1_atom_validated, fp.gate2_within_tol
        d.atom_citation = fp.atom_citation
        d.promoted, d.promotion_reason = fp.promoted, fp.reason
        # The distinction the ratified rule collapses but the re-run needs.
        if fp.gate3_cross_engine:
            d.gate3_state = GATE3_OK
        elif rec is None:
            d.gate3_state = GATE3_UNEVALUABLE
        else:
            d.gate3_state = GATE3_FAILED
        d.blocker = _blocker_for(d, reason)
        dispositions.append(d)

    stale = detect_stale_inputs(phase_c_rows, two_engine)

    # If the two-engine artifact is demonstrably behind the live channel ANYWHERE, then
    # every gate-3 number read from it is provisional — including on elements where no
    # contradiction could be detected, because an element whose live value is null (every
    # owed-HELD row) offers nothing to compare against. Silence there is absence of
    # evidence, not evidence of freshness. Any promotion resting on that delta is
    # therefore PROVISIONAL-on-the-re-run, and says so on its own row.
    provisional = bool(stale)
    for d in dispositions:
        if d.promoted and d.cross_engine_delta is not None and provisional:
            d.promotion_reason += (
                "  [PROVISIONAL: gate 3 read a cross-engine delta from the two-engine "
                "REVIEW artifact, which is demonstrably behind the live channel on "
                f"{len(stale)} element(s). Confirm on the RYA-527 re-run before freezing.]")

    # RYA-674 §2C: this report ADOPTS a value per element (`_value_for` will take the
    # two-engine reported number when phase_c has none), so it is an emission path even
    # though it changes nothing. That adoption is exactly where the Li 1.409 and Cr II
    # 5.676 leaks would re-enter.
    assert_ratified_constraints_satisfied(
        dispositions, "per-element disposition report (RYA-663)")

    return {
        "ticket": ticket,
        "thresholds": dict(FLOOR_PROMOTION),
        "phase_c_summary": summary,
        "inputs": [git_provenance(p) for p in
                   (phase_c_path or PHASE_C_PATH, two_engine_path or TWO_ENGINE_PATH,
                    GOLD_HELD_PATH, gold_current_path())],
        "stale_input_evidence": stale,
        "gate3_provisional": provisional,
        "value_disagreements": value_disagreements(phase_c_rows, two_engine),
        "can_flip_now": [d.element for d in dispositions if d.promoted],
        "dispositions": [d.as_dict() for d in dispositions],
    }


def render_markdown(report: dict) -> str:
    """The human view. Same data as the JSON — never a second, hand-kept narrative."""
    L: list[str] = []
    s = report["phase_c_summary"]
    L.append("# Per-element disposition — pre-527 straggler sweep (RYA-663)\n")
    L.append("**GENERATED — do not hand-edit.** Regenerate with "
             "`python scripts/gen_element_disposition.py`.\n")
    counts = " · ".join(f"{k} **{v}**" for k, v in (s.get("counts") or {}).items())
    L.append(f"Live channel: {counts} over {s.get('n_elements')} elements, "
             f"phase_c generated {s.get('generated')}.\n")
    L.append(f"Gates (ratified RYA-561, applied via `engine_selection."
             f"evaluate_floor_promotion`): tolerance **{report['thresholds']['tol_pass_dex']}** dex, "
             f"cross-engine **{report['thresholds']['cross_engine_dex']}** dex. "
             "Gate 3 is STRICT — a missing delta fails.\n")

    flip = report["can_flip_now"]
    L.append("## Can flip to PASS now\n")
    L.append(f"**{', '.join(flip) if flip else 'none'}**"
             + ("  — and this rests on a stale input; see below."
                if report["gate3_provisional"] and flip else ""))

    L.append("\n## Inputs read\n")
    L.append("| artifact | commit | committed |")
    L.append("|---|---|---|")
    for i in report["inputs"]:
        L.append(f"| `{i['path']}` | `{i['commit']}` | {i['committed']} |")

    if report["stale_input_evidence"]:
        L.append("\n## ⚠ Stale-input evidence\n")
        L.append("The two-engine record disagrees with the live verdict channel on these "
                 "elements, so it demonstrably predates them. **Every gate-3 number in this "
                 "report is read from that artifact.**\n")
        for e in report["stale_input_evidence"]:
            L.append(f"- {e}")

    if report["value_disagreements"]:
        L.append("\n## Value disagreements across artifacts (the cleanup list)\n")
        L.append("Invisible to the RYA-632 ledger guard, which compares verdicts and counts "
                 "but never values. Each must be explained or corrected before the v3 freeze.\n")
        L.append("| element | spread | values |")
        L.append("|---|---|---|")
        for v in report["value_disagreements"]:
            vals = "; ".join(f"{k} **{x:.3f}**" for k, x in v["values"].items())
            L.append(f"| {v['element']} | {v['spread']:.3f} | {vals} |")

    L.append("\n## Per-element\n")
    L.append("| El | verdict | bucket | A(X) | ref | dCE | g1 | g2 | g3 | blocker |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for d in report["dispositions"]:
        fmt = lambda x: "—" if x is None else f"{x:.3f}"      # noqa: E731
        tick = lambda b: "—" if b is None else ("✓" if b else "✗")   # noqa: E731
        # Ti's gate-1 citation is a paragraph-long RCA. Full text stays in the JSON;
        # the table keeps one scannable line per element.
        blocker = d["blocker"].replace("\n", " ")
        if len(blocker) > _BLOCKER_MAX:
            blocker = blocker[:_BLOCKER_MAX].rstrip() + f"… *(full text in {JSON_NAME})*"
        L.append(f"| {d['element']} | {d['verdict']} | {d['owed_reason']} | "
                 f"{fmt(d['value'])} | {fmt(d['reference'])} | {fmt(d['cross_engine_delta'])} | "
                 f"{tick(d['gate1'])} | {tick(d['gate2'])} | {d['gate3_state']} | "
                 f"{blocker} |")

    L.append("\n## What each bucket means\n")
    for reason, plan in OWED_REASON_PLAN.items():
        members = [d["element"] for d in report["dispositions"] if d["owed_reason"] == reason]
        L.append(f"- **{reason}** ({', '.join(members) if members else 'none'}) — {plan}")
    return "\n".join(L) + "\n"


def _blocker_for(d: ElementDisposition, reason: str) -> str:
    """The ONE thing standing between this element and PASS, in plain terms."""
    if d.promoted:
        return "none — promotes under the ratified three gates"
    if reason == OwedReason.OWED_BLANK and d.value is None:
        return ("no value exists (zero graded survivors) — needs a real line, "
                "not a decision")
    if reason == OwedReason.OWED_BLANK:
        # Mg's shape: the EW pool is blank, but a value exists off another channel.
        # Saying "no value exists" here while printing one would be the report
        # contradicting itself on a single row.
        return (f"EW pool is blank (0 graded survivors) but a value exists off another "
                f"channel [{d.value_source}] — reconcile the channels before promoting; "
                + "; ".join(_gate_failures(d)))
    return "; ".join(_gate_failures(d)) or d.promotion_reason


def _gate_failures(d: ElementDisposition) -> list[str]:
    """Each failing gate, with the number that failed it."""
    bits = []
    if d.gate1 is False:
        bits.append(f"gate 1: Engine-B atom not RYA-534-validated ({d.atom_citation})")
    if d.gate2 is False:
        bits.append(
            f"gate 2: |A - ref| = "
            f"{abs(d.value - d.reference):.3f} > {FLOOR_PROMOTION['tol_pass_dex']}"
            if d.value is not None and d.reference is not None else
            "gate 2: no value to test")
    if d.gate3_state == GATE3_UNEVALUABLE:
        bits.append("gate 3: UNEVALUABLE — no two-engine record (the RYA-527 re-run "
                    "is what produces one)")
    elif d.gate3_state == GATE3_FAILED:
        bits.append(
            f"gate 3: |dCE| = {abs(d.cross_engine_delta):.3f} > "
            f"{FLOOR_PROMOTION['cross_engine_dex']}"
            if d.cross_engine_delta is not None else
            "gate 3: single-engine record — zero independent confirmation")
    return bits

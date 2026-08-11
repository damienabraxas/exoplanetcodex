#!/usr/bin/env python3
"""
Element status tracker GENERATOR (RYA-654) -- RYA-436 Move B, results side.
==========================================================================

`data/audit/element_status_tracker.csv` used to be hand-maintained: the RYA-632 audit
found it had NO generator at all, and nothing but docs/CONVENTIONS.md even referenced
it. That is how it came to disagree with the live verdict channel on Co (still carrying
the DEMOTED blue-edge +1.188) and on N (still NLTE-OWED after RYA-556 cleared the debt),
which is what the RYA-632 count check trips on.

RYA-436's rule for exactly this shape of problem is GENERATE, DO NOT HAND-SYNC. So:

  * the STATUS half of every row is derived from its single source, every time;
  * the ANALYST half lives in a hand-authored sidecar the generator reads;
  * `--check` fails if the committed CSV differs from a fresh regeneration, which makes
    hand-editing the CSV the thing that breaks the build.

WHERE EACH COLUMN COMES FROM (one source per column, no exceptions)
-------------------------------------------------------------------
  element, verdict, verdict_value, sigma, n_lines, delta_vs_asplund, method
      -> data/audit/cno_synthesis/solar_phase_c_verdict.json. RATIFIED CANONICAL for
         the tracker's status columns (RYA-654 §1). `method` is the phase_c `channel`
         string VERBATIM -- the canonical source already states the method, so
         re-deriving a short label here would be inventing a second vocabulary.
  ion, regime_verdict
      -> config/physics_regime_rya400.yaml (RYA-400), the single source for an
         element's prescribed ion and its routing state.
  tier
      -> scripts/build_solar_reference_v2_rya522.confidence_of(), the RYA-522 ratified
         freeze tiers (gold / gf_floor / upper_limit / owed). NOT re-implemented here.
  refinement_debt
      -> data/audit/element_refinement_registry.csv via pipeline/refinement_debt_join.py
         (RYA-676). The tracker already said an element was owed; it never said whether
         a ticket to fix it existed, which is how RYA-581/585/565 sat in Backlog through
         eight architecture tickets. Registry membership is a human judgement against
         that file's admission rule; the RENDERING is generated, so it cannot be
         hand-patched into agreement.
  engine_a/b_model_vintage, classification, action_needed, source_tickets,
  editorial_updated, notes
      -> data/audit/element_status_tracker_editorial.yaml, hand-authored. These have no
         machine source; inventing one would be worse than admitting they are editorial.

WHAT IS DELIBERATELY *NOT* A SOURCE
-----------------------------------
The frozen gold reference. Per RYA-654 §1, gold vN owns frozen reference VALUES, not
status -- and gold v2 currently lags the live channel on Co and N (RYA-653 has the
corrected content ready as a candidate, pending the RYA-527 v3 re-freeze). Sourcing
status from a frozen snapshot is what made the tracker stale in the first place.

THE `verdict` COLUMN IS NEW
---------------------------
The old tracker had no verdict column and its `tier` column mixed two vocabularies
("gold PASS" next to "gf_floor"), so RYA-632's guard had to infer a verdict from a tier.
They are different things -- a verdict is live status, a tier is a freeze decision -- so
each now has its own column and its own source.

SIRIUS-ONLY PROVENANCE
----------------------
The tracker must be generated from the ratified run, never a local Mac canary. The
phase_c artifact carries no host/run-id field to assert on (flagged as owed below), so
the enforceable invariant is: generate only from the artifact as COMMITTED. If the
working-tree phase_c differs from its committed blob, this refuses to run -- an
uncommitted local re-run cannot leak into the ledger. The emitted header pins the phase_c
blob SHA and its last commit, so any row can be traced back to the exact run.

USAGE
-----
    python scripts/generate_element_status_tracker_rya654.py            # write
    python scripts/generate_element_status_tracker_rya654.py --check    # verify, exit 1
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import model_attempt_ledger as mal   # RYA-695 `models_tried`
from pipeline import state_surfaces  # noqa: E402

# Every file this generator reads or writes is a registered state surface (RYA-659).
# Taken from the registry's named handles, never re-spelled here, so a rename there is
# a one-line edit that cannot leave this generator reading a stale path.
TRACKER_REL = state_surfaces.TRACKER
PHASE_C_REL = state_surfaces.PHASE_C_VERDICT_JSON
REGIME_REL = state_surfaces.PHYSICS_REGIME
EDITORIAL_REL = state_surfaces.TRACKER_EDITORIAL

TRACKER_PATH = ROOT / TRACKER_REL
PHASE_C_PATH = ROOT / PHASE_C_REL
REGIME_PATH = ROOT / REGIME_REL
EDITORIAL_PATH = ROOT / EDITORIAL_REL

#: The Fe II row is an ARBITER/diagnostic, not an element verdict -- phase_c carries one
#: Fe row ("62 Fe I + 3 Fe II"). It is emitted from the ratified arbiter constant so the
#: tracker keeps carrying it without anyone hand-typing the number.
FE_II_KEY = "Fe_II"

COLUMNS = [
    "element", "ion", "verdict", "verdict_value", "sigma", "n_lines",
    "delta_vs_asplund", "tier", "method", "regime_verdict",
    "two_engine_wiring_status", "chosen_engine", "selection_reason", "models_tried",
    "refinement_debt", "engine_reach",
    "engine_a_model_vintage", "engine_b_model_vintage", "classification",
    "action_needed", "source_tickets", "editorial_updated", "notes",
]

GENERATED_COLUMNS = frozenset({
    "element", "ion", "verdict", "verdict_value", "sigma", "n_lines",
    "delta_vs_asplund", "tier", "method", "regime_verdict",
    "two_engine_wiring_status", "chosen_engine", "selection_reason", "models_tried",
    "refinement_debt", "engine_reach",
})

#: RYA-695. The tracker said which engines COVER an element
#: (`two_engine_wiring_status`) and never which one was CHOSEN, why, or what else had
#: been tried. Ryan: "it should show which element was chosen, where we got the damn
#: thing and why we chose it. But it should also reference any model we tried. ... we
#: keep repeating work."
#:
#: `chosen_engine` / `selection_reason` come from the two-engine record, which since
#: RYA-695 carries the per-line `reason` and `regime` the RYA-525 selector always
#: computed and the emitter always discarded. `models_tried` is a JOIN over records
#: that already exist (grid provenance `supersedes` chains, the Engine-B deck ledger,
#: the availability table) — see pipeline/model_attempt_ledger.py. No column here is
#: hand-filled; an element with nothing on record says "none on record" rather than
#: emitting a blank that reads as "never investigated".
TWO_ENGINE_REL = "data/audit/rya527_phase3/solar_two_engine_records.json"
TWO_ENGINE_PATH = ROOT / TWO_ENGINE_REL

#: RYA-673. The tracker already says what each element's verdict IS; it never said how
#: many engines that verdict rests on. Those are different facts, and the second one is
#: what Beta's "best of abilities on all engines" bar is actually about — an element
#: reporting a clean PASS on one engine with no cross-check looked identical here to one
#: confirmed on both.
WIRING_REL = "data/audit/two_engine_wiring_audit.csv"
WIRING_PATH = ROOT / WIRING_REL

#: RYA-676. The tracker said an element was `owed`; Linear said a ticket to fix it was
#: `Backlog`; nothing joined the two, so RYA-581/585/565 sat unfired through eight
#: architecture tickets. `refinement_debt` is that join, generated from
#: `data/audit/element_refinement_registry.csv` — hand-maintained, because whether a debt
#: EXISTS is a judgement against that file's admission rule, not something derivable.
#: It is a GENERATED column here for the same reason every other one is: so it cannot be
#: hand-patched into agreement.


class SourceError(RuntimeError):
    """A source is missing, uncommitted or inconsistent. Never degraded to a default."""


# ─────────────────────────────────────────────────────────────────────────────
#  git provenance -- loud on every failure (a failed lookup is not a clean run)
# ─────────────────────────────────────────────────────────────────────────────
def _git(*args: str) -> str:
    try:
        out = subprocess.run(["git", *args], cwd=ROOT, text=True,
                             capture_output=True, check=True)
    except FileNotFoundError as exc:
        raise SourceError(f"git executable not found: {exc}") from exc
    except subprocess.CalledProcessError as exc:
        raise SourceError(
            f"git {' '.join(args)} failed (exit {exc.returncode}): "
            f"{(exc.stderr or '').strip()}") from exc
    return out.stdout.strip()


def phase_c_provenance() -> dict:
    """Blob SHA + last commit of the COMMITTED phase_c artifact, and a hard refusal to
    generate from a working tree that has diverged from it.

    This is the Sirius-only lever we can actually enforce today: the ratified run is the
    one that was committed and reviewed. A Mac canary re-run sitting uncommitted in the
    working tree is exactly the thing that must not reach the ledger.
    """
    if _git("status", "--porcelain", "--", PHASE_C_REL):
        raise SourceError(
            f"{PHASE_C_REL} has uncommitted changes. The tracker is generated from the "
            f"RATIFIED (committed) phase_c run, never from a local re-run -- commit the "
            f"verdict artifact first, or `git checkout -- {PHASE_C_REL}` to drop it.")
    blob = _git("rev-parse", f"HEAD:{PHASE_C_REL}")
    commit = _git("log", "-1", "--format=%H", "--", PHASE_C_REL)
    when = _git("log", "-1", "--format=%cI", "--", PHASE_C_REL)
    if not commit:
        raise SourceError(f"{PHASE_C_REL} has no commit history -- it is not a ratified run")
    return {"blob": blob, "commit": commit, "committed_at": when}


# ─────────────────────────────────────────────────────────────────────────────
#  sources
# ─────────────────────────────────────────────────────────────────────────────
def load_phase_c() -> dict:
    if not PHASE_C_PATH.exists():
        raise SourceError(f"phase_c verdict artifact not found at {PHASE_C_REL} -- "
                          f"regenerate with scripts/phase_c_verdict_rya371.py")
    return json.loads(PHASE_C_PATH.read_text(encoding="utf-8"))


def load_regime() -> dict:
    import yaml
    if not REGIME_PATH.exists():
        raise SourceError(f"physics regime map not found at {REGIME_REL}")
    return (yaml.safe_load(REGIME_PATH.read_text(encoding="utf-8")) or {}).get("elements") or {}


def load_editorial() -> dict:
    import yaml
    if not EDITORIAL_PATH.exists():
        raise SourceError(
            f"editorial sidecar not found at {EDITORIAL_REL} -- it holds the columns "
            f"that have no machine source; the generator will not invent them")
    return yaml.safe_load(EDITORIAL_PATH.read_text(encoding="utf-8")) or {}


def load_wiring() -> dict:
    """(element, ion) -> two-engine wiring status, from the RYA-673 audit.

    Loud on absence like every other source here: a blank wiring column would read as
    "no problem" on precisely the elements the audit exists to flag.
    """
    if not WIRING_PATH.exists():
        raise SourceError(
            f"two-engine wiring audit not found at {WIRING_REL} -- regenerate it with "
            f"scripts/rya673_two_engine_wiring_audit.py (needs Sirius: it drives both "
            f"engines over real solar data)")
    import csv as _csv
    with WIRING_PATH.open(encoding="utf-8") as fh:
        return {(r["element"], r["ion"]): r["wiring_status"]
                for r in _csv.DictReader(fh)}


def load_wiring_rows() -> dict:
    """(element, ion) -> the FULL audit row (RYA-695).

    `load_wiring` returns only the status string, which is all the
    `two_engine_wiring_status` column needs. `selection_reason` additionally needs the
    audit's `single_engine_reason` — the cited evidence that a missing engine is
    impossible rather than unfinished — so the whole row is read here rather than
    widening the existing loader's contract.
    """
    if not WIRING_PATH.exists():
        raise SourceError(f"two-engine wiring audit not found at {WIRING_REL}")
    import csv as _csv
    with WIRING_PATH.open(encoding="utf-8") as fh:
        return {(r["element"], r["ion"]): r for r in _csv.DictReader(fh)}


def tier_of(element: str) -> str:
    """RYA-522 ratified freeze tier, read from the gold builder -- not re-implemented."""
    from build_solar_reference_v2_rya522 import confidence_of
    return confidence_of(element)


# ─────────────────────────────────────────────────────────────────────────────
#  row assembly
# ─────────────────────────────────────────────────────────────────────────────
def _fmt(value, spec: str = "") -> str:
    if value is None:
        return ""
    if spec and isinstance(value, (int, float)):
        return format(value, spec)
    return str(value)


def _wiring_for(wiring: dict, element: str, ion: str) -> str:
    """Wiring status for one row, matched on the ion the tracker reports.

    The audit is keyed per SPECIES because wiring is per species — Cr I and Cr II do not
    share an Engine-B atom. A tracker ion that finds no audit row is a real disagreement
    between the two about which species the Codex reports, so it raises rather than
    emitting a blank that would read as "fine".
    """
    if (element, ion) in wiring:
        return wiring[(element, ion)]

    # A COMBINED ion token ("Fe" -> "I+II"): the tracker carries one row for an element
    # analysed on both stages, while wiring is per species because it IS per species --
    # Fe I and Fe II do not share an Engine-B atom. Report each stage rather than picking
    # one, so a row that is wired on one stage and not the other cannot read as clean.
    if "+" in ion:
        parts = [p.strip() for p in ion.split("+") if p.strip()]
        statuses = [wiring.get((element, p)) for p in parts]
        if all(statuses):
            uniq = set(statuses)
            return statuses[0] if len(uniq) == 1 else "; ".join(
                f"{p}:{s}" for p, s in zip(parts, statuses))

    same_element = {k[1] for k in wiring if k[0] == element}
    raise SourceError(
        f"{element} {ion}: no row in {WIRING_REL} (it has {element} on "
        f"{sorted(same_element) or 'no ion'}). The tracker and the wiring audit "
        f"disagree about which species is reported -- reconcile them, do not blank it")


def _editorial_row(editorial: dict, key: str, element: str) -> dict:
    if key not in editorial:
        raise SourceError(
            f"no editorial entry for {key!r} in {EDITORIAL_REL}. Every row needs one: "
            f"add it (classification / action_needed / vintages / source_tickets / "
            f"notes) rather than letting {element} emit blank analyst columns")
    row = editorial[key]
    missing = [f for f in ("engine_a_model_vintage", "engine_b_model_vintage",
                           "classification", "action_needed", "source_tickets",
                           "editorial_updated", "notes") if f not in row]
    if missing:
        raise SourceError(f"editorial entry {key!r} is missing field(s) {missing}")
    return row


def load_refinement_debt() -> dict:
    """element -> rendered `refinement_debt` cell (RYA-676).

    Loud on absence, like every other source here: a blank column would read as
    "nothing owed" on precisely the rows the registry exists to make visible.
    """
    from pipeline import refinement_debt_join
    try:
        return refinement_debt_join.by_element()
    except refinement_debt_join.RegistryError as exc:
        raise SourceError(str(exc)) from exc


def load_two_engine() -> dict:
    """element -> the two-engine record for its REPORTED ion (RYA-695).

    Loud on absence, like every other source: a blank `chosen_engine` column would
    read as "no engine was chosen", which is a much stronger claim than "the record
    was not regenerated".
    """
    if not TWO_ENGINE_PATH.exists():
        raise SourceError(
            f"two-engine record not found at {TWO_ENGINE_REL} -- it holds the columns "
            f"`chosen_engine` / `selection_reason`. Regenerate it on Sirius with "
            f"scripts/rya527_two_engine_run.py --out-dir data/audit/rya527_phase3")
    doc = json.loads(TWO_ENGINE_PATH.read_text(encoding="utf-8"))
    by: dict[str, list[dict]] = {}
    for r in doc.get("records", []):
        by.setdefault(str(r["element"]), []).append(r)
    return by


_ENGINE_LABEL = {"engineA_1dnlte": "A", "engineB_synth": "B"}


def _engine_cells(two_engine: dict, wiring: dict, element: str, ion: str) -> tuple:
    """(chosen_engine, selection_reason) for one species.

    Four outcomes, all of them meaningful and none of them blank:
      * both engines won lines        -> "both (aggregated)" + the per-reason split
      * one engine won every line     -> that engine + why it won
      * no two-engine record           -> the wiring audit's ratified single-engine
                                          reason, or an explicit "not in the record"
    """
    recs = [r for r in two_engine.get(element, []) if str(r.get("ion")) == ion] \
        or two_engine.get(element, [])
    if not recs:
        w = wiring.get((element, ion)) or {}
        single = str(w.get("single_engine_reason") or "").strip()
        if single:
            return "single-engine (ratified)", single
        return "none", ("element produces no two-engine record -- it is absent from "
                        "both engines' coverage; see two_engine_wiring_audit.csv")
    rec = max(recs, key=lambda r: r.get("n_lines") or 0)
    engines = [_ENGINE_LABEL.get(e, e) for e in (rec.get("selected_engines") or [])]
    if not engines:
        return "none", "two-engine record carries no selected engine"
    chosen = "both (aggregated)" if len(engines) > 1 else engines[0]
    bits = []
    for sr in (rec.get("selection_reasons") or []):
        bits.append(f"{_ENGINE_LABEL.get(sr.get('engine'), sr.get('engine'))}"
                    f" x{sr.get('n_lines')}: {sr.get('reason')}")
    reason = " | ".join(bits) if bits else "no per-line reason recorded"
    if rec.get("mix_flagged"):
        reason = (f"MIX FLAGGED (mean cross-engine delta "
                  f"{rec.get('mean_cross_engine_delta')} exceeds the gate -- adjudicate, "
                  f"do not trust the mean): " + reason)
    elif rec.get("cross_engine_mix"):
        reason = (f"cross-engine mix, within gate (mean delta "
                  f"{rec.get('mean_cross_engine_delta')}): " + reason)
    if rec.get("engineB_source"):
        reason += f" || Engine-B source: {rec['engineB_source']}"
    return chosen, reason


def build_rows() -> list[dict]:
    phase_c = load_phase_c()
    regime = load_regime()
    editorial = load_editorial()
    wiring = load_wiring()
    debt = load_refinement_debt()
    two_engine = load_two_engine()
    wiring_rows = load_wiring_rows()
    attempts = mal.attempts_by_element()
    reach = _engine_reach_table()

    rows = []
    for rec in phase_c["verdicts"]:
        element = str(rec["element"]).strip()
        spec = regime.get(element)
        if spec is None:
            raise SourceError(
                f"{element} is in the phase_c verdict but not in {REGIME_REL} -- the two "
                f"element sets must agree; map it there rather than emitting a blank ion")
        ed = _editorial_row(editorial, element, element)
        ion = str(spec.get("ion", "")).strip()
        rows.append({
            "element": element,
            "ion": ion,
            "verdict": _fmt(rec.get("verdict")),
            "verdict_value": _fmt(rec.get("A_measured"), ".3f"),
            "sigma": _fmt(rec.get("sigma"), ".3f"),
            "n_lines": _fmt(rec.get("n_lines")),
            "delta_vs_asplund": _fmt(rec.get("delta_vs_asplund"), "+.3f"),
            "tier": tier_of(element),
            "method": _fmt(rec.get("channel")),
            "regime_verdict": str(spec.get("verdict", "")).strip(),
            "two_engine_wiring_status": _wiring_for(wiring, element, ion),
            **dict(zip(("chosen_engine", "selection_reason"),
                       _engine_cells(two_engine, wiring_rows, element, ion))),
            "models_tried": mal.render_cell(element, attempts),
            "refinement_debt": _debt_cell(debt, element),
            "engine_reach": _engine_reach_cell(reach, element, ion),
            **{k: str(ed[k]) for k in (
                "engine_a_model_vintage", "engine_b_model_vintage", "classification",
                "action_needed", "source_tickets", "editorial_updated", "notes")},
        })

    rows.append(_fe_ii_row(editorial, wiring, debt, two_engine, wiring_rows))

    # A registry row naming an element no tracker row carries would render nowhere --
    # the same silent drop this column exists to stop, one level up. Checked against the
    # rows actually emitted, so it cannot pass on a stale element set.
    from pipeline import refinement_debt_join
    try:
        refinement_debt_join.assert_elements_known({r["element"] for r in rows})
    except refinement_debt_join.RegistryError as exc:
        raise SourceError(str(exc)) from exc
    return rows


def _debt_cell(debt: dict, element: str) -> str:
    from pipeline import refinement_debt_join
    return refinement_debt_join.debt_cell(element, debt)


#: RYA-776. The tracker names WHICH grid and WHAT STATE per engine and has never carried
#: the WAVELENGTH REACH, so "do we have Engine A on Fe in the IR?" was unanswerable here
#: and got re-derived by hand every time it came up (RYA-763 was that re-derivation).
#:
#: This is a JOIN, not an absorption. The keyed (element, ion, engine, grid, band) table
#: is the SIBLING file `data/catalog/engine_coverage.csv`; the tracker carries one compact
#: cell pointing into it -- bands the engine SERVES, then bands it only REACHES marked
#: `?`. Putting the full table in a tracker row would bloat it and, worse, fork the
#: coverage answer into two places that drift apart.
def _engine_reach_table():
    """The generated engine-coverage rows, or None where the reference is not built.

    Returns None rather than raising. The reach table is generated on SIRIUS (the grids
    are never on the Mac), so a Mac regeneration of the tracker must not fail on its
    absence -- the tracker's own sources are all committed and its verdict logic does not
    depend on this column. A missing table prints as "(not generated)", which is visibly
    different from an engine that genuinely reaches nothing.
    """
    from pipeline import coverage
    try:
        return coverage.load_engine_coverage()
    except coverage.CoverageError:
        return None


def _engine_reach_cell(table, element: str, ion: str) -> str:
    if table is None:
        return "(not generated — run scripts/generate_engine_coverage_rya776.py on Sirius)"
    from pipeline import coverage
    return coverage.engine_summary(element, ion, table)


def _fe_ii_row(editorial: dict, wiring: dict, debt: dict,
               two_engine: dict, wiring_rows: dict) -> dict:
    """The Fe II ionization-arbiter row, from the ratified arbiter constant (RYA-305/341/405).

    It asserts no element verdict -- it is the diagnostic the Fe I anchor is gated on --
    so `verdict` is blank and `tier` is `diagnostic`. RYA-632's reader drops it for the
    count tally on exactly that basis.
    """
    from config.constants import FE_IONIZATION_SYNTH_ARBITER
    solar = FE_IONIZATION_SYNTH_ARBITER.get("solar")
    if not solar or solar.get("fe2_synth") is None:
        raise SourceError(
            "FE_IONIZATION_SYNTH_ARBITER['solar'] carries no fe2_synth -- the Fe II "
            "arbiter row has no source. Fix the constant; do not hand-type the value.")
    ed = _editorial_row(editorial, FE_II_KEY, "Fe II")
    return {
        "element": "Fe",
        "ion": "II",
        "verdict": "",
        "verdict_value": _fmt(solar["fe2_synth"], ".3f"),
        "sigma": _fmt(solar.get("fe2_uncertainty"), ".3f"),
        "n_lines": "",
        "delta_vs_asplund": "",
        "tier": "diagnostic",
        "method": f"arbiter synthesis -- {solar.get('provenance', '')}"
                  f" (dFe {_fmt(solar.get('dFe'), '+.3f')})",
        "regime_verdict": "",
        "two_engine_wiring_status": _wiring_for(wiring, "Fe", "II"),
        **dict(zip(("chosen_engine", "selection_reason"),
                   _engine_cells(two_engine, wiring_rows, "Fe", "II"))),
        "models_tried": mal.render_cell("Fe"),
        # Fe carries no registry row; this resolves to "" and says so honestly rather
        # than being special-cased blank.
        "refinement_debt": _debt_cell(debt, "Fe"),
        "engine_reach": _engine_reach_cell(_engine_reach_table(), "Fe", "II"),
        **{k: str(ed[k]) for k in (
            "engine_a_model_vintage", "engine_b_model_vintage", "classification",
            "action_needed", "source_tickets", "editorial_updated", "notes")},
    }


# ─────────────────────────────────────────────────────────────────────────────
#  rendering
# ─────────────────────────────────────────────────────────────────────────────
def render(rows: list[dict], prov: dict, phase_c_summary: dict) -> str:
    counts = phase_c_summary.get("counts", {})
    counts_txt = " ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    head = f"""\
# ELEMENT STATUS TRACKER — **GENERATED from phase_c @ {prov['commit'][:12]} — do not hand-edit**
# =====================================================================================
# GENERATED by scripts/generate_element_status_tracker_rya654.py (RYA-654, RYA-436 Move B).
# Hand-editing this file is a BUILD BREAK, not a shortcut: the generator's --check mode
# compares the committed file against a fresh regeneration and exits 1 on any difference.
#
#   to change a STATUS column   -> it is derived; fix its source and re-run the generator
#   to change an ANALYST column -> edit data/audit/element_status_tracker_editorial.yaml
#   then, either way            -> python scripts/generate_element_status_tracker_rya654.py
#
# PROVENANCE OF THE STATUS COLUMNS
#   source      : {PHASE_C_REL}
#   phase_c run : {phase_c_summary.get('ticket', '?')} | star={phase_c_summary.get('star', '?')} | generated={phase_c_summary.get('generated', '?')}
#   counts      : {counts_txt} (n_elements={phase_c_summary.get('n_elements', '?')})
#   git blob    : {prov['blob']}
#   git commit  : {prov['commit']} ({prov['committed_at']})
# The generator refuses to run while that artifact is dirty in the working tree, so these
# rows always trace to the RATIFIED, COMMITTED run — never to a local Mac canary re-run.
# CAVEAT, owed: the phase_c artifact records no host/run-id, so "ran on Sirius" is
# enforced upstream (RYA-567) and pinned here only by commit identity, not asserted.
#
# COLUMN SOURCES (one source per column)
#   element verdict verdict_value sigma n_lines delta_vs_asplund method
#                                        <- phase_c (RATIFIED CANONICAL for status, RYA-654 §1)
#   ion regime_verdict                   <- config/physics_regime_rya400.yaml (RYA-400)
#   tier                                 <- RYA-522 ratified freeze tiers, read from
#                                           scripts/build_solar_reference_v2_rya522.confidence_of
#   refinement_debt                      <- data/audit/element_refinement_registry.csv
#                                           (RYA-676) — which owed rows have a resolving
#                                           ticket, and which have NONE. A cell reading
#                                           "TBD - no resolving ticket" means a ticket is
#                                           owed; an EMPTY cell means no known refinement
#                                           path, which is not the same as "nothing owed".
#   engine_reach                         <- data/catalog/engine_coverage.csv (RYA-776) —
#                                           a JOIN, not a copy. `A:VIS · B:VIS,red-optical?`
#                                           reads "Engine A serves the optical; Engine B
#                                           serves the optical and REACHES the red without
#                                           an extract (`?`)". The keyed per-band table is
#                                           the sibling file; read it via
#                                           pipeline.coverage.engine_reach(). A cell
#                                           reading "(not generated)" means the reference
#                                           has not been built on Sirius — it is NOT a
#                                           statement that an engine reaches nothing.
#   engine_a/b_model_vintage classification action_needed source_tickets
#   editorial_updated notes              <- data/audit/element_status_tracker_editorial.yaml (hand)
#
# `verdict` (live status) and `tier` (freeze decision) are DIFFERENT THINGS and now have
# separate columns. The frozen gold reference is deliberately NOT a source here: it owns
# frozen VALUES, not status (RYA-654 §1), and gold v2 currently lags the live channel on
# Co and N pending the RYA-527 v3 re-freeze.
#
# The Fe II row is the ionization ARBITER (diagnostic, asserts no verdict), emitted from
# config.constants.FE_IONIZATION_SYNTH_ARBITER. Cross-artifact consistency is enforced by
# `python -m pipeline.ledger_consistency_guard` (RYA-632).
# Read with pandas.read_csv(path, comment='#').
"""
    buf = io.StringIO()
    buf.write(head)
    writer = csv.DictWriter(buf, fieldnames=COLUMNS, lineterminator="\n",
                            quoting=csv.QUOTE_MINIMAL)
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    return buf.getvalue()


def generate() -> str:
    prov = phase_c_provenance()
    phase_c = load_phase_c()
    return render(build_rows(), prov, phase_c.get("summary", {}))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--check", action="store_true",
                    help="verify the committed tracker equals a fresh regeneration; "
                         "exit 1 on any difference (hand-edit detector)")
    args = ap.parse_args(argv)

    try:
        fresh = generate()
    except SourceError as exc:
        print(f"TRACKER GENERATOR FAILED: {exc}", file=sys.stderr)
        return 1

    if args.check:
        if not TRACKER_PATH.exists():
            print(f"TRACKER GENERATOR --check FAILED: {TRACKER_REL} does not exist",
                  file=sys.stderr)
            return 1
        current = TRACKER_PATH.read_text(encoding="utf-8")
        if current != fresh:
            import difflib
            diff = list(difflib.unified_diff(
                current.splitlines(), fresh.splitlines(),
                fromfile=f"committed {TRACKER_REL}", tofile="regenerated", lineterm="", n=1))
            print("TRACKER GENERATOR --check FAILED: the committed tracker is NOT what the "
                  "generator produces. It was hand-edited, or a source moved and nobody "
                  "re-ran the generator.", file=sys.stderr)
            for line in diff[:60]:
                print(f"  {line}", file=sys.stderr)
            if len(diff) > 60:
                print(f"  ... {len(diff) - 60} more diff line(s)", file=sys.stderr)
            print(f"\n  Fix: python {Path(__file__).relative_to(ROOT)}", file=sys.stderr)
            return 1
        print(f"Element status tracker OK: committed {TRACKER_REL} == regeneration.")
        return 0

    TRACKER_PATH.write_text(fresh, encoding="utf-8")
    n = sum(1 for line in fresh.splitlines()
            if line and not line.startswith("#")) - 1        # minus the column header
    print(f"Wrote {TRACKER_REL} ({n} data rows) from {PHASE_C_REL}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

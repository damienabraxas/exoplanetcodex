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
    "two_engine_wiring_status",
    "engine_a_model_vintage", "engine_b_model_vintage", "classification",
    "action_needed", "source_tickets", "editorial_updated", "notes",
]

GENERATED_COLUMNS = frozenset({
    "element", "ion", "verdict", "verdict_value", "sigma", "n_lines",
    "delta_vs_asplund", "tier", "method", "regime_verdict",
    "two_engine_wiring_status",
})

#: RYA-673. The tracker already says what each element's verdict IS; it never said how
#: many engines that verdict rests on. Those are different facts, and the second one is
#: what Beta's "best of abilities on all engines" bar is actually about — an element
#: reporting a clean PASS on one engine with no cross-check looked identical here to one
#: confirmed on both.
WIRING_REL = "data/audit/two_engine_wiring_audit.csv"
WIRING_PATH = ROOT / WIRING_REL


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


def build_rows() -> list[dict]:
    phase_c = load_phase_c()
    regime = load_regime()
    editorial = load_editorial()
    wiring = load_wiring()

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
            **{k: str(ed[k]) for k in (
                "engine_a_model_vintage", "engine_b_model_vintage", "classification",
                "action_needed", "source_tickets", "editorial_updated", "notes")},
        })

    rows.append(_fe_ii_row(editorial, wiring))
    return rows


def _fe_ii_row(editorial: dict, wiring: dict) -> dict:
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

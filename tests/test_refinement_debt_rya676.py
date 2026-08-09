"""
tests/test_refinement_debt_rya676.py
====================================
RYA-676 — the refinement-debt architecture.

The defect class: RYA-524's audit spawned children, the architectural ones were
executed, and the science-refinement ones (RYA-581/585/565) sat in Backlog through
eight architecture tickets. Nothing was broken; no surface carried the sentence
"this owed row has a ticket, and it has not fired."

What these tests defend, in order of what would actually rot first:

  * the registry's ADMISSION RULE has teeth — a row without a provenance ticket
    is speculation, and speculation is what makes a debt registry useless;
  * the join is by ELEMENT and cannot silently drop a row whose ion differs from
    the ion the tracker reports (RYA-475 is about Y I; Y is reported as Y II);
  * a registry row naming an element the tracker does not carry is a build break,
    not an invisible row;
  * the report is INFORMATIONAL BY CONSTRUCTION — exit 0 regardless of content —
    and cannot change the ledger guard's exit code either;
  * `--phase-close` genuinely escalates, so the escalation is real and not a
    comment;
  * the two open-debt buckets stay SEPARATE: "has a ticket, unfired" and "has no
    ticket" are different problems and folding them would make the first
    un-actionable.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import ledger_consistency_guard as guard   # noqa: E402
from pipeline import refinement_debt_join as rdj         # noqa: E402
from pipeline import state_surfaces as ss                # noqa: E402

TRACKER = ROOT / ss.TRACKER
REGISTRY = ROOT / ss.REFINEMENT_REGISTRY


def _tracker():
    return pd.read_csv(TRACKER, comment="#")


def _row(**over):
    base = dict(element="Ba", ion="II", situation="engine-B-unwired",
                resolving_ticket="RYA-680", ticket_state="Backlog",
                provenance_ticket="RYA-673", phase="Solar Beta",
                short_label="Ba engine-B harness", notes="n")
    base.update(over)
    return rdj.DebtRow(**base)


# ── the registry itself ──────────────────────────────────────────────────────
def test_registry_is_well_formed_and_populated():
    rows = rdj.load_registry()
    assert rows, "the registry is the SSOT; an empty one is not a valid state"
    assert set(rdj.REQUIRED_COLUMNS) <= set(pd.read_csv(REGISTRY, comment="#").columns)


def test_every_row_cites_the_ticket_that_established_the_debt(tmp_path):
    """The admission rule, enforced. Without this the registry drifts into a wishlist."""
    text = REGISTRY.read_text(encoding="utf-8")
    # strip the provenance ticket from the FIRST data row only
    lines = text.splitlines(keepends=True)
    idx = next(i for i, ln in enumerate(lines)
               if ln.startswith("Ba,II,engine-B-unwired"))
    lines[idx] = lines[idx].replace("RYA-673", "", 1)
    bad = tmp_path / "registry.csv"
    bad.write_text("".join(lines), encoding="utf-8")

    with pytest.raises(rdj.RegistryError, match="provenance"):
        rdj.load_registry(bad)


def test_unknown_ticket_state_is_not_guessed_at():
    with pytest.raises(rdj.RegistryError, match="ticket_state"):
        _row(ticket_state="Shipped").render_class


def test_missing_registry_raises_rather_than_rendering_everything_clean(tmp_path):
    """A blank column reads as 'nothing owed' on exactly the rows this exists to show."""
    with pytest.raises(rdj.RegistryError, match="not found"):
        rdj.load_registry(tmp_path / "absent.csv")


# ── rendering (RYA-676 §2B value forms) ──────────────────────────────────────
@pytest.mark.parametrize("state,expected", [
    ("Backlog", "Backlog: RYA-680 (Ba engine-B harness)"),
    ("Todo", "Backlog: RYA-680 (Ba engine-B harness)"),
    ("In Progress", "In Progress: RYA-680 (Ba engine-B harness)"),
    ("Done", "Done: RYA-680 (Ba engine-B harness)"),
])
def test_value_forms(state, expected):
    assert _row(ticket_state=state).render() == expected


def test_unticketed_row_renders_the_fixed_literal():
    r = _row(resolving_ticket="TBD", ticket_state="none",
             short_label="second clean Ba II line")
    assert r.render().startswith(rdj.TBD_TEXT)
    assert r.is_unticketed


def test_an_element_with_two_debts_renders_both_most_actionable_first():
    """Ba is the worked example: an unwired Engine B and an unconfirmed single line
    are different fixes, and neither substitutes for the other.

    The expected ORDER moved when RYA-680 landed (RYA-705): Ba's engine-B row is now
    `Done`, and `Done` sits last in `_RENDER_ORDER` precisely because a discharged debt
    is the least actionable thing in the cell. So the unticketed row leads. Both rows
    are still present and still distinct, which is what this test is for -- discharging
    one debt must not swallow the other."""
    cell = rdj.debt_cell("Ba", rdj.by_element())
    assert cell.startswith(rdj.TBD_TEXT), cell
    assert "Done: RYA-680" in cell, cell
    assert cell.index(rdj.TBD_TEXT) < cell.index("Done: RYA-680"), cell


# ── the join ─────────────────────────────────────────────────────────────────
def test_join_is_by_element_so_an_ion_mismatch_cannot_hide_a_row():
    """RYA-475 is about Y I lines; the Codex reports Y as Y II (RYA-683). Joining on
    the (element, ion) pair would drop it silently -- the exact failure this ticket
    is about, one level up."""
    rows = rdj.load_registry()
    y_ions = {r.ion for r in rows if r.element == "Y"}
    assert y_ions == {"I", "II"}, "test premise moved; re-read the registry"

    tracker = _tracker()
    reported_ion = tracker.loc[tracker.element == "Y", "ion"].iloc[0]
    assert reported_ion == "II"
    assert "RYA-475" in tracker.loc[tracker.element == "Y", "refinement_debt"].iloc[0]


def test_registry_element_absent_from_the_tracker_is_a_build_break():
    with pytest.raises(rdj.RegistryError, match="joins to nothing|does not carry"):
        rdj.assert_elements_known({"Fe"})


def test_tracker_carries_the_generated_column_and_it_is_not_hand_edited():
    df = _tracker()
    assert "refinement_debt" in df.columns
    r = subprocess.run(
        [sys.executable, "scripts/generate_element_status_tracker_rya654.py", "--check"],
        cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_blank_cell_means_no_known_path_not_nothing_owed():
    """Si/Ni/Na are owed and deliberately carry no row -- no landed artifact or filed
    ticket names a specific refinement for them. If that ever silently becomes 'they
    have no debt', the registry has started lying by omission."""
    df = _tracker()
    owed_blank = df[(df.tier == "owed") & (df.refinement_debt.isna())]
    assert not owed_blank.empty, (
        "no owed element lacks a registry row any more -- good, but update this test "
        "rather than deleting it: the distinction it pins is still real")


# ── the report: informational BY CONSTRUCTION ────────────────────────────────
def test_report_always_exits_zero():
    assert rdj.report() == 0


def test_report_cannot_be_made_blocking_by_content(monkeypatch):
    """Not 'it happens to pass today' -- exit 0 must not depend on the debt count."""
    monkeypatch.setattr(rdj, "open_debt", lambda *a, **k: {
        "refinement_debt_open": [{"element": "X", "ion": "I", "tier": "owed",
                                  "situation": "s", "refinement_debt": "Backlog: RYA-1 (x)",
                                  "provenance_ticket": "RYA-0"}] * 99,
        "refinement_debt_unticketed": [],
    })
    assert rdj.report() == 0


def test_ledger_guard_json_carries_both_buckets_without_changing_its_exit_code():
    r = subprocess.run([sys.executable, "pipeline/ledger_consistency_guard.py", "--json"],
                       cwd=ROOT, capture_output=True, text=True)
    payload = json.loads(r.stdout)
    assert "refinement_debt_open" in payload
    assert "refinement_debt_unticketed" in payload
    # the exit code is the CONSISTENCY verdict alone -- debt never contributes
    assert r.returncode == (1 if payload["failures"] else 0)


def test_the_two_buckets_stay_separate():
    """'has a ticket, unfired' is fireable; 'has no ticket' is not. Folding them makes
    the first count un-actionable, which is how a report stops being read."""
    rows = rdj.load_registry()
    tiers = rdj.tracker_tiers()
    buckets = rdj.open_debt(tiers, rows)
    assert all(r["refinement_debt"].split(":")[0] in rdj._OPEN_CLASSES
               for r in buckets["refinement_debt_open"])
    assert all(r["refinement_debt"].startswith(rdj.TBD_TEXT)
               for r in buckets["refinement_debt_unticketed"])
    assert all(r["refinement_debt"].startswith("Done:")
               for r in buckets["refinement_debt_discharged"])
    seen = [set(map(str, buckets[k])) for k in
            ("refinement_debt_open", "refinement_debt_unticketed",
             "refinement_debt_discharged")]
    assert not (seen[0] & seen[1]) and not (seen[0] & seen[2]) and not (seen[1] & seen[2])


def test_in_progress_debt_is_OPEN_debt(monkeypatch):
    """RYA-705: the selector was `rendered.startswith("Backlog:")`, so an `In Progress`
    row landed in NO bucket -- invisible to the report and to the phase-close gate. That
    inverts the intent: `_RENDER_ORDER` ranks `In Progress` FIRST as the most actionable
    class. A phase could have closed over debt that was mid-flight. Ca is the live case
    (RYA-561, promotion ratified but not applied)."""
    row = _row(ticket_state="In Progress", resolving_ticket="RYA-561")
    buckets = rdj.open_debt({row.element: "owed"}, [row], row.phase)
    assert len(buckets["refinement_debt_open"]) == 1, buckets
    assert not buckets["refinement_debt_discharged"]


def test_a_discharged_row_is_reported_even_though_it_is_not_counted(capsys):
    """RYA-705: `Done` rows used to fall out of every bucket and vanish. A debt leaving
    the report in silence is how a debt gets lost -- 'the ticket closed' is a claim that
    should be auditable against the element's actual state."""
    rdj.report()
    out = capsys.readouterr().out
    assert "discharged (ticket Done, kept for provenance)" in out
    assert "NOT counted as debt" in out
    assert "Done: RYA-680" in out


def test_phase_scoping_excludes_a_deferred_row():
    """Sr II's INASAN pull is explicitly deferred post-Beta by RYA-672 §2, so it must
    not count against the Solar Beta gate."""
    rows = rdj.load_registry()
    sr = [r for r in rows if r.element == "Sr"]
    assert sr and all(r.phase == "Post-Beta" for r in sr)
    tiers = rdj.tracker_tiers()
    assert not any(r["element"] == "Sr"
                   for r in rdj.open_debt(tiers, rows)["refinement_debt_unticketed"])


# ── the escalation is real ───────────────────────────────────────────────────
def test_phase_close_escalates_to_a_failure():
    r = subprocess.run(
        [sys.executable, "-m", "pipeline.refinement_debt_join", "--phase-close"],
        cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 1, "phase-close must not pass while debt is open"
    assert "PHASE-CLOSE ESCALATION" in r.stderr


def test_plain_report_invocation_stays_green():
    r = subprocess.run(
        [sys.executable, "-m", "pipeline.refinement_debt_join", "--report"],
        cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


# ── read-set membership (RYA-676 §2C) ────────────────────────────────────────
def test_disposition_report_is_a_read_set_member_and_named_in_ledgers():
    assert ss.DISPOSITION_REPORT_MD in ss.LEDGER_PATHS
    assert ss.DISPOSITION_REPORT_MD in (ROOT / ss.LEDGERS_INDEX).read_text(encoding="utf-8")


def test_registry_is_a_state_surface_but_not_a_read_set_member():
    """You read the debt through the tracker column it generates, not the file."""
    assert ss.REFINEMENT_REGISTRY in ss.SURFACE_PATHS
    assert ss.REFINEMENT_REGISTRY not in ss.LEDGER_PATHS

"""
tests/test_harness_provenance_rya873.py — RYA-873
=================================================
`error_budget.harness_term` printed a FIXED sentence beside whatever number it was handed:

    f"{handler} MEASURED against the known optical answer, not assumed zero"

For `ProfileFitHandler` and its 0.0129 that is true. For `SynthesisHandler` and its 0.0000
it is false in the most specific way available — the sentence explicitly disclaims the
thing the number is. Every synthesis-route budget in the repo carried that contradiction.

WHY THE NUMBER DID NOT CHANGE
-----------------------------
The obvious repair is to charge the handler's banked control offset, 0.0100. RYA-873
tested the one available explanation for why that control is untrustworthy — its EW angle
disagrees by +0.1562 dex while its abundance angle agrees to −0.0100, which the control
itself calls "compensating errors ... it looks like success" — and **REFUTED** it
(`data/results/rya873/`: 15 of 18 lines carry no predicted contamination and still scatter
0.50–2.07). So 0.0100 is the abundance angle of an unexplained compensating-error pass and
charging it would launder "looks like success" into a published systematic.

Neither number is established, so **nothing is charged and the budget says so**. RYA-875
is the RCA that resolves the split. These tests pin the honesty property, not a value:

1. an UNCHARGED residual is never described as measured;
2. a MEASURED one still is, so the fix is not a blanket hedge;
3. an unstated provenance says it is unstated rather than inheriting the claim;
4. the uncharged prose NAMES the reason and the ticket, so the absence is auditable;
5. provenance travels with the value and the label — one decision, one mapping;
6. the divergence declaration still holds, and still points at a real artifact.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import harness_residual as hr                        # noqa: E402
from pipeline.error_budget import build, harness_term              # noqa: E402


#: RYA-875 — `SynthesisHandler` left `UNCHARGED_CONTROL_RESIDUAL_DEX` because its numbers
#: agreed, so no real handler is uncharged today. The MECHANISM still has to be tested, and
#: deleting these tests because the fixture went away is how a guard quietly stops guarding.
#: So the uncharged path is exercised against an INJECTED declaration instead.
_FAKE = "AHandlerNobodyControlled"


@pytest.fixture
def uncharged(monkeypatch):
    monkeypatch.setitem(hr.HANDLER_RESIDUAL_DEX, _FAKE, 0.0)
    monkeypatch.setitem(hr.UNCHARGED_CONTROL_RESIDUAL_DEX, _FAKE, {
        "control_dex": 0.0250,
        "control_artifact": "data/audit/synthesis_control/control_FeI.json",
        "ticket": "RYA-999",
        "why_not_charged": "a test fixture standing in for a real divergence",
    })
    return hr.for_handler(_FAKE)


def _harness_line(**kw) -> str:
    b = build("Fe", 5000.0, 40, scatter_dex=0.10, gf_graded=False, **kw)
    lines = [x for x in b.describe().splitlines() if "harness residual" in x]
    assert len(lines) == 1
    return lines[0]


# ── 1/2. the claim is made only where it is true ─────────────────────────────────────

def test_an_uncharged_residual_is_never_described_as_measured(uncharged):
    """THE DEFECT, in one assertion. RED on the pre-RYA-873 code by construction: the
    sentence was a constant, so this line read 'MEASURED ... not assumed zero' beside a
    zero that was exactly that."""
    line = _harness_line(**uncharged.budget_kwargs())
    assert "MEASURED" not in line, line
    assert "NO RESIDUAL IS CHARGED" in line
    assert "ABSENCE, not a measured zero" in line


def test_a_measured_residual_still_says_measured():
    """VERIFY THE FIX DISCRIMINATES (RYA-805/845). A repair that stopped claiming
    'measured' everywhere would pass test 1 while destroying the information — the
    profile fitter's 0.0129 IS measured and must keep saying so."""
    line = _harness_line(**hr.for_handler("ProfileFitHandler").budget_kwargs())
    assert "MEASURED against the known optical answer" in line
    assert "NO RESIDUAL IS CHARGED" not in line


def test_an_unstated_provenance_does_not_inherit_the_claim():
    """A caller that says nothing gets a term that says nothing — not the old default,
    which was to assert measurement on the caller's behalf."""
    line = _harness_line(harness_residual_dex=0.0129, handler="ProfileFitHandler")
    assert "NOT STATED" in line
    assert "MEASURED against" not in line


# ── 3. the absence is auditable, not merely stated ───────────────────────────────────

def test_the_uncharged_prose_names_its_reason_its_ticket_and_its_artifact(uncharged):
    """An absence with no reason is indistinguishable from an omission (RYA-833). The
    budget is an artifact a reader may see without the code, so the reason travels IN it.
    """
    line = _harness_line(**uncharged.budget_kwargs())
    note = hr.UNCHARGED_CONTROL_RESIDUAL_DEX[_FAKE]
    assert note["ticket"] in line, "the budget must name the ticket that owes the RCA"
    assert note["control_artifact"] in line, "and the artifact a reader can go check"
    assert note["why_not_charged"] in line, "and the stated reason, verbatim"


def test_the_refuted_contamination_hypothesis_is_still_recorded():
    """RYA-873's content is a NEGATIVE result, and a negative result nobody wrote down
    gets re-derived. It no longer belongs in `harness_residual.py` — RYA-875 showed the
    number it was doubting (−0.0100) was a line-set artifact, so the contamination
    refutation is now a finding about the control's EW ANGLE, which is still unexplained.
    It has to survive somewhere checkable, and that somewhere is its own artifact."""
    import json
    p = ROOT / "data" / "results" / "rya873" / "rya873_crowding_summary.json"
    assert p.exists(), "the refutation artifact is gone"
    bulk = json.loads(p.read_text())["bulk_test"]
    assert bulk["verdict"].startswith("REFUTED")
    assert bulk["n_bulk"] == 15 and bulk["bulk_observed_ratio_max"] > 2.0


# ── 4. provenance travels with the pair ──────────────────────────────────────────────

def test_provenance_travels_with_the_value_and_the_label():
    for handler in hr.HANDLER_RESIDUAL_DEX:
        kw = hr.for_handler(handler).budget_kwargs()
        assert set(kw) == {"harness_residual_dex", "handler", "harness_provenance"}
        assert kw["harness_provenance"] in (hr.PROV_MEASURED, hr.PROV_UNCHARGED)


def test_a_declared_divergence_is_exactly_the_uncharged_set():
    """The two must not be able to disagree: a handler is described as uncharged if and
    only if it is declared as one. Deriving the prose from anything else would let a
    budget claim a measurement the declaration denies."""
    declared = set(hr.UNCHARGED_CONTROL_RESIDUAL_DEX)
    uncharged = {h for h in hr.HANDLER_RESIDUAL_DEX
                 if hr.for_handler(h).provenance == hr.PROV_UNCHARGED}
    assert declared == uncharged, (declared, uncharged)


# ── 5. the interim is an interim ─────────────────────────────────────────────────────

def test_the_synthesis_handler_is_now_charged_a_MEASURED_zero():
    """RYA-873 landed 'charge nothing' and said this test was the one to update
    deliberately if a later ticket charged a number. RYA-875 did: the −0.0100 it was
    waiting on turned out to be a LINE-SET artifact, and paired per line the offset is
    0.0000. So the charged VALUE is unchanged and its provenance flipped from assumed to
    measured — and the declaration went with it, **because the numbers agree**."""
    assert hr.HANDLER_RESIDUAL_DEX["SynthesisHandler"] == 0.0
    assert "SynthesisHandler" not in hr.UNCHARGED_CONTROL_RESIDUAL_DEX
    assert hr.for_handler("SynthesisHandler").provenance == hr.PROV_MEASURED
    line = _harness_line(**hr.for_handler("SynthesisHandler").budget_kwargs())
    assert "MEASURED against the known optical answer" in line


def test_the_uncharged_note_is_empty_for_a_handler_that_is_charged(uncharged):
    assert hr.uncharged_note("ProfileFitHandler") == ""
    assert hr.uncharged_note("SynthesisHandler") == ""      # charged since RYA-875
    assert hr.uncharged_note(_FAKE)                          # and the machinery still works


def test_one_function_owns_the_harness_term():
    """ONE HOME for the term, which is the invariant that matters.

    The first cut of this guard grepped for the sentence "MEASURED against the known
    optical answer" and failed on two files that merely QUOTE it in a docstring — a
    guard that cannot tell a citation from a claim. What must be unique is the place the
    Term is CONSTRUCTED: a second constructor is how a caller that "knows" its number is
    measured would reintroduce the fixed sentence.
    """
    hits = [f.relative_to(ROOT)
            for f in sorted((ROOT / "pipeline").rglob("*.py"))
            + sorted((ROOT / "scripts").rglob("*.py"))
            if 'Term("harness residual"' in f.read_text()]
    assert hits == [Path("pipeline/error_budget.py")], hits


def test_the_prose_is_reachable_only_through_the_provenance_branch():
    """And the sentence itself must sit inside `harness_term`, not at module scope where
    another function could reach it."""
    import ast
    tree = ast.parse((ROOT / "pipeline" / "error_budget.py").read_text())
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "harness_term")
    claims = [n for n in ast.walk(fn) if isinstance(n, ast.Constant)
              and isinstance(n.value, str)
              and "MEASURED against the known optical answer" in n.value]
    assert len(claims) == 1, f"the measured claim appears {len(claims)}x in harness_term"


def test_the_term_value_is_unchanged_by_the_prose_fix():
    """RYA-873 moves NO number. The term's dex is whatever it was; only its description
    changed. Asserted so a future edit to the prose cannot quietly become an edit to the
    bar."""
    for handler, dex in hr.HANDLER_RESIDUAL_DEX.items():
        t = harness_term(dex, handler, hr.for_handler(handler).provenance)
        assert t.dex == dex
        assert t.name == "harness residual"
        assert t.averages_down is False

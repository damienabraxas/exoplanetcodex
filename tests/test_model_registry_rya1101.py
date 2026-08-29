"""
tests/test_model_registry_rya1101.py — RYA-1101
===============================================
The 8-model roster used to live only in chat memory — the nicknames "Frankenstein"
(models 5+6) and "The Bride" (model 8) were in no ledger, so the roster got re-litigated
every session. This pins the registry AND the property that makes it trustworthy: **every
live token is read from the code, so the file cannot launder a value supplied from memory.**
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import band_products, model_registry  # noqa: E402
from pipeline.model_registry import ModelRegistryError  # noqa: E402

NOTE = ROOT / "docs" / "catalog" / "model_registry_notes.md"


@pytest.fixture(scope="module")
def rows():
    return model_registry.load()


def test_every_live_token_exists_verbatim_in_the_code(rows):
    """🔴 THE ANTI-MEMORY GUARD. A token in the registry that is not in
    band_products.TREATMENTS means it was typed from memory rather than read from source —
    which is the exact failure this registry exists to end."""
    code = set(band_products.TREATMENTS)
    for r in rows:
        if r["status"] != "live":
            continue
        assert r["stored_token"], f"model {r['model_id']} is live with no token"
        assert r["stored_token"] in code, (
            f"model {r['model_id']} token {r['stored_token']!r} is not in the code")


def test_the_in_dev_model_has_no_invented_token(rows):
    """Model 8 (VoronoiRT / 'The Bride') has no token in the code yet. The registry leaves
    it BLANK — an absence recorded is not the same as a name made up."""
    bride = next(r for r in rows if r["model_id"] == "8")
    assert bride["stored_token"] == ""
    assert bride["status"] == "in-dev"
    assert bride["nickname"] == "The Bride"


def test_frankensteins_dog_is_a_row_with_no_invented_token(rows):
    """RYA-1115: Ryan ratified the Dog as a DISTINCT model, so it is in the roster. But
    `treatment_axes.AXIS_NATIVE` registers no 1D-LTE leg on a mean-3D reference and
    `ATMOSPHERES` is a closed set with no such atmosphere, so every code-sourced field is
    left BLANK — the model 8 rule applied to model 9."""
    dog = next(r for r in rows if r["model_id"] == "9")
    assert dog["nickname"] == "Frankenstein's Dog"
    assert dog["stored_token"] == ""
    assert dog["atmosphere"] == ""
    assert dog["status"] == "not-emitted"


def test_the_dog_is_not_in_the_frankenstein_pair(rows):
    """It would be the ATMOSPHERE-rung comparand, not a member of the mandatory 5/6 NLTE
    pair. Putting it in the pair group would let a consumer difference it as if it were."""
    dog = next(r for r in rows if r["model_id"] == "9")
    assert dog["pair_group"] != "frankenstein"


def test_not_emitted_is_a_distinct_status_from_in_dev(rows):
    """The Bride is being BUILT; the Dog is a leg nothing emits. One token for two
    situations is the defect this whole registry exists to end."""
    assert "not-emitted" in model_registry.STATUSES
    assert {r["status"] for r in rows} <= set(model_registry.STATUSES)
    assert next(r for r in rows if r["model_id"] == "8")["status"] == "in-dev"
    assert next(r for r in rows if r["model_id"] == "9")["status"] == "not-emitted"


def test_the_line_set_hook_is_open_with_a_closed_vocabulary(rows):
    """RYA-1115 opens the column for RYA-1111. It is model-scoped nowhere yet, so every row
    reads '-', but the vocabulary exists NOW so the first real value cannot be typed from
    memory."""
    assert "line_set" in rows[0]
    assert {r["line_set"] for r in rows} == {"-"}
    assert model_registry.check_line_sets(rows) == []
    # ⚠️ `asplund`, NOT `asplund-graded`. RYA-1101 opened this column before the axis had
    # an owner and I guessed the spelling; RYA-1111 owns the axis and names the values
    # {asplund, gbs, our-graded, our-deep-graded}. Corrected there, asserted here so the
    # old guess cannot come back.
    for pool in ("asplund", "gbs", "our-graded", "our-deep-graded"):
        assert pool in model_registry.LINE_SETS
    assert "asplund-graded" not in model_registry.LINE_SETS


def test_the_line_set_guard_refuses_a_value_outside_the_vocabulary():
    """The control: the guard must actually FIRE. `consistent` is the realistic mistake —
    RYA-1105 retired the tier but `--lines-tier consistent` is still live in the code."""
    bad = model_registry.check_line_sets([{"model_id": "1", "line_set": "consistent"}])
    assert len(bad) == 1 and "consistent" in bad[0]


def test_bare_ENGINE_B_refuses_to_bind(rows):
    """🔴 THE COLLISION. `ENGINE-B` resolves to model 1 today but historically meant model 4
    (now ENGINE-B-NLTE) — same string, two physics, and a pre-RYA-906 artifact does not say
    which. It must raise, naming both, not pick one."""
    with pytest.raises(ModelRegistryError) as e:
        model_registry.resolve("ENGINE-B")
    msg = str(e.value)
    assert "1D-LTE" in msg and "ENGINE-B-NLTE" in msg, "the error must name BOTH candidates"


def test_the_collision_is_documented_not_just_guarded():
    """A guard nobody can read is a trap. The hazard, both candidates, and the fact that the
    canonical binding is Ryan's call all have to be written down."""
    txt = NOTE.read_text()
    assert "ENGINE-B" in txt and "ENGINE-B-NLTE" in txt
    assert "Ryan" in txt, "the binding decision must be recorded as unresolved, not settled"


def test_the_canonical_tokens_still_resolve(rows):
    """The guard must not be so broad that it refuses the real thing."""
    assert model_registry.resolve("1D-LTE")["model_id"] == "1"
    assert model_registry.resolve("ENGINE-B-NLTE")["model_id"] == "4"


def test_an_unknown_token_refuses_rather_than_returning_none():
    with pytest.raises(ModelRegistryError):
        model_registry.resolve("ENGINE-Z")


def test_the_frankenstein_pair_is_marked_on_both_rows(rows):
    """RYA-542: the NLTE effect is (6 − 5) on ONE atmosphere. If either row could be read as
    a standalone product, someone would difference it against 1D-LTE and report the
    1D→mean-3D ATMOSPHERE shift as non-LTE physics."""
    frank = {r["model_id"] for r in rows if r["pair_group"] == "frankenstein"}
    assert frank == {"5", "6"}
    for r in rows:
        if r["model_id"] in ("5", "6"):
            assert r["nickname"] == "Frankenstein"
            assert "pair" in r["notes"].lower()


def test_no_live_count_column_exists(rows):
    """Live counts are volatile and owned by RYA-1015. A copy here would be a second source
    of truth for a number that changes every run."""
    for col in rows[0]:
        assert not any(h in col.lower() for h in ("count", "n_products", "n_lines"))


def test_the_mean3d_rows_carry_the_axes_the_CODE_stores(rows):
    """Not the ticket's candidate spelling. The candidate said scale 'mean3D-LTE'; the code
    registers '<3D>-LTE'. The code wins, and this pins that it did."""
    from pipeline import treatment_axes as ta
    for r in rows:
        tok = r["stored_token"]
        if tok in ta.AXIS_NATIVE:
            ax = ta.AXIS_NATIVE[tok]
            assert r["scale"] == ax.scale
            assert r["atmosphere"] == ax.atmos
            assert r["route"] == ax.route


def test_every_row_names_a_source_ticket(rows):
    """No value from memory — each non-code field traces to a decision ticket."""
    for r in rows:
        assert r["source_ticket"].strip() and "RYA-" in r["source_ticket"]


def test_the_verifier_passes_and_reports_the_out_of_roster_tokens():
    """The smoke test the ticket asks for, run as a test. It must also SURFACE the live code
    tokens that are not models here (ENGINE-B, 1D-LTE-LABGF) rather than dropping them."""
    r = subprocess.run([sys.executable, "scripts/verify_model_registry.py"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert ("MODEL_REGISTRY OK: 9 models (7 live, 1 in-dev, 1 not-emitted)"
            in r.stdout)
    assert "1D-LTE-LABGF" in r.stdout and "ENGINE-B" in r.stdout


def test_the_registry_is_in_the_ledger_read_set():
    """RYA-659: if it is not in LEDGERS.md it is not read at session start, and the roster
    goes back to living in chat memory."""
    assert "data/catalog/model_registry.csv" in (ROOT / "LEDGERS.md").read_text()
    reg = (ROOT / "CODEX_STATE_REGISTER.md").read_text()
    assert "model_registry.csv" in reg
    assert "3D-NLTE Amarsi" not in reg, "the register must POINT at the roster, not copy it"

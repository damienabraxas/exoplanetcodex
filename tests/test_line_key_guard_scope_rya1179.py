"""RYA-1179 — the wavelength-only-join guard must scope its EP test to the DECISION.

RYA-1142 check A2c: `_enclosing_has_ep()` tested the whole enclosing FunctionDef, so a
single `ep` bound anywhere laundered every wavelength-only comparison in that function.
The guard reported clean on the exact defect it exists to catch, and never looked at the
RYA-1136 CNO intake or the RYA-1132 Al intake at all.

`test_an_unrelated_ep_in_the_same_function_no_longer_launders` is the QA's positive
control, promoted to a regression test.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import audit_line_keys_rya1037 as A  # noqa: E402


def kinds(src: str) -> list[str]:
    """Run the shipped guard over a source string; return the finding kinds."""
    tree = ast.parse(src)
    v = A.Visitor("t.py", src, A._key_positions(tree), tree)
    v.visit(tree)
    return [f.kind for f in v.found]


# ── the three the spec asks for ──────────────────────────────────────────────
def test_a_wavelength_only_comparison_is_flagged():
    assert "WAVE_ONLY_TOL" in kinds(
        "def f(a, b, tol):\n"
        "    return abs(a.wavelength_air_A - b.wavelength_air_A) < tol\n")


def test_an_unrelated_ep_in_the_same_function_no_longer_launders():
    """🔴 THE REGRESSION. Identical comparison, plus an `ep` bound earlier for an
    unrelated purpose. The old function-scoped test went silent here; that silence is
    what hid the CNO and Al intakes."""
    src = ("def f(a, b, tol, rows):\n"
           "    ep = rows[0].excitation_potential_eV\n"
           "    note = f'ep={ep}'\n"
           "    return abs(a.wavelength_air_A - b.wavelength_air_A) < tol, note\n")
    assert "WAVE_ONLY_TOL" in kinds(src), (
        "an unrelated `ep` in the enclosing function laundered a wavelength-only "
        "comparison — RYA-1179 has regressed")


def test_a_genuine_dual_key_conjunction_passes():
    assert "WAVE_ONLY_TOL" not in kinds(
        "def f(a, b, tol, eptol):\n"
        "    return (abs(a.wavelength_air_A - b.wavelength_air_A) < tol\n"
        "            and abs(a.excitation_potential_eV - b.excitation_potential_eV) < eptol)\n")


# ── the shapes real code actually uses ───────────────────────────────────────
def test_the_two_line_augmented_form_passes():
    """`ok = <lambda test>` / `ok &= <EP test>` is one decision written over two lines,
    and the builder's own matcher is written that way."""
    assert "WAVE_ONLY_TOL" not in kinds(
        "def f(a, b, tol, eptol):\n"
        "    ok = abs(a.wavelength_air_A - b.wavelength_air_A) < tol\n"
        "    ok &= abs(a.excitation_potential_eV - b.excitation_potential_eV) < eptol\n"
        "    return ok\n")


def test_a_comprehension_filter_with_both_terms_passes():
    """The shape `nearest_canonical()` uses — and the reason the CNO intake is clean."""
    assert "WAVE_ONLY_TOL" not in kinds(
        "def f(rows, w, ep):\n"
        "    return [r for r in rows\n"
        "            if abs(float(r['wavelength_air_A']) - w) <= 0.03\n"
        "            and abs(float(r['excitation_potential_eV']) - ep) <= 0.002]\n")


def test_a_mention_is_not_a_constraint():
    """RYA-1141's rule. Naming `ep` while WRITING it into an output must not license a
    wavelength-only filter — and the substring-tolerant EP pattern would otherwise be
    satisfied by a variable called `rows_ep`."""
    assert "WAVE_ONLY_TOL" in kinds(
        "def f(rows, w, tol):\n"
        "    rows_ep = [r for r in rows if abs(r['wavelength_air_A'] - w) < tol]\n"
        "    return rows_ep\n")
    assert "WAVE_ONLY_TOL" in kinds(
        "def f(rows, w, tol, ep):\n"
        "    hit = [r for r in rows if abs(r['wavelength_air_A'] - w) < tol]\n"
        "    return {'match': hit, 'ep_reported': ep}\n")


# ── the sites the ticket says must be reachable ──────────────────────────────
def test_the_cno_n_i_join_is_dual_key_and_clean():
    """Spec item 3. 🔴 ALREADY FIXED — by RYA-1143, after this ticket was filed. The
    N I loop calls `nearest_canonical(..., published_ep, ...)` and that function gates on
    wavelength AND EP in one conjunction, so the fixed guard correctly does NOT flag it.
    Asserted here so the repair cannot silently regress."""
    src = (ROOT / "scripts/build_cno_intake_rya1136.py").read_text()
    tree = ast.parse(src)
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "nearest_canonical")
    conj = [n for n in ast.walk(fn) if isinstance(n, ast.BoolOp) and isinstance(n.op, ast.And)]
    assert any(A._mentions_wave(c) and A._mentions_ep(c) for c in conj), (
        "nearest_canonical no longer gates on wavelength AND EP in one conjunction")
    assert "WAVE_ONLY_TOL" not in kinds(src)

    census = next(n for n in ast.walk(tree)
                  if isinstance(n, ast.FunctionDef) and n.name == "atomic_census")
    assert "N_I_ADOPTED_SET" in A._expr_names(census), "the N I loop is gone"


def test_the_guard_is_no_longer_function_scoped():
    """Structural: the decision test must not consult the function stack."""
    src = (ROOT / "scripts/audit_line_keys_rya1037.py").read_text()
    fn = next(n for n in ast.walk(ast.parse(src))
              if isinstance(n, ast.FunctionDef) and n.name == "_decision_has_ep")
    used = A._expr_names(fn)
    assert "_fn_stack" not in used, "_decision_has_ep still reads the enclosing function"
    assert "_statement_of" in used


def test_the_repo_scan_stays_green_and_every_waiver_is_owned():
    """The 19 sites the hole was hiding are waived, not suppressed: each needs a ticket
    and a reason, and the count refuses a twentieth."""
    import yaml
    doc = yaml.safe_load((ROOT / "config/line_key_waivers.yaml").read_text())
    new = [w for w in doc["waivers"] if w.get("ticket") == "RYA-1179"]
    assert sum(w["count"] for w in new) == 19, "the RYA-1179 waiver set changed size"
    for w in new:
        assert w.get("reason", "").strip(), w
        assert "[JOIN]" in w["reason"] or "[NON-JOIN]" in w["reason"], (
            f"{w['file']}: every RYA-1179 waiver must classify itself for triage")
    joins = sum(w["count"] for w in new if "[JOIN]" in w["reason"])
    assert joins == 11, f"the genuine-join count moved: {joins}"

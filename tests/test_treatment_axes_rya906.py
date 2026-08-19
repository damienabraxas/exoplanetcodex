"""RYA-906 — the agree-guard, landing WITH the mapping module rather than after it.

This is the guard that would have caught RYA-869. That defect was a string compare
(`treatment == "ENGINE-B"`) against a label set that had silently grown a variant, and it
reached publication because nothing asserted that the label and the artifact still agreed
about what engine a product was.

So: for every committed cell that carries both a legacy label and a route witness, the
axes derived here must agree with the witness — and the ONE known disagreement (the
near-UV Fe I cell, whose `1D-LTE` label is simply false about its route) is asserted as an
expected exception rather than papered over.
"""
from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import pytest

from pipeline.treatment_axes import (
    Axes, LEGACY, ROUTE_BY_HANDLER, UnknownTreatment, axes_for, display_for, resolve_route)

ROOT = Path(__file__).resolve().parents[1]


def _committed_csvs():
    out = subprocess.run(["git", "ls-files", "*.csv"], cwd=ROOT,
                         capture_output=True, text=True).stdout.split()
    return [ROOT / p for p in out]


def _treatment_rows():
    """Every committed row carrying a treatment, with whatever route evidence it has."""
    for f in _committed_csvs():
        try:
            head = f.open(encoding="utf-8", errors="replace").readline()
        except OSError:
            continue
        if "treatment" not in head:
            continue
        try:
            rows = list(csv.DictReader(f.open(encoding="utf-8", errors="replace")))
        except Exception:
            continue
        for r in rows:
            t = (r.get("treatment") or "").strip()
            if t:
                yield f, t, (r.get("handler") or "").strip(), r.get("ew_inversion")


# ── the vocabulary is closed, and closing it is the point ────────────────────

def test_every_committed_treatment_label_is_known_to_the_module():
    """An unknown label must be a HARD failure, not a default.

    RYA-869 published four wrong systematics because a new variant met an `==` that did
    not know about it and quietly did nothing. The module raises instead.
    """
    seen = {t for _, t, _, _ in _treatment_rows()}
    assert seen, "no committed treatment rows found — this guard would pass vacuously"
    unknown = seen - set(LEGACY)
    assert not unknown, (
        f"treatment labels on disk that pipeline.treatment_axes has never heard of: "
        f"{sorted(unknown)}. Add them to LEGACY and _ROUTE_BY_LABEL.")


def test_an_unknown_label_raises_rather_than_defaulting():
    """Positive control for the guard above: the module must actually refuse."""
    with pytest.raises(UnknownTreatment):
        axes_for("ENGINE-C-SUPER-NLTE")


# ── the agree-guard proper ───────────────────────────────────────────────────

def test_derived_route_agrees_with_the_handler_on_every_committed_cell():
    """The handler is the authoritative witness. Where a row has one, the axes must match."""
    checked = 0
    for f, t, handler, _ in _treatment_rows():
        if not handler:
            continue
        expected = ROUTE_BY_HANDLER.get(handler)
        if expected is None:
            pytest.fail(f"{f}: handler {handler!r} is not in ROUTE_BY_HANDLER — a new "
                        f"route vocabulary appeared and nothing mapped it")
        assert axes_for(t, handler=handler).route == expected, f"{f}: {t} x {handler}"
        checked += 1
    assert checked, "no row carried a handler — the agree-guard would pass vacuously"


def test_the_engine_labels_really_do_pin_one_route_each():
    """`_ROUTE_BY_LABEL` claims ENGINE-A* is EW and ENGINE-B* is synthesis. RE-MEASURE it
    against the tree rather than trusting the comment that says so — that claim is the
    only place the module still reads a route off a label, so it has to keep earning it."""
    pairs: dict[str, set[str]] = {}
    for _, t, handler, _ in _treatment_rows():
        if handler and t.startswith("ENGINE-"):
            pairs.setdefault(t, set()).add(ROUTE_BY_HANDLER.get(handler, handler))
    assert pairs, "no ENGINE-* row carried a handler — nothing was actually verified"
    for label, routes in pairs.items():
        assert len(routes) == 1, (
            f"{label} appears with MORE THAN ONE route {sorted(routes)} — it no longer "
            f"pins a route, so _ROUTE_BY_LABEL must stop claiming it does")


def test_the_1D_LTE_label_is_genuinely_ambiguous_and_that_is_why_route_is_not_read_from_it():
    """The exception the whole module is built around, asserted as a POSITIVE control.

    If this ever fails it means `1D-LTE` stopped being used by both routes — at which
    point the design rationale changes and someone must revisit it deliberately.
    """
    routes = set()
    for _, t, handler, _ in _treatment_rows():
        if t == "1D-LTE" and handler:
            routes.add(ROUTE_BY_HANDLER.get(handler, handler))
    assert routes == {"ew", "synth"}, (
        f"`1D-LTE` was expected to appear on BOTH routes (that is the RYA-906 premise); "
        f"found {sorted(routes)}")


# ── the canary ───────────────────────────────────────────────────────────────

def test_CANARY_near_uv_fe_i_renders_synth_not_ew():
    """🔴 THE ACCEPTANCE CANARY (RYA-906 §2).

    The near-UV Fe I cell is a RYA-759 synthesis flux fit whose legacy label is `1D-LTE`,
    identical to the VIS EW cells. If the migration ever keys on the label instead of the
    route, this cell renders `EW` and is WRONG. It is the one cell that proves the rename
    was done correctly.
    """
    assert display_for("1D-LTE", handler="SynthesisHandler") == "Synth · 1D-LTE"
    # ...and on rows written before `handler` existed, the per-line bool must carry it.
    assert display_for("1D-LTE", ew_inversion=False) == "Synth · 1D-LTE"
    # The control: the SAME label on the EW route must NOT say Synth.
    assert display_for("1D-LTE", handler="ProfileFitHandler") == "EW · 1D-LTE"


def test_CANARY_the_real_near_uv_artifact_on_disk_resolves_to_synth():
    """The canary against the committed artifact, not a hand-built argument."""
    hits = []
    for f in _committed_csvs():
        if "3000_3780" not in f.name or "products" not in f.name:
            continue
        for r in csv.DictReader(f.open(encoding="utf-8", errors="replace")):
            if (r.get("treatment") or "").strip() != "1D-LTE":
                continue
            ax = axes_for("1D-LTE", handler=(r.get("handler") or "").strip() or None)
            hits.append((f.name, ax.display, ax.route_basis))
    assert hits, "no committed near-UV Fe I 1D-LTE product row found to check"
    for name, disp, basis in hits:
        assert disp == "Synth · 1D-LTE", f"{name}: rendered {disp!r} (basis {basis})"


# ── the display name is derived, and the mapping matches the ratified table ──

@pytest.mark.parametrize("treatment,handler,expected", [
    ("1D-LTE",          "ProfileFitHandler", "EW · 1D-LTE"),
    ("1D-LTE-LABGF",    "ProfileFitHandler", "EW · 1D-LTE · lab-gf"),
    ("ENGINE-A",        None,                "EW · 1D-NLTE · Bergemann"),
    ("ENGINE-A-3DNLTE", None,                "EW · 3D-NLTE · Amarsi"),
    ("ENGINE-B",        None,                "Synth · 1D-LTE"),
    ("ENGINE-B-NLTE",   None,                "Synth · 1D-NLTE · Gerber"),
    ("1D-LTE",          "SynthesisHandler",  "Synth · 1D-LTE"),
])
def test_display_names_match_the_ratified_table(treatment, handler, expected):
    """RYA-906 §2's mapping table, asserted verbatim. The table is the contract."""
    assert display_for(treatment, handler=handler) == expected


# ── the whole point: a FIELD test replaces a spelling list ───────────────────

def test_everything_nlte_is_a_field_test_not_a_list_of_spellings():
    """`scale != '1D-LTE'` cannot be outgrown by a new variant. A spelling list can, and
    that is precisely how RYA-869 happened."""
    nlte = {t for t in LEGACY if axes_for(t).is_nlte}
    assert nlte == {"ENGINE-A", "ENGINE-A-3DNLTE", "ENGINE-B-NLTE"}
    lte = {t for t in LEGACY if not axes_for(t).is_nlte}
    assert lte == {"1D-LTE", "1D-LTE-LABGF", "ENGINE-B"}
    # 🔴 The RYA-869 shape itself: the naive compare misses the NLTE synthesis variant,
    # the field test does not. This asserts the RELATIONSHIP, not a constant (RYA-845).
    naive = {t for t in LEGACY if t == "ENGINE-B"}
    field = {t for t in LEGACY if axes_for(t).route == "synth"}
    assert "ENGINE-B-NLTE" in field and "ENGINE-B-NLTE" not in naive


def test_atmos_is_stored_not_derived_when_the_caller_knows_it():
    """Ryan ratified an independent field: a derived atmosphere is correct today and
    silently wrong the first day a model runs on a different one."""
    assert axes_for("ENGINE-B-NLTE").atmos == "marcs-ges"      # legacy inference
    assert axes_for("ENGINE-A").atmos == "atlas9"
    # a caller that KNOWS overrides the inference — that is the whole point of storing it
    assert axes_for("ENGINE-B-NLTE", atmos="atlas9").atmos == "atlas9"


def test_route_basis_records_which_witness_settled_it():
    """An inference and a reading must never look alike in an artifact."""
    assert axes_for("1D-LTE", handler="SynthesisHandler").route_basis == "handler"
    assert axes_for("1D-LTE", ew_inversion=False).route_basis == "ew_inversion"
    assert axes_for("ENGINE-B").route_basis == "label-family"
    unknown = axes_for("1D-LTE")
    assert unknown.route is None and unknown.route_basis == "unknown"
    # and an unknown route must not silently render as a real one
    assert unknown.display.startswith("route?")


def test_handler_outranks_ew_inversion_when_both_are_present():
    """Declared precedence, asserted — otherwise 'strongest witness' is just prose."""
    route, basis = resolve_route("1D-LTE", handler="SynthesisHandler", ew_inversion=True)
    assert (route, basis) == ("synth", "handler")


# ── RYA-850 reconciliation: graded-ness becomes an AXIS, not a suffix ────────

def test_the_gf_axis_and_the_LABGF_suffix_agree_on_every_legacy_label():
    """RYA-906 §6 requires `1D-LTE-LABGF` to become `gf:lab` (an axis value) WITHOUT
    breaking RYA-850's graded-primary reporting. This is the equivalence proof: over the
    whole legacy vocabulary, the axis test and the string test give the same answer, so
    migrating the call site cannot move a value."""
    from pipeline.graded_reporting import is_graded
    for t in LEGACY:
        assert is_graded(t) == is_graded(t, gf=axes_for(t).gf), t


def test_the_gf_axis_OUTRANKS_the_suffix_when_both_are_present():
    """The axis is what the emitter measured; the suffix is what someone typed. When they
    disagree the measurement wins — otherwise the fallback is not a fallback."""
    from pipeline.graded_reporting import is_graded
    assert is_graded("1D-LTE", gf="lab") is True
    assert is_graded("1D-LTE-LABGF", gf="kurucz") is False


def test_LABGF_is_not_assumed_to_be_the_EW_route():
    """🔴 A SECOND CANARY, and it corrects RYA-906's own mapping table.

    The ticket's table maps `1D-LTE-LABGF` -> `EW · 1D-LTE · lab-gf`. On the committed
    near-UV cell that is WRONG: `rya836_nearuv_lab_gf_subpool.py` builds it with
    `handler="SynthesisHandler"` and says so in as many words — "both pools are the
    RYA-759 near-UV SYNTHESIS route (flux fit, no EW anywhere)". The label carries the
    same lie as the `1D-LTE` cell, one variant over.

    So the module must render it from the handler, both ways, and never from the name.
    """
    assert display_for("1D-LTE-LABGF", handler="SynthesisHandler") == "Synth · 1D-LTE · lab-gf"
    assert display_for("1D-LTE-LABGF", handler="ProfileFitHandler") == "EW · 1D-LTE · lab-gf"


def test_LABGF_really_does_appear_on_both_routes_in_the_committed_tree():
    """Measured, not asserted from the comment above."""
    routes = set()
    for _, t, handler, _ in _treatment_rows():
        if t == "1D-LTE-LABGF" and handler:
            routes.add(ROUTE_BY_HANDLER.get(handler, handler))
    assert routes == {"ew", "synth"}, (
        f"`1D-LTE-LABGF` was expected on BOTH routes; found {sorted(routes)}")

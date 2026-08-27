"""RYA-1069: the preflight readiness conductor never defaults an unknown to READY.

Every CRITICAL failure condition in the ticket gets a POSITIVE CONTROL here -- a case
constructed so the test would FAIL if the guard were removed. A test that only exercises
the happy path cannot tell a working gate from an absent one (RYA-905's lesson, and the
reason `codex-data-audit` demands positive controls by name).

The conductor is loaded BY PATH. Importing it must not require the measurement harness,
whose Kitt Peak resolution `SystemExit`s at import wherever the atlas is not staged
(RYA-1064 friction #3) -- that is itself one of the behaviours under test.
"""
import csv
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location(
        "_preflight_readiness_test", ROOT / "scripts" / "preflight_readiness.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_preflight_readiness_test"] = mod
    spec.loader.exec_module(mod)
    return mod


pr = _load()


class _Spec:
    """A stand-in for `measure_band_ew.HoldingSpec` -- same `covers` contract."""

    def __init__(self, holding_id="h", span_A=None):
        self.holding_id, self.span_A = holding_id, span_A

    def covers(self, centre, pad):
        if self.span_A is None:
            return True
        return self.span_A[0] <= centre - pad and centre + pad <= self.span_A[1]


# ── CRITICAL: telluric_satisfied on an IR band that is not-applied / unknown ──

@pytest.mark.parametrize("applied", ["not-applied", "unknown"])
def test_ir_band_is_never_satisfied_without_an_applied_correction(applied):
    """THE IR RULE, no exceptions. The positive control is the `applied` case below."""
    lo, hi = pr.TELLURIC_DOMAIN_MAX_A - 100.0, pr.TELLURIC_DOMAIN_MAX_A + 1000.0
    ok, tag, why = pr.telluric_satisfied(
        "crires_plus", "NIR", lo, hi, "correction_required", applied, "verified")
    assert ok is False
    assert tag == "ir-band-uncorrected"
    assert f"{pr.TELLURIC_DOMAIN_MAX_A:.0f}" in why


def test_ir_band_IS_satisfied_once_the_correction_is_applied_and_verified():
    """Positive control: the same band, same instrument, only `applied` differs.

    Without this the test above would also pass on a function that always returned False,
    which would gate every IR product forever and prove nothing about the rule.
    """
    lo, hi = pr.TELLURIC_DOMAIN_MAX_A - 100.0, pr.TELLURIC_DOMAIN_MAX_A + 1000.0
    ok, tag, _why = pr.telluric_satisfied(
        "crires_plus", "NIR", lo, hi, "correction_required", "applied", "verified")
    assert ok is True and tag == ""


def test_optical_line_selection_is_satisfied_but_only_inside_the_enumerated_domain():
    """HARPS optical runs on the RYA-460/786 clean-line basis; past the enumeration it does not."""
    inside = pr.telluric_satisfied(
        "harps", "VIS", 3780.0, 6910.0, "line_selection", "not-applied", "verified")
    assert inside[0] is True
    outside = pr.telluric_satisfied(
        "kpno_solar_atlas", "NIR", 9199.0, 13000.0, "line_selection", "not-applied",
        "verified")
    assert outside[0] is False and outside[1] == "ir-band-uncorrected"


def test_a_declared_ambiguity_is_not_a_declared_basis():
    """`mode_dependent` must not reach satisfied. See the run_bug_ledger entry."""
    ok, tag, _why = pr.telluric_satisfied(
        "uves", "VIS", 3780.0, 6910.0, "mode_dependent", "not-applied", "verified")
    assert ok is False and tag == "basis-mode_dependent"


def test_applied_on_unverified_evidence_is_not_satisfied():
    ok, tag, _ = pr.telluric_satisfied(
        "harps", "VIS", 3780.0, 6910.0, "correction_required", "applied", "audited")
    assert ok is False and tag == "applied-unverified"


def test_no_atmosphere_means_the_telluric_axis_does_not_gate():
    """RYA-806: HST is above the atmosphere, so `not-applied` is literally true and inert."""
    ok, _tag, why = pr.telluric_satisfied(
        "hst_stis", "VIS", 3780.0, 6910.0, "not_applicable", "not-applied", "verified")
    assert ok is True and "not_applicable" in why


# ── CRITICAL: an unwired holding is an honest NO-GO, never a crash ────────────

def test_an_unwired_holding_reports_reader_no_go_and_does_not_crash():
    centre, wired = pr._reader_reach(None, 3780.0, 6910.0, 0.62)
    assert wired is False and centre == pytest.approx(5345.0)


def test_reader_reach_is_asked_at_the_readers_own_overlap_not_the_band_centre():
    """The RYA-794 Y arm (10280-10680 A) SERVES the NIR band; it does not span it.

    Positive control for the bug this replaced: asking `covers()` at the band centre
    (11249.5 A) returns False for a holding that demonstrably measures Fe I there, and
    made `reader_wired` contradict `line_pool_reachable` on the same reader.
    """
    spec = _Spec("solar_crires_plus_y_rya794", span_A=(10280.0, 10680.0))
    assert spec.covers(0.5 * (9199.0 + 13000.0), 1.4) is False       # the old question
    _centre, wired = pr._reader_reach(spec, 9199.0, 13000.0, 1.4)    # the right one
    assert wired is True


def test_a_reader_that_does_not_reach_the_band_at_all_is_out_of_band():
    spec = _Spec("solar_iag_reiners2016", span_A=(4047.46, 5001.10))
    _centre, wired = pr._reader_reach(spec, 6910.0, 9199.0, 1.10)
    assert wired is False


def test_no_reader_means_zero_reachable_lines_not_an_empty_pool():
    """RYA-833: a wiring gap and an empty linelist must not collapse into one verdict."""
    ok, in_band, servable = pr.pool_reachable(None, 3780.0, 6910.0, 0.62, "Fe", "I")
    assert ok is False and servable == 0
    assert in_band > 0, "the graded Fe I pool is not empty in the VIS -- see canonical_gf"


# ── the gates are read from catalogs, never written down here ─────────────────

def test_band_edges_come_from_config_and_are_clipped_to_the_instrument():
    from config import synth_bands
    bands = dict((n, (lo, hi)) for n, lo, hi in pr.bands_for("harps"))
    assert set(bands) == {"VIS"}, "HARPS reaches 378-691 nm and nothing else"
    assert bands["VIS"][1] == pytest.approx(synth_bands.SYNTH_BANDS["VIS"].hi_A)
    stis = dict((n, (lo, hi)) for n, lo, hi in pr.bands_for("hst_stis"))
    assert stis["NIR"][1] == pytest.approx(10270.0), "clipped to the catalog reach"


def test_no_band_edge_or_absolute_path_is_hardcoded_in_the_conductor():
    """The numbers this file may contain are indices and formatting, not physics."""
    import re
    src = (ROOT / "scripts" / "preflight_readiness.py").read_text(encoding="utf-8")
    code = "\n".join(ln for ln in src.splitlines() if not ln.lstrip().startswith("#"))
    for edge in ("3780", "6910", "9199", "13000", "3000.0", "11560", "4500"):
        assert edge not in code, (
            f"{edge} is a band edge or a detector threshold and must be READ from "
            f"config/synth_bands.yaml, instrument_catalog.csv or the owning module")
    assert not re.search(r"""['"](/mnt/|/srv/|/home/|/Users/)""", code)


def test_the_graded_tier_set_is_imported_not_re_enumerated():
    from pipeline.gf_empirical import GRADED_TIERS
    assert pr.GRADED_TIERS is GRADED_TIERS


def test_the_telluric_domain_edge_is_derived_from_the_enumerated_bands():
    from pipeline.telluric_policy import TELLURIC_BANDS
    assert pr.TELLURIC_DOMAIN_MAX_A == max(hi for _lo, hi, _n in TELLURIC_BANDS)


# ── product selection: the RYA-1064 gate ─────────────────────────────────────

def test_an_n_file_inventory_is_not_a_selected_product():
    """The RYA-1064 STOP, reproduced: 205 loose ADPs are not one approved run product."""
    ok, tag, why = pr.product_selection(
        "alpha_cen_a_harps", "alpha_cen_a", "harps",
        "data/audit/alpha_cen_optical_ir/acen_optical_ir_manifest_rya479.csv")
    assert ok is False and tag.startswith("n-file-inventory")
    assert "RYA-1064" in why


def test_a_shared_manifest_is_attributed_per_holding_not_counted_whole():
    """One file serves six alpha Cen holdings; each must see only its own rows."""
    _ok, harps, _ = pr.product_selection(
        "alpha_cen_a_harps", "alpha_cen_a", "harps",
        "data/audit/alpha_cen_optical_ir/acen_optical_ir_manifest_rya479.csv")
    _ok, crires, _ = pr.product_selection(
        "alpha_cen_a_crires_plus", "alpha_cen_a", "crires_plus",
        "data/audit/alpha_cen_optical_ir/acen_optical_ir_manifest_rya479.csv")
    assert harps != crires, "the two arms of one manifest cannot have the same count"
    assert "(88)" in harps and "(6)" in crires


def test_the_facility_name_for_an_instrument_resolves_to_its_registered_id():
    """The alpha Cen manifest writes `CRIRES` for products registered `crires_plus`."""
    assert pr._instrument_alias("CRIRES") == "crires_plus"
    assert pr._instrument_alias("HARPS") == "harps", "an exact match never widens"
    assert pr._instrument_alias("zzz") is None


def test_a_one_row_catalogue_entry_IS_a_selected_product():
    """Positive control: without it, `product` could be a gate that never passes."""
    ok, tag, _why = pr.product_selection(
        "solar_harps", "solar", "harps",
        "data/catalog/solar_reference_holdings_rya708.csv")
    assert ok is True and tag == "selected"


def test_a_missing_manifest_is_not_selected():
    ok, tag, _why = pr.product_selection("x", "y", "harps", "data/does/not/exist.csv")
    assert ok is False and tag == "manifest-missing"


# ── the artifact and the refusal ─────────────────────────────────────────────

def test_an_unregistered_system_is_refused_rather_than_guessed():
    with pytest.raises(SystemExit) as exc:
        pr.holdings_for_system("no_such_star", None)
    assert "refuse, do not guess" in str(exc.value)


def test_the_readiness_schema_is_exactly_the_ticket_spec():
    from dataclasses import fields
    assert [f.name for f in fields(pr.ReadinessRow)] == [
        "system_id", "holding_id", "instrument_id", "band", "coverage_A",
        "evidence_state", "product_selected", "normalization_state", "telluric_basis",
        "telluric_applied", "telluric_satisfied", "reader_wired", "line_pool_reachable",
        "measurement_ready", "blocking_gate", "source_issue_ids"]


@pytest.mark.parametrize("system", ["alpha_cen_a", "solar"])
def test_no_committed_readiness_row_is_GO_with_an_unknown_gate(system):
    """The artifact itself is the control: a GO beside an `unknown` is the CRITICAL case."""
    path = ROOT / "data" / "audit" / "readiness" / f"{system}_readiness.csv"
    if not path.exists():
        pytest.skip(f"{path.name} not generated in this checkout")
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    assert rows, f"{path.name} is empty"
    for r in rows:
        if r["measurement_ready"] != "GO":
            assert r["blocking_gate"], "a NO-GO must name its blocking gate"
            continue
        assert not r["blocking_gate"], "a GO carries no blocking gate"
        assert r["normalization_state"] != "unknown"
        assert r["evidence_state"] == "verified"
        assert r["telluric_satisfied"] == "True" and r["reader_wired"] == "True"
        assert r["product_selected"] == "True" and r["line_pool_reachable"] == "True"
        if r["telluric_basis"] != "not_applicable":
            assert r["telluric_applied"] != "unknown"


def test_alpha_cen_a_harps_still_names_the_three_RYA_1064_gaps():
    """The acceptance case. If any of these three ever passes, it was CONDITIONED --
    or the gate broke. Either way this test is the place that says so."""
    path = ROOT / "data" / "audit" / "readiness" / "alpha_cen_a_readiness.csv"
    if not path.exists():
        pytest.skip("alpha_cen_a_readiness.csv not generated in this checkout")
    rows = [r for r in csv.DictReader(path.open(encoding="utf-8"))
            if r["holding_id"] == "alpha_cen_a_harps"]
    assert rows, "alpha_cen_a_harps is registered and must appear"
    for r in rows:
        assert r["measurement_ready"] == "NO-GO"
        gates = {g.split(":", 1)[0] for g in r["blocking_gate"].split(";")}
        assert {"product", "normalization", "reader"} <= gates

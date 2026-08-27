"""RYA-1079: correction is the standard; only a saturated core is genuinely lost.

The policy under test: a line's telluric state is DECIDED FROM MEASURED DEPTH over its own
window, not asserted from `instrument_catalog.telluric_basis`; observability runs BEFORE
grading; and a line lost to a saturated core is recorded UNMEASURED-WITH-REASON, never
"ungraded".

BOTH POSITIVE CONTROLS THE TICKET DEMANDS ARE HERE, and they are the same test's two
halves: one clearly-recoverable line that the RYA-460/786 band rule would have skipped and
this classifier admits, and one clearly-saturated line it rejects. A classifier shown only
to reject proves nothing -- it could be a function that always says no.

The measured half runs off the COMMITTED census artifact so it runs on any machine; the
live measurement needs the Kitt Peak atlas staged and is skipped where it is not (a
skipped test is not a passing test -- the artifact tests are what carry the claim).
"""
import csv
from pathlib import Path

import numpy as np
import pytest

from pipeline import telluric_observability as OBS

ROOT = Path(__file__).resolve().parents[1]
CENSUS = (ROOT / "data" / "audit" / "rya1079_observability"
          / "solar_kpno_red-optical_observability.csv")
PC = (ROOT / "data" / "audit" / "rya1079_observability"
      / "solar_kpno_red-optical_problem_children.csv")



# ── source introspection: CODE only, never the module's own prose ────────────
#
# The module docstring necessarily NAMES `telluric_basis` and `line_selection` -- it is
# explaining what this ticket retires. A raw grep cannot tell an explanation from a
# decision and fails on the module's own documentation, so these read the parsed tree.

def _tree():
    import ast
    return ast.parse((ROOT / "pipeline" / "telluric_observability.py")
                     .read_text(encoding="utf-8"))


def _docstring_ids():
    import ast
    out = set()
    for node in ast.walk(_tree()):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)):
            body = getattr(node, "body", None)
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                out.add(id(body[0].value))
    return out


def _code_strings():
    import ast
    skip = _docstring_ids()
    return {n.value for n in ast.walk(_tree())
            if isinstance(n, ast.Constant) and isinstance(n.value, str)
            and id(n) not in skip}


def _code_literals():
    import ast
    return {repr(n.value) for n in ast.walk(_tree())
            if isinstance(n, ast.Constant) and isinstance(n.value, (int, float))} | {
            str(n.value) for n in ast.walk(_tree())
            if isinstance(n, ast.Constant) and isinstance(n.value, (int, float))}


def _code_names():
    import ast
    return {n.id for n in ast.walk(_tree()) if isinstance(n, ast.Name)} | {
            n.attr for n in ast.walk(_tree()) if isinstance(n, ast.Attribute)}


def _census():
    if not CENSUS.exists():
        pytest.skip("census artifact not generated in this checkout")
    return list(csv.DictReader(CENSUS.open(encoding="utf-8")))


# ── the threshold is DERIVED, not borrowed (RYA-161) ─────────────────────────

def test_both_edges_fall_out_of_one_photon_noise_argument():
    """No magic number: both edges come from S0 and PIPELINE['snr_min_science']."""
    th = OBS.thresholds(1000.0)
    assert th.snr_min_science == OBS.snr_floor()
    assert th.clean_max_depth == pytest.approx(1.0 / 1000.0)
    assert th.saturated_min_depth == pytest.approx(1.0 - (OBS.snr_floor() / 1000.0) ** 2)


def test_a_noisier_spectrum_can_recover_LESS():
    """The physics the derivation encodes, asserted as a direction rather than a value.

    Measured on our own arms: Kitt Peak S0~1382 in the red-optical tolerates depth to
    0.979; HARPS S0~458 only to 0.809. A single global cut would be wrong for both.
    """
    deep, shallow = OBS.thresholds(2000.0), OBS.thresholds(400.0)
    assert deep.saturated_min_depth > shallow.saturated_min_depth
    assert deep.clean_max_depth < shallow.clean_max_depth


def test_no_threshold_literal_is_written_into_the_module():
    """A borrowed constant would be a CRITICAL failure condition on this ticket."""
    import ast
    import inspect
    # Assert the PROPERTY, not the absence of digits: the verdict may only be decided by
    # comparing the measured depth against the two DERIVED edges. An earlier draft banned
    # numeric literals outright and flagged `(i + 0.5)`, a bin-centring offset in a probe
    # loop -- geometry, not a threshold. Banning digits is not the same as banning magic
    # numbers.
    fn = ast.parse(inspect.getsource(OBS.observe).lstrip()).body[0]
    compared = set()
    for node in ast.walk(fn):
        if isinstance(node, ast.Compare):
            for side in [node.left, *node.comparators]:
                if isinstance(side, ast.Constant) and isinstance(side.value, (int, float)):
                    compared.add(side.value)
                if isinstance(side, ast.Attribute):
                    compared.add(side.attr)
    assert "clean_max_depth" in compared and "saturated_min_depth" in compared, compared
    numeric = {c for c in compared if not isinstance(c, str)}
    assert not numeric, (
        f"the verdict compares the measured depth against bare number(s) {numeric} "
        f"instead of the derived edges -- every threshold on this ticket must be derived "
        f"(RYA-161)")


def test_the_snr_floor_is_read_from_the_repo_constant():
    from config.constants import PIPELINE
    assert OBS.snr_floor() == float(PIPELINE["snr_min_science"])


# ── the three states, and the dispositions they imply ────────────────────────

def test_continuum_snr_is_measured_from_pixel_scale_scatter():
    """Differencing must read the NOISE, not the solar structure sitting under it."""
    rng = np.random.default_rng(1079)
    # Structure that is genuinely smooth on the pixel scale. An earlier draft used a
    # cosine over 6 radians in 4000 points, whose per-pixel step (~1.5e-3) is the SAME
    # SIZE as the noise at S/N 1000 -- so the test was measuring its own fixture.
    smooth = np.cos(np.linspace(0, 0.5, 20000))
    for true_snr in (200.0, 1000.0, 2500.0):
        flux = smooth + rng.normal(0, 1.0 / true_snr, smooth.size)
        assert OBS.continuum_snr(flux) == pytest.approx(true_snr, rel=0.10), (
            "differencing must measure the NOISE, not the solar structure")


def test_the_estimator_errs_CONSERVATIVELY_when_structure_is_pixel_scale():
    """The bias is stated in the module and asserted here: unresolved structure inflates
    the measured noise, which LOWERS S0, which lowers `saturated_min_depth` -- so more
    lines are called saturated. That costs coverage; it never admits an unmeasurable line."""
    rng = np.random.default_rng(11)
    sharp = np.cos(np.linspace(0, 60, 20000))          # varies on the pixel scale
    flux = sharp + rng.normal(0, 1.0 / 2000.0, sharp.size)
    assert OBS.continuum_snr(flux) < 2000.0


def test_a_flat_noiseless_window_refuses_rather_than_claiming_infinite_snr():
    with pytest.raises(OBS.TransmissionUnavailable):
        OBS.continuum_snr(np.ones(500))


def test_an_unmeasured_depth_is_never_defaulted_either_way():
    """🔴 The RYA-833/1072 rule, on this axis. Not clean, not saturated -- unmeasured."""
    assert issubclass(OBS.TransmissionUnavailable, RuntimeError)
    src = (ROOT / "pipeline" / "telluric_observability.py").read_text(encoding="utf-8")
    assert "NEVER DEFAULTED" in src


def test_a_low_snr_spectrum_is_not_blamed_on_the_atmosphere():
    """🔴 `saturated_min_depth` goes NEGATIVE when S0 < SNR_min, which would call every
    line -- including ones under a clear sky -- telluric_saturated_core. That is a false
    accusation against the atmosphere for a property of the DATA, so it is a separate
    exception type."""
    assert OBS.thresholds(50.0).saturated_min_depth < 0
    assert OBS.ContinuumBelowScienceFloor is not OBS.TransmissionUnavailable


# ── OBSERVABILITY BEFORE TIER (spec 3) ───────────────────────────────────────

def test_a_saturated_line_is_recorded_unmeasured_NOT_ungraded():
    """The label-vs-reality rule. `Partition.unmeasured` exists and `ungraded` does not."""
    p = OBS.Partition()
    assert hasattr(p, "unmeasured") and not hasattr(p, "ungraded")
    src = (ROOT / "pipeline" / "telluric_observability.py").read_text(encoding="utf-8")
    assert "unmeasured, with a reason" in src or "unmeasured" in src


def test_no_saturated_line_appears_in_the_measurable_pool():
    """The ordering, asserted on the real census: the two sets must not intersect."""
    rows = _census()
    sat = {r["wavelength_air_A"] for r in rows if r["verdict"] == OBS.SATURATED}
    keep = {r["wavelength_air_A"] for r in rows
            if r["verdict"] in (OBS.CLEAN, OBS.RECOVERABLE)}
    assert sat, "the census must contain at least one saturated line to be a control"
    assert not (sat & keep), f"saturated lines leaked into the measurable pool: {sat & keep}"


def test_every_saturated_line_carries_its_reason_and_its_measured_depth():
    for r in _census():
        if r["verdict"] != OBS.SATURATED:
            continue
        assert r["reason"] == OBS.SATURATED_CORE_REASON
        assert float(r["depth"]) >= float(r["saturated_min_depth"])
        assert r["disposition"] == OBS.EXCLUDE


# ── THE TWO POSITIVE CONTROLS, in one place ──────────────────────────────────

def test_the_classifier_both_ADMITS_and_REJECTS_inside_one_telluric_complex():
    """🔴 THE TICKET'S REQUIRED CONTROL PAIR, and deliberately drawn from the SAME
    H2O 8100-8400 complex so band membership cannot be what separates them.

        Fe I 8198.921  depth 0.53  -> RECOVERABLE, rescued (RYA-460/786 would skip it)
        Fe I 8232.316  depth 0.98  -> SATURATED, genuinely lost

    Both sit inside a registered telluric band. Under the old rule both were excluded.
    Under this one, the measured depth adjudicates and one of them comes back.
    """
    rows = {r["wavelength_air_A"]: r for r in _census()}
    rescued = rows.get("8198.9210")
    lost = rows.get("8232.3160")
    assert rescued is not None and lost is not None, sorted(rows)[:5]

    assert rescued["verdict"] == OBS.RECOVERABLE
    assert rescued["disposition"] == OBS.NEEDS_CORRECTION
    assert 0.4 < float(rescued["depth"]) < 0.7

    assert lost["verdict"] == OBS.SATURATED
    assert lost["disposition"] == OBS.EXCLUDE
    assert float(lost["depth"]) > float(rescued["depth"])

    from pipeline.telluric_policy import TELLURIC_BANDS
    for r in (rescued, lost):
        w = float(r["wavelength_air_A"])
        assert any(lo <= w <= hi for lo, hi, _n in TELLURIC_BANDS), (
            "both controls must be INSIDE an enumerated band, or they prove nothing "
            "about superseding the band rule")


def test_the_rescued_class_is_not_empty_on_a_real_band():
    """If this ever hits zero, the policy has stopped rescuing anything."""
    rows = _census()
    rec = [r for r in rows if r["verdict"] == OBS.RECOVERABLE]
    assert len(rec) >= 5, f"only {len(rec)} recoverable lines -- the rescue is vacuous"
    assert all(r["disposition"] == OBS.NEEDS_CORRECTION for r in rec)


def test_clean_lines_are_the_majority_and_are_measured_as_is():
    rows = _census()
    clean = [r for r in rows if r["verdict"] == OBS.CLEAN]
    assert len(clean) > len(rows) / 2
    assert all(r["disposition"] == OBS.MEASURE for r in clean)


# ── spec 4: `line_selection` is retired as a DECISION INPUT ──────────────────

def test_the_verdict_never_consults_telluric_basis():
    """🔴 CRITICAL failure condition: the state must be MEASURED, not asserted from the
    instrument column. The module may not read the basis axis at all."""
    for forbidden in ("telluric_basis", "line_selection"):
        assert forbidden not in _code_strings(), (
            f"{forbidden!r} appears in the verdict CODE; RYA-1079 retires the basis "
            f"column as a decision input -- the measured depth decides")


def test_the_registered_telluric_bands_are_not_used_to_exclude_anything():
    """The band list stays useful as an enumeration; it is no longer the DECISION."""
    assert "TELLURIC_BANDS" not in _code_names()


# ── the problem_children sink carries the measurement (spec 2) ───────────────

def test_the_saturated_row_is_a_valid_problem_children_row_carrying_its_depth():
    from pipeline.problem_children import SCHEMA_COLUMNS, TREATMENTS
    o = OBS.Observability(
        holding_id="solar_kpno", wavelength_A=7620.513, verdict=OBS.SATURATED,
        disposition=OBS.EXCLUDE, depth=1.0, transmission_min=0.0,
        thresholds=OBS.thresholds(1382.4), corrected=False,
        reason=OBS.SATURATED_CORE_REASON, evidence="measured")
    row = OBS.problem_child_row(o, "Fe I")
    assert list(row) == SCHEMA_COLUMNS
    assert row["required_treatment"] == "exclude" and "exclude" in TREATMENTS
    assert OBS.SATURATED_CORE_REASON in row["notes"]
    assert "T=0.0000" in row["notes"], "the measured depth must travel with the exclusion"
    assert "never graded" in row["notes"]


def test_the_emitted_problem_children_file_matches_the_saturated_census():
    if not PC.exists():
        pytest.skip("problem_children artifact not generated in this checkout")
    emitted = {r["lambda_or_scope"] for r in csv.DictReader(PC.open(encoding="utf-8"))}
    sat = {f"{float(r['wavelength_air_A']):.3f}" for r in _census()
           if r["verdict"] == OBS.SATURATED}
    assert emitted == sat


# ── spec 5: the conductor routes RECOVERABLE-uncorrected as WORK ─────────────

def test_the_conductor_distinguishes_needs_correction_from_a_terminal_no_go():
    path = ROOT / "data" / "audit" / "readiness" / "solar_readiness.csv"
    if not path.exists():
        pytest.skip("solar readiness artifact not generated in this checkout")
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    verdicts = {r["measurement_ready"] for r in rows}
    assert "NEEDS-CORRECTION" in verdicts, (
        "no actionable telluric cell in the artifact -- spec 5 has no positive control")
    for r in rows:
        if r["measurement_ready"] != "NEEDS-CORRECTION":
            continue
        gates = {g.split(":", 1)[0] for g in r["blocking_gate"].split(";")}
        assert gates == {"telluric"}, (
            "NEEDS-CORRECTION must mean telluric is the ONLY thing owed; anything else "
            "failing makes it a NO-GO with work still to do elsewhere")
    assert "GO" in verdicts and "NO-GO" in verdicts, "all three verdicts must be reachable"

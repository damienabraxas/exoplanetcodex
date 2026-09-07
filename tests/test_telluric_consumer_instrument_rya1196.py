"""RYA-1196 — no production measurement path may ask the telluric policy without saying
WHICH data it is asking about.

`measure_band_ew.py` called `telluric_reason(wave)` with no instrument. That is not the
neutral setting: with nothing to resolve, `telluric_policy.exclusion` sees basis
`unspecified` and EVERY registered band fires unconditionally. In a measurement path that
silently over-excludes real lines, and it is the same axis-collapse trap RYA-1193 found in
the GBS line set and RYA-1194 fixed in `exclusion()` itself.

The guard below is an AST scan, not a grep: it reads the CODE, so a mention in a docstring
neither trips it nor satisfies it (RYA-1079's lesson).
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

#: Functions that answer "is this line telluric-excluded?" and therefore need to know
#: whose data is being asked about.
_RESOLVERS = {"telluric_reason", "exclusion", "_telluric_exclusion",
              "telluric_exclusion", "in_telluric_band"}

#: The paths that decide what gets MEASURED or what lands in a line set. A wrong answer
#: here changes a number; elsewhere it changes a report.
PRODUCTION_PATHS = (
    "pipeline/measure/synthesis.py",
    "scripts/measure_band_ew.py",
    "scripts/measure_band_profilefit.py",
    "scripts/derive_band_products.py",
    "scripts/rya1110_build_gbs_fe_lineset.py",
)

#: Deliberate no-argument callers, exempted BY NAME and never by pattern -- a pattern
#: exemption would grow to cover the next real defect (RYA-686's lesson, where the audit
#: script matched its own grep). Each entry must still CONTAIN such a call, asserted
#: below, so an exemption cannot outlive the thing it excuses.
DELIBERATE_POLICY_ONLY = {
    "scripts/rya1192_telluric_verification.py":
        "an AUDIT of the policy itself: it asks what the ENUMERATION says about a "
        "wavelength, with no holding, because that is precisely the question. It measures "
        "the policy rather than measuring a line with it.",
}


def _calls(path: Path):
    """(lineno, func_name, n_args, has_instrument_kw) for every resolver call in CODE."""
    for n in ast.walk(ast.parse(path.read_text())):
        if not isinstance(n, ast.Call):
            continue
        fn = getattr(n.func, "id", None) or getattr(n.func, "attr", None)
        if fn not in _RESOLVERS:
            continue
        kw = {k.arg for k in n.keywords}
        yield n.lineno, fn, len(n.args) + len(n.keywords), bool(
            kw & {"instrument", "holding_id", "holding"})


@pytest.mark.parametrize("rel", PRODUCTION_PATHS)
def test_no_production_call_omits_the_instrument(rel):
    """The spec's guard. A one-argument call is the most-aggressive setting."""
    bad = [(ln, fn) for ln, fn, n, kw in _calls(ROOT / rel) if n < 2 and not kw]
    assert not bad, (
        f"{rel}: telluric resolver called with only a wavelength at line(s) "
        f"{[ln for ln, _ in bad]}. That is telluric_basis=unspecified, which fires every "
        f"band regardless of the data. Pass the instrument, and the holding where the "
        f"caller has one (RYA-1194).")


def test_the_rya1196_target_line_passes_the_instrument():
    """Named explicitly so the fix cannot be reverted without this failing by name."""
    src = (ROOT / "scripts/measure_band_ew.py").read_text()
    assert "telluric_reason(r.wave_air_A, a.instrument)" in src
    assert "telluric_reason(r.wave_air_A)" not in src


def test_the_wrapper_can_carry_a_holding():
    """RYA-1194's third axis has to be reachable from this consumer, or the next caller
    that DOES know its holding still cannot say so."""
    import inspect, sys
    sys.path.insert(0, str(ROOT / "scripts"))
    import measure_band_ew
    p = inspect.signature(measure_band_ew.telluric_reason).parameters
    assert "instrument" in p and "holding_id" in p


def test_every_exemption_is_real_and_still_needed():
    """An exemption list that outlives its reason is how the next defect gets waved
    through. Each exempt file must still contain the call it is exempt for."""
    for rel, why in DELIBERATE_POLICY_ONLY.items():
        p = ROOT / rel
        assert p.exists(), f"{rel} is exempted but does not exist — stale exemption"
        bare = [ln for ln, _fn, n, kw in _calls(p) if n < 2 and not kw]
        assert bare, (
            f"{rel} no longer makes a no-argument policy call, so its exemption is dead "
            f"and must be deleted. Reason on file: {why}")
        assert rel not in PRODUCTION_PATHS, f"{rel} cannot be both production and exempt"


def test_the_guard_would_fire_on_the_defect_it_was_written_for():
    """Mutation control: feed the scanner the exact shape RYA-1196 fixed and it must
    object. Without this the parametrised test above passes just as well when `_calls`
    silently returns nothing."""
    import tempfile
    src = "def f(w):\n    return telluric_reason(w)\n"
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(src)
        tmp = Path(fh.name)
    try:
        bad = [(ln, fn) for ln, fn, n, kw in _calls(tmp) if n < 2 and not kw]
        assert bad, "the scanner does not detect a bare telluric_reason(wave) call"
    finally:
        tmp.unlink()

    # ...and it must NOT object to the fixed shape.
    src2 = "def f(w, a):\n    return telluric_reason(w, a.instrument)\n"
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(src2)
        tmp2 = Path(fh.name)
    try:
        assert not [1 for _l, _f, n, kw in _calls(tmp2) if n < 2 and not kw]
    finally:
        tmp2.unlink()


def test_the_fix_moves_no_value_on_this_driver():
    """🔴 THE VALUE CLAIM, MEASURED. This driver pins the instrument to kpno_solar_atlas,
    whose telluric_basis is `line_selection` — the same verdict `unspecified` produces. So
    the corrected call returns the same DECISION at every wavelength and no number moves.

    The defect was the reasoning, not the number: point the driver at a corrected holding
    and the bare call would have thrown those lines away."""
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from measure_band_ew import telluric_reason
    from pipeline.telluric_policy import TELLURIC_BANDS, basis

    assert basis("kpno_solar_atlas") == "line_selection", "premise gone"
    probes = [3000.0, 5000.0, 6271.0, 6875.0, 7200.0, 7600.0, 9400.0, 30000.0]
    for lo, hi, _n in TELLURIC_BANDS:
        probes += [lo - 0.01, lo, (lo + hi) / 2, hi, hi + 0.01]
    for w in probes:
        assert bool(telluric_reason(w)) == bool(
            telluric_reason(w, "kpno_solar_atlas")), w

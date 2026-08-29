"""RYA-1089 — the solar Type B artifact must agree with the CODE that defines delta_p.

🔴 THE DEFECT THIS EXISTS TO CATCH. RYA-1093 moved the solar microturbulence allowance to
0.2912 km/s in `pipeline/uncertainty_stack.py`. `solar_uncertainty_rya158.json` kept
RYA-158's retired 0.05. Both sat on main together, and all three existing readers of the
artifact passed, because not one of them compared the file to the module.

The stale value is also the FLATTERING one — 0.05 puts sigma_B_vmic at 0.0120 dex, inside
the 0.05 dex solar gate, while the real 0.2912 puts it at 0.0699 dex, which fails. An
artifact drifting toward passing a gate is the RYA-161 hazard arriving by accident, so the
agreement is asserted rather than assumed.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import uncertainty_stack  # noqa: E402

ARTIFACT = ROOT / "data" / "audit" / "uncertainty" / "solar_uncertainty_rya158.json"
GATE_DEX = 0.05


@pytest.fixture(scope="module")
def doc():
    return json.loads(ARTIFACT.read_text())


@pytest.fixture(scope="module")
def fe(doc):
    return next(r for r in doc["per_element"]
                if r["element"] == "Fe" and r["ion"] == "I")


def test_artifact_delta_p_equals_the_code(doc):
    """THE GUARD. The artifact may not carry its own copy of a value the code defines."""
    _, deltas = uncertainty_stack.params_and_deltas("solar")
    art = doc["delta_p_solar"]
    assert art["vmic"] == pytest.approx(deltas["vturb_kms"]), (
        f"artifact delta_p vmic={art['vmic']} but uncertainty_stack says "
        f"{deltas['vturb_kms']} — the file has drifted from the module again")
    assert art["Teff"] == pytest.approx(deltas["teff_K"])
    assert art["logg"] == pytest.approx(deltas["logg"])
    assert art["FeH"] == pytest.approx(deltas["feh"])


def test_the_retired_005_is_gone(doc, fe):
    """RYA-158's uncited 0.05 and the sigma_B it produced must not be reachable."""
    assert doc["delta_p_solar"]["vmic"] != 0.05
    assert fe["sigma_B_vmic"] != 0.012


def test_the_borrowed_numbers_were_not_adopted(doc, fe):
    """Neither Jofre's 0.18 nor RYA-311's 0.0588 formal error is the stamped allowance.

    0.0588 is the one that matters: it CLEARS the gate, and RYA-1093 2E names choosing it
    for that reason as the RYA-161 failure."""
    assert doc["delta_p_solar"]["vmic"] not in (0.18, 0.0588)


def test_the_xi_term_alone_fails_the_gate_and_is_reported_anyway(fe):
    """Gate-as-FLAG. The honest bar exceeds the target and is published, not tuned down."""
    assert fe["sigma_B_vmic"] > GATE_DEX
    assert fe["sigma_B_vmic"] == pytest.approx(0.0699, abs=5e-5)


def test_sigma_reported_is_recomputed_not_stale(fe):
    """sigma_reported must move with sigma_B, or the budget understates the bar."""
    se = fe["raw_sigma"] / math.sqrt(fe["n_lines"])
    sp = math.sqrt(fe["sigma_B_Teff"] ** 2 + fe["sigma_B_vmic"] ** 2)
    assert fe["sigma_reported"] == pytest.approx(math.sqrt(se ** 2 + sp ** 2), abs=5e-5)
    assert fe["sigma_reported"] > fe["sigma_solar"] - 1e-9


def test_the_derivatives_are_in_three_states_and_the_file_says_so(doc):
    """Only Fe I is a RESULT. Eight elements are NaN (unmeasured, RYA-907). S I is a third
    state — an EXACT 0.0 on both derivatives at n=2 — which is implausible for a real line
    set and is flagged rather than trusted.

    Asserted so that measuring one later is a deliberate change, and so the caveat cannot
    quietly stop being true while the zeros stay."""
    nan_rows, zero_rows, measured = [], [], []
    for r in doc["per_element"]:
        v = r.get("dA_dvmic_per_kms", float("nan"))
        t = r.get("dA_dTeff_per100K", float("nan"))
        if math.isnan(v):
            nan_rows.append(r["element"])
        elif v == 0.0 and t == 0.0:
            zero_rows.append(r["element"])
        else:
            measured.append(r["element"])
    assert measured == ["Fe"]
    assert zero_rows == ["S"]
    assert nan_rows == ["Ba", "C", "Co", "Cu", "Mn", "N", "O", "V"]
    assert any("RYA-907" in c for c in doc["caveats"])
    assert any("S I" in c for c in doc["caveats"])


def test_no_element_but_Fe_gained_a_sigma_B_from_the_restamp(doc):
    """The re-stamp changes delta_p for the whole star, but only a row with a real
    derivative may move. If a NaN or an implausible 0.0 started producing a non-zero
    term, the generator would be manufacturing a budget out of an absence."""
    for r in doc["per_element"]:
        if r["element"] != "Fe":
            assert r["sigma_B_vmic"] == 0.0, r["element"]
            assert r["sigma_B_Teff"] == 0.0, r["element"]


def test_the_generator_reproduces_the_artifact_from_its_own_delta_p():
    """THE CONTROL, run as a test. A generator that cannot reproduce the file it is about
    to rewrite cannot tell 'the input changed' from 'my arithmetic is wrong'."""
    r = subprocess.run(
        [sys.executable, "scripts/rya1089_stamp_honest_delta_xi.py", "--control"],
        cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "CONTROL REPRODUCED" in r.stdout


def test_the_control_is_not_vacuous(tmp_path):
    """The control must FAIL on a perturbed artifact. A control that passes on anything
    proves nothing about the file it blessed."""
    bad = json.loads(ARTIFACT.read_text())
    for r in bad["per_element"]:
        if r["element"] == "Fe":
            r["sigma_B_vmic"] = round(r["sigma_B_vmic"] + 0.001, 4)
    p = tmp_path / "perturbed.json"
    p.write_text(json.dumps(bad, indent=2))
    r = subprocess.run(
        [sys.executable, "scripts/rya1089_stamp_honest_delta_xi.py",
         "--control", "--artifact", str(p)],
        cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 1
    assert "CONTROL FAILED" in r.stdout

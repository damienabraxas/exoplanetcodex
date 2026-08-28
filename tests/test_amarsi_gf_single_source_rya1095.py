"""
RYA-1095 — the Amarsi 3D-NLTE leg must run on the SAME oscillator strengths as its base.

🔴 THE DEFECT, MEASURED. `scripts/rya817_run_3dnlte_bands.py` reads `a_1dlte` straight out
of a band product built with canonical (laboratory) gf applied — and then resolved the
Elo/Eup/log gf beside it from VALD. So the published product paired a LAB-gf abundance
with a VALD-gf correction. `log gf` is an INPUT FEATURE of the Amarsi MLP, not a label, so
the mismatch propagated into the correction itself.

On the 2026-08-24 redo (Fe I VIS, 50 lines, every one `gf_tier=LAB` in canonical_gf) the
used-minus-canonical log gf had mean −0.166 dex, max |Δ| 2.199, and ZERO exact matches.

⚠️ The first thing to rule out was a misidentified line, and it is ruled out: 49 of the 50
had exactly ONE VALD candidate in the window, and EP agreed with canonical to <0.001 eV on
every line including the three worst. The transitions were right; the gf came from the
wrong table.
"""
from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "rya817_run_3dnlte_bands.py"
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"


def test_the_runner_canonicalises_its_gf():
    """The fix, pinned by structure rather than by a string: `attach_atomic` must hand its
    result through the canonicaliser, and the canonicaliser must call `resolve_df_gf`."""
    tree = ast.parse(RUNNER.read_text())
    fns = {n.name: n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    assert "_canonicalise_gf" in fns, "the gf canonicaliser is gone"
    src = ast.get_source_segment(RUNNER.read_text(), fns["_canonicalise_gf"])
    assert "resolve_df_gf" in src, "must resolve through pipeline.gf_resolver"
    assert "keep_unresolved=False" in src, (
        "an authoritative use gets NO silent fallback -- a half-canonical pool is harder "
        "to reason about than either honest alternative")
    attach = ast.get_source_segment(RUNNER.read_text(), fns["attach_atomic"])
    assert "_canonicalise_gf(" in attach, (
        "attach_atomic must route its output through the canonicaliser, or the VALD gf "
        "reaches the network unchanged")


def test_the_upper_level_still_comes_from_VALD():
    """Only the gf moves. canonical_gf carries no upper level and the network needs the
    transition energy, so Eup must keep its VALD provenance -- replacing it would be a
    different (and unsourced) change riding along with this one."""
    tree = ast.parse(RUNNER.read_text())
    fns = {n.name: n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    src = ast.get_source_segment(RUNNER.read_text(), fns["_canonicalise_gf"])
    assert "eup_eV" not in src.split('"""')[-1], (
        "the canonicaliser must not touch Eup -- canonical_gf has no upper level")
    assert "excitation_potential_eV" not in src or "ep_col" in src


@pytest.mark.skipif(not CANON.exists(), reason="canonical_gf absent")
def test_canonical_resolution_covers_the_pool_it_will_be_asked_for():
    """`keep_unresolved=False` RAISES on a line canonical_gf does not carry, so the fix
    only works if the graded Fe I VIS pool is actually covered. Asserted rather than
    assumed: a refusal that fires on every run is not a guard, it is an outage."""
    from pipeline.gf_resolver import GfResolutionError, resolve
    from pipeline.species import species_key

    cg = pd.read_csv(CANON, low_memory=False)
    lab = cg[(cg.species == "Fe I")
             & cg.gf_tier.astype(str).str.contains("LAB", na=False)
             & cg.wavelength_air_A.between(4200, 6910)]
    assert len(lab) > 50, f"only {len(lab)} LAB Fe I lines in VIS -- pool too thin to test"
    key = species_key("Fe", "I")
    unresolved = []
    for r in lab.itertuples():
        try:
            resolve(key, float(r.wavelength_air_A), float(r.excitation_potential_eV))
        except GfResolutionError:
            unresolved.append(float(r.wavelength_air_A))
    assert not unresolved, (
        f"{len(unresolved)} LAB Fe I VIS line(s) do not resolve against canonical_gf "
        f"itself: {unresolved[:6]}")


def test_the_measured_mismatch_is_recorded_where_someone_will_find_it():
    """The numbers that justify the change live beside the code that makes it. A fix whose
    evidence is only in a ticket comment is a fix nobody can re-derive."""
    src = RUNNER.read_text()
    for token in ("-0.166", "2.199", "gf_tier=LAB", "RYA-353"):
        assert token in src, f"the canonicaliser's docstring lost {token!r}"

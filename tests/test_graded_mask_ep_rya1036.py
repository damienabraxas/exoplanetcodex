"""
tests/test_graded_mask_ep_rya1036.py — RYA-1036
===============================================
`_graded_mask` decides which measured lines carry a PRIMARY LABORATORY gf. It used to key on
wavelength alone within 5 mA, so a LAB line whose canonical row sits further away was
published UNGRADED — tier disagreeing with provenance. That is the product-path twin of
RYA-871 and the mirror of RYA-1034 (there an ungraded line wore a GRADED tier by λ
coincidence; here a lab line wore an UNGRADED one).

The sequence in the ticket is non-negotiable and these tests pin it: **EP first, widening
second.** Widening to 30 mA on wavelength alone would swallow the 340 Fe I clusters that sit
within 5 mA of each other AND disagree on gf — e.g. 6065.4820 (EP 2.609, NIST A, −1.530)
against 6065.4850 (EP 4.956, K07, −3.471), 1.94 dex apart. That is the coin flip RYA-1033
killed and it must not come back through a widened window.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import line_match  # noqa: E402
from pipeline.line_match import LineMatchError  # noqa: E402


@pytest.fixture(scope="module")
def dbp():
    spec = importlib.util.spec_from_file_location(
        "dbp_rya1036", ROOT / "scripts" / "derive_band_products.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_the_mask_requires_EP_and_will_not_take_wavelengths_alone(dbp):
    """The whole point. A caller without EP must thread it through, not ask for a looser
    key — the old signature took wavelengths and could not do anything else."""
    import inspect
    sig = inspect.signature(dbp._graded_mask)
    assert list(sig.parameters) == ["waves", "eps"]
    with pytest.raises(TypeError):
        dbp._graded_mask(np.array([6705.1169]))          # no EP -> cannot even be called


def test_the_widening_is_paired_with_EP_not_standalone(dbp):
    """Sequence guard. 30 mA is only safe BECAUSE EP is required with it."""
    assert dbp._GRADED_MASK_TOL_A == 0.030
    assert dbp._GRADED_MASK_EP_TOL_EV <= 0.05
    src = (ROOT / "scripts" / "derive_band_products.py").read_text()
    i = src.index("def _graded_mask")
    body = src[i:i + 2000]
    assert "require_ep=True" in body, "the widened window lost its EP requirement"
    assert "line_match" in body, "must route through the one canonical matcher (RYA-1037)"


# ── the recovery, and what it must NOT recover ───────────────────────────────

@pytest.fixture(scope="module")
def lab():
    cg = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    return cg[cg.gf_tier.astype(str).str.contains("LAB", na=False)].reset_index(drop=True)


def _match(lab, wl, ep, tol=0.030, ep_tol=0.05):
    return line_match.match(np.array([wl]), lab.wavelength_air_A.astype(float).values,
                            want_ep=np.array([ep]),
                            src_ep=lab.excitation_potential_eV.astype(float).values,
                            tol_A=tol, ep_tol_eV=ep_tol, require_ep=True)


def test_fe1_6705_1169_recovers_and_it_is_a_ruffoni_line(lab):
    """The ticket's headline. 15.9 mA from its canonical row with the EP agreeing to
    0.0000 eV — same-transition evidence, not a looser window."""
    r = _match(lab, 6705.1169, 4.6070)
    assert r.resolved[0], "6705.1169 no longer recovers"
    row = lab.iloc[int(r.index[0])]
    assert "Ruffoni" in str(row.loggf_reference)
    assert abs(float(row.excitation_potential_eV) - 4.6070) < 1e-3
    assert 0.005 < abs(float(row.wavelength_air_A) - 6705.1169) <= 0.030, (
        "the whole point is that it is FURTHER than the old 5 mA window")


def test_it_is_NOT_recovered_at_the_old_5mA_window(lab):
    """The defect, pinned: at 5 mA this line has no candidate at all, which is why it was
    published ungraded."""
    assert not _match(lab, 6705.1169, 4.6070, tol=0.005).resolved[0]


@pytest.mark.parametrize("wl,ep,d_ep", [(6858.1396, 3.6025, 1.006),
                                        (8713.1976, 2.9488, 2.039),
                                        (8876.0059, 4.5844, 0.436)])
def test_a_lambda_coincidence_is_REFUSED_however_close_in_wavelength(lab, wl, ep, d_ep):
    """🔴 THESE ARE THE ONES THE WIDENING MUST NOT LET IN. Each sits inside 30 mA of a LAB
    row and would have been "recovered" by a wavelength-only widening — while its EP is
    0.4-2.0 eV away, i.e. a different transition. Measured: they DID come through until the
    matcher's lone-candidate hole was closed."""
    r = _match(lab, wl, ep)
    assert not r.resolved[0], (
        f"{wl} was accepted against a row {d_ep} eV away — the EP test is not being applied")


def test_the_matchers_lone_candidate_hole_stays_closed():
    """🔴 THE DEFECT FOUND WHILE DOING THIS TICKET. `match()` applied the EP filter only when
    there was MORE THAN ONE candidate, so a lone row in the wavelength window was accepted
    with no EP check at all — necessary-but-not-sufficient. That is the Fe I 6065.490 shape
    exactly: one candidate, wrong level, no ambiguity flag to warn anyone."""
    r = line_match.match([6705.1169], [6705.1010], want_ep=[4.6070], src_ep=[9.9999],
                         tol_A=0.030, ep_tol_eV=0.05, require_ep=True)
    assert not r.resolved[0], "a lone candidate is again being accepted without an EP check"
    assert r.unresolved, "and the miss must be REPORTED, not silently dropped"


def test_without_strict_mode_existing_callers_are_unchanged():
    """Additive. RYA-1033's callers must not have their match sets moved by this ticket."""
    r = line_match.match([6705.1169], [6705.1010], want_ep=[4.6070], src_ep=[9.9999],
                         tol_A=0.030, ep_tol_eV=0.05)
    assert r.resolved[0], "default behaviour changed — existing products would move"


def test_the_5mA_clusters_that_disagree_on_gf_are_separated_by_EP(lab):
    """The ticket's guard, measured on the real table rather than asserted: Fe I lines
    within 5 mA of each other that disagree on gf are separable by EP, which is what makes
    the 30 mA window safe."""
    cg = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    fe = cg[cg.species.astype(str) == "Fe I"].sort_values("wavelength_air_A")
    w = fe.wavelength_air_A.astype(float).values
    e = fe.excitation_potential_eV.astype(float).values
    g = fe.log_gf.astype(float).values
    clusters = separable = 0
    i = 0
    while i < len(w) - 1:
        j = i
        while j + 1 < len(w) and w[j + 1] - w[j] <= 0.005:
            j += 1
        if j > i:
            clusters += 1
            if np.nanmax(g[i:j+1]) - np.nanmin(g[i:j+1]) > 0.1:
                if np.nanmax(e[i:j+1]) - np.nanmin(e[i:j+1]) > 0.05:
                    separable += 1
            i = j + 1
        else:
            i += 1
    assert clusters > 300
    assert separable > 300, (
        "the clusters that disagree on gf are no longer separable by EP — the 30 mA "
        "window's safety argument has changed and must be re-derived")

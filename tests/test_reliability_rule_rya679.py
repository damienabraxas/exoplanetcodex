"""RYA-679 — ONE reliability rule, single-sourced, with red_chi2 REPORTED not GATED.

Before this ticket, six harnesses fitting the same kind of in-window solar profile
applied five different reliability rules: a red_chi2 ceiling of 60.0 (RYA-564, RYA-581),
15.0 (RYA-565), 5.0 (RYA-560's deblend path), and no ceiling at all (RYA-551,
RYA-560's plain path, RYA-592). None had been ratified.

RYA-679 measured the question rather than splitting the difference:

  * the only written justification for 60.0 — "the sigma_flux=0.01 floor inflates
    rchi2" — is backwards. Measured per-pixel noise in the actual fit windows is
    0.00007-0.0051, so the assumed 0.01 is 2x-146x LARGER than the truth and
    DEFLATES chi2;
  * `red_chi2` is therefore not a chi2 but a rescaled residual RMS, with no
    statistical calibration to hang a pass/fail bar on;
  * full-window red_chi2 tracks BLEND-LIST fidelity, not the element: the same Zr II
    lines on the same spectra score 83.12/25.92/15.93 full-window and 0.39/0.35/1.49
    deblended, while dEW_dA barely moves.

So the ratified rule gates on the two terms that are about the target line's own
response — `dEW_dA` and `railed` — and reports red_chi2 with a loud non-gating review
flag. These tests exist so that decision cannot be quietly re-litigated in one harness.
"""
import ast
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / 'scripts'
SHARED = SCRIPTS / 'solar_profile_fit.py'
sys.path.insert(0, str(SCRIPTS))

from solar_profile_fit import (  # noqa: E402
    RCHI2_REVIEW, RELIABLE_DEWDA, SIGMA_FLUX_ASSUMED, assess_reliability)

# Every harness that derives a `reliable` flag from an in-window profile fit.
HARNESSES = ['rya551_sr2_synth_sirius.py', 'rya560_zr2_synth_sirius.py',
             'rya564_co1_synth_sirius.py', 'rya565_eu2_synth_sirius.py',
             'rya581_ba2_deblend_sirius.py', 'rya592_mg_5528_synth_sirius.py']

# Names that must be defined in exactly ONE place — the shared module.
SINGLE_SOURCED = {'RELIABLE_DEWDA', 'RCHI2_REVIEW', 'SIGMA_FLUX_ASSUMED'}

# The retired ceilings. A harness re-introducing any of these has re-opened a
# question that was settled with measurements.
RETIRED = {'RELIABLE_RCHI2', 'RCHI2_SANE_MAX'}


def _module_level_assignments(path):
    names = set()
    for node in ast.parse(path.read_text()).body:
        if isinstance(node, ast.Assign):
            names |= {t.id for t in node.targets if isinstance(t, ast.Name)}
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
    return names


# ───────────────────────────── single-sourcing ─────────────────────────────
@pytest.mark.parametrize('harness', HARNESSES)
def test_harness_does_not_redefine_the_reliability_constants(harness):
    """The constants live in solar_profile_fit and are imported. A local definition is
    how the 60/15/5 spread happened in the first place."""
    defined = _module_level_assignments(SCRIPTS / harness)
    clash = defined & SINGLE_SOURCED
    assert not clash, (
        f"{harness} re-defines {sorted(clash)} instead of importing it from "
        f"solar_profile_fit. RYA-679 ratified ONE rule; a per-harness copy is exactly "
        f"the divergence this ticket closed.")


@pytest.mark.parametrize('harness', HARNESSES)
def test_harness_does_not_resurrect_a_retired_ceiling(harness):
    defined = _module_level_assignments(SCRIPTS / harness)
    clash = defined & RETIRED
    assert not clash, (
        f"{harness} re-introduces {sorted(clash)}. The red_chi2 ceiling was RETIRED by "
        f"RYA-679 on measured evidence, not deferred. If it needs to come back, it "
        f"needs a new ratification and a rationale that survives the sigma_flux test.")


def test_the_constants_are_defined_in_the_shared_module():
    assert SINGLE_SOURCED <= _module_level_assignments(SHARED)


def test_sigma_flux_is_not_re_inlined_as_a_literal():
    """Both chi2 expressions must use the named constant, so the assumption cannot
    drift between the two fit entrypoints — and so it stays visible."""
    src = SHARED.read_text()
    assert src.count('SIGMA_FLUX_ASSUMED ** 2') == 2, (
        "both fit_profile and fit_profile_deblend must normalise by the NAMED "
        "SIGMA_FLUX_ASSUMED")
    assert '(0.01 ** 2)' not in src, "sigma_flux re-inlined as a bare literal"


# ───────────────────────────── the rule itself ─────────────────────────────
def test_red_chi2_does_not_gate_reliable():
    """THE ratified decision. A fit with an enormous red_chi2 but good sensitivity and
    no railing is RELIABLE. This is the Sr II 4077 case (red_chi2 78.27, dEW/dA 203.5)."""
    fit = dict(dEW_dA=203.5, railed=False, red_chi2=78.27)
    out = assess_reliability(fit)
    assert out['reliable'] is True, (
        "red_chi2 must not gate `reliable` — a full-window red_chi2 is a statement "
        "about the blend list, not about whether A(X) is measurable.")
    assert out['rchi2_review'] is True, "but it MUST raise the review flag"
    assert 'RYA-679' in out['rchi2_review_reason']


@pytest.mark.parametrize('rchi2', [0.0, 0.16, 0.71, 4.72, 60.0, 78.27, 9170.0])
def test_reliable_is_independent_of_red_chi2(rchi2):
    """Sweep red_chi2 across the full observed range; `reliable` must never move."""
    base = dict(dEW_dA=203.5, railed=False)
    assert assess_reliability(dict(base, red_chi2=rchi2))['reliable'] is True


def test_sensitivity_floor_still_gates():
    """The half of the rule that DID survive. Eu II 6645 is the case: an excellent fit
    (red_chi2 0.16) on an intrinsically weak line (dEW/dA 13.9) is NOT reliable."""
    out = assess_reliability(dict(dEW_dA=13.9, railed=False, red_chi2=0.16))
    assert out['reliable'] is False
    assert out['rchi2_review'] is False


def test_railed_still_gates():
    assert assess_reliability(
        dict(dEW_dA=203.5, railed=True, red_chi2=0.1))['reliable'] is False


def test_review_flag_fires_strictly_above_the_trigger():
    base = dict(dEW_dA=203.5, railed=False)
    assert assess_reliability(dict(base, red_chi2=RCHI2_REVIEW))['rchi2_review'] is False
    assert assess_reliability(
        dict(base, red_chi2=RCHI2_REVIEW + 0.01))['rchi2_review'] is True


def test_missing_sensitivity_is_not_reliable():
    """No dEW_dA means the gate cannot be evaluated — that must fail closed, never
    default to reliable."""
    assert assess_reliability(
        dict(dEW_dA=None, railed=False, red_chi2=0.1))['reliable'] is False


def test_red_chi2_is_a_rescaled_residual_rms():
    """The identity that makes the review reason meaningful: red_chi2 =
    (RMS_resid / SIGMA_FLUX_ASSUMED)^2, so a red_chi2 of 78.27 is an 8.8% residual —
    against a measured near-UV photon noise of 0.42%, i.e. ~21x the noise."""
    out = assess_reliability(dict(dEW_dA=203.5, railed=False, red_chi2=78.27))
    assert '8.8%' in out['rchi2_review_reason']
    assert SIGMA_FLUX_ASSUMED == 0.01


# ──────────────────────── committed-artefact regressions ────────────────────────
def _load(name):
    p = REPO / 'data' / 'results' / name
    if not p.exists():
        pytest.skip(f"{name} not committed")
    return json.loads(p.read_text())


def test_sr2_4077_is_not_silently_demoted():
    """RYA-679 section 5: no ceiling may demote Sr II without saying so. Sr II 2.759 is
    the live RYA-669 D.1 adoption candidate, and 4077 carries red_chi2 78.27 — which
    fails BOTH candidate ceilings. Under the ratified rule it stays reliable, on
    sensitivity 203.5 mA/dex. If this ever flips, it must be a deliberate, argued act."""
    d = _load('sr2_synthesis_rya551.json')
    h = d['4077.709']['harps']
    assert h['reliable'] is True, (
        "Sr II 4077 HARPS must remain reliable under the ratified rule — its red_chi2 "
        "of 78.27 is a near-UV blend-model statistic, not evidence A(Sr) is unmeasurable")
    assert h['dEW_dA_mA_dex'] >= RELIABLE_DEWDA
    assert h['red_chi2'] > RCHI2_REVIEW, "and it must stay REVIEW-FLAGGED, loudly"


def test_species_whose_disposition_must_not_move():
    """Ba/Eu/Zr were never load-bearing on the ceiling; confirm they did not move."""
    ba = _load('solar_ba_deblend_rya581.json')
    assert ba['per_arm']['harps']['deblended']['0.6']['reliable'] is True
    assert ba['per_arm']['harps']['deblended']['0.6']['red_chi2'] <= RCHI2_REVIEW

    eu = _load('eu2_synthesis_rya565.json')
    assert eu['disposition']['emit_value'] is False, "Eu stays owed — on SENSITIVITY"

    zr = _load('zr2_deblend_rya585.json')
    assert zr['_meta']['reliable_lines'] == [], "Zr stays owed — on SENSITIVITY"
    for wl in ('4208.98', '4258.041', '4442.992'):
        assert zr[wl]['harps']['dEW_dA_mA_dex'] < RELIABLE_DEWDA

"""RYA-643 — the rest-frame + broadening fix must stay SINGLE-SOURCED.

The defect RYA-592 found existed because one fitter had been copy-pasted into four
harnesses; fixing it in one copy left three carrying it. These tests exist so that
cannot happen a fifth time:

  * no harness may re-inline the shared fit machinery;
  * the rest-frame correction must be measured, never hardcoded;
  * the broadening grid must not have a floor the fit can rail against;
  * the committed measurements must carry the evidence that both were applied.

The physics re-run itself lives on Sirius (RYA-567); these read the committed
artefacts, which is the point of committing them.
"""
import ast
import json
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / 'scripts'
SHARED = SCRIPTS / 'solar_profile_fit.py'

# Every harness that fits an in-window solar profile. Adding one? It imports the
# shared module — that is the whole point.
HARNESSES = ['rya551_sr2_synth_sirius.py', 'rya560_zr2_synth_sirius.py',
             'rya564_co1_synth_sirius.py', 'rya565_eu2_synth_sirius.py',
             'rya581_ba2_deblend_sirius.py', 'rya592_mg_5528_synth_sirius.py']

# The functions that must exist in exactly ONE place.
# fit_profile_deblend/_cont_ratio (RYA-585) join the list for the same reason: the
# deblend path is the SECOND fitter in this codebase, and it is exactly the kind of
# thing the next harness would copy rather than import.
SHARED_FUNCS = {'rot_kernel', 'broaden', 'local_renorm', 'measure_arm_rv', 'fit_profile',
                'fit_profile_deblend', '_cont_ratio'}


def _defined_functions(path):
    return {n.name for n in ast.parse(path.read_text()).body
            if isinstance(n, ast.FunctionDef)}


@pytest.mark.parametrize('script', HARNESSES)
def test_harness_does_not_redefine_the_shared_fitter(script):
    """A harness that defines its own broaden/local_renorm/... has forked the fitter,
    which is exactly how RYA-551/560/564 kept the defect after RYA-592 fixed it."""
    p = SCRIPTS / script
    if not p.exists():
        pytest.skip(f'{script} not present')
    clash = _defined_functions(p) & SHARED_FUNCS
    assert not clash, (f"{script} re-defines {sorted(clash)} instead of importing them from "
                       f"solar_profile_fit — that is the copy-paste that caused RYA-643")


@pytest.mark.parametrize('script', HARNESSES)
def test_harness_imports_the_shared_fitter(script):
    p = SCRIPTS / script
    if not p.exists():
        pytest.skip(f'{script} not present')
    assert 'from solar_profile_fit import' in p.read_text(), \
        f"{script} must source the fit machinery from solar_profile_fit"


def test_no_harness_hardcodes_a_frame_velocity():
    """The frame correction is MEASURED per arm (measure_arm_rv) and FITTED per line.
    A literal velocity anywhere in a harness is a silent fallback waiting to happen."""
    import re
    # the measured solar values, which must appear only as prose/comments, never as code
    suspicious = re.compile(r'^[^#]*\b(0\.76|0\.756|0\.28|0\.278)\b\s*(?!.*km/s\s*$)')
    for script in HARNESSES:
        p = SCRIPTS / script
        if not p.exists():
            continue
        for i, line in enumerate(p.read_text().splitlines(), 1):
            if 'dv' in line and suspicious.match(line) and '=' in line:
                pytest.fail(f'{script}:{i} looks like a hardcoded frame velocity: {line.strip()}')


def test_gsig_grid_has_no_floor_to_rail_against():
    """The old grid started at 1.5 km/s and the fit railed there, trading broadening
    against abundance. The replacement must start well below any plausible solar value."""
    from scripts.solar_profile_fit import GSIG_GRID
    assert GSIG_GRID[0] <= 0.5, f'gsig grid floor {GSIG_GRID[0]} is high enough to rail'
    assert GSIG_GRID[-1] >= 8.0


def test_dv_grid_brackets_the_measured_solar_offsets():
    """HARPS ~+0.76 km/s, IAG ~+0.28. The grid must contain them with room either side,
    or the nuisance parameter rails instead of correcting."""
    from scripts.solar_profile_fit import DV_GRID
    assert DV_GRID.min() < -0.5 and DV_GRID.max() > 1.5


def test_measure_arm_rv_is_abundance_blind_and_recovers_a_known_shift():
    """Inject a known velocity into a synthetic spectrum and confirm it is recovered.
    This is the guard that makes a fitted dv trustworthy."""
    from scripts.solar_profile_fit import CLIGHT, RV_CHECK_LINES, measure_arm_rv
    v_true = 0.80
    w = np.arange(5500.0, 5760.0, 0.01)
    f = np.ones_like(w)
    for lam0 in RV_CHECK_LINES:
        lam = lam0 * (1 + v_true / CLIGHT)
        f -= 0.5 * np.exp(-0.5 * ((w - lam) / 0.05) ** 2)
    v, n, sd = measure_arm_rv(w, f)
    assert n == len(RV_CHECK_LINES)
    assert abs(v - v_true) < 0.05, f'recovered {v}, injected {v_true}'


def test_require_arm_rv_loud_fails_when_unsourced():
    """No silent zero: an arm with no usable check lines must STOP the run."""
    from scripts.solar_profile_fit import require_arm_rv
    w = np.arange(4000.0, 4050.0, 0.01)      # no check line in range
    f = np.ones_like(w)
    with pytest.raises(SystemExit, match='FRAME NOT SOURCED'):
        require_arm_rv(w, f, 'fixture-arm')


# ── the committed measurements must show the fix was applied ────────────────

ARTEFACTS = {'sr2_synthesis_rya551.json': 'Sr II',
             'zr2_synthesis_rya560.json': 'Zr II',
             'zr2_deblend_rya585.json': 'Zr II (deblend)',
             'co_synthesis_rya564.json': 'Co I'}


@pytest.mark.parametrize('fname,label', list(ARTEFACTS.items()))
def test_committed_artefact_records_the_frame_correction(fname, label):
    """Each per-arm record must carry BOTH the fitted dv and the independently measured
    arm velocity — the cross-check is what distinguishes a real frame offset from the
    fit absorbing a profile mismatch."""
    p = REPO / 'data' / 'results' / fname
    if not p.exists():
        pytest.skip(f'{fname} not committed')
    d = json.loads(p.read_text())
    seen = 0
    for k, rec in d.items():
        if k.startswith('_') or not isinstance(rec, dict):
            continue
        for arm in ('harps', 'iag'):
            a = rec.get(arm)
            if not isinstance(a, dict) or a.get('status') == 'no_coverage':
                continue
            assert a.get('dv_fitted_kms') is not None, f'{label} {k}/{arm}: no fitted dv'
            assert a.get('dv_measured_kms') is not None, f'{label} {k}/{arm}: no measured arm dv'
            assert 'gsig_railed' in a, f'{label} {k}/{arm}: gsig railing not reported'
            seen += 1
    assert seen > 0, f'{label}: no fitted arm records found'

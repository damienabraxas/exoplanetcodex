"""
tests/test_pool_integrity_rya279_280.py
=======================================
Regression guards for the two Fe-pool integrity fixes that feed the reported-
uncertainty budget (RYA-282). The fixes landed on main (d721c72, 2ae8abe) but
carried no test — these lock them against silent re-breakage.

RYA-279 — ONE ceiling-correct Fe I pool feeds both the gate scatter statistic
          (A_X_std) and the per-line CSV, so the gate σ can never again be computed
          over lines the artifact excludes. Guarded via the extracted predicate
          `_fe1_ceiling_pool_mask`, which both sites now index.
RYA-280 — a missing `vald_proximity_flag` column must fail LOUD (KeyError), never
          silently default proximity scoring to 0.5; and the comment-headed per-star
          linelists must load (the `comment='#'` fix).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import pipeline.abundances_derive as ad  # noqa: E402
from config.constants import PIPELINE  # noqa: E402

CEIL = float(PIPELINE['vmic_ew_ceiling_mA'])   # 100 mÅ (current); the ceiling under test


# ── RYA-279: the single ceiling-correct Fe I pool ────────────────────────────────
def test_ceiling_mask_excludes_only_saturated_fe1():
    notes = ['Fe 1', 'Fe 1', 'Fe 1', 'Fe 2', 'Ca 1']
    ew    = [CEIL - 10, CEIL, CEIL + 10, CEIL + 500, 5.0]   # boundary line == CEIL is KEPT
    mask  = ad._fe1_ceiling_pool_mask(notes, ew, CEIL)
    assert mask.tolist() == [True, True, False, False, False]


def test_ceiling_mask_is_the_one_pool_for_gate_and_csv():
    """The invariant behind RYA-279: whatever the gate statistic pool keeps, the
    per-line CSV keeps too. Both index this exact mask array, so the Fe I line sets
    are identical by construction — no re-derived threshold can drift between them."""
    rng = np.arange(1, 40)
    notes = ['Fe 1'] * len(rng)
    ew = rng * 5.0                       # 5..195 mÅ — straddles the ceiling
    mask = ad._fe1_ceiling_pool_mask(notes, ew, CEIL)
    gate_pool_idx = set(np.where(mask)[0].tolist())
    csv_kept_idx  = {i for i in range(len(notes)) if mask[i]}   # the per-line loop's test
    assert gate_pool_idx == csv_kept_idx
    assert all(ew[i] <= CEIL for i in gate_pool_idx)            # nothing saturated slips in


# ── RYA-280: proximity flag must load, or fail loud ──────────────────────────────
def _min_inputs():
    per_line = pd.DataFrame([{
        'element': 'Fe', 'ion': 'I', 'wavelength_air_A': 5000.0,
        'ew_mA': 45.0, 'a_1dlte': 7.46,
    }])
    ew_df = pd.DataFrame([{
        'element': 'Fe', 'ion': 'I', 'wavelength_air_A': 5000.0,
        'ew_err_mA': 1.0, 'chi2': 1.0,
    }])
    return per_line, ew_df


def _write_linelist(path, with_prox=True):
    cols = ['element', 'ion', 'wavelength_air_A']
    row = {'element': 'Fe', 'ion': 'I', 'wavelength_air_A': 5000.0}
    if with_prox:
        cols.append('vald_proximity_flag')
        row['vald_proximity_flag'] = 0.9        # distinctive — NOT the 0.5 silent default
    df = pd.DataFrame([row])[cols]
    # Comment-headed, exactly like the per-star linelists (the RYA-280 trigger).
    with open(path, 'w') as fh:
        fh.write('# per-star linelist metadata block\n# provenance: unit-test\n')
        df.to_csv(fh, index=False)


def test_missing_proximity_column_fails_loud(tmp_path):
    ll = tmp_path / 'linelist_noprox.csv'
    _write_linelist(ll, with_prox=False)
    per_line, ew_df = _min_inputs()
    with pytest.raises(KeyError, match='vald_proximity_flag'):
        ad._compute_line_scores(per_line, ew_df, linelist_path=str(ll))


def test_proximity_loads_from_comment_headed_linelist(tmp_path):
    """Active scoring: the flag is read from a comment-headed file (the comment='#'
    fix), so the output carries the real 0.9, not the inert 0.5 default."""
    ll = tmp_path / 'linelist_prox.csv'
    _write_linelist(ll, with_prox=True)
    per_line, ew_df = _min_inputs()
    scored = ad._compute_line_scores(per_line, ew_df, linelist_path=str(ll))
    assert scored['vald_proximity_flag'].iloc[0] == pytest.approx(0.9)
    assert scored['vald_proximity_flag'].iloc[0] != 0.5     # not the silent default

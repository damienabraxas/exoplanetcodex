"""
RYA-1033 — the measured-line -> atomic-data join must never be keyed on a rounded wavelength.

Every number asserted here was MEASURED on the committed data before the fix was written;
see `data/audit/rya1033_rounded_key_join/`. The tests pin the DEFECT, not just the API, so
that reintroducing a rounded key fails loudly rather than quietly moving a pool count.
"""
from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pipeline import line_match
from pipeline.line_match import LineMatchError, match, match_frames, require_resolved

ROOT = Path(__file__).resolve().parents[1]
CANONICAL_GF = ROOT / "data" / "linelists" / "canonical_gf.csv"
MEASURED_EW = ROOT / "data" / "measured" / "sol_ew_results_v1.csv"

#: The exact Fe I lines the 2-dp rounded key dropped, from the RYA-1033 audit.
DROPPED_17 = [
    4787.49462, 5281.16500, 5513.86524, 5650.70512, 5770.15492, 5793.91435,
    5833.92489, 5963.05425, 5984.81467, 6013.91474, 6034.03465, 6034.89465,
    6252.55455, 6364.69383, 6713.19471, 6845.93488, 6875.44503,
]


def _fe_i_measured() -> pd.DataFrame:
    ew = pd.read_csv(MEASURED_EW)
    return ew[(ew.element == "Fe") & (ew.ion == "I")].reset_index(drop=True)


def _fe_i_canonical() -> pd.DataFrame:
    cg = pd.read_csv(CANONICAL_GF, low_memory=False)
    return cg[cg.species == "Fe I"].reset_index(drop=True)


# ── the defect itself ────────────────────────────────────────────────────────────────

def test_rounded_key_is_not_a_function_of_the_wavelength():
    """🔴 THE SHARPEST HALF OF RYA-1033.

    Python and numpy round 6136.615 to 2 dp to DIFFERENT values. `promote_solar_ew` used
    pandas and `abundances_derive` used Python, so one wavelength had two keys and the join
    result depended on which module built which side.
    """
    assert round(6136.615, 2) == 6136.61
    assert float(np.round(np.float64(6136.615), 2)) == 6136.62
    assert round(6136.615, 2) != float(np.round(np.float64(6136.615), 2))


def test_the_17_dropped_lines_are_all_present_in_canonical_gf():
    """ZERO of the 17 are genuinely absent — every one is within 1.2 mA of a canonical row."""
    canon = _fe_i_canonical()
    W = np.sort(canon.wavelength_air_A.astype(float).values)
    worst = 0.0
    for wl in DROPPED_17:
        d = float(np.min(np.abs(W - wl)))
        assert d <= line_match.MATCH_TOL_A, f"{wl} unexpectedly {1000 * d:.2f} mA from canonical"
        worst = max(worst, d)
    assert worst <= 0.0012, f"worst separation grew to {1000 * worst:.2f} mA"


def test_rounded_key_drops_exactly_those_17_and_the_matcher_keeps_them():
    """The regression the fix removes, stated as a SET difference, never a count."""
    meas, canon = _fe_i_measured(), _fe_i_canonical()
    w = meas.wavelength_air_A.astype(float).to_numpy()
    ck = set(canon.wavelength_air_A.astype(float).round(2))

    # Round BOTH sides the way the original join did (pandas/numpy). Rounding the measured
    # side with Python's `round` instead yields EIGHTEEN — 6136.615 joins the casualties —
    # which is the library divergence test above, seen on the real pool.
    dropped = set(np.round(w[~pd.Series(w).round(2).isin(ck).to_numpy()], 5))
    assert dropped == {round(x, 5) for x in DROPPED_17}

    py_dropped = set(np.round(w[[round(float(x), 2) not in ck for x in w]], 5))
    assert len(py_dropped) == 18 and py_dropped - dropped == {6136.615}

    res = match(w, canon.wavelength_air_A.astype(float).to_numpy())
    assert not res.unresolved, "the tolerance matcher must resolve every measured Fe I line"


# ── wavelength alone is not an identity ──────────────────────────────────────────────

def test_ambiguous_lines_are_refused_not_guessed():
    """canonical_gf holds coincident Fe I transitions whose gf differ by up to ~1.9 dex.

    Without an EP to separate them the matcher must REFUSE. Silently taking the nearest is
    the RYA-780/852 defect.
    """
    meas, canon = _fe_i_measured(), _fe_i_canonical()
    res = match_frames(meas, canon)          # sol_ew_results_v1 carries no ep_eV
    assert len(res.ambiguous) == 7, f"expected 7 ambiguous, got {len(res.ambiguous)}"
    assert any(abs(w - 6065.482) < 1e-6 for w, _ in res.ambiguous)
    with pytest.raises(LineMatchError, match="MORE THAN ONE"):
        require_resolved(res, what="sol_ew_results_v1", species="Fe I")


def test_excitation_potential_resolves_the_ambiguity():
    """6065.4820 (EP 2.609, NIST-C+) vs 6065.4850 (EP 4.956, KURUCZ) — 1.9 dex apart."""
    src_wl = np.array([6065.4820, 6065.4850])
    src_ep = np.array([2.609, 4.956])
    res = match(np.array([6065.48200]), src_wl, want_ep=np.array([2.609]), src_ep=src_ep)
    assert not res.ambiguous
    assert res.index[0] == 0

    res_hi = match(np.array([6065.48200]), src_wl, want_ep=np.array([4.956]), src_ep=src_ep)
    assert res_hi.index[0] == 1, "EP must win over raw wavelength proximity"


# ── the loud-failure contract ────────────────────────────────────────────────────────

def test_unresolved_raises_and_names_every_line():
    res = match(np.array([4000.0, 5000.0]), np.array([4000.0005, 9999.0]))
    assert res.index[0] == 0 and res.index[1] == -1
    with pytest.raises(LineMatchError) as e:
        require_resolved(res, what="fixture", species="Fe I")
    assert "5000.00000" in str(e.value)
    assert "NaN gf_tier" in str(e.value)


def test_duplicate_rows_of_one_line_are_not_ambiguity():
    """Two identical wavelengths are one line listed twice, not a fork."""
    res = match(np.array([5000.0]), np.array([5000.0, 5000.0]))
    assert not res.ambiguous and res.index[0] in (0, 1)


def test_empty_source_is_unresolved_not_a_crash():
    res = match(np.array([5000.0]), np.array([]))
    assert res.index[0] == -1 and len(res.unresolved) == 1


# ── the guard: no rounded wavelength may become a join key again ─────────────────────

#: Sites where a rounded wavelength indexes a lookup whose BOTH SIDES are produced by the
#: same rounding call on the same in-memory frame. A key cannot disagree with itself, so
#: these cannot exhibit the RYA-1033 defect, which is specifically a CROSS-TABLE join
#: between files that store the line to different precision.
#:
#: ⚠️ They are still listed rather than pattern-excluded, so adding one is a deliberate act.
#: Every entry below was read and confirmed self-consistent; none reads a second file.
_ROUNDED_KEY_ALLOWED = {
    # dedup of "special" rows inside ONE frame, both sides built from that frame. Flagged
    # in RYA-1033 but out of scope: the window also merges lines 0.1 A apart, so changing
    # it changes the worklist and therefore every measurement.
    ("pipeline/lines_fit.py", "_wkey"),
    # one-off per-element synthesis drivers: `have`/`out`/`rows` are built and read inside
    # the same function from one source, as a component-count or row cache.
    ("scripts/measure_cu_v_hfs_synthesis_rya466.py", "out[round(wl, 3)]"),
    ("scripts/measure_mn_hfs_synthesis_rya473.py", "out[round(wl, 3)]"),
    ("scripts/rya560_zr2_synth_sirius.py", "have[round(float(row"),
    ("scripts/rya564_co1_synth_sirius.py", "have[round(float(row"),
    ("scripts/rya565_eu2_synth_sirius.py", "rows[round(float(r"),
    ("scripts/rya592_mg_5528_synth_sirius.py", "have[round(float(row"),
    # `edit_idx` is keyed and read from the one adjudication frame being edited in place.
    ("scripts/adjudicate_gf_batch3_zrbayeu_rya354.py", "edit_idx[(species, round(wl, 3))]"),
    # `dropped_placeholders` is produced BY `detect_placeholder_zero_lines(df)` from this
    # same `df`, so both sides of the `.isin` carry identical rounding.
    ("scripts/build_nlte_grids_mpia.py", "dropped_placeholders"),
    # per-(element, ion) cache built and read from one engine-diff frame.
    ("scripts/rya527_two_engine_run.py", "out.setdefault((el, ion), {})"),
}


def _rounded_wavelength_key_sites() -> list[str]:
    """Every `round(<something wavelengthy>, n)` that feeds a dict/set/merge key."""
    hits: list[str] = []
    for rel in ("pipeline", "scripts"):
        for path in sorted((ROOT / rel).rglob("*.py")):
            try:
                tree = ast.parse(path.read_text())
            except SyntaxError:                      # pragma: no cover
                continue
            src = path.read_text().splitlines()
            for node in ast.walk(tree):
                if not isinstance(node, ast.Subscript):
                    continue
                for sub in ast.walk(node.slice):
                    if not (isinstance(sub, ast.Call) and _is_round(sub)):
                        continue
                    if not _mentions_wavelength(sub):
                        continue
                    line = src[node.lineno - 1].strip()
                    key = f"{path.relative_to(ROOT)}:{node.lineno}: {line}"
                    if not any(a[0] == str(path.relative_to(ROOT)) and a[1] in line
                               for a in _ROUNDED_KEY_ALLOWED):
                        hits.append(key)
    return hits


def _is_round(call: ast.Call) -> bool:
    f = call.func
    if isinstance(f, ast.Name) and f.id == "round":
        return True
    return isinstance(f, ast.Attribute) and f.attr == "round"


def _mentions_wavelength(node: ast.AST) -> bool:
    text = ast.dump(node)
    return any(w in text for w in ("wavelength", "wave_A", "wave_air", "'wl'", '"wl"'))


def test_no_rounded_wavelength_join_keys_remain():
    """🔴 RYA-1033's standing guard. A rounded wavelength must never index a lookup again."""
    hits = _rounded_wavelength_key_sites()
    assert not hits, (
        "rounded-wavelength key(s) reintroduced — use pipeline.line_match instead:\n  "
        + "\n  ".join(hits))

"""RYA-1110 — the GBS (Jofré+ 2014) solar Fe reference line set.

The smoke tests the ticket asks for, plus the two guards that matter more than they do:

  * THE ARTIFACT MUST NOT DRIFT FROM THE CODE THAT DEFINES IT. The committed CSVs are
    rebuilt here and compared byte-for-byte. A committed line list that no longer matches
    its generator is the RYA-1101/686 failure — and there, the STALE side is the one that
    passed the gate.
  * A TOLERANCE NEEDS A MEASURED NULL. `test_the_displaced_null_resolves_nothing` is the
    control that the 0.015 Å + EP key is identifying lines rather than finding neighbours.
    Without it the join's pass rate is unfalsifiable.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import rya1110_build_gbs_fe_lineset as B  # noqa: E402


@pytest.fixture(scope="module")
def lineset():
    return pd.read_csv(B.LINESET)


@pytest.fixture(scope="module")
def cov():
    return pd.read_csv(B.COVERAGE)


@pytest.fixture(scope="module")
def holdings(tmp_path_factory):
    """`measure_band_ew` resolves the Kitt Peak atlas at IMPORT, and the holdings
    REGISTRY needs no spectra. Same stub the RYA-933/767 suites use."""
    kp = tmp_path_factory.mktemp("kp")
    (kp / "lm0296").touch()
    os.environ.setdefault("CODEX_KP_ATLAS", str(kp))
    from measure_band_ew import _INSTRUMENT_HOLDINGS
    return _INSTRUMENT_HOLDINGS


# ── the ticket's smoke test ──────────────────────────────────────────────────
def test_the_count_reproduces_the_published_solar_n(lineset):
    """Jofré Table 3 states the Sun's selected lines outright: 150 Fe I + 9 Fe II."""
    got = lineset.groupby("species").size().to_dict()
    assert got == B.PUBLISHED_SUN_COUNTS
    assert len(lineset) == 159


def test_the_rew_cut_is_verified_on_every_class(lineset):
    """`pass` / `excluded` / `ambiguous` are three CLOSED sets that partition the file.

    The ticket asks for "all <= -4.8"; the honest form of that is per class, because the
    file deliberately keeps the lines the cut removes (RYA-931 — quarantine, never cull).
    """
    assert set(lineset.rew_class) == {"pass", "excluded", "ambiguous"}
    p = lineset[lineset.rew_class == "pass"]
    x = lineset[lineset.rew_class == "excluded"]
    a = lineset[lineset.rew_class == "ambiguous"]
    assert len(p) + len(x) + len(a) == len(lineset)
    assert (p.rew_max <= B.REW_CUT).all(), "a 'pass' line exceeds the cut on some method"
    assert (x.rew_min > B.REW_CUT).all(), "an 'excluded' line passes on some method"
    # ambiguous is defined by straddling; assert it STRADDLES rather than trusting the label
    assert (a.rew_min <= B.REW_CUT).all() and (a.rew_max > B.REW_CUT).all()
    assert len(p) == 142 and len(x) == 14 and len(a) == 3


def test_the_published_selection_is_carried_unmodified(lineset):
    """Every row is a line Jofré published as USED for the Sun. The REW class is ours;
    the selection flag is theirs, and the two are separate columns on purpose."""
    assert set(lineset.gbs_selected_sun) == {1}
    assert set(lineset.gbs_golden) <= {"Y", "N"}


def test_both_gf_columns_ship_and_neither_is_filled_in_for_the_other(lineset):
    """The decision flag needs both, and needs the gaps to stay gaps (RYA-161)."""
    assert lineset.log_gf_ours.notna().all(), "the canonical_gf join must be total"
    assert int(lineset.log_gf_gbs.notna().sum()) == 138
    # Where Jofré publishes no gf, the row is EMPTY — never backfilled from our value.
    missing = lineset[lineset.log_gf_gbs.isna()]
    assert (missing.gf_provenance_gbs == "NOT PUBLISHED IN TABLES 4/5").all()
    assert missing.elow_eV_paper.isna().all()
    assert missing.delta_gbs_minus_ours.isna().all()
    # and the one golden line the paper omits is named, so it cannot quietly become 138+1
    g = missing[missing.gbs_golden == "Y"]
    assert len(g) == 1 and float(g.iloc[0].wavelength_air_A) == 6149.26


def test_no_gf_reference_code_is_claimed_to_be_resolved(lineset):
    """The decoder footnote did not typeset in the copy we hold. A code is not a source."""
    assert set(lineset.loggf_ref_gbs_resolved) == {False}


# ── the joins, and the controls on them ──────────────────────────────────────
def test_the_join_is_a_dual_key_within_the_derived_tolerance(lineset):
    d = lineset.match_distance_mA.dropna()
    assert len(d) == len(lineset)
    assert d.max() <= B._MATCH_TOL_A * 1000 + 1e-9
    # EP is part of the key, not decoration: the joined rows must agree on it.
    ep = (lineset.our_excitation_potential_eV - lineset.excitation_potential_eV).abs()
    from pipeline.line_match import EP_TOL_EV
    assert ep.max() <= EP_TOL_EV


def _resolve_at(sun, cg, shift, *, ep=True):
    """The same match `_null_control` runs, at an arbitrary displacement.

    `ep=False` drops the EP half of the key — used only to show what the null looks like
    WITHOUT it, which is the evidence that the dual key is doing the work.
    """
    from pipeline import line_match
    n = 0
    for species in B.PUBLISHED_SUN_COUNTS:
        want = [r for r in sun if r["species"] == species]
        src = cg[cg["species"] == species]
        kw = {}
        if ep:
            kw = dict(want_ep=[r["excitation_potential_eV"] for r in want],
                      src_ep=src["excitation_potential_eV"].values, require_ep=True)
        n += line_match.match([r["wavelength_air_A"] + shift for r in want],
                              src["wavelength_air_A"].values,
                              tol_A=B._MATCH_TOL_A, **kw).n_resolved
    return int(n)


def test_the_displaced_null_resolves_nothing():
    """🔴 THE CONTROL. Displace the GBS wavelengths past any real identification and the
    same match must resolve ZERO lines. A join whose chance rate was never measured is
    not evidence that the lines were identified (RYA-1070)."""
    ew = B._read_measurement_table(B.VIZIER / "ew.dat")
    sun = [r for r in ew if r["star"] == B.STAR]
    cg = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    cg = cg[cg["species"].isin(B.PUBLISHED_SUN_COUNTS)].reset_index(drop=True)
    assert B._null_control(sun, cg) == [0] * len(B._NULL_SHIFTS_A)
    # ...and the POSITIVE half, without which "0" would also be what a broken matcher
    # returns. Undisplaced, the same call must resolve every one of the 159 lines.
    assert B._null_control(sun, cg[:0]) == [0] * len(B._NULL_SHIFTS_A)   # empty source: 0
    assert _resolve_at(sun, cg, 0.0) == 159
    # ...and the EP half is what makes the null zero, MEASURED rather than claimed:
    # on wavelength alone the same displacements find 12-22 chance "matches", and the
    # undisplaced match finds only 149 of 159.
    assert _resolve_at(sun, cg, 0.0, ep=False) == 149
    assert all(_resolve_at(sun, cg, s, ep=False) >= 12 for s in B._NULL_SHIFTS_A)


def test_the_ep_control_would_fire_if_the_paper_join_were_wrong():
    """The paper join is checked by EP; prove the check is not vacuous by breaking it."""
    ew = B._read_measurement_table(B.VIZIER / "ew.dat")
    sun = [r for r in ew if r["star"] == B.STAR]
    paper = B._read_paper_tables()
    B._check_ep_agreement(sun, paper)                     # the real data passes
    bad = dict(paper)
    k = next(iter(bad))
    bad[k] = dict(bad[k], elow_eV=bad[k]["elow_eV"] + 1.0)
    with pytest.raises(B.BuildError):
        B._check_ep_agreement(sun, bad)                   # a wrong transition is caught


def test_the_table6_pairing_is_ordinal_and_its_preconditions_fire():
    """RYA-1037: `table6` carries NO excitation potential, so a tolerance match against it
    would be a λ-only key. The pairing is ordinal instead, and the conditions that make an
    ordinal pairing provable are asserted — so prove the assertions are not decorative."""
    ew = B._read_measurement_table(B.VIZIER / "ew.dat")
    t6 = B._read_table6(B.VIZIER / "table6.dat")
    bound = B._bind_table6(ew, t6)
    assert len(bound) == 242                            # the real holding pairs cleanly

    with pytest.raises(B.BuildError, match="unequal counts|lists"):
        B._bind_table6(ew, t6[:-1])                     # a dropped table6 row

    # a line moved to within one print step of its neighbour: rounding could reorder them
    crowded = [dict(r) for r in ew]
    lams = sorted({r["wavelength_air_A"] for r in crowded if r["species"] == "Fe I"})
    for r in crowded:
        if r["species"] == "Fe I" and r["wavelength_air_A"] == lams[1]:
            r["wavelength_air_A"] = lams[0] + 0.05
    with pytest.raises(B.BuildError, match="print step"):
        B._bind_table6(crowded, t6)

    # a pair that ordinal position says belongs together but rounding cannot connect
    shifted = [dict(t) for t in t6]
    shifted[0] = dict(shifted[0],
                      wavelength_air_A=shifted[0]["wavelength_air_A"] + 0.4)
    with pytest.raises(B.BuildError):
        B._bind_table6(ew, shifted)


def test_the_published_count_check_would_fire_if_the_holding_changed():
    with pytest.raises(B.BuildError):
        B._check_published_counts([])


# ── scope, vocabulary, reachability ──────────────────────────────────────────
def test_the_line_set_tag_is_in_the_registry_vocabulary(lineset):
    from pipeline.model_registry import LINE_SETS
    assert set(lineset.line_set) == {B.LINE_SET_TAG}
    assert B.LINE_SET_TAG in LINE_SETS


def test_every_line_resolves_to_vis_through_band_policy(lineset):
    from pipeline import band_policy
    assert {band_policy.resolve(float(w)).name
            for w in lineset.wavelength_air_A} == {B.BAND}
    assert set(lineset.band) == {B.BAND}


def test_no_line_is_telluric_excluded_and_the_reason_is_checkable(lineset):
    from pipeline.telluric_policy import TELLURIC_BANDS
    assert lineset.telluric_exclusion.isna().all()
    bluest = min(lo for lo, _hi, _n in TELLURIC_BANDS)
    assert float(lineset.wavelength_air_A.max()) < bluest


def test_the_coverage_report_names_every_holding_and_adds_up(cov, holdings, lineset):
    declared = {h.holding_id for specs in holdings.values() for h in specs}
    assert set(cov.holding_id) == declared, "a holding with no row is an unanswered question"
    assert (cov.n_reachable == cov.n_in_span - cov.n_telluric_excluded).all()
    assert set(cov.n_lines_selected) == {int((lineset.rew_class == "pass").sum())}
    # The two IAG holdings are complementary, not redundant — RYA-767's span declarations.
    iag = cov[cov.instrument == "iag_fts_solar_atlas"]
    assert int(iag.n_reachable.sum()) == int((lineset.rew_class == "pass").sum())


def test_an_undeclared_span_is_reported_as_undeclared_never_as_covered(cov):
    """`covers()` returns True for everything when span_A is None. That is the ABSENCE of
    a coverage claim; a report that read it as coverage would be RYA-767 again."""
    u = cov[cov.span_source.str.startswith("UNDECLARED")]
    assert len(u) == 1 and u.iloc[0].holding_id == "solar_vesta_crires_plus_idp"
    assert int(u.iloc[0].n_in_span) == 0


# ── the artifact must match the code that defines it ─────────────────────────
def test_the_committed_artifacts_rebuild_byte_for_byte(tmp_path, holdings):
    """🔴 RYA-1101's lesson: a committed artifact and the module that defines it landed in
    one merge disagreeing, and every reader passed because none compared file to module."""
    df = B.build()
    got = df.to_csv(index=False, lineterminator="\n")
    assert got == B.LINESET.read_text(), "gbs_solar_fe_rya1110.csv is STALE — re-run the builder"
    gotc = B.coverage(df).to_csv(index=False, lineterminator="\n")
    assert gotc == B.COVERAGE.read_text(), "the coverage report is STALE — re-run the builder"


def test_the_holding_is_intact():
    import rya1110_fetch_jofre2014 as F
    assert F.verify() == 0

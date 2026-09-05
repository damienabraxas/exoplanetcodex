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

import numpy as np
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


def test_the_gf_reference_codes_are_now_decoded(lineset):
    """SUPERSEDED PREMISE, kept as the record of what changed.

    The first RYA-1110 pass asserted `set(loggf_ref_gbs_resolved) == {False}` — true then,
    because the arXiv copy we held (1309.1099v2) printed the footnote as
    `References: 102: ????????.` The PUBLISHED A&A PDF typesets it, so the codes decode and
    the flag varies. The test is inverted rather than deleted: a guard that recorded a real
    limitation should show when the limitation lifted.
    """
    assert set(lineset.loggf_ref_gbs_resolved) == {True, False}
    assert int(lineset.loggf_ref_gbs_resolved.sum()) == 137
    # every row that carries a code decodes, except the one the paper never defined
    coded = lineset[lineset.loggf_ref_code_gbs.notna()]
    assert int((~coded.loggf_ref_gbs_resolved).sum()) == 1


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


#: 🔴 RYA-1193 ADDED THE O2 GAMMA-BAND AND IT REACHES THREE GBS LINES — A VALUE MOVE THAT
#: TICKET FORBADE ITSELF FROM MAKING.
#:
#: RYA-940's molecfit fitted SEVEN bands; `telluric_policy.TELLURIC_BANDS` listed six.
#: RYA-1193 added the missing one, `O2 gamma 6270-6300 A` (spec item 3), sourced from
#: RYA-940's own `o2gamma/fit_manifest.json`. That is correct on its own terms.
#:
#: But three GBS reference lines sit inside it — 6270.22, 6271.28, 6297.79, all currently
#: `rew_class == pass` — so a rebuild stamps QUARANTINED-TELLURIC on them and the usable
#: GBS Fe set goes 142 -> 139. RYA-1193 is explicitly "config + doc only, no value moved,
#: feed/linelists diff-clean", and `data/linelists/reference_sets/` is a linelist. The
#: spec and the constraint cannot both be honoured, so the POLICY is corrected and the
#: ARTIFACT is left exactly as committed, with the divergence recorded here rather than
#: resolved by whichever of the two I preferred.
#:
#: ⚠️ AND THE UNDERLYING QUESTION IS RYAN'S, NOT A REGENERATION. The exclusion fires at
#: `telluric_basis=unspecified`; on a holding whose basis is `corrected` it would not fire
#: at all. Whether an instrument-AGNOSTIC reference line set should carry an
#: instrument-DEPENDENT exclusion is a science decision, and rebuilding the artifact would
#: have answered it silently in one direction.
#:
#: strict=True so these flip to a hard failure the moment the artifact is rebuilt — the
#: marker cannot outlive the divergence (the RYA-853 idiom).
_GBS_O2GAMMA_REASON = (
    "RYA-1193: O2 gamma-band 6270-6300 A newly enumerated; 3 GBS lines (6270.22, "
    "6271.28, 6297.79) would become telluric-excluded, 142 -> 139 usable. Artifact "
    "deliberately NOT rebuilt — RYA-1193 is no-value-moved. Ryan's call.")


@pytest.mark.xfail(strict=True, reason=_GBS_O2GAMMA_REASON)
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


# ── gf PROVENANCE DECODE (RYA-1110 second pass) ──────────────────────────────
def test_every_line_carries_a_decoded_gf_provenance(lineset):
    """The column is never empty and never says something the basis does not support."""
    assert lineset.gf_source_per_line.notna().all()
    assert (lineset.gf_source_per_line.str.len() > 0).all()
    assert set(lineset.gf_source_basis) == {
        "heiter2021-exact", "jofre2014-footnote", "no-gbs-value", "unresolved"}
    counts = lineset.gf_source_basis.value_counts().to_dict()
    assert counts == {"heiter2021-exact": 93, "jofre2014-footnote": 44,
                      "no-gbs-value": 21, "unresolved": 1}


def test_the_heiter_route_is_used_ONLY_where_the_values_agree(lineset):
    """🔴 THE RULE THAT KEEPS THIS FROM FABRICATING A PEDIGREE.

    Heiter+2021 is GES v6; Jofré used v3. Where the versions carry different log gf,
    Heiter's per-line `r_loggf` is the source of a DIFFERENT NUMBER, and citing it for the
    GBS value would make the line look sourced to a paper whose value it does not carry —
    the `gf_grades` SCALE-MISMATCH defect. So the basis must track the value agreement
    exactly, in both directions.
    """
    d = (lineset.log_gf_gbs - lineset.heiter2021_log_gf).abs()
    exact = lineset.gf_source_basis == "heiter2021-exact"
    assert (d[exact] <= B._GF_SAME_DEX).all(), "a heiter2021-exact row has a different value"
    foot = lineset.gf_source_basis == "jofre2014-footnote"
    assert (d[foot] > B._GF_SAME_DEX).all(), "a footnote row could have used the finer route"
    # and the split is real, not one class swallowing everything
    assert int(exact.sum()) > 0 and int(foot.sum()) > 0


def test_no_gbs_value_rows_claim_no_attribution(lineset):
    n = lineset[lineset.gf_source_basis == "no-gbs-value"]
    assert n.log_gf_gbs.isna().all()
    assert n.gf_source_per_line.str.startswith("NO GBS gf PUBLISHED").all()
    # Heiter still answers for these — the v6 value and its source are recorded.
    assert n.heiter2021_log_gf.notna().all() and n.heiter2021_r_loggf.notna().all()


def test_the_one_undecodable_code_is_named_and_not_guessed(lineset):
    """Jofré's Table 4 body uses code 190; the published footnote never defines it."""
    u = lineset[lineset.gf_source_basis == "unresolved"]
    assert len(u) == 1
    r = u.iloc[0]
    assert (r.species, float(r.wavelength_air_A), int(r.loggf_ref_code_gbs)) == (
        "Fe I", 4985.55, 190)
    assert r.rew_class == "excluded"          # so it is NOT in the 142-line set
    assert "UNRESOLVED" in r.gf_source_per_line
    assert 190 not in B._read_jofre_codes()


def test_the_two_decoders_are_cross_checked_not_merely_stacked(lineset):
    """Where BOTH decoders answer, do they name the same people?

    This is the control on the whole decode. If the λ+EP join to Heiter were wrong, or the
    footnote transcription were wrong, the author sets would disagree at random. They do
    not: on every row where Heiter's value equals the GBS value, the two independently
    sourced author sets overlap.
    """
    codes = B._read_jofre_codes()
    refs = B._read_ges_refs()
    agree = disagree = 0
    for _, r in lineset.iterrows():
        if r.gf_source_basis != "heiter2021-exact" or pd.isna(r.loggf_ref_code_gbs):
            continue
        j = codes.get(int(r.loggf_ref_code_gbs))
        if j is None:
            continue
        _, names = B._decode_ges_code(r.heiter2021_r_loggf, refs)
        if set(names) & set(j[2]):
            agree += 1
        else:
            disagree += 1
    assert agree >= 88, f"only {agree} rows cross-check"
    assert disagree == 0, f"{disagree} rows where the value agrees but the sources do not"


def test_the_firewall_check_fires_on_the_lines_it_should(lineset):
    """🔴 Decoding the provenance is what makes RYA-161 checkable, so the check ships here.

    Jofré code 158 is Meléndez & Barbuy (2009) — `melendez2009` in the bibliography,
    *"partly solar-fitted, so it must never referee a solar abundance"*. Three GBS Fe II
    lines carry it, and all three are in the 142-line replication set.
    """
    f = lineset[lineset.gf_source_firewalled.notna()]
    assert set(f.species) == {"Fe II"}
    assert sorted(f.wavelength_air_A.astype(float)) == [5414.07, 5425.26, 6432.68]
    assert (f.rew_class == "pass").all()
    assert (f.loggf_ref_code_gbs == 158).all()
    assert f.gf_source_firewalled.str.contains("RYA-161").all()
    # non-vacuous: an unfirewalled line must come back clean
    assert lineset[lineset.gf_source_basis == "heiter2021-exact"].gf_source_firewalled.isna().any()


def test_the_ratified_disposition_is_flag_and_keep(lineset):
    """🔴 Ryan, RYA-1110, 2026-08-29. The three solar-fitted Fe II lines STAY.

    Replication fidelity: we replicate Jofré's PUBLISHED set and flag its properties — the
    same principle already applied to the −4.8 quirk. A later pass that "cleans" them out
    would be undoing a decision, so the count is pinned here.
    """
    circ = B.solar_circular_lines(lineset)
    assert len(circ) == 3
    assert (circ.rew_class == "pass").all(), "a flagged line was dropped from the 142"
    assert int((lineset.rew_class == "pass").sum()) == 142, "the set must not shrink"
    assert set(lineset.gbs_solar_validity) <= set(B.SOLAR_VALIDITY)
    assert lineset.gbs_solar_validity.value_counts().to_dict() == {
        "not-flagged": 156, "method-reproduction-only": 3}


def test_the_flag_carries_the_MEANING_not_just_the_paper_name(lineset):
    """The point of Ryan's decision: a flag that only names Meléndez & Barbuy leaves the
    reader to rediscover the consequence. It must say what it means for a solar number."""
    circ = B.solar_circular_lines(lineset)
    for _, r in circ.iterrows():
        t = r.gf_source_firewalled
        assert "METHOD-REPRODUCTION CHECK, NOT AN INDEPENDENT VALIDATION" in t
        assert "FLAG-AND-KEEP" in t and "RYA-1110" in t
        assert "STAYS in the replication set" in t


def test_whether_our_gf_escapes_the_circularity_is_stated_PER_LINE(lineset):
    """Ryan: the circularity does not escape via "use our gf" either — on two of the three
    our adopted value is the same MB09 number. That is a per-line fact, so it is derived
    per line rather than asserted once in prose no filter can see."""
    circ = B.solar_circular_lines(lineset)
    circ = circ.set_index(circ.wavelength_air_A.astype(float))
    same = circ.loc[[5414.07, 5425.26]]
    assert (same.log_gf_gbs == same.log_gf_ours).all()
    assert same.gf_source_firewalled.str.contains("does NOT escape").all()
    diff = circ.loc[[6432.68]]
    assert (diff.log_gf_gbs != diff.log_gf_ours).all()
    assert diff.gf_source_firewalled.str.contains("not circular here").all()
    assert not diff.gf_source_firewalled.str.contains("does NOT escape").any()


def test_the_coverage_report_carries_the_flag_to_every_holding(cov):
    """"Wherever these lines feed a solar number" — the coverage report is the only such
    surface that exists today (RYA-1111 is the measurement path and is not built)."""
    assert "n_reachable_solar_circular" in cov.columns
    assert (cov.n_reachable_solar_circular <= cov.n_reachable).all()
    reach = cov[cov.n_reachable > 0]
    # every holding that reaches the red-optical end of the set reaches all three
    assert set(reach.n_reachable_solar_circular) == {0, 3}
    # and the one that reaches none is the blue-capped IAG arm, not an accident
    zero = reach[reach.n_reachable_solar_circular == 0]
    assert list(zero.holding_id) == ["solar_iag_reiners2016"]


def test_the_resolved_flag_is_derived_not_asserted(lineset):
    """It was hardcoded False before the published footnote was in hand. It must now say
    whether THIS row's code actually decodes — a status column that cannot change is the
    defect it was recording."""
    codes = set(B._read_jofre_codes())
    for _, r in lineset.iterrows():
        want = (not pd.isna(r.loggf_ref_code_gbs)) and int(r.loggf_ref_code_gbs) in codes
        assert bool(r.loggf_ref_gbs_resolved) == want, r.wavelength_air_A
    assert lineset.loggf_ref_gbs_resolved.nunique() == 2, "the flag never varies — vacuous"


def test_the_provenance_join_is_a_dual_key_with_a_measured_null():
    """Same discipline as the canonical_gf join: λ+EP, and the null is asserted, not hoped."""
    ges = B._read_ges_lines()
    ew = B._read_measurement_table(B.VIZIER / "ew.dat")
    sun = [r for r in ew if r["star"] == B.STAR]
    src = ges.rename(columns={"lambda": "wavelength_air_A", "Elow": "excitation_potential_eV"})
    src = src.assign(species=np.where(ges["Ion"] == 1, "Fe I", "Fe II"))
    assert B._null_control(sun, src) == [0] * len(B._NULL_SHIFTS_A)
    assert _resolve_at(sun, src, 0.0) == 159


def test_no_gf_VALUE_was_changed_by_the_provenance_pass(lineset):
    """The brief's hard constraint. The decode adds columns; it must move no number."""
    for col in ("log_gf_gbs", "log_gf_ours", "gf_synth_ges", "delta_gbs_minus_ours",
                "delta_gbs_minus_ges", "rew_class", "gbs_selected_sun"):
        assert col in lineset.columns
    # Heiter's value is recorded in its OWN column and never substituted into the GBS one.
    assert "heiter2021_log_gf" in lineset.columns
    d = lineset.gf_source_basis == "jofre2014-footnote"
    assert (lineset.log_gf_gbs[d] != lineset.heiter2021_log_gf[d]).all()


def test_the_decoder_holdings_are_intact():
    import rya1110_stage_heiter2021 as S
    assert S.verify() == 0


# ── the artifact must match the code that defines it ─────────────────────────
@pytest.mark.xfail(strict=True, reason=_GBS_O2GAMMA_REASON)
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

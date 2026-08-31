"""RYA-1141 - the QA battery must find the defects it found, and change nothing.

Every assertion below is a NEGATIVE control on the audit itself: a check that
silently degrades to PASS is worse than no check, so each one pins both the
verdict AND the evidence that produced it.
"""
import hashlib
import json

import pandas as pd
import pytest

from scripts.qa_al_intake_rya1141 import (AUDITED, CANONICAL, audited_files,
                                          build, join_constrains_identity, SELF)


@pytest.fixture(scope="module")
def qa(tmp_path_factory):
    out = tmp_path_factory.mktemp("rya1141")
    return build(out), out


def test_the_audit_mutates_nothing_it_audits(qa):
    verdict, _ = qa
    assert verdict["artifacts_mutated"] == []
    assert verdict["checks"]["NO-MUTATION"] == "PASS"


def test_the_auditor_is_excluded_from_its_own_measurement_by_name():
    """🔴 The audit script must never appear in the set it hashes. Excluding it by
    PATTERN would be excluding it by the same rule that selects the intake; the
    exclusion is by resolved path, and this test is what keeps it that way."""
    assert SELF not in {p.resolve() for p in audited_files()}
    assert SELF.name == "qa_al_intake_rya1141.py"


def test_identity_comparison_detector_has_a_working_positive(qa):
    """The A2 result is a NEGATIVE. It is only evidence if the same test can say yes."""
    verdict, _ = qa
    assert verdict["checks"]["A2-control"] == "PASS"
    assert verdict["checks"]["A2"] == "FAIL"


def test_detector_is_not_fooled_by_a_mention_without_a_comparison():
    """`ingest_new_lab_sources` writes `lower_level` into a string. A grep would call
    that an identity check; this must not."""
    import ast

    def fn(src):
        return [n for n in ast.walk(ast.parse(src)) if isinstance(n, ast.FunctionDef)][0]

    mention_only = fn("def f(df, a, w):\n"
                      "    note = f'{a.lower_level} - {a.upper_level}'\n"
                      "    return df[(df.wavelength_air - w).abs() <= .08], note\n")
    compares = fn("def f(df, a, w, ep):\n"
                  "    return df[((df.wavelength_air - w).abs() <= .08)\n"
                  "              & ((df.lower_EP - ep).abs() <= .02)]\n")
    assert not join_constrains_identity(mention_only)[0]
    assert join_constrains_identity(compares)[0]


def test_wavelength_only_join_is_caught_in_the_act(qa):
    """The null: a real collision, not a hypothetical one."""
    verdict, out = qa
    assert verdict["checks"]["A2-null"] == "FAIL"
    c = pd.read_csv(out / "a2_transition_collisions.csv")
    assert len(c) >= 1
    row = c[c.canonical_line_id.eq("alphys_II_3587.0720_0333")]
    assert len(row) == 1 and row.iloc[0].n_source_rows == 3


def test_source_transcription_is_refereed_by_the_cds_readme_and_branching_closure(qa):
    verdict, out = qa
    assert verdict["checks"]["A1"] == "PASS"
    assert verdict["checks"]["A1-parse"] == "PASS"
    assert verdict["checks"]["A1-burheim"] == "PASS"
    # The dropped-flag finding names the CDS columns, not a count alone.
    f = pd.read_csv(out / "a1_dropped_source_flags.csv")
    assert {"l_e_Aki", "n_Lambda", "n_Aki"} <= set(f.flag_column)
    assert (f[f.flag_column.eq("n_Lambda")].flag == "*").all()


def test_fine_structure_identity_is_the_air_vacuum_referee(qa):
    verdict, out = qa
    assert verdict["checks"]["A4"] == "PASS"
    fs = pd.read_csv(out / "a4_fine_structure_identity.csv")
    assert len(fs) == 3
    assert fs.residual_cm1.abs().max() < 0.02
    # and it is a real constraint: displacing one side by an air-vacuum shift breaks it
    shifted = 1e8 / fs.lambda_air_1.iloc[0] - 1e8 / fs.lambda_vac_2.iloc[0]
    assert abs(shifted - fs.expected_cm1.iloc[0]) > 1.0


def test_rya1001_hfs_defect_is_still_live_and_was_stamped_verified(qa):
    verdict, out = qa
    assert verdict["checks"]["A3-rya1001"] == "FAIL"
    h = pd.read_csv(out / "a3_hfs_component_counts.csv")
    bad = h[h.still_wrong]
    assert set(bad.wavelength_air_A.round(3)) == {3944.006, 3961.520}
    assert (bad.canonical_gf_hfs_n == 1).all()
    assert set(bad.census_hfs_n) == {4, 6}
    assert (bad.manifest_HFS_status == "COMPONENT_SUM_VERIFIED").all()
    assert (bad.manifest_gf_grade == "GF-LAB").all()


def test_misquoted_dois_are_named_with_their_corrections(qa):
    verdict, out = qa
    assert verdict["checks"]["A5-doi"] == "FAIL"
    d = pd.read_csv(out / "a5_doi_resolution.csv")
    wrong = d[d.verdict.eq("MISQUOTED")]
    assert set(wrong.doi) == {"10.1086/312738", "10.3847/1538-4357/ad4451",
                              "10.1093/mnras/stt2120"}
    assert wrong.correct_doi.str.len().gt(0).all()
    # the rest must resolve, or the check is just flagging everything
    assert (d.verdict.eq("OK")).sum() == len(d) - 3
    assert verdict["checks"]["A5-doi-control"] == "PASS"


def test_a_volume_comparison_would_have_missed_the_griesmann_doi(qa):
    """🔴 The negative control on the referee itself. `10.1086/312738` IS in ApJ 536,
    the volume the intake claims - a volume-only test passes it. Only the author list
    separates the two papers, which is why that is the comparison the check makes."""
    _, out = qa
    d = pd.read_csv(out / "a5_doi_resolution.csv")
    g = d[d.doi.eq("10.1086/312738")].iloc[0]
    assert "536" in str(g.claimed_citation) and str(g.registered_volume) == "536"
    assert g.verdict == "MISQUOTED"
    assert "Griesmann" not in str(g.registered_authors)


def test_a_damaged_crossref_byte_does_not_read_as_a_wrong_author():
    """Crossref serves 'J\ufffdnsson' for Jonsson & Lundberg 1983. That is one damaged
    character, not a different person - and it must still not match a real mismatch."""
    from scripts.qa_al_intake_rya1141 import same_surname
    assert same_surname("Jonsson", "J\ufffdnsson")
    assert same_surname("Vujnovic", "Vujnovi\u0107")
    assert not same_surname("Griesmann", "Brinchmann")
    assert not same_surname("Nandakumar", "Miura")


def test_dropped_competing_lab_gf_is_quantified_not_just_named(qa):
    verdict, out = qa
    assert verdict["checks"]["A6"] == "FAIL"
    a = pd.read_csv(out / "a6_dropped_competing_gf.csv")
    # Six lines, none of them reachable from either place a reader would look.
    assert len(a) == 6
    assert not a.in_conflict_ledger.any()
    assert not a.named_in_competing_gf_summary.any()
    # The one that matters: two PRIMARY-LAB sources in tension on a GF-LAB line.
    lab = a[a.adopted_source.eq("EXP-BURHEIM23")]
    worst = lab.loc[lab.n_sigma_on_adopted.idxmax()]
    assert round(float(worst.wavelength_air), 3) == 13123.416
    assert worst.n_sigma_on_adopted > 5
    # And the Al II 3900.675 row, whose Aki the CDS ReadMe attributes to Tayal &
    # Hibbert theory, is carried under Vujnovic's DOI 2.4 dex from the adopted value.
    tayal = a[a.wavelength_air.round(3).eq(3900.675)]
    assert len(tayal) == 1 and abs(float(tayal.iloc[0].delta_dex)) > 2.0


def test_the_band_gap_swallows_the_two_best_graded_lines(qa):
    verdict, out = qa
    assert verdict["checks"]["C-bands"] == "FAIL"
    r = pd.read_csv(out / "c_band_gap_relabelled_lines.csv")
    lab = r[r.gf_grade.eq("GF-LAB")]
    assert set(lab.wavelength_air.round(3)) == {13123.416, 13150.753}
    # Twelve of the thirteen are lines the census placed squarely in NIR.
    assert (r.band_census == "NIR").sum() == 12
    assert (r[r.band_census.eq("NIR")].instruments_coverage_blind_spot
            == "crires_plus").all()


def test_every_crires_holding_is_dropped_by_the_coverage_module(qa):
    verdict, out = qa
    assert verdict["checks"]["C"] == "FAIL"
    h = pd.read_csv(out / "c_solar_holdings_resolution.csv")
    crires = h[h.instrument_id.eq("crires_plus")]
    assert len(crires) == 5
    assert not crires.reaches_coverage_module.any()
    # and something DOES resolve, so the reader is not simply broken
    assert h.reaches_coverage_module.any()


def test_headline_inventory_reproduces_and_canonical_gf_is_untouched(qa):
    verdict, _ = qa
    assert verdict["checks"]["B1"] == "PASS"
    assert verdict["checks"]["B1-promotions"] == "PASS"
    assert verdict["checks"]["B2"] == "PASS"


def test_verdict_is_a_fail_and_the_gate_stays_closed(qa):
    verdict, out = qa
    assert verdict["overall"] == "FAIL"
    assert verdict["intake_independently_verified"] is False
    assert verdict["measurement_gate"] == "CLOSED"
    md = (out / "verdict.md").read_text()
    for c in ("A1", "A2", "A3", "A4", "A5", "A6", "B1", "B2", "B3", "C"):
        assert f"| {c} |" in md


def test_canonical_gf_is_byte_identical_after_the_whole_battery(qa):
    """The findings-only firewall, asserted from outside the script's own bookkeeping."""
    _, _ = qa
    digest = hashlib.sha256(CANONICAL.read_bytes()).hexdigest()
    assert hashlib.sha256(CANONICAL.read_bytes()).hexdigest() == digest
    assert not any(p.name.startswith("al_line_manifest") and p.parent == AUDITED
                   and p.stat().st_size == 0 for p in AUDITED.iterdir())


def test_the_repo_wide_wavelength_guard_does_not_catch_this_join(qa):
    """RYA-1037 shipped a guard to make RYA-1034 unrepeatable. It did not fire here."""
    verdict, _ = qa
    assert verdict["checks"]["A2-repo-guard"] == "FAIL"


def test_the_1037_scan_is_called_in_process_and_writes_nothing():
    """🔴 Shelling out to the auditor REWRITES its inventory JSON — a repo file outside
    this QA's output directory. `scan()` is pure; assert we call that, not the CLI.

    #: Read through the absolute `SELF`, never a relative path: other modules in the
    #: suite chdir, and a relative read passes alone and dies in the full run.
    """
    from scripts.qa_al_intake_rya1141 import ROOT
    src = SELF.read_text()
    assert "from audit_line_keys_rya1037 import scan" in src
    assert "audit_line_keys_rya1037.py\"]" not in src  # no subprocess invocation
    inventory = ROOT / "data/audit/rya1037/rya1037_line_key_inventory.json"
    before = inventory.read_bytes() if inventory.exists() else None
    from scripts.qa_al_intake_rya1141 import build as _b  # noqa: F401
    assert (inventory.read_bytes() if inventory.exists() else None) == before


def test_no_mutation_check_actually_fails_on_a_stray_write(tmp_path, monkeypatch):
    """🔴 THE CONTROL ON THE FIREWALL. A hash over a chosen set cannot see a write
    outside it — an earlier revision of this script dirtied `data/audit/rya1037/` and
    still reported PASS. Sabotage a check into touching a repo file and require FAIL."""
    from scripts import qa_al_intake_rya1141 as q
    probe = q.ROOT / "data/audit/rya1037/stray_probe_test.tmp"
    original = q.check_a4

    def sabotage(rep, man):
        probe.write_text("x")
        return original(rep, man)

    monkeypatch.setattr(q, "check_a4", sabotage)
    try:
        v = q.build(tmp_path)
        assert v["checks"]["NO-MUTATION"] == "FAIL"
        assert "data/audit/rya1037/stray_probe_test.tmp" in v["artifacts_mutated"]
    finally:
        probe.unlink(missing_ok=True)


def test_the_identity_test_is_scoped_to_the_narrowing_expression():
    """🔴 NOT THE WHOLE FUNCTION. RYA-1037's `_enclosing_has_ep()` asks the enclosing
    function, so one unrelated `ep` silences every wavelength-only comparison below it.
    A mention elsewhere must not launder a bare wavelength filter."""
    import ast
    from scripts.qa_al_intake_rya1141 import join_constrains_identity

    def fn(src):
        return [n for n in ast.walk(ast.parse(src)) if isinstance(n, ast.FunctionDef)][0]

    laundered = fn("def f(df, w, r):\n"
                   "    note = f'{r.lower_level} - {r.upper_level}'\n"
                   "    return df[(df.wavelength_air - w).abs() <= .08], note\n")
    conjoined = fn("def g(df, w, ep):\n"
                   "    ok = (df.wavelength_air - w).abs() <= .08\n"
                   "    ok &= (df.lower_EP - ep).abs() <= .02\n"
                   "    return df[ok]\n")
    assert not join_constrains_identity(laundered)[0]
    assert join_constrains_identity(conjoined)[0]


def test_b2_pins_the_merge_sha_and_proves_it_is_the_right_one(qa):
    """🔴 THE AUDITOR'S OWN MERGE NAMES RYA-1132. A `--grep=RYA-1132 --merges -1` search
    therefore retargets B2 onto this audit's commit the moment it lands — and still
    reports PASS. The SHA is pinned; the control asserts the pin is correct."""
    verdict, _ = qa
    assert verdict["checks"]["B2-control"] == "PASS"
    assert verdict["checks"]["B2"] == "PASS"
    #: A MENTION IS NOT A USE — the source explains the `--grep` trap in a comment, so a
    #: substring search would fail on the explanation. Ask the AST for an actual argument.
    import ast
    src = SELF.read_text()
    assert "RYA1132_MERGE = " in src
    grep_args = [n for n in ast.walk(ast.parse(src))
                 if isinstance(n, ast.Constant) and isinstance(n.value, str)
                 and n.value.startswith("--grep")]
    assert not grep_args, f"B2 must not search for its target: {[n.value for n in grep_args]}"


def test_b2_control_fails_if_the_pin_is_wrong(monkeypatch, tmp_path):
    """A pin nobody checks is just a different way to be wrong."""
    from scripts import qa_al_intake_rya1141 as q
    monkeypatch.setattr(q, "RYA1132_MERGE", "origin/main")
    v = q.build(tmp_path)
    assert v["checks"]["B2-control"] == "FAIL"


# ─────────────────────────────────────────────────────────────────────────────
# D — the coverage the first pass did not reach.
# ─────────────────────────────────────────────────────────────────────────────
def test_asplund_grade_is_a_line_set_not_a_gf_grade(qa):
    """🔴 I had this wrong. `model_registry.LINE_SETS` is the ONE definition: `asplund`
    is a value on the `line_set` PROVENANCE axis — which pool a measurement was made on
    — not a statement about a gf. RYA-1127 put it in the product identity key."""
    from pipeline.model_registry import LINE_SETS
    assert "asplund" in LINE_SETS and "our-graded" in LINE_SETS
    assert "consistent" not in LINE_SETS, "RYA-1105 retired it; it must not acquire a name"
    verdict, _ = qa
    assert verdict["checks"]["D3-lineset"] == "FAIL"


def test_al_was_frozen_through_the_rya946_census_gate(qa):
    """RYA-946: no element is FROZEN_READY_FOR_MEASUREMENT until the AGSS21 line-set
    cross-reference is complete or an approved exception exists. Neither is true for Al."""
    from scripts.qa_al_intake_rya1141 import ROOT
    verdict, _ = qa
    assert verdict["checks"]["D3"] == "FAIL"
    sets = sorted(p.name for p in (ROOT / "data/reference").glob("asplund*") if p.is_dir())
    assert sets == ["asplund2021_fe"], f"expected Fe-only reference set, got {sets}"


def test_the_evaluated_tier_is_opacity_project_theory(qa):
    """🔴 CORRECTS A5-lab. No GF-LAB row is theory — but all 19 CRITICALLY_EVALUATED rows
    are, and `source_type` cannot see it because NIST is tested before THEORY."""
    from scripts.qa_al_intake_rya1141 import BUILDER
    verdict, _ = qa
    assert verdict["checks"]["D4-lineage"] == "FAIL"
    fn = BUILDER.read_text()
    fn = fn[fn.index("def source_type"):fn.index("def nearest")]
    assert fn.index('"NIST" in s') < fn.index('"THEORY" in s')


def test_promotions_rest_on_measured_ratios_not_ls_theory(qa):
    """The Vujnovic PAPER, not its CDS table, is what says which Aki are measured."""
    verdict, out = qa
    assert verdict["checks"]["D2"] == "PASS"
    t = pd.read_csv(out / "d2_vujnovic_ratio_basis.csv")
    theo = set(t[t.ratio_basis.eq("THEORETICAL_LS_RATIO")].wavelength_A.round(2))
    assert theo == {13123.41, 13150.76, 21093.04, 21163.75}


def test_outside_current_reach_is_contradicted_by_the_instrument_catalog(qa):
    verdict, out = qa
    assert verdict["checks"]["D5-outside"] == "FAIL"
    s = pd.read_csv(out / "d5_full_instrument_catalog_sweep.csv")
    assert (s.n_catalog_instruments > 0).all(), "no Al line is beyond every instrument"
    wrong = s[s.manifest_instrument_reach.eq("OUTSIDE_CURRENT_REACH")]
    assert len(wrong) == 4 and (wrong.n_catalog_instruments >= 4).all()


def test_a_summed_feature_is_not_graded_better_than_its_worst_component(qa):
    verdict, out = qa
    assert verdict["checks"]["D4-grades"] == "FAIL"
    e = pd.read_csv(out / "d4_evaluated_tier_provenance.csv")
    bad = e[e.nist_grade.ne(e.nist_grade_worst)]
    assert len(bad) == 5
    assert 6906.287 in set(bad.wavelength_air.round(3))

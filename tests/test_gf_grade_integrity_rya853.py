"""
tests/test_gf_grade_integrity_rya853.py — RYA-853
=================================================
The "Type B uncertainty anchor" disagrees with the source it cites on 70% of its verifiable
grades. RYA-850 keys `graded_gf_term` on that metadata, so a stored `A` on a line NIST grades
`D` publishes 0.013 dex where the source says 0.176.

The most valuable test here needs NO NETWORK: the two extracts are meant to describe the same
lines, so where they disagree with EACH OTHER at least one is wrong. That check would have
caught RYA-592's half-applied fix — `nist_reference` corrected, `nist_crosscheck` left stale
— a month before this audit did.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

NIST_REF = ROOT / "data" / "linelists" / "nist_reference.csv"
NIST_XC = ROOT / "data" / "linelists" / "nist_crosscheck.csv"
SUMMARY = ROOT / "data" / "results" / "rya853" / "rya853_summary.json"

VALUE_TOL_DEX = 0.02
WAVE_TOL_A = 0.05


@pytest.fixture(scope="module")
def summary():
    if not SUMMARY.exists():
        pytest.skip("RYA-853 summary absent (Sirius artifact)")
    return json.loads(SUMMARY.read_text())


def _pairs():
    """Lines carried by BOTH extracts, matched on species + wavelength."""
    ref = pd.read_csv(NIST_REF, comment="#")
    xc = pd.read_csv(NIST_XC, comment="#")
    out = []
    for _, a in xc.iterrows():
        m = ref[(ref.element == a.element) & (ref.ion == a.ion)
                & (np.abs(ref.wavelength_air_A - a.wavelength_air_A) <= WAVE_TOL_A)]
        if len(m):
            out.append((a, m.iloc[0]))
    return out


# ── the offline guard — no network, runs anywhere, catches half-applied fixes ──

@pytest.mark.xfail(strict=True, reason=
    "KNOWN DEFECT, RYA-853. The two extracts disagree and the files have NOT been\n     corrected here: these rows feed published anchors (O I 6300, Li I 6707,\n     Ba II 5853, Ni I 6300.336), so applying 31 grade and 14 value corrections is a\n     reviewable change with a real blast radius, not an audit side effect.\n     strict=True on purpose — when the files ARE corrected this flips to an\n     unexpected PASS and fails, forcing the marker to be removed deliberately.\n     The corrections are tabulated in data/results/rya853/rya853_corrections.csv.")
def test_the_two_extracts_agree_on_grade():
    """🔴 CURRENTLY FAILING BY DESIGN IS NOT THE INTENT — this documents the defect and
    must be made to pass by CORRECTING the files, not by loosening the check.

    Mg I 5711.090 is the smoking gun: RYA-592 corrected it in nist_reference and left
    nist_crosscheck stale, so the same line carries two different grades depending on which
    file you read.
    """
    bad = [(a.element, a.ion, float(a.wavelength_air_A),
            str(b.nist_grade).strip(), str(a.nist_grade).strip())
           for a, b in _pairs()
           if str(a.nist_grade).strip() != str(b.nist_grade).strip()]
    assert not bad, (
        "the two NIST extracts disagree on the GRADE for lines both carry — at least one "
        f"is wrong, and no network is needed to see it: {bad}")


@pytest.mark.xfail(strict=True, reason=
    "KNOWN DEFECT, RYA-853. The two extracts disagree and the files have NOT been\n     corrected here: these rows feed published anchors (O I 6300, Li I 6707,\n     Ba II 5853, Ni I 6300.336), so applying 31 grade and 14 value corrections is a\n     reviewable change with a real blast radius, not an audit side effect.\n     strict=True on purpose — when the files ARE corrected this flips to an\n     unexpected PASS and fails, forcing the marker to be removed deliberately.\n     The corrections are tabulated in data/results/rya853/rya853_corrections.csv.")
def test_the_two_extracts_agree_on_log_gf():
    """S I 6052.670 differs by 0.92 dex between our own two files; Ni I 6300.336 by 0.20,
    and that Ni line is the blend under [O I] 6300 that gates the solar oxygen abundance
    (RYA-365)."""
    bad = [(a.element, a.ion, float(a.wavelength_air_A),
            float(b.log_gf), float(a.log_gf))
           for a, b in _pairs()
           if abs(float(a.log_gf) - float(b.log_gf)) > VALUE_TOL_DEX]
    assert not bad, (
        f"the two NIST extracts disagree on log gf for lines both carry: {bad}")


# ── the audited result, pinned ────────────────────────────────────────────────

def test_the_stored_grades_disagree_with_nist_at_scale(summary):
    """The finding. If this ever passes trivially because the mismatch count dropped to
    zero, the files were corrected — update this test in the same commit."""
    assert summary["n_grade_mismatch"] > 0, (
        "no grade mismatches remain — if the extracts were corrected, replace this test "
        "with one asserting they STAY correct")
    frac = summary["n_grade_mismatch"] / max(summary["n_uniquely_matched"], 1)
    assert frac > 0.5, "the mismatch rate collapsed — re-derive before trusting it"


def test_only_uniquely_matched_rows_are_judged(summary):
    """A wavelength+EP window is not a unique line identifier: several transitions share a
    wavelength and a lower level. Taking the first match manufactured 12-dex 'defects'
    (Mg I 5183.604 stored +0.180 against a NIST row at -11.908). Those rows must be
    reported AMBIGUOUS, not judged."""
    assert summary["n_ambiguous"] > 0, (
        "no ambiguous rows — the guard against multi-candidate matches may have been lost")
    assert (summary["n_uniquely_matched"] + summary["n_ambiguous"]
            + summary["n_unmatched"]) >= summary["n_stored_rows"] - 1
    for v in summary["value_mismatches"]:
        assert abs(v["delta"]) < 1.0, (
            f"a >1 dex 'mismatch' survived into the judged set ({v}) — that is a "
            f"matching artifact, not a transcription defect")


def test_the_understatement_is_not_one_directional(summary):
    """Li I is stored WIDER than NIST (A/A+ against AAA) while most rows are stored
    tighter. A one-directional error would suggest simple optimism; both directions means
    transcription drift, which is what the file's own header describes."""
    ratios = [g["understated_by"] for g in summary["grade_mismatches"]
              if g.get("understated_by")]
    assert any(r > 1 for r in ratios) and any(r < 1 for r in ratios)


# ── scope 4: the cross-match guards ───────────────────────────────────────────

def test_the_live_nist_sites_carry_both_guards(summary):
    """astroquery.nist defaults to VACUUM (+1.71 A at 6150), and EP must be matched on both
    sides. Every site that queries NIST live must have both."""
    for s in summary["cross_match_guards"]["sites"]:
        if s.get("air_guard") is None:
            continue                      # local-file matcher, no query to guard
        assert s["air_guard"], f"{s['file']} lost its wavelength_type='vac+air'"
        assert s["ep_guard"], f"{s['file']} lost its EP guard"


def test_the_wavelength_only_matcher_is_recorded(summary):
    """rya347 matches on wavelength alone within 0.1 A. It did NOT cause the Fe II defect
    — the EP there is correct and the right line was matched — but it is a latent risk and
    is carried so it does not get rediscovered as a cause."""
    loose = [s for s in summary["cross_match_guards"]["sites"]
             if s.get("wavelength_only_window_A")]
    assert loose, "the wavelength-only matcher vanished from the audit"
    assert "corrected_diagnosis" in summary
    assert "NOT" in summary["corrected_diagnosis"]


def test_the_dh19_referee_is_marked_owed_not_answered(summary):
    """Scope 3 decides whether the ionization balance is an independent check or a circular
    one, and it is NOT done. It must not read as settled."""
    txt = json.dumps(summary).lower()
    assert "n_grade_mismatch" in summary
    # the audit must not claim the scale offset was characterised
    assert "scale_offset_verdict" not in summary


# ── scope 3: the DH19 referee ─────────────────────────────────────────────────
#
# 🔴 WHAT WENT WRONG HERE, AND WHAT THESE TESTS NOW PIN
# [RYA-945] ingested Den Hartog 2019 into canonical_gf. The referee read "ours" from
# canonical_gf, so by 2026-08-27 it was comparing DH19 against DH19: ours - DH = +0.000,
# CI [+0.000, +0.000], sd 0.000 — and printed "LEGITIMATE ... solar-fitting REFUTED".
# The verdict flipped from INCONCLUSIVE without anyone touching RYA-853.
#
# The bootstrap CI was added to stop the median being overread, and it did not catch this:
# a zero-width CI passes every "is it well determined?" test there is. So the tests below
# check the two properties a self-match cannot have — an independent source, and a spread.

REFEREE = ROOT / "data" / "results" / "rya853" / "rya853_dh19_referee.json"
REFEREE_LIVE = ROOT / "data" / "results" / "rya853" / "rya853_dh19_referee_live.json"
FROZEN_SNAPSHOT = (ROOT / "data" / "reference" / "fe_gf_lab"
                   / "fe2_pre945_scale_snapshot.csv")


@pytest.fixture(scope="module")
def referee():
    if not REFEREE.exists():
        pytest.skip("DH19 referee artifact absent")
    return json.loads(REFEREE.read_text())


def test_the_referee_is_pure_lab(referee):
    """A referee that touched the solar spectrum could not settle a question about solar
    fitting. DH19 is branching fractions x LIF lifetimes — no solar normalisation."""
    assert "pure lab" in referee["referee"].lower()
    assert "BF" in referee["referee"] or "LIF" in referee["referee"]


def test_the_referee_is_not_on_both_sides_of_its_own_comparison(referee):
    """🔴 THE REGRESSION TEST. 'ours' must be the pre-945 snapshot, and no scored line may
    cite the referee. Reading canonical_gf here is what produced the false REFUTED."""
    assert referee["ours_source"] == "frozen", (
        "the referee is reading live canonical_gf again — post-RYA-945 that file CONTAINS "
        "Den Hartog 2019, so the comparison scores the referee against itself")
    assert "pre945" in referee["ours_source_file"]
    assert referee["n_excluded_self_referential"] == 0, (
        "a scored line cites the referee — the frozen snapshot has been contaminated")


def test_the_estimator_did_not_collapse(referee):
    """A spread of exactly zero is not agreement, it is the same number on both sides. This
    is the numeric signature the false verdict had: sd 0.000, CI width 0.000."""
    st = referee["ours_minus_dh"]
    assert st["std"] > 0.0, "zero spread across the pool — this is a self-match, not a test"
    assert st["ci_width"] > 1e-6
    assert st["n"] >= 5, "too few lines to carry a verdict"


def test_a_verdict_is_only_claimed_when_the_ci_separates_the_readings(referee):
    """The two readings predict ~0.0 (pure lab) and ~+0.13 (solar-fitted). Claiming either
    requires a CI that EXCLUDES the other — the original bug was thresholding |median| on
    its own, which a CI covering both still passes."""
    lo, hi = referee["ours_minus_dh"]["ci95"]
    v = referee["verdict"]
    if v.startswith("LEGITIMATE"):
        assert not (lo <= 0.13 <= hi), (
            "LEGITIMATE claimed while the CI still covers the solar-fitted prediction")
    elif v.startswith("HYPOTHESIS LIVES"):
        assert not (lo <= 0.0 <= hi), (
            "HYPOTHESIS LIVES claimed while the CI still covers the pure-lab prediction")
    elif v.startswith("INCONCLUSIVE"):
        assert lo <= 0.0 <= hi and lo <= 0.13 <= hi
    else:
        assert v.startswith("UNUSABLE")


def test_the_verdict_rests_on_the_near_uv_lines_not_the_blue_overlap(referee):
    """WHERE THE ANSWER COMES FROM. The blue overlap alone is the same ten lines that were
    INCONCLUSIVE in August and still are; the twelve near-UV lines RYA-945's full table
    added are what tighten it. Stated so nobody reads this as the blue arm having changed
    its mind."""
    per_band = referee["ours_minus_dh_per_band"]
    blue = per_band["blue"]
    lo, hi = blue["ci95"]
    assert lo <= 0.0 <= hi and lo <= 0.13 <= hi, (
        "the blue arm now separates the readings on its own — re-read the verdict, it is "
        "no longer carried by the near-UV lines")
    assert per_band["near-UV"]["n"] >= 10


def test_the_referee_is_reproducible_from_the_repo(referee):
    """The previous run transcribed ten values out of the DH19 PDF and said so as a real
    limitation. RYA-945 vendored the machine-readable Table 6, so the referee is now read
    from a committed file — 131 lines, not 10."""
    assert referee["referee_source"].endswith("fe2_lab_loggf_dh19.csv")
    assert referee["referee_n_lines"] == 131
    assert (ROOT / referee["referee_source"]).exists()


def test_both_band_arms_are_measured_not_carried_as_a_literal(referee):
    """RYA-852's +0.106 was taken before 852/877/945 touched the pool. Quoting it beside a
    freshly measured blue arm compares two pools on two dates and calls the difference
    band-dependence. Both arms must carry their own n."""
    b = referee["band_dependence"]
    assert b is not None
    assert b["blue_n"] >= 5 and b["red_n"] >= 5
    assert b["sign_flips"] is True
    assert abs(b["swing_dex"]) > 0.1


def test_the_arbiter_lines_are_not_claimed_to_be_refereed(referee):
    """DH19 stops at 4584 A; the three arbiter lines are redward of it."""
    assert "arbiter" in referee["caveat_arbiter_lines"].lower()
    assert referee["scored_span_A"][1] < 6000.0


def test_the_independence_caveat_travels_with_the_verdict(referee):
    """GUARD 1 checks the reference STRING. The near-UV lines carrying the verdict are
    plain 'VALD3', and VALD3 is an aggregator — that label does not prove independence
    from Den Hartog, it only proves we did not adopt him directly."""
    c = referee["caveat_independence_is_string_deep_only"].lower()
    assert "vald3" in c and "aggregator" in c


def test_the_balance_caveat_records_that_the_measurement_leg_is_pre_continuum_fix(referee):
    """🔴 A CLEARED gf SCALE DOES NOT CLEAR THE BALANCE. Fe I 7.586 / Fe II 7.568 were
    measured 2026-08-16..18; every continuum fix (RYA-911/913, 1000/1006, 1026, 1030)
    landed 2026-08-19 or later, and the Fe II VIS cell has never been re-run."""
    c = referee["caveat_balance_is_pre_continuum_fix"]
    for t in ("1026", "1030", "911", "re-run"):
        assert t in c, f"the balance caveat lost {t!r}"


def test_the_frozen_snapshot_does_not_contain_the_referee():
    """The snapshot is the whole defence. If a future ingest is ever frozen into it, the
    referee silently grades its own homework again."""
    if not FROZEN_SNAPSHOT.exists():
        pytest.skip("frozen pre-945 snapshot absent")
    d = pd.read_csv(FROZEN_SNAPSHOT, low_memory=False)
    ref = d.loggf_reference.astype(str).str.lower()
    hits = ref.str.contains("denhartog2019|den hartog 2019|dh19|2019apjs..243", regex=True,
                            na=False)
    assert int(hits.sum()) == 0, (
        f"{int(hits.sum())} snapshot rows cite the referee — this is no longer a "
        f"pre-adoption scale")


def test_the_guard_refuses_the_degenerate_run_if_it_was_recorded():
    """`--source live` reproduces the broken comparison on purpose. When that artifact is
    present it must carry a refusal, never a verdict."""
    if not REFEREE_LIVE.exists():
        pytest.skip("live-source control run not recorded")
    live = json.loads(REFEREE_LIVE.read_text())
    assert live["verdict"].startswith("UNUSABLE")
    assert live["guard_refusal"]
    assert live["n_excluded_self_referential"] > 0
    assert live["n_scored"] == 0


def test_the_degeneracy_threshold_actually_catches_a_constant_series():
    """Unit-level, no network: the exact input that fooled the old code. stat() on a
    constant series must report zero spread, so the guard has something to fire on."""
    import importlib
    m = importlib.import_module("rya853_dh19_scale_referee")
    st = m.stat(pd.Series([0.0] * 10))
    assert st["std"] == 0.0
    assert st["ci_width"] < m.DEGENERATE_CI_WIDTH


# ── scope 2: the Fe I lab pool (the pool RYA-850 promotes) ────────────────────

LABPOOL = ROOT / "data" / "results" / "rya853" / "rya853_fe1_labpool_referee.json"


@pytest.fixture(scope="module")
def labpool():
    if not LABPOOL.exists():
        pytest.skip("Fe I lab-pool referee artifact absent")
    return json.loads(LABPOOL.read_text())


def test_the_fe1_lab_pool_is_materially_clean(labpool):
    """THE RESULT THAT CLEARS RYA-850's POOL. 207 of 217 refereed lines reproduce their
    source paper's log gf — the opposite of the NIST extracts' 70% failure. If this rate
    ever drops, 850's graded bars are back in question."""
    assert labpool["coverage_frac"] > 0.95, (
        "coverage collapsed — the CDS machine-readable tables may have stopped parsing")
    clean = labpool["n_refereed"] - labpool["n_loggf_mismatch"]
    frac = clean / max(labpool["n_refereed"], 1)
    assert frac > 0.95, (
        f"only {frac:.1%} of refereed lab lines reproduce their source — the pool "
        f"RYA-850 promotes is no longer clean")


def test_the_one_real_outlier_is_named(labpool):
    """Belmonte 3935.307: ours -2.199 vs the paper's -1.820, and its sigma is wrong too
    (0.070 vs 0.180). Both axes disagree on one line — a bad row, not rounding."""
    big = [m for m in labpool["loggf_mismatches"] if abs(m["delta"]) > 0.15]
    assert len(big) == 1
    assert abs(big[0]["wavelength_air_A"] - 3935.307) < 0.01


def test_the_unrefereed_remainder_is_not_called_clean(labpool):
    """Coverage was 47% while only the PDF excerpts were available; vendoring the CDS
    machine-readable tables took it to 99.6%. What remains unrefereed is 2 Belmonte lines,
    and absent from the check is still not the same as verified (RYA-833)."""
    unrefereed = labpool["n_pool_rows"] - labpool["n_refereed"]
    assert unrefereed >= 0
    assert labpool["coverage_frac"] > 0.95, "coverage regressed — did the CDS tables move?"
    assert "UNVERIFIED" in labpool["caveat"] or "unverified" in labpool["caveat"].lower()


def test_the_encoding_traps_are_recorded(labpool):
    """Each produced a confident wrong answer: the unicode minus turned 95 of 99 lines into
    'mismatches' of exactly twice their value, and PUA padding made Belmonte look like it
    covered none of our lines."""
    traps = " ".join(labpool["traps"]).lower()
    assert "nanometre" in traps or "nanometres" in traps
    assert "pdftotext" in traps
    src = (ROOT / "scripts" / "rya853_fe1_labpool_referee.py").read_text()
    assert "u2212" in src.lower(), "the unicode-minus guard is gone"
    assert "uf8ff" in src.lower(), "the private-use-area guard is gone"


def test_ruffoni_and_denhartog_are_perfect(labpool):
    """345 lines across the two Wisconsin sources, refereed against their CDS machine-
    readable tables, with ZERO defects on value or cited sigma. Every surviving mismatch
    is Belmonte, which is the one source refereed from a typeset PDF."""
    bad_sources = {m["source"] for m in labpool["loggf_mismatches"]}
    assert bad_sources <= {"Belmonte2017"}, (
        f"a Wisconsin source now has mismatches ({bad_sources}) — that would be a real "
        f"defect in the pool RYA-850 promotes, not a PDF-extraction artifact")


def test_belmonte_coverage_did_not_silently_collapse(labpool):
    """Belmonte covers 98.3% of its own pool lines. A drop to zero means the PUA/± parse
    broke again, not that the paper stopped covering them."""
    assert labpool["papers"]["Belmonte2017"]["rows_parsed"] > 100

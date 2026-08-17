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

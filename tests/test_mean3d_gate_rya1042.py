"""RYA-1042: the RYA-534 gate, re-expressed for a ⟨3D⟩ deck.

The old gate could not run on one at all — it drives `interpol_modeles_nlte`, and no vendor
binary reads a ⟨3D⟩ deck (RYA-821). The LOGIC is unchanged; only the read path is.

What these tests protect is the gate's ability to FAIL. A validation gate that cannot
return NOT_VALIDATED is decoration, and several of the ways this one could quietly stop
failing are subtle: a sign flip invisible to a `|difference|` test, a thin match rate
reading as a solid result, a vturb bracket whose two ends disagree, and an anchor
accidentally sourced from the machinery under test.
"""
import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
GATE = ROOT / "scripts" / "rya1042_mean3d_gate.py"
ANCHOR = ROOT / "data" / "nlte_grids" / "amarsi2016_fe" / "amarsi2016_mean3d_solar_anchor.csv"
SRC = GATE.read_text()


def _lines_csv(path: Path, rows):
    """A minimal per-line product CSV in the schema the band route emits."""
    cols = ["element", "ion", "wavelength_air_A", "abundance", "in_aggregate", "ep_eV"]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for wave, ab, ep in rows:
            w.writerow(dict(element="Fe", ion="I", wavelength_air_A=wave,
                            abundance=ab, in_aggregate=True, ep_eV=ep))


def _run(tmp, nlte_rows, lte_rows):
    n, l = tmp / "n.csv", tmp / "l.csv"
    _lines_csv(n, nlte_rows); _lines_csv(l, lte_rows)
    out = tmp / "v.json"
    r = subprocess.run([sys.executable, str(GATE), "--nlte", str(n), "--lte", str(l),
                        "--out", str(out)], capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr
    return json.loads(out.read_text()), r.stdout


@pytest.fixture(scope="module")
def anchor_waves():
    """Real anchor wavelengths, so the match path is exercised against the real file."""
    with ANCHOR.open() as fh:
        return [float(r["lambda_air_nm"]) * 10.0 for r in csv.DictReader(fh)
                if r["species"] == "Fe1" and r["vturb_kms"] == "1.50"
                and r["clean"] == "yes"][:12]


# ── the gate must be able to FAIL ───────────────────────────────────────────

def test_a_deck_that_disagrees_with_the_anchor_FAILS(tmp_path, anchor_waves):
    """The whole point. A differential far from the anchor must not pass."""
    rows_n = [(w, 7.50 + 0.5, 3.0) for w in anchor_waves]      # +0.5 dex differential
    rows_l = [(w, 7.50, 3.0) for w in anchor_waves]
    doc, _ = _run(tmp_path, rows_n, rows_l)
    assert doc["verdict"] == "NOT_VALIDATED"
    assert any(v.get("verdict") == "FAIL" for v in doc["leg1_anchor_agreement"].values())


def test_a_sign_flip_is_reported_even_when_the_magnitude_passes(tmp_path, anchor_waves):
    """🔴 A `|difference|` TEST CANNOT SEE A SIGN FLIP. The anchor's solar Fe I differential
    is negative (≈ −0.03); a deck returning the same size with the opposite sign can sit
    inside a 0.05 tolerance while describing the opposite physics. Stated separately."""
    rows_n = [(w, 7.50 + 0.030, 3.0) for w in anchor_waves]    # +0.030 vs anchor -0.03
    rows_l = [(w, 7.50, 3.0) for w in anchor_waves]
    doc, out = _run(tmp_path, rows_n, rows_l)
    assert "SIGN DISAGREEMENT" in out
    assert any(v.get("sign_disagreement") for v in doc["leg1_anchor_agreement"].values())


def test_a_thin_match_rate_is_a_headline_not_a_footnote(tmp_path, anchor_waves):
    """A gate run over a handful of matched lines is weak however clean its arithmetic
    looks. Burying that in a notes list is how a thin result reads as a solid one."""
    rows_n = [(anchor_waves[0], 7.47, 3.0), (9999.999, 7.47, 3.0)]
    rows_l = [(anchor_waves[0], 7.50, 3.0), (9999.999, 7.50, 3.0)]
    doc, out = _run(tmp_path, rows_n, rows_l)
    assert "THIN" in out
    for v in doc["leg1_anchor_agreement"].values():
        if v.get("n"):
            assert v["match_rate"] < 1.0 and "n_ours" in v


def test_an_unmatchable_product_cannot_silently_pass(tmp_path):
    """No overlap with the anchor is NOT agreement. It must not read as a pass."""
    rows_n = [(9000.0 + i, 7.47, 3.0) for i in range(5)]
    rows_l = [(9000.0 + i, 7.50, 3.0) for i in range(5)]
    doc, out = _run(tmp_path, rows_n, rows_l)
    assert doc["verdict"] == "NOT_VALIDATED"
    assert "NO LINES MATCHED" in out


def test_no_differential_exists_when_the_legs_share_no_line(tmp_path, anchor_waves):
    """The NLTE effect is a per-line DIFFERENCE. A line measured in one leg and excluded
    in the other cannot contribute to it, and zero shared lines is not a zero effect."""
    n, l = tmp_path / "n.csv", tmp_path / "l.csv"
    _lines_csv(n, [(anchor_waves[0], 7.47, 3.0)])
    _lines_csv(l, [(anchor_waves[1], 7.50, 3.0)])
    r = subprocess.run([sys.executable, str(GATE), "--nlte", str(n), "--lte", str(l)],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode != 0
    assert "no line is in-aggregate in BOTH" in r.stdout + r.stderr


# ── the anchor stays external, and its caveats stay visible ─────────────────

def test_the_anchor_is_the_scale_matched_file_and_never_our_own_output():
    """🔴 RYA-161/1035. The reference must never come from the machinery under test —
    that is the closed loop RYA-1035 found in this deck's own abundance record."""
    assert "amarsi2016_mean3d_solar_anchor.csv" in SRC
    assert "gerber" not in SRC.split('"""')[2].lower() or "read_deck_node" not in SRC.split('"""')[2]


def test_the_sentinel_is_excluded_from_the_anchor():
    """One Fe I row in the release is exactly −4.0000 — a floor value, not a correction."""
    assert "SENTINEL = -4.0" in SRC
    assert "d == SENTINEL" in SRC


def test_the_vturb_bracket_is_carried_not_interpolated():
    """The Sun is vturb 1.0 and the grid has no such node. Interpolating would invent one;
    picking the end that passes would be tuning. Both ends are carried to the verdict."""
    assert 'VTURB_BRACKET = ("0.75", "1.50")' in SRC
    assert "STRADDLES THE TOLERANCE" in SRC


def test_an_ambiguous_wavelength_match_is_refused():
    """A rounded number is not an identity — two anchor lines inside the window is a
    refusal, never resolved by taking the nearer."""
    assert "AMBIGUOUS" in SRC
    assert "refused rather than resolved by proximity" in SRC


def test_the_tolerance_is_declared_as_inherited():
    """⚠️ A borrowed threshold is not a control (RYA-847). 0.05 comes from RYA-534 because
    the ticket says the gate LOGIC is unchanged; it is NOT re-derived on this data, and the
    measured margin is printed beside every verdict so the cut can be judged too."""
    assert "INHERITED FROM RYA-534, NOT RE-DERIVED" in SRC
    assert "tolerance_basis" in SRC


def test_both_legs_are_required():
    """Leg 1 can be satisfied by two codes sharing an error; leg 2 by a deck that is
    self-consistently wrong. Only both together is a validation."""
    assert "leg1_anchor_agreement" in SRC and "leg2_structure" in SRC
    assert 'doc["verdict"] = "VALIDATED" if both else "NOT_VALIDATED"' in SRC


def test_leg2_declares_its_envelope_as_declared():
    """It bounds the absurd; it does not certify the value. Saying so is the difference
    between a sanity check and a fabricated threshold."""
    assert "DECLARED, not derived" in SRC
    assert "bounds the absurd" in SRC

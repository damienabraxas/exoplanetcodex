"""RYA-1189 — the continuum RCA: is the UV/IR continuum placed too high?

DIAGNOSTIC. These tests pin what was MEASURED and, more importantly, pin the two ways
this measurement could be misread — because both are easy and both point at a redo that
would make the numbers worse.

🔴 MISREAD 1: "the local continuum sits below unity, so the current one is too high."
    EW = integral(1 - f/cont). The 95th percentile of a real absorbed spectrum is below
    unity essentially everywhere, so cont_local < 1 = cont_current is FORCED, and so is
    dA < 0. The direction of travel is arithmetic. It is not evidence.

🔴 MISREAD 2: "the near-UV shift is the biggest, so fix the near-UV first."
    The near-UV shift is the biggest AND the least meaningful: its side-band carries 1.65
    of catalogued central depth over 59 lines/A, and 0 of 10 side-bands pass
    SIDEBAND_CLEAN_MIN. `cont_local` there is a PSEUDO-continuum. Placing the continuum on
    it removes real line flux and biases EW LOW — the failure `band_products` already
    documents from measurement (Fe I 6910-9199 lines whose side-bands sat at 0.90/0.94
    lost 71% and 60% of their EW to re-normalisation).
"""
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
PERLINE = ROOT / "data/results/rya1189/rya1189_per_line.csv"
PERBAND = ROOT / "data/results/rya1189/rya1189_per_band.csv"
DOC = ROOT / "data/results/rya1189/rya1189_continuum_rca.json"


@pytest.fixture(scope="module")
def doc():
    if not DOC.exists():
        pytest.skip("RYA-1189 artifact absent")
    return json.loads(DOC.read_text())


@pytest.fixture(scope="module")
def band():
    if not PERBAND.exists():
        pytest.skip("RYA-1189 artifact absent")
    return pd.read_csv(PERBAND).set_index("band")


@pytest.fixture(scope="module")
def lines():
    if not PERLINE.exists():
        pytest.skip("RYA-1189 artifact absent")
    return pd.read_csv(PERLINE)


def test_the_control_band_barely_moves(band):
    """VIS is the control: we trust it, and the ticket says it should barely move. If VIS
    ever shows a material shift, the METHOD is suspect and no other band's number is
    readable."""
    vis = band.loc["VIS"]
    assert vis["n"] >= 8
    assert abs(vis["mean_delta_pct"]) < 3.0, "the control moved — method suspect"
    assert abs(vis["mean_d_A_dex"]) < 0.05
    assert vis["n_sideband_clean"] == vis["n"], (
        "a VIS side-band stopped passing SIDEBAND_CLEAN_MIN — the control's side-bands "
        "being clean is what makes it a control")


def test_the_near_uv_shift_is_the_largest_AND_the_least_meaningful(band):
    """Both halves, together, or the row invites exactly the wrong redo."""
    uv = band.loc["near-UV"]
    assert uv["mean_delta_pct"] > 10.0, "the headline shift"
    assert uv["n_sideband_clean"] == 0, "not one near-UV side-band is a continuum"
    assert uv["mean_blanket_frac"] == 1.0, "every side-band pixel is absorbed"
    assert uv["mean_sideband_catalogued_absorption"] > 1.0, (
        "the catalogued absorption that EXPLAINS the shift is gone — re-read the verdict")
    assert "BLEND-DRIVEN" in json.loads(DOC.read_text())["verdict_per_band"]["near-UV"]


def test_the_near_uv_has_no_isolated_lines_at_all(lines):
    """🔴 THE FINDING BEHIND THE VERDICT, and it is why the ticket's own spec ("6-10 clean,
    isolated Fe I lines in EACH band") is not satisfiable in the near-UV.

    Measured over the lab-graded Fe I set, the deepest catalogued neighbour within 0.6 A:
        near-UV  n=59   min 0.790  median 0.966
        VIS      n=176  min 0.001  median 0.540
    Not one near-UV line has a neighbour shallower than 0.79. A per-band continuum needs
    LINE-FREE WINDOWS; the near-UV has none. The first cut of the script used isolation as
    a FILTER and returned zero near-UV lines — an absence manufactured by the cut."""
    uv = lines[(lines.band == "near-UV") & (lines.status == "ok")]
    vis = lines[(lines.band == "VIS") & (lines.status == "ok")]
    assert len(uv) >= 8, "the near-UV must be MEASURED, never filtered to empty"
    # the least-crowded near-UV line is still more blended than the WORST VIS one
    assert uv.blend_frac.min() > vis.blend_frac.max(), (
        "the bands stopped separating on blending — re-read the census")
    # ⚠️ On DENSITY the claim is distributional, not min-vs-max: one VIS line sits in a
    # locally crowded spot (22.2 lines/A) above the least-crowded near-UV one (16.7), so
    # a min > 5*max assertion is simply false. Medians are what separate, ~9x.
    assert (uv.sideband_line_density_per_A.median()
            > 5 * vis.sideband_line_density_per_A.median())
    assert uv.sideband_line_density_per_A.min() > vis.sideband_line_density_per_A.median()


def test_the_direction_of_dA_is_not_treated_as_evidence(doc):
    """MISREAD 1, pinned in the artifact's own words. cont_local < 1 forces dA < 0 in every
    band, so "a per-band continuum moves the UV toward 7.466" is arithmetic, not support."""
    note = doc["anchor_check"]["note"]
    assert "arithmetic" in note and "NOT evidence" in note
    assert doc["anchor_check"]["anchor"] == 7.466
    assert "used for nothing" in doc["anchor_check"]["what"]


def test_not_blend_driven_is_not_reported_as_a_misplaced_continuum(doc):
    """⚠️ Excluding catalogued blending is not the same as establishing a continuum error.
    Both flagged bands are telluric-heavy, so uncatalogued absorption stays live."""
    for b in ("red-optical", "NIR-H"):
        v = doc["verdict_per_band"][b]
        assert v.startswith("NOT BLEND-DRIVEN")
        assert "does NOT establish a misplaced continuum" in v
        assert "telluric residual" in v and "molecular opacity" in v


def test_the_abundance_shift_is_labelled_a_lower_bound_where_saturated(lines):
    """dA = log10(EW ratio) is EXACT on the linear curve of growth and an UNDERESTIMATE of
    |dA| on the saturated part. Every near-UV line here is saturated (REW > -5), so the
    near-UV dA is a lower bound — which matters only for how the number is quoted, since
    the band's verdict is blend-driven anyway."""
    uv = lines[(lines.band == "near-UV") & (lines.status == "ok")]
    assert uv.saturated.all(), "the near-UV pool stopped being saturated — requote dA"
    src = (ROOT / "scripts/rya1189_continuum_rca.py").read_text()
    assert "LOWER BOUND ON |dA|" in src
    assert "MOOGSILENT" in src, "the owed true-inversion route is no longer named"


def test_the_diagnostic_writes_ONLY_into_its_own_results_directory():
    """The ticket is explicit: measure, report, change nothing.

    ⚠️ The first version of this test banned the STRING "data/products" and duly failed
    the moment the script started READING the feed for the closing anchor check -- which
    is exactly what the ticket asks for. Reading a published value is the job; writing one
    is the thing to forbid. So the guard is now on the WRITE SITES, parsed rather than
    grepped: every `to_csv` / `write_text` target must resolve through `out_dir`.
    """
    import ast as _ast
    src = (ROOT / "scripts/rya1189_continuum_rca.py").read_text()
    tree = _ast.parse(src)
    # ⚠️ The TARGET is a different node per method, and getting that wrong is how this
    # guard first failed on its own correct code: `to_csv(path)` takes the path as its
    # ARGUMENT, while `write_text(content)` takes the CONTENT and the path is the
    # RECEIVER. Reading args[0] for both asks whether the JSON payload is under out_dir.
    PATH_IS_ARG = {"to_csv", "to_json", "savefig"}
    PATH_IS_RECEIVER = {"write_text", "write_bytes", "open"}
    writes = []
    for node in _ast.walk(tree):
        if not (isinstance(node, _ast.Call) and isinstance(node.func, _ast.Attribute)):
            continue
        if node.func.attr in PATH_IS_ARG and node.args:
            writes.append(_ast.unparse(node.args[0]))
        elif node.func.attr in PATH_IS_RECEIVER:
            writes.append(_ast.unparse(node.func.value))
    assert writes, "no write sites found — the parser stopped seeing them"
    for w in writes:
        assert "out_dir" in w, f"writes outside its own results directory: {w}"

    # and it must never place a new continuum on a holding
    assert "fit_continuum" not in src, "this diagnostic must never place a new continuum"
    assert "spectra_normalize" not in src


def test_every_band_in_scope_is_pre_normalised_or_the_baseline_is_wrong(doc):
    """`cont_current = 1.0` is only what the pipeline does on a PRE-NORMALISED holding. On
    any other holding that baseline would be fiction, so the script refuses rather than
    measuring a shift against a continuum the holding never had."""
    src = (ROOT / "scripts/rya1189_continuum_rca.py").read_text()
    assert "is_pre_normalised" in src and "Refusing rather than measuring a shift" in src
    assert doc["placements"]["current"].startswith("cont = 1.0")

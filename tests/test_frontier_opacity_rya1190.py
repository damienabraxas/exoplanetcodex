"""RYA-1190 — frontier-band uncatalogued opacity: what is resolvable, and what is not.

DIAGNOSTIC. These tests pin the measurement and the three ways it could be misread. All
three are live: the ticket is written on the premise that completing the line list is the
lever, and the measurement says it is not.
"""
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "data/results/rya1190/rya1190_frontier_opacity.json"
PERBAND = ROOT / "data/results/rya1190/rya1190_per_band.csv"
SRC = ROOT / "scripts/rya1190_frontier_opacity.py"


@pytest.fixture(scope="module")
def doc():
    if not DOC.exists():
        pytest.skip("RYA-1190 artifact absent")
    return json.loads(DOC.read_text())


def test_the_production_list_IS_the_vald_extract(doc):
    """🔴 PART A'S PREMISE, REFUTED BY MEASUREMENT. The ticket asks which blends are
    "resolvable by completing the line list" from a VALD3 extract-all. Our near-UV list
    ALREADY IS that extract: 0 lines in VALD only, 0 in ours only, over 62,072 rows. There
    is nothing to complete at the extraction's threshold."""
    uv = doc["part_A_vald_completeness"]["near-UV"]
    assert uv["available"] is True
    assert uv["in_vald_only"] == 0, (
        f"{uv['in_vald_only']} VALD lines are missing from production — the premise is "
        f"back and the payoff must be re-estimated")
    assert uv["in_production_only"] == 0
    assert uv["vald_rows"] == uv["production_rows"] > 50000
    assert uv["extraction_depth_threshold"] == 0.001, (
        "the extraction threshold changed — the weak-line-haze bound is keyed to it")


def test_every_band_in_scope_is_already_vald_complete(doc):
    """Not just the near-UV. If a band ever shows VALD-only lines, that band DOES have a
    completeness lever and this ticket's answer does not apply to it."""
    for band, v in doc["part_A_vald_completeness"].items():
        if not v.get("available"):
            continue
        assert v["in_vald_only"] == 0, f"{band} has {v['in_vald_only']} un-ingested VALD lines"


def test_the_molecular_surplus_is_the_same_subthreshold_tail_not_a_new_source(doc):
    """⚠️ MISREAD 1: "the molecular lists reach below 3780 A, so add them."

    They do reach there — CH MoLLIST down to 2751 A — and reading the filenames alone
    (`12C14N_400-450.bsyn`, nanometres) would have said the opposite. But the sub-3780
    content is extreme vibrational overtone (`v12-0`, `v19-0`) at log gf -10 to -14 for
    OH and NH, and the ELECTRONIC region-split files all start at 4200 A. What CH does
    offer is material VALD itself evaluated (its reference footer credits "Masseron 2014
    obs: CH") and dropped below the same 0.001 threshold. Same population, not a new one.
    """
    mol = doc["part_A_molecular_coverage"]
    assert mol["lowest_wavelength_A"] < 3780.0, (
        "the coverage claim is now the other way round — re-read the finding rather than "
        "inverting the conclusion")
    src = SRC.read_text()
    assert "NANOMETRES" in src, "the filename/unit trap is no longer recorded"


def test_the_completeness_payoff_is_bounded_and_small(doc):
    """The number the recommendation rests on. Both factors are upper bounds and the
    artifact must say so, because a payoff quoted without that reads as a plan."""
    p = doc["part_A_payoff_estimate"]
    assert p["max_catalogueable_fraction"] < 0.5, (
        "the haze can now explain most of the deficit — the recommendation flips")
    assert abs(p["estimated_payoff_dex"]) < 0.05
    assert p["published_gap_to_anchor_dex"] == [0.13, 0.19]
    assert "UPPER BOUND" in p["reading"] and "ORDER, not a number to plan on" in p["reading"]


def test_the_red_optical_verdict_names_UNCORRECTED_not_residual_telluric(doc):
    """🔴 The distinction decides the fix. "Residual after molecfit" presumes a correction
    happened; in 7 of 10 red-optical windows the corrected product is BYTE-IDENTICAL to
    the raw atlas, and those untouched windows are the MORE depressed ones. The lever is
    extending the correction, not re-placing a continuum."""
    v = doc["verdict_per_band"]["red-optical"]
    assert v.startswith("OPACITY")
    assert "UNCORRECTED telluric, not residual-after-molecfit" in v
    b = pd.read_csv(PERBAND).set_index("band").loc["red-optical"]
    assert b["n_windows_uncorrected"] > b["n_windows_telluric_corrected"]
    assert b["residual_where_untouched"] < b["residual_where_corrected"], (
        "the untouched windows are no longer the more depressed ones — the attribution "
        "rested on that ordering")


def test_the_nir_verdict_refuses_to_classify_a_sigma_artefact(doc):
    """⚠️ MISREAD 2: "NIR-H shows the biggest IR deficit, so fix the IR continuum."

    Its deficit CHANGES SIGN across the width sweep (-0.057 at the physical sigma, +0.099
    at sigma x2), because CRIRES+ at R~70k has an instrumental width comparable to the
    Doppler one and the accounting goes insensitive. UNDETERMINED is the honest answer."""
    v = doc["verdict_per_band"]["NIR-H"]
    assert v.startswith("UNDETERMINED")
    assert "CHANGES SIGN" in v
    assert "reading the sigma, not the data" in v


def test_the_width_carries_the_instrument_not_just_the_doppler_term(doc):
    """🔴 THE ERROR THAT NEARLY MISCLASSIFIED THE IR, pinned. A first cut used the Doppler
    width alone for every band — a 2% correction at Kitt Peak (R 300-500k) and a ~46%
    one at CRIRES+ (R 50-100k). Under-broadening makes the catalogue absorb less than it
    should and manufactures a deficit."""
    src = SRC.read_text()
    assert "sigma_inst" in src and "resolving_power" in src
    assert "OMITTING IT NEARLY COST THE NIR VERDICT" in src
    rp = doc["method"]["resolving_powers_used"]
    assert rp["kpno_solar_atlas"] > 3 * rp["crires_plus"], (
        "the two instruments no longer differ in resolution the way this correction "
        "assumes — re-derive it")


def test_the_control_floor_is_subtracted_and_not_explained_away(doc):
    """The VIS control's side-bands are line-free by selection and it STILL sits ~1.1%
    below the accounting. This cannot tell a method floor from a small real continuum
    offset, so it subtracts the floor and says so rather than picking one."""
    floor = doc["method"]["control_floor_subtracted"]
    assert -0.03 < floor < 0.0
    assert "NOT radiative transfer" in SRC.read_text() or "NOT a synthesis differential" in SRC.read_text()
    b = pd.read_csv(PERBAND).set_index("band")
    assert abs(b.loc["VIS", "excess_over_control"]) < 1e-9, "the control must define zero"


def test_it_is_not_presented_as_a_synthesis(doc):
    """⚠️ MISREAD 3: quoting this as "the synthetic". It is a catalogued-opacity
    accounting with no atmosphere and no radiative transfer, and the reason a real
    differential was not run must travel with it."""
    m = doc["method"]
    assert "NOT a synthesis differential" in m["what"]
    assert "bsyn_lu is absent" in m["why"] and "Linux ELF" in m["why"]
    assert "owed as a Sirius run" in m["why"]


def test_the_diagnostic_writes_ONLY_into_its_own_results_directory():
    """Diagnostic only: identify, classify, count, estimate. Guarded on the WRITE SITES,
    parsed rather than grepped, with the target node taken per method — `to_csv` takes the
    path as an argument, `write_text` takes the content and the path is the receiver."""
    import ast as _ast
    tree = _ast.parse(SRC.read_text())
    PATH_IS_ARG = {"to_csv", "to_json", "savefig"}
    PATH_IS_RECEIVER = {"write_text", "write_bytes"}
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

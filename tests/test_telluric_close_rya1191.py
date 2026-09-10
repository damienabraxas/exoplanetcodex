"""RYA-1191 — the telluric close for Fe: what was verified, what was fixed, what was not.

VERIFICATION + FIX. These pin the four things that would otherwise be misread:

  * a line SET re-derived rather than re-measured — worthless unless it is controlled
    against the artifacts that exist, so that control is asserted here first;
  * a "raw-like" band read as a correction gap when the correction had nothing to remove,
    or read as clean when the test was blind to a shallow one;
  * a quarantine treated as conservative when it is throwing away measurable lines;
  * a per-window CRIRES+ verdict carried onto lines the window cannot reach.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "data/results/rya1191"
DOC_A = R / "rya1191_graded_line_telluric.json"
CSV_A = R / "rya1191_graded_line_telluric.csv"
DOC_B = R / "rya1191_telluric_residual.json"
DOC_C = R / "rya1191_crires_h_lines.json"
EVID = ROOT / "data/catalog/telluric_correction_evidence.csv"


def _load(p):
    if not p.exists():
        pytest.skip(f"{p.name} absent")
    return json.loads(p.read_text()) if p.suffix == ".json" else pd.read_csv(p)


@pytest.fixture(scope="module")
def doc_a(): return _load(DOC_A)
@pytest.fixture(scope="module")
def lines(): return _load(CSV_A)
@pytest.fixture(scope="module")
def doc_b(): return _load(DOC_B)
@pytest.fixture(scope="module")
def doc_c(): return _load(DOC_C)


# ── (A) the control that licenses re-deriving a line set ─────────────────────────────
def test_the_line_selection_is_controlled_against_every_artifact_that_exists(doc_a):
    """🔴 THE LINCHPIN. The red-optical and NIR line sets are re-derived by calling the
    production selectors, because no `_lines.csv` exists for those products. That is only
    trustworthy because the same call reproduces every graded artifact that DOES exist —
    line for line, not in count. If this ever drops below 100%, every verdict downstream
    is a guess about which lines were measured."""
    c = doc_a["selection_control"]
    assert c["checked"] >= 30, f"only {c['checked']} artifacts checked — too thin a control"
    assert c["exact"] == c["checked"], f"MISMATCHED: {c['mismatched']}"
    assert not c["mismatched"]


def test_every_graded_line_has_a_verdict_and_none_is_unverifiable_for_want_of_an_artifact(lines):
    """The ticket's smoke test. RYA-1192 could not speak for red-optical/NIR because the
    per-line artifacts do not exist; every graded line is judged here, and the few that
    are UNVERIFIABLE are so for a NAMED reason about the data, never a missing file."""
    assert len(lines) > 1000
    assert lines.flux_state.notna().all()
    u = lines[lines.flux_state == "UNVERIFIABLE"]
    assert (u.flux_basis.str.contains("DIRECT-COMPARISON")).all(), (
        "an UNVERIFIABLE line must name the measurement that failed, not a missing artifact")
    for band in ("red-optical", "NIR"):
        assert (lines.band == band).any(), f"{band} is still unrepresented"


def test_no_graded_line_carries_an_abundance_on_verified_raw_telluric_flux(doc_a):
    """🔴 THE HEADLINE, AND IT IS THE REASSURING DIRECTION. Across every band, not one
    graded Fe line that production KEPT sits inside a registered telluric band on flux
    verified NOT corrected. The KP label defect is a provenance defect; it never reached
    an abundance."""
    assert doc_a["disagreements"]["kept_though_the_flux_is_RAW"]["n"] == 0


# ── (B) telluric residual vs uncatalogued blend ──────────────────────────────────────
def test_the_o2gamma_verdict_is_WITHDRAWN_and_says_why(doc_b):
    """🔴 I REPORTED THIS ONE AND IT WAS WRONG. Leg B called o2gamma "harps essentially
    UNCORRECTED, +12.0 sigma". The template it was measured against is the ratio of two
    SEPARATELY REDUCED IAG products, and that ratio is not flat where there is no
    telluric: in the clean 6400-6500 A window it reaches 0.476 minimum transmission with
    5.7% of pixels below 0.95. o2gamma's own template is 11.8% — not clearly above it.

    molecfit settles it independently: fitted on HARPS's own night with real GDAS, its
    PHYSICAL model puts the deepest O2 gamma pixel at 1.78% absorption against this
    template's 82.76%, a factor of 47, with only r = +0.21 shared. The +12 sigma was
    HARPS correlating with a Baker-vs-Reiners REDUCTION difference."""
    b = next(x for x in doc_b["bands"] if x["lo_A"] == 6270.0)
    assert b["template_clears_the_clean_window_null"] is False
    assert b["verdict"].startswith("NO TEMPLATE")
    assert "NOT a clean verdict" in b["verdict"]
    assert not b["gaps"], "a withdrawn band must claim no gaps"


def test_the_clean_window_null_exists_and_is_not_flat(doc_b):
    """⚠️ THE NULL LEG B DID NOT HAVE. A displaced null controls for REGISTRATION — it
    cannot notice that the template is not telluric to begin with. This one measures the
    template where no telluric complex exists, and it is emphatically not zero."""
    c = doc_b["clean_window_null"]
    assert c["worst_frac"] > 0.02, (
        "if the template were now flat in clean windows, the o2gamma withdrawal needs "
        "re-deriving")
    assert len(c["windows"]) >= 3
    assert "NOT flat" in c["note"]


def test_the_deep_bands_are_unaffected_by_the_withdrawal(doc_b):
    """The damage is bounded. O2 B, H2O 7160-7340 and H2O 9280-9600 sit one to two orders
    above the clean-window level, so their verdicts rest on real correction structure and
    stand."""
    for lo in (6867.0, 7160.0, 9280.0):
        b = next(x for x in doc_b["bands"] if x["lo_A"] == lo)
        assert b["template_clears_the_clean_window_null"] is True, lo


def test_a_difference_cannot_tell_a_rescale_from_a_line_removal(doc_b):
    """🔴 THE OTHER HALF OF THE SAME MISTAKE. RYA-1192 returned VERIFIED-CORRECTED for
    kpno_molecfit at o2gamma on a mean |diff| of 0.030 against a zero floor. Measured
    against the telluric line positions it still carries a residual — the 0.030 is a
    featureless rescale, and the O2 gamma lines are where they were."""
    b = next(x for x in doc_b["bands"] if x["lo_A"] == 6270.0)
    assert b["holdings"]["solar_kpno_molecfit_corrected"]["residual_z"] >= 3.0
    s = doc_b["supersedes"]
    assert any("o2gamma" in f.lower() or "O2 gamma" in f for f in s["findings_rejudged"])
    assert "rescale" in s["why_the_earlier_tests_missed_it"] or \
           any("rescale" in f for f in s["findings_rejudged"])


def test_kurucz2005_is_clean_in_H2O_7160_and_rya1192_was_wrong(doc_b):
    """🔴 THIS LEG REFUTES MY OWN EARLIER FINDING. RYA-1192 called
    solar_kpno_kurucz2005_corrected RAW-LIKE in H2O 7160-7340 from a depth ratio of 4.2
    against a 5x bar — the wrong comparand. Against the telluric template it is at +1.0
    sigma while the two uncorrected references are at +28 and +32."""
    b = next(x for x in doc_b["bands"] if x["lo_A"] == 7160.0)
    assert b["holdings"]["solar_kpno_kurucz2005_corrected"]["residual_z"] < 3.0
    for ref in ("solar_kpno", "solar_iag_reiners2016"):
        assert b["holdings"][ref]["residual_z"] > 20.0, "the positive control must be loud"
    assert not b["gaps"]


def test_which_test_answered_is_recorded_per_band(doc_b):
    """⚠️ Neither test is the better one. The template loses power where a band saturates
    (O2 A/B), depth is blind where it is shallow (o2gamma). Every band must say which one
    answered it, or two incomparable numbers get read as one map."""
    for b in doc_b["bands"]:
        if b.get("state") != "MEASURED":
            continue
        assert "judged_by" in b, b["band"]
        assert b["judged_by"] in (None, "TEMPLATE (line positions)",
                                  "DEPTH (template saturated)"), b["judged_by"]
    sat = next(x for x in doc_b["bands"] if x["lo_A"] == 7594.0)
    assert sat["judged_by"] == "DEPTH (template saturated)", (
        "the O2 A-band is saturated — a shape correlation cannot judge it")


def test_the_sign_convention_is_stated_because_it_is_easy_to_read_backwards(doc_b):
    """The template is raw/corrected, so it dips BELOW 1 at telluric lines and a holding
    that still carries them correlates NEGATIVELY. `residual_r` is pre-flipped; if that
    ever stops being documented, every verdict inverts silently."""
    assert "LARGER POSITIVE = MORE RESIDUAL" in doc_b["sign_convention"]
    b = next(x for x in doc_b["bands"] if x["lo_A"] == 7160.0)
    assert b["holdings"]["solar_kpno"]["residual_r"] > 0, (
        "an uncorrected reference must score POSITIVE after the flip")


# ── the fix ──────────────────────────────────────────────────────────────────────────
def test_the_evidence_table_is_generated_from_the_measurement_not_typed():
    """RYA-686: a hand-kept copy of a measurement drifts, and the stale side is the one
    that passes. The table is regenerable and `--check` diffs it."""
    if not EVID.exists():
        pytest.skip("evidence table absent")
    head = EVID.read_text().splitlines()[0]
    assert "GENERATED" in head
    assert "rya1191_build_telluric_evidence.py" in EVID.read_text()
    d = pd.read_csv(EVID, comment="#")
    assert {"holding_id", "band_name", "lo_A", "hi_A", "state"} <= set(d.columns)
    assert set(d.state) <= {"clean", "partial", "uncorrected", "undetermined"}


def test_partial_and_undetermined_both_refuse_and_are_distinguishable():
    """⚠️ A partial correction leaves a MEASURED residual; an undetermined band is one no
    test had power in. Neither grants corrected treatment, and only the first is a
    statement about the flux (RYA-833) — so they must not collapse into one word."""
    from pipeline import telluric_policy as tp
    if not EVID.exists():
        pytest.skip("evidence table absent")
    # o2gamma is now `undetermined` on every holding (its template was withdrawn), and
    # H2O 9280-9600 is `partial` on the molecfit holding — the two states this test exists
    # to keep apart, on real rows.
    st, prov = tp.verified_band_state("solar_harps_molecfit_corrected", 6285.0)
    assert st == "undetermined" and prov
    st2, _ = tp.verified_band_state("solar_kpno_molecfit_corrected", 9400.0)
    assert st2 == "partial"
    assert st != st2
    src = (ROOT / "pipeline/telluric_policy.py").read_text()
    assert "BOTH REFUSALS, AND THEY ARE NOT THE SAME REFUSAL" in src


def test_the_lift_is_granted_on_evidence_and_still_refuses_the_real_gap():
    """🔴 THE FIX, AND THE THING IT MUST NOT BREAK. The quarantine lift used to speak only
    for the Kitt Peak readers, so a provider-corrected holding could never be heard and 25
    graded red-optical lines were thrown away on flux measured clean. It now also hears a
    measurement — but H2O 7160-7340 on the molecfit holding got NO admissible RYA-940 fit
    and is NOT measured clean, so that quarantine must still stand."""
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    if not EVID.exists():
        pytest.skip("evidence table absent")
    try:
        from measure_band_ew import serves_corrected_flux
    except SystemExit as e:      # the Kitt Peak atlas is not staged on every machine
        pytest.skip(f"measure_band_ew unavailable here: {str(e)[:60]}")
    assert serves_corrected_flux("solar_kpno_kurucz2005_corrected", 7164.448), (
        "kurucz2005 is measured clean in H2O 7160-7340 — the lift must be granted")
    assert not serves_corrected_flux("solar_kpno_molecfit_corrected", 7164.448), (
        "🔴 the H2O 7160-7340 gap on the molecfit holding is real and must keep refusing")


def test_the_quarantine_no_longer_throws_away_measurably_clean_lines(doc_a):
    """After the fix, no graded line may be refused before fitting on a holding whose flux
    at that wavelength is verified corrected."""
    d = doc_a["disagreements"]["quarantined_though_the_flux_IS_corrected"]
    assert d["n"] == 0, d["detail"]


# ── (C) CRIRES+ H ────────────────────────────────────────────────────────────────────
def test_the_2p84_sigma_window_is_resolved_and_not_kept(doc_c):
    """Ryan's ruling: +2.84 sigma may not be carried as 'acceptable'. It is resolved by
    UNIT — the 15700-15800 A window contains no graded Fe line at all, so the statistic
    cannot reach any abundance."""
    w = doc_c["the_2p84_sigma_window"]
    assert w["n_graded_lines_in_window"] == 0
    assert "RESOLVED BY UNIT" in w["verdict"]


def test_a_100A_window_verdict_does_not_survive_contact_with_its_lines(doc_c):
    """🔴 THE UNIT WAS WRONG. An abundance is fitted over 1.89 A; a residual 40 A away
    cannot reach it. Two of the seven graded lines inside RYA-1192's flagged windows have
    NO telluric within 5 A of them at all."""
    s = doc_c["the_seven_lines_in_rya1192_flagged_windows"]
    assert s["n"] == 7
    assert s["by_cause"].get("NO TELLURIC HERE", 0) >= 2
    assert "PER LINE, NOT PER 100 A WINDOW" in doc_c["unit_note"]


def test_an_underpowered_per_line_test_says_so_rather_than_reporting_clean(doc_c):
    """⚠️ The per-line window buys the right unit and pays in power: the RAW frame — the
    positive control — clears 3 sigma at only a few lines. Those are UNDETERMINED, and
    must never be counted as clean (RYA-833). What settles the ticket is that no CRIRES+ H
    Fe product is live, so none of these lines carries a published abundance."""
    assert doc_c["by_cause"].get("UNDETERMINED", 0) > 0
    assert doc_c["product_is_live"] is False
    assert "PRECONDITION" in doc_c["power_note"]
    for r in doc_c["lines"]:
        if r.get("cause") == "UNDETERMINED":
            assert r.get("why"), "an UNDETERMINED line must name why"


# ── (E) provenance housekeeping ──────────────────────────────────────────────────────
def test_the_stale_zero_vesta_fits_claim_is_gone():
    """🔴 A SENTENCE ABOUT A FILESYSTEM IS TRUE ON A DATE. `normalize_vesta_ir.py` asserted
    that a sweep found ZERO Vesta FITS on Sirius, so molecfit 'had nothing to run on'.
    All 18 are staged, and every CRIRES+ telluric verdict inherited from that sentence."""
    src = (ROOT / "scripts/normalize_vesta_ir.py").read_text()
    assert "STALE" in src and "All 18 CRIRES+ Vesta IDPs are staged" in src
    assert "60.A-9051(A)" in src


def test_the_crires_evidence_is_the_measurement_not_a_band_it_cannot_reach():
    """🔴 The registry justified `telluric_applied=applied` with an O2 A-band statistic
    (7594-7685 A) for an arm that starts at 9796.5 A. The claim is true; that was never
    the evidence for it."""
    reg = pd.read_csv(ROOT / "data/catalog/holdings_manifest_registry.csv")
    note = str(reg.set_index("holding_id").loc["elgueta2026_vizier", "notes"])
    assert "cannot have come from these spectra" in note
    assert "60.A-9051(A)" in note and "MTRANS" in note
    src = (ROOT / "scripts/normalize_vesta_ir.py").read_text()
    assert "O2 A-band is 7594-7685 A" in src and "9796.5 A" in src


# ── the HARPS O2 gamma gap, and the fix's two near-misses ────────────────────────────
def _include_regions():
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        from rya931_molecfit_model import include_regions_A
    except SystemExit as e:                      # esorex absent on this machine
        pytest.skip(f"rya931 runner unavailable: {str(e)[:60]}")
    return include_regions_A


def test_the_harps_fit_window_is_derived_from_the_band_registry():
    """🔴 THE ROOT CAUSE, AND IT IS NOT A MOLECFIT FAILURE. RYA-931's intake manifest
    records `registered_telluric_bands_in_instrument_range: ['O2 B 6867-6884']` — one
    band, so one hardcoded window, 6857-6911 A. RYA-1193 later added o2gamma 6270-6300 to
    TELLURIC_BANDS and nothing re-scoped the fit, so HARPS carries the O2 gamma lines at
    +12 sigma IDENTICALLY in its raw and corrected products. A hardcoded window cannot
    notice the registry grew; a derived one follows it."""
    regions = _include_regions()()
    assert any(a <= 6270.0 and 6300.0 <= b for a, b in regions), (
        f"o2gamma 6270-6300 is not covered by the derived fit window: {regions}")
    # the hardcoded pair must be GONE as a live constant — a comment naming it is fine,
    # and is in fact where the defect is explained
    import ast as _ast
    tree = _ast.parse((ROOT / "scripts/rya931_molecfit_model.py").read_text())
    assigned = {t.id for n in tree.body if isinstance(n, _ast.Assign)
                for t in n.targets if isinstance(t, _ast.Name)}
    assert "INCLUDE_LO_A" not in assigned and "INCLUDE_HI_A" not in assigned, (
        "the hardcoded fit window is still a live module constant")
    assert "TELLURIC_BANDS" in (ROOT / "scripts/rya931_molecfit_model.py").read_text()


def test_the_derivation_never_shrinks_the_window_that_already_worked():
    """🔴 THE REGRESSION I ALMOST SHIPPED. TELLURIC_BANDS is a LINE-SELECTION registry —
    its edges say where to refuse a line, not where to constrain a fit. Deriving from it
    alone gave O2 B as 6857-6894, SHRINKING RYA-931's measured 6857-6911 by 17 A, because
    the registered 6867-6884 understates a band that runs past the S1D cutoff. The derived
    regions are unioned with the empirical ones, never substituted for them."""
    regions = _include_regions()()
    assert any(a <= 6857.0 and 6911.0 <= b for a, b in regions), (
        f"RYA-931's measured O2 B window 6857-6911 is no longer covered: {regions}")


def test_the_molecule_test_is_a_token_not_a_substring():
    """🔴 `"O2" in "CO2"` IS TRUE, and the substring version pulled the CO2 15700-16100
    band into a HARPS O2 fit — a band 9000 A outside the instrument. The test is on the
    band name's first token."""
    regions = _include_regions()()
    assert not any(a > 12000.0 for a, _ in regions), (
        f"an IR band leaked into the O2 fit window: {regions}")
    src = (ROOT / "scripts/rya931_molecfit_model.py").read_text()
    assert "FIRST TOKEN, NOT A SUBSTRING" in src


def test_the_fit_regions_are_clipped_to_what_the_exposure_covers():
    """A band outside the instrument's range must be dropped before molecfit sees it,
    not handed over as an empty interval — and the runner refuses outright rather than
    emitting a no-op 'correction' if nothing is left."""
    f = _include_regions()
    harps = f(3782.6, 6912.0)
    assert harps and all(a >= 3782.6 and b <= 6912.0 for a, b in harps)
    assert not any(a >= 7000.0 for a, _ in harps), harps
    src = (ROOT / "scripts/rya931_molecfit_model.py").read_text()
    assert "refusing to run a fit" in src


# ── the blend audit and the curated exclusion it earned ──────────────────────────────
@pytest.fixture(scope="module")
def blend(): return _load(R / "rya1191_blend_audit.json")


def test_the_blend_audit_rediscovers_the_known_artifact_before_judging_anything(blend):
    """🔴 THE CONTROL, AND IT HAS TO REDISCOVER *WHY*, NOT JUST FLAG. RYA-515 identified
    Fe II 4303.170 independently: A = 10.0, 36 of 65 transitions CH. This audit counts 66
    transitions in the window and names CH I as the top contaminant at 43.8% of the summed
    depth. A verdict of TARGET-NOT-DOMINANT alone would NOT be a pass — 164 of 353 lines
    carry it — so the control demands the molecular share and the CH identification."""
    c = blend["control"]
    assert c["found"] and c["control_passes"]
    assert str(c["top_contaminant"]).startswith("CH")
    assert c["molecular_share"] > 0.3, c
    assert 60 <= c["n_transitions"] <= 72, (
        f"window population {c['n_transitions']} no longer matches RYA-515's 65")


def test_the_synthesis_list_carries_no_molecules_and_that_is_the_root_cause(blend):
    """🔴 THE FINDING UNDER THE FINDING. The VIS/red-optical synthesis list is
    GESv6_ATOM_hfs_iso — 206 species, none molecular — while the stellar catalogue carries
    CH, CN, OH, C2, NH and MgH. In the G band CH is the LARGEST species present. So the
    fit could not represent the dominant absorber and had to put it in Fe. That makes
    4303.170 an artifact of the LINE LIST, which is why removing it is curation."""
    import pandas as _pd
    cat = _pd.read_csv(ROOT / "data/linelists/linelist_solar.csv", low_memory=False,
                       usecols=["element", "wavelength_air_A"])
    sp = set(cat.element.astype(str))
    assert {"CH", "CN", "OH", "C2"} <= sp, "the catalogue no longer carries molecules"
    g = cat[cat.wavelength_air_A.between(4290, 4315)].element.astype(str).value_counts()
    assert g.index[0] == "CH", f"CH is no longer the dominant G-band species: {g.head(3)}"
    b = blend["synthesis_blind_to_molecular"]
    assert b["n"] >= 1 and "ATOM-ONLY" in b["meaning"]


def test_a_worse_line_than_the_flagged_one_is_reported_not_quietly_dropped(blend):
    """⚠️ Fe I red-optical 8432.174 scores WORSE than the line RYA-515 flagged: molecular
    share 1.0000, all CN, and no catalogued Fe absorption in its fit window at all. It is
    REPORTED for adjudication, not dropped — zero catalogued Fe could equally mean the
    catalogue is incomplete there, and that is a different problem from a blend."""
    worst = blend["synthesis_blind_to_molecular"]["lines"][0]
    assert worst["wavelength_air_A"] == pytest.approx(8432.174, abs=0.01)
    assert worst["molecular_share"] > 0.9
    import sys
    sys.path.insert(0, str(ROOT))
    from pipeline.line_curation import excluded
    assert not excluded("Fe I", 8432.174), (
        "8432.174 needs adjudication before it is dropped, not a silent removal")


def test_the_curated_exclusion_is_evidenced_and_is_not_a_threshold(blend):
    """🔴 164 of 353 graded lines are TARGET-NOT-DOMINANT. Auto-dropping on that statistic
    would delete a third of the pool on a number nobody ratified — a quota dressed as a
    cut. Exclusions are individual calls, each carrying its evidence and the ticket that
    ruled, and the loader refuses a row without them."""
    import sys
    sys.path.insert(0, str(ROOT))
    from pipeline.line_curation import excluded, excluded_wavelengths, REGISTRY
    assert excluded("Fe II", 4303.170).startswith("CURATED-EXCLUSION (RYA-1191)")
    assert not excluded("Fe I", 4303.170), "the exclusion must be species-specific"
    n_dropped = len(excluded_wavelengths("Fe I")) + len(excluded_wavelengths("Fe II"))
    assert n_dropped < 10, f"{n_dropped} curated drops — that is a cut, not curation"
    assert blend["by_verdict"].get("TARGET-NOT-DOMINANT", 0) > n_dropped * 10, (
        "the audit must rank far more lines than the registry drops")
    src = REGISTRY.read_text()
    assert "NOT A THRESHOLD" in src


def test_the_curated_exclusion_reaches_the_selection():
    """A registry nothing consults is a note. The synthesis route drops curated lines
    beside the telluric quarantine, and says which kind of exclusion each was."""
    src = (ROOT / "scripts/derive_band_products.py").read_text()
    assert "from pipeline.line_curation import excluded as _curated" in src
    assert "CURATED OUT (blend/artifact, not telluric)" in src
    assert "SEPARATE FROM THE TELLURIC QUARANTINE" in src


# ── the twin test: telluric offset vs gf-limited width ───────────────────────────────
@pytest.fixture(scope="module")
def twin(): return _load(R / "rya1191_twin_compare.json")


def test_a_wide_bar_is_attributed_by_changing_the_SPECTRUM_not_the_line_list(twin):
    """🔴 THE DISCRIMINATOR. `solar_kpno_molecfit_corrected` NIR ENGINE-A ships sigma_stat
    0.497 with a note reading "gf-limited, IRREDUCIBLE, reducible only by laboratory gf".
    If that were true it would be a property of the LINE LIST and identical on both
    holdings — same lines, same gf, same atmosphere. The verified-clean twin is 7x
    tighter. A statistic that moves when you swap the spectrum and hold the line list
    fixed is a statistic about the spectrum."""
    src = (ROOT / "scripts/rya1191_twin_compare.py").read_text()
    assert "statistic about the SPECTRUM" in src
    assert twin["gate"]["max_offset_dex"] == 0.03
    assert twin["gate"]["max_scatter_ratio"] == 2.0


def test_the_shared_width_is_blamed_on_neither_holding(twin):
    """⚠️ RYAN'S OWN CORRECTION, AND IT MATTERS AS MUCH AS THE FINDING. The NIR width IS
    genuinely gf-limited — the clean twin is wide there too (0.159 on verified-clean
    flux). So a failing pair must separate the width the two SHARE, which belongs to the
    line list, from the EXCESS and the OFFSET, which are the only things a raw holding can
    be blamed for. Replacing one blanket story with another is the error here."""
    for p in twin["pairs"]:
        if p.get("state") != "MEASURED" or "FAILS" not in p["verdict"]:
            continue
        assert "shared_scatter" in p and "excess_scatter" in p
        assert p["excess_scatter"] <= p["scatter_kp"] + 1e-9
        assert "belongs to neither holding" in p["reading"]


def test_only_lines_in_BOTH_pools_are_differenced(twin):
    """⚠️ The two holdings exclude different lines, so an unmatched mean difference mixes a
    per-line offset with a different line SET. The comparison is the intersection, and the
    pool sizes are reported so the restriction is visible rather than assumed."""
    for p in twin["pairs"]:
        if p.get("state") != "MEASURED":
            continue
        assert p["n_matched"] <= min(p["n_kp"], p["n_iag"])
        assert p["n_matched"] >= 3


def test_sigma_stat_and_the_printed_scatter_are_never_compared(twin):
    """⚠️ `sigma_stat` is the printed line scatter over sqrt(n). Comparing one against the
    other manufactures a 5-8x discrepancy out of arithmetic alone (RYA-1084), which on
    this ticket would have looked exactly like the telluric finding being hunted."""
    assert "divided by\n" in twin["statistic_note"] or "sqrt(n)" in twin["statistic_note"]
    for p in twin["pairs"]:
        if p.get("state") != "MEASURED":
            continue
        for h in ("kp", "iag"):
            n = p["n_matched"]
            # ⚠️ Tolerance is set by the artifact's 4-decimal STORAGE, not by taste:
            # scatter is rounded before division, so the two disagree in the 5th decimal
            # by construction. The drift this guards against is a factor of sqrt(n) —
            # 5-8x here — which no rounding tolerance could hide.
            assert p[f"sigma_stat_{h}"] == pytest.approx(
                p[f"scatter_{h}"] / np.sqrt(n), abs=2e-4), (
                f"{h}: sigma_stat is not scatter/sqrt(n) — the two statistics have drifted")


def test_guarding_first_moves_red_optical_inside_the_gate(twin):
    """🔴 THE ORDER RYAN RULED ON (2026-09-07 20:01): guard first, THEN assess telluric —
    because a railed fit contaminates the offset and the width alike. It changes the
    answer. Unguarded, both red-optical pairs failed on offset (+0.119, +0.169) and one on
    width (x3.51). With `pipeline.fit_validity` applied to BOTH pools before matching they
    are +0.025 and +0.029 against a 0.03 gate, and every width ratio drops to <= 1.47.
    Two of the four rows Ryan scanned as garbage are not garbage."""
    assert "guard" in twin and "guard first" in twin["guard"]
    ro = [p for p in twin["pairs"]
          if p["band"] == "red-optical" and p.get("state") == "MEASURED"]
    assert len(ro) == 2
    for p in ro:
        assert abs(p["delta_A"]) <= 0.03, p
        assert p["scatter_ratio"] <= 2.0, p
        assert "PASSES" in p["verdict"], p["verdict"]
        assert p["n_dropped_kp"] + p["n_dropped_iag"] > 0, (
            "if the guard drops nothing here, this result is not about guarding")


def test_the_NIR_offset_survives_the_guard_and_is_the_real_signal(twin):
    """⚠️ AND THE WIDTH STORY DOES NOT SURVIVE IT. What is left after guarding is a value
    OFFSET in the NIR only — +0.049 dex (1D-LTE) and -0.112 (ENGINE-A), the latter
    matching Ryan's own post-guard figure exactly — with every width ratio at or below
    1.47. The "sigma 0.497, 6.9x IAG, therefore uncorrected telluric" reading was one
    non-convergent line on each side."""
    nir = [p for p in twin["pairs"] if p["band"] == "NIR" and p.get("state") == "MEASURED"]
    assert len(nir) == 2
    assert all(abs(p["delta_A"]) > 0.03 for p in nir), [p["delta_A"] for p in nir]
    ea = next(p for p in nir if p["treatment"] == "ENGINE-A")
    assert ea["delta_A"] == pytest.approx(-0.112, abs=0.005)
    assert ea["A_iag"] == pytest.approx(7.599, abs=0.002), (
        "the clean twin must land on its shipped value once guarded")
    for p in twin["pairs"]:
        if p.get("state") == "MEASURED":
            assert p["scatter_ratio"] <= 2.0, (
                f"{p['band']} {p['treatment']}: a width blow-up survived the guard")


def test_a_contaminated_line_outside_every_declared_band_is_reported(twin):
    """🔴 THE HOLE IN MY OWN MAP. 9012.075 lies in NO registered telluric band — it is
    between H2O 8100-8400 and 9280-9600 — and RYA-1191's per-line map therefore recorded
    it as "no telluric band declared, no correction expected". The twin test says the flux
    there is contaminated by 3.2 dex. The earlier headline (zero graded lines kept on
    verified-raw flux INSIDE a telluric band) was literally true and materially
    incomplete: the weak link is the REGISTRY, not the check (RYA-833)."""
    import sys
    sys.path.insert(0, str(ROOT))
    from pipeline.telluric_policy import TELLURIC_BANDS
    w = 9012.075
    assert not any(lo <= w <= hi for lo, hi, _ in TELLURIC_BANDS), (
        "9012.075 is now inside a declared band — re-derive this reasoning")
    # 🔴 AND IT IS NO LONGER IN THE TWIN OUTLIERS, BECAUSE THE GUARD REMOVED IT FIRST.
    # 9012.075 returned A = 10.166/10.180 with the pool-maximum red_chi2 of 530 — a
    # non-convergent fit, caught by `pipeline.fit_validity`, not a telluric offset. Its
    # 3.2 dex "discrepancy" was never evidence about the registry, and asserting it here
    # would have kept a refuted reading alive.
    assert not any(round(x["wavelength_air_A"], 3) == w for x in twin["worst_lines"]), (
        "9012.075 is back in the twin outliers — the guard is no longer removing it")
    val = _load(R / "rya1191_validity_impact.json")
    assert any(round(r["wavelength_air_A"], 3) == w for r in val["dropped_lines"])
    # ⚠️ AND THE TWO CASES ARE DIFFERENT, WHICH IS THE POINT. 9437.793 — the largest
    # discrepancy of all — IS inside a declared band (H2O 9280-9600, where RYA-1191
    # measured kpno_molecfit PARTIAL). So the pool splits into "declared and imperfectly
    # corrected" and "not declared at all", and only the second is a registry gap.
    assert any(lo <= 9437.793 <= hi for lo, hi, _ in TELLURIC_BANDS), (
        "9437.793 is no longer inside a declared band — the two cases have merged")


# ── the fit-validity bound ───────────────────────────────────────────────────────────
@pytest.fixture(scope="module")
def validity(): return _load(R / "rya1191_validity_impact.json")


def test_a_fit_that_returned_an_impossible_abundance_is_not_a_measurement():
    """🔴 RYAN, ON FINDING A = 4.53 IN A LIVE POOL: "that seems like a bad line". It is.
    9437.793 returns A = 4.539 on solar_iag and 10.988 on solar_kpno_molecfit_corrected —
    same line, same gf, 6.4 dex apart, red_chi2 152 and 246, `excluded_reason` blank on
    both. The optimiser railed and nothing stopped it entering the aggregate."""
    import sys
    sys.path.insert(0, str(ROOT))
    from pipeline.fit_validity import fit_is_physical
    assert not fit_is_physical(4.539) and not fit_is_physical(10.988)
    assert fit_is_physical(7.46) and fit_is_physical(7.0) and fit_is_physical(8.0)
    assert fit_is_physical(None) and fit_is_physical(float("nan")), (
        "absence is a different problem from an impossible value (RYA-833)")


def test_the_bound_is_too_wide_to_shape_a_result(validity):
    """⚠️ NOT AN OUTLIER CUT, AND THAT DISTINCTION IS THE WHOLE JUSTIFICATION. Dropping
    points for being far from the others is the RYA-981 error, and RYA-515 says so
    directly for Fe II. This bound spans a FACTOR OF ~1000 in iron against a line-to-line
    scatter of ~0.2 dex — it can catch non-convergence and nothing else."""
    b = validity["bound"]
    assert b["half_width_dex"] >= 1.0, "a tighter bound would start shaping results"
    assert b["admits"] == [b["centre_A"] - b["half_width_dex"],
                           b["centre_A"] + b["half_width_dex"]]
    assert validity["fraction_dropped"] < 0.03, validity["fraction_dropped"]


def test_two_guards_were_tried_first_and_are_recorded_as_having_failed():
    """The bound is on the abundance because the two obvious fit-quality guards do not
    work on this data: red_chi2 > 10 flags 1262 of 1366 in-aggregate lines (the chi2 is
    uncalibrated — HARPS's ERR column is all NaN), and missing fit diagnostics flags 334
    of which only 3 are bad. Recording the dead ends stops the next reader re-proposing
    them."""
    src = (ROOT / "pipeline/fit_validity.py").read_text()
    assert "1262 of 1366" in src and "334 lines of which only 3" in src
    assert "NOT AN OUTLIER CUT" in src


def test_the_existing_constraint_gate_is_described_accurately():
    """⚠️ A CORRECTION TO MY OWN FIRST ACCOUNT. I wrote that RYA-992's cut "was never wired
    in". It was: `synth_gof_cut` and `ARM_SCALE` are called from `pipeline/gf_empirical.py`.
    What RYA-992 deliberately did not do is threshold the band-product route — RYA-847's
    sweep refuted every candidate and `constraint_gate` documents `SYNTH_CONSTRAINT` as
    None PERMANENTLY. The real gap is that both arms of that gate ask whether the fit was
    CONSTRAINED and neither asks whether the answer is POSSIBLE, and a line with NaN
    `frac_rise_weaker` slips both."""
    cg = (ROOT / "pipeline/constraint_gate.py").read_text()
    assert "PERMANENTLY rather than pending" in cg, (
        "the permanence of SYNTH_CONSTRAINT=None is the fact this reasoning rests on")
    src = (ROOT / "pipeline/fit_validity.py").read_text()
    assert "gf_empirical" in src, "the correction must name where RYA-992's cut IS called"
    assert "DIFFERENT QUESTIONS" in src
    # arm 1 cannot fire without a frac_rise, which is how the bad line slipped through
    import sys
    sys.path.insert(0, str(ROOT))
    from pipeline.constraint_gate import verdict
    assert verdict({"frac_rise_weaker": None}).ok
    assert verdict({"frac_rise_weaker": float("nan")}).ok
    assert not verdict({"frac_rise_weaker": -0.1}).ok


def test_the_bound_rediscovers_the_line_ryan_named_and_the_blend_audit_flagged(validity):
    """🔴 THREE INDEPENDENT ROUTES, ONE LINE. Fe II 4303.170 was named by RYA-515 from its
    abundance, found by RYA-1191's blend audit from the CH content of its window, and is
    caught here by a bound that knows nothing about either — on a red_chi2 of 1212. The
    blend audit's 4286.863 and 8432.174 come out too."""
    got = {round(float(r["wavelength_air_A"]), 3) for r in validity["dropped_lines"]}
    assert 4303.170 in got, "the bound no longer catches the known artifact"
    assert {4286.863, 8432.174} <= got, "the blend-flagged lines are no longer caught"


def test_the_bound_restores_a_shipped_number_the_current_code_does_not(validity):
    """⚠️ THE CHECK THAT MAKES THIS A FIX RATHER THAN A PREFERENCE. solar_iag NIR ENGINE-A
    ships A = 7.599 / sigma_stat 0.072 / n = 6. The current code reproduces 7.581 / 0.433 /
    n = 7 because it keeps 9437.793 at A = 4.539. With the bound the product returns to the
    committed value exactly."""
    hit = [r for r in validity["products_affected"]
           if "9199_11081_iag" in r["artifact"] and r["n_before"] == 7]
    assert hit, "the IAG NIR ENGINE-A product is no longer in the affected set"
    r = hit[0]
    assert r["n_after"] == 6
    assert r["A_after"] == pytest.approx(7.599, abs=0.002)
    assert r["sigma_stat_after"] == pytest.approx(0.072, abs=0.002)
    assert "reproduces_a_shipped_number" in validity


def test_the_impact_is_reported_for_every_product_including_the_big_move(validity):
    """⚠️ THE MOVE IS NOT UNIFORMLY SMALL, AND SAYING SO MATTERS. On the graded Fe I pools
    the medians shift by at most 0.021 dex; on the 3-line Fe II VIS ENGINE-A pool one
    dropped line moves A by 0.130 dex and leaves n = 2, exactly the RYA-1031 floor.
    Quoting the Fe I figure alone would have understated it by 6x."""
    assert validity["max_abs_delta_A"] >= 0.1, (
        "if the max move is now small, re-derive — it was 0.140 dex")
    big = [r for r in validity["products_affected"] if abs(r["delta_A"]) > 0.1]
    assert big, "the large-move product must stay visible in the artifact"
    assert any(r["n_after"] <= 2 for r in big), (
        "a pool driven to the n=2 floor must be visible, not averaged away")


def test_the_guard_reaches_the_aggregation_and_says_what_it_is_not():
    """A bound nothing calls is a note. It runs after the status and constraint guards —
    both of which 9437.793 passed — and its excluded_reason states that it is rejecting a
    non-convergent fit, not an outlier."""
    src = (ROOT / "scripts/derive_band_products.py").read_text()
    assert "fit_is_physical(lm.abundance, a.element)" in src
    assert "NOT AN OUTLIER CUT" in src
    from pipeline.fit_validity import rejection_reason
    r = rejection_reason(4.539)
    assert r.startswith("FIT-NOT-PHYSICAL") and "NON-CONVERGENT FIT, not an outlier" in r


# ── the provenance note that cited a bad fit as its evidence ─────────────────────────
def test_the_irreducible_dispersion_note_no_longer_rests_on_a_bad_fit():
    """🔴 THE NOTE NAMED ITS OWN COUNTER-EXAMPLE AS PROOF. It read "the widest bars carry
    the LARGEST n (KP ENGINE-A 1.315 at n=7 is the exception that proves it — CRIRES+
    curated lab-gf is ~6x tighter at n=5)". That 1.315 is sigma_stat 0.497 x sqrt(7), and
    it is ONE non-convergent fit: 9437.793 at A = 10.988 here against 4.539 on solar_iag.
    Guarded it is 0.144 at n=6 — indistinguishable from CRIRES+'s 0.138, so the ~6x ratio
    was that single line."""
    import json as _json
    feed = _json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    raw = (ROOT / "data/products/solar/Fe.json").read_text()
    assert "IRREDUCIBLE under the current NIR line list. It tracks gf QUALITY" not in raw, (
        "the refuted claim is back in the feed")
    # 🔴 COUNTED OVER THE LIVE PRODUCTS, NOT THE RAW FILE, AND NOW 10 NOT 6.
    # RYA-1208 added four NIR `synth-1D-LTE-gerber` cells, each carrying an
    # irreducible_dispersion block like every other NIR product. A raw string count also
    # sweeps in the ARCHIVE, which grows whenever a product is superseded and has nothing
    # to do with the claim under test. `carrying == len(live)` is what actually asserts
    # "every live block carries the correction"; the literal is pinned alongside so a
    # block VANISHING still fails rather than passing on a shrunken set.
    live = [p for p in feed["products"] if p.get("irreducible_dispersion")]
    carrying = sum(1 for p in live
                   if "RYA-1191 RE-DERIVED THIS NOTE"
                   in str(p["irreducible_dispersion"].get("note", "")))
    assert carrying == len(live) == 10, (
        f"every LIVE irreducible_dispersion block must carry the correction "
        f"({carrying} of {len(live)} do)")
    assert "NOT more lines\"" not in raw.replace(
        "PARTLY laboratory gf for NIR Fe I; a LARGE PART was non-convergent fits and has "
        "already been removed (RYA-1191)", ""), "a bare 'NOT more lines' claim survives"
    assert feed, "the feed must still parse"


def test_the_note_keeps_the_part_of_the_claim_that_survived():
    """⚠️ Ryan's instruction was "keep the true one, fix the false one" — not to replace
    one blanket story with another. Re-measured on guarded pools the many-line NIR bars
    ARE still the widest (IAG 0.307 at n=22, KP 0.386 at n=23) against the few-line
    curated ones (CRIRES+ 0.138 at n=5), so "tracks gf quality, not line count" stands.
    What does not stand is IRREDUCIBLE: it more than halved."""
    raw = (ROOT / "data/products/solar/Fe.json").read_text()
    assert "still tracks gf QUALITY rather than line count" in raw
    assert "would still widen it" in raw
    assert "it is NOT IRREDUCIBLE" in raw
    assert "0.895 -> 0.386" in raw and "0.949 -> 0.307" in raw


def test_the_stale_dispersion_number_is_flagged_not_silently_left():
    """⚠️ `dispersion_dex` is derived from the product's own sigma_stat and cannot be
    corrected by editing prose — it is regenerated when the product is. Saying so is the
    difference between a known-stale number and a wrong one (RYA-686)."""
    raw = (ROOT / "data/products/solar/Fe.json").read_text()
    assert "dispersion_dex above is the PRE-GUARD value" in raw

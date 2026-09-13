"""RYA-1214 — `cno_synthesis` is region-aware for the first time, and VIS must not move.

Ryan, 2026-09-12: run the molecular indicators in THIS ticket, near-UV and VIS first. The
VIS route was wired and validated; the near-UV region did not exist, and `run_cno` could
not have expressed one — `--region` chose the `RegionConfig` while three things stayed
pinned to the optical:

    diagnostics = VIS_DIAGNOSTICS                       (whatever the region said)
    obs_w, obs_f = _load_observed_spectrum(star_id)      (the HARPS loader)
    primary = {'C': by_key['CH_Gband'], ...}             (three VIS keys)

So a second region would have fitted VIS windows against a near-UV spectrum and reported
it under the new region's name — or, for the primaries, raised KeyError. These tests pin
the derivation and, more importantly, pin that VIS is UNCHANGED by it: every assertion
about VIS below is a statement that the generalisation is behaviour-preserving.

⚠️ NO SYNTHESIS RUNS HERE. iSpec is Sirius-only and a band fit is minutes; these are
structural, which is exactly the layer a refactor like this can break silently.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


@pytest.fixture(scope="module")
def cs():
    """Import `cno_synthesis`, skipping ONLY if iSpec genuinely is not installed.

    ⚠️ `importlib.util.find_spec("ispec")` IS THE WRONG TEST AND IT SKIPPED ALL TWELVE.
    iSpec is not on `sys.path` until a pipeline module puts it there —
    `abundances_derive.py:95` does `sys.path.insert(0, str(ISPEC_DIR))` before its own
    `import ispec` — so `find_spec` answers None on the very machine where iSpec works,
    and every test here reported `s` instead of running. A skipped test is not a passing
    test, and twelve of them silently not running is worse than not writing them.

    So the import is ATTEMPTED and only a real ImportError skips.
    """
    try:
        from pipeline import cno_synthesis
    except ImportError as exc:                  # pragma: no cover - env dependent
        pytest.skip(f"cno_synthesis not importable here: {exc}")
    return cno_synthesis


# ── the region registry ──────────────────────────────────────────────────────
def test_every_region_has_its_own_diagnostics(cs):
    """A region in REGIONS with no REGION_DIAGNOSTICS entry would KeyError at run time —
    after the atmosphere and line lists are loaded, which is minutes in."""
    assert set(cs.REGIONS) == set(cs.REGION_DIAGNOSTICS), (
        f"REGIONS {sorted(cs.REGIONS)} and REGION_DIAGNOSTICS "
        f"{sorted(cs.REGION_DIAGNOSTICS)} disagree")
    assert "nearuv" in cs.REGIONS and "vis" in cs.REGIONS


def test_the_nearuv_region_is_the_kitt_peak_near_uv_band(cs):
    r = cs.REGIONS["nearuv"]
    assert r.instrument == "kpno_solar_atlas"
    assert (r.wave_min_A, r.wave_max_A) == (3000.0, 3780.0)
    # The band this region claims must be the band_policy near-UV band, not a private one.
    from pipeline.band_policy import POLICIES
    nuv = next(p for p in POLICIES if p.name == "near-UV")
    assert r.wave_min_A == nuv.lo_A and r.wave_max_A <= nuv.hi_A


def test_vis_region_and_diagnostics_are_unchanged(cs):
    """The generalisation must be invisible to VIS."""
    assert cs.REGION_DIAGNOSTICS["vis"] is cs.VIS_DIAGNOSTICS
    keys = [d.key for d in cs.VIS_DIAGNOSTICS]
    # ⚠️ SIX, not seven. `OI_777` is a real Diagnostic but it belongs to
    # ESPRESSO_DIAGNOSTICS — the 777 triplet is not in the HARPS VIS band at all
    # (7771-7775 A is redward of 6910). Asserted from the tuple rather than from memory.
    assert keys == ["CH_Gband", "CI_5052", "CI_5380", "C2_Swan", "CN_red",
                    "OI_6300"], keys
    assert "OI_777" in [d.key for d in cs.ESPRESSO_DIAGNOSTICS]


# ── the derived primary ──────────────────────────────────────────────────────
def test_vis_primaries_derive_to_exactly_the_keys_that_were_hardcoded(cs):
    """The three names the old code named. If this fails the refactor moved VIS."""
    prim = cs.primary_by_element(cs.VIS_DIAGNOSTICS)
    assert {e: d.key for e, d in prim.items()} == {
        "C": "CH_Gband", "N": "CN_red", "O": "OI_6300"}


def test_a_region_may_measure_fewer_than_three_elements(cs):
    prim = cs.primary_by_element(cs.NEARUV_DIAGNOSTICS)
    assert {e: d.key for e, d in prim.items()} == {"N": "NH_AX", "O": "OH_AX"}
    assert "C" not in prim, (
        "the near-UV has no carbon primary — CH A-X (3870-4320 A) is redward of the "
        "3780 A edge and IS the VIS G-band")


def test_two_primaries_for_one_element_is_refused_not_resolved(cs):
    """Taking the first would make the region's own choice silently (RYA-780)."""
    d = cs.NEARUV_DIAGNOSTICS[0]
    import dataclasses
    twin = dataclasses.replace(d, key=d.key + "_twin")
    with pytest.raises(ValueError, match="two primary diagnostics"):
        cs.primary_by_element(cs.NEARUV_DIAGNOSTICS + (twin,))


# ── the near-UV windows are the measured ones ────────────────────────────────
def test_the_nearuv_windows_are_the_dominance_measured_ones(cs):
    """Chosen by the target molecule's dominance, and the two CH-dominated OH candidates
    were REFUSED — 3144-3147 (48 CH vs 19 OH) and 3180-3183 (26 CH vs 10 OH). Those hold
    MORE total molecular lines than either adopted OH window, so a wider-is-better rule
    would have taken exactly the wrong two."""
    by = {d.key: d for d in cs.NEARUV_DIAGNOSTICS}
    assert by["OH_AX"].windows_A == ((3063.0, 3066.0), (3122.0, 3125.0))
    assert by["NH_AX"].windows_A == ((3358.0, 3361.0), (3370.0, 3373.0))
    for w in by["OH_AX"].windows_A + by["NH_AX"].windows_A:
        assert not (3144.0 <= w[0] <= 3147.0), "3144-3147 is CH-dominated and was refused"
        assert not (3180.0 <= w[0] <= 3183.0), "3180-3183 is CH-dominated and was refused"


def test_every_nearuv_window_is_inside_the_region_band(cs):
    r = cs.REGIONS["nearuv"]
    for d in cs.NEARUV_DIAGNOSTICS:
        for lo, hi in d.windows_A:
            assert r.wave_min_A <= lo < hi <= r.wave_max_A, f"{d.key} window {lo}-{hi}"


def test_the_nearuv_bands_are_molecular_and_declare_their_coupling(cs):
    for d in cs.NEARUV_DIAGNOSTICS:
        assert d.kind == "molecular_band" and d.use_molecules is True
        assert "C" in d.depends_on, (
            f"{d.key}: carbon competes for both N and O through CN and CO, so A(C) must be "
            f"fixed before this band is fit")
        assert d.nlte_flag == "lte_molecular_band"


def test_the_nearuv_products_are_labelled_upper_bounds(cs):
    """The bias direction is known in advance: the near-UV synthesis under-corrects
    opacity (RYA-1204/1207/1189) and a fit compensates by RAISING the abundance. A number
    from this region that does not say so would be read as a measurement."""
    assert "UPPER BOUND" in cs.REGIONS["nearuv"].notes
    for d in cs.NEARUV_DIAGNOSTICS:
        assert "UPPER BOUND" in d.reference, f"{d.key} does not declare the bound"


# ── the loader dispatch ──────────────────────────────────────────────────────
def test_an_unwired_region_instrument_fails_loud(cs):
    """It must not fall through to the HARPS loader — that is what produced "the near-UV
    has no data" from a band whose atlas is right there."""
    bad = cs.RegionConfig(name="x", instrument="nonesuch", R=1e5, wave_min_A=1.0,
                          wave_max_A=2.0, telluric_correction_required=False,
                          nlte_backend="lte_by_design")
    with pytest.raises(cs.ArmNotWired, match="no spectrum loader"):
        cs._load_region_spectrum("solar", bad)


def test_the_kitt_peak_loader_is_solar_only(cs):
    with pytest.raises(cs.ArmNotWired, match="SOLAR flux atlas"):
        cs._load_kpno_atlas_arm("procyon", cs.REGIONS["nearuv"])


# ── the IR CN regions: where AGSS21's nitrogen actually is ────────────────────
def test_the_ir_cn_regions_exist_and_cover_the_agss21_band(cs):
    """AGSS21's CN indicator is the A-X band at 10872-13204 A. The optical `CN_red`
    diagnostic fits 6125/6195 A and ZERO of their 59 lines fall in it — same molecule,
    different band. These two regions are the band they actually used."""
    for name, lo, hi in (("nir_cn_kp", 10872.0, 13205.0),
                         ("nir_cn_iag", 10872.0, 11083.0)):
        r = cs.REGIONS[name]
        assert (r.wave_min_A, r.wave_max_A) == (lo, hi), name
        assert r.telluric_correction_required is True, (
            f"{name}: H2O 11120-11560 sits inside this band; the gate must be armed")
    # IAG stops at ITS OWN red edge, not a round number.
    assert cs.REGIONS["nir_cn_iag"].wave_max_A <= 11083.46


def test_the_ir_windows_are_agss21_own_line_positions(cs):
    """Not positions chosen by us. Clustering AGSS21's published CN wavelengths is what
    makes this a replication rather than a new selection (the `--lines-from-set` principle
    applied to a band)."""
    import pandas as pd
    src = (ROOT / "data/reference/amarsi2021_cno/derived"
           / "amarsi2021_cno_molecular_lines.csv")
    a = pd.read_csv(src)
    cn = a[(a.element_parameter == "logepsN") & (a.species == "CN")]
    lam = cn.wavelength_vac_nm * 10.0
    s2 = (1e4 / lam) ** 2
    air = lam / (1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2))
    air = sorted(float(x) for x in air)
    for name in ("nir_cn_kp", "nir_cn_iag"):
        d = cs.REGION_DIAGNOSTICS[name][0]
        for lo, hi in d.windows_A:
            assert any(lo <= w <= hi for w in air), (
                f"{name} window {lo}-{hi} A contains no AGSS21 CN line — the windows must "
                f"be their line positions, not ours")


def test_no_ir_fit_window_sits_inside_the_h2o_band(cs):
    """H2O 11120-11560 A swallows 5 of the 59 lines. EXCLUDED, not corrected: it is 440 A
    of the band's 2333 and the rest is enumerated-clean, so there is nothing to gain by
    leaning on a correction there (RYA-1193: reachability is per-observation)."""
    from pipeline.telluric_policy import TELLURIC_BANDS
    h2o = [(lo, hi) for lo, hi, n in TELLURIC_BANDS if lo == 11120.0]
    assert h2o, "the H2O 11120-11560 band is gone from TELLURIC_BANDS — re-derive this test"
    lo_h, hi_h = h2o[0]
    for name in ("nir_cn_kp", "nir_cn_iag"):
        for lo, hi in cs.REGION_DIAGNOSTICS[name][0].windows_A:
            assert not (hi > lo_h and lo < hi_h), f"{name} window {lo}-{hi} overlaps H2O"


def test_the_holding_resolver_is_one_decision_for_loader_and_gate(cs):
    """Two copies of this choice is how a gate clears one spectrum while the fit measures
    another (RYA-845). And the two Kitt Peak holdings are NOT interchangeable: the residual
    atlas stops at 10010 A, which is blueward of every AGSS21 CN line."""
    assert cs.holding_for_region(cs.REGIONS["nearuv"]) == "solar_kpno_kurucz2005_corrected"
    assert cs.holding_for_region(cs.REGIONS["nir_cn_kp"]) == "solar_kpno_molecfit_corrected"
    assert cs.holding_for_region(cs.REGIONS["nir_cn_iag"]) == "solar_iag"
    with pytest.raises(cs.ArmNotWired, match="no holding declared"):
        cs.holding_for_region(cs.RegionConfig(
            name="x", instrument="nonesuch", R=1e5, wave_min_A=1.0, wave_max_A=2.0,
            telluric_correction_required=False, nlte_backend="lte_by_design"))


def test_every_ir_holding_is_a_corrected_science_basis(cs):
    """The gate reads the holding's VERIFIED state, never a caller's flag (RYA-1026)."""
    from pipeline import telluric_display_policy as tdp
    for name in ("nir_cn_kp", "nir_cn_iag"):
        h = cs.holding_for_region(cs.REGIONS[name])
        assert tdp.display_state(h) in ("CLEAN", "CLEAN_WITH_ANOMALY"), (
            f"{name} -> {h} is not a corrected science basis")
    # solar_kpno reaches the same band and must NOT be the one chosen: CONTROL_ONLY.
    assert tdp.display_state("solar_kpno") == "CONTROL_ONLY"


def test_the_ir_cn_line_list_exists_and_covers_the_band():
    """Built from Brooke+2014; every 12C14N list we held stops at 9200 A."""
    d = ROOT / "data/linelists/molecular/turbospectrum/CN"
    ir = sorted(d.glob("12C14N_1087-*.bsyn"))
    assert ir, "the IR CN list is missing — AGSS21's CN band has no line list without it"
    head = ir[0].read_text().splitlines()[:2]
    assert "0607.012014" in head[0], f"wrong species code: {head[0]}"
    assert "Brooke" in head[1], f"source line does not name Brooke: {head[1]}"
    n = int(head[0].split()[-1])
    assert n > 5000, f"only {n} lines — too few to cover 10872-13205 A"

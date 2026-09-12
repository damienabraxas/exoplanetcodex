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

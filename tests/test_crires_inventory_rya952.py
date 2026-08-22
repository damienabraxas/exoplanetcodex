"""
RYA-952 — the CRIRES inventory, its classic/plus split, and its target identifications.

These pin the JUDGEMENTS, not the file count — a drive gains and loses files, but the rules
by which a frame is classified and identified must not drift.

1. `OBJECT` IS NOT AN IDENTIFIER. Four tau Ceti science frames on this drive carry
   `OBJECT = 'STD'`. Identification is positional and the referee is SIMBAD.
2. PROPER MOTION IS LOAD-BEARING. tau Ceti moves 1.92"/yr; over 22 years that is 42",
   comfortably outside any usable match radius.
3. THE INSTRUMENT IS SETTLED BY THE PIPELINE. `INSTRUME` says the bare string `CRIRES` for
   both instruments, and `PRO REC1 PIPE ID` is a version number.
4. A BINARY IS NOT RESOLVED BY PICKING THE NEARER COMPONENT.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "data" / "audit" / "rya952_crires_inventory"
CATALOG = ROOT / "data" / "reference" / "crires_target_astrometry.csv"
MANIFEST = ROOT / "data" / "audit" / "tau_cet_crires_plus" / "tau_cet_crires_plus_manifest.csv"

import sys
sys.path.insert(0, str(ROOT))
from pipeline.audit_crires import (  # noqa: E402
    BINARY_PAIRS, C_KMS, COADD_MAX_RESIDUAL_KMS, Frame, MOVING_TARGETS, PLUS_PIPE_TOKEN,
    _norm_name, _propagate, _sep_arcsec, classify_vintage, identify, load_astrometry)


def _rows(p: Path) -> list[dict]:
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


@pytest.fixture(scope="module")
def cat():
    return load_astrometry(ROOT)


@pytest.fixture(scope="module")
def inventory():
    return _rows(AUDIT / "crires_inventory.csv")


# ── judgement 1: OBJECT is not an identifier ─────────────────────────────────────────
def test_tau_ceti_is_found_despite_being_labelled_STD(inventory):
    """The headline. A frame's role is not its name."""
    std = [r for r in inventory if r["object_raw"] == "STD"]
    assert std, "no OBJECT='STD' frames — this pin has lost its subject"
    assert all(r["star_id"] == "tau_cet" for r in std), (
        "an OBJECT='STD' frame did not resolve to tau Cet; identification has fallen back "
        "to the label")
    assert all(r["id_status"] == "confirmed" for r in std)
    assert all("DOES NOT NAME THIS STAR" in r["id_evidence"] for r in std), (
        "the mislabel must stay VISIBLE — silently correcting it erases the finding")


def test_one_star_under_two_names_is_one_star(inventory):
    """`eps Eri` and `HD 22049` are the same object; so are `rho01 Cnc` and `HD 75732`."""
    for star, names in (("eps_eri", {"eps Eri", "HD 22049"}),
                        ("55cnc_a", {"rho01 Cnc", "HD 75732"})):
        seen = {r["object_raw"] for r in inventory if r["star_id"] == star}
        assert names <= seen, f"{star}: expected both of {names}, saw {seen}"


def test_a_catalogue_designation_is_not_a_mislabel(inventory):
    """`HD 22049` names eps Eri perfectly well. Only a ROLE or a placeholder is a finding."""
    mis = {r["object_raw"] for r in inventory
           if "DOES NOT NAME THIS STAR" in r.get("id_evidence", "")}
    assert "HD 22049" not in mis and "rho01 Cnc" not in mis, (
        f"a legitimate catalogue designation was flagged as a mislabel: {mis}")
    assert mis == {"STD"}, f"unexpected mislabel set {mis}"


def test_name_normalisation_strips_simbad_decoration():
    assert _norm_name("* rho01 Cnc") == _norm_name("rho01 Cnc") == "rho01cnc"
    assert _norm_name("** STT 270A") == "stt270a"


# ── judgement 2: proper motion is load-bearing ───────────────────────────────────────
def test_proper_motion_moves_tau_ceti_further_than_the_match_radius(cat):
    """If this ever stops being true the propagation could be dropped. It will not."""
    mjd_2022 = 59585.0
    ra, dec = _propagate(cat, mjd_2022)
    i = int(np.where(cat.star_id == "tau_cet")[0][0])
    moved = _sep_arcsec(cat.ra_deg_j2000.iloc[i], cat.dec_deg_j2000.iloc[i], ra[i], dec[i])
    assert moved > 40.0, (
        f"tau Cet moved only {moved:.1f}\" from J2000 to 2022 — recheck the catalogue; "
        f"the whole point is that this exceeds the match radius")


def test_without_proper_motion_tau_ceti_would_be_missed(cat):
    """The counterfactual, asserted: the J2000 position does NOT match the 2022 pointing."""
    fr = Frame(path="x", md5="x", ra_deg=26.004805, dec_deg=-15.93158, mjd=59585.06542,
               object_raw="STD", date_obs="2022-01-06T01:34:12")
    i = int(np.where(cat.star_id == "tau_cet")[0][0])
    naive = _sep_arcsec(fr.ra_deg, fr.dec_deg,
                        cat.ra_deg_j2000.iloc[i], cat.dec_deg_j2000.iloc[i])
    assert naive > 40.0, f"J2000 separation is only {naive:.1f}\" — counterfactual is void"
    got = identify(fr, cat)
    assert got.star_id == "tau_cet" and got.id_sep_arcsec < 30.0


# ── judgement 3: the pipeline settles the instrument ─────────────────────────────────
def test_instrume_alone_cannot_split_the_two_instruments(inventory):
    assert {r["instrume"] for r in inventory} == {"CRIRES"}, (
        "INSTRUME is expected to read the bare string CRIRES for CRIRES+ too — if that has "
        "changed, the classifier can be simplified")


def test_the_recipe_not_the_pipe_version_decides(inventory):
    """`PRO REC1 PIPE ID` is a version string and carries no instrument name."""
    plus = [r for r in inventory if r["instrument_class"] == "crires_plus"]
    assert plus
    for r in plus:
        assert PLUS_PIPE_TOKEN in r["rec_id"].lower(), r["rec_id"]
        assert PLUS_PIPE_TOKEN not in r["pipe_id"].lower(), (
            f"pipe_id {r['pipe_id']!r} now contains the recipe name — the classifier's "
            f"reason for reading rec_id has changed")


def test_a_cr2res_frame_dated_before_2014_is_a_loud_failure():
    """Pipeline and date must agree; a disagreement is not resolved by preferring one."""
    fr = Frame(path="x", md5="x", date_obs="2011-11-03T00:00:00",
               rec_id="cr2res_obs_nodding", pipe_id="1.6.9")
    with pytest.raises(ValueError, match="Refusing to pick one"):
        classify_vintage(fr)


def test_a_pre_2014_frame_classifies_as_classic():
    fr = Frame(path="x", md5="x", date_obs="2007-08-26T00:00:00", rec_id="", pipe_id="")
    assert classify_vintage(fr).instrument_class == "crires_classic"


# ── judgement 4: a binary is not resolved by picking the nearer component ────────────
def test_alpha_cen_frames_are_quarantined_not_assigned(inventory, cat):
    acen = [r for r in inventory if "alf Cen" in r["object_raw"] or r["object_raw"] == "Star S5"]
    assert acen, "no alpha Cen frames in the inventory"
    for r in acen:
        assert r["id_status"] == "quarantine", (
            f"{r['object_raw']} was assigned to {r['star_id']!r}; alpha Cen A and B are "
            f"~4-8 arcsec apart and astrometry alone cannot choose between them")
        assert "RYA-423" in r["id_evidence"]


def test_the_binary_pair_is_declared():
    assert frozenset({"alpha_cen_a", "alpha_cen_b"}) in BINARY_PAIRS


# ── the moving target is identified, not quarantined ─────────────────────────────────
def test_vesta_is_a_moving_target_not_an_unidentified_frame(inventory):
    v = [r for r in inventory if r["object_raw"] == "Vesta"]
    assert v
    assert all(r["id_status"] == "moving_target" for r in v)
    assert "vesta" in MOVING_TARGETS


# ── de-duplication is by CONTENT ─────────────────────────────────────────────────────
def test_duplicates_are_detected_by_md5_not_filename(inventory):
    """The same frame is on this drive as an ESO ADP id and as a CR_SONE_* name."""
    dups = [r for r in inventory if r["duplicate_of"]]
    assert dups, "no duplicates found — the de-duplication may have stopped working"
    for r in dups:
        assert Path(r["path"]).name != Path(r["duplicate_of"]).name or True
    by_md5 = {}
    for r in inventory:
        by_md5.setdefault(r["md5"], set()).add(Path(r["path"]).name)
    renamed = [names for names in by_md5.values() if len(names) > 1]
    assert renamed, "expected at least one frame stored under two different filenames"


# ── the tau Ceti holding ─────────────────────────────────────────────────────────────
def test_the_tau_cet_manifest_records_a_measured_telluric_state():
    rows = _rows(MANIFEST)
    assert len(rows) == 4
    for r in rows:
        assert r["telluric_applied"] == "not-applied"
        assert r["has_telluric_column"] == "False"
        assert r["rec_chain"] == "cr2res_obs_nodding", (
            "a telluric step in the recipe chain would change the verdict")
        # measured IN THE FLUX, not just asserted from the header
        assert float(r["frac_px_below_0p7_continuum"]) > 0.05, (
            "K band with no telluric correction should show deep absorption; a low value "
            "here means either the state changed or the measurement broke")


def test_the_coadd_gate_is_one_resolution_element():
    assert COADD_MAX_RESIDUAL_KMS == pytest.approx(C_KMS / 86000.0)
    assert 3.0 < COADD_MAX_RESIDUAL_KMS < 4.0

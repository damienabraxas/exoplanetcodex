"""
RYA-964 — `resolve_star`: one alias lookup at intake.

The judgements pinned here:

1. IT REFUSES RATHER THAN GUESSES. `STD`, `CAL_*` and `Star S5` are the labels RYA-952 found
   on real frames, and they must land in quarantine, not on a star.
2. ONE STAR, MANY SPELLINGS. RYA-952 found eps Eri's 236 HARPS files split across four
   different `OBJECT` strings. All four resolve to one id.
3. THE ALIAS SET LIVES IN THE CATALOGUE. Nothing hardcodes a name map, so the module must
   read `system_catalog.csv` and a new spelling must be a CSV edit.
4. AN AMBIGUOUS ALIAS IS A CATALOGUE BUG, not a tie to break.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))
from pipeline.star_id import (  # noqa: E402
    CATALOG, UNRESOLVED, StarLabelUnresolved, _alias_index, known_aliases, normalize_label,
    resolve_star)


# ── judgement 1: junk labels quarantine, loudly ──────────────────────────────────────
@pytest.mark.parametrize("junk", [
    "STD",                                   # the tau Ceti frames (RYA-952)
    "CAL_TauCet-RV_AugDec-NoAO_correct",     # their OBS NAME
    "Star S5",                               # the alpha Cen B junk frame (RYA-423)
    "", "   ", "OBJECT", "SKY", "DARK",
])
def test_a_role_or_placeholder_is_never_a_star(junk):
    assert resolve_star(junk) == UNRESOLVED


def test_strict_mode_raises_and_says_what_to_do():
    with pytest.raises(StarLabelUnresolved) as e:
        resolve_star("Star S5", strict=True)
    msg = str(e.value)
    assert "system_catalog.csv" in msg, "the error must name the file a human would edit"
    assert "quarantine" in msg.lower()


# ── judgement 2: one star, many spellings ────────────────────────────────────────────
@pytest.mark.parametrize("label,expected", [
    # the ticket's own smoke test
    ("alf Cen A", "alpha_cen_a"), ("HD128620", "alpha_cen_a"), ("tau_cet", "tau_ceti"),
    # the four eps Eri spellings found in one HARPS pile (RYA-952)
    ("eps Eri", "eps_eri"), ("HD 22049", "eps_eri"), ("Epsilon Eridani", "eps_eri"),
    ("Epsilon-Eridani", "eps_eri"), ("EPSERI", "eps_eri"),
    # separator and case folding
    ("alf_Cen_A", "alpha_cen_a"), ("ALPHA-CEN-A", "alpha_cen_a"), ("alpha cen a", "alpha_cen_a"),
    ("HD 128621", "alpha_cen_b"), ("alf Cen B", "alpha_cen_b"),
    ("rho01 Cnc", "55cnc_a"), ("HD75732", "55cnc_a"), ("55Cnc", "55cnc_a"),
    ("tau Ceti", "tau_ceti"), ("HD 10700", "tau_ceti"), ("HIP8102", "tau_ceti"),
    ("Procyon", "procyon"), ("SUN", "solar"), ("sol", "solar"),
])
def test_known_labels_resolve(label, expected):
    assert resolve_star(label) == expected


def test_hd128620_and_hd128621_do_not_bleed_into_each_other():
    """The old inline test needed an explicit guard for exactly this pair."""
    assert resolve_star("HD128620") == "alpha_cen_a"
    assert resolve_star("HD128621") == "alpha_cen_b"


def test_normalisation_folds_separators_and_simbad_decoration():
    assert normalize_label("* tau Cet") == normalize_label("tau_cet") == "taucet"
    assert normalize_label("ALF-CEN-A") == normalize_label("alf cen a") == "alfcena"


# ── judgement 3: the catalogue is the only source ────────────────────────────────────
def test_the_alias_set_comes_from_the_catalogue_file():
    with open(CATALOG, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert "aliases" in rows[0], "system_catalog.csv must carry an `aliases` column"
    keyed = [r for r in rows if (r.get("star_params_key") or "").strip()]
    assert keyed, "no system carries a star_params_key"
    for r in keyed:
        assert (r.get("aliases") or "").strip(), (
            f"{r['system_name']} has params but no aliases — intake cannot name it")


def test_a_system_with_no_params_is_not_resolvable():
    """tau Boo has CRIRES data but no parameters (RYA-957), so it is not yet a target."""
    assert resolve_star("tau Boo") == UNRESOLVED
    assert resolve_star("HD120136") == UNRESOLVED


def test_every_system_id_resolves_to_itself():
    for sid in known_aliases():
        assert resolve_star(sid) == sid


# ── judgement 4: ambiguity is a bug ──────────────────────────────────────────────────
def test_an_ambiguous_alias_is_refused(tmp_path):
    src = list(csv.DictReader(open(CATALOG, encoding="utf-8")))
    cols = list(src[0].keys())
    a = next(r for r in src if r["star_params_key"] == "alpha_cen_a")
    b = next(r for r in src if r["star_params_key"] == "alpha_cen_b")
    b["aliases"] = (b["aliases"] or "") + "|" + a["aliases"].split("|")[0]
    bad = tmp_path / "system_catalog.csv"
    with bad.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(src)
    _alias_index.cache_clear()
    with pytest.raises(ValueError, match="cannot mean two stars"):
        _alias_index(str(bad))
    _alias_index.cache_clear()


# ── the cited backfill ───────────────────────────────────────────────────────────────
def test_hd_and_hip_are_backfilled_for_the_stars_we_hold():
    with open(CATALOG, encoding="utf-8") as fh:
        rows = {r["star_params_key"]: r for r in csv.DictReader(fh)
                if (r.get("star_params_key") or "").strip()}
    for sid in ("procyon", "alpha_cen_a", "alpha_cen_b", "55cnc_a", "tau_ceti", "eps_eri"):
        assert rows[sid]["hd_id"].strip(), f"{sid} has no hd_id"
        assert rows[sid]["hip_id"].strip(), f"{sid} has no hip_id"
    # The Sun has neither, correctly — it is not in the HD or HIP catalogues.
    assert not rows["solar"]["hd_id"].strip()

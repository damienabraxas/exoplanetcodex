"""RYA-1181 — the intake discarded identity it had parsed, and called held data unpublished.

RYA-1142 A3/A9. Two separate ways of understating what the repo actually holds:

  A3  J'' is parsed for every primary transition and thrown away (it survived only folded
      inside gf = f*(2J''+1)); `system` and the vibrational band reach
      primary_molecular_crossmatch and are dropped in the merge to the artifact everything
      downstream reads. Missing rotational identity is the intake's OWN stated blocker.
  A9  NH A-X, OH A-X and CN B-X are recorded `count=NOT_PUBLISHED`, "individual list not
      published" — while the transitions sit in this repo, unread.

No gf value changes here; every artifact gains columns and moves none.
"""
from __future__ import annotations

import csv
import io
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

AUDIT = ROOT / "data/audit/rya1136_cno_intake"
PHYS = AUDIT / "molecular_physical_crossmatch.csv"
PRIM = AUDIT / "primary_molecular_crossmatch.csv"
LEDGER = AUDIT / "rejected_indicator_ledger.csv"
HELD = AUDIT / "held_uv_transitions_rya1181.csv"
SCOPE = AUDIT / "delivered_scope_rya1181.json"


def _rows(p: Path) -> list[dict]:
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh))


# ── A3: identity carried ─────────────────────────────────────────────────────
def test_system_and_vibrational_band_reach_the_physical_crossmatch():
    """They were published, parsed, and dropped in the merge. `source_band` is the
    REPORTING bin (VIS/NIR/IR) and is not a substitute for the vibrational band."""
    rows = _rows(PHYS)
    assert all(r["system"].strip() for r in rows), "system missing on some row"
    assert all(r["vibrational_band"].strip() for r in rows), "vibrational band missing"
    assert {r["system"] for r in rows} >= {"A-X", "Swan", "X-X"}
    assert {r["vibrational_band"] for r in rows} != {r["source_band"] for r in rows}


def test_j_lower_is_carried_rather_than_folded_into_gf():
    """🔴 The intake's own blocker is missing rotational identity, and J'' was parsed for
    every primary transition and discarded."""
    prim = _rows(PRIM)
    assert "primary_j_lower" in prim[0]
    matched = [r for r in prim if r["join_status"] != "UNMATCHED"]
    assert sum(1 for r in matched if r["primary_j_lower"].strip()) == len(matched)

    rows = _rows(PHYS)
    with_j = [r for r in rows
              if r["primary_j_lower"].strip() or r["j_lower_label"].strip()]
    assert len(with_j) >= 400, f"only {len(with_j)}/{len(rows)} rows carry a J''"


def test_the_co_and_ch_label_remainders_are_split():
    """Both 'transition label' columns were unsplit parser leftovers: the CO Turbospectrum
    record tail (which hides J', J'', v', v'') and CH's branch glued to its o-c residual."""
    rows = _rows(PHYS)
    co = [r for r in rows if r["species"] == "12C16O"]
    ch = [r for r in rows if r["species"] == "CH"]
    assert co and ch
    assert all(r["j_upper"].strip() and r["j_lower_label"].strip()
               and r["v_upper"].strip() and r["v_lower"].strip() for r in co)
    assert all(r["branch"].strip() and r["o_minus_c"].strip() for r in ch)
    # and a field a species does not publish stays EMPTY rather than invented
    assert all(not r["branch"].strip() for r in co)


def test_no_value_moved_only_columns_added():
    """🔴 THE CONSTRAINT. Every artifact gains identity columns and moves no number."""
    for path, keys in ((PHYS, ("published_loggf", "lower_energy_eV",
                               "wavelength_vac_nm", "join_status", "species")),
                       (PRIM, ("published_loggf", "summed_loggf", "primary_loggfs",
                               "join_status", "primary_energies_eV"))):
        before = subprocess.run(["git", "show", f"HEAD:{path.relative_to(ROOT)}"],
                                cwd=ROOT, capture_output=True, text=True)
        if before.returncode != 0 or not before.stdout:
            pytest.skip("pre-change revision unavailable")
        old = list(csv.DictReader(io.StringIO(before.stdout)))
        new = _rows(path)
        if set(new[0]) == set(old[0]):
            continue                      # HEAD already contains the change
        assert len(old) == len(new)
        assert set(old[0]) < set(new[0]), "a column was removed, not added"
        for k in keys:
            assert all(a[k] == b[k] for a, b in zip(old, new)), f"{path.name}: {k} moved"


# ── A9: the held UV ──────────────────────────────────────────────────────────
def test_the_held_uv_lists_are_read_and_staged():
    """The two loose A-X CSVs sit beside the archives the ingest opens; `parse_brooke_xx`
    reads only the X-X member INSIDE each zip, so these were never opened."""
    assert HELD.exists()
    rows = _rows(HELD)
    nh = [r for r in rows if r["species"] == "NH"]
    oh = [r for r in rows if r["species"] == "OH"]
    assert len(nh) > 12000 and len(oh) > 9000
    # the identity the intake says it lacks is present on every staged NH row
    assert all(r["j_lower"].strip() and r["j_upper"].strip() for r in nh)
    assert all(r["branch"].strip() for r in nh)
    assert {r["system"] for r in rows} == {"A-X"}


def test_the_negative_is_restated_precisely():
    """🔴 'individual list not published' reads as NOT AVAILABLE. What is unpublished is
    the SELECTION; the transitions are in the repo."""
    rows = {(r["species"], r["system"]): r for r in _rows(LEDGER)}
    assert len(rows) == 4, "the ledger lost rows — DictWriter takes fieldnames from row 0"
    for key in (("NH", "A-X"), ("OH", "A-X"), ("CN", "B-X")):
        r = rows[key]
        assert r["count"] != "NOT_PUBLISHED" and r["count"].isdigit() and int(r["count"]) > 0
        assert r["use_status"] == "NOT_SELECTED_BY_SOURCE"
        assert "SELECTION is unpublished" in r["reason"]
        held = r["transitions_held_at"]
        assert held and (ROOT / held).exists(), f"{key} claims a path that is not there"


def test_the_parsed_but_unqueried_systems_are_counted():
    """CN B-X is parsed off disk into memory and never looked up, because the index is
    keyed on (species, system) and Amarsi's Table 2 has no B-X row. Not lost — but nothing
    counted it, which is how ~46k held lines coexist with 'not published'."""
    import rya1181_stage_held_uv as S
    unq = S.parsed_but_unqueried()
    assert unq.get("CN B-X", 0) > 40000
    assert sum(unq.values()) > 80000
    doc = json.loads((AUDIT / "held_uv_summary_rya1181.json").read_text())
    assert doc["total_held_uv_not_previously_counted"] > 100000


def test_staging_did_not_rewire_the_match_index():
    """⚠️ Wiring the A-X lists into inventory() would make them candidates for the A-X
    targets and move join_status on live rows. This ticket forbids a value change, so the
    read is a STAGING step and the match set is untouched."""
    import ingest_cno_molecular_primary_rya1136 as M
    src = M.inventory.__doc__ or ""
    body = Path(M.__file__).read_text()
    assert "NH-A-X-linelist.csv" not in body, "the A-X list was wired into the ingest"
    assert "OH-A-X-linelist-final.csv" not in body


# ── scope truth ──────────────────────────────────────────────────────────────
def test_the_delivered_scope_is_stated_and_is_vis_to_ir():
    d = json.loads(SCOPE.read_text())
    assert d["delivered"] == "VIS to IR"
    assert d["fuv_matched"] == 0 and d["nuv_matched"] == 0
    lo, hi = d["atomic_span_A"]
    assert lo > 5000 and hi < 11000, (lo, hi)
    assert "FUV/NUV/IR" in d["claimed_by_ticket_titles"]
    assert "FOLLOW-ON" in d["why_the_uv_bins_are_empty"]


def test_the_empty_uv_bins_are_explained_not_left_bare():
    """An empty bin in a coverage matrix reads as 'looked and found nothing'. Here nothing
    was looked at, and the record now says which."""
    d = json.loads(SCOPE.read_text())
    why = d["why_the_uv_bins_are_empty"]
    assert "Not a null result" in why
    assert "held" in why and "rejected_indicator_ledger" in why

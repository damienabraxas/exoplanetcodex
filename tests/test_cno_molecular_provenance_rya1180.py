"""RYA-1180 — the CNO molecular intake must not overstate what it acquired.

RYA-1142 checks A4/A5 found two provenance defects, both claiming more than the disk holds:

  A4  `Li2015_CO` was labelled a primary "Li et al. 2015, ApJS 216, 15" acquisition. The
      asset's own second line reads `'ExoMol Li2015'`, and it is the SOLE source behind
      all 80 12C16O PHYSICAL_TUPLE_MATCH rows — the intake's entire clean-match class.
  A5  the constants ledger asserted `PRIMARY_TABLES_ACQUIRED` with
      `partition_function_source='Barklem & Collet 2016'` for table6/table7, which are not
      on disk, while recording none of the six De values in the table that IS.

No gf, energy or constant value is changed by this ticket.
"""
from __future__ import annotations

import csv
import hashlib
import statistics
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

AUDIT = ROOT / "data/audit/rya1136_cno_intake"
BIB = AUDIT / "source_bibliography.csv"
LEDGER = AUDIT / "molecular_constants_ledger.csv"
XMATCH = AUDIT / "molecular_physical_crossmatch.csv"
BARKLEM = ROOT / "data/reference/cno_molecular_primary/constants_barklem2016"
CO_DAT = ROOT / "data/linelists/molecular/turbospectrum/CO/CO_IR_Li2015.dat"


def _rows(p: Path) -> list[dict]:
    return list(csv.DictReader(p.open()))


# ── A4 — the CO redistribution ───────────────────────────────────────────────
def test_the_co_asset_says_exomol_on_its_own_second_line():
    """The premise, read from the file rather than taken from the ticket."""
    with CO_DAT.open() as fh:
        fh.readline()
        assert "ExoMol Li2015" in fh.readline()


def test_no_co_primary_table_was_ever_acquired():
    primary = ROOT / "data/reference/cno_molecular_primary"
    assert primary.is_dir()
    assert not [d for d in primary.iterdir() if d.is_dir() and "co" == d.name.split("_")[0]], \
        "a CO primary directory now exists — re-run the join and re-grade the 80 matches"


def test_li2015_co_reads_as_a_redistribution():
    row = next(r for r in _rows(BIB) if r["source_id"] == "Li2015_CO")
    assert row["source_type"] == "redistribution"
    assert row["status"] == "ACQUIRED_AS_REDISTRIBUTION"
    note = row["provenance_note"]
    for phrase in ("NOT a primary", "ExoMol", "RYA-236", "EXTERNAL to this repo",
                   "no CO directory"):
        assert phrase in note, phrase


def test_every_bibliography_row_is_typed():
    rows = _rows(BIB)
    assert rows and all(r["source_type"] for r in rows), "an untyped source is unclassified"
    assert {r["source_type"] for r in rows} <= {
        "article", "primary_published_table", "derived_grid", "redistribution"}
    # exactly one redistribution today; a second must be a deliberate decision
    assert [r["source_id"] for r in rows if r["source_type"] == "redistribution"] == ["Li2015_CO"]


def test_all_eighty_co_matches_carry_the_redistribution_grade():
    """🔴 The consequence, on the rows themselves. The intake's only clean-match class
    rests entirely on a twice-derived list, and the row now says so."""
    rows = _rows(XMATCH)
    co = [r for r in rows if r["species"] == "12C16O"]
    assert len(co) == 80
    assert {r["join_status"] for r in co} == {"PHYSICAL_TUPLE_MATCH"}
    assert {r["source_provenance_grade"] for r in co} == {
        "REDISTRIBUTION_VENDORED_SYNTHESIS_LIST"}
    # and the grade is not applied indiscriminately: the primary-matched rows differ
    other = {r["source_provenance_grade"] for r in rows
             if r["species"] != "12C16O" and r["matched_file"]}
    assert other == {"PRIMARY_PUBLISHED_TABLE"}, other


def test_an_unrecognised_source_is_unclassified_not_primary():
    """RYA-1072: an allow-list for the true case must not launder the case it does not
    recognise. A new vendor path must never default to 'primary'."""
    from build_cno_intake_rya1136 import provenance_grade
    assert provenance_grade("") == ""
    assert provenance_grade("data/somewhere/new_vendor/foo.dat") == "UNCLASSIFIED_SOURCE"
    assert provenance_grade("data/reference/cno_molecular_primary/x") == "PRIMARY_PUBLISHED_TABLE"


# ── A5 — the constants ledger ────────────────────────────────────────────────
#: Adopted De from Barklem & Collet 2016 table1.dat, bytes 84-93.
EXPECTED_DE = {"C2": "6.371000", "CH": "3.469600", "CN": "7.737000",
               "CO": "11.117000", "NH": "3.419000", "OH": "4.417100"}


def test_the_ledger_records_the_six_real_de_values_with_table_id_and_checksum():
    rows = {r["molecule"]: r for r in _rows(LEDGER)}
    assert set(rows) == set(EXPECTED_DE)
    sha = hashlib.sha256((BARKLEM / "table1.dat").read_bytes()).hexdigest()
    for mol, want in EXPECTED_DE.items():
        r = rows[mol]
        assert r["dissociation_energy_eV"] == want, mol
        assert "table1.dat" in r["dissociation_energy_source"]
        assert "84-93" in r["dissociation_energy_source"], "the adopted column must be named"
        assert r["de_source_sha256"] == sha, mol
        assert r["de_source_row_label"] == mol


def test_the_ledger_no_longer_claims_tables_that_are_not_on_disk():
    for absent in ("table6.dat", "table7.dat", "table2"):
        assert not (BARKLEM / absent).exists(), f"{absent} now exists — re-state the ledger"
    for r in _rows(LEDGER):
        assert r["partition_function_status"] == "NOT_ACQUIRED"
        assert r["equilibrium_constant_status"].startswith("NOT_ACQUIRED")
        assert "PRIMARY_TABLES_ACQUIRED" not in r["verdict"]


def test_the_ledger_is_no_longer_byte_identical_boilerplate():
    """Six identical rows carry no per-molecule information, which is how a ledger can be
    wrong about all six at once."""
    rows = _rows(LEDGER)
    assert len({r["dissociation_energy_eV"] for r in rows}) == 6


def test_the_adopted_column_is_read_positionally_not_by_order():
    """🔴 RYA-853's defect in advance. table1.dat has four De columns and they DISAGREE:
    NH is HH 3.470, Luo 3.470, G2 3.378 and ADOPTED 3.419 — a value equal to none of its
    inputs. Whitespace order would have taken a comparison column."""
    from build_cno_intake_rya1136 import barklem_adopted_de
    got = barklem_adopted_de(BARKLEM / "table1.dat")
    assert got["NH"][0] == "3.419000"
    nh = next(l for l in (BARKLEM / "table1.dat").read_text().splitlines()
              if l[0:5].strip() == "NH")
    assert nh.split()[3] == "3.470000", "the first numeric column is NOT the adopted value"
    assert got["NH"][0] != nh.split()[3]


def test_the_bibliography_row_no_longer_names_the_absent_tables_as_its_role():
    row = next(r for r in _rows(BIB) if r["source_id"] == "BarklemCollet2016")
    assert "table1.dat" in row["asset"]
    assert "DISSOCIATION ENERGIES" in row["role"]
    assert "NOT acquired" in row["role"]
    assert row["status"] == "ACQUIRED_PARTIAL"


# ── A5b — the C2 origin constant ─────────────────────────────────────────────
def test_the_c2_origin_offset_is_derived_without_the_circular_energy_cut():
    """🔴 The constant is used INSIDE the join's energy filter, so deriving it from the
    matched rows would be circular. Candidates here are selected on wavenumber and
    vibrational band ONLY."""
    import bisect
    import ingest_cno_molecular_primary_rya1136 as M

    src = sorted(M.parse_c2(), key=lambda t: 1e8 / t.wavelength_vac_A)
    wn = [1e8 / t.wavelength_vac_A for t in src]
    offsets = []
    for target in csv.DictReader(M.AMARSI.open()):
        if target["species"] != "C2":
            continue
        twn = 1e8 / (float(target["wavelength_vac_nm"]) * 10)
        vp, vl = (int(x) for x in target["band"].strip("()").split("-"))
        lo, hi = bisect.bisect_left(wn, twn - 2.0), bisect.bisect_right(wn, twn + 2.0)
        near = [x for x in src[lo:hi] if x.vp == vp and x.vl == vl]
        if not near:
            continue
        best = min(near, key=lambda t: abs(1e8 / t.wavelength_vac_A - twn))
        offsets.append(float(target["lower_energy_eV"])
                       - (best.lower_energy_eV - M.C2_LOWER_ORIGIN_EV))

    assert len(offsets) == 39
    median = statistics.median(offsets)
    cluster = [o for o in offsets if 0.075 <= o <= 0.077]
    assert len(cluster) >= 34, f"the offset cluster dispersed: {len(cluster)}/39"
    assert abs(median - M.C2_LOWER_ORIGIN_EV) < 0.001, (
        f"derived median {median:.6f} no longer reproduces the constant "
        f"{M.C2_LOWER_ORIGIN_EV}; it must be re-derived or cited, not left standing")


def test_the_c2_constant_declares_that_it_is_derived_and_not_cited():
    import ingest_cno_molecular_primary_rya1136 as M
    p = M.C2_LOWER_ORIGIN_PROVENANCE
    assert p["status"] == "DERIVED_NOT_CITED"
    assert p["rows_used"] == 39 and p["rows_in_cluster"] >= 34
    assert "NO energy cut" in p["selection"]
    assert abs(p["derived_median_eV"] - M.C2_LOWER_ORIGIN_EV) < 0.001

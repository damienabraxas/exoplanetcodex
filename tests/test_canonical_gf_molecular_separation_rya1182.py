"""RYA-1182 — molecular transition provenance must not live in the atomic store.

RYA-1130 exists to keep it out. It was not being kept out: RYA-1142 check B2 found 8,977
molecular rows inside `canonical_gf.csv`, every one seeded `linelist(VALD3)` — the exact
redistribution the RYA-1136 CNO intake spent five primary archives avoiding, one join from
being picked up by the next CNO measurement that reached for canonical_gf.

`test_a_reinserted_molecular_row_is_caught` is the smoke test's deliberate reinsertion.
"""
from __future__ import annotations

import collections
import csv
import io
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline import canonical_gf_invariants as INV  # noqa: E402

QUARANTINE = INV.QUARANTINE

#: What RYA-1182 moved, enumerated independently of the ticket's own count.
EXPECTED = {"CN": 2766, "C2": 2654, "CH": 2265, "MgH": 594, "SiH": 583, "NH": 114, "OH": 1}
TOTAL = 8977


@pytest.fixture(scope="module")
def rows():
    return INV.read_rows()


# ── the invariant ────────────────────────────────────────────────────────────
def test_the_atomic_store_holds_only_atomic_transitions(rows):
    bad = INV.non_atomic(rows)
    assert not bad, (
        "canonical_gf.csv holds non-atomic rows: "
        f"{dict(collections.Counter(r['species'].strip() for r in bad))}")
    assert INV.check() == 0


def test_a_reinserted_molecular_row_is_caught():
    """🔴 The smoke test's deliberate reinsertion. A guard that has only ever seen a clean
    store has not been shown to fail."""
    row = {"key_z": "CH", "ion": "", "species": "CH", "seed_source": "linelist(VALD3)",
           "loggf_reference": "VALD3", "gf_tier": "VALD3"}
    assert not INV.is_atomic(row)
    assert INV.non_atomic([row]) == [row]
    assert INV.violations([row]) == [row]


def test_the_rule_is_positive_not_a_molecule_blacklist():
    """A recognised set only catches what somebody already thought of. TiO and FeH are in
    no list here and must still be refused — that is the whole design."""
    for sp in ("TiO", "FeH", "CaH", "ZrO", "H2O", "some_new_molecule"):
        assert not INV.is_atomic({"key_z": sp, "ion": "", "species": sp}), sp
    # and the contract accepts ordinary atomic rows
    assert INV.is_atomic({"key_z": "26", "ion": "I", "species": "Fe I"})
    assert INV.is_atomic({"key_z": "8", "ion": "II", "species": "O II"})
    # a numeric key_z is not enough on its own
    assert not INV.is_atomic({"key_z": "26", "ion": "", "species": "Fe"})
    assert not INV.is_atomic({"key_z": "0", "ion": "I", "species": "?"})
    assert not INV.is_atomic({"key_z": "93", "ion": "I", "species": "Np I"})


def test_a_vetted_molecular_row_is_distinguishable_from_redistribution():
    """The spec's escape hatch, and its narrowness. Being molecular is enough to fail the
    strict contract; `violations()` additionally requires the row to be REDISTRIBUTION, so
    a vetted measurement can be told apart from a VALD3 grab."""
    vetted = {"key_z": "CH", "ion": "", "species": "CH",
              "seed_source": "Masseron2014 MoLLIST", "loggf_reference": "Masseron2014",
              "gf_tier": "LAB"}
    assert INV.violations([vetted]) == [], "a vetted molecular row is not redistribution"
    assert INV.non_atomic([vetted]) == [vetted], "it is still not an atomic transition"


# ── the relocation ───────────────────────────────────────────────────────────
def test_every_relocated_row_is_present_and_unaltered():
    """RYA-946 keep-not-delete: the rows are moved, not dropped, and a reader can still
    find them. Byte-identical, under the same header, in the original order."""
    assert QUARANTINE.exists(), "the quarantine file is missing — rows were deleted"
    text = QUARANTINE.read_text()
    kept = list(csv.DictReader(io.StringIO(text)))
    assert len(kept) == TOTAL
    assert dict(collections.Counter(r["species"].strip() for r in kept)) == EXPECTED

    # same 26-column header as the store they came from
    assert next(csv.reader(io.StringIO(text))) == \
        next(csv.reader(io.StringIO(INV.CANONICAL_GF.read_text())))

    # provenance preserved: all of it is the VALD3 redistribution, unedited
    assert {r["seed_source"] for r in kept} == {"linelist(VALD3)"}
    assert {r["loggf_reference"] for r in kept} == {"VALD3"}
    assert {r["gf_tier"] for r in kept} == {"VALD3"}


def test_the_split_reproduces_the_pre_relocation_store_exactly():
    """🔴 THE VALUE CLAIM. Recombining the atomic store with the quarantine must give back
    the committed pre-RYA-1182 file, line for line — that is what proves no atomic row was
    rewritten and no gf value moved. Skipped where git history is unavailable."""
    before = subprocess.run(
        ["git", "show", f"HEAD:{INV.CANONICAL_GF.relative_to(ROOT)}"],
        cwd=ROOT, capture_output=True, text=True)
    if before.returncode != 0 or not before.stdout:
        pytest.skip("pre-relocation revision not reachable in this checkout")

    old = before.stdout.rstrip("\n").split("\n")
    new = INV.CANONICAL_GF.read_text().rstrip("\n").split("\n")
    quar = QUARANTINE.read_text().rstrip("\n").split("\n")
    if len(old) == len(new):
        pytest.skip("HEAD already contains the relocation")

    assert old[0] == new[0] == quar[0]
    cols = next(csv.reader(io.StringIO(old[0])))
    iz = cols.index("key_z")

    def numeric(v):
        try:
            int(v.strip())
            return True
        except ValueError:
            return False

    exp_atomic, exp_mol = [], []
    for line in old[1:]:
        r = next(csv.reader(io.StringIO(line)))
        (exp_atomic if numeric(r[iz]) else exp_mol).append(line)
    assert exp_atomic == new[1:], "an atomic row was rewritten or reordered"
    assert exp_mol == quar[1:], "a relocated row was rewritten or reordered"


def test_the_quarantine_says_it_is_not_adopted():
    """A file full of unvetted gf values with no label is a trap for the next reader."""
    readme = QUARANTINE.parent / "README.md"
    assert readme.exists()
    text = readme.read_text()
    for phrase in ("NOT ADOPTED", "RYA-1130", "keep-not-delete", "VALD3"):
        assert phrase in text, phrase

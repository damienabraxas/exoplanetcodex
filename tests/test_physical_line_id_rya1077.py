"""RYA-1077 — a canonical_gf id must survive a rebuild of canonical_gf.

🔴 `line_id` did not. It is `gf_000000`, `gf_000001`, … assigned by ROW POSITION, so a
block removal or replacement shifts every id after it. Measured on the committed
artifacts: **1,739 of 6,979 rows (25%) cite an id that now points at a different line** —
`gf_177842` recorded against 21553.299 A now resolves to 21816.566 A, and `gf_177933`
recorded as Fe I is now Ni I.

Appending is NOT the culprit and the tests should not imply it is: the extenders append a
sorted block and leave earlier ids alone. The damage came from a block being replaced
(RYA-1052 found RYA-1047's lab lines unselectable; RYA-1053 re-extended), which is exactly
why rya1060 — generated afterwards — is 0.1% rotted while rya1059 is 66%.
"""
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"

from pipeline.physical_line_id import (  # noqa: E402
    PhysicalKeyError, add_physical_ids, physical_id, physical_key)


def test_the_id_is_a_function_of_the_physics_only():
    """Same transition, same id -- that is the entire contract."""
    a = physical_id("Fe I", 5250.2084, 0.121)
    b = physical_id(" Fe I ", 5250.2084, 0.121)     # whitespace is not physics
    assert a == b == physical_id("Fe I", float("5250.2084"), 0.121)


def test_lines_differing_in_the_FOURTH_decimal_get_DIFFERENT_ids():
    """🔴 THE ROUNDING TRAP, REFUSED BY CONSTRUCTION. canonical_gf's wavelengths carry a
    median of 3 decimals but 13.3% carry four or more, so ANY fixed precision fuses real
    lines. RYA-1033 paid for this once: a rounded number is not an identity."""
    assert physical_id("Fe I", 5250.2084, 0.121) != physical_id("Fe I", 5250.2085, 0.121)
    assert physical_id("Fe I", 5250.2084, 0.1210) == physical_id("Fe I", 5250.2084, 0.121)


def test_the_key_carries_no_rounding_at_all():
    """`repr` is the shortest string that round-trips a float64 exactly, so the precision
    is DERIVED from the value. A `format`-based key would invent precision on some rows and
    discard it on others, in the same file."""
    k = physical_key("Fe I", 5250.2084, 0.121)
    assert "5250.2084" in k and "0.121" in k
    assert "%" not in k and ":.3f" not in k


def test_species_is_part_of_the_identity():
    """gf_177933 was recorded as Fe I and now resolves to Ni I. Wavelength alone is not an
    identity either."""
    assert physical_id("Fe I", 5250.2084, 0.121) != physical_id("Ni I", 5250.2084, 0.121)


def test_a_row_with_no_physical_key_is_REFUSED():
    """A row that cannot be identified must not be given an id that looks like one."""
    with pytest.raises(PhysicalKeyError):
        physical_id("Fe I", float("nan"), 0.121)
    with pytest.raises(PhysicalKeyError):
        physical_id("Fe I", "not-a-number", 0.121)


def test_a_COLLISION_is_refused_never_tiebroken():
    """Two rows sharing a key are a duplicate or a key too coarse -- both want a human.
    Keeping the first silently is how a join starts lying (RYA-1033)."""
    df = pd.DataFrame([{"species": "Fe I", "wavelength_air_A": 5250.2084,
                        "excitation_potential_eV": 0.121, "log_gf": -4.9},
                       {"species": "Fe I", "wavelength_air_A": 5250.2084,
                        "excitation_potential_eV": 0.121, "log_gf": -4.8}])
    with pytest.raises(PhysicalKeyError) as e:
        add_physical_ids(df)
    assert "REFUSING" in str(e.value)


def test_the_live_table_has_a_stable_id_on_every_row_and_no_collisions():
    """Asserted on the ARTIFACT, not a fixture -- a fixture cannot rot."""
    cg = pd.read_csv(CANON, low_memory=False)
    assert "physical_id" in cg.columns, "canonical_gf must carry the stable id"
    assert cg.physical_id.notna().all()
    assert cg.physical_id.nunique() == len(cg), (
        f"{len(cg) - cg.physical_id.nunique()} colliding physical_id(s)")
    assert cg.physical_id.str.startswith("pk_").all()


def test_the_stable_id_REGENERATES_identically_from_the_physics():
    """The property `line_id` lacked: rebuild it and it must come out the same."""
    cg = pd.read_csv(CANON, low_memory=False).head(2000)
    again = [physical_id(r.species, r.wavelength_air_A, r.excitation_potential_eV)
             for _, r in cg.iterrows()]
    assert again == list(cg.physical_id)


def test_line_id_is_KEPT_so_nothing_breaks_on_this_commit():
    """The positional id is not deleted here. 18 artifacts still cite it, and removing the
    column would turn a wrong reference into a missing one -- which is worse, because a
    missing reference cannot be re-resolved by physical key later."""
    cg = pd.read_csv(CANON, low_memory=False, nrows=5)
    assert "line_id" in cg.columns

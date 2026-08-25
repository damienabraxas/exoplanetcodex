"""RYA-1035: which published ⟨3D⟩ aux tables carry the metallicity-zeroing defect.

RYA-1035 found the defect in Fe and decided, for Fe alone, which of the two shipped aux
files to register against. That decision is needed once per element — and discovering it
one element at a time is the pattern this project keeps paying for (RYA-710's own standing
note calls it the third occurrence of "check the source, not the extract"; the Fe deck was
the fourth). So it was swept across all seventeen at once.

These tests pin the SWEEP'S CONCLUSIONS, not the fetch: the artifact is committed and the
generator re-derives it live. What must not drift is the decision each element's row
implies, because that decision is a registry line somebody will write later.
"""
import csv
import json
from pathlib import Path

import pytest

RES = Path(__file__).resolve().parents[1] / "data" / "results" / "rya1035"
JSON = RES / "mean3d_aux_defect_sweep.json"
CSV = RES / "mean3d_aux_defect_sweep.csv"


@pytest.fixture(scope="module")
def sweep():
    return json.loads(JSON.read_text())


@pytest.fixture(scope="module")
def rows():
    with CSV.open() as fh:
        return {r["element"]: r for r in csv.DictReader(fh)}


def test_all_seventeen_published_decks_were_scored(sweep):
    """The sweep is only worth having if it is exhaustive — a partial one recreates the
    element-at-a-time discovery it exists to replace."""
    assert len(sweep["elements"]) == 17
    assert not [el for el, d in sweep["elements"].items() if "error" in d]


def test_exactly_Fe_and_Mn_are_defective(sweep):
    """🔴 TWO OF SEVENTEEN. Both are May-2021 vintage (Fe May-21, Mn May-17), and both put
    the defect in exactly the same place: the seven Teff=5777 STAGGER solar members."""
    defective = sorted(el for el, d in sweep["elements"].items()
                       if d["plain"]["n_feh_overridden"])
    assert defective == ["Fe", "Mn"]


def test_the_defect_is_confined_to_the_stagger_solar_member(sweep):
    """Every overridden row sits at (5777, 4.44) on both decks. That is what makes the
    model NAME a usable referee: the other thousands of rows agree with their names
    exactly, so the column is wrong and the name is right — not the other way round."""
    for el in ("Fe", "Mn"):
        nodes = sweep["elements"][el]["plain"]["overridden_teff_logg"]
        assert nodes == [[5777.0, 4.44]], el


def test_the_row_count_scales_with_the_decks_abundance_axis(sweep):
    """🔴 THE SAME DEFECT LOOKS DIFFERENT SIZED ON DIFFERENT DECKS. Six metallicities are
    wrong on each; Fe resolves ONE abundance per node so that is 6 rows, Mn resolves 25 so
    it is 150. Counting rows rather than nodes would make Mn look like a worse problem and
    Fe like a rounding error."""
    assert sweep["elements"]["Fe"]["plain"]["n_feh_overridden"] == 6
    assert sweep["elements"]["Mn"]["plain"]["n_feh_overridden"] == 150
    assert sweep["elements"]["Mn"]["plain"]["solar_node_distinct_abundances"] == 25
    assert sweep["elements"]["Fe"]["plain"]["solar_node_distinct_abundances"] == 1


def test_the_conversion_is_unrecoverable_on_both_defective_decks(sweep):
    """🔴 `convert_3d_grid_to_marcs_names.py` BUILDS THE NAME FROM THE [Fe/H] COLUMN, so on
    a defective deck it propagates the zeroing into the name and seven distinct atmospheres
    become one byte-identical string. After that nothing can referee it: name and column
    agree, and both are wrong."""
    for el in ("Fe", "Mn"):
        conv = sweep["elements"][el]["marcs_names"]
        assert conv["n_distinct_names_at_5777"] == 1, el
        assert conv["n_rows_at_5777"] > 1, el
        assert conv["n_feh_overridden"] == 0, f"{el}: nothing left to referee with"


def test_Al_is_the_positive_control(sweep):
    """The conversion is faithful when the column it reads is right — 217 rows at Teff=5777
    keep all SEVEN names. Without this the collapse above could be blamed on the converter
    in general rather than on the defective input it was given."""
    al = sweep["elements"]["Al"]
    assert al["plain"]["n_feh_overridden"] == 0
    assert al["marcs_names"]["n_distinct_names_at_5777"] == 7
    assert al["plain"]["n_distinct_names_at_5777"] == 7


def test_every_element_has_an_addressable_solar_node_after_refereeing(sweep):
    """The referee's whole job: one unambiguous row at the solar node, for all seventeen.
    Before it, the defective decks tie seven ways and file order picks [Fe/H] = −1.0."""
    for el, d in sweep["elements"].items():
        assert d["plain"]["solar_node_distinct_ids"] == 1, el


def test_no_deck_carries_an_unparseable_model_name(sweep):
    """`read_deck_node` refuses a record it cannot identify, so an unparseable name is an
    unusable node. This is what caught the short-Teff form (`p50g25m40`) on Fe, which is
    182 of its 189 rows."""
    for el, d in sweep["elements"].items():
        assert d["plain"]["unparseable_names"] == 0, el


def test_the_csv_carries_the_registry_decision_for_each_element(rows):
    """The row a future wiring ticket actually reads. `plain` is a REQUIREMENT, `either`
    is permission — Al is registered against `_marcs_names` and must stay valid."""
    assert rows["Fe"]["register_against"] == "plain"
    assert rows["Mn"]["register_against"] == "plain"
    assert rows["Al"]["register_against"] == "either"
    assert rows["Cr"]["register_against"] == "either"
    assert rows["Eu"]["register_against"] == "either"
    assert rows["Y"]["register_against"] == "either"
    assert rows["Fe"]["marcs_names_collapsed"] == "True"
    assert rows["Al"]["marcs_names_collapsed"] == "False"


def test_the_registered_Fe_deck_obeys_the_sweeps_own_verdict(rows):
    """The sweep is only useful if the registry follows it. Fe is the one ⟨3D⟩ deck wired
    so far whose verdict is `plain` — RYA-710 registered it that way."""
    from pipeline import gerber_nlte as G
    assert rows["Fe"]["register_against"] == "plain"
    assert "_marcs_names" not in G.DECKS["Fe@mean3D"]["aux"]
    # ...and Al's `either` verdict is why its `_marcs_names` registration stays valid
    assert rows["Al"]["register_against"] == "either"
    assert "_marcs_names" in G.DECKS["Al@mean3D"]["aux"]

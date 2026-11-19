"""RYA-470: a published-set line is graded from ITS OWN canonical row, told apart from a second
transition within 0.01 A by the lower level -- never left ungraded by a neighbour, never guessed."""
import pandas as pd

from pipeline import line_sets as ls
from pipeline.gf_grades import is_gf_graded


def test_every_set_line_with_one_level_matched_row_takes_that_rows_grade():
    canon = pd.read_csv(ls.CANONICAL, low_memory=False)
    g = is_gf_graded(canon)
    reg = pd.read_csv(ls.REGISTRY)
    for r in reg[reg.element == "V"].itertuples():
        t = pd.read_csv(ls.ROOT / r.csv)
        for x in t.itertuples():
            m = ((canon.species == x.species) & ((canon.wavelength_air_A - x.wavelength_air_A).abs() < 0.01)
                 & ((canon.excitation_potential_eV - x.ep_eV).abs() <= 0.02))
            if m.sum() == 1:
                assert bool(x.gf_graded) == bool(g[m].iloc[0]), (r.set_name, x.wavelength_air_A)


def test_a_neighbour_within_001_A_no_longer_ungrades_the_asplund_line():
    t = pd.read_csv(ls.ROOT / "data/reference/line_sets/v_agss21_VI.csv")
    row = t[(t.wavelength_air_A - 5727.6521).abs() < 1e-4]
    assert len(row) == 1 and bool(row.gf_graded.iloc[0])

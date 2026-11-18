"""RYA-470: a gf adopted into canonical_gf with its source's own per-line sigma is GRADED on
that sigma for EVERY species -- not only the Fe/Al/Si species that have a LAB_TABLES entry.
The first V I VIS Codex pool (Lawler 2014 + NIST classes, every line carrying a published
sigma) was graded '152 x K07' and stamped NOT PUBLISHABLE before this."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import gf_grades as gg, gf_rung  # noqa: E402


def _rows(species, tier_has):
    c = gg.canonical_species(species)
    sig = pd.to_numeric(c.gf_sigma_dex, errors="coerce")
    return c[c.gf_tier.astype(str).str.contains(tier_has) & sig.notna()]


def test_an_adopted_lab_row_grades_gf_lab_on_its_own_sigma():
    c = _rows("V I", "LAB")
    assert len(c) > 100 and "V I" not in gg.LAB_TABLES      # graded WITHOUT a lab table
    r = c.iloc[0]
    v = gg.grade_line(r.wavelength_air_A, r.excitation_potential_eV, r.log_gf, species="V I")
    assert v.gf_grade == gg.GRADE_LAB and v.gf_sigma_dex == float(r.gf_sigma_dex)


def test_the_sigma_describes_only_the_value_it_was_published_with():
    r = _rows("V I", "LAB").iloc[0]
    v = gg.grade_line(r.wavelength_air_A, r.excitation_potential_eV, r.log_gf + 0.3,
                      species="V I")
    assert v.gf_grade == gg.GRADE_MISMATCH and not v.has_cited_sigma


def test_a_lab_plus_nist_pool_grades_rung_3_and_one_ungraded_line_still_binds_rung_1():
    lab, nist = _rows("V I", "LAB").head(8), _rows("V I", "NIST").head(3)
    pool = pd.DataFrame({"wavelength_air_A": list(lab.wavelength_air_A) + list(nist.wavelength_air_A),
                         "ep_eV": list(lab.excitation_potential_eV) + list(nist.excitation_potential_eV),
                         "log_gf": list(lab.log_gf) + list(nist.log_gf)})
    r = gf_rung.decide("V", "I", pool)
    assert r.rung == 3 and r.n_graded == len(pool) and np.isfinite(r.cited_sigma_dex)
    k = gg.canonical_species("V I")
    k = k[k.gf_tier.astype(str).eq("KURUCZ")].iloc[0]
    mixed = pd.concat([pool, pd.DataFrame([{"wavelength_air_A": k.wavelength_air_A,
                                            "ep_eV": k.excitation_potential_eV,
                                            "log_gf": k.log_gf}])], ignore_index=True)
    assert gf_rung.decide("V", "I", mixed).rung == 1

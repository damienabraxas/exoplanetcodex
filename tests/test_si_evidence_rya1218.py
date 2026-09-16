import pandas as pd
import pytest

from pipeline import si_evidence as si
from pipeline import uncertainty_contract as uc


def lines(values=(7.0, 7.1, 8.0)):
    return pd.DataFrame([dict(element="Si", ion="I", wavelength_air_A=5000+i,
                              ep_eV=4.0, abundance=v, in_aggregate=True,
                              line_id=f"si:{i}") for i, v in enumerate(values)])


def test_estimator_mismatch_is_exposed_without_changing_the_abundance():
    result = si.diagnostic_statistics(lines(), {"n_lines": 3, "A": 7.1, "stat_dex": .3})
    assert result["median"] == 7.1
    assert result["mean"] != result["median"]
    assert result["statistical_admission"].startswith("HOLD")


def test_rejected_line_cannot_enter_statistics_or_pairing():
    a = lines(); a.loc[2, "in_aggregate"] = False
    result = si.paired_difference(a, lines((7.1, 7.2, 99)))
    assert result["n"] == 2
    assert result["mean_delta"] == pytest.approx(-.1)
    assert result["right_only"] == ["si:2"]
    assert result["sigma_differential"] is None


def test_identity_cannot_cross_ion_or_excitation():
    a = lines((7.0,))
    for species, ep in [("Si II", 4.), ("Si I", 8.)]:
        c = pd.DataFrame([dict(species=species, wavelength_air_A=5000.,
                              excitation_potential_eV=ep, line_id="foreign",
                              log_gf=-1., loggf_reference="synthetic")])
        with pytest.raises(si.line_match.LineMatchError):
            si.attach_identity(a, c)


def test_duplicate_physical_lines_are_not_independent_evidence():
    a = lines(); a.loc[1, "line_id"] = a.loc[0, "line_id"]
    with pytest.raises(ValueError, match="repeats"):
        si.accepted(a)


def test_legacy_budget_maps_to_canonical_holds_not_false_zeroes():
    p = dict(element="Si", ion="I", instrument="harps", holding="solar_harps_molecfit_corrected",
             band="VIS", tier="ALL", selector="FROMEW", route="SYNTH", treatment="1D-LTE",
             line_set="si-agss21")
    b = si.hold_budget(p, ["si:0"], source="synthetic archived budget", statistics={})
    assert b["schema"] == uc.SCHEMA
    assert set(b["holds"]) == set(uc.COMPONENTS)
    assert b["sigma_reported"] is None
    assert b["sigma_measured_partial"] is None
    with pytest.raises(uc.UncertaintyError, match="unresolved"):
        uc.validate(b, scope=b["scope"])


def test_si_ii_single_line_results_are_explicitly_diagnostic():
    expected = {"harps": 7.562, "iag": 7.618, "kpno": 7.578}
    root = __import__("pathlib").Path(__file__).parents[1] / "data/results/rya1218/recovered_20260916/si_ii_synth"
    for holding, abundance in expected.items():
        row = __import__("pandas").read_csv(next((root / holding).glob("*products.csv"))).iloc[0]
        assert row.A == pytest.approx(abundance)
        assert row.n_lines == 1
        assert "UNMEASURED" in row.stat_basis
        assert row.syst_dex == pytest.approx(0.17)

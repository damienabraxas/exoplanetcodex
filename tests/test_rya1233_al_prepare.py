"""RYA-1233: Al through prepare --apply -- the pieces Si needed by hand, now orchestrator steps.

  - a published solar-literature line missing everywhere (Al I 10768.363) is ADDED: its VALD
    components into the band's synthesis list, one canonical row at the authors' gf
  - an adopted gf carries the authors' reference, not the replaced value's
  - every gf source of Al's graded lines is classified in the correlation registry
    (NIST-graded lines without a lab tag map to NIST_ASD)
  - a cell's identity hashes canonical_gf and the synthesis list ONLY in its own window
"""
import pandas as pd
import pytest

from pipeline import gf_error_model, published_line_add as pla, run_matrix as rm
from pipeline.element_prepare import adopted_reference

ROOT = rm.ROOT


def test_source_key_maps_untagged_nist_rows():
    assert gf_error_model.source_key("Burheim2023", "x", "") == "Burheim2023"
    assert gf_error_model.source_key("", "10.1063/1.2734566", "C+") == gf_error_model.NIST_KEY
    assert gf_error_model.source_key("", "NIST ASD v5.11; Kelleher", "") == gf_error_model.NIST_KEY
    assert gf_error_model.source_key("", "K75", "") == ""


def test_al_gf_sources_are_classified():
    model = gf_error_model.load()
    for tag in ("AL_AGSS21", "Burheim2023", gf_error_model.NIST_KEY):
        assert tag in model, tag
    # TOPbase is one calculation measured 0.11 dex off experiment on both 5p lines: correlated
    assert model["AL_AGSS21"]["correlation"] == "fully_correlated"
    # the budget passes the KEYED tag (source_key with the row's NIST grade)
    tags = [gf_error_model.source_key("", "10.1063/1.2734566", "C+"),
            gf_error_model.source_key("Burheim2023", "x", "")]
    _, note, unreviewed = gf_error_model.covariance([0.06, 0.05], ["10.1063/1.2734566", "x"], tags)
    assert unreviewed == [] and "UNREVIEWED" not in note


def test_adopted_reference_names_the_set_and_the_per_line_source():
    cite = {"AL_AGSS21": {"citation": "Nordlander & Lind 2017"}}
    assert adopted_reference({"set": "AL_AGSS21", "authors_gf_source": "TOPbase"}, cite) == \
        "Nordlander & Lind 2017 (gf: TOPbase)"
    assert adopted_reference({"set": "AL_AGSS21"}, cite) == "Nordlander & Lind 2017"


def test_al_10768_is_in_canonical_at_the_authors_gf():
    c = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    r = c[(c.species == "Al I") & ((c.wavelength_air_A - 10768.363).abs() < 0.005)]
    assert len(r) == 1
    r = r.iloc[0]
    assert (r.log_gf, r.gf_sigma_dex, r.gf_tier) == (-2.02, 0.06, "PUBLISHED-SET")
    assert "Nordlander & Lind 2017" in r.loggf_reference and "TOPbase" in r.loggf_reference
    t = pd.read_csv(ROOT / "data/linelists/ispec_ir_9200_13000/atomic_lines.tsv", sep="\t",
                    usecols=["element", "wave_A"])
    assert ((t.element == "Al 1") & ((t.wave_A - 10768.36).abs() < 0.05)).sum() == 6


def test_al_adoptions_cite_the_authors():
    c = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    a = c[c.adjudication_status.astype(str).str.endswith("rya1233_al_agss21")]
    assert len(a) == 6
    assert a.loggf_reference.str.startswith("Asplund 2021 -> Nordlander & Lind 2017").all()


def _comps(w, gfs, eup=5.2363):
    return [{"element": "Al", "ion": "I", "wavelength": w + 0.004 * i, "log_gf": g,
             "e_low_eV": 4.0853, "e_up_eV": eup} for i, g in enumerate(gfs)]


def test_plan_takes_the_band_list_and_refuses_on_a_clash(monkeypatch, tmp_path):
    monkeypatch.setattr(pla, "components", lambda sp, w, ep: (_comps(w, [-2.4, -2.3, -1.9]),
                                                               tmp_path / "vald_x.txt"))
    monkeypatch.setattr(pla, "clashes", lambda el, waves: [])
    p = pla.plan({"species": "Al I", "wavelength_A": 10768.363, "ep_eV": 4.085, "set": "S",
                  "authors_loggf": -2.02, "authors_sigma_dex": 0.06})
    assert p["synthesis_list"].endswith("ispec_ir_9200_13000/atomic_lines.tsv")
    assert p["n_components"] == 3 and p["clashes"] == []
    monkeypatch.setattr(pla, "clashes", lambda el, waves: ["SiI 10768.900 (x)"])
    p = pla.plan({"species": "Al I", "wavelength_A": 10768.363, "ep_eV": 4.085, "set": "S",
                  "authors_loggf": -2.02, "authors_sigma_dex": 0.06})
    with pytest.raises(pla.LineAddError, match="would move"):
        pla.apply([p], citation={})


def test_cell_identity_hashes_only_its_own_window(tmp_path):
    f = tmp_path / "t.csv"
    f.write_text("species,wavelength_air_A,log_gf\nAl I,6696.0,-1.5\nAl I,10768.4,-2.0\n")
    a = rm._rows_fingerprint(f, "wavelength_air_A", 6600, 6800)
    f.write_text("species,wavelength_air_A,log_gf\nAl I,6696.0,-1.5\nAl I,10768.4,-1.9\n")
    assert rm._rows_fingerprint(f, "wavelength_air_A", 6600, 6800) == a      # outside: same
    f.write_text("species,wavelength_air_A,log_gf\nAl I,6696.0,-1.6\nAl I,10768.4,-1.9\n")
    assert rm._rows_fingerprint(f, "wavelength_air_A", 6600, 6800) != a      # inside: moves


def test_synthesis_lists_overlap_by_window():
    names = [p.parent.name for p in rm.synthesis_lists_for(10700, 10800)]
    assert names == ["ispec_ir_9200_13000"]
    assert rm.synthesis_lists_for(6600, 6800) == []


def test_cog_inversion_never_rails():
    """np.interp clamps: an NLTE EW outside the LTE curve read as exactly the bracket edge
    (37 Si I lines at -0.200/+0.200). Out of range or non-monotonic is NaN."""
    import math
    from pipeline.pysme_nlte import cog_invert
    A = [7.31, 7.41, 7.51, 7.61, 7.71]
    cog = [10.0, 12.0, 14.0, 16.0, 18.0]
    assert cog_invert(15.0, cog, A) == pytest.approx(7.56)
    assert math.isnan(cog_invert(19.0, cog, A))
    assert math.isnan(cog_invert(9.0, cog, A))
    assert math.isnan(cog_invert(15.0, [10, 12, 11, 16, 18], A))

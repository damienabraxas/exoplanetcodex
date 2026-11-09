"""RYA-1250 fix: the gf the synthesis READS for solar Mg is the literature's.

`gf_resolver.apply_to_synth_array` scales each synthesis cluster to the nearest canonical row, so
these tests resolve the canonical table the way the synthesis does, cluster by cluster, instead of
reading canonical rows in isolation (PR #585's rows were right in the table and wrong in the
synthesis).
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pipeline import gf_resolver as gfr
from pipeline import published_line_add as pla

ROOT = Path(__file__).resolve().parents[1]
LISTS = {"ir": ROOT / "data/linelists/ispec_ir_9200_13000/atomic_lines.tsv",
         "h": ROOT / "data/linelists/ispec_h_15007_17494/atomic_lines.tsv",
         "k": ROOT / "data/linelists/ispec_k_19452_24846/atomic_lines.tsv"}


def _synth_totals(lst: str, tag: str) -> dict:
    """centroid -> canonical total the synthesis applies to that Mg cluster."""
    d = pd.read_csv(LISTS[lst], sep="\t", usecols=["element", "wave_A", "loggf", "lower_state_eV"],
                    low_memory=False)
    d = d[d.element == tag].reset_index(drop=True)
    key = (12, int(tag.split()[1]))
    wl, ep, gf = (d.wave_A.to_numpy(float), d.lower_state_eV.to_numpy(float),
                  d.loggf.to_numpy(float))
    out = {}
    for cl in gfr.cluster_physical_lines([key] * len(d), wl, ep):
        w = 10 ** gf[cl]
        cen = float(np.sum(wl[cl] * w) / w.sum())
        try:
            out[round(cen, 3)] = gfr.resolve(key, cen, float(np.mean(ep[cl])))
        except gfr.GfResolutionError:
            out[round(cen, 3)] = None
    return out


@pytest.mark.parametrize("lst,tag,centroid,expected", [
    ("ir", "Mg 1", 12083.282, -1.116),   # 3d1D2-4f3F pair .278+.346 (PR17 exp + Bergemann 2015)
    ("ir", "Mg 1", 12083.649, 0.377),    # 3d1D2-4f1F, PR17 exp
    ("ir", "Mg 1", 10965.445, -0.184),   # Bergemann's three components summed
    ("ir", "Mg 1", 10811.112, 0.401),    # Bergemann's six components summed, not VALD's 7
    ("ir", "Mg 1", 9983.188, -2.177),    # AGSS21 line now in the synthesis list
    ("h", "Mg 1", 15765.825, 0.471),
    ("h", "Mg 1", 15886.206, -1.340),
    ("h", "Mg 1", 15879.558, -1.226),    # AGSS21 feature total
    ("k", "Mg 1", 19723.532, -0.739),    # NIST grades the dominant component; total = VALD sum
    ("k", "Mg 1", 19733.868, -0.517),
])
def test_synthesis_reads_the_published_feature_total(lst, tag, centroid, expected):
    got = _synth_totals(lst, tag).get(centroid)
    assert got is not None and got == pytest.approx(expected, abs=0.002)


def test_mg2_10092_is_two_components_summing_to_agss21():
    t = _synth_totals("ir", "Mg 2")
    a, b = t.get(10092.095), t.get(10092.217)
    assert a is not None and b is not None
    assert np.log10(10 ** a + 10 ** b) == pytest.approx(1.283, abs=0.002)


def test_5183_carries_the_lab_gf_not_a_dropped_sign():
    c = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    r = c[c.line_id == "gf_087247"].iloc[0]
    assert r.log_gf == pytest.approx(-0.168) and r.gf_tier == "LAB"


def test_no_added_row_shadows_an_existing_line():
    """A row the line-add path wrote must not sit within the resolver tolerance of another row:
    the resolver reads only the nearest, so one of the two is dead weight (PR #585 left 14)."""
    c = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    added = c[c.adjudication_status.astype(str).str.startswith("added_")]
    for _, r in added.iterrows():
        near = c[(c.species == r.species) & (c.line_id != r.line_id)
                 & ((c.wavelength_air_A - r.wavelength_air_A).abs() <= gfr.WTOL)
                 & ((c.excitation_potential_eV - r.excitation_potential_eV).abs() <= gfr.EPTOL)]
        assert near.empty, f"{r.line_id} {r.species} {r.wavelength_air_A} shadows {list(near.line_id)}"


def test_line_add_refuses_a_feature_spanning_two_synthesis_clusters(monkeypatch):
    comps = [{"wavelength": 10092.095, "log_gf": 0.91}, {"wavelength": 10092.217, "log_gf": 1.02}]
    monkeypatch.setattr(pla, "components", lambda sp, w, ep: (comps, Path("x.txt")))
    with pytest.raises(pla.LineAddError, match="more than one synthesis cluster"):
        pla.plan({"species": "Mg II", "wavelength_A": 10092.2, "ep_eV": 11.6297, "set": "S",
                  "authors_loggf": 1.283, "authors_sigma_dex": 0.01})

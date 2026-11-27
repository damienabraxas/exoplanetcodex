"""RYA-731 -- model_atmosphere for a species Scott+2015 does not tabulate (Ni II), priced from
the unit's own ATLAS9 / MARCS.GES deck-LTE pair instead of holding the product."""
import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import rya1230_cno_budget as B  # noqa: E402

EVID = ROOT / "data/audit/rya731_ni2_clear2026/nominal/sirius"
STEM = "NiII_3000_3780_kpno_solar_atlas_solar_kpno_kurucz2005_corrected_SYNTH_SET-NI_CLEAR2026_{}_lines.csv"


def _price(nominal_dir: Path, token: str):
    stem = STEM.format(token)
    acc = B._acc(pd.read_csv(nominal_dir / stem))
    comps, notes = [], []
    B._grid_pair_atmosphere(stem, nominal_dir, acc, len(acc), "median", "x", comps, notes)
    return comps, notes, acc


def _pair_median(marcs: Path, atlas: Path) -> float:
    m, a = B._acc(pd.read_csv(marcs)), B._acc(pd.read_csv(atlas))
    m, a = m.assign(wavelength_air_A=m.wavelength_air_A.round(3)), a.assign(wavelength_air_A=a.wavelength_air_A.round(3))
    j = m.merge(a, on="wavelength_air_A", suffixes=("_m", "_a"))
    return float((j.abundance_m - j.abundance_a).median())


@pytest.fixture
def legs(tmp_path):
    for u in ("unit15", "unit16", "unit17"):
        shutil.copytree(EVID / u / "nominal", tmp_path / u / "nominal")
    return tmp_path


def test_marcs_and_atlas9_products_carry_the_same_paired_grid_term(legs):
    want = abs(_pair_median(legs / "unit15/nominal" / STEM.format("synth-1D-LTE-gerber"),
                            legs / "unit15/nominal" / STEM.format("1D-LTE")))
    got = {}
    for unit, tok in (("unit15", "synth-1D-LTE-gerber"), ("unit15", "1D-LTE"),
                      ("unit17", "1D-LTE"), ("unit17", "ENGINE-B"), ("unit16", "1D-LTE")):
        comps, notes, acc = _price(legs / unit / "nominal", tok)
        assert comps, (unit, tok, notes)
        assert comps[0]["evidence"]["n_paired"] == len(acc)
        got[(unit, tok)] = comps[0]["sigma_dex"]
    assert all(v == pytest.approx(want, abs=1e-4) for v in got.values())
    # the ATLAS9 product in a unit without a MARCS comparand found the peer unit's
    comps, _, _ = _price(legs / "unit17" / "nominal", "1D-LTE")
    assert comps[0]["evidence"]["comparand"].startswith("unit15/")


def test_a_mean3d_product_gets_half_its_own_3d_effect_with_the_1d_grid_term(legs):
    stem = STEM.format("synth-mean3D-LTE-gerber-stagger")
    acc = B._acc(pd.read_csv(legs / "unit16/nominal" / stem))
    comps, notes = [], []
    B._grid_pair_atmosphere(stem, legs / "unit16/nominal", acc, len(acc), "median", "x", comps, notes)
    assert comps, notes
    ev = comps[0]["evidence"]
    one = legs / "unit15/nominal" / STEM.format("synth-1D-LTE-gerber")
    assert ev["half_mean3d_effect_dex"] == pytest.approx(
        abs(_pair_median(legs / "unit16/nominal" / stem, one)) / 2, abs=1e-6)
    grid = abs(_pair_median(one, legs / "unit15/nominal" / STEM.format("1D-LTE")))
    assert ev["one_d_grid_dex"] == pytest.approx(grid, abs=1e-4)
    assert comps[0]["sigma_dex"] == pytest.approx((grid ** 2 + ev["half_mean3d_effect_dex"] ** 2) ** 0.5, abs=1e-4)


def test_a_peer_with_a_different_atlas9_pool_is_refused(legs):
    # unit17 has no MARCS sibling; perturb unit15's ATLAS9 pool so the peer no longer matches
    f = legs / "unit15/nominal" / STEM.format("1D-LTE")
    d = pd.read_csv(f)
    i = d.index[d.in_aggregate == True][0]  # noqa: E712
    d.loc[i, "abundance"] = float(d.loc[i, "abundance"]) + 0.01
    d.to_csv(f, index=False)
    comps, notes, _ = _price(legs / "unit17" / "nominal", "1D-LTE")
    assert not comps
    assert any("no ATLAS9/MARCS.GES pair" in n for n in notes)

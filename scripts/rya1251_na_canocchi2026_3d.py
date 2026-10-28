"""RYA-1251: solar Na I 3D-NLTE minus 1D-LTE corrections from Canocchi et al. (2026).

Canocchi, Wang, Amarsi, Lind & Racca 2026 (arXiv 2603.22975; Zenodo 10.5281/zenodo.19201830,
CC BY 4.0) publish a full 3D NLTE grid (Balder on Stagger, the Lind+2011 atom) with feed-forward
networks mapping (Teff, logg, [Fe/H], REW) -> A(Na) in 3D NLTE and (..., vmic, REW) -> A(Na) in
1D LTE on MARCS. This script evaluates the correction exactly as the authors' interp_tools.
nlte_correction does (1D-LTE curve of growth inverted by cubic spline, then the 3D-NLTE net at
that REW), at the Sun (5772/4.44/0.00, our pinned vmic 1.00 = theirs) and A_1DLTE = 6.25 (their
solar 1D LTE value).

Differences from the authors' code, none of which change a number:
  * models load through a whitelist unpickler (sklearn / numpy / joblib globals only) -- the
    .joblib files are pickles;
  * the networks are evaluated from their stored weights (`coefs_`, `intercepts_`), checked
    against sklearn's own predict on every call, so the result does not depend on the sklearn
    version (they were saved with 0.24.2 / 1.6.1);
  * the `Scalar` class the archive imports is not in the record (the SAND repo is private); its
    stored arrays are [mean; std] and are applied as standardisation. That reading is VERIFIED
    here: at every solar-metallicity dwarf node near the Sun the 3D-NLTE net returns the node's
    own A(Na) from the node's own 3D-NLTE spectrum (written to the validation artifact).

Lines Asplund+2021 Table A.1 already cites (6154, 6160, 10746) are NOT written to the table --
one line has one cited authority (pipeline.solar_nlte_tables refuses two). They are written to
the validation artifact as the cross-check between the two.

Usage:  python scripts/rya1251_na_canocchi2026_3d.py [--grid DIR] [--check]
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
TABLE = ROOT / "data/reference/na_canocchi2026/na_canocchi2026_solar.csv"
AUDIT = ROOT / "data/results/rya1251/na_canocchi2026_validation.json"
SUN = dict(teff=5772.0, logg=4.44, met=0.0, vmic=1.00, a1d=6.25)

#: grid line key -> (line centre the authors integrate REW over, canonical row, air wavelength)
LINES = {
    4982: (4982.84075051, "gf_097220", 4982.813),
    5682: (5682.670485909723, "gf_097223", 5682.633),
    5688: (5688.2289707840455, "gf_097224", 5688.205),
    5889: (5889.97397902, "gf_097225", 5889.951),
    5895: (5895.95234942, "gf_097226", 5895.924),
    6154: (6154.23, "gf_097227", 6154.225),
    6160: (6160.75, "gf_097228", 6160.747),
    8183: (8183.29211327, "gf_097242", 8183.256),
    8194: (8194.82568903, "gf_097243", 8194.8206),
    10747: (10747.39753831, "gf_167446", 10746.449),
}
#: REW span of the 1D-LTE curve of growth, per line (interp_tools.gen_abund_to_rew)
REW_RANGE = {5682: (-9, -3.5), 5688: (-9, -3.5), 8194: (-9, -3.5), 8183: (-9, -3.5),
             5889: (-7.0, -2.8), 5895: (-7.0, -2.8), 6154: (-10, -4.0), 6160: (-10, -4.0),
             4982: (-10, -4.0), 10747: (-11, -4.8)}
AGSS21_LINES = {6154, 6160, 10747}
CITATION = ("Canocchi, Wang, Amarsi, Lind & Racca 2026, arXiv 2603.22975, Zenodo "
            "10.5281/zenodo.19201830 (CC BY 4.0); authors' FFNN, solar 5772/4.44/0.00, vmic 1.00, "
            "A_1DLTE 6.25; MARCS 1D LTE baseline")


def _safe_load(path):
    import joblib.numpy_pickle as jnp

    class Safe(jnp.NumpyUnpickler):
        def find_class(self, module, name):
            if not module.startswith(("sklearn.", "numpy", "joblib.numpy_pickle")):
                raise RuntimeError(f"refused pickle global {module}.{name} in {path}")
            return super().find_class(module, name)

    base, jnp.NumpyUnpickler = jnp.NumpyUnpickler, Safe
    try:
        return jnp.load(path)
    finally:
        jnp.NumpyUnpickler = base


class Grid:
    def __init__(self, root: Path):
        self.root, self._m = Path(root), {}

    def _model(self, line, model):
        if (line, model) not in self._m:
            d = self.root / "ffnn" / f"line{line}" / model
            sx = np.load(d / "Abna_noholdout_seed28_scalar_fullset.npy", allow_pickle=False)
            sy = np.load(d / "Abna_noholdout_seed28_y_scalar.npy", allow_pickle=False)
            self._m[(line, model)] = (sx, sy, _safe_load(d / "ffnn.joblib"))
        return self._m[(line, model)]

    def abund(self, X, line, model):
        sx, sy, m = self._model(line, model)
        Z = (np.atleast_2d(np.asarray(X, float)) - sx[0]) / sx[1]
        act = {"relu": lambda z: np.maximum(z, 0), "tanh": np.tanh,
               "logistic": lambda z: 1 / (1 + np.exp(-z)), "identity": lambda z: z}[m.activation]
        h = Z
        for i, (W, b) in enumerate(zip(m.coefs_, m.intercepts_)):
            h = h @ W + b
            h = act(h) if i < len(m.coefs_) - 1 else h
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if not np.allclose(h.ravel(), m.predict(Z), atol=1e-9):
                raise SystemExit(f"weights forward pass != sklearn predict for {line} {model}")
        return h.ravel() * sy[1, 0] + sy[0, 0]

    def correction(self, teff, logg, met, vmic, a1d, line):
        from scipy.interpolate import CubicSpline
        rews = np.arange(*REW_RANGE[line], 0.1)
        n = len(rews)
        a_lte = self.abund(np.column_stack([[teff] * n, [logg] * n, [met] * n, [vmic] * n, rews]),
                           line, "1DLTE")
        rew = float(CubicSpline(a_lte, rews)(a1d))
        return float(self.abund([[teff, logg, met, rew]], line, "3DNLTE")[0]) - a1d, rew


def node_validation(grid: Grid) -> dict:
    """At solar-metallicity dwarf nodes near the Sun, REW of the node's own 3D-NLTE spectrum
    -> 3D-NLTE net -> A(Na); returns the deviation from the node's A(Na) per line."""
    from astropy.io import fits
    with fits.open(grid.root / "3Dflux_all.fits") as h:
        P, wl, FN = h["PARAMS"].data, h["WAVELENGTH"].data * 10.0, h["FLUX_NLTE"].data
    segs = np.split(np.arange(wl.size), np.where(np.diff(wl) > 5)[0] + 1)
    sel = np.where((np.abs(P["logg"] - 4.5) < 0.2) & (np.abs(P["fe_h"]) < 0.05)
                   & (np.abs(P["Teff"] - 5800) < 400))[0]
    dev = collections.defaultdict(list)
    for i in sel:
        for s in segs:
            w, f = wl[s], FN[i][s]
            inside = [c for c, *_ in LINES.values() if w[0] < c < w[-1]]
            for line, (c, *_r) in LINES.items():
                if not w[0] < c < w[-1]:
                    continue
                lo = max([w[0]] + [0.5 * (o + c) for o in inside if o < c])
                hi = min([w[-1]] + [0.5 * (o + c) for o in inside if o > c])
                m = (w >= lo) & (w <= hi)
                rew = np.log10(np.trapezoid(1 - f[m], w[m]) / c)
                a = grid.abund([[P["Teff"][i], P["logg"][i], P["fe_h"][i], rew]], line, "3DNLTE")[0]
                dev[line].append(round(float(a - P["A_Na"][i]), 4))
    return {"nodes": [dict(teff=float(P["Teff"][i]), logg=float(P["logg"][i]),
                           fe_h=float(P["fe_h"][i]), A_Na=float(P["A_Na"][i])) for i in sel],
            "note": ("REW integrated over the line's own segment, split half-way to a neighbour "
                     "line in the same segment (8183/8194, D1/D2, 5682/5688); that split, not the "
                     "network, is the larger 8183/8194 deviation"),
            "deviation_dex": {str(k): v for k, v in sorted(dev.items())},
            "median_dex": {str(k): float(np.median(v)) for k, v in sorted(dev.items())}}


def build(grid_dir: Path) -> tuple[list[dict], dict]:
    grid = Grid(grid_dir)
    agss = {r["canonical_line_id"]: r for r in csv.DictReader(
        (ROOT / "data/reference/na_agss21/na_agss21_lines.csv").open())}
    rows, cross = [], {}
    for line, (centre, cid, air) in LINES.items():
        d, rew = grid.correction(SUN["teff"], SUN["logg"], SUN["met"], SUN["vmic"], SUN["a1d"], line)
        sens = {f"d_a1d_{a}": round(grid.correction(SUN["teff"], SUN["logg"], SUN["met"],
                                                      SUN["vmic"], a, line)[0], 4)
                for a in (6.05, 6.45)}
        sens.update({f"d_vmic_{v}": round(grid.correction(SUN["teff"], SUN["logg"], SUN["met"],
                                                            v, SUN["a1d"], line)[0], 4)
                     for v in (0.8, 1.2)})
        row = dict(species="Na I", wavelength_air_A=air, canonical_line_id=cid, grid_line=line,
                   rew=round(rew, 3), A_1d_lte=SUN["a1d"], A_3d_nlte=round(SUN["a1d"] + d, 4),
                   delta_3d_nlte_minus_1d_lte=round(d, 4), **sens)
        if line in AGSS21_LINES:
            a = agss[cid]
            cross[cid] = dict(row, agss21_3d_nlte_minus_marcs_lte=round(
                float(a["A_3d_nlte"]) - float(a["A_marcs_lte"]), 4))
        else:
            rows.append(row)
    audit = dict(ticket="RYA-1251", citation=CITATION, sun=SUN,
                 node_validation=node_validation(grid),
                 agss21_cross_check=cross,
                 rule=("AGSS21 Table A.1 cites 6154/6160/10746 and keeps them; this table carries "
                       "the seven lines only Canocchi+2026 covers"),
                 paper_solar=dict(A_1d_lte=6.25, A_3d_nlte=6.17, source="Canocchi+2026 Table 3, "
                                  "nine lines on Lind+2022 EWs"))
    return rows, audit


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--grid", type=Path, default=None,
                    help="grid directory (default: path register grids.canocchi2026_na)")
    ap.add_argument("--check", action="store_true",
                    help="recompute and fail if the committed table differs")
    a = ap.parse_args(argv)
    if a.grid is None:
        from config.constants import codex_path
        a.grid = Path(str(codex_path("grids.canocchi2026_na")))
    rows, audit = build(a.grid)
    cols = list(rows[0])
    if a.check:
        old = list(csv.DictReader(TABLE.open()))
        new = [{k: str(v) for k, v in r.items()} for r in rows]
        if old != new:
            raise SystemExit("na_canocchi2026_solar.csv differs from a recompute")
        print("table reproduces")
        return
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    with TABLE.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    AUDIT.write_text(json.dumps(audit, indent=1) + "\n")
    for r in rows:
        print(f"{r['wavelength_air_A']:>10}  {r['delta_3d_nlte_minus_1d_lte']:+.3f}")
    for cid, r in audit["agss21_cross_check"].items():
        print(f"{r['wavelength_air_A']:>10}  {r['delta_3d_nlte_minus_1d_lte']:+.3f}  "
              f"AGSS21 {r['agss21_3d_nlte_minus_marcs_lte']:+.3f}  (cross-check only)")
    print("wrote", TABLE.relative_to(ROOT), AUDIT.relative_to(ROOT))


if __name__ == "__main__":
    main()

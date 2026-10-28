#!/usr/bin/env python3
"""RYA-1250: measure the `hfs_isotopes` budget term for the strong, split Mg lines.

The budget zeroes the term only in the optically thin limit (every split component <= 0.10
deep); a stronger split line has no measured term and HOLDs its product -- 58 Mg products,
all of H and K. This measures a BOUND per line: the synthesis as run (components split, as
in the list) against the same physical line collapsed to ONE component at the gf-weighted
centroid with the total gf. The abundance that makes the collapsed line's synthetic EW equal
the split line's at A0 is A'; |A' - A0| is the largest error any mistake in the splitting
can cause, since collapse is the extreme of every splitting error.

Run on Sirius (iSpec venv); writes data/results/rya1250/mg_hfs_collapse.csv, which
`scripts/rya1230_cno_budget.py` reads.

    ~/venv_rya1103/bin/python scripts/rya1250_mg_hfs_term.py
"""
from __future__ import annotations

import glob
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import gf_resolver as gfr  # noqa: E402
from pipeline.element_prepare import SYNTH_LISTS  # noqa: E402

OUT = ROOT / "data/results/rya1250/mg_hfs_collapse.csv"
A0 = 7.56            # Amarsi & Grevesse 2026; the bound is evaluated at the solar value
HALF_WIDTH_A = 0.3   # fit half-width handed to synth_ew (integration is >= 1 A either side)
SLICE_A = 6.0        # list rows kept around the line; integration reaches 1 A


def pool_lines() -> pd.DataFrame:
    """Every in-aggregate Mg line of a built product whose canonical row is split."""
    from pipeline.gf_resolver import EPTOL, WTOL
    g = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    g = g[g.species.isin(["Mg I", "Mg II"]) & (g.hfs_n_components > 1)]
    seen = {}
    for f in glob.glob(str(ROOT / "data/results/band_products/Mg*_lines.csv")):
        d = pd.read_csv(f)
        d = d[d.in_aggregate.astype(str) == "True"]
        for _, l in d.iterrows():
            sp = f"Mg {str(l.ion).strip()}"
            m = g[(g.species == sp) & ((g.wavelength_air_A - l.wavelength_air_A).abs() <= WTOL)
                  & ((g.excitation_potential_eV - l.ep_eV).abs() <= EPTOL)]
            if len(m):
                r = m.loc[(m.wavelength_air_A - l.wavelength_air_A).abs().idxmin()]
                seen[r.physical_id] = r
    return pd.DataFrame(list(seen.values()))


def _list_for(w: float):
    path = None
    for lo, f in SYNTH_LISTS:
        if w >= lo:
            path = f
    return None if path is None else str(path)


def collapse(ll: np.ndarray, tag: str, w: float, ep: float) -> tuple[np.ndarray, int]:
    """The list with the physical line at (w, ep) collapsed to one component."""
    near = np.where((ll["element"] == tag) & (np.abs(ll["wave_A"] - w) < 0.5)
                    & (np.abs(ll["lower_state_eV"] - ep) < 0.05))[0]
    cl = gfr.cluster_physical_lines([0] * len(near), ll["wave_A"][near], ll["lower_state_eV"][near])
    cl = min(cl, key=lambda c: np.min(np.abs(ll["wave_A"][near[c]] - w)))
    idx = near[cl]
    lin = 10.0 ** ll["loggf"][idx]
    keep = idx[np.argmax(lin)]
    out = ll.copy()
    out["wave_A"][keep] = float(np.sum(ll["wave_A"][idx] * lin) / lin.sum())
    if "wave_nm" in (out.dtype.names or ()):
        out["wave_nm"][keep] = out["wave_A"][keep] / 10.0
    out["loggf"][keep] = float(np.log10(lin.sum()))
    return np.delete(out, [i for i in idx if i != keep]), len(idx)


def measure(r) -> dict:
    from pipeline.abundances_derive import _synth_flux_at_abund
    from pipeline.nearuv_synth import build_solar_context
    from pipeline.synth_ew import synthetic_ew_mA
    w, ep = float(r.wavelength_air_A), float(r.excitation_potential_eV)
    ctx = build_solar_context("Mg", 100000.0, linelist_file=_list_for(w))
    ll = ctx["linelist"]
    ll = ll[np.abs(ll["wave_A"] - w) < SLICE_A]
    tag = f"Mg {1 if r.species.endswith(' I') else 2}"
    coll, n = collapse(ll, tag, w, ep)
    kw = {k: ctx[k] for k in ("atmosphere", "teff", "logg", "feh", "vturb", "isotopes",
                              "solar_abund", "atom_code")}
    kw.update(element="Mg", R=ctx["resolving_power"], macroturbulence=ctx["macroturbulence"],
              vsini=ctx["vsini"], tmp_dir=tempfile.mkdtemp(prefix="rya1250_hfs_"))

    def ew(lst, a):
        return synthetic_ew_mA(_synth_flux_at_abund, centre_A=w, abundance=a,
                               fit_half_width_A=HALF_WIDTH_A,
                               synth_kwargs=dict(kw, linelist=lst))
    target = ew(ll, A0)
    # local slope of ln EW in A for the collapsed line (the shift is small), then one
    # correction step and a check at the solution
    e0, e1 = ew(coll, A0), ew(coll, A0 + 0.05)
    slope = (np.log(e1) - np.log(e0)) / 0.05
    a = A0 + (np.log(target) - np.log(e0)) / slope
    ea = ew(coll, a)
    a = a + (np.log(target) - np.log(ea)) / slope
    ea = ew(coll, a)
    # what the solve leaves unmatched (synthesis/trapezoid noise, ~0.3% in EW) is added to
    # the bound rather than discarding the line
    resid = float(np.log(ea / target))
    ok = bool(np.isfinite(ea) and np.isfinite(slope) and slope > 0)
    a2 = a
    return dict(species=r.species, wavelength_air_A=w, ep_eV=ep, physical_id=r.physical_id,
                n_components_in_synthesis=n, ew_split_mA=round(target, 3),
                ew_collapsed_at_A0_mA=round(e0, 3), dlnEW_dA=round(slope, 4),
                resid_ln_ew=round(resid, 5),
                A_collapsed=round(a2, 4) if ok else np.nan,
                dA_collapse_dex=(round(float(np.sign(a2 - A0) * (abs(a2 - A0) + abs(resid) / slope)), 4)
                                 if ok else np.nan))


def main() -> int:
    lines = pool_lines()
    rows = []
    for _, r in lines.sort_values("wavelength_air_A").iterrows():
        m = measure(r)
        print(f"{m['species']} {m['wavelength_air_A']:.3f} n={m['n_components_in_synthesis']} "
              f"EW {m['ew_split_mA']} mA  dA(collapse) {m['dA_collapse_dex']}", flush=True)
        rows.append(m)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"wrote {OUT.relative_to(ROOT)}: {len(rows)} split, non-thin Mg lines")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

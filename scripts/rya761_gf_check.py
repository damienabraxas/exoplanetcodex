#!/usr/bin/env python3
"""RYA-761 closing question: of the recovered lines, how many could reach a product?

RYA-761: *"A recovered line only reaches the product if its gf is also laboratory-grade
— currently 29 of 101 clear that bar, so a naive scaling would suggest ~80 new product
lines, but that assumes the recoverable set has the same gf provenance mix as the clean
set, which is exactly the kind of assumption RYA-760 just falsified for a different cut.
Measure it."*

So this measures it: the gf provenance of the RECOVERED lines against the provenance of
the CLEAN (already in-aggregate) lines, from the synthesis line list's own
`reference_code`. It does not project.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, os.environ.get("ISPEC_DIR", "/srv/codex/engines/ispec_src"))

from pipeline.abundances_derive import _load_synth_resources  # noqa: E402

REC = ROOT / "data/results/rya761/FeI_6910_9199_kpno_solar_atlas_PROFILEFIT_blend_recovery.csv"
EW = Path("/mnt/codex-data/codex/rya713/data/measured/band_ew/"
          "FeI_6910_9199_kpno_solar_atlas_PROFILEFIT_ew.csv")

# Semi-empirical / theoretical source families. Anything else we treat as measured;
# the point of the comparison is the MIX, not a grade letter (RYA-711 warns our grade
# vocabulary collides with NIST's).
SEMI_PREFIXES = ("K75", "K88", "KUR", "K07", "K14")


def gf_ref(ll, w_A, elem_prefix="FE", tol=0.05):
    w = np.asarray(ll["wave_A"], dtype=float)
    el = np.asarray([str(x) for x in ll["element"]])
    gf = np.asarray(ll["loggf"], dtype=float)
    idx = [i for i in np.where(np.abs(w - w_A) < tol)[0]
           if el[i].strip().upper().startswith(elem_prefix)]
    if not idx:
        return None
    i = max(idx, key=lambda k: gf[k])
    return str(ll["reference_code"][i]).strip()


def classify(ref):
    if ref is None:
        return "NOT-IN-LINELIST"
    return ("semi-empirical" if any(ref.upper().startswith(p) for p in SEMI_PREFIXES)
            else "measured/other")


def main() -> None:
    ll, _, _ = _load_synth_resources()
    rec = pd.read_csv(REC)
    ew = pd.read_csv(EW)

    recovered = rec[rec.recovered.fillna(False)]
    clean = ew[ew.in_aggregate.fillna(False)]

    out = {}
    for name, waves in (("RECOVERED", recovered.wavelength_air_A.values),
                        ("CLEAN (already in-aggregate)", clean.wavelength_air_A.values)):
        refs = [gf_ref(ll, float(x)) for x in waves]
        cls = pd.Series([classify(r) for r in refs])
        out[name] = (len(waves), cls.value_counts(), pd.Series(refs).value_counts())

    print("\n" + "=" * 78)
    print("gf PROVENANCE — recovered vs clean. Measured, not assumed.")
    print("=" * 78)
    for name, (n, counts, refs) in out.items():
        print(f"\n{name}  (n = {n})")
        for k, v in counts.items():
            print(f"    {k:<20} {v:4d}   ({100.0*v/max(n,1):.1f}%)")
        print(f"    top sources: "
              f"{', '.join(f'{a}:{b}' for a, b in list(refs.items())[:6])}")

    n_rec = out["RECOVERED"][0]
    m_rec = out["RECOVERED"][1].get("measured/other", 0)
    n_cln = out["CLEAN (already in-aggregate)"][0]
    m_cln = out["CLEAN (already in-aggregate)"][1].get("measured/other", 0)
    print("\n" + "-" * 78)
    print(f"  recovered lines on a measured gf source: {m_rec} of {n_rec} "
          f"({100.0*m_rec/max(n_rec,1):.1f}%)")
    print(f"  clean lines on a measured gf source:     {m_cln} of {n_cln} "
          f"({100.0*m_cln/max(n_cln,1):.1f}%)")
    print("\n  => the naive assumption RYA-761 flags is that these two rates are equal.")
    print(f"     They are {'COMPARABLE' if abs(m_rec/max(n_rec,1) - m_cln/max(n_cln,1)) < 0.10 else 'DIFFERENT'}"
          f" — so the projection is "
          f"{'defensible' if abs(m_rec/max(n_rec,1) - m_cln/max(n_cln,1)) < 0.10 else 'NOT defensible'}.")


if __name__ == "__main__":
    main()

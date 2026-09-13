#!/usr/bin/env python3
"""RYA-1214 — `cno_synthesis` output -> the band-product schema, so it can be PUBLISHED.

`pipeline.cno_synthesis` writes `data/audit/cno_synthesis/<star>_<region>_cno_*.csv`, which
is its own shape and is not `data/products/solar/<El>.json`. That directory is exactly
where a result goes to be forgotten: `scripts/publish_product.py` exists because "a real,
current, four-arm 3D-NLTE result was produced on Sirius and then could not be found", and
the CNO molecular numbers have been sitting in that audit directory since 2026-06-27 under
a STALE.md nobody could act on. So they get converted and offered to the publisher, which
is the only thing allowed to write the store (RYA-1034).

HOW A MOLECULAR BAND IS KEYED, AND WHY `ion` IS NOT THE ANSWER
--------------------------------------------------------------
The product identity is (element, ion, band, instrument, holding, tier, selector, route,
treatment, line_set) — RYA-1127. A molecular band measures an ELEMENTAL abundance: A(C)
from CH, A(N) from CN, A(O) from OH. There is no ionisation stage involved, and writing
`ion='I'` would be the closest available lie rather than a description.

So the DIAGNOSTIC travels in `selector`, which is a full identity field (RYA-984), as
`MOL-CH_Gband` / `MOL-CN_red` / `ATOM-CI_5052` / `FORB-OI_6300`. `ion` carries the
element's neutral stage because the schema requires a value and the abundance is the
neutral-scale A(X) every other product reports — stated here so it is read as a convention
and not as a claim that CH is a carbon ion. Two diagnostics of one element in one band are
then two products on two keys, which is what they are.

⚠️ n_lines IS THE NUMBER OF LINES IN THE BAND, NOT n_pix. `cno_synthesis` reports `n_pix`,
the pixel count of the fit windows, and calling that `n_lines` would put a 600 in a column
every other product fills with a line count — a fictitious pool size, and a `sigma/sqrt(n)`
downstream would then divide by the wrong thing. The molecular line count in each window is
counted FROM THE .bsyn LISTS and `n_pix` is carried separately under its own name.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CNOSYNTH = ROOT / "data" / "audit" / "cno_synthesis"
OUT = ROOT / "data" / "results" / "band_products"

#: Our diagnostic -> (selector prefix, the band name the product is filed under, the
#: molecule whose lines carry the window, or None for an atomic diagnostic).
DIAG = {
    "CH_Gband": ("MOL", "VIS", "12CH"),
    "C2_Swan":  ("MOL", "VIS", "12C12C"),
    "CN_red":   ("MOL", "VIS", "12C14N"),
    "CI_5052":  ("ATOM", "VIS", None),
    "CI_5380":  ("ATOM", "VIS", None),
    "OI_6300":  ("FORB", "VIS", None),
    "NH_AX":    ("MOL", "near-UV", "NH"),
    "OH_AX":    ("MOL", "near-UV", "OH"),
}
#: Which holding each region's spectrum came from. Read from the region, never guessed:
#: two holdings of one instrument are two different PRODUCTS (RYA-1026).
REGION_HOLDING = {
    "vis": ("harps", "solar_harps_molecfit_corrected"),
    "nearuv": ("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected"),
}


def _windows_for(key: str) -> tuple:
    from pipeline import cno_synthesis as cs
    for tup in (cs.VIS_DIAGNOSTICS, cs.NEARUV_DIAGNOSTICS):
        for d in tup:
            if d.key == key:
                return d.windows_A
    return ()


def _count_molecular_lines(molecule: str, windows) -> int | None:
    """Lines of `molecule` inside the fit windows, counted from the vendored .bsyn lists."""
    base = ROOT / "data" / "linelists" / "molecular" / "turbospectrum"
    cand = list(base.rglob(f"{molecule}_*.bsyn")) + list(base.rglob(f"*{molecule}*rovib*.bsyn"))
    if not cand:
        return None
    total = 0
    for f in cand:
        m = re.match(r".*_(\d+)-(\d+)\.bsyn", f.name)
        if m:
            lo_nm, hi_nm = float(m.group(1)), float(m.group(2))
            if not any(lo_nm <= lo / 10 <= hi_nm or lo_nm <= hi / 10 <= hi_nm
                       for lo, hi in windows):
                continue
        with f.open() as fh:
            for ln in fh:
                p = ln.split()
                if not p:
                    continue
                try:
                    w = float(p[0])
                except ValueError:
                    continue
                if any(lo <= w <= hi for lo, hi in windows):
                    total += 1
    return total


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--star", default="solar")
    ap.add_argument("--region", action="append", required=True)
    a = ap.parse_args()
    sys.path.insert(0, str(ROOT))
    OUT.mkdir(parents=True, exist_ok=True)

    written = []
    for region in a.region:
        pb = CNOSYNTH / f"{a.star}_{region}_cno_per_band.csv"
        prod = CNOSYNTH / f"{a.star}_{region}_cno_product.csv"
        if not pb.exists():
            print(f"  {region}: no per-band artifact at {pb.name} — skipped")
            continue
        bands = pd.read_csv(pb)
        elem_unc = pd.read_csv(prod).set_index("element")
        inst, holding = REGION_HOLDING[region]
        for _, r in bands.iterrows():
            key = str(r.key)
            if key not in DIAG:
                print(f"  {region}/{key}: no product mapping — skipped, not guessed")
                continue
            kind, band, molecule = DIAG[key]
            el = str(r.element)
            unc = elem_unc.loc[el] if el in elem_unc.index else None
            # sigma_stat is the DIAGNOSTIC's own curvature; the element-level number in
            # `_cno_product.csv` is a scatter over several diagnostics and belongs to the
            # element, not to this band (RYA-712: a product carries its OWN sigma).
            s_stat = float(r.sigma_fit) if pd.notna(r.sigma_fit) else float("nan")
            s_syst = (float(unc["sigma_sys"]) if unc is not None
                      and pd.notna(unc["sigma_sys"]) else float("nan"))
            windows = _windows_for(key)
            n_mol = _count_molecular_lines(molecule, windows) if molecule else None
            row = {
                "element": el, "ion": "I", "band": band, "instrument": inst,
                "treatment": "1D-LTE", "handler": "CNOSynthesis(RYA-237)",
                "A": round(float(r.A_X), 3),
                "n_lines": (int(n_mol) if n_mol else 1),
                "n_excluded": 0,
                "stat_dex": (round(s_stat, 4) if np.isfinite(s_stat) else ""),
                "syst_dex": (round(s_syst, 4) if np.isfinite(s_syst) else ""),
                "stat_basis": ("measured — 1 sigma from the chi2 curvature of THIS band's "
                               "fit (pipeline.fit_constraint.curvature_sigma, rescaled to "
                               "red_chi2 = 1); NaN where the curvature was not measurable, "
                               "never 0 and never the retired 1.000 clip"),
                "dominant": "", "route": "SYNTH", "scale": "1D-LTE", "model": "none",
                "atmos": "atlas9", "gf": "canonical", "route_basis": "handler",
                "deck": "none",
                "selector": f"{kind}-{key}",
                "n_pix": int(r.n_pix),
                "red_chi2": round(float(r.red_chi2), 3),
                "constrained": bool(r.get("constrained", True)),
                "diagnostic_kind": kind,
                "molecule": molecule or "",
            }
            stem = (f"{el}I_{kind}_{key}_{inst}_{holding}_SYNTH_{kind}-{key}")
            f = OUT / f"{stem}_products.csv"
            pd.DataFrame([row]).to_csv(f, index=False)
            written.append({"region": region, "key": key, "file": f.name,
                            "A": row["A"], "selector": row["selector"],
                            "n_lines": row["n_lines"], "molecule": molecule})
            print(f"  {region:7s} {key:10s} {kind:4s} A({el})={row['A']:.3f}  "
                  f"n_lines={row['n_lines']}  n_pix={row['n_pix']}  -> {f.name}")

    (CNOSYNTH / "rya1214_product_shim.prov.json").write_text(json.dumps({
        "ticket": "RYA-1214",
        "what": "cno_synthesis per-band output mapped onto the band-product schema so "
                "publish_product.py can write it into data/products/solar/<El>.json",
        "keying": "the DIAGNOSTIC travels in `selector` (MOL-/ATOM-/FORB- prefix), which "
                  "is a full identity field (RYA-984). `ion='I'` is a schema convention, "
                  "not a claim: a molecular band measures an ELEMENTAL abundance and has "
                  "no ionisation stage.",
        "n_lines": "counted from the vendored .bsyn lists inside the fit windows. NOT "
                   "n_pix — calling a 600-pixel window a 600-line pool would put a "
                   "fictitious pool size in a column every other product fills with a "
                   "line count, and a sigma/sqrt(n) downstream would divide by it.",
        "rows": written,
    }, indent=2) + "\n")
    print(f"\n  {len(written)} product artifact(s) written to "
          f"{OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

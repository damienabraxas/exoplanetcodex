#!/usr/bin/env python3
"""RYA-1214 — the CO line list for AGSS21's carbon indicator in the CRIRES+ K band, in AIR.

AGSS21's CO first-overtone indicators fall in K: 45 lines of (2-0) and (3-1) in the
conditioned CRIRES+ K product's extent. The vendored `CO_IR_Li2015.dat` (RYA-1182: the
ExoMol Li 2015 redistribution that carries AGSS21's CO matches) is measured here against
those 45 lines:

  * log gf: IDENTICAL — median d log gf 0.0000, median |dEP| 0.0003 eV, 45/45 matched.
    Its column 5 is 2J'+1 (J'=9 -> 19.0), the RYA-1207 writer's convention.
  * wavelengths: VACUUM — |d lambda| 0.002 A in the vacuum frame, ~6 A in air.

`pipeline.ace_co_feasibility` synthesised in the VACUUM frame to match ACE's vacuum column,
which was right for that comparison. Every spectrum this repo fits on the band-product and
cno_synthesis routes is AIR, including the RYA-1214 CRIRES+ K product, so staging the .dat
as-is would put every CO line ~83 km/s off. This writes the band subset converted to air
through `pipeline.wavelength_util` (RYA-501), gf and every other column untouched.

GATE (refuses to write): after conversion, every AGSS21 CO line in the band must have a
list line within 0.05 A in AIR with |d log gf| <= 0.02 dex.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.wavelength_util import air_to_vac, vac_to_air  # noqa: E402

SRC = ROOT / "data/linelists/molecular/turbospectrum/CO/CO_IR_Li2015.dat"
AGSS21 = ROOT / "data/audit/rya1136_cno_intake/molecular_physical_crossmatch.csv"
OUT_DIR = ROOT / "data/linelists/molecular/turbospectrum/CO"
AUDIT = ROOT / "data/audit/rya1214_crires_jk/co_k_bsyn_build.json"
LO_AIR, HI_AIR = 19452.42, 24845.62          # config/synth_bands.yaml K
MATCH_TOL_A, GF_TOL = 0.05, 0.02


def main() -> int:
    lines = SRC.read_text().splitlines()
    header, body = lines[:2], lines[2:]
    lo_vac, hi_vac = float(air_to_vac(LO_AIR)), float(air_to_vac(HI_AIR))
    keep, w_vac = [], []
    for ln in body:
        p = ln.split()
        try:
            w = float(p[0])
        except (ValueError, IndexError):
            continue
        if lo_vac <= w <= hi_vac:
            keep.append(ln)
            w_vac.append(w)
    if not keep:
        raise SystemExit("no CO line in the K band — refusing to write an empty list")
    w_air = vac_to_air(np.asarray(w_vac, float))
    out_rows = [f"{wa:10.3f}{ln[len(ln.split()[0]) + (len(ln) - len(ln.lstrip())):]}"
                for wa, ln in zip(w_air, keep)]
    tab = pd.DataFrame({"w_air": w_air,
                        "ep": [float(l.split()[1]) for l in keep],
                        "loggf": [float(l.split()[2]) for l in keep]})

    m = pd.read_csv(AGSS21)
    a = m[m.species.astype(str).str.contains("12C16O")].copy()
    a["air"] = vac_to_air(a.wavelength_vac_nm.to_numpy(float) * 10.0)
    a = a[a.air.between(LO_AIR, HI_AIR)]
    miss, bad = [], []
    for _, r in a.iterrows():
        c = tab[(tab.w_air - r.air).abs() <= MATCH_TOL_A]
        if c.empty:
            miss.append(round(float(r.air), 3))
            continue
        if (c.loggf - r.published_loggf).abs().min() > GF_TOL:
            bad.append(round(float(r.air), 3))
    if miss or bad or len(a) == 0:
        raise SystemExit(f"CO validation FAILED: {len(a)} AGSS21 lines in band, unmatched "
                         f"{miss[:5]}, gf-off {bad[:5]} — refusing to write")

    n0 = header[0].rsplit(None, 1)[0]
    fn = OUT_DIR / f"16O12C_{int(LO_AIR // 10)}-{int(HI_AIR // 10) + 1}.bsyn"
    with fn.open("w") as fh:
        fh.write(f"{n0}   {len(out_rows)}\n")
        fh.write("'ExoMol Li2015 12C16O, vacuum->air (RYA-1214 K band)'\n")
        fh.write("\n".join(out_rows) + "\n")
    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    AUDIT.write_text(json.dumps({
        "ticket": "RYA-1214", "source": str(SRC.relative_to(ROOT)),
        "written": str(fn.relative_to(ROOT)), "n_lines": len(out_rows),
        "air_span_A": [round(float(w_air.min()), 3), round(float(w_air.max()), 3)],
        "conversion": "vacuum -> air, pipeline.wavelength_util.vac_to_air; gf unchanged",
        "agss21_validation": {"lines_in_band": int(len(a)), "unmatched": 0, "gf_off": 0,
                              "match_tol_A_air": MATCH_TOL_A, "gf_tol_dex": GF_TOL},
    }, indent=2) + "\n")
    print(f"wrote {fn.relative_to(ROOT)}: {len(out_rows)} lines, "
          f"{w_air.min():.2f}-{w_air.max():.2f} A air; AGSS21 {len(a)}/{len(a)} validated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

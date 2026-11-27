#!/usr/bin/env python3
"""RYA-731 -- Ni I lab gf from the two furnace papers, transcribed by hand from the ADS scans.

* Blackwell, Booth, Petford & Laming 1989 (MNRAS 236, 235), Table 2 (p. 240-241): 75 Ni I lines
  298-397 nm, EP 0-0.42 eV, "Log gf absolute value". Oxford furnace absorption; relative values
  0.7 per cent (0.003 dex), absolute scale from the Becker et al. 1974 lifetimes, which scatter
  6.7 per cent (0.028 dex) about the Oxford lifetimes (Sect. 4). sigma = 0.03 dex for every line.
  VALD's code for this paper is BBPL. Wood et al. 2014 add +0.024 dex to it; NOT applied here --
  the paper's own value column is transcribed, and Wood wins wherever both measured a line.
* Huber & Sandeman 1980 (A&A 86, 95), Table 1 (p. 98-99): hook + absorption in a furnace,
  1964-4094 A. Errors are printed as 90 per cent confidence limits on the ABSOLUTE scale
  (footnote d); sigma = that / 1.645. Where a line has hook, absorption and recommended values,
  the recommended one (footnote g) is taken. 3548.185 (footnote k, a blend with the multiplet-3
  line at the same wavelength) is not transcribed: it measures two lines at once.

Only 3000 A and redder is kept (nothing we hold reaches bluer). Neither table prints the upper
level, so the levels come from VALD: the row must match exactly ONE distinct VALD transition
(lower EP, upper EP, J) within VALD_TOL_A and the printed lower EP (Huber prints 2 decimals,
Blackwell 4) and, for Huber, the printed lower J. VALD repeats a Ni line once per isotope, so
the check counts distinct transitions, not rows. A is derived from log gf for the schema only.

Blackwell is written before Huber, and rya1268_adopt_lab_gf.py takes the first lab table that
claims a row: on 3669.24 and 3739.23 the two disagree by 0.41 dex (Huber: one short-separation
hook measurement each, footnote f) and Blackwell's 0.7 per cent furnace value is the one adopted.
Usage: rya731_transcribe_ni1_furnace.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
VALD = ROOT / "data" / "audit" / "rya1268_line_coverage" / "vald_solar_gf_refs.csv"
OUTDIR = ROOT / "data" / "reference" / "lab_gf_tables" / "transcribed"
COLS = ["paper", "species", "wavelength_air_A", "ep_eV", "j_lower", "ep_upper_eV",
        "A_1e6_s", "dA_1e6_s", "loggf", "loggf_sigma_dex", "doi", "vizier"]
LO_A = 3000.0
VALD_TOL_A = 0.03

BLACKWELL_DOI = "10.1093/mnras/236.2.235"
BLACKWELL_SIGMA = 0.03
# Table 2: (air wavelength nm, EP eV, log gf absolute). Rows bluer than 300 nm kept for the record.
BLACKWELL = [
    (298.165, 0.1091, -0.997), (298.413, 0.0000, -1.507), (299.260, 0.0254, -1.217),
    (299.446, 0.0254, -0.969), (300.248, 0.0254, -0.083),
    (300.363, 0.1091, -0.323), (301.200, 0.4228, +0.008), (301.914, 0.0000, -1.363),
    (303.187, 0.0000, -1.883), (303.794, 0.0254, -0.536),
    (304.501, 0.1652, -1.822), (305.432, 0.1091, -0.621), (305.764, 0.2124, -0.376),
    (306.462, 0.1091, -1.187), (308.076, 0.2124, -1.216),
    (309.712, 0.1652, -1.529), (310.155, 0.1091, -0.166), (310.188, 0.4228, -0.246),
    (310.547, 0.2748, -1.526), (312.931, 0.2748, -2.555),
    (313.411, 0.2124, -0.273), (314.572, 0.1652, -2.376), (318.437, 0.2748, -2.260),
    (320.042, 0.0254, -2.992), (322.165, 0.0000, -1.920),
    (322.502, 0.4228, -1.404), (323.296, 0.0000, -1.028), (324.306, 0.0254, -1.311),
    (325.074, 0.4228, -1.845), (327.112, 0.1091, -2.349),
    (328.695, 0.0254, -2.570), (331.566, 0.1091, -1.240), (332.026, 0.1652, -1.432),
    (332.231, 0.4228, -1.403), (336.156, 0.1091, -1.458),
    (336.281, 0.2124, -3.570), (336.617, 0.1652, -1.437), (336.789, 0.0254, -2.865),
    (339.105, 0.0000, -1.061), (340.958, 0.0000, -2.467),
    (341.477, 0.0254, -0.029), (342.074, 0.2748, -3.305), (343.356, 0.0254, -0.683),
    (343.728, 0.0000, -1.207), (344.626, 0.1091, -0.407),
    (345.847, 0.2124, -0.238), (346.165, 0.0254, -0.362), (346.949, 0.2748, -1.835),
    (349.296, 0.1091, -0.265), (350.085, 0.1652, -1.294),
    (351.505, 0.1091, -0.226), (351.977, 0.2748, -1.422), (352.454, 0.0254, -0.007),
    (355.153, 0.1652, -2.700), (356.637, 0.4228, -0.251),
    (357.187, 0.1652, -1.154), (359.771, 0.2124, -1.115), (360.228, 0.1652, -2.192),
    (361.046, 0.1091, -1.164), (361.274, 0.2748, -1.423),
    (361.939, 0.4228, +0.020), (362.473, 0.0000, -2.447), (366.410, 0.2748, -2.130),
    (366.924, 0.1652, -2.692), (367.043, 0.1652, -2.399),
    (368.842, 0.2748, -2.346), (373.923, 0.1652, -2.755), (374.905, 0.0000, -3.347),
    (377.557, 0.4228, -1.408), (379.361, 0.2748, -2.989),
    (380.714, 0.4228, -1.220), (383.169, 0.4228, -2.272), (383.287, 0.1652, -3.633),
    (385.830, 0.4228, -0.951), (397.356, 0.4228, -2.620),
]

HUBER_DOI = ""   # A&A 86, 95 (1980) has no DOI; ADS 1980A&A....86...95H
CL90 = 1.645
# Table 1, air, >= 3000 A: (wavelength A, EP eV (2 dp), J lower, log gf, 90% abs error dex).
HUBER = [
    (3002.484, 0.03, 3, -0.06, 0.22), (3003.622, 0.11, 2, -0.33, 0.22), (3012.004, 0.42, 2, 0.00, 0.11),
    (3019.143, 0.00, 4, -1.34, 0.11), (3031.870, 0.00, 4, -1.83, 0.11), (3037.935, 0.03, 3, -0.51, 0.11),
    (3045.006, 0.16, 3, -1.77, 0.22), (3050.819, 0.03, 3, -0.13, 0.22), (3054.316, 0.11, 2, -0.62, 0.22),
    (3057.638, 0.21, 1, -0.37, 0.22), (3064.619, 0.11, 2, -1.13, 0.11), (3080.754, 0.21, 1, -1.18, 0.10),
    (3097.118, 0.16, 3, -1.47, 0.11), (3099.115, 0.16, 3, -1.58, 0.27), (3101.554, 0.11, 2, -0.14, 0.11),
    (3101.879, 0.42, 2, -0.31, 0.11), (3105.469, 0.27, 2, -1.53, 0.22), (3114.124, 0.11, 2, -1.52, 0.10),
    (3129.314, 0.27, 2, -2.41, 0.27), (3134.108, 0.21, 1, -0.28, 0.10), (3145.719, 0.16, 3, -2.13, 0.17),
    (3164.367, 0.27, 2, -1.99, 0.22), (3195.573, 0.27, 2, -2.36, 0.27), (3197.113, 0.21, 1, -1.74, 0.11),
    (3200.423, 0.03, 3, -2.61, 0.27), (3221.652, 0.00, 4, -1.90, 0.10), (3225.020, 0.42, 2, -1.30, 0.11),
    (3226.984, 0.00, 4, -2.50, 0.17), (3232.963, 0.00, 4, -1.04, 0.10), (3234.649, 0.11, 2, -1.66, 0.10),
    (3243.058, 0.03, 3, -1.27, 0.10), (3248.457, 0.03, 3, -2.17, 0.11), (3249.440, 0.27, 2, -2.40, 0.22),
    (3250.743, 0.42, 2, -1.70, 0.17), (3271.116, 0.11, 2, -2.24, 0.11), (3282.696, 0.16, 3, -2.30, 0.17),
    (3286.946, 0.03, 3, -2.50, 0.10), (3315.663, 0.11, 2, -1.21, 0.09), (3320.257, 0.16, 3, -1.40, 0.10),
    (3322.310, 0.42, 2, -1.30, 0.17), (3361.556, 0.11, 2, -1.42, 0.10), (3365.766, 0.42, 2, -1.24, 0.10),
    (3366.168, 0.16, 3, -1.40, 0.13), (3367.892, 0.03, 3, -2.65, 0.27), (3369.573, 0.00, 4, -0.69, 0.10),
    (3371.993, 0.16, 3, -1.40, 0.11), (3374.221, 0.03, 3, -1.76, 0.10), (3380.574, 0.42, 2, -0.22, 0.10),
    (3380.885, 0.27, 2, -1.34, 0.11), (3391.050, 0.00, 4, -1.07, 0.10), (3392.992, 0.03, 3, -0.54, 0.09),
    (3409.578, 0.00, 4, -2.41, 0.10), (3413.478, 0.16, 3, -1.48, 0.10), (3414.765, 0.03, 3, -0.06, 0.10),
    (3423.711, 0.21, 1, -0.73, 0.10), (3433.558, 0.03, 3, -0.67, 0.09), (3437.280, 0.00, 4, -1.19, 0.10),
    (3446.263, 0.11, 2, -0.41, 0.09), (3452.890, 0.11, 2, -0.91, 0.09), (3458.474, 0.21, 1, -0.26, 0.09),
    (3461.652, 0.03, 3, -0.36, 0.09), (3467.502, 0.16, 3, -2.00, 0.10), (3469.486, 0.27, 2, -1.80, 0.10),
    (3472.545, 0.11, 2, -0.81, 0.09), (3483.774, 0.27, 2, -1.11, 0.10), (3485.888, 0.21, 1, -2.20, 0.11),
    (3492.956, 0.11, 2, -0.27, 0.09), (3498.19, 0.00, 4, -4.84, 0.12), (3500.852, 0.16, 3, -1.29, 0.10),
    (3502.595, 0.00, 4, -2.76, 0.17), (3507.694, 0.16, 3, -2.46, 0.17), (3510.338, 0.21, 1, -0.67, 0.09),
    (3513.933, 0.21, 1, -2.01, 0.12), (3515.054, 0.11, 2, -0.24, 0.10), (3519.766, 0.27, 2, -1.42, 0.11),
    (3523.444, 0.03, 3, -2.37, 0.27), (3524.541, 0.03, 3, -0.03, 0.10), (3527.982, 0.16, 3, -2.25, 0.22),
    (3561.751, 0.00, 4, -2.69, 0.12), (3566.372, 0.42, 2, -0.27, 0.10), (3571.869, 0.16, 3, -1.16, 0.09),
    (3587.931, 0.03, 3, -2.38, 0.11), (3597.705, 0.21, 1, -1.08, 0.10), (3602.281, 0.16, 3, -2.15, 0.17),
    (3609.314, 0.11, 2, -2.13, 0.10), (3610.462, 0.11, 2, -1.15, 0.09), (3612.741, 0.27, 2, -1.39, 0.10),
    (3619.392, 0.42, 2, 0.00, 0.10), (3624.733, 0.00, 4, -2.48, 0.12), (3664.095, 0.27, 2, -1.95, 0.22),
    (3669.241, 0.16, 3, -2.28, 0.27), (3670.427, 0.16, 3, -2.23, 0.27), (3674.060, 0.03, 3, -2.57, 0.22),
    (3674.150, 0.42, 2, -2.07, 0.27), (3688.415, 0.27, 2, -2.11, 0.27), (3722.484, 0.21, 1, -2.21, 0.27),
    (3736.813, 0.42, 2, -1.84, 0.27), (3739.229, 0.16, 3, -2.35, 0.27), (3749.045, 0.00, 4, -3.40, 0.12),
    (3775.572, 0.42, 2, -1.35, 0.11), (3783.530, 0.42, 2, -1.31, 0.10), (3807.144, 0.42, 2, -1.18, 0.10),
    (3831.690, 0.42, 2, -2.00, 0.11), (3832.873, 0.16, 3, -3.66, 0.12), (3858.301, 0.42, 2, -0.97, 0.11),
    (3885.87, 0.27, 2, -4.56, 0.10), (3946.18, 0.16, 3, -5.47, 0.12), (3973.562, 0.42, 2, -2.06, 0.27),
    (4093.62, 0.16, 3, -5.69, 0.11),
]


def _transitions(vald: pd.DataFrame, w: float, ep: float, ep_tol: float, j_lo: float | None) -> pd.DataFrame:
    h = vald[((vald.wavelength_air_A - w).abs() <= VALD_TOL_A) & ((vald.ep_eV - ep).abs() <= ep_tol)]
    if j_lo is not None:
        h = h[(h.j_lower - j_lo).abs() < 0.01]
    return h.assign(_up=h.ep_upper_eV.round(3)).drop_duplicates(["ep_eV", "_up", "j_lower", "j_upper"])


def _row(paper, v, w, loggf, sig, doi):
    g_up = 2.0 * float(v.j_upper) + 1.0
    A = 10 ** loggf / (1.49919e-16 * w ** 2 * g_up) / 1e6
    return {"paper": paper, "species": "Ni I", "wavelength_air_A": round(w, 4), "ep_eV": float(v.ep_eV),
            "j_lower": float(v.j_lower), "ep_upper_eV": float(v.ep_upper_eV), "A_1e6_s": round(A, 4),
            "dA_1e6_s": round(A * (10 ** sig - 1), 4), "loggf": loggf, "loggf_sigma_dex": round(sig, 4),
            "doi": doi, "vizier": ""}


def build(vald: pd.DataFrame) -> dict[str, tuple[pd.DataFrame, dict]]:
    out = {}
    rows, c = [], {"transcribed": 0, "bluer_than_3000": 0, "no_single_vald_transition": 0}
    for nm, ep, gf in BLACKWELL:
        c["transcribed"] += 1
        w = nm * 10.0
        if w < LO_A:
            c["bluer_than_3000"] += 1
            continue
        t = _transitions(vald, w, ep, 0.002, None)
        if len(t) != 1:
            c["no_single_vald_transition"] += 1
            print(f"  blackwell1989 {w:.2f}: {len(t)} VALD transitions -- dropped")
            continue
        rows.append(_row("blackwell1989_ni1", t.iloc[0], w, gf, BLACKWELL_SIGMA, BLACKWELL_DOI))
    c["written"] = len(rows)
    out["blackwell1989_ni1"] = (pd.DataFrame(rows, columns=COLS), c)

    rows, c = [], {"transcribed": 0, "no_single_vald_transition": 0}
    for w, ep, jl, gf, e90 in HUBER:
        c["transcribed"] += 1
        t = _transitions(vald, w, ep, 0.011, jl)
        if len(t) != 1:
            c["no_single_vald_transition"] += 1
            print(f"  huber1980 {w:.3f}: {len(t)} VALD transitions -- dropped")
            continue
        rows.append(_row("huber1980_ni1", t.iloc[0], w, gf, e90 / CL90, HUBER_DOI))
    c["written"] = len(rows)
    out["huber1980_ni1"] = (pd.DataFrame(rows, columns=COLS), c)
    return out


def main() -> None:
    vald = pd.read_csv(VALD)
    for paper, (df, counts) in build(vald[vald.species == "Ni I"]).items():
        df.to_csv(OUTDIR / f"{paper}.csv", index=False)
        print(paper, counts)


if __name__ == "__main__":
    main()

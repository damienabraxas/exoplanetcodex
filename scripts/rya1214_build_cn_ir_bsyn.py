#!/usr/bin/env python3
"""RYA-1214 — the CN line list AGSS21's nitrogen actually needs: 1.087-1.108 um.

Ryan: *"visible N is going to be hard, hence the IR."* This is the IR.

WHY THIS LIST DOES NOT EXIST YET. Our vendored `12C14N_*.bsyn` set runs 4200-9200 A and
stops. AGSS21's nitrogen rests on CN at **10872-13204 A** (Amarsi et al. 2021 Table 2, all
59 CN lines in band "(0-0)" of the A-X system) and on NH at 2.9-15 um. So every CN line
list we hold is blueward of every CN line AGSS21 uses, and the optical CN A-X red system we
fitted in `CN_red` shares a molecule with theirs and nothing else — ZERO of their 59 lines
fall in our 6125-6130 / 6195-6200 A windows.

WHICH HOLDINGS REACH IT — MEASURED BY PROBING EACH ONE, NOT BY READING DECLARED SPANS.
Three holdings report no span at all (`span=None`, the reader inventories its own segments),
so a declared-span answer would have been wrong in both directions:

    solar_kpno_molecfit_corrected   CLEAN         the WHOLE band 10872-13204 A   all 59 lines
    solar_iag                       CLEAN         10872-11083 A (its red edge)   16 lines
    solar_kpno                      CONTROL_ONLY  same reach — not a science basis (RYA-1026)
    CRIRES+ Y / Y-wide              CLEAN         end at 10796 A — 76 A SHORT of the band
    CRIRES+ H                       CLEAN         starts at 15007 A — 1803 A PAST it
    solar_vesta_crires_plus_idp     BLOCKED       RestFrameNotConditioned

🔴 SO CRIRES+ CANNOT MEASURE AGSS21's NITROGEN, AND THAT IS A MEASUREMENT, NOT A SKIP. Its Y
arm stops 76 A blueward of the band's first line and its H arm starts 1803 A redward of its
last; the raw Vesta IDPs reach further but are telluric-uncorrected AND rest-frame
unconditioned, so they are BLOCKED rather than merely uncorrected. The instrument catalogue
says CRIRES+ spans 9500-53000 A, which is why this had to be probed per HOLDING: the
instrument's reach and the holdings' reach differ by the whole band.

⚠️ ONE TELLURIC BAND OVERLAPS: H2O at 11120-11560 A (`telluric_policy.TELLURIC_BANDS`). Fit
windows AVOID it rather than lean on a correction — it is 440 A of the band's 2333, leaving
1893 A of enumerated-clean spectrum, so there is nothing to gain by fitting inside it.

THE SOURCE IS PRIMARY AND THE PARSE IS VALIDATED, NOT ASSUMED
-------------------------------------------------------------
Brooke et al. 2014 (ApJS 210, 23; CDS J/ApJS/210/23) Table 4 — 195,120 CN lines across the
A-X, B-X and X-X systems — is vendored at
`data/reference/cno_molecular_primary/cn_brooke2014/table4.dat.gz`. It yields **11,290 lines
in 10872-13205 A**, 1,092 of them in the A-X (0-0) band AGSS21 uses.

Read by BYTE POSITION from the ReadMe's own byte-by-byte block, never by whitespace
splitting: the fields are fixed-width and `Obs` (39-49) is nullable, so a split-on-space
parse shifts every column after it on those rows. A first pass did exactly that and put the
wavelength span at 862069-100000000 A, i.e. it was reading the O-C residual as a wavenumber.

🔴 E'' IS THE TOTAL LOWER-STATE ENERGY, AND THAT HAD TO BE CHECKED. The ReadMe calls bytes
71-80 "Lower state energy relative to v''=0", which can be read as rotational-only. It is
not: measured on the data, the minimum E'' per lower vibrational level is 0.00 (v''=0),
2042.42 (v''=1), 8011.77 (v''=4), 11859.33 (v''=6) and 21013.29 cm-1 (v''=11) — i.e. it
steps by CN's omega_e ~ 2068 cm-1 and therefore INCLUDES the vibrational term. Had it been
rotational-only, every v''>0 line would have carried an excitation potential too low by up
to 2.6 eV.

THE VALIDATION GATE. `log gf = log10((2J'+1) * f)` is checked against Amarsi 2021 Table 2's
own published `log gf` on his 59 CN lines before anything is written: 58 of 59 match within
0.05 A, median |d log gf| = 0.0077 dex. A UNIQUE wavelength match is required — 13 of the
58 sit in windows with 2-3 candidates and nearest-lambda picks wrong there, which is what
produced a 1.99 dex outlier in the first pass. Ambiguity is refused, not resolved.

The `.bsyn` writing convention (species code, `gu = 2J+1`, fdamp and gamma_rad left 0.0 so
Turbospectrum applies its own defaults) is taken from `pipeline.nearuv_linelist.
write_molecular_bsyn` — RYA-1207's writer, whose output is the validated near-UV molecular
opacity. ⚠️ It is NOT taken from the vendored files: `12C14N_900-950.bsyn` writes 3.5 in
column 5 where `CO_IR_Li2015.dat` writes 19.0 for J=9, so the two disagree about whether
that column is J or 2J+1, and following the wrong one scales every line strength.
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SRC = ROOT / "data/reference/cno_molecular_primary/cn_brooke2014/table4.dat.gz"
AMARSI = ROOT / "data/reference/amarsi2021_cno/derived/amarsi2021_cno_molecular_lines.csv"
OUT_DIR = ROOT / "data/linelists/molecular/turbospectrum/CN"
AUDIT = ROOT / "data/audit/rya1214_cno_products"

#: AGSS21's CN band, END TO END. The list covers the whole thing; the REGION's fit windows
#: decide what any one holding uses, which is the right split — a line list is not a claim
#: about coverage, and writing one list per holding would be three copies of one parse.
#:
#: Measured coverage per holding (`allow_uncorrected` probe, not declared spans):
#:   solar_kpno_molecfit_corrected  CLEAN         the WHOLE band, 10872-13204 A
#:   solar_iag                      CLEAN         10872-11083 only (its own red edge)
#:   solar_kpno                     CONTROL_ONLY  same reach, not a science basis (RYA-1026)
#:   CRIRES+ Y / Y-wide             CLEAN         end at 10796 A — 76 A SHORT of the band
#:   CRIRES+ H                      CLEAN         starts at 15007 A — 1803 A PAST it
#:   solar_vesta_crires_plus_idp    BLOCKED       RestFrameNotConditioned
LO_A, HI_A = 10872.0, 13205.0
#: The ONE telluric band inside it — `telluric_policy.TELLURIC_BANDS`. Fit windows avoid it
#: rather than relying on a correction: the band is 440 A wide and the rest of the CN band
#: is 1893 A of enumerated-clean spectrum, so there is nothing to gain by fitting inside it.
H2O_BAND_A = (11120.0, 11560.0)
#: `12C14N`, the dominant isotopologue — the same code RYA-1207's writer uses for CN.
SPECIES_CODE = "0607.012014"
#: Validation thresholds. The match must be UNIQUE, and the gf must reproduce Amarsi's.
MATCH_TOL_A = 0.05
GF_TOL_DEX = 0.02
MIN_VALIDATED = 20


def parse_brooke() -> pd.DataFrame:
    """Byte-position parse of Brooke 2014 Table 4. See the module docstring for why."""
    rows = []
    with gzip.open(SRC, "rt") as fh:
        for ln in fh:
            if len(ln.rstrip("\n")) < 106:
                continue
            try:
                eu, el = ln[0], ln[2]
                vu, vl = int(ln[4:6]), int(ln[7:9])
                ju, jl = float(ln[10:15]), float(ln[16:21])
                obs, cal = ln[38:49].strip(), ln[50:60].strip()
                e_low = float(ln[70:80])
                a_ki = float(ln[81:93])
                f_val = float(ln[94:106])
            except ValueError:
                continue
            wn = float(obs) if obs else (float(cal) if cal else 0.0)
            if wn <= 0:
                continue
            rows.append((eu, el, vu, vl, ju, jl, wn, e_low, a_ki, f_val))
    d = pd.DataFrame(rows, columns=["eu", "el", "vu", "vl", "Ju", "Jl", "wn",
                                    "Ecm", "A", "f"])
    d["lam_vac"] = 1e8 / d.wn
    s2 = (1e4 / d.lam_vac) ** 2
    d["lam_air"] = d.lam_vac / (1 + 0.0000834254 + 0.02406147 / (130 - s2)
                                + 0.00015998 / (38.9 - s2))
    d["gu"] = 2.0 * d.Ju + 1.0
    d["loggf"] = np.log10(d.gu * d.f)
    d["Elow_eV"] = d.Ecm / 8065.54429
    return d


def assert_e_includes_vibration(d: pd.DataFrame) -> dict:
    """E'' must step by ~omega_e per lower vibrational level, or it is rotational-only."""
    mins = {int(v): round(float(g.Ecm.min()), 2) for v, g in d.groupby("vl") if v <= 11}
    steps = [mins[v] - mins[v - 1] for v in sorted(mins) if v >= 1 and (v - 1) in mins]
    med = float(np.median(steps)) if steps else 0.0
    if not (1500.0 < med < 2600.0):
        raise SystemExit(
            f"E'' does not step by CN's omega_e (~2068 cm-1): median step {med:.1f}. It is "
            f"either rotational-only or a different column — refusing to write excitation "
            f"potentials from it.")
    return {"min_Ecm_by_vlow": mins, "median_step_cm1": round(med, 1),
            "verdict": "E'' INCLUDES the vibrational term; Elow_eV = E''/8065.54 is the "
                       "total lower-state energy"}


def validate_against_amarsi(d: pd.DataFrame) -> dict:
    """Reproduce Amarsi 2021 Table 2's own published CN log gf, or refuse to write."""
    a = pd.read_csv(AMARSI)
    cn = a[(a.element_parameter == "logepsN") & (a.species == "CN")].copy()
    cn["lam_vac"] = cn.wavelength_vac_nm * 10.0
    hits, amb, miss = [], 0, 0
    for _, r in cn.iterrows():
        w = float(r.lam_vac)
        m = d[(d.lam_vac - w).abs() <= MATCH_TOL_A]
        if m.empty:
            miss += 1
            continue
        if len(m) > 1:
            amb += 1                      # refused, never resolved by proximity
            continue
        hits.append(float(m.iloc[0].loggf) - float(r.published_loggf))
    hits = np.asarray(hits, float)
    n_ok = int((np.abs(hits) <= GF_TOL_DEX).sum())
    rep = {"amarsi_cn_lines": int(len(cn)), "uniquely_matched": int(len(hits)),
           "ambiguous_refused": amb, "unmatched": miss,
           "median_abs_dloggf": round(float(np.median(np.abs(hits))), 4) if len(hits) else None,
           "max_abs_dloggf": round(float(np.abs(hits).max()), 4) if len(hits) else None,
           "within_0p02_dex": n_ok,
           "match_tol_A": MATCH_TOL_A}
    if len(hits) < MIN_VALIDATED or n_ok < 0.8 * len(hits):
        raise SystemExit(
            f"gf validation FAILED: {n_ok}/{len(hits)} within {GF_TOL_DEX} dex of Amarsi's "
            f"published values (need >=80% of >={MIN_VALIDATED}). Refusing to write a line "
            f"list whose gf does not reproduce the source AGSS21 used.")
    return rep


def main() -> int:
    AUDIT.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    d = parse_brooke()
    print(f"parsed {len(d)} Brooke 2014 CN lines, {d.lam_air.min():.0f}-{d.lam_air.max():.0f} A")
    e_rep = assert_e_includes_vibration(d)
    print(f"  E'' check: median step {e_rep['median_step_cm1']} cm-1 — includes vibration")
    v_rep = validate_against_amarsi(d)
    print(f"  gf validation: {v_rep['within_0p02_dex']}/{v_rep['uniquely_matched']} within "
          f"{GF_TOL_DEX} dex of Amarsi 2021 (median |d| {v_rep['median_abs_dloggf']}), "
          f"{v_rep['ambiguous_refused']} ambiguous REFUSED")

    band = d[(d.lam_air >= LO_A) & (d.lam_air <= HI_A)].sort_values("lam_air")
    if band.empty:
        raise SystemExit(f"no CN line in {LO_A}-{HI_A} A — refusing to write an empty list")
    fn = OUT_DIR / f"12C14N_{int(LO_A//10)}-{int(HI_A//10)}.bsyn"
    with fn.open("w") as fh:
        fh.write(f"'         {SPECIES_CODE}'    1   {len(band)}\n")
        fh.write("'Brooke+2014 ApJS 210 23 CN A-X/B-X/X-X (RYA-1214); 12C14N'\n")
        for r in band.itertuples():
            fh.write(f"{r.lam_air:12.3f} {r.Elow_eV:9.5f} {r.loggf:7.3f}    0.000 "
                     f"{r.gu:6.1f}  0.00E+00 'X' 'X'   0.0    0.0\n")
    sysc = band.groupby([band.eu + "-" + band.el, band.vu.astype(str) + "-"
                         + band.vl.astype(str)]).size().nlargest(5)
    print(f"\nwrote {fn.relative_to(ROOT)}: {len(band)} lines, "
          f"{band.lam_air.min():.2f}-{band.lam_air.max():.2f} A")
    for k, v in sysc.items():
        print(f"    {k[0]} ({k[1]}) {v}")

    (AUDIT / "cn_ir_bsyn_build.prov.json").write_text(json.dumps({
        "ticket": "RYA-1214",
        "why": "AGSS21's nitrogen rests on CN at 10872-13204 A; every CN list we held runs "
               "4200-9200 A and stops. The optical CN_red band we fitted shares a molecule "
               "with AGSS21's and nothing else — 0 of their 59 lines are in our windows.",
        "source": "Brooke, Ram, Western, Li, Schwenke & Bernath 2014, ApJS 210, 23 "
                  "(CDS J/ApJS/210/23) Table 4 — 195,120 CN lines, A-X / B-X / X-X",
        "source_file": str(SRC.relative_to(ROOT)),
        "parse": "BYTE POSITION from the ReadMe's byte-by-byte block. `Obs` (39-49) is "
                 "nullable, so whitespace splitting shifts every later column on those "
                 "rows — a first pass did that and read the O-C residual as a wavenumber.",
        "e_low_check": e_rep,
        "gf_validation": v_rep,
        "window_A": [LO_A, HI_A],
        "window_basis": "AGSS21's CN band end to end. MEASURED coverage per holding, by "
                        "probing each one rather than reading declared spans: "
                        "solar_kpno_molecfit_corrected (CLEAN) reaches the WHOLE band; "
                        "solar_iag (CLEAN) reaches 10872-11083 only; solar_kpno is "
                        "CONTROL_ONLY; CRIRES+ Y ends at 10796 A (76 A short) and H starts "
                        "at 15007 A (1803 A past), so CRIRES+ CANNOT serve this band; the "
                        "raw Vesta IDPs are BLOCKED (RestFrameNotConditioned).",
        "telluric": "one band overlaps — H2O 11120-11560 A (telluric_policy."
                    "TELLURIC_BANDS). Fit windows AVOID it rather than relying on a "
                    "correction: it is 440 A of the band's 2333, and the remaining 1893 A "
                    "is enumerated-clean.",
        "h2o_band_A": list(H2O_BAND_A),
        "n_lines_written": int(len(band)),
        "bsyn_convention": "pipeline.nearuv_linelist.write_molecular_bsyn (RYA-1207): "
                           "gu = 2J+1, fdamp and gamma_rad 0.0 for Turbospectrum defaults. "
                           "NOT copied from the vendored files, which disagree with each "
                           "other about whether column 5 is J or 2J+1.",
    }, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

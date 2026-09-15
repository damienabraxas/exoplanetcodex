#!/usr/bin/env python3
"""RYA-1214 — the OH line list for AGSS21's oxygen indicators in the CRIRES+ H band.

AGSS21's OH vibration-rotation indicators include 14 lines at 15283-16914 A (vacuum), in the
(2-0), (3-1) and (4-2) bands, which sit inside the CRIRES+ H holding (15007-17494 A, CLEAN,
RYA-1191). No synthesis-ready OH list covers them.

🔴 THE VENDORED MYTHOS FILE IS NOT USABLE AS-IS, AND THAT WAS MEASURED.
`data/linelists/molecular/turbospectrum/OH/16O-1H__MYTHOS_rovib.bsyn` carries 4,188 lines in
the band, but against AGSS21's own 14 lines:
  * its wavelengths are VACUUM (|d lambda| <= 0.022 A in the vacuum frame, up to 1.5 A in
    air), and Turbospectrum reads AIR;
  * its log gf is +0.30 to +0.36 dex too strong on all 14, because its `gu` is 2(2J+1) —
    ExoMol folds in the H nuclear-spin degeneracy (J=8.5 -> gu 36) — so every line would be
    ~2x too strong.

THE SOURCE IS THE ONE AGSS21 USED. Brooke et al. 2016 (JQSRT 168, 142) X-X line list,
vendored in `data/reference/cno_molecular_primary/oh_brooke2016/OH-Supplementary.zip`. Its
transition labels are AGSS21's labels verbatim (`pP1e(9.5)`).

THE GATES, BOTH OF WHICH REFUSE TO WRITE:
  * E'' must step by OH's vibrational quantum (~3570 cm-1) per lower level, i.e. include the
    vibrational term. Measured: 0.0 / 3568.5 / 6971.3 / 10210.5 cm-1 for v''=0..3.
  * log gf = log10((2J''+1) f) must reproduce AGSS21's published values on the 14 H lines,
    matched by LABEL and wavelength together, never by proximity alone. Measured median
    |d| 0.0004 dex. The upper-state g gives 0.041 dex and is refused.

Positions: the OBSERVED wavenumber where Brooke has one, else PGOPHER's calculated one;
vacuum -> air through `pipeline.wavelength_util` (RYA-501). Written in RYA-1207's `.bsyn`
convention (gu = 2J'+1, damping 0 for Turbospectrum defaults), as the IR CN list was.
"""
from __future__ import annotations

import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.wavelength_util import vac_to_air  # noqa: E402  RYA-501 single source

SRC = ROOT / "data/reference/cno_molecular_primary/oh_brooke2016/OH-Supplementary.zip"
MEMBER = "OH-Supplementary/OH-XX-Line_list.txt"
AGSS21 = ROOT / "data/audit/rya1136_cno_intake/molecular_physical_crossmatch.csv"
OUT_DIR = ROOT / "data/linelists/molecular/turbospectrum/OH"
AUDIT = ROOT / "data/audit/rya1214_cno_products"

#: The CRIRES+ H holding's band (config.synth_bands 'H'), air Angstrom.
LO_A, HI_A = 15007.11, 17493.69
#: 16O-1H — the code RYA-1207's writer uses for OH.
SPECIES_CODE = "0108.000016"
MATCH_TOL_A = 0.05
GF_TOL_DEX = 0.02
MIN_VALIDATED = 10


def parse_brooke() -> pd.DataFrame:
    """Whitespace parse is safe here: absent fields are written as `...`, never blank."""
    lines = zipfile.ZipFile(SRC).read(MEMBER).decode("latin-1").splitlines()
    head = next(i for i, ln in enumerate(lines) if ln.strip().startswith("v' v''"))
    rows = []
    for ln in lines[head + 1:]:
        p = ln.split()
        if len(p) != 17:
            continue
        obs, calc = p[8], p[9]
        wn = float(obs) if obs != "..." else float(calc)
        rows.append(dict(vu=int(p[0]), vl=int(p[1]), Ju=float(p[2]), Jl=float(p[3]),
                         wn=wn, observed=obs != "...", Ecm=float(p[11]), f=float(p[13]),
                         label=p[14]))
    d = pd.DataFrame(rows)
    d = d[d.wn > 0].copy()
    d["lam_vac"] = 1e8 / d.wn
    d["lam_air"] = vac_to_air(d.lam_vac.to_numpy(float))
    d["loggf"] = np.log10((2.0 * d.Jl + 1.0) * d.f)
    d["gu"] = 2.0 * d.Ju + 1.0
    d["Elow_eV"] = d.Ecm / 8065.54429
    return d


def assert_e_includes_vibration(d: pd.DataFrame) -> dict:
    mins = {int(v): round(float(g.Ecm.min()), 1) for v, g in d.groupby("vl") if v <= 5}
    steps = [mins[v] - mins[v - 1] for v in sorted(mins) if v >= 1 and (v - 1) in mins]
    med = float(np.median(steps)) if steps else 0.0
    if not (2800.0 < med < 3800.0):
        raise SystemExit(f"E'' does not step by OH's vibrational quantum (~3570 cm-1): "
                         f"median {med:.1f}. Refusing to write excitation potentials.")
    return {"min_Ecm_by_vlow": mins, "median_step_cm1": round(med, 1)}


def validate_against_agss21(d: pd.DataFrame) -> dict:
    m = pd.read_csv(AGSS21)
    oh = m[(m.species == "OH") & (m.wavelength_vac_nm * 10 >= 15000)
           & (m.wavelength_vac_nm * 10 <= 17600)]
    dgf, dep, amb, miss = [], [], 0, 0
    for _, r in oh.iterrows():
        c = d[(d.label == str(r.raw_transition_label).strip())
              & ((d.lam_vac - r.wavelength_vac_nm * 10).abs() <= MATCH_TOL_A)]
        if c.empty:
            miss += 1
            continue
        if len(c) > 1:
            amb += 1
            continue
        dgf.append(float(c.iloc[0].loggf) - float(r.published_loggf))
        dep.append(float(c.iloc[0].Elow_eV) - float(r.lower_energy_eV))
    dgf = np.asarray(dgf)
    n_ok = int((np.abs(dgf) <= GF_TOL_DEX).sum())
    rep = {"agss21_oh_lines_in_band": int(len(oh)), "matched_by_label_and_wavelength":
           int(len(dgf)), "ambiguous_refused": amb, "unmatched": miss,
           "median_abs_dloggf": round(float(np.median(np.abs(dgf))), 4) if len(dgf) else None,
           "max_abs_dloggf": round(float(np.abs(dgf).max()), 4) if len(dgf) else None,
           "median_abs_dEP_eV": round(float(np.median(np.abs(dep))), 4) if dep else None,
           "within_tol": n_ok, "gf_tol_dex": GF_TOL_DEX, "match_tol_A": MATCH_TOL_A}
    if len(dgf) < MIN_VALIDATED or n_ok < 0.8 * len(dgf):
        raise SystemExit(f"gf validation FAILED: {n_ok}/{len(dgf)} within {GF_TOL_DEX} dex "
                         f"of AGSS21. Refusing to write.")
    return rep


def main() -> int:
    AUDIT.mkdir(parents=True, exist_ok=True)
    d = parse_brooke()
    e_rep = assert_e_includes_vibration(d)
    v_rep = validate_against_agss21(d)
    print(f"parsed {len(d)} Brooke 2016 OH X-X lines; E'' median step "
          f"{e_rep['median_step_cm1']} cm-1; gf {v_rep['within_tol']}/"
          f"{v_rep['matched_by_label_and_wavelength']} within {GF_TOL_DEX} dex of AGSS21 "
          f"(median |d| {v_rep['median_abs_dloggf']})")
    band = d[(d.lam_air >= LO_A) & (d.lam_air <= HI_A)].sort_values("lam_air")
    if band.empty:
        raise SystemExit("no OH line in the H band — refusing to write an empty list")
    fn = OUT_DIR / f"16OH_{int(LO_A // 10)}-{int(HI_A // 10)}.bsyn"
    with fn.open("w") as fh:
        fh.write(f"'         {SPECIES_CODE}'    1   {len(band)}\n")
        fh.write("'Brooke+2016 JQSRT 168 142 OH X-X (RYA-1214); 16OH'\n")
        for r in band.itertuples():
            fh.write(f"{r.lam_air:12.3f} {r.Elow_eV:9.5f} {r.loggf:7.3f}    0.000 "
                     f"{r.gu:6.1f}  0.00E+00 'X' 'X'   0.0    0.0\n")
    print(f"wrote {fn.relative_to(ROOT)}: {len(band)} lines, "
          f"{band.lam_air.min():.2f}-{band.lam_air.max():.2f} A "
          f"({int(band.observed.sum())} at observed positions)")
    (AUDIT / "oh_h_bsyn_build.prov.json").write_text(json.dumps({
        "ticket": "RYA-1214",
        "source": "Brooke, Bernath, Western et al. 2016, JQSRT 168, 142 — OH X-X line list",
        "source_file": f"{SRC.relative_to(ROOT)}::{MEMBER}",
        "why_not_mythos": "vendored 16O-1H__MYTHOS_rovib.bsyn is in VACUUM and its gu "
                          "includes nuclear-spin degeneracy: +0.30..+0.36 dex vs AGSS21 on "
                          "all 14 H-band lines",
        "position": "observed wavenumber where available, else PGOPHER calculated; "
                    "vacuum -> air via pipeline.wavelength_util",
        "loggf": "log10((2J''+1) f) — validated; the upper-state g is refused (0.041 dex)",
        "e_low_check": e_rep,
        "gf_validation": v_rep,
        "window_A": [LO_A, HI_A],
        "n_lines_written": int(len(band)),
        "n_at_observed_positions": int(band.observed.sum()),
    }, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

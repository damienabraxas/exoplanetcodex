#!/usr/bin/env python3
"""Vendor the primary-laboratory Al I log gf measurements — RYA-1002.

    python3 scripts/rya1002_fetch_al_gf_lab.py

WHY THIS ONE, AND WHY PRIMARY
-----------------------------
Al had NO lab-gf table in this repo, and `pipeline.gf_grades.GF-LAB` was hard-wired to
the Fe file, so no Al line could reach the graded state however good its provenance.

  Burheim, Hartman & Nilsson 2023, A&A 672 A197 (arXiv:2309.06273)
      EXPERIMENTAL Al I oscillator strengths: branching fractions from Fourier-transform
      spectrometry of a hollow-cathode lamp, combined with radiative lifetimes. 12 lines,
      670-4200 nm, 2-11 % accuracy, each with the paper's own per-line Unc_gf. That
      per-line uncertainty is the quantity a compilation tag cannot give.

Deliberately NOT treated as independent sources, and NOT in this table:
  * Nordlander & Lind 2017 — its Al uncertainties trace to NIST/Kelleher-Podobedova 2008,
    so counting it would double-count the NIST rung (RYA-835).
  * `1995JPhB..` — RYA-1001 identified it as Mendoza, Eissner, Le Dourneuf & Zeippen 1995,
    J.Phys.B 28 3485, the Opacity Project close-coupling calculation. THEORY, not a lab
    measurement, despite the bibcode reading like one.
  * Papoulia 2019 — independent THEORY (MCDHF). Carried as a comparison column only.

⚠️ THE COLUMN TRAP THIS SCRIPT EXISTS TO NOT REPEAT
---------------------------------------------------
Burheim's tables print `sigma [cm^-1]` immediately beside `lambda_vac [A]`, and over this
range the two are the same order of magnitude. RYA-835 read the union of both columns as
a wavelength list; every "nearest Burheim line" it then computed was meaningless, and it
concluded the paper covered none of our lines when it covers eight. This script parses the
columns positionally out of the LaTeX and ASSERTS `sigma * lambda_vac == 1e8` row by row,
which is the identity that makes confusing them impossible.

WHAT IS EMITTED
---------------
`data/reference/al_gf_lab/al1_lab_loggf.csv` — one row per measured line, in the SAME
schema as `fe_gf_lab/fe1_lab_loggf.csv` so `gf_grades` consumes both with no special case,
plus the provenance JSON. All 12 lines are carried, not only the 8 inside our instrument
reach: a subset is not a line list, and the out-of-reach four are the control that the
transcription is complete.

CONTROLS, asserted — the script refuses to write if any fails:
  * exactly 12 derived-log gf rows (the count the paper states);
  * lambda_vac spans the paper's stated 670-4200 nm;
  * the sigma and lambda column SETS are disjoint, and sigma*lambda == 1e8 per row;
  * every row carries a finite uncertainty in dex — the entire point of using this source;
  * levels come from NIST ASD, an INDEPENDENT source, and `Ek - Ei` there reproduces
    Burheim's sigma. That is what proves the level columns and the wavelength column
    describe the SAME transition rather than having been mis-associated while parsing.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
import tarfile
import urllib.request
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.wavelength_util import vac_to_air        # noqa: E402  SSOT converter

ARXIV_ID = "2309.06273"
EPRINT = f"https://arxiv.org/e-print/{ARXIV_ID}"
CITATION = "Burheim, Hartman & Nilsson 2023, A&A 672, A197"
DOI = "10.1051/0004-6361/202245394"
SOURCE_TAG = "Burheim2023"

NIST_TSV = ROOT / "data/linelists/primary_gf/nist_asd_AlI_6600_42000.tsv"
OUT_DIR = ROOT / "data" / "reference" / "al_gf_lab"
OUT_CSV = OUT_DIR / "al1_lab_loggf.csv"
OUT_PROV = OUT_DIR / "al1_lab_loggf.prov.json"

N_EXPECTED = 12                 # the paper's own count
CM1_TO_EV = 1.0 / 8065.543937   # CODATA
#: How closely NIST's (Ek - Ei) must reproduce Burheim's sigma for the two sources to be
#: describing the same transition. 0.5 cm^-1 is far tighter than any level spacing here
#: (the closest distinct lower levels, 3d 2D3/2 and 2D5/2, are 1.34 cm^-1 apart).
SIGMA_MATCH_CM1 = 0.5
NIST_WAVE_TOL_A = 0.05


def fetch_source(cache: Path | None = None) -> tuple[str, str]:
    """(LaTeX text, md5 of the tarball). The e-print, never a rendered PDF — the table is
    the artifact and a PDF's column order is a rendering of it, not the thing itself."""
    if cache and cache.exists():
        blob = cache.read_bytes()
    else:
        req = urllib.request.Request(EPRINT, headers={"User-Agent": "exoplanetcodex/1.0"})
        with urllib.request.urlopen(req, timeout=120) as fh:
            blob = fh.read()
        if cache:
            cache.write_bytes(blob)
    md5 = hashlib.md5(blob).hexdigest()
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tf:
        texs = [m for m in tf.getmembers() if m.name.endswith(".tex")]
        if len(texs) != 1:
            raise SystemExit(f"expected exactly one .tex in the e-print, got "
                             f"{[t.name for t in texs]}")
        return tf.extractfile(texs[0]).read().decode("utf8", "replace"), md5


def _cells(line: str) -> list[str]:
    return [c.strip() for c in line.rstrip("\\ ").split("&")]


def _num(x: str) -> float:
    return float(x.replace("$", "").replace("{", "").replace("}", "").strip())


def parse_table3(tex: str) -> pd.DataFrame:
    """Table 3 `tab:loggf_comp` — the lines with a DERIVED log gf. This is the table that
    grades a line. Table 2 `tab:BRtable` carries branching fractions with NO log gf, so a
    line appearing only there is NOT graded by Burheim."""
    i = tex.index("\\label{tab:loggf_comp}")
    body = tex[i:tex.index("\\end{table}", i)]
    rows = []
    for ln in body.splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("%") or "\\hline" in ln or "&" not in ln:
            continue
        c = _cells(ln)
        if len(c) != 8:
            continue
        try:
            rows.append(dict(lower_level=c[0], upper_level=c[1], sigma_cm1=_num(c[2]),
                             lam_vac_A=_num(c[3]), loggf=_num(c[4]),
                             loggf_papoulia19=_num(c[5]), loggf_kurucz95=_num(c[6]),
                             loggf_topbase00=_num(c[7])))
        except ValueError:
            continue                      # a header row, not data
    return pd.DataFrame(rows)


def parse_table1_uncertainties(tex: str) -> dict[float, float]:
    """Table 1 `tab:BFtable` is where the per-line `Unc_gf` [%] lives."""
    i = tex.index("\\label{tab:BFtable}")
    body = tex[i:tex.index("\\end{table}", i)]
    unc: dict[float, float] = {}
    for ln in body.splitlines():
        ln = ln.strip()
        if not ln or "\\hline" in ln or "&" not in ln or "Residual" in ln:
            continue
        c = _cells(ln)
        if len(c) != 9:
            continue
        try:
            lam = _num(c[3])
            unc[round(lam, 3)] = _num(c[8].replace("$*$", "").replace("*", ""))
        except ValueError:
            continue
    return unc


def nist_levels() -> pd.DataFrame:
    """Ei/Ek per line from NIST ASD — an INDEPENDENT source for the levels, which is what
    makes the sigma control below a real test rather than an identity."""
    if not NIST_TSV.exists():
        raise SystemExit(
            f"NIST pull missing at {NIST_TSV}. Regenerate with:\n"
            f"  venv_ci/bin/python scripts/rya822_pull_nist_nearuv.py "
            f"--species 'Al I' --lo-A 6600 --hi-A 42000 --step-A 500\n"
            f"Refusing to continue: deriving the levels from Burheim's own sigma would "
            f"make the consistency check circular.")
    t = pd.read_csv(NIST_TSV, sep="\t")
    raw = t["ei_ek_raw"].astype(str).str.replace(r"[\[\]a-zA-Z+?]", "", regex=True)
    parts = raw.str.split("-", n=1, expand=True)
    t["ei_eV_p"] = pd.to_numeric(parts[0], errors="coerce")
    t["ek_eV_p"] = pd.to_numeric(parts[1], errors="coerce")
    t["wl"] = pd.to_numeric(t["wavelength_A"], errors="coerce")
    return t.dropna(subset=["wl", "ei_eV_p", "ek_eV_p"])


def build(cache: Path | None) -> tuple[pd.DataFrame, dict]:
    tex, md5 = fetch_source(cache)
    t3 = parse_table3(tex)
    unc = parse_table1_uncertainties(tex)

    # ── CONTROL 1: the count the paper states ────────────────────────────────
    assert len(t3) == N_EXPECTED, (
        f"expected {N_EXPECTED} derived-log gf rows in tab:loggf_comp, parsed {len(t3)}. "
        f"The paper states 12; a different count means the parse drifted.")
    # ── CONTROL 2: the stated range, 670-4200 nm ─────────────────────────────
    assert 6600 < t3.lam_vac_A.min() and t3.lam_vac_A.max() < 42000, (
        f"lambda_vac {t3.lam_vac_A.min():.1f}-{t3.lam_vac_A.max():.1f} A is outside the "
        f"paper's stated 670-4200 nm")
    # ── CONTROL 3: THE COLUMN TRAP (RYA-835) ─────────────────────────────────
    overlap = set(t3.sigma_cm1.round(3)) & set(t3.lam_vac_A.round(3))
    assert not overlap, f"sigma[cm-1] and lambda_vac[A] column sets overlap: {overlap}"
    prod = t3.sigma_cm1 * t3.lam_vac_A
    assert np.allclose(prod, 1e8, rtol=1e-4), (
        f"sigma * lambda != 1e8 — the two columns are NOT the same transitions in two "
        f"units, so the parse mis-associated them: {prod.tolist()}")

    t3["unc_pct"] = t3.lam_vac_A.round(3).map(unc)
    # ── CONTROL 4: every row has a finite uncertainty ────────────────────────
    missing = t3[t3.unc_pct.isna()].lam_vac_A.tolist()
    assert not missing, (
        f"no Unc_gf found in tab:BFtable for {missing} — a lab table whose whole value is "
        f"the per-line sigma must not carry a blank one")
    t3["e_loggf_dex"] = np.log10(1.0 + t3.unc_pct / 100.0)

    t3["wavelength_air_A"] = vac_to_air(t3.lam_vac_A.values)

    # ── CONTROL 5: levels from NIST, and Ek-Ei must reproduce Burheim's sigma ─
    nl = nist_levels()
    elo, eup, matched = [], [], 0
    for _, r in t3.iterrows():
        cand = nl[(nl.wl - r.wavelength_air_A).abs() <= NIST_WAVE_TOL_A]
        hit = None
        for _, c in cand.iterrows():
            dsig = abs((c.ek_eV_p - c.ei_eV_p) / CM1_TO_EV - r.sigma_cm1)
            if dsig <= SIGMA_MATCH_CM1:
                hit = c
                break
        if hit is None:
            elo.append(np.nan); eup.append(np.nan)
        else:
            elo.append(float(hit.ei_eV_p)); eup.append(float(hit.ek_eV_p)); matched += 1
    t3["elo_eV"], t3["eup_eV"] = elo, eup
    assert matched >= 8, (
        f"only {matched}/12 Burheim lines could be tied to a NIST level pair whose "
        f"Ek-Ei reproduces the paper's sigma to {SIGMA_MATCH_CM1} cm^-1. Below 8 the "
        f"level assignment is not corroborated and the table must not be written.")
    t3["elo_cm1"] = t3.elo_eV / CM1_TO_EV
    t3["eup_cm1"] = t3.eup_eV / CM1_TO_EV
    t3["source"] = SOURCE_TAG

    cols = ["source", "wavelength_air_A", "elo_cm1", "eup_cm1", "elo_eV", "eup_eV",
            "loggf", "e_loggf_dex",
            # Al-specific provenance beyond the Fe schema. Additive: gf_grades reads the
            # shared columns by name and ignores these.
            "lam_vac_A", "sigma_cm1", "unc_pct", "lower_level", "upper_level",
            "loggf_papoulia19", "loggf_kurucz95", "loggf_topbase00"]
    out = t3[cols].sort_values("wavelength_air_A").reset_index(drop=True)

    prov = dict(
        ticket="RYA-1002", source=CITATION, doi=DOI, arxiv=ARXIV_ID,
        eprint_url=EPRINT, eprint_md5=md5, table="Table 3 (tab:loggf_comp)",
        uncertainties_from="Table 1 (tab:BFtable), column Unc_gf [%]",
        n_lines=int(len(out)), n_with_nist_levels=int(matched),
        pulled_utc=date.today().isoformat(),
        wavelength_frame="lambda_vac as published; wavelength_air_A via "
                         "pipeline.wavelength_util.vac_to_air (Birch & Downs 1994, SSOT)",
        levels_source=f"NIST ASD via {NIST_TSV.name} — independent of Burheim; "
                      f"Ek-Ei reproduces the paper's sigma to {SIGMA_MATCH_CM1} cm^-1",
        not_included=[
            "Table 2 (tab:BRtable): branching fractions with NO derived log gf. A line "
            "appearing only there is NOT graded by this paper — which is what makes the "
            "RYA-835 claim that Burheim covers 10875.953/16723.541/21098.84 false.",
            "Nordlander & Lind 2017: uncertainties trace to NIST/Kelleher-Podobedova "
            "2008, so it is not an independent rung.",
        ],
        comparison_columns="loggf_papoulia19 / loggf_kurucz95 / loggf_topbase00 are the "
                           "paper's own P19 / K95 / M00 comparison columns, carried as "
                           "evidence. THEORY, never a lab value; not gradeable.",
        firewall="RYA-161 — this table is a SOURCE. Nothing here adopts a gf into a "
                 "product; that is a pool rebuild (RYA-839).",
    )
    return out, prov


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cache", default="", help="local copy of the e-print tarball")
    a = ap.parse_args(argv)
    df, prov = build(Path(a.cache) if a.cache else None)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    OUT_PROV.write_text(json.dumps(prov, indent=2) + "\n")
    print(f"{CITATION}\n  e-print md5 {prov['eprint_md5']}")
    print(f"  {len(df)} lines, {prov['n_with_nist_levels']} tied to NIST levels")
    print(df[["wavelength_air_A", "lam_vac_A", "loggf", "unc_pct", "e_loggf_dex",
              "elo_eV"]].to_string(index=False))
    print(f"\nwrote {OUT_CSV}\nwrote {OUT_PROV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

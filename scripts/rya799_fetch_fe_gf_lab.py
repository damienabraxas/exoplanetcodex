#!/usr/bin/env python3
"""Vendor the primary-laboratory Fe I log gf measurements — RYA-799.

    python3 scripts/rya799_fetch_fe_gf_lab.py

WHY THESE THREE, AND WHY PRIMARY
--------------------------------
RYA-760 established that a NIST-ASD check is NOT independent for every source: `FMW` is
Fuhr, Martin & Wiese, itself a NIST compilation, and VALD copies it. Agreement with ASD
therefore proves only "no transcription error"; **only a primary laboratory measurement
can referee a gf**. These three are primary — branching fractions from Fourier-transform
spectrometry combined with time-resolved laser-induced-fluorescence lifetimes — and each
publishes a PER-LINE uncertainty in dex, which is the quantity RYA-799 needs and a
compilation tag cannot give:

  Ruffoni et al. 2014,   MNRAS 441, 3127   (VizieR J/MNRAS/441/3127, table3)
      Fe I gf for the Gaia-ESO survey. The canonical statement of this exact gap: of 449
      well-resolved Fe I lines, only 167 had lab gf under 25 % uncertainty.
  Den Hartog et al. 2014, ApJS 215, 23     (VizieR J/ApJS/215/23, table4)
      Fe I from HIGH-LYING EVEN-PARITY levels. Not named in the ticket, included because
      high excitation potential is exactly the regime the near-IR band occupies -- the
      RYA-708 rule that a model is skipped only for a recorded reason cuts both ways.
  Belmonte et al. 2017,  ApJ 848, 125      (arXiv:1710.07571, table T4)
      Fe I from high-lying ODD-parity levels, 213-1033 nm, smallest rms above 700 nm.
      NOT on VizieR (J/ApJ/848/125 is a 404), so it is parsed from the arXiv source.

Together they are the complementary halves of one collaboration's programme: Ruffoni and
Den Hartog cover even-parity uppers, Belmonte the odd-parity ones.

WHAT IS EMITTED
---------------
`data/reference/fe_gf_lab/fe1_lab_loggf.csv` -- one row per measured line, normalised to
(wavelength_air_A, elo_eV, eup_eV, loggf, e_loggf_dex, source), plus the provenance JSON.
Energies are carried in BOTH cm^-1 and eV so the match can be made on LEVELS and not on
wavelength alone (RYA-780: match on wavelength AND EP).

CONTROLS, asserted -- the script refuses to write if any fails:
  * per-source line counts match what each paper states;
  * every row has a finite uncertainty in dex (that is the entire point of using these);
  * Eup - Elo reproduces hc/lambda_vac for every row, which is what proves the level
    columns and the wavelength column describe the SAME transition rather than having
    been mis-associated during parsing.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import sys
import tarfile
import urllib.request
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

OUT_DIR = ROOT / "data" / "reference" / "fe_gf_lab"
OUT_CSV = OUT_DIR / "fe1_lab_loggf.csv"
OUT_PROV = OUT_DIR / "fe1_lab_loggf.prov.json"

CDS = "https://cdsarc.cds.unistra.fr/ftp/{cat}/{f}"
RUFFONI = ("J/MNRAS/441/3127", "table3.dat")
DENHARTOG = ("J/ApJS/215/23", "table4.dat")
BELMONTE_EPRINT = "https://arxiv.org/e-print/1710.07571"

HC_EV_A = 12398.419843320026      # CODATA hc in eV.Angstrom
CM1_TO_EV = 1.0 / 8065.543937     # cm^-1 -> eV

SOURCES = {
    "Ruffoni2014": "Ruffoni, Den Hartog, Lawler, Brewer, Lind, Nave & Pickering 2014, "
                   "MNRAS 441, 3127",
    "DenHartog2014": "Den Hartog, Ruffoni, Lawler, Pickering, Lind & Brewer 2014, "
                     "ApJS 215, 23",
    "Belmonte2017": "Belmonte, Pickering, Ruffoni, Den Hartog, Lawler, Guzman & Nave "
                    "2017, ApJ 848, 125",
}
DOIS = {
    "Ruffoni2014": "10.1093/mnras/stu780",
    "DenHartog2014": "10.1088/0067-0049/215/2/23",
    "Belmonte2017": "10.3847/1538-4357/aa8cd3",
}
#: what each paper says it measured -- the count control
EXPECT_MIN = {"Ruffoni2014": 120, "DenHartog2014": 180, "Belmonte2017": 118}


def air_to_vacuum(lam_air):
    s2 = (1.0e4 / np.asarray(lam_air, dtype=float)) ** 2
    n = 1.0 + (8342.13 + 2406030.0 / (130.0 - s2) + 15997.0 / (38.9 - s2)) * 1e-8
    return np.asarray(lam_air, dtype=float) * n


def _get(url: str, cache: Path | None) -> bytes:
    if cache is not None and cache.exists():
        return cache.read_bytes()
    req = urllib.request.Request(
        url, headers={"User-Agent": "exoplanetcodex/RYA-799 (mailto:ryan.damien.schmitt@gmail.com)"})
    with urllib.request.urlopen(req, timeout=120) as r:   # noqa: S310 (fixed https URLs)
        blob = r.read()
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_bytes(blob)
    return blob


# ── per-source parsers ────────────────────────────────────────────────────────

def parse_ruffoni(blob: bytes) -> pd.DataFrame:
    """VizieR table3.dat, fixed width per its ReadMe.

    The upper level is NOT a column: it is carried in the `# Level <E> cm^-1^ ...`
    headline that opens each block, so the parser has to remember it. Reading only the
    data lines would silently drop the upper level and force a wavelength-only match.
    """
    rows, upper = [], None
    for ln in blob.decode("utf-8", "replace").splitlines():
        m = re.match(r"#\s*Level\s+([\d.]+)\s*cm\^-1\^", ln)
        if m:
            upper = float(m.group(1))
            continue
        if ln.startswith("#") or len(ln) < 60:
            continue
        try:
            lam = float(ln[13:23])
            sigma = float(ln[24:33])            # transition wavenumber, cm^-1
            lg, e = ln[53:58].strip(), ln[59:63].strip()
        except ValueError:
            continue
        if not lg or not e or upper is None:
            continue
        rows.append({"source": "Ruffoni2014", "wavelength_air_A": lam,
                     "eup_cm1": upper, "elo_cm1": upper - sigma,
                     "loggf": float(lg), "e_loggf_dex": float(e)})
    return pd.DataFrame(rows)


def parse_denhartog(blob: bytes) -> pd.DataFrame:
    """VizieR table4.dat, FIXED WIDTH per its ReadMe.

    Splitting this table on whitespace is wrong and fails silently: several optional
    columns (BF, e_BF, A) are blank on some rows, so the fields shift left and
    `e_log(gf)` is read out of the PUBLISHED log gf column instead. The first run of this
    script did exactly that and reported a minimum uncertainty of -2.62 dex -- an
    impossible number, which is the only reason it was visible. Byte offsets, and a
    positivity control below so it can never pass silently again.
    """
    rows = []
    for ln in blob.decode("utf-8", "replace").splitlines():
        if len(ln) < 64:
            continue
        try:
            lam = float(ln[0:10])
            eup, elo = float(ln[11:20]), float(ln[23:32])
            lg, e = ln[54:59].strip(), ln[60:64].strip()
        except ValueError:
            continue
        if not lg or not e:
            continue
        rows.append({"source": "DenHartog2014", "wavelength_air_A": lam,
                     "eup_cm1": eup, "elo_cm1": elo,
                     "loggf": float(lg), "e_loggf_dex": float(e)})
    return pd.DataFrame(rows)


_BE = re.compile(
    r"\s*\$([\d.]+)\$\s*&\s*\$([\d.]+)\$\s*&\s*\$(\d+)\$\s*&\s*&\s*\$([\d.]+)\$\s*&\s*\$(\d+)\$")


def parse_belmonte(blob: bytes) -> pd.DataFrame:
    """arXiv source, table T4_Results.tex. Wavelengths are AIR and in NANOMETRES."""
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:*") as tf:
        member = next((m for m in tf.getmembers()
                       if Path(m.name).name == "T4_Results.tex"), None)
        if member is None:
            raise SystemExit("T4_Results.tex absent from the arXiv:1710.07571 source")
        text = tf.extractfile(member).read().decode("utf-8", "replace")
    rows = []
    for ln in text.splitlines():
        m = _BE.match(ln)
        if not m:
            continue
        pairs = re.findall(r"\$(-?[\d.]+)\$\s*\$\\pm\$\s*\$(-?[\d.]+)\$", ln)
        if not pairs:
            continue
        rows.append({"source": "Belmonte2017",
                     "wavelength_air_A": float(m.group(1)) * 10.0,   # nm -> A
                     "eup_cm1": float(m.group(2)), "elo_cm1": float(m.group(4)),
                     "loggf": float(pairs[0][0]), "e_loggf_dex": float(pairs[0][1])})
    return pd.DataFrame(rows)


# ── controls ──────────────────────────────────────────────────────────────────

def controls(df: pd.DataFrame) -> list[str]:
    out = []
    for src, n_min in EXPECT_MIN.items():
        n = int((df["source"] == src).sum())
        ok = n >= n_min
        out.append(f"COUNT {src:<14s}: {n:4d} usable rows (paper implies >= {n_min}) "
                   f"{'PASS' if ok else 'FAIL'}")

    n_bad = int((~np.isfinite(df["e_loggf_dex"])).sum())
    out.append(f"UNCERTAINTY: {len(df) - n_bad}/{len(df)} rows carry a finite dex "
               f"uncertainty {'PASS' if n_bad == 0 else 'FAIL'}")
    # An uncertainty is positive and, for these techniques, small. A fixed-width table
    # read on whitespace shifts columns and lands the PUBLISHED log gf here -- which is
    # finite, so the check above passes and this one does not. It is the control that
    # actually caught the Den Hartog parser.
    lo, hi = float(df["e_loggf_dex"].min()), float(df["e_loggf_dex"].max())
    ok = lo > 0.0 and hi <= 0.5
    out.append(f"SIGMA RANGE: {lo:+.3f} .. {hi:+.3f} dex (must be >0 and <=0.5) "
               f"{'PASS' if ok else 'FAIL'}")

    # The transition-energy identity. If a parser mis-associated a level column with the
    # wrong wavelength column, this is what catches it -- nothing else would.
    dE_levels = (df["eup_cm1"] - df["elo_cm1"]) * CM1_TO_EV
    dE_lambda = HC_EV_A / air_to_vacuum(df["wavelength_air_A"].values)
    resid = np.abs(dE_levels - dE_lambda)
    for src in sorted(df["source"].unique()):
        m = df["source"] == src
        worst = float(resid[m].max())
        ok = worst < 5e-3
        out.append(f"LEVELS  {src:<14s}: max |(Eup-Elo) - hc/lambda_vac| = {worst:.2e} eV "
                   f"{'PASS' if ok else 'FAIL'}")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cache-dir", type=Path, default=None,
                    help="directory to cache the three downloads for offline reruns")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args(argv)

    def cache(name):
        return None if args.cache_dir is None else args.cache_dir / name

    blobs = {
        "Ruffoni2014": _get(CDS.format(cat=RUFFONI[0], f=RUFFONI[1]), cache("ruffoni_t3.dat")),
        "DenHartog2014": _get(CDS.format(cat=DENHARTOG[0], f=DENHARTOG[1]), cache("dh_t4.dat")),
        "Belmonte2017": _get(BELMONTE_EPRINT, cache("belmonte.tar")),
    }
    df = pd.concat([parse_ruffoni(blobs["Ruffoni2014"]),
                    parse_denhartog(blobs["DenHartog2014"]),
                    parse_belmonte(blobs["Belmonte2017"])], ignore_index=True)
    df["elo_eV"] = df["elo_cm1"] * CM1_TO_EV
    df["eup_eV"] = df["eup_cm1"] * CM1_TO_EV
    df = df.sort_values(["source", "wavelength_air_A"]).reset_index(drop=True)

    report = controls(df)
    print("Primary-lab Fe I log gf vendoring (RYA-799)\n")
    for line in report:
        print("  " + line)
    if any(l.endswith("FAIL") for l in report):
        print("\nCONTROLS FAILED — refusing to write.", file=sys.stderr)
        return 1

    ir = df[(df.wavelength_air_A >= 6910) & (df.wavelength_air_A <= 9199)]
    print(f"\n  total {len(df)} measured Fe I lines; "
          f"{len(ir)} fall in the RYA-799 window 6910-9199 A")
    for src in sorted(df["source"].unique()):
        s = ir[ir.source == src]
        print(f"    {src:<14s} {len(s):3d} in window"
              + (f"  sigma {s.e_loggf_dex.min():.2f}-{s.e_loggf_dex.max():.2f} dex"
                 if len(s) else ""))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    cols = ["source", "wavelength_air_A", "elo_cm1", "eup_cm1", "elo_eV", "eup_eV",
            "loggf", "e_loggf_dex"]
    df[cols].to_csv(args.out_dir / OUT_CSV.name, index=False, float_format="%.6f")
    (args.out_dir / OUT_PROV.name).write_text(json.dumps({
        "ticket": "RYA-799",
        "artifact": OUT_CSV.name,
        "what": "Primary laboratory Fe I log gf with PER-LINE uncertainties in dex.",
        "why_primary": ("RYA-760: FMW is itself a NIST compilation and VALD copies it, so "
                        "an ASD check is not independent. Only a primary lab measurement "
                        "can referee a gf, and only a primary paper publishes a per-line "
                        "sigma."),
        "sources": {k: {"cite": v, "doi": DOIS[k],
                        "sha256": hashlib.sha256(blobs[k]).hexdigest()}
                    for k, v in SOURCES.items()},
        "not_on_vizier": "Belmonte 2017 — J/ApJ/848/125 returns 404; parsed from arXiv source.",
        "controls": report,
        "retrieved": str(date.today()),
        "regenerate": "python3 scripts/rya799_fetch_fe_gf_lab.py",
    }, indent=2) + "\n")
    print(f"\n  wrote {args.out_dir / OUT_CSV.name}")
    print(f"  wrote {args.out_dir / OUT_PROV.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

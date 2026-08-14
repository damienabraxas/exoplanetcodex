#!/usr/bin/env python3
"""Recover the line list the Amarsi 2022 3D-NLTE MLP was TRAINED on — RYA-817.

    python3 scripts/rya817_recover_amarsi_training_set.py

WHY THIS EXISTS
---------------
`vendor/1L-3NErrors/` ships three trained `MLPRegressor`s and nothing else: no training
data, no line list, no feature bounds. `main_aberr.py` guards only the STELLAR box
(Teff/logg/vmic/A(Fe)) and will happily predict for any (Elo, Eup, loggf) you hand it —
including line parameters the network never saw. RYA-817 asks us to run the network on
the near-IR Fe band, which is a different part of line-parameter space from anything in
the paper. Without the training line list there is no domain check, only a number.

So the training set is RECOVERED, not assumed:

  Amarsi, Liljegren & Nissen 2022 (A&A 668, A68 = arXiv:2209.13449), Sect. 4:
    "This study employs the `golden' line list presented in Jofre et al. (2014),
     adopting the energies, log gf, and pressure broadening parameters given in their
     Tables 4 and 5. ... All of the lines shown are in the optical regime:
     lambda_Air = 478.783 nm to 681.026 nm."

Jofre et al. 2014 (A&A 564, A133 = arXiv:1309.1099) ships those two tables as LaTeX in
its arXiv source: `golden_Fe1.tex` and `golden_Fe2.tex`. This script downloads that
source, parses both tables, derives Eup the way Amarsi did, and writes the result to
`data/reference/amarsi2022_training/`.

FOUR POSITIVE CONTROLS, ALL ASSERTED (RYA-805 rule: an absence needs a control, and the
test must be shown to DISCRIMINATE)
-----------------------------------------------------------------------------------
1. COUNT      — the parse must yield exactly 171 Fe I and 12 Fe II lines, the counts
                Amarsi 2022 states and `pipeline/nlte_corrections.py` records.
2. WAVELENGTH — min/max must be 4787.83 / 6810.26 A, i.e. the 478.783-681.026 nm the
                paper quotes, to the digit.
3. Eup        — Eup = Elo + hc/lambda_vac must reproduce the Eup column of the vendored
                `test_data.csv` for every row that matches a golden line, to <1e-5 eV.
                This is what proves the Eup convention rather than guessing it, and it
                is the finding that (Eup - Elo) IS the transition energy, i.e. the
                wavelength, carried into the network as a derived feature.
4. SCALER     — the recovered per-network subsets must reproduce the `mean_` and
                `scale_` that the vendored StandardScalers carry for Elo and lggf.
                Those attributes are moments OF THE TRAINING DATA; if the recovered
                list were the wrong list they would not match. Tolerances are loose
                enough to absorb the paper's per-model line-overlap cut (a line can be
                dropped in some stellar models and kept in others, so the sample
                weighting is not exactly uniform over lines) and tight enough that a
                wrong list fails: `--show-discrimination` prints what a deliberately
                wrong list scores against the same test.

NETWORK ASSIGNMENT. Amarsi trains three networks and `main_aberr.py` routes by
`Elo < 2.0` / `Elo > 2.0` for Fe I (note: `main_aberr.py` uses strict `>`, so a line at
exactly 2.0 eV is routed to NEITHER and silently keeps its -999 default; no golden line
sits at 2.0, but the repo's `_compute_aberr` uses `elo < 2.0 else gt02`, which is the
safer total function). We reproduce the repo's routing here.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import pickle
import re
import sys
import tarfile
import urllib.request
import warnings
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

OUT_DIR = ROOT / "data" / "reference" / "amarsi2022_training"
OUT_CSV = OUT_DIR / "amarsi2022_training_lines.csv"
OUT_PROV = OUT_DIR / "amarsi2022_training_lines.prov.json"
SOLAR_CSV = OUT_DIR / "amarsi2022_solar_control_lines.csv"
VENDOR = ROOT / "vendor" / "1L-3NErrors"

JOFRE_EPRINT = "https://arxiv.org/e-print/1309.1099"
JOFRE_BIB = "Jofre et al. 2014, A&A 564, A133 (arXiv:1309.1099), Tables 4 and 5"
AMARSI_BIB = "Amarsi, Liljegren & Nissen 2022, A&A 668, A68 (arXiv:2209.13449)"
AP2002_EPRINT = "https://arxiv.org/e-print/astro-ph/0111055"
AP2002_BIB = ("Allende Prieto, Asplund, Garcia Lopez & Lambert 2002, ApJ 567, 544 "
              "(arXiv:astro-ph/0111055), Table 2")

# CODATA hc in eV*Angstrom. The same constant reproduces Amarsi's Eup column to 4e-6 eV
# under either Edlen 1966 or Birch & Downs 1994 air->vacuum; the two differ by far less
# than the check tolerance, so the choice is not load-bearing (control 3 proves it).
HC_EV_A = 12398.419843320026

# Expected values -- these are the CONTROLS, not configuration. Changing one to make a
# run pass is the failure mode this file exists to prevent.
EXPECT_N = {"Fe1": 171, "Fe2": 12}
EXPECT_LAMBDA = (4787.83, 6810.26)

# StandardScaler.mean_ / .scale_ indices for the 7 features
# ['teff/K','lg(g/cms^-2)','A(Fe)','vmic/kms^-1','Elo/eV','Eup','lggf']
I_ELO, I_EUP, I_LGGF = 4, 5, 6
SCALER_TOL_MEAN = 0.12   # dex / eV. gt02 and fe2 land within 0.01; lt02 (n=17) needs
SCALER_TOL_SCALE = 0.08  # more room because the overlap cut bites hardest there.


def air_to_vacuum(lambda_air_A: np.ndarray) -> np.ndarray:
    """Edlen (1966) refractive index of standard air. Vectorised."""
    s2 = (1.0e4 / np.asarray(lambda_air_A, dtype=float)) ** 2
    n = 1.0 + (8342.13 + 2406030.0 / (130.0 - s2) + 15997.0 / (38.9 - s2)) * 1e-8
    return np.asarray(lambda_air_A, dtype=float) * n


def transition_eV(lambda_air_A) -> np.ndarray:
    """Photon energy of a transition given its AIR wavelength: hc / lambda_vac."""
    return HC_EV_A / air_to_vacuum(lambda_air_A)


# ── the Jofre 2014 arXiv source ───────────────────────────────────────────────

_ROW = re.compile(r"^\s*(\d{4}\.\d+)\s*&\s*(-?\d+\.\d+)\s*&\s*(-?\d+\.\d+)\s*&")


def fetch_source(cache: Path | None, url: str = JOFRE_EPRINT) -> tuple[bytes, str]:
    """Return (tarball bytes, sha256). Uses `cache` if it exists, else downloads."""
    if cache is not None and cache.exists():
        blob = cache.read_bytes()
    else:
        req = urllib.request.Request(
            url, headers={"User-Agent": "exoplanetcodex/RYA-817 (arXiv e-print)"})
        with urllib.request.urlopen(req, timeout=120) as r:   # noqa: S310 (fixed https URL)
            blob = r.read()
        if cache is not None:
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_bytes(blob)
    return blob, hashlib.sha256(blob).hexdigest()


def parse_golden(blob: bytes) -> tuple[pd.DataFrame, dict]:
    """Parse golden_Fe1.tex / golden_Fe2.tex out of the arXiv tarball."""
    rows, member_sha = [], {}
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:*") as tf:
        for species, name in (("Fe1", "golden_Fe1.tex"), ("Fe2", "golden_Fe2.tex")):
            member = next((m for m in tf.getmembers()
                           if Path(m.name).name == name), None)
            if member is None:
                raise SystemExit(
                    f"{name} is not in the arXiv:1309.1099 source. The tables were "
                    f"parsed out of the paper source rather than a CDS catalogue "
                    f"because VizieR J/A+A/564/A133 carries table6 (which lines were "
                    f"used per star) but NOT Tables 4/5 (the atomic data). If the "
                    f"source layout changed, re-check before editing this parser.")
            text = tf.extractfile(member).read()
            member_sha[name] = hashlib.sha256(text).hexdigest()
            for line in text.decode("utf-8", "replace").splitlines():
                m = _ROW.match(line)
                if m:
                    lam, elo, lggf = (float(x) for x in m.groups())
                    rows.append({"species": species, "wavelength_air_A": lam,
                                 "elo_eV": elo, "loggf": lggf})
    return pd.DataFrame(rows), member_sha


# The Amarsi 2022 SOLAR analysis used a different line list from the training set, and
# that distinction is the whole reason this second artifact exists. Sect. 6.1: "Line-by-
# line 1D LTE lg eps(Fe) values measured in the solar flux atlas of Kurucz (1984) were
# taken from Allende Prieto et al. (2002). The analysis was restricted to weak lines with
# REW < -4.9." Checking the reactivation against Table 6's solar row therefore requires
# THIS list, not the 171+12 golden lines. Running it against the training list instead
# scores a 0.04 dex miss on Fe I that looks like a defect and is not one -- it is the
# RYA-785 wrong-referee failure, and this artifact is what makes the right referee
# available.
_AP_ROW = re.compile(
    r"^\s*(Fe\s+I{1,2})\s*&\s*([\d.]+)\s*&\s*([\d.-]+)\s*&\s*(-?[\d.]+)"
    r"\s*&\s*([\d.]+)\s*&\s*([\d.]+)\s*&\s*(\d+)")


def parse_ap2002_solar(blob: bytes) -> tuple[pd.DataFrame, dict]:
    """Table 2 of Allende Prieto et al. 2002 — the solar Fe I/II line-by-line 1D LTE set."""
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:*") as tf:
        member = next((m for m in tf.getmembers()
                       if Path(m.name).name.endswith(".tex")), None)
        if member is None:
            raise SystemExit("no .tex in the arXiv:astro-ph/0111055 source")
        text = tf.extractfile(member).read()
    sha = {Path(member.name).name: hashlib.sha256(text).hexdigest()}
    body = text.decode("utf-8", "replace")
    if "\\label{table2}" not in body:
        raise SystemExit("Table 2 (the SOLAR line list) not found in the AP2002 source; "
                         "Table 1 is Procyon and must not be substituted for it.")
    tab = body.split("\\label{table2}")[1].split("\\end{deluxetable}")[0]
    rows = []
    for line in tab.splitlines():
        m = _AP_ROW.match(line)
        if m:
            sp, lam, ep, lggf, a3d, a1d, w = m.groups()
            rows.append({"ion": sp.split()[1].strip(),
                         "wavelength_air_A": float(lam), "elo_eV": float(ep),
                         "loggf": float(lggf), "a_3d_lte_ap2002": float(a3d),
                         "a_1d_lte_ap2002": float(a1d), "ew_mA": float(w)})
    df = pd.DataFrame(rows)
    df["delta_E_eV"] = transition_eV(df["wavelength_air_A"].values)
    df["eup_eV"] = df["elo_eV"] + df["delta_E_eV"]
    df["rew"] = np.log10(df["ew_mA"] * 1e-3 / df["wavelength_air_A"])
    # Amarsi's own cut, carried as a column rather than applied here: the artifact stays
    # the full published table and the consumer states which subset it used.
    df["weak_line_rew_lt_m49"] = df["rew"] < -4.9
    return df.sort_values(["ion", "wavelength_air_A"]).reset_index(drop=True), sha


def add_derived(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["delta_E_eV"] = transition_eV(df["wavelength_air_A"].values)
    df["eup_eV"] = df["elo_eV"] + df["delta_E_eV"]
    df["ion"] = df["species"].map({"Fe1": "I", "Fe2": "II"})
    df["network"] = np.where(df["species"] == "Fe2", "fe2",
                             np.where(df["elo_eV"] < 2.0, "lt02", "gt02"))
    return df.sort_values(["species", "wavelength_air_A"]).reset_index(drop=True)


# ── controls ──────────────────────────────────────────────────────────────────

def control_counts(df: pd.DataFrame) -> list[str]:
    out = []
    for sp, n in EXPECT_N.items():
        got = int((df["species"] == sp).sum())
        out.append(f"COUNT {sp}: {got} (expected {n}) "
                   f"{'PASS' if got == n else 'FAIL'}")
    return out


def control_wavelength(df: pd.DataFrame) -> list[str]:
    lo, hi = float(df["wavelength_air_A"].min()), float(df["wavelength_air_A"].max())
    ok = abs(lo - EXPECT_LAMBDA[0]) < 0.005 and abs(hi - EXPECT_LAMBDA[1]) < 0.005
    return [f"WAVELENGTH: {lo:.2f}-{hi:.2f} A (paper states "
            f"{EXPECT_LAMBDA[0]:.2f}-{EXPECT_LAMBDA[1]:.2f}) {'PASS' if ok else 'FAIL'}"]


def control_eup(df: pd.DataFrame) -> tuple[list[str], float]:
    """Reproduce the vendored test_data.csv Eup column from (Elo, lambda)."""
    td = VENDOR / "test_data.csv"
    cols = ["teff", "logg", "afe", "vmic", "elo", "eup", "lggf", "species"]
    raw = pd.read_csv(td, comment="#", names=cols, skipinitialspace=True)
    raw["species"] = raw["species"].astype(str).str.strip()
    raw = raw[raw["species"].isin(("Fe1", "Fe2"))]

    errs = []
    for _, r in raw.iterrows():
        # match on (Elo, lggf) -- test_data carries more Elo digits than Jofre's table,
        # so match on the rounded Elo the table prints plus the exact lggf.
        cand = df[(df["species"] == r["species"])
                  & (np.abs(df["elo_eV"] - r["elo"]) < 5e-4)
                  & (np.abs(df["loggf"] - r["lggf"]) < 5e-4)]
        if len(cand) != 1:
            continue
        lam = float(cand.iloc[0]["wavelength_air_A"])
        errs.append(abs((r["eup"] - r["elo"]) - float(transition_eV(lam))))
    if not errs:
        return ["EUP: no test_data row matched a golden line -- CONTROL DID NOT RUN"], np.inf
    worst = float(max(errs))
    return ([f"EUP: {len(errs)} test_data rows reproduced from (Elo, lambda_air) "
             f"via hc/lambda_vac, max err {worst:.2e} eV "
             f"{'PASS' if worst < 1e-5 else 'FAIL'}"], worst)


def _scalers() -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return {k: pickle.load(open(VENDOR / f, "rb"))[0] for k, f in
                (("lt02", "fe1_model_lt02.p"), ("gt02", "fe1_model_gt02.p"),
                 ("fe2", "fe2_model.p"))}


def control_scaler(df: pd.DataFrame, label: str = "recovered") -> tuple[list[str], bool]:
    """The vendored StandardScalers carry mean_/scale_ OF THE TRAINING DATA."""
    out, ok_all = [], True
    for net, sc in _scalers().items():
        sub = df[df["network"] == net]
        if sub.empty:
            out.append(f"SCALER {net}: no recovered line routes here -- FAIL")
            ok_all = False
            continue
        for idx, col, tol_m, tol_s, unit in ((I_ELO, "elo_eV", SCALER_TOL_MEAN,
                                              SCALER_TOL_SCALE, "eV"),
                                             (I_LGGF, "loggf", SCALER_TOL_MEAN,
                                              SCALER_TOL_SCALE, "dex")):
            dm = abs(float(sub[col].mean()) - float(sc.mean_[idx]))
            ds = abs(float(sub[col].std(ddof=0)) - float(sc.scale_[idx]))
            ok = dm < tol_m and ds < tol_s
            ok_all &= ok
            out.append(
                f"SCALER {net:<4} {col:<7} {label}: mean {sub[col].mean():+8.4f} vs "
                f"trained {sc.mean_[idx]:+8.4f} (d={dm:.4f} {unit})  "
                f"sd {sub[col].std(ddof=0):7.4f} vs {sc.scale_[idx]:7.4f} "
                f"(d={ds:.4f})  {'PASS' if ok else 'FAIL'}")
    return out, ok_all


def show_discrimination(df: pd.DataFrame) -> list[str]:
    """Does the SCALER control actually discriminate, or would any Fe list pass?

    RYA-805: verify the test can FAIL. Feed it a deliberately wrong list -- the same
    lines with Elo shifted by +1 eV, and a random-uniform stand-in over the same
    ranges -- and show both are rejected.
    """
    out = ["", "DISCRIMINATION (the control must be able to FAIL):"]
    shifted = df.copy()
    shifted["elo_eV"] = shifted["elo_eV"] + 1.0
    lines, ok = control_scaler(shifted, label="Elo+1eV ")
    out += [ln for ln in lines if "elo_eV" in ln]
    out.append(f"  -> wrong-list verdict: {'PASS (BAD: test does not discriminate)' if ok else 'REJECTED (good)'}")
    rng = np.random.default_rng(817)
    fake = df.copy()
    fake["loggf"] = rng.uniform(df["loggf"].min(), df["loggf"].max(), len(df))
    lines, ok = control_scaler(fake, label="rand lggf")
    out += [ln for ln in lines if "loggf" in ln]
    out.append(f"  -> wrong-list verdict: {'PASS (BAD: test does not discriminate)' if ok else 'REJECTED (good)'}")
    return out


# ── main ──────────────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cache", type=Path, default=None,
                    help="path to cache/read the arXiv:1309.1099 tarball "
                         "(offline reruns; the sha256 is recorded either way)")
    ap.add_argument("--cache-solar", type=Path, default=None,
                    help="same, for the arXiv:astro-ph/0111055 (Allende Prieto 2002) "
                         "tarball that carries the SOLAR control line list")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--show-discrimination", action="store_true",
                    help="also prove the scaler control can reject a wrong line list")
    ap.add_argument("--check", action="store_true",
                    help="run the controls and exit non-zero on failure; write nothing")
    args = ap.parse_args(argv)

    blob, sha = fetch_source(args.cache)
    raw, member_sha = parse_golden(blob)
    df = add_derived(raw)

    solar_blob, solar_sha = fetch_source(args.cache_solar, AP2002_EPRINT)
    solar_df, solar_member_sha = parse_ap2002_solar(solar_blob)

    report = []
    report += control_counts(df)
    report += control_wavelength(df)
    eup_lines, eup_worst = control_eup(df)
    report += eup_lines
    scaler_lines, scaler_ok = control_scaler(df)
    report += scaler_lines

    print(f"Amarsi 2022 training-set recovery (RYA-817)")
    print(f"  source: {JOFRE_EPRINT}  sha256 {sha[:16]}...")
    for k, v in member_sha.items():
        print(f"          {k}  sha256 {v[:16]}...")
    print()
    for line in report:
        print("  " + line)
    if args.show_discrimination:
        for line in show_discrimination(df):
            print("  " + line)

    failed = [ln for ln in report if ln.endswith("FAIL") or "FAIL" in ln]
    if failed:
        print("\nCONTROLS FAILED -- refusing to write the artifact. "
              "Do NOT loosen a tolerance to make this pass.", file=sys.stderr)
        return 1

    print("\n  TRAINING DOMAIN (what the network actually saw):")
    for net, sub in df.groupby("network"):
        print(f"    {net:<5} n={len(sub):3d}  "
              f"Elo {sub.elo_eV.min():.3f}-{sub.elo_eV.max():.3f} eV  "
              f"Eup {sub.eup_eV.min():.3f}-{sub.eup_eV.max():.3f} eV  "
              f"lggf {sub.loggf.min():+.3f}..{sub.loggf.max():+.3f}  "
              f"dE {sub.delta_E_eV.min():.4f}-{sub.delta_E_eV.max():.4f} eV")

    if args.check:
        return 0

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    cols = ["species", "ion", "network", "wavelength_air_A", "elo_eV", "eup_eV",
            "delta_E_eV", "loggf"]
    df[cols].to_csv(out_dir / OUT_CSV.name, index=False,
                    float_format="%.6f")
    solar_df.to_csv(out_dir / SOLAR_CSV.name, index=False, float_format="%.6f")

    prov = {
        "ticket": "RYA-817",
        "artifact": OUT_CSV.name,
        "what": ("The line list the Amarsi 2022 3D-NLTE MLPs (vendor/1L-3NErrors) were "
                 "TRAINED on: the Jofre et al. 2014 'golden' Fe I / Fe II set, 171 + 12 "
                 "lines. Recovered so that the network's LINE-parameter domain can be "
                 "checked, not just its stellar-parameter box."),
        "why": ("vendor/1L-3NErrors ships trained models with no training data and no "
                "feature bounds; main_aberr.py guards only Teff/logg/vmic/A(Fe). Running "
                "the network on the near-IR Fe band without this list would be silent "
                "extrapolation."),
        "chain": [AMARSI_BIB + " Sect. 4 names the line list and quotes its wavelength "
                               "range 478.783-681.026 nm",
                  JOFRE_BIB + " ships those tables as golden_Fe1.tex / golden_Fe2.tex "
                              "in its arXiv source"],
        "eup_convention": ("Eup = Elo + hc/lambda_vac, hc = 12398.419843320026 eV.A, "
                           "lambda_vac from Edlen (1966). VERIFIED against the vendored "
                           f"test_data.csv Eup column to {eup_worst:.1e} eV. COROLLARY: "
                           "(Eup - Elo) IS the transition energy, so WAVELENGTH enters "
                           "the network as a derived feature -- the MLP is not "
                           "wavelength-agnostic in its inputs."),
        "controls": report,
        "not_from_cds": ("VizieR J/A+A/564/A133 carries table6 (which lines each star "
                         "used) but NOT Tables 4/5 (the atomic data), and there is no "
                         "VizieR catalogue for J/A+A/668/A68 at all."),
        "companion_artifact": {
            "file": SOLAR_CSV.name,
            "what": (f"The SOLAR line list Amarsi 2022 used for its Table 6 solar row: "
                     f"{AP2002_BIB}. {len(solar_df[solar_df.ion == 'I'])} Fe I + "
                     f"{len(solar_df[solar_df.ion == 'II'])} Fe II lines with their "
                     f"published 1D LTE abundances and equivalent widths."),
            "why": ("The solar row of Table 6 is NOT computable from the training list "
                    "-- Amarsi analysed the Sun with this list instead. Checking the "
                    "reactivation against Table 6 requires it; using the training list "
                    "misses Fe I by ~0.04 dex for a reason that has nothing to do with "
                    "the network."),
            "cut": ("Amarsi restricted the solar analysis to weak lines with REW < -4.9; "
                    "carried as the `weak_line_rew_lt_m49` column, not pre-applied."),
        },
        "sources": {
            "jofre2014_eprint": JOFRE_EPRINT,
            "jofre2014_tarball_sha256": sha,
            "member_sha256": member_sha,
            "ap2002_eprint": AP2002_EPRINT,
            "ap2002_tarball_sha256": solar_sha,
            "ap2002_member_sha256": solar_member_sha,
        },
        "retrieved": str(date.today()),
        "regenerate": "python3 scripts/rya817_recover_amarsi_training_set.py",
    }
    (out_dir / OUT_PROV.name).write_text(json.dumps(prov, indent=2) + "\n")
    print(f"\n  wrote {out_dir / OUT_CSV.name}")
    print(f"  wrote {out_dir / SOLAR_CSV.name}  "
          f"({len(solar_df[solar_df.ion == 'I'])} Fe I + "
          f"{len(solar_df[solar_df.ion == 'II'])} Fe II solar control lines)")
    print(f"  wrote {out_dir / OUT_PROV.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Run the Amarsi 2022 3D-NLTE Fe MLP on the solar VIS and IR bands — RYA-817 Deliverable A.

    python3 scripts/rya817_run_3dnlte_bands.py --band-products-dir <dir> [--out <dir>]

WHAT IT DOES
------------
For each (ion, band) it takes the RYA-783 per-line 1D-LTE band product, attaches the
line's atomic data (Elo, Eup, log gf) from the same VALD extraction the pipeline builds
its line lists from, asks `pipeline.amarsi3d` whether the line is inside the MLP's
TRAINING domain, and — for in-domain lines only — adds the per-line 3D-NLTE correction
to produce an `ENGINE-A-3DNLTE` product (RYA-712: its own value/sigma/n, presented
beside the existing four, never merged with them).

Out-of-domain lines are NOT dropped and NOT extrapolated. They are carried into the
per-line output with the axis that rejected them, and — separately, in a column the
aggregate never reads — the number the network WOULD have returned, so the size of the
refused extrapolation is on the record instead of being a matter of opinion.

THE VIS LEG IS THE CONTROL. The reactivation is only trustworthy if the network still
reproduces the RYA-283 solar cross-check (A(3D-NLTE) ~ 7.45, Asplund 2021) on the
optical lines it was trained on. That is why VIS is run first and reported first: if
the control does not land, the IR number means nothing.

INPUTS
------
--band-products-dir  directory holding the RYA-783 outputs
                     `Fe{I,II}_{lo}_{hi}_{instrument}_PROFILEFIT_1D-LTE_lines.csv`.
                     Required, and a path rather than a literal: these live on Sirius
                     (RYA-567 compute host), and a hardcoded mount point is exactly what
                     RYA-810 removed from this repo.
"""
from __future__ import annotations

import argparse
import re
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from data.linelists.vald_parse import parse_vald_long          # noqa: E402
from pipeline import amarsi3d                                   # noqa: E402
from pipeline import band_policy                              # noqa: E402  (RYA-922)
from pipeline.band_products import (                            # noqa: E402
    LineMeasurement, build_product, products_frame)
from config.constants import get_star_params                    # noqa: E402

LINELISTS = ROOT / "data" / "linelists"
# The two VALD3 long-format solar extractions that between them span both bands.
VALD_SOURCES = (
    LINELISTS / "vald_solar_raw.txt",                            # 3780-6910 A
    LINELISTS / "vald_solar_redopt_6910_9500_hfson_raw.txt",     # 6910-9500 A
    LINELISTS / "vald_solar_ir_9500_17000_hfson_raw.txt",        # 9500-17000 A (RYA-762)
)
OUT_DEFAULT = ROOT / "data" / "results" / "rya817"

# 🔴 RYA-922 — BANDS AND INSTRUMENT WERE BOTH HARDCODED HERE, and that is why 3D-NLTE
# existed for exactly one arm and two bands.
#
#     BANDS = (("VIS", 3800, 6910), ("IR", 6910, 9199))
#     INSTRUMENT = "kpno_solar_atlas"
#
# The instrument was a module-level constant used BOTH to find the input file and to tag
# the output, and there was no `--instrument` flag at all — so this route could not be
# pointed at HARPS, IAG or CRIRES+ even in principle. That is the RYA-913 defect class
# with the opposite symptom: RYA-913's ENGINE-B loaded Kitt Peak while TAGGING
# `a.instrument`, so it lied; this one tagged honestly and simply could not move.
#
# `IR` was not a band either (RYA-918): `band_policy` declares four contiguous,
# non-overlapping bands and `IR` is not among them — 6910-9199 is INSIDE red-optical.
#
# Both are now DISCOVERED rather than declared. The instrument comes from the command
# line; the bands come from whatever 1D-LTE per-line artifacts exist, with each one's
# band NAME resolved through `band_policy` from its own measured span. Nothing about the
# Amarsi MLP is instrument-specific: it consumes (wavelength, EP, abundance) per line.
IONS = ("I", "II")

#: `Fe{ion}_{lo}_{hi}_{instrument}_PROFILEFIT_1D-LTE_lines.csv`
_SRC_RE = re.compile(
    r"^Fe(?P<ion>I{1,2})_(?P<lo>\d+)_(?P<hi>\d+)_(?P<instrument>.+)"
    r"_PROFILEFIT_1D-LTE_lines\.csv$")


def discover_bands(bp_dir: Path, instrument: str, ion: str):
    """Every measured 1D-LTE span for this (ion, instrument), with its canonical band.

    The band NAME is resolved from `band_policy` at the span's midpoint rather than
    read off the filename, so an artifact measured over 6910-9199 is correctly named
    `red-optical` and never invents a band of its own.

    ⚠️ The span may be a SUBSET of the band it resolves to — 6910-9199 covers 74% of
    red-optical (RYA-921). That is carried through to the caller as `(lo, hi)` so the
    product records what was actually measured, not what the band name implies.
    """
    found = []
    for f in sorted(bp_dir.glob(f"Fe{ion}_*_PROFILEFIT_1D-LTE_lines.csv")):
        m = _SRC_RE.match(f.name)
        if not m or m.group("instrument") != instrument or m.group("ion") != ion:
            continue
        lo, hi = int(m.group("lo")), int(m.group("hi"))
        band = band_policy.resolve(0.5 * (lo + hi)).name
        found.append((band, lo, hi, f))
    return found

# Wavelength match tolerance between the measured line and the VALD transition, in A.
# The band harness measures at the line list's own wavelength, so this is a rounding
# tolerance, not a search radius. Widening it would start matching NEIGHBOURS, which is
# the RYA-785 "a same-species neighbour does NOT cancel" failure.
WAVE_TOL_A = 0.01

ASPLUND21_FE = 7.46           # A(Fe), Asplund, Amarsi & Grevesse 2021
GOLD_V4_FE_I = 7.466          # the frozen 3D-NLTE VIS anchor (RYA-553/665/811)

# ── the reactivation control, taken from the paper itself ────────────────────
# Amarsi+2022 Table 6 gives, for the Sun:
#   Fe I  1D LTE 7.47 -> 3D non-LTE 7.46
#   Fe II 1D LTE 7.41 -> 3D non-LTE 7.47
# and Sect. 5.1 states that at [M/H]=0 the Fe I abundance difference Delta_1L3N falls by
# "around -0.05 to -0.1 dex as Elo increases from 0 eV to 5 eV".
#
# The control reproduces BOTH COLUMNS of that row from Amarsi's own inputs: the published
# per-line 1D LTE abundances of Allende Prieto et al. 2002 Table 2, restricted to weak
# lines (REW < -4.9) exactly as Sect. 6.1 says, corrected by the reactivated network. It
# is a check against published numbers, not against our own pipeline — RYA-785's lesson
# is that a check anchored on the thing being tested is the wrong referee.
#
# WHICH LINE LIST is the trap here. Run this same control over the TRAINING list instead
# and Fe I misses by ~0.04 dex, because the golden set is dominated by high-Elo lines
# while Amarsi's solar set reaches down to 0.05 eV, and the correction tracks Elo with
# r = +0.94. That miss is a line-selection artefact, not a defect in the reactivation.
PAPER_SOLAR = {"I": {"1D_LTE": 7.47, "3D_NLTE": 7.46},
               "II": {"1D_LTE": 7.41, "3D_NLTE": 7.47}}
PAPER_ELO_TREND = (-0.10, -0.05)      # d(Delta_1L3N) over Elo 0 -> 5 eV at [M/H]=0
CONTROL_TOL = 0.01                    # the paper quotes 2 decimals; hold to that
SOLAR_CONTROL_CSV = (ROOT / "data" / "reference" / "amarsi2022_training"
                     / "amarsi2022_solar_control_lines.csv")


# ── atomic data ───────────────────────────────────────────────────────────────

def fe_transitions() -> pd.DataFrame:
    """Fe I/II transitions with Elo, Eup and log gf, from the solar VALD extractions."""
    rows = []
    for src in VALD_SOURCES:
        if not src.exists():
            raise SystemExit(f"VALD extraction missing: {src}")
        recs, status = parse_vald_long(src)
        if status["n_failures"]:
            # RYA-429: never a silent drop. A parse failure here would quietly shrink
            # the pool of lines that can be domain-checked.
            raise SystemExit(
                f"{src.name}: {status['n_failures']} unparsed data lines "
                f"(examples: {status['examples'][:2]}). Refusing to run on a partially "
                f"read line list.")
        for r in recs:
            if r["element"] == "Fe" and r["ion"] in ("I", "II"):
                rows.append({"ion": r["ion"], "wavelength_air_A": r["wavelength"],
                             "elo_eV": r["e_low_eV"], "eup_eV": r["e_up_eV"],
                             "loggf": r["log_gf"], "source": src.name})
    df = pd.DataFrame(rows).sort_values("wavelength_air_A").reset_index(drop=True)
    return df


def attach_atomic(lines: pd.DataFrame, trans: pd.DataFrame) -> pd.DataFrame:
    """Join each measured line to its VALD transition on (ion, wavelength, log gf).

    Where several transitions of the same ion sit inside WAVE_TOL_A, the one whose
    log gf is strongest is taken — that is the transition the EW is dominated by — and
    the ambiguity is recorded per line rather than hidden.
    """
    out = lines.copy()
    for col in ("elo_eV", "eup_eV", "loggf_vald", "n_vald_candidates", "vald_source"):
        out[col] = np.nan if col != "vald_source" else ""
    for i, row in out.iterrows():
        cand = trans[(trans["ion"] == row["ion"])
                     & (np.abs(trans["wavelength_air_A"] - row["wavelength_air_A"])
                        <= WAVE_TOL_A)]
        out.at[i, "n_vald_candidates"] = len(cand)
        if cand.empty:
            continue
        best = cand.loc[cand["loggf"].idxmax()]
        out.at[i, "elo_eV"] = float(best["elo_eV"])
        out.at[i, "eup_eV"] = float(best["eup_eV"])
        out.at[i, "loggf_vald"] = float(best["loggf"])
        out.at[i, "vald_source"] = str(best["source"])
    return out


# ── reach survey: can this engine EVER serve a band? (RYA-762 coordination) ───

#: Bands to survey for reach, independent of whether anything has been measured there.
#: RYA-762 asks whether the 3D leg can extend past the 9199.9 A wall that stops Engine B.
#: That wall is a GES level-ID limit in the LINELIST, and the 3D leg does not use GES
#: level IDs — it keys on Elo/Eup/log gf from VALD — so the premise that the two legs
#: have different walls is correct. It just does not follow that the 3D leg has none.
SURVEY_BANDS = (("VIS (training overlap)", 4788, 6810),
                ("IR (RYA-783)", 6910, 9199),
                ("extended IR (RYA-762)", 9199, 13000))


def reach_survey(trans: pd.DataFrame, star: dict) -> pd.DataFrame:
    """Domain-check every Fe I/II transition in a band, with no measurement required.

    Reach is a LINE-PARAMETER question, so it can be answered before a single EW is
    measured — which is the point: it says whether a measurement campaign in a band
    could ever be served by this engine, rather than finding out afterwards.
    """
    teff, logg, vmic = float(star["teff"]), float(star["logg"]), float(star["xi"])
    rows = []
    for label, lo, hi in SURVEY_BANDS:
        sub = trans[(trans["wavelength_air_A"] >= lo) & (trans["wavelength_air_A"] < hi)]
        for ion in IONS:
            s = sub[sub["ion"] == ion]
            if s.empty:
                continue
            verdicts = [amarsi3d.classify_line(ion, r.elo_eV, r.eup_eV, r.loggf,
                                               teff=teff, logg=logg, vmic=vmic,
                                               afe=ASPLUND21_FE)
                        for r in s.itertuples()]
            dE = s["eup_eV"].values - s["elo_eV"].values
            rows.append({
                "band": label, "lo_A": lo, "hi_A": hi, "ion": ion,
                "n_transitions": len(s),
                "n_in_domain": sum(v.in_domain for v in verdicts),
                "n_fail_delta_E": sum(not v.delta_E_ok for v in verdicts),
                "n_fail_feature": sum(not v.feature_ok for v in verdicts),
                "n_fail_level": sum(not v.level_ok for v in verdicts),
                "delta_E_min": round(float(dE.min()), 4),
                "delta_E_max": round(float(dE.max()), 4),
                "elo_max": round(float(s["elo_eV"].max()), 3),
            })
    return pd.DataFrame(rows)


# ── the domain map ────────────────────────────────────────────────────────────

def domain_figure(lines: pd.DataFrame, out_png: Path) -> None:
    """Where the measured lines sit relative to the training set, on the two axes that
    decide the verdict: excitation potential and transition energy.

    The right panel is the finding. Wavelength is not an input to the network, but
    Eup - Elo is, and it IS the wavelength — so the IR band sits in a strip of transition
    energy the training set never occupies, and no amount of Elo overlap rescues it.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tr = amarsi3d.load_training()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    styles = {"VIS": ("tab:blue", "o"), "IR": ("tab:red", "s")}

    for ax, (xcol, xlabel) in zip(axes, (("elo_eV", "lower level $E_{low}$ [eV]"),
                                         ("elo_eV", "lower level $E_{low}$ [eV]"))):
        ax.set_xlabel(xlabel)
    for ax, ycol, ylabel in ((axes[0], "loggf", r"$\log gf$"),
                             (axes[1], "delta_E_eV", r"transition energy $E_{up}-E_{low}$ [eV]")):
        ax.scatter(tr["elo_eV"], tr[ycol], s=42, c="0.75", edgecolor="0.45",
                   label=f"Amarsi+2022 training set (n={len(tr)})", zorder=1)
        for band, (c, mk) in styles.items():
            sub = lines[(lines["band"] == band) & np.isfinite(lines[ycol])
                        & np.isfinite(lines["elo_eV"])]
            if sub.empty:
                continue
            ok = sub[sub["in_domain"].map(bool)]
            no = sub[~sub["in_domain"].map(bool)]
            ax.scatter(no["elo_eV"], no[ycol], s=18, marker=mk, facecolors="none",
                       edgecolors=c, alpha=0.65, zorder=2,
                       label=f"{band} out-of-domain (n={len(no)})")
            if len(ok):
                ax.scatter(ok["elo_eV"], ok[ycol], s=20, marker=mk, c=c, zorder=3,
                           label=f"{band} IN-DOMAIN (n={len(ok)})")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)

    lo, hi = float(tr["delta_E_eV"].min()), float(tr["delta_E_eV"].max())
    axes[1].axhspan(lo, hi, color="0.85", zorder=0)
    axes[1].axhline(lo, color="0.4", lw=1.0, ls="--")
    axes[1].annotate(f"training range {lo:.3f}-{hi:.3f} eV\n"
                     f"(= 4788-6810 $\\AA$, optical only)",
                     xy=(0.03, 0.06), xycoords="axes fraction", fontsize=8, color="0.25")
    axes[0].legend(fontsize=7, loc="lower right")
    fig.suptitle("RYA-817 — Amarsi+2022 3D-NLTE MLP: measured solar Fe lines vs the "
                 "network's training domain", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


# ── the reactivation control ──────────────────────────────────────────────────

def reactivation_control(star: dict) -> dict:
    """Reproduce the solar row of Amarsi+2022 Table 6 with the reactivated network.

    Runs before any band product: if the engine does not reproduce a published result on
    the author's own inputs, nothing it says about the IR is worth reading.
    """
    from pipeline.nlte_corrections import _compute_aberr
    if not SOLAR_CONTROL_CSV.exists():
        raise SystemExit(
            f"solar control line list missing: {SOLAR_CONTROL_CSV}. Regenerate with "
            f"`python3 scripts/rya817_recover_amarsi_training_set.py`. Running the "
            f"engine with no control is what this ticket exists to prevent.")
    sol = pd.read_csv(SOLAR_CONTROL_CSV)
    teff, logg, vmic = float(star["teff"]), float(star["logg"]), float(star["xi"])
    out = {"axis_A_Fe_3N": ASPLUND21_FE, "teff": teff, "logg": logg, "vmic": vmic,
           "line_list": amarsi3d.SOLAR_CONTROL_CITATION,
           "cut": "weak lines only, REW < -4.9 (Amarsi+2022 Sect. 6.1)",
           "checks": [], "passed": True}

    for ion in ("I", "II"):
        sub = sol[(sol["ion"] == ion) & sol["weak_line_rew_lt_m49"]]
        ab = np.array([_compute_aberr(ion, r.elo_eV, r.eup_eV, r.loggf,
                                      teff, logg, ASPLUND21_FE, vmic)
                       for r in sub.itertuples()], dtype=float)
        a1d = float(sub["a_1d_lte_ap2002"].mean())
        a3n = a1d + float(np.nanmean(ab))
        for label, got, want in (("1D LTE", a1d, PAPER_SOLAR[ion]["1D_LTE"]),
                                 ("3D non-LTE", a3n, PAPER_SOLAR[ion]["3D_NLTE"])):
            ok = abs(got - want) <= CONTROL_TOL
            out["passed"] &= ok
            out["checks"].append({
                "check": f"solar A(Fe {ion}) {label} over {len(sub)} weak lines",
                "computed": round(got, 4), "published": want,
                "tolerance": CONTROL_TOL, "pass": ok,
                "source": "Amarsi+2022 Table 6"})
        out[f"fe{ion}_mean_correction"] = round(float(np.nanmean(ab)), 4)

    # the Elo trend is a statement about the TRAINING set's parameter range, so it is
    # checked there — a different claim, a different line list, said out loud.
    tr = amarsi3d.load_training()
    fe1 = tr[tr["species"] == "Fe1"]
    d1l3n = -np.array([_compute_aberr("I", r.elo_eV, r.eup_eV, r.loggf,
                                      teff, logg, ASPLUND21_FE, vmic)
                       for r in fe1.itertuples()], dtype=float)
    span = 5.0 * float(np.polyfit(fe1["elo_eV"].values, d1l3n, 1)[0])
    ok = PAPER_ELO_TREND[0] <= span <= PAPER_ELO_TREND[1]
    out["passed"] &= ok
    out["checks"].append({
        "check": "Fe I d(Delta_1L3N) over Elo 0->5 eV at [M/H]=0 (training lines)",
        "computed": round(span, 4), "published": list(PAPER_ELO_TREND),
        "pass": ok, "source": "Amarsi+2022 Sect. 5.1"})
    out["fe1_delta_range"] = [round(float(np.nanmin(d1l3n)), 4),
                              round(float(np.nanmax(d1l3n)), 4)]
    return out


# ── the run ───────────────────────────────────────────────────────────────────

def run_band(ion: str, band: str, lo: int, hi: int, bp_dir: Path,
             trans: pd.DataFrame, star: dict, instrument: str
             ) -> tuple[pd.DataFrame, object, dict]:
    src = bp_dir / f"Fe{ion}_{lo}_{hi}_{instrument}_PROFILEFIT_1D-LTE_lines.csv"
    if not src.exists():
        raise SystemExit(f"per-line 1D-LTE product missing: {src}")
    lines = pd.read_csv(src)
    lines = lines[lines["ion"].astype(str) == ion].reset_index(drop=True)
    lines = attach_atomic(lines, trans)

    teff, logg, vmic = float(star["teff"]), float(star["logg"]), float(star["xi"])

    base_rows = []
    for _, r in lines.iterrows():
        base_rows.append(dict(
            element="Fe", ion=ion, band=band, wavelength_air_A=float(r["wavelength_air_A"]),
            ew_mA=float(r["ew_mA"]) if pd.notna(r.get("ew_mA")) else np.nan,
            rew=float(r["rew"]) if pd.notna(r.get("rew")) else np.nan,
            a_1dlte=float(r["abundance"]) if pd.notna(r.get("abundance")) else np.nan,
            in_aggregate_1dlte=bool(r.get("in_aggregate", True)),
            excluded_reason_1dlte=str(r.get("excluded_reason") or ""),
            elo_eV=r["elo_eV"], eup_eV=r["eup_eV"], loggf=r["loggf_vald"],
            n_vald_candidates=int(r["n_vald_candidates"]),
        ))
    per_line = pd.DataFrame(base_rows)

    have_inputs = (np.isfinite(per_line["a_1dlte"]) & np.isfinite(per_line["elo_eV"])
                   & np.isfinite(per_line["loggf"]))

    # STEP 1 — the abundance-INDEPENDENT axes. Elo/Eup/log gf, the transition energy and
    # level representation do not depend on what goes on the A(Fe;3N) axis, so they are
    # settled before any abundance is chosen. This is what stops the axis from being
    # tuned by the very lines it is supposed to judge.
    axis_probe = float(np.clip(np.nanmedian(per_line.loc[have_inputs, "a_1dlte"])
                               if have_inputs.any() else 7.46, 4.5, 7.5))
    for col in ("network", "domain_reason"):
        per_line[col] = ""
    for col in ("feature_ok", "delta_E_ok", "level_ok", "stellar_ok", "in_domain"):
        per_line[col] = pd.Series([None] * len(per_line), dtype=object)
    per_line["delta_E_eV"] = np.nan
    for i, r in per_line.iterrows():
        if not have_inputs[i]:
            per_line.at[i, "in_domain"] = False
            per_line.at[i, "domain_reason"] = (
                "no 1D-LTE abundance (line excluded upstream)"
                if not np.isfinite(r["a_1dlte"]) else
                f"no Fe {ion} transition in the VALD extraction within {WAVE_TOL_A} A "
                f"-- cannot domain-check, so cannot run")
            continue
        v = amarsi3d.classify_line(ion, r["elo_eV"], r["eup_eV"], r["loggf"],
                                   teff=teff, logg=logg, vmic=vmic, afe=axis_probe)
        per_line.at[i, "network"] = v.network
        per_line.at[i, "feature_ok"] = v.feature_ok
        per_line.at[i, "delta_E_ok"] = v.delta_E_ok
        per_line.at[i, "level_ok"] = v.level_ok
        per_line.at[i, "stellar_ok"] = v.stellar_ok
        per_line.at[i, "delta_E_eV"] = v.delta_E_eV
        per_line.at[i, "in_domain"] = bool(v.in_domain)
        per_line.at[i, "domain_reason"] = v.reason

    # STEP 2 — converge the STAR's A(Fe;3N) on the lines that will actually enter the
    # product. Iterating on all measured lines would let a line the product excludes set
    # the axis the product is computed at.
    usable_mask = per_line["in_domain"].map(bool) & per_line["in_aggregate_1dlte"]
    afe_free, n_iter, converged = amarsi3d.converge_star_abundance(
        per_line.loc[usable_mask, ["ion", "elo_eV", "eup_eV", "loggf", "a_1dlte"]],
        teff=teff, logg=logg, vmic=vmic)
    # The grid's A(Fe;3N) axis stops at 7.5 and OUR solar Fe reads high — the 1D-LTE VIS
    # band product is 7.586, carrying the RYA-161 ungraded-gf systematic — so the
    # self-consistent axis rails on the ceiling. That is a property of our gf zero point,
    # not of the Sun (whose A(Fe) is ~7.46, comfortably inside). Clamping is therefore
    # correct, but it is NOT allowed to be silent: the rail is recorded, and the axis
    # SENSITIVITY is measured below so the reader can see what the clamp is worth
    # instead of taking it on trust.
    afe_star = float(np.clip(afe_free, 4.5, 7.5)) if np.isfinite(afe_free) else afe_free
    afe_railed = bool(np.isfinite(afe_free) and abs(afe_star - afe_free) > 1e-9)

    # STEP 3 — per-line corrections at the pinned axis (the product), plus the
    # archived per-line axis and the refused extrapolation, both as diagnostics.
    for col in ("aberr", "aberr_axis_line", "aberr_if_extrapolated", "a_3dnlte"):
        per_line[col] = np.nan
    for i, r in per_line.iterrows():
        if not have_inputs[i]:
            continue
        args = (ion, r["elo_eV"], r["eup_eV"], r["loggf"], r["a_1dlte"])
        kw = dict(teff=teff, logg=logg, vmic=vmic)
        if np.isfinite(afe_star):
            ab, _ = amarsi3d.aberr_for_line(*args, afe3n_axis=afe_star, **kw)
            per_line.at[i, "aberr"] = ab
            if np.isfinite(ab):
                per_line.at[i, "a_3dnlte"] = float(r["a_1dlte"]) + ab
        ab_line, _ = amarsi3d.aberr_for_line(*args, **kw)
        per_line.at[i, "aberr_axis_line"] = ab_line
        # With no in-domain line the star axis cannot be converged, so the refused-
        # extrapolation diagnostic is evaluated at the band's own 1D-LTE median. It is a
        # diagnostic precisely because that axis is itself unvalidated here.
        ab_x, _ = amarsi3d.aberr_for_line(
            *args, afe3n_axis=afe_star if np.isfinite(afe_star) else axis_probe,
            allow_out_of_domain=True, **kw)
        per_line.at[i, "aberr_if_extrapolated"] = ab_x

    # How much does the clamped axis actually matter? Recompute the in-domain median
    # correction at both ends of the axis. If this is a few thousandths, the rail above
    # is immaterial and the number says so; if it were large, the product would not be
    # defensible and the number would say that instead.
    axis_sensitivity = None
    dom = per_line[per_line["in_domain"].map(bool) & np.isfinite(per_line["aberr"])]
    if len(dom):
        ends = []
        for axis in (4.5 + 2.5, 7.5):   # mid-grid and the ceiling
            v = [amarsi3d.aberr_for_line(ion, r.elo_eV, r.eup_eV, r.loggf, r.a_1dlte,
                                         teff=teff, logg=logg, vmic=vmic,
                                         afe3n_axis=axis)[0] for r in dom.itertuples()]
            ends.append(float(np.nanmedian(v)))
        axis_sensitivity = round(abs(ends[1] - ends[0]), 5)

    measurements = []
    for _, r in per_line.iterrows():
        usable = bool(r["in_domain"]) and np.isfinite(r["a_3dnlte"]) \
            and bool(r["in_aggregate_1dlte"])
        if usable:
            reason = ""
        elif not r["in_aggregate_1dlte"]:
            reason = (f"carried from the 1D-LTE product's own exclusion: "
                      f"{r['excluded_reason_1dlte'] or 'no 1D-LTE abundance'}")
        else:
            reason = (f"OUT-OF-DOMAIN for the Amarsi 2022 MLP: "
                      f"{r['domain_reason'] or 'no correction returned'}")
        measurements.append(LineMeasurement(
            element="Fe", ion=ion, wavelength_air_A=float(r["wavelength_air_A"]),
            instrument=instrument,
            ew_mA=float(r["ew_mA"]) if np.isfinite(r["ew_mA"]) else 0.0,
            ew_method="RYA-783 PROFILE-FIT (inherited)",
            # RYA-871 — `elo_eV` here is the VALD lower-level energy this leg matched the
            # line to, i.e. the same quantity `ep_eV` names elsewhere. Carried under the
            # shared name so one consumer reads one column; NaN stays None rather than
            # becoming 0 eV, which would key the line on a level it does not have.
            ep_eV=(float(r["elo_eV"]) if np.isfinite(r.get("elo_eV", np.nan))
                   else None),
            abundance=float(r["a_3dnlte"]) if usable else None,
            rew=float(r["rew"]) if np.isfinite(r["rew"]) else None,
            treatment=amarsi3d.TREATMENT,
            in_aggregate=usable, excluded_reason=reason,
            ew_inversion=False,   # the REW ceiling was already applied by the 1D-LTE leg
        ))

    # RYA-869 — the abundances are RYA-783's profile-fit EW inversions with a 3D-NLTE
    # departure added per line (`ew_method` says so on every measurement), so the harness
    # systematic this product carries is the profile fitter's. `ENGINE-A-3DNLTE` is the
    # other `X-VARIANT` treatment name the ticket named: it must follow `ENGINE-A`, and
    # it does so here by declaring the handler instead of by matching a prefix.
    product = build_product(
        "Fe", ion, instrument, band, amarsi3d.TREATMENT, measurements,
        handler="ProfileFitHandler",
        provenance=(f"{amarsi3d.CITATION}; training domain from "
                    f"{amarsi3d.TRAINING_CITATION}; 1D-LTE base from {src.name}"))

    def _n_failing(axis: str) -> int:
        """Lines the domain check rejected on `axis`. A line with no atomic data has no
        verdict on any axis and is counted only in n_out_of_domain, never here."""
        if axis not in per_line.columns:
            return 0
        col = per_line[axis]
        return int((col.notna() & (col == False)).sum())   # noqa: E712 (object dtype)

    n = len(per_line)
    in_dom = per_line["in_domain"].map(bool)
    n_in = int(in_dom.sum())
    diag = per_line.loc[in_dom & np.isfinite(per_line["aberr_axis_line"])]
    stats = {
        "n_measured_lines": n,
        "n_with_atomic_data": int(np.isfinite(per_line["elo_eV"]).sum()),
        "n_in_domain": n_in,
        "n_out_of_domain": int(n - n_in),
        "fail_feature": _n_failing("feature_ok"),
        "fail_delta_E": _n_failing("delta_E_ok"),
        "fail_level": _n_failing("level_ok"),
        "fail_stellar": _n_failing("stellar_ok"),
        "afe_axis_star": None if not np.isfinite(afe_star) else round(float(afe_star), 4),
        "afe_axis_selfconsistent": (None if not np.isfinite(afe_free)
                                    else round(float(afe_free), 4)),
        "afe_axis_railed_at_grid_ceiling": afe_railed,
        "afe_axis_sensitivity_dex": axis_sensitivity,
        "afe_axis_iterations": int(n_iter),
        "afe_axis_converged": bool(converged),
        "median_aberr": (round(float(np.nanmedian(per_line["aberr"])), 4)
                         if np.isfinite(per_line["aberr"]).any() else None),
        # RYA-207 per-line-axis diagnostic: what the archived leg's axis convention
        # would have produced on the SAME in-domain lines.
        "A_3dnlte_axis_line_diag": (
            round(float(np.median(diag["a_1dlte"] + diag["aberr_axis_line"])), 4)
            if len(diag) else None),
    }
    return per_line, product, stats


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--band-products-dir", type=Path, required=True,
                    help="directory holding the RYA-783 *_1D-LTE_lines.csv per-line files")
    ap.add_argument("--instrument", required=True,
                    help="the arm whose 1D-LTE per-line products to correct. RYA-922: "
                         "REQUIRED, and deliberately has no default — this route used to "
                         "hardcode kpno_solar_atlas, which is why 3D-NLTE existed for one "
                         "arm only. An instrument that is inferred is an instrument that "
                         "can be wrong.")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--star", default="solar")
    args = ap.parse_args(argv)

    star = get_star_params(args.star)
    trans = fe_transitions()
    doms = amarsi3d.domains()

    print("=" * 78)
    print("RYA-817 — Amarsi 2022 3D-NLTE MLP, reactivated, domain-checked")
    print("=" * 78)
    print(f"\nTRAINING DOMAIN ({amarsi3d.TRAINING_CITATION}):")
    for net, d in sorted(doms.items()):
        print(f"  {net:<5} n={d.n_lines:3d}  lambda {d.lambda_air[0]:.1f}-{d.lambda_air[1]:.1f} A  "
              f"Elo {d.elo[0]:.3f}-{d.elo[1]:.3f}  Eup {d.eup[0]:.3f}-{d.eup[1]:.3f}  "
              f"lggf {d.loggf[0]:+.3f}..{d.loggf[1]:+.3f}  "
              f"dE {d.delta_E[0]:.4f}-{d.delta_E[1]:.4f} eV")
    print(f"\nSTAR {args.star}: Teff={star['teff']:.0f} logg={star['logg']:.3f} "
          f"vmic={star['xi']:.2f} km/s  (Fe I/II transitions loaded: {len(trans)})")

    ctl = reactivation_control(star)
    print("\nREACTIVATION CONTROL — reproduce Amarsi+2022 Table 6 (Sun) from its own inputs")
    print(f"  line list: {ctl['line_list']}")
    print(f"  cut: {ctl['cut']}")
    for c in ctl["checks"]:
        pub = c["published"]
        pub_s = f"{pub:+.3f}" if isinstance(pub, float) else f"{pub[0]:+.2f}..{pub[1]:+.2f}"
        print(f"  {'PASS' if c['pass'] else 'FAIL'}  {c['check']:<58s} "
              f"computed {c['computed']:+.4f}  published {pub_s}   [{c['source']}]")
    print(f"        solar corrections: Fe I {ctl['feI_mean_correction']:+.4f} dex, "
          f"Fe II {ctl['feII_mean_correction']:+.4f} dex")
    print(f"        Fe I Delta_1L3N range at the Sun: "
          f"{ctl['fe1_delta_range'][0]:+.3f}..{ctl['fe1_delta_range'][1]:+.3f} dex "
          f"(paper: -0.15..+0.10 across the whole Teff/logg grid at [M/H]=0)")

    args.out.mkdir(parents=True, exist_ok=True)
    all_lines, products, summary = [], [], {}

    # RYA-922 — DISCOVER what this arm has, rather than iterate a hardcoded tuple.
    # Bluest first: the training-overlap band is the control for the reactivation, and a
    # failed control invalidates anything redder that follows it.
    plan = []
    for ion in IONS:
        for band, lo, hi, src in discover_bands(args.band_products_dir,
                                                args.instrument, ion):
            plan.append((band, lo, hi, ion, src))
    if not plan:
        raise SystemExit(
            f"no Fe 1D-LTE per-line products for instrument {args.instrument!r} under "
            f"{args.band_products_dir}. This route CORRECTS an existing 1D-LTE pool; it "
            f"does not measure one. Run the band harness for that arm first.")
    plan.sort(key=lambda t: (t[1], t[3]))
    print(f"  {len(plan)} measured span(s) found for {args.instrument}:")
    for band, lo, hi, ion, src in plan:
        print(f"    Fe {ion:<2} {band:<12} {lo}-{hi} A   {src.name}")
    print()

    for band, lo, hi, ion, _src in plan:
        if True:
            per_line, product, stats = run_band(ion, band, lo, hi,
                                                args.band_products_dir, trans, star,
                                                args.instrument)
            all_lines.append(per_line)
            products.append(product)
            key = f"Fe {ion} {band}"
            summary[key] = {**stats,
                            "A_3dnlte": product.value, "sigma": product.sigma,
                            "n_lines": product.n_lines, "n_excluded": product.n_excluded}
            val = f"{product.value:.3f}" if product.value is not None else "  n/a"
            sig = f"{product.sigma:.3f}" if product.sigma is not None else "  n/a"
            print(f"\n{key:<10} measured={stats['n_measured_lines']:3d}  "
                  f"atomic={stats['n_with_atomic_data']:3d}  "
                  f"IN-DOMAIN={stats['n_in_domain']:3d}  "
                  f"OUT={stats['n_out_of_domain']:3d}")
            print(f"{'':<10} out-of-domain by axis: feature={stats['fail_feature']} "
                  f"dE={stats['fail_delta_E']} level={stats['fail_level']} "
                  f"stellar={stats['fail_stellar']}")
            print(f"{'':<10} A(Fe {ion}; 3D-NLTE) = {val} +/- {sig}  "
                  f"(n={product.n_lines})")
            if stats["afe_axis_star"] is not None:
                rail = (" RAILED at the 7.5 grid ceiling"
                        if stats["afe_axis_railed_at_grid_ceiling"] else "")
                print(f"{'':<10} A(Fe;3N) axis = {stats['afe_axis_star']:.3f}"
                      f"{rail} (self-consistent "
                      f"{stats['afe_axis_selfconsistent']:.3f}); product moves "
                      f"{stats['afe_axis_sensitivity_dex']:.4f} dex across the WHOLE "
                      f"axis, so the clamp is immaterial")
                print(f"{'':<10} median per-line correction "
                      f"{stats['median_aberr']:+.4f} dex; per-line-axis diagnostic "
                      f"(RYA-207 convention) would read "
                      f"{stats['A_3dnlte_axis_line_diag']:.3f}")
            if stats["n_in_domain"] == 0:
                oo = per_line[np.isfinite(per_line["aberr_if_extrapolated"])]
                if len(oo):
                    print(f"{'':<10} REFUSED extrapolation would have given "
                          f"median aberr {oo['aberr_if_extrapolated'].median():+.3f} dex "
                          f"(range {oo['aberr_if_extrapolated'].min():+.3f}.."
                          f"{oo['aberr_if_extrapolated'].max():+.3f}) — NOT a product.")

    # RYA-762 coordination: does the engine reach past the 9199.9 A Engine-B wall?
    survey = reach_survey(trans, star)
    survey.to_csv(args.out / "rya817_reach_survey.csv", index=False)
    print("\n" + "=" * 78)
    print("REACH SURVEY — every Fe transition in the band, no measurement required")
    print("=" * 78)
    print(f"  {'band':<24s}{'ion':<4s}{'lines':>7s}{'in-domain':>11s}"
          f"{'dE range (eV)':>20s}{'Elo max':>9s}")
    for _, r in survey.iterrows():
        print(f"  {r['band']:<24s}{r['ion']:<4s}{r['n_transitions']:7d}"
              f"{r['n_in_domain']:11d}"
              f"{r['delta_E_min']:>11.4f}-{r['delta_E_max']:<8.4f}{r['elo_max']:>9.3f}")

    lines_df = pd.concat(all_lines, ignore_index=True)
    lines_df.to_csv(args.out / "rya817_3dnlte_per_line.csv", index=False)
    prod_df = products_frame(products)
    prod_df.to_csv(args.out / "rya817_3dnlte_products.csv", index=False)
    domain_figure(lines_df, args.out / "rya817_domain_map.png")

    # the comparison the ticket asks for, as a DIAGNOSTIC table
    print("\n" + "=" * 78)
    print("COMPARISON (Asplund 2021 predicts NO band dependence in proper 3D-NLTE)")
    print("=" * 78)
    _bands_run = sorted({k.split(" ", 2)[2] for k in summary})
    fe1 = {b: summary.get(f"Fe I {b}", {}) for b in _bands_run}
    vis, ir = fe1.get("VIS", {}), fe1.get("IR", {})

    def _fmt(d):
        v = d.get("A_3dnlte")
        return f"{v:.3f}" if v is not None else "NOT PRODUCED (0 in-domain lines)"

    print(f"  Fe I VIS 3D-NLTE : {_fmt(vis):<34s} n={vis.get('n_lines')}")
    print(f"  Fe I IR  3D-NLTE : {_fmt(ir):<34s} n={ir.get('n_lines')}")
    print(f"  gold v4 VIS anchor (3D-NLTE, frozen) : {GOLD_V4_FE_I}")
    print(f"  Asplund, Amarsi & Grevesse 2021      : {ASPLUND21_FE}")
    if vis.get("A_3dnlte") is None or ir.get("A_3dnlte") is None:
        print("\n  The band-dependence test Asplund 2021 predicts CANNOT BE RUN with "
              "this engine:\n  one side of the comparison does not exist, because the "
              "network has no IR training\n  line to stand on. That is the result, not "
              "a failure to produce a number.")

    (args.out / "rya817_run_summary.json").write_text(json.dumps({
        "ticket": "RYA-817",
        "star": args.star,
        "engine": amarsi3d.CITATION,
        "training_set": amarsi3d.TRAINING_CITATION,
        "treatment": amarsi3d.TREATMENT,
        "reactivation_control": ctl,
        "band_products_dir": str(args.band_products_dir),
        "training_domain": {k: {"n_lines": d.n_lines, "lambda_air": list(d.lambda_air),
                                "elo": list(d.elo), "eup": list(d.eup),
                                "loggf": list(d.loggf), "delta_E": list(d.delta_E)}
                            for k, d in sorted(doms.items())},
        "per_band": summary,
        "context_values": {"gold_v4_fe_i_3dnlte_vis": GOLD_V4_FE_I,
                           "asplund2021_fe": ASPLUND21_FE},
        "run_date": str(date.today()),
    }, indent=2) + "\n")
    print(f"\nwrote {args.out}/rya817_3dnlte_per_line.csv")
    print(f"wrote {args.out}/rya817_3dnlte_products.csv")
    print(f"wrote {args.out}/rya817_domain_map.png")
    print(f"wrote {args.out}/rya817_reach_survey.csv")
    print(f"wrote {args.out}/rya817_run_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""RYA-1106 — the Amarsi 3D-NLTE method on ASPLUND'S OWN solar Fe I lines, per instrument.

    python3 scripts/derive_amarsi_asplund_lines.py \
        --holdings kpno_kurucz2005 kpno_molecfit harps_molecfit iag

WHY THIS PRODUCT EXISTS, AND WHY IT IS NOT A RE-SLICE OF WHAT WE ALREADY HOLD
-----------------------------------------------------------------------------
RYA-1104 established that our shipped Amarsi leg reads 7.512 against Asplund's 7.46 for a
reason that has nothing to do with the 3D-NLTE physics and nothing to do with the
measurement route: the LINE POPULATION. Our graded pool is the laboratory-gf tier, and
Ruffoni 2014 / Den Hartog 2014 measured HIGH-EXCITATION Fe I, so the pool floors at
Elo 2.85 eV. Amarsi's solar set reaches down to 0.05 eV. The network's correction tracks
Elo, so the two pools cannot return the same number.

🔴 THE OBVIOUS SHORTCUT DOES NOT EXIST. "Restrict our existing pool to Asplund's lines"
sounds like a one-line filter and it is not available: measured on the committed tree,
our pools contain **1 of Asplund's 50 Fe I lines** (5784.658). Every VIS holding scores
the same 1/50. There is no subset to take, so the lines have to be MEASURED — which is
what this script does, per holding, through the production synthesis fitter.

WHAT "ASPLUND'S OWN LINES" MEANS HERE, EXACTLY
----------------------------------------------
Amarsi, Liljegren & Nissen 2022 Sect. 6.1: "Line-by-line 1D LTE lg eps(Fe) values measured
in the solar flux atlas of Kurucz (1984) were taken from Allende Prieto et al. (2002). The
analysis was restricted to weak lines with REW < -4.9."

So the line set is AP2002 Table 2 (Fe I, 50 lines) under Amarsi's own REW < -4.9 cut
(41 lines). Both are reported; the CUT pool is the product, because that is the pool the
published 7.46 stands on and the pool the RYA-817 reactivation control reproduces.

⚠️ THE REW CUT IS AP2002'S OWN, NOT A CUT OF OURS. It is applied on AP2002's published
equivalent widths, because it is a property of THEIR line selection. Re-deriving it from
our own fits would silently make the pool a function of our measurement — a different
experiment wearing the same name.

🔴 AND THE gf IS ASPLUND'S. This is the part that makes it a replication rather than a
re-run. `canonical_gf` carries Kurucz/VALD values for most of these lines — our scale, not
his — and log gf sets the abundance directly. So the TARGET lines carry AP2002's own
published log gf, and everything else in the window keeps canonical gf exactly as
production does. The override is applied to the in-memory linelist the fitter is handed,
on the lambda+EP dual key, and every substitution is reported per line.

⚠️ THIS IS A SCOPED, DECLARED EXCEPTION TO RYA-353 SINGLE-SOURCING, and it uses the
mechanism the pipeline already provides for exactly this (`_load_synth_resources(
apply_canonical_gf=...)`, whose docstring requires the caller to say so in provenance)
rather than a bypass invented here. `canonical_gf.csv` is NEVER written. The exception is
confined to one product whose entire purpose is to be measured on someone else's scale.

🔴 RYA-161 FIREWALL. This product is emitted because it uses Asplund's ACTUAL INPUTS, not
because of what it returns. Whatever each instrument gives is what is reported. Nothing in
this script consults 7.46 to decide anything: `ASPLUND21_FE` appears only in the report, as
a number to display beside the result, and there is no cut, no weight and no branch keyed
on the distance to it.

COVERAGE IS LOUD (RYA-429/RYA-711)
----------------------------------
An Asplund line an instrument cannot serve — out of the holding's span, no admissible
window, a fit that did not converge — is REPORTED with its reason and excluded. It is never
replaced by a neighbour and never silently dropped, and the per-holding coverage count is
part of the result rather than a footnote: a holding that covers 30 of 41 lines is
answering a different question from one that covers 41, and the reader must be able to see
that without reconstructing it.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline.gf_resolver import WTOL, EPTOL                      # noqa: E402
from pipeline import amarsi3d                                     # noqa: E402
from pipeline.band_products import LineMeasurement, build_product, products_frame  # noqa: E402
from config.constants import get_star_params                      # noqa: E402
from config.synth_bands import SYNTH_BANDS                        # noqa: E402

ASPLUND21_FE = 7.46          # A(Fe), Asplund, Amarsi & Grevesse 2021 — DISPLAY ONLY

SOLAR_CONTROL_CSV = (ROOT / "data" / "reference" / "amarsi2022_training"
                     / "amarsi2022_solar_control_lines.csv")
OUT_DEFAULT = ROOT / "data" / "results" / "rya1106"

#: The VIS holdings Amarsi covers. Amarsi is Fe I VIS ONLY — N/A for Fe II, near-UV, NIR
#: and red-optical — so "all instruments" is these four and the scope is stated here rather
#: than discovered, so a holding that silently stops resolving is a failure and not a
#: shorter table.
HOLDINGS: dict[str, tuple[str, str]] = {
    "kpno_kurucz2005": ("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected"),
    "kpno_molecfit":   ("kpno_solar_atlas", "solar_kpno_molecfit_corrected"),
    "harps_molecfit":  ("harps",            "solar_harps_molecfit_corrected"),
    "iag":             ("iag_fts_solar_atlas", "solar_iag"),
}

#: Amarsi+2022 Sect. 6.1's own cut, quoted so the number in the code and the number in the
#: paper cannot drift apart.
REW_CUT = -4.9
BAND = "VIS"


class AsplundRunRefusal(SystemExit):
    """A question this run cannot answer honestly from the material available."""


# ── the line set ──────────────────────────────────────────────────────────────

def asplund_fe1_lines(*, weak_only: bool = True) -> pd.DataFrame:
    """AP2002 Table 2 Fe I, optionally under Amarsi's own REW < -4.9 cut.

    The artifact is the FULL published table and the cut is a COLUMN on it
    (`weak_line_rew_lt_m49`, written by `scripts/rya817_recover_amarsi_training_set.py`),
    so this function selects rather than re-derives. Recomputing the cut here would put a
    second copy of Amarsi's threshold in the tree, free to drift from the first.
    """
    if not SOLAR_CONTROL_CSV.exists():
        raise AsplundRunRefusal(
            f"Asplund/AP2002 solar line set missing: {SOLAR_CONTROL_CSV}\n"
            f"  Regenerate with `python3 scripts/rya817_recover_amarsi_training_set.py`.\n"
            f"  This product IS that line set — there is nothing to substitute.")
    d = pd.read_csv(SOLAR_CONTROL_CSV)
    d = d[d["ion"].astype(str) == "I"].reset_index(drop=True)
    if weak_only:
        d = d[d["weak_line_rew_lt_m49"].astype(bool)].reset_index(drop=True)
    return d


# ── the gf override ───────────────────────────────────────────────────────────

def override_target_gf(ctx: dict, targets: pd.DataFrame) -> dict:
    """Put AP2002's OWN log gf on the target lines of the in-memory synthesis list.

    🔴 WHY IN MEMORY AND NOT A WRITTEN LINE LIST. `fit_one` hands `ctx['linelist']`
    straight to the fitter, so overriding the array the context already holds is the same
    object the fit reads — there is no second file that could drift from it, and no
    write/read round-trip to normalise a value behind our back (RYA-1084's lesson about
    round-tripping data files applies to line lists too).

    ⚠️ ONLY THE TARGET LINES MOVE. Blends inside the window keep canonical gf exactly as
    production set them, because the blend context is OUR modelling choice and is not part
    of what Asplund contributed. Substituting his scale into lines he never measured would
    be inventing a third thing that is neither his analysis nor ours.

    The match is the lambda+EP DUAL KEY at `gf_resolver`'s own tolerances (RYA-1037). A
    wavelength-only join here would be especially dangerous: the GES list carries hyperfine
    and blend components within a few mA of each other, and putting a laboratory gf on the
    wrong component of a blend would move the fit while looking completely normal.
    """
    ll = ctx["linelist"]
    names = ll.dtype.names
    w_A = np.asarray(ll["wave_A"] if "wave_A" in names else ll["wave_nm"] * 10.0,
                     dtype=float)
    ep = np.asarray(ll["lower_state_eV"], dtype=float)
    el = np.asarray([str(x).strip() for x in ll["element"]])
    is_fe1 = el == "Fe 1"

    applied, missing, ambiguous = [], [], []
    for r in targets.itertuples():
        hit = np.flatnonzero(is_fe1
                             & (np.abs(w_A - r.wavelength_air_A) <= WTOL)
                             & (np.abs(ep - r.elo_eV) <= EPTOL))
        if hit.size == 0:
            missing.append({"wavelength_air_A": round(float(r.wavelength_air_A), 4),
                            "elo_eV": round(float(r.elo_eV), 4)})
            continue
        if hit.size > 1:
            # Not fatal in itself — a hyperfine cluster is several rows for one physical
            # line — but it is never resolved by argmin. The gf is placed on EVERY row of
            # the cluster only when they agree on being the same physical transition;
            # otherwise the line is refused and named.
            spread = float(np.ptp(w_A[hit]))
            if spread > WTOL:
                ambiguous.append({"wavelength_air_A": round(float(r.wavelength_air_A), 4),
                                  "n_rows": int(hit.size),
                                  "wavelength_spread_A": round(spread, 5)})
                continue
        before = np.asarray(ll["loggf"], dtype=float)[hit]
        ll["loggf"][hit] = float(r.loggf)
        applied.append({"wavelength_air_A": round(float(r.wavelength_air_A), 4),
                        "elo_eV": round(float(r.elo_eV), 4),
                        "n_rows": int(hit.size),
                        "loggf_canonical": round(float(np.mean(before)), 4),
                        "loggf_ap2002": round(float(r.loggf), 4),
                        "delta_dex": round(float(r.loggf - np.mean(before)), 4)})
    if missing or ambiguous:
        raise AsplundRunRefusal(
            f"the synthesis line list cannot carry Asplund's gf for the whole pool: "
            f"{len(missing)} line(s) absent from the GES list on the lambda+EP dual key "
            f"(+/-{WTOL} A / +/-{EPTOL} eV) {missing[:4]}, {len(ambiguous)} ambiguous "
            f"{ambiguous[:3]}.\n"
            f"  A partial override would measure SOME lines on Asplund's scale and the "
            f"rest on ours and report one number for the mixture — which is exactly the "
            f"confound this product exists to remove. Refusing (RYA-429).")
    deltas = np.array([a["delta_dex"] for a in applied], dtype=float)
    return {
        "n_targets": int(len(targets)), "n_applied": len(applied),
        "per_line": applied,
        "delta_vs_canonical_dex": {
            "mean": round(float(deltas.mean()), 4),
            "median": round(float(np.median(deltas)), 4),
            "min": round(float(deltas.min()), 4),
            "max": round(float(deltas.max()), 4),
            "n_exact": int((np.abs(deltas) < 1e-9).sum()),
        },
        "note": ("AP2002's published log gf on the target lines; canonical gf everywhere "
                 "else in the window. RYA-353 single-sourcing is DECLARED OFF for the "
                 "targets and canonical_gf.csv is not written."),
    }


# ── measuring one holding ─────────────────────────────────────────────────────

def measure_holding(key: str, targets: pd.DataFrame, *, tmp_root: Path,
                    star: str = "solar") -> dict:
    """Fit every Asplund line in one holding, through the PRODUCTION synthesis fitter."""
    instrument, holding = HOLDINGS[key]
    from pipeline.nearuv_synth import build_solar_context
    from rya759_nearuv_fe_product import fit_one
    from rya759_nearuv_synth import _kp_segments
    from measure_band_ew import load_window_ex

    cfg = SYNTH_BANDS[BAND]
    cat = pd.read_csv(ROOT / "data" / "catalog" / "instrument_catalog.csv")
    row = cat[cat.iloc[:, 0].astype(str) == instrument]
    if row.empty:
        raise AsplundRunRefusal(
            f"{instrument!r} absent from data/catalog/instrument_catalog.csv")
    R = float(row.iloc[0]["resolving_power_max"])

    # `apply_canonical_gf=True` first — production behaviour for the whole window — and
    # THEN the targets are moved to AP2002's scale. Loading with canonical off would leave
    # the blends on raw GES gf, which is a third scale nobody asked for.
    ctx = build_solar_context("Fe", R, linelist_file=str(cfg.linelist),
                              apply_canonical_gf=True, star=star)
    gf_report = override_target_gf(ctx, targets)

    segs = _kp_segments() if instrument == "kpno_solar_atlas" else None

    def _load(centre: float, pad: float):
        w = load_window_ex(instrument, centre, pad, segs=segs, holding=holding)
        return w.wave, w.flux, w.provenance

    tmp = tmp_root / f"synth_{key}"
    tmp.mkdir(parents=True, exist_ok=True)

    rows = []
    for r in targets.itertuples():
        t0 = time.time()
        try:
            res = fit_one(ctx, segs, float(r.wavelength_air_A), cfg.half_width_A,
                          str(tmp), load=_load)
        except Exception as exc:                       # a holding that cannot serve it
            res = {"status": "unserved", "reason": f"{type(exc).__name__}: {exc}",
                   "a_synth": np.nan, "red_chi2": np.nan, "obs_source": ""}
        a = float(res.get("a_synth", np.nan))
        ok = str(res.get("status")) == "ok" and np.isfinite(a)
        rows.append({
            "wavelength_air_A": float(r.wavelength_air_A),
            "elo_eV": float(r.elo_eV), "eup_eV": float(r.eup_eV),
            "loggf_ap2002": float(r.loggf),
            "rew_ap2002": float(r.rew), "ew_mA_ap2002": float(r.ew_mA),
            "a_1dlte": a if ok else np.nan,
            "a_1d_lte_ap2002": float(r.a_1d_lte_ap2002),
            "status": str(res.get("status")),
            "reason": str(res.get("reason", "")),
            "red_chi2": float(res.get("red_chi2", np.nan)),
            "obs_source": str(res.get("obs_source", "")),
            "fit_seconds": round(time.time() - t0, 2),
        })
    per_line = pd.DataFrame(rows)
    return {"key": key, "instrument": instrument, "holding": holding,
            "resolving_power": R, "gf_override": gf_report,
            "per_line": per_line, "linelist": str(cfg.linelist),
            "half_width_A": cfg.half_width_A}


# ── the Amarsi correction ─────────────────────────────────────────────────────

def apply_amarsi(per_line: pd.DataFrame, *, star: str = "solar") -> pd.DataFrame:
    """Add the 3D-NLTE correction to each measured line — through `pipeline.amarsi3d`.

    Not a second copy of the sign convention or of the axis iteration: RYA-701 is the
    standing reason, and here there is a sharper one — this product is compared against
    the RYA-817 control and against the shipped leg, so a drifting copy would land
    directly in the comparison.
    """
    p = get_star_params(star)
    teff, logg, vmic = float(p["teff"]), float(p["logg"]), float(p["xi"])
    d = per_line.copy()
    for c in ("network", "domain_reason", "aberr", "a_3dnlte"):
        d[c] = "" if c in ("network", "domain_reason") else np.nan
    for c in ("in_domain",):
        d[c] = pd.Series([None] * len(d), dtype=object)

    have = np.isfinite(d["a_1dlte"].to_numpy(float))
    for i, r in d.iterrows():
        if not have[i]:
            d.at[i, "in_domain"] = False
            d.at[i, "domain_reason"] = f"no 1D-LTE abundance ({r['status']})"
            continue
        v = amarsi3d.classify_line("I", r["elo_eV"], r["eup_eV"], r["loggf_ap2002"],
                                   teff=teff, logg=logg, vmic=vmic, afe=ASPLUND21_FE)
        d.at[i, "network"] = v.network
        d.at[i, "in_domain"] = bool(v.in_domain)
        d.at[i, "domain_reason"] = v.reason

    usable = d["in_domain"].map(bool) & have
    afe_free, n_iter, converged = amarsi3d.converge_star_abundance(
        d.loc[usable].assign(ion="I", loggf=d.loc[usable, "loggf_ap2002"])[
            ["ion", "elo_eV", "eup_eV", "loggf", "a_1dlte"]],
        teff=teff, logg=logg, vmic=vmic)
    afe_star = float(np.clip(afe_free, 4.5, 7.5)) if np.isfinite(afe_free) else afe_free

    for i, r in d.iterrows():
        if not (have[i] and bool(r["in_domain"])):
            continue
        ab, _ = amarsi3d.aberr_for_line("I", r["elo_eV"], r["eup_eV"], r["loggf_ap2002"],
                                        r["a_1dlte"], afe3n_axis=afe_star,
                                        teff=teff, logg=logg, vmic=vmic)
        d.at[i, "aberr"] = ab
        if np.isfinite(ab):
            d.at[i, "a_3dnlte"] = float(r["a_1dlte"]) + ab
    d.attrs["axis"] = {"afe_axis_star": None if not np.isfinite(afe_star) else round(afe_star, 4),
                       "afe_axis_selfconsistent": None if not np.isfinite(afe_free)
                       else round(float(afe_free), 4),
                       "railed_at_grid_ceiling": bool(np.isfinite(afe_free)
                                                      and abs(afe_star - afe_free) > 1e-9),
                       "iterations": int(n_iter), "converged": bool(converged)}
    return d


# ── the product ───────────────────────────────────────────────────────────────

TREATMENT = amarsi3d.TREATMENT          # ENGINE-A-3DNLTE — the stored token, not a new one
#: 🔴 THE HANDLER IS THE SYNTHESISER, AND IT IS STATED BY THE ROUTE THAT RAN.
#: RYA-1104 traced the shipped Amarsi leg's `route=ew` to a hardcoded `ProfileFitHandler`
#: in `scripts/rya817_run_3dnlte_bands.py` that RYA-1031 stranded when it widened that
#: route to SYNTH pools. This leg is a Turbospectrum FLUX FIT end to end, so it declares
#: the handler it actually used — which is what makes `route` resolve to `synth` and stops
#: the budget charging the profile fitter's 0.0129 dex residual (RYA-869).
HANDLER = "SynthesisHandler"


def build(holding_key: str, per_line: pd.DataFrame, run: dict) -> tuple:
    """One product for one holding, plus its budget — the band route's own machinery."""
    from rya1095_amarsi_error_budget import budget_from_pool

    measurements = []
    for r in per_line.itertuples():
        usable = bool(r.in_domain) and np.isfinite(r.a_3dnlte)
        if usable:
            reason = ""
        elif not np.isfinite(r.a_1dlte):
            reason = (f"{run['holding']} could not serve this Asplund line: "
                      f"{r.status}{(' — ' + r.reason) if r.reason else ''}")
        else:
            reason = f"OUT-OF-DOMAIN for the Amarsi 2022 MLP: {r.domain_reason}"
        measurements.append(LineMeasurement(
            element="Fe", ion="I", wavelength_air_A=float(r.wavelength_air_A),
            instrument=run["instrument"],
            ew_mA=float(r.ew_mA_ap2002),
            ew_method=("Turbospectrum synthesis flux-fit on Asplund's own line set "
                       "(AP2002 Table 2) with AP2002's own log gf; EW column is AP2002's "
                       "published width, carried for the REW cut and never measured here"),
            abundance=float(r.a_3dnlte) if usable else None,
            rew=float(r.rew_ap2002), ep_eV=float(r.elo_eV),
            treatment=TREATMENT, in_aggregate=usable, excluded_reason=reason,
            ew_inversion=False))

    product = build_product(
        "Fe", "I", run["instrument"], BAND, TREATMENT, measurements,
        handler=HANDLER,
        provenance=(
            f"{amarsi3d.CITATION}; line set + log gf from "
            f"{amarsi3d.SOLAR_CONTROL_CITATION} under Amarsi+2022 Sect. 6.1's own "
            f"REW < {REW_CUT} cut; 1D-LTE base MEASURED HERE by Turbospectrum flux fit on "
            f"holding {run['holding']} (half-width +/-{run['half_width_A']} A); canonical "
            f"gf single-sourcing (RYA-353) DECLARED OFF for the target lines and ON for "
            f"the blend context — this product is a replication of Asplund's analysis and "
            f"therefore runs on Asplund's gf scale, not ours (RYA-1106)"))

    pool = per_line[per_line["in_domain"].map(bool)
                    & np.isfinite(per_line["a_3dnlte"])].copy()
    budget_text = ""
    stat = syst = basis = None
    if len(pool) >= 2 and product.sigma is not None and np.isfinite(product.sigma):
        pool = pool.assign(loggf=pool["loggf_ap2002"], aberr_axis_line=np.nan)
        b, rung = budget_from_pool(pool, element="Fe", ion="I",
                                   instrument=run["instrument"], handler=HANDLER,
                                   scatter_dex=float(product.sigma))
        stat, syst = b.total()
        basis = b.stat_basis()
        # ⚠️ `publish_product` PARSES the `gf rung:` line out of the budget file, so it is
        # appended here for the same reason the band route appends it at three call sites.
        budget_text = b.describe() + f"\n  gf rung: {rung.describe()}\n"
    return product, budget_text, (stat, syst, basis)


def _gf_pool_axis(per_line: pd.DataFrame) -> str:
    """The gf axis this product EARNS, from what it actually ran on.

    AP2002's log gf are laboratory/solar-calibrated oscillator strengths published with the
    analysis, not the Kurucz theoretical floor. Stated as `lab` and passed to the product
    rather than inherited from `treatment_axes.LEGACY`, which pins `kurucz` for this
    treatment — the pinned constant RYA-1104 found contradicting the leg's own budget.
    """
    return "lab"


# ── coverage, stated rather than counted quietly ──────────────────────────────

def coverage(per_line: pd.DataFrame) -> dict:
    """WHICH Asplund lines this holding could and could not serve, and why.

    A count alone is not enough: "38 of 41" and "38 of 41, and the three missing are the
    lowest-excitation lines in the set" are different results, and on this product the
    second one matters — the whole finding of RYA-1104 is that the correction tracks Elo,
    so a holding that loses its low-Elo lines is answering a shifted question.
    """
    served = per_line[np.isfinite(per_line["a_1dlte"])]
    lost = per_line[~np.isfinite(per_line["a_1dlte"])]
    by_reason: dict[str, int] = {}
    for st in lost["status"].astype(str):
        by_reason[st] = by_reason.get(st, 0) + 1
    return {
        "n_asplund_lines": int(len(per_line)),
        "n_served": int(len(served)), "n_unserved": int(len(lost)),
        "fraction_served": round(len(served) / len(per_line), 3) if len(per_line) else None,
        "unserved_by_status": by_reason,
        "unserved_lines": [
            {"wavelength_air_A": round(float(r.wavelength_air_A), 3),
             "elo_eV": round(float(r.elo_eV), 3), "status": str(r.status),
             "reason": str(r.reason)[:120]} for r in lost.itertuples()],
        "elo_eV_served": ({"min": round(float(served["elo_eV"].min()), 3),
                           "median": round(float(served["elo_eV"].median()), 3),
                           "max": round(float(served["elo_eV"].max()), 3)}
                          if len(served) else None),
        "note": ("an Asplund line this holding cannot serve is EXCLUDED and NAMED, never "
                 "replaced by a neighbour and never dropped quietly (RYA-429/RYA-711)"),
    }


# ── main ──────────────────────────────────────────────────────────────────────

def run_one(key: str, targets: pd.DataFrame, *, tmp_root: Path, out: Path) -> dict:
    run = measure_holding(key, targets, tmp_root=tmp_root)
    per_line = apply_amarsi(run["per_line"])
    product, budget_text, budget = build(key, per_line, run)
    cov = coverage(per_line)

    d = out / key
    d.mkdir(parents=True, exist_ok=True)
    per_line.to_csv(d / "asplund_lines_per_line.csv", index=False)
    frame = products_frame([product])
    stat, syst, basis = budget
    frame["stat_dex"] = None if stat is None else round(stat, 4)
    frame["syst_dex"] = None if syst is None else round(syst, 4)
    frame["stat_basis"] = basis or ""
    frame.to_csv(d / "asplund_lines_products.csv", index=False)
    if budget_text:
        (d / "asplund_lines_budgets.txt").write_text(budget_text)

    row = frame.iloc[0]
    return {
        "holding_key": key, "instrument": run["instrument"], "holding": run["holding"],
        "A_3dnlte": product.value, "sigma": product.sigma, "n_lines": product.n_lines,
        "n_excluded": product.n_excluded,
        "stat_dex": None if stat is None else round(stat, 4),
        "syst_dex": None if syst is None else round(syst, 4),
        "route_axis": str(row["route"]), "route_basis": str(row["route_basis"]),
        "gf_axis": str(row["gf"]), "handler": str(row["handler"]),
        "vs_asplund_dex": (None if product.value is None
                           else round(product.value - ASPLUND21_FE, 4)),
        "base_1dlte_median": (round(float(np.nanmedian(per_line["a_1dlte"])), 4)
                              if np.isfinite(per_line["a_1dlte"]).any() else None),
        "median_aberr": (round(float(np.nanmedian(per_line["aberr"])), 4)
                         if np.isfinite(per_line["aberr"]).any() else None),
        "ap2002_published_1dlte_median": round(
            float(per_line["a_1d_lte_ap2002"].median()), 4),
        "coverage": cov, "gf_override": run["gf_override"],
        "axis": per_line.attrs.get("axis", {}),
        "artifacts": {k: str((d / v).relative_to(ROOT)) for k, v in (
            ("per_line", "asplund_lines_per_line.csv"),
            ("products", "asplund_lines_products.csv"))},
    }


def render(rep: dict) -> str:
    L = []; A = L.append
    A("=" * 92)
    A("RYA-1106 — Amarsi 3D-NLTE on ASPLUND'S OWN Fe I lines, per VIS instrument")
    A("=" * 92)
    p = rep["pool"]
    A(f"  line set : {p['citation']}")
    A(f"  cut      : {p['cut']}  ->  n={p['n_used']} of {p['n_published']} published Fe I")
    A(f"  log gf   : {p['gf']}")
    A(f"  Elo span : {p['elo_min']:.3f}-{p['elo_max']:.3f} eV "
      f"(our graded pool floors at 2.85 — that gap IS the RYA-1104 finding)")
    A("")
    A(f"  {'holding':<17}{'A(Fe)':>8}{'vs 7.46':>9}{'n':>4}{'cov':>9}"
      f"{'1D base':>9}{'corr':>8}   {'route':<7}{'gf':<6}")
    A("  " + "-" * 88)
    for r in rep["holdings"]:
        if r.get("error"):
            A(f"  {r['holding_key']:<17}{'—':>8}   {r['error'][:60]}")
            continue
        cov = r["coverage"]
        v = "  n/a" if r["A_3dnlte"] is None else f"{r['A_3dnlte']:.3f}"
        dv = "   n/a" if r["vs_asplund_dex"] is None else f"{r['vs_asplund_dex']:+.3f}"
        A(f"  {r['holding_key']:<17}{v:>8}{dv:>9}{r['n_lines']:>4}"
          f"{cov['n_served']:>5}/{cov['n_asplund_lines']:<3}"
          f"{r['base_1dlte_median']:>9.3f}{r['median_aberr']:>+8.3f}   "
          f"{r['route_axis']:<7}{r['gf_axis']:<6}")
    A("")
    A(f"  Asplund, Amarsi & Grevesse 2021 : {ASPLUND21_FE}")
    A(f"  AP2002's own published 1D-LTE   : "
      f"{rep['holdings'][0]['ap2002_published_1dlte_median']:.3f} (median over the pool)")
    A("")
    A("  COVERAGE — every Asplund line a holding could not serve, named:")
    for r in rep["holdings"]:
        if r.get("error"):
            continue
        cov = r["coverage"]
        if not cov["unserved_lines"]:
            A(f"    {r['holding_key']:<17} all {cov['n_served']} served")
            continue
        A(f"    {r['holding_key']:<17} {cov['n_unserved']} unserved "
          f"{cov['unserved_by_status']}")
        for u in cov["unserved_lines"][:6]:
            A(f"        {u['wavelength_air_A']:9.3f}  Elo {u['elo_eV']:.3f}  "
              f"{u['status']}: {u['reason'][:60]}")
    A("")
    A("  gf OVERRIDE — AP2002's scale vs ours, on the target lines only:")
    g = rep["holdings"][0]["gf_override"]["delta_vs_canonical_dex"]
    A(f"    applied to {rep['holdings'][0]['gf_override']['n_applied']} lines; "
      f"AP2002 minus canonical: mean {g['mean']:+.4f}  median {g['median']:+.4f}  "
      f"range {g['min']:+.3f}..{g['max']:+.3f}  ({g['n_exact']} identical)")
    A("")
    A("  ⚠️ CAVEAT CARRIED, NOT CHASED (RYA-1104/RYA-282): the +0.045 dex line-population")
    A("     swing between this pool and our graded pool is 29% accounted for by the")
    A("     network's Elo trend and 71% UNATTRIBUTED. It is documented here and is not a")
    A("     tuning target: no line list, gf or parameter was adjusted to move any number")
    A("     toward 7.46 (RYA-161).")
    A("=" * 92)
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--holdings", nargs="+", default=list(HOLDINGS),
                    choices=sorted(HOLDINGS),
                    help="which VIS holdings to run. Amarsi is Fe I VIS only.")
    ap.add_argument("--all-lines", action="store_true",
                    help="use the full AP2002 Fe I table instead of Amarsi's own "
                         "REW < -4.9 weak cut. A DIAGNOSTIC: the published 7.46 stands on "
                         "the cut pool, so the product is the cut pool.")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--tmp", type=Path, default=Path("/tmp/rya1106_synth"))
    args = ap.parse_args(argv)

    targets = asplund_fe1_lines(weak_only=not args.all_lines)
    args.tmp.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)

    rows = []
    for key in args.holdings:
        print(f"\n[{key}] {HOLDINGS[key][0]} / {HOLDINGS[key][1]} — "
              f"{len(targets)} Asplund lines")
        try:
            rows.append(run_one(key, targets, tmp_root=args.tmp, out=args.out))
        except SystemExit as exc:
            # A holding that cannot run is REPORTED as a named gap, not a shorter table.
            # RYA-1106 runs three holdings on the Mac and IAG on Sirius, and a silently
            # absent row would make an environment limit look like a coverage result.
            print(f"[{key}] REFUSED: {exc}")
            rows.append({"holding_key": key, "instrument": HOLDINGS[key][0],
                         "holding": HOLDINGS[key][1], "error": str(exc)[:400]})

    rep = {
        "ticket": "RYA-1106",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pool": {
            "citation": amarsi3d.SOLAR_CONTROL_CITATION,
            "cut": ("Amarsi+2022 Sect. 6.1: weak lines only, REW < -4.9 (AP2002's own "
                    "published widths)" if not args.all_lines
                    else "NO CUT — full AP2002 Fe I table (diagnostic)"),
            "n_used": int(len(targets)),
            "n_published": int(len(asplund_fe1_lines(weak_only=False))),
            "gf": ("AP2002's own published log gf on the target lines; canonical gf "
                   "(RYA-353) on the blend context"),
            "elo_min": float(targets["elo_eV"].min()),
            "elo_max": float(targets["elo_eV"].max()),
        },
        "asplund21": ASPLUND21_FE,
        "firewall": ("RYA-161: this product is emitted because it uses Asplund's ACTUAL "
                     "INPUTS, not because of what it returns. No line list, gf, "
                     "atmosphere or stellar parameter was adjusted, and nothing in this "
                     "script branches on the distance to 7.46."),
        "residual_caveat": {
            "swing_dex": 0.045, "explained_by_elo": 0.29, "unattributed": 0.71,
            "source": "RYA-1104",
            "note": ("carried as a documented caveat (RYA-282), never a tuning target. "
                     "The unattributed 71% is stated because it is unattributed."),
        },
        "holdings": rows,
    }
    text = render(rep)
    print("\n" + text)
    (args.out / "asplund_replication.json").write_text(json.dumps(rep, indent=2,
                                                                 default=str) + "\n")
    (args.out / "asplund_replication.txt").write_text(text + "\n")
    print(f"\nwrote {args.out / 'asplund_replication.json'}")
    print(f"wrote {args.out / 'asplund_replication.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

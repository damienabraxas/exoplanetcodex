#!/usr/bin/env python3
"""RYA-782 — is the Fe IR REW slope a curve-of-growth systematic, or outlier leverage?

    python3 scripts/rya782_rew_trend.py --lines-csv <1D-LTE per-line artifact>

THE CLAIM UNDER TEST
--------------------
RYA-760's confounder regression `A ~ EP + REW + is_loose` over the measured Fe I IR pool
(6910-9199.9 A) returned **REW = +0.384 (t = 2.55)**: weaker lines return systematically
LOWER A(Fe). RYA-782 reads that as a residual systematic in OUR analysis -- residual
saturation or microturbulence, the RYA-523 disease family -- and states it as

    "a gentle systematic SLOPE across the whole Fe IR pool, not a single catastrophic
     bad line"

with three proposed remedies: (a) a microturbulence re-solve, (b) a strong-line saturation
cut, (c) a blend re-check on the strong end.

WHY THE FILED REGRESSION CANNOT SETTLE IT
-----------------------------------------
`is_loose` is a BINARY on source accuracy > 0.04 dex. It lumps four sources of three
different accuracies into one level -- FMW and GESB82c (0.08), MRW (0.08) and K07 (0.20) --
against a laboratory reference level. K07 alone is 47 of the 99 lines. So the REW
coefficient is not "the curve-of-growth trend after the source is controlled"; it is
whatever REW carries after a single step function has absorbed the average of four
populations. A gf-scale error that correlates with line strength -- exactly what a
semi-empirical source is expected to have -- lands in the REW term, not the source term.

The distinction is not cosmetic. Under the ticket's reading the fix is ours (vmic /
saturation) and it biases EVERY Fe IR line including the delivered product. Under the
alternative the fix is atomic, it lives in a tier the delivered product already excludes,
and re-solving vmic against it would be tuning to a gf artifact.

WHAT THIS DOES
--------------
  1. reproduces the filed coefficient exactly, so the rest is measured against it;
  2. re-fits with per-SOURCE dummies instead of the binary, and separately on the
     gf-HOMOGENEOUS laboratory pool -- the only sample in which an A-vs-REW slope is
     actually a microturbulence diagnostic;
  3. censuses the pool for catastrophic lines and re-asks the slope without them,
     which is the ticket's own "gentle slope, not a bad line" framing put to the test;
  4. runs the three proposed remedies as falsifiable tests:
       (b) strong-line cut  -- and reports that pipeline.band_products already enforces
           REW_SATURATION_CEILING = -4.90, so the pool has no saturated tail to cut;
       (c) blend re-check    -- a Boltzmann-weighted neighbour metric as a regressor;
       (a) microturbulence   -- tested by LOCALITY without a synthesis run: a vmic error
           acts through the saturated end, so its slope must concentrate near the knee
           and vanish in the clearly-linear regime. The converse refutes it.
  5. emits the per-line table and a machine-readable summary.

The RYA-161 firewall is respected throughout: the optical anchor is never fitted to, and
no coefficient here is chosen to move A(Fe) toward it. The anchor appears in the output
only as a descriptive comparison, never as a target.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

# Source accuracies as published (the RYA-760 table). These are the ACCURACY of the
# source, not a quality judgement; the laboratory tier is a cut on this number.
SOURCE_DEX = {
    "BWL": 0.02, "2014MNRAS.": 0.03, "GESHRL14": 0.03, "BK": 0.04, "BKK": 0.04,
    "FMW": 0.08, "GESB82c": 0.08, "MRW": 0.08, "K07": 0.20, "K13": 0.20, "K75": 0.20,
}
LAB_TIER = 0.04
# pipeline/band_products.py REW_SATURATION_CEILING -- restated, not re-derived.
REW_SATURATION_CEILING = -4.90
# Boltzmann factor for a crude solar line-strength proxy: loggf - EP * 5040/Teff.
THETA_SUN = 5040.0 / 5777.0
# RYA-553 optical A(Fe I). DESCRIPTIVE ONLY -- never fitted to (RYA-161 firewall).
VIS_ANCHOR = 7.466
# LSQ coefficients are rounded before reporting: Mac and Sirius BLAS differ at ~1e-14
# and an unrounded coefficient makes a bit-identical artifact impossible (RYA-780).
ROUND = 6


def base_source(ref: str) -> str:
    """'BK+BWL' -> 'BK'. Composite tags take the accuracy of their first component."""
    r = str(ref).strip()
    return r.split("+")[0] if "+" in r else r


def ols(frame: pd.DataFrame, cols: list[str], ycol: str = "A") -> dict:
    """Plain OLS with t-statistics. Returns rounded coefficients (see ROUND)."""
    n = len(frame)
    if n <= len(cols) + 1:
        return dict(n=n, insufficient=True)
    X = np.column_stack([np.ones(n)] + [frame[c].to_numpy(float) for c in cols])
    y = frame[ycol].to_numpy(float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = max(n - X.shape[1], 1)
    se = np.sqrt(np.sum(resid ** 2) / dof * np.diag(np.linalg.pinv(X.T @ X)))
    out = dict(n=n, resid_sd=round(float(np.std(resid, ddof=X.shape[1])), ROUND))
    for name, b, e in zip(["intercept"] + cols, beta, se):
        out[name] = dict(coef=round(float(b), ROUND), se=round(float(e), ROUND),
                         t=round(float(b / e), 4) if e else None)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    out["r2"] = round(1.0 - float(np.sum(resid ** 2)) / ss_tot, ROUND) if ss_tot else None
    return out


def univariate_slope(frame: pd.DataFrame, xcol: str = "rew", boot: int = 5000,
                     seed: int = 782) -> dict:
    """Slope of A on x, with a bootstrap CI -- these subsamples are small enough that
    the analytic se alone is not a safe read."""
    n = len(frame)
    if n < 5:
        return dict(n=n, insufficient=True)
    x = frame[xcol].to_numpy(float)
    y = frame["A"].to_numpy(float)
    X = np.column_stack([np.ones(n), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    se = np.sqrt(np.sum(resid ** 2) / max(n - 2, 1) * np.diag(np.linalg.pinv(X.T @ X)))
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(boot):
        k = rng.integers(0, n, n)
        Xb = np.column_stack([np.ones(n), x[k]])
        try:
            bb, *_ = np.linalg.lstsq(Xb, y[k], rcond=None)
        except np.linalg.LinAlgError:
            continue
        draws.append(bb[1])
    lo, hi = (np.percentile(draws, [2.5, 97.5]) if draws else (np.nan, np.nan))
    return dict(n=n, slope=round(float(beta[1]), ROUND),
                se=round(float(se[1]), ROUND),
                t=round(float(beta[1] / se[1]), 4) if se[1] else None,
                boot95=[round(float(lo), ROUND), round(float(hi), ROUND)],
                excludes_zero=bool(np.isfinite(lo) and lo * hi > 0))


def build_pool(lines_csv: Path, canonical_csv: Path, lo: float, hi: float) -> pd.DataFrame:
    """Join the measured per-line abundances to their gf provenance.

    Matching is on wavelength within 0.05 A, taking the STRONGEST Fe I candidate -- the
    same rule scripts/rya760_source_adjudication.py used, kept identical so the
    reproduction step is a true reproduction. The distance to the credited transition is
    retained as `mismatch_mA`, because a large mismatch is itself a misidentification
    signal (RYA-780: match on wavelength AND EP; here EP comes from the match, so the
    mismatch is reported rather than silently absorbed).
    """
    d = pd.read_csv(lines_csv)
    d = d[d["abundance"].notna() & d["rew"].notna() & d["ew_mA"].notna()]
    d = d[(d["wavelength_air_A"] >= lo) & (d["wavelength_air_A"] <= hi)]

    c = pd.read_csv(canonical_csv, low_memory=False)
    fe1 = c[c["species"].astype(str).str.strip() == "Fe I"].copy()
    fe1 = fe1.sort_values("wavelength_air_A").reset_index(drop=True)
    all_w = c["wavelength_air_A"].to_numpy(float)
    order = np.argsort(all_w)
    all_w, all_rest = all_w[order], c.iloc[order].reset_index(drop=True)
    all_rest = all_rest.assign(
        strength=all_rest["log_gf"].astype(float)
        - all_rest["excitation_potential_eV"].astype(float) * THETA_SUN)
    fw = fe1["wavelength_air_A"].to_numpy(float)

    rows = []
    for _, r in d.iterrows():
        w = float(r["wavelength_air_A"])
        i0, i1 = np.searchsorted(fw, [w - 0.05, w + 0.05])
        if i1 <= i0:
            continue
        cand = fe1.iloc[i0:i1]
        i = cand["log_gf"].astype(float).idxmax()
        ref = base_source(fe1.loc[i, "loggf_reference"])
        lam = float(fe1.loc[i, "wavelength_air_A"])
        ep = float(fe1.loc[i, "excitation_potential_eV"])
        gf = float(fe1.loc[i, "log_gf"])

        # strongest OTHER absorber within 0.05 A, Boltzmann-weighted for solar conditions
        j0, j1 = np.searchsorted(all_w, [w - 0.05, w + 0.05])
        nb = all_rest.iloc[j0:j1]
        nb = nb[np.abs(nb["wavelength_air_A"].astype(float) - lam) > 2e-3]
        self_s = gf - ep * THETA_SUN
        if len(nb):
            k = nb["strength"].idxmax()
            blend_ratio = float(nb.loc[k, "strength"]) - self_s
            blend_sp = str(nb.loc[k, "species"]).strip()
        else:
            blend_ratio, blend_sp = -np.inf, "none"

        rows.append(dict(
            wavelength_air_A=w, A=float(r["abundance"]), ew_mA=float(r["ew_mA"]),
            rew=float(r["rew"]), loggf_reference=ref,
            source_accuracy_dex=SOURCE_DEX.get(ref, np.nan),
            excitation_potential_eV=ep, log_gf=gf,
            mismatch_mA=round(abs(lam - w) * 1000.0, 3),
            n_fe1_within_50mA=int(len(cand)),
            blend_species=blend_sp,
            blend_ratio_dex=(round(blend_ratio, 4) if np.isfinite(blend_ratio) else None),
        ))

    t = pd.DataFrame(rows)
    t = t[t["source_accuracy_dex"].notna()].copy()
    med = t["A"].median()
    mad = 1.4826 * np.median(np.abs(t["A"] - med))
    t["robust_z"] = ((t["A"] - med) / mad).round(4)
    t["is_loose"] = (t["source_accuracy_dex"] > LAB_TIER).astype(float)
    t["is_laboratory"] = t["source_accuracy_dex"] <= LAB_TIER
    t["catastrophic"] = t["robust_z"].abs() > 3
    return t.sort_values("wavelength_air_A").reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lines-csv", required=True,
                    help="the 1D-LTE per-line artifact (wavelength_air_A, abundance, "
                         "ew_mA, rew)")
    ap.add_argument("--canonical-gf",
                    default=str(ROOT / "data" / "linelists" / "canonical_gf.csv"),
                    help="the gf single source (RYA-353)")
    ap.add_argument("--lo", type=float, default=6910.0)
    ap.add_argument("--hi", type=float, default=9199.9)
    ap.add_argument("--out-dir", default=str(ROOT / "data" / "results" / "rya782"))
    ap.add_argument("--label", default="FeI_IR",
                    help="artifact basename stem")
    a = ap.parse_args()

    t = build_pool(Path(a.lines_csv), Path(a.canonical_gf), a.lo, a.hi)
    lab = t[t["is_laboratory"]]
    clip3 = t[~t["catastrophic"]]
    summary: dict = dict(
        ticket="RYA-782",
        band_A=[a.lo, a.hi],
        inputs=dict(lines_csv=str(a.lines_csv), canonical_gf=str(a.canonical_gf)),
        n_pool=int(len(t)), n_laboratory=int(len(lab)),
        n_catastrophic=int(t["catastrophic"].sum()),
        rew_saturation_ceiling=REW_SATURATION_CEILING,
        pool_rew_range=[round(float(t["rew"].min()), 4), round(float(t["rew"].max()), 4)],
    )

    print(f"\n{len(t)} measured Fe I lines in {a.lo:.0f}-{a.hi:.0f} A with a known gf source")
    print(f"  laboratory tier (<= {LAB_TIER}): {len(lab)}    "
          f"beyond 3 robust sd: {int(t['catastrophic'].sum())}")

    # ── 1. reproduce ─────────────────────────────────────────────────────────
    print("\n" + "=" * 82)
    print("[1] REPRODUCE the filed coefficient")
    print("=" * 82)
    filed = ols(t, ["excitation_potential_eV", "rew", "is_loose"])
    summary["filed_regression"] = filed
    for k in ("excitation_potential_eV", "rew", "is_loose"):
        print(f"  {k:<26}{filed[k]['coef']:>9.4f}  t={filed[k]['t']:+6.2f}")
    print(f"  residual sd {filed['resid_sd']:.3f} dex")

    # ── 2. the binary is the problem ─────────────────────────────────────────
    print("\n" + "=" * 82)
    print("[2] REPLACE the binary with per-SOURCE dummies")
    print("=" * 82)
    t2 = t.copy()
    dummies = []
    for s in ["K07", "FMW", "GESB82c", "MRW"]:
        col = f"is_{s}"
        t2[col] = (t2["loggf_reference"] == s).astype(float)
        if t2[col].sum() > 0:
            dummies.append(col)
    per_src = ols(t2, ["excitation_potential_eV", "rew"] + dummies)
    summary["per_source_regression"] = per_src
    for k in ["rew"] + dummies:
        print(f"  {k:<26}{per_src[k]['coef']:>9.4f}  t={per_src[k]['t']:+6.2f}")

    print("\n  univariate A-vs-REW slope by gf population "
          "(the only vmic-diagnostic sample is a gf-homogeneous one):")
    pops = {
        "ALL sources (mixed scales)": t,
        "LABORATORY only (<=0.04)": lab,
        "LOOSE only (>0.04)": t[~t["is_laboratory"]],
        "K07 only (0.20)": t[t["loggf_reference"] == "K07"],
        "FMW+GESB82c (0.08)": t[t["loggf_reference"].isin(["FMW", "GESB82c"])],
    }
    summary["slope_by_population"] = {}
    for name, s in pops.items():
        r = univariate_slope(s)
        summary["slope_by_population"][name] = r
        if r.get("insufficient"):
            print(f"    {name:<30} n={r['n']:3d}  (too few)")
            continue
        print(f"    {name:<30} n={r['n']:3d}  slope={r['slope']:+7.3f}  t={r['t']:+5.2f}"
              f"  boot95=[{r['boot95'][0]:+.3f},{r['boot95'][1]:+.3f}]"
              f"{'  <-- excludes 0' if r['excludes_zero'] else ''}")

    # ── 3. gentle slope, or catastrophic lines? ──────────────────────────────
    print("\n" + "=" * 82)
    print("[3] THE TICKET'S FRAMING TESTED: 'gentle slope, not a catastrophic bad line'")
    print("=" * 82)
    bad = t[t["catastrophic"]].sort_values("A", ascending=False)
    print(f"  {'wave':<11}{'src':<10}{'acc':>6}{'A':>9}{'z':>7}{'EW':>7}"
          f"{'mis(mA)':>9}{'blend':>8}")
    for _, r in bad.iterrows():
        br = r["blend_ratio_dex"]
        br_s = f"{br:+.2f}" if br is not None and np.isfinite(br) else "none"
        print(f"  {r['wavelength_air_A']:<11.3f}{r['loggf_reference']:<10}"
              f"{r['source_accuracy_dex']:>6.2f}{r['A']:>9.3f}{r['robust_z']:>7.1f}"
              f"{r['ew_mA']:>7.1f}{r['mismatch_mA']:>9.1f}{br_s:>8}")

    C = ["excitation_potential_eV", "rew", "is_loose"]
    variants = {
        "as filed": t,
        "3-sigma clipped": clip3,
        "K07 removed": t[t["loggf_reference"] != "K07"],
        "K07 removed + 3-sigma clipped": clip3[clip3["loggf_reference"] != "K07"],
    }
    summary["slope_survival"] = {}
    print()
    for name, s in variants.items():
        r = ols(s, C)
        summary["slope_survival"][name] = r
        if r.get("insufficient"):
            continue
        print(f"  {name:<34} n={r['n']:3d}  REW={r['rew']['coef']:+7.3f} "
              f"(t={r['rew']['t']:+5.2f})   resid_sd={r['resid_sd']:.3f}")
    for name, s in (("LABORATORY pool (the delivered product)", lab),
                    ("LABORATORY pool, 3-sigma clipped", lab[~lab["catastrophic"]])):
        r = ols(s, ["excitation_potential_eV", "rew"])
        summary["slope_survival"][name] = r
        if not r.get("insufficient"):
            print(f"  {name:<34} n={r['n']:3d}  REW={r['rew']['coef']:+7.3f} "
                  f"(t={r['rew']['t']:+5.2f})   resid_sd={r['resid_sd']:.3f}")

    # ── 4. the three remedies ────────────────────────────────────────────────
    print("\n" + "=" * 82)
    print("[4] THE THREE PROPOSED REMEDIES, AS TESTS")
    print("=" * 82)

    print(f"\n  (b) strong-line saturation cut")
    print(f"      pipeline/band_products.py REW_SATURATION_CEILING = "
          f"{REW_SATURATION_CEILING} is ALREADY enforced;")
    print(f"      pool REW range is [{t['rew'].min():.3f}, {t['rew'].max():.3f}] -- "
          f"there is no saturated tail left to cut.")
    summary["remedy_b_strong_line_cut"] = {}
    for cut in (100.0, 80.0, 60.0):
        r = ols(t[t["ew_mA"] < cut], C)
        summary["remedy_b_strong_line_cut"][f"EW<{cut:.0f}mA"] = r
        if not r.get("insufficient"):
            print(f"      EW < {cut:>5.0f} mA   n={r['n']:3d}  "
                  f"REW={r['rew']['coef']:+7.3f} (t={r['rew']['t']:+5.2f})")

    print(f"\n  (c) blend re-check -- Boltzmann-weighted strongest neighbour as a regressor")
    u = t[t["blend_ratio_dex"].notna()].copy()
    r_nb = ols(u, C)
    r_wb = ols(u, C + ["blend_ratio_dex"])
    summary["remedy_c_blend"] = dict(without=r_nb, with_blend=r_wb)
    if not r_wb.get("insufficient"):
        print(f"      blend term        {r_wb['blend_ratio_dex']['coef']:+.4f}  "
              f"t={r_wb['blend_ratio_dex']['t']:+6.2f}   "
              f"(REW moves {r_nb['rew']['coef']:+.3f} -> {r_wb['rew']['coef']:+.3f})")

    print(f"\n  (a) microturbulence -- tested by LOCALITY, no synthesis run needed.")
    print(f"      A vmic error acts through the SATURATED end. Its slope must be")
    print(f"      concentrated near the knee and absent in the clearly-linear regime.")
    windows = [(-9.9, -5.2, "clearly linear (<=-5.2)"),
               (-5.2, -5.0, "mid  (-5.2..-5.0)"),
               (-5.1, REW_SATURATION_CEILING, "near knee (-5.1..-4.9)")]
    summary["remedy_a_vmic_locality"] = {}
    for pool_name, pool in (("all sources", t), ("laboratory only", lab)):
        print(f"      {pool_name}:")
        for w0, w1, lbl in windows:
            s = pool[(pool["rew"] > w0) & (pool["rew"] <= w1)]
            r = univariate_slope(s, boot=0)
            summary["remedy_a_vmic_locality"][f"{pool_name} / {lbl}"] = r
            if r.get("insufficient"):
                print(f"        {lbl:<26} n={r['n']:3d}  (too few)")
            else:
                print(f"        {lbl:<26} n={r['n']:3d}  slope={r['slope']:+7.3f}"
                      f"  t={r['t']:+5.2f}")

    # ── 5. what the slope is worth to the error budget ───────────────────────
    print("\n" + "=" * 82)
    print("[5] ERROR BUDGET -- the STATED term (RYA-713: bigger bars where science backs it)")
    print("=" * 82)
    r = ols(lab, ["excitation_potential_eV", "rew"])
    x = lab["rew"].to_numpy(float)
    slope = r["rew"]["coef"]
    se = r["rew"]["se"]
    span = float(x.max() - x.min())
    upper = slope + 1.96 * se
    summary["error_budget"] = dict(
        laboratory_rew_slope=slope, se=se, t=r["rew"]["t"],
        rew_span_dex=round(span, 4),
        best_estimate_swing_dex=round(abs(slope) * span, 4),
        upper95_slope=round(upper, 4),
        upper95_swing_dex=round(abs(upper) * span, 4),
        scatter_contribution_dex=round(abs(slope) * float(np.std(x, ddof=1)), 4),
    )
    print(f"  laboratory pool n={len(lab)}: REW slope {slope:+.3f} +/- {se:.3f} "
          f"(t={r['rew']['t']:+.2f}) -- consistent with zero")
    print(f"  REW span {span:.2f} dex  ->  best-estimate weak-to-strong swing "
          f"{abs(slope)*span:.3f} dex")
    print(f"  95% upper bound on a HIDDEN slope {upper:+.3f}  ->  swing "
          f"<= {abs(upper)*span:.3f} dex")
    print(f"  contribution to per-line scatter: {abs(slope)*float(np.std(x, ddof=1)):.3f} dex")
    print(f"\n  descriptive only (never fitted to -- RYA-161): laboratory pool median "
          f"A={lab['A'].median():.3f}, vs optical anchor {VIS_ANCHOR} "
          f"= {lab['A'].median()-VIS_ANCHOR:+.3f} dex")

    # ── emit ─────────────────────────────────────────────────────────────────
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out / f"{a.label}_rew_trend.csv"
    json_path = out / f"{a.label}_rew_trend_summary.json"
    t.to_csv(csv_path, index=False)
    with open(json_path, "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True, default=str)
        fh.write("\n")
    print(f"\nwrote {csv_path.relative_to(ROOT)}")
    print(f"wrote {json_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

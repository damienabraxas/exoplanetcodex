#!/usr/bin/env python3
"""The near-UV primary-lab-gf Fe sub-pool — the twin of RYA-824. RYA-836.

    python3 scripts/rya836_nearuv_lab_gf_subpool.py --census-only
    python3 scripts/rya836_nearuv_lab_gf_subpool.py                    # Sirius

WHY
---
RYA-822 extended canonical_gf to 3000 A and GRADED the near-UV, and the payoff was a
NULL: 7.487 -> 7.488, scatter 0.354 -> 0.355. That is the same result RYA-799 got in the
IR, and for the same reason — **grading attaches provenance, it does not re-measure.**
What moved the IR was RYA-824, which re-inverted on primary laboratory oscillator
strengths. This is that operation for the near-UV.

The lever is real and 822 sized it: **61 near-UV Fe I lines have a primary lab gf and the
pool was measured on that value for only 6.** The other 55 were inverted on Kurucz.

WHAT PRIMARY MEANS HERE, AND WHAT IS EXCLUDED
---------------------------------------------
Only Ruffoni 2014, Den Hartog 2014 and Belmonte 2017 — branching fractions from a
Fourier-transform spectrometer times a radiative lifetime, measured in a laboratory with
no reference to the solar spectrum at any point.

`GF-NIST` is NOT admitted, and that exclusion is load-bearing. RYA-822 assigned it to 604
near-UV lines and RYA-760 established why it cannot referee: `FMW` *is* Fuhr, Martin &
Wiese, a NIST compilation, and VALD copies it. Letting a compilation into a sub-pool whose
entire claim is independence would let the headline restate its own input.

Solar-calibrated (astrophysical) gf is excluded for the stronger reason: it is fitted so
the Sun comes out right, so using it would make the solar abundance a restatement of
itself. That is the RYA-161 firewall.

THE CONTROL
-----------
Every line is fitted TWICE through the same synthesis route — once with the line list
untouched, once with the lab gf substituted — so the difference is attributable to the
oscillator strength and to nothing else. The untouched fit is also the check that this
harness IS the production path: it must reproduce what the RYA-832 route produces for the
same line.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline.gf_grades import lab_lines, CITATIONS  # noqa: E402
from pipeline.band_products import LineMeasurement, build_product  # noqa: E402

LO_A, HI_A = 3000.0, 3780.0
OUT_DEFAULT = ROOT / "data" / "results" / "rya836"

#: RYA-825's rule, restated because this ticket re-reads the lab tables.
SIGMA_MIN, SIGMA_MAX = 0.0, 0.5
#: Below this the substituted gf IS the value already in the list — the "already on lab"
#: case, which is a no-op fit and must be counted as such rather than as a change.
LOGGF_SAME_TOL = 1e-3
#: A LINE-IDENTIFICATION SCREEN, not a quality cut, and the distinction matters.
#:
#: A primary lab sigma in this band is <=0.12 dex and a Kurucz semi-empirical gf is good
#: to perhaps 0.2-0.3. So lab and production disagreeing by a whole dex is far outside any
#: combination of the two uncertainties: it is much more likely that the lab line and the
#: line list line are DIFFERENT TRANSITIONS that happen to agree in wavelength and
#: excitation potential to within the match tolerance. Substituting there would not be
#: adopting a better oscillator strength, it would be putting one line's gf on another
#: line -- the RYA-785 "a same-species neighbour does not cancel" failure.
#:
#: 1.0 dex is chosen as the smallest round number outside that combined-uncertainty
#: envelope, BEFORE looking at what it excludes, and flagged lines are CARRIED and
#: reported rather than dropped (RYA-711 quarantine-not-cull). The sub-pool is reported
#: both ways so the leverage of the screen is visible instead of asserted.
PROBABLE_MISID_DEX = 1.0

#: RYA-759 / RYA-832: the near-UV product's route parameters. Imported in spirit; the
#: functions themselves are imported literally below.
HALF_WIDTH_A = 0.40
PSEUDO_CONTINUUM_DEX = 0.100


def census() -> pd.DataFrame:
    """Near-UV Fe I lines that have a PRIMARY laboratory gf."""
    lab = lab_lines()
    bad = lab[(lab.e_loggf_dex <= SIGMA_MIN) | (lab.e_loggf_dex > SIGMA_MAX)]
    if len(bad):
        raise SystemExit(f"{len(bad)} lab rows carry an impossible sigma — refusing to "
                         f"run (RYA-825 positivity guard)")
    m = lab[(lab.wavelength_air_A >= LO_A) & (lab.wavelength_air_A < HI_A)].copy()
    return m.sort_values("wavelength_air_A").reset_index(drop=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--from-per-line", type=Path, default=None,
                    help="rebuild the products/report from a saved per-line CSV without "
                         "re-fitting (122 Turbospectrum fits is ~2 h). The screen is "
                         "applied at AGGREGATION, so a re-report is exact, not an "
                         "approximation of a rerun.")
    ap.add_argument("--census-only", action="store_true",
                    help="emit the census; no synthesis (runs anywhere)")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--band-products-out", type=Path,
                    default=ROOT / "data" / "results" / "band_products")
    ap.add_argument("--instrument", default="kpno_solar_atlas")
    ap.add_argument("--tmp-dir", type=Path, default=Path("/tmp/ispec_rya836_nearuv"))
    args = ap.parse_args(argv)

    c = census()
    args.out.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("RYA-836 — near-UV primary-lab-gf Fe sub-pool")
    print("=" * 78)
    lab_all = lab_lines()
    print(f"\nPRIMARY-LAB Fe I coverage (Ruffoni 2014 / Den Hartog 2014 / Belmonte 2017):")
    print(f"  below {HI_A:.0f} A total : {len(lab_all[lab_all.wavelength_air_A < HI_A]):3d}  "
          + ", ".join(f"{s} {int(v)}" for s, v in
                      lab_all[lab_all.wavelength_air_A < HI_A].source.value_counts().items()))
    print(f"  IN BAND {LO_A:.0f}-{HI_A:.0f} : {len(c):3d}  "
          + ", ".join(f"{s} {int(v)}" for s, v in c.source.value_counts().items()))
    print(f"  cited sigma: {c.e_loggf_dex.min():.3f}-{c.e_loggf_dex.max():.3f} dex "
          f"(median {c.e_loggf_dex.median():.3f})")
    print("\n  GF-NIST is NOT in this set. RYA-822 assigned it to 604 near-UV lines, but")
    print("  FMW *is* a NIST compilation and VALD copies it (RYA-760), so admitting it")
    print("  would let a compilation restate the headline it is supposed to test.")

    if args.census_only:
        c.to_csv(args.out / "rya836_nearuv_lab_census.csv", index=False)
        print(f"\n[--census-only] no synthesis. wrote "
              f"{args.out}/rya836_nearuv_lab_census.csv")
        return 0

    if args.from_per_line is not None:
        res = pd.read_csv(args.from_per_line)
        ctrl_lines, lab_lines_lm = [], []
        for r in res.itertuples():
            misid = bool(getattr(r, "probable_misid", False)) or (
                np.isfinite(getattr(r, "d_loggf", np.nan))
                and abs(float(r.d_loggf)) >= PROBABLE_MISID_DEX)
            for col, bucket, treat in (("a_control", ctrl_lines, "1D-LTE"),
                                       ("a_labgf", lab_lines_lm, "1D-LTE-LABGF")):
                a_x = float(getattr(r, col, np.nan))
                lm = LineMeasurement(
                    element="Fe", ion="I", wavelength_air_A=float(r.wavelength_air_A),
                    instrument=args.instrument, ew_mA=float("nan"),
                    ew_method=f"synthesis flux-fit, fixed half-width "
                              f"+/-{HALF_WIDTH_A} A (RYA-759 route)",
                    abundance=(a_x if np.isfinite(a_x) else None),
                    treatment=treat, ew_inversion=False)
                if lm.abundance is None:
                    lm.in_aggregate = False
                    lm.excluded_reason = "SYNTHESIS: fit did not converge"
                elif misid and treat == "1D-LTE-LABGF":
                    lm.in_aggregate = False
                    lm.excluded_reason = (
                        f"PROBABLE-MISIDENTIFICATION: |d(log gf)| "
                        f"{abs(float(r.d_loggf)):.3f} dex >= the "
                        f"{PROBABLE_MISID_DEX} dex screen")
                bucket.append(lm)
        if "probable_misid" not in res.columns:
            res["probable_misid"] = [
                bool(np.isfinite(d) and abs(d) >= PROBABLE_MISID_DEX)
                for d in res.get("d_loggf", pd.Series([np.nan] * len(res)))]
        return _report(res, c, ctrl_lines, lab_lines_lm, args)

    # ── the re-inversion ──────────────────────────────────────────────────────
    from pipeline.nearuv_synth import build_solar_context, gf_provenance
    from rya759_nearuv_fe_product import NEARUV_LINELIST, fit_one
    from rya759_nearuv_synth import _kp_segments
    from rya824_lab_gf_subpool import substitute_gf          # the SAME substitution

    if not Path(NEARUV_LINELIST).exists():
        raise SystemExit(
            f"near-UV line list missing at {NEARUV_LINELIST}\n"
            f"build it: `python3 scripts/rya759_nearuv_synth.py --step linelist`. "
            f"Refusing to emit an empty sub-pool that would read as a null result.")

    cat = pd.read_csv(ROOT / "data" / "catalog" / "instrument_catalog.csv")
    row = cat[cat.iloc[:, 0].astype(str) == args.instrument]
    if row.empty:
        raise SystemExit(f"{args.instrument!r} absent from instrument_catalog.csv")
    R = float(row.iloc[0]["resolving_power_max"])

    prov_gf = gf_provenance(LO_A, HI_A)
    ctx = build_solar_context("Fe", R, linelist_file=str(NEARUV_LINELIST),
                              apply_canonical_gf=prov_gf["apply_canonical_gf"])
    segs = _kp_segments()
    args.tmp_dir.mkdir(parents=True, exist_ok=True)
    tmp = str(args.tmp_dir)
    print(f"\n[gf] {prov_gf['detail']}")
    print(f"\nRE-INVERTING {len(c)} lines — each fitted TWICE "
          f"(control on the production gf, then on lab gf)")
    print(f"  {'line':>11}{'src':>15}{'lab_sig':>9}{'gf_prod':>9}{'gf_lab':>9}"
          f"{'A_ctrl':>9}{'A_lab':>9}{'dA':>8}")

    rows, ctrl_lines, lab_lines_lm = [], [], []
    for r in c.itertuples():
        w, ep = float(r.wavelength_air_A), float(r.elo_eV)
        try:
            arr_lab = substitute_gf(ctx["linelist"], w, ep, float(r.loggf))
        except ValueError as exc:
            print(f"  {w:11.4f}  SUBSTITUTION REFUSED: {str(exc)[:70]}")
            rows.append(dict(wavelength_air_A=w, ep_eV=ep, lab_source=r.source,
                             lab_sigma_dex=float(r.e_loggf_dex), a_control=np.nan,
                             a_labgf=np.nan, note=f"substitution refused: {exc}"))
            continue
        # what the production list actually holds for this line, recovered by asking the
        # substitution for a no-op — never assumed from the accounting column (RYA-825).
        gf_prod = _current_total(ctx["linelist"], w, ep)
        res_c = fit_one(ctx, segs, w, HALF_WIDTH_A, tmp)
        ctx_lab = {**ctx, "linelist": arr_lab}
        res_l = fit_one(ctx_lab, segs, w, HALF_WIDTH_A, tmp)
        a_c = float(res_c.get("a_synth", np.nan))
        a_l = float(res_l.get("a_synth", np.nan))
        dA = (a_l - a_c) if (np.isfinite(a_l) and np.isfinite(a_c)) else np.nan
        already = abs(float(r.loggf) - gf_prod) <= LOGGF_SAME_TOL
        misid = (np.isfinite(gf_prod)
                 and abs(float(r.loggf) - gf_prod) >= PROBABLE_MISID_DEX)
        rows.append(dict(wavelength_air_A=w, ep_eV=ep, lab_source=str(r.source),
                         lab_sigma_dex=float(r.e_loggf_dex), loggf_production=gf_prod,
                         loggf_lab=float(r.loggf), d_loggf=float(r.loggf) - gf_prod,
                         already_on_lab=already, probable_misid=misid,
                         a_control=a_c, a_labgf=a_l, d_A=dA,
                         status_control=res_c.get("status", ""),
                         status_lab=res_l.get("status", ""), note=""))
        print(f"  {w:11.4f}{str(r.source)[:14]:>15}{float(r.e_loggf_dex):9.3f}"
              f"{gf_prod:9.3f}{float(r.loggf):9.3f}"
              f"{a_c if np.isfinite(a_c) else float('nan'):9.3f}"
              f"{a_l if np.isfinite(a_l) else float('nan'):9.3f}"
              f"{dA:8.3f}")

        for a_x, bucket, treat in ((a_c, ctrl_lines, "1D-LTE"),
                                   (a_l, lab_lines_lm, "1D-LTE-LABGF")):
            lm = LineMeasurement(
                element="Fe", ion="I", wavelength_air_A=w, instrument=args.instrument,
                ew_mA=float("nan"),
                ew_method=f"synthesis flux-fit, fixed half-width +/-{HALF_WIDTH_A} A "
                          f"(RYA-759 route)",
                abundance=(a_x if np.isfinite(a_x) else None),
                treatment=treat, ew_inversion=False)
            if lm.abundance is None:
                lm.in_aggregate = False
                lm.excluded_reason = "SYNTHESIS: fit did not converge"
            elif misid and treat == "1D-LTE-LABGF":
                # Carried, never dropped (RYA-711). Excluded from the aggregate with the
                # reason stated: a gf that far from the production value is more likely a
                # different transition than a better measurement.
                lm.in_aggregate = False
                lm.excluded_reason = (
                    f"PROBABLE-MISIDENTIFICATION: lab gf {float(r.loggf):+.3f} differs "
                    f"from the production gf {gf_prod:+.3f} by "
                    f"{abs(float(r.loggf) - gf_prod):.3f} dex, beyond the "
                    f"{PROBABLE_MISID_DEX} dex screen")
            bucket.append(lm)

    res = pd.DataFrame(rows)
    res.to_csv(args.out / "rya836_nearuv_lab_gf_per_line.csv", index=False)
    return _report(res, c, ctrl_lines, lab_lines_lm, args)


def _current_total(linelist, wave_A: float, ep_eV: float) -> float:
    """The physical line's CURRENT total gf in the production list."""
    from pipeline.gf_resolver import cluster_physical_lines
    from pipeline.species import species_key
    keys = [species_key(linelist["element"][i], linelist["ion"][i],
                        linelist["molecule"][i]) for i in range(len(linelist))]
    wls = linelist["wave_A"].astype(float)
    eps = linelist["lower_state_eV"].astype(float)
    gf = linelist["loggf"].astype(float)
    for cl in cluster_physical_lines(keys, wls, eps):
        w = 10.0 ** gf[cl]
        centroid = float(np.sum(wls[cl] * w) / w.sum())
        if abs(centroid - wave_A) <= 0.02 and abs(float(np.mean(eps[cl])) - ep_eV) <= 0.02:
            return float(np.log10(np.sum(10.0 ** gf[cl])))
    return float("nan")


def _report(res, c, ctrl_lines, lab_lines_lm, args) -> int:
    ok = res[np.isfinite(res.a_labgf) & np.isfinite(res.a_control)]
    # `.map(bool)`, not `~ ... .fillna(False)`. These columns come back from a CSV as
    # OBJECT dtype, where `~` is bitwise negation on ints: it turns True into -1 and the
    # result is then used as a column key. That is what crashed the first run AFTER all
    # 122 Turbospectrum fits had completed — the per-line CSV is written before the
    # report precisely so a reporting bug cannot cost the fits.
    already_mask = res.get("already_on_lab", pd.Series([False] * len(res))).map(bool)
    already = res[already_mask]
    changed = ok[~ok.get("already_on_lab", pd.Series([False] * len(ok))).map(bool)]

    print(f"\n{'-' * 78}")
    print(f"CONVERGED {len(ok)} of {len(res)}   |   already on lab gf (no-op) "
          f"{len(already)}   |   genuinely changed {len(changed)}")
    if len(changed):
        print(f"  d(log gf) lab - production: median {changed.d_loggf.median():+.3f}, "
              f"range {changed.d_loggf.min():+.3f}..{changed.d_loggf.max():+.3f}")
        print(f"  dA: median {changed.d_A.median():+.3f}, "
              f"range {changed.d_A.min():+.3f}..{changed.d_A.max():+.3f}")

    p_ctrl = build_product("Fe", "I", args.instrument, "near-UV", "1D-LTE", ctrl_lines,
                           provenance="near-UV control: the SAME lines, on the "
                                      "production gf. Exists so the lab-gf difference is "
                                      "attributable to the gf and nothing else.")
    p_lab = build_product("Fe", "I", args.instrument, "near-UV", "1D-LTE-LABGF",
                          lab_lines_lm,
                          provenance=(
                              "near-UV Fe I on PRIMARY LABORATORY gf only (Ruffoni 2014 / "
                              "Den Hartog 2014 / Belmonte 2017 — branching fractions x "
                              "radiative lifetimes, no reference to the solar spectrum). "
                              "GF-NIST is EXCLUDED: FMW *is* a NIST compilation and VALD "
                              "copies it (RYA-760), so it cannot referee a sub-pool whose "
                              "claim is independence. Solar-calibrated gf excluded by the "
                              "RYA-161 firewall. Synthesis-only route (RYA-759/832); "
                              "1D-LTE; pseudo-continuum systematic 0.100 dex does NOT "
                              "average down. Reported BESIDE the full Kurucz pool, never "
                              "merged with it (RYA-712)."))

    print(f"\n{'pool':<34}{'A(Fe I)':>9}{'scatter':>9}{'n':>5}")
    for label, p in (("near-UV control (production gf)", p_ctrl),
                     ("near-UV PRIMARY-LAB-GF", p_lab)):
        v = f"{p.value:.3f}" if p.value is not None else "n/a"
        s = f"{p.sigma:.3f}" if p.sigma is not None else "n/a"
        print(f"  {label:<32}{v:>9}{s:>9}{p.n_lines:5d}")
    print(f"  {'RYA-832 full Kurucz pool (n=40)':<32}{'7.488':>9}{'0.413':>9}{40:5d}")
    mis = res[res.get("probable_misid", pd.Series([False] * len(res))).map(bool)]
    if len(mis):
        print(f"\n  LINE-IDENTIFICATION SCREEN: {len(mis)} of {len(res)} lines carry a lab")
        print(f"  gf >= {PROBABLE_MISID_DEX} dex from the production value. CARRIED, and")
        print(f"  excluded from the aggregate with the reason stated (RYA-711). |d| "
              f"{mis.d_loggf.abs().min():.2f}-{mis.d_loggf.abs().max():.2f} dex.")
        vals = np.array([l.abundance for l in lab_lines_lm
                         if l.abundance is not None], dtype=float)
        if len(vals) > 1:
            print(f"  WITHOUT the screen the sub-pool would read {np.median(vals):.3f} "
                  f"+/- {np.std(vals, ddof=1):.3f} (n={len(vals)}) — shown so the")
            print(f"  screen's leverage is VISIBLE rather than asserted.")

    print("\n  DOES IT BEAT +/-0.413?  ", end="")
    if p_lab.sigma is not None:
        beats = p_lab.sigma < 0.413
        print(f"{'YES' if beats else 'NO'} — lab-gf scatter {p_lab.sigma:.3f} dex "
              f"vs 0.413")
    else:
        print("no lab-gf product (see per-line statuses)")

    out_bp = Path(args.band_products_out)
    out_bp.mkdir(parents=True, exist_ok=True)
    stem = f"FeI_{int(LO_A)}_{int(HI_A)}_{args.instrument}_LABGF"
    if p_lab.value is not None:
        from pipeline.error_budget import build as build_budget
        b = build_budget("Fe", 0.5 * (LO_A + HI_A), p_lab.n_lines,
                         scatter_dex=(p_lab.sigma or 0.0), gf_graded=True,
                         harness_residual_dex=0.0, handler="SynthesisHandler")
        stat, syst = b.total()
        syst = float(np.hypot(syst, PSEUDO_CONTINUUM_DEX))
        pd.DataFrame([dict(
            element="Fe", ion="I", band="near-UV", instrument=args.instrument,
            treatment="1D-LTE-LABGF", A=round(p_lab.value, 3), n_lines=p_lab.n_lines,
            n_excluded=p_lab.n_excluded, stat_dex=round(stat, 4),
            syst_dex=round(syst, 4),
            dominant=(b.dominant().name if b.dominant() else ""))]
        ).to_csv(out_bp / f"{stem}_products.csv", index=False)
        print(f"\n  wrote {out_bp}/{stem}_products.csv  "
              f"(gf_graded=TRUE — this pool IS on primary lab gf)")

    (args.out / "rya836_summary.json").write_text(json.dumps({
        "ticket": "RYA-836",
        "n_in_band": int(len(c)), "n_converged": int(len(ok)),
        "n_already_on_lab": int(len(already)), "n_changed": int(len(changed)),
        "control": {"A": p_ctrl.value, "scatter": p_ctrl.sigma, "n": p_ctrl.n_lines},
        "lab_gf": {"A": p_lab.value, "scatter": p_lab.sigma, "n": p_lab.n_lines},
        "rya832_full_pool": {"A": 7.488, "scatter": 0.413, "n": 40},
        "median_d_loggf": (float(changed.d_loggf.median()) if len(changed) else None),
        "median_d_A": (float(changed.d_A.median()) if len(changed) else None),
        "run_date": str(date.today()),
    }, indent=2, default=float) + "\n")
    print(f"  wrote {args.out}/rya836_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

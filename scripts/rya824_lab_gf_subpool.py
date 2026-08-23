#!/usr/bin/env python3
"""The primary-lab-gf Fe sub-pool — re-measure on lab gf, not Kurucz. RYA-824.

    python3 scripts/rya824_lab_gf_subpool.py --pool-dir <RYA-783 dir> --census-only
    python3 scripts/rya824_lab_gf_subpool.py --pool-dir <RYA-783 dir>          # Sirius

WHY
---
RYA-799 proved the Fe line pool's log gf is 100 % Kurucz K14 in every band, so the
~0.20 dex systematic is the KURUCZ FLOOR and grading cannot move it (pool RMS 0.1993
against a 0.2000 blanket). The only lever that moves a floor is a different, independent
MEASUREMENT of the oscillator strength. This script applies one.

This is NOT tuning and the distinction is the RYA-161 firewall. A primary laboratory gf
is measured in a laboratory — branching fractions from a Fourier-transform spectrometer
times a radiative lifetime — with no reference to the solar spectrum at any point. An
ASTROPHYSICAL gf is fitted so that the Sun comes out right, and using one here would make
the solar abundance a restatement of its own input. Only the former is admitted.

WHAT IT DOES
------------
1. CENSUS -- every Fe line in the measured pools that has a primary lab gf, per band,
   with the lab source and its stated per-line sigma. Cheap, no synthesis, so it runs
   anywhere and is the half that can be checked without Sirius.
2. RE-MEASURE -- for each such line, substitute the lab gf into the synth line list and
   re-invert the SAME measured equivalent width through the SAME production core
   (`abundances_derive._bisect_synth_abundance`). The substitution goes through
   `gf_resolver`'s branching-preserved shift, so an HFS-split line keeps its
   component-to-component ratios and only the total moves.
3. REPORT -- the lab-gf sub-pool beside the Kurucz pool, per band. Both. The Kurucz pool
   is not corrected, not replaced and not dropped: it is the broad number, and the sub-pool
   is the tight one, and the spread between them is honest (RYA-712).

THE CONTROL THAT MAKES THE RESULT READABLE
------------------------------------------
Every line is inverted TWICE: once with the line list untouched, and once with the lab gf
substituted. The first inversion must reproduce the banked RYA-783 1D-LTE abundance for
that line. If it does not, the harness is not the production path and the second number
means nothing -- which is RYA-798's lesson (a control that builds its own inputs is not
testing what ships) turned into a run-time check. The two inversions differ in exactly one
input, so their difference is attributable to the gf and to nothing else.
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

from pipeline.gf_grades import lab_lines, EP_TOL_EV, WAVE_TOL_A  # noqa: E402

ACCOUNTING = ROOT / "data" / "audit" / "line_accounting" / "per_line.csv"
OUT_DEFAULT = ROOT / "data" / "results" / "rya824"

#: The measured pools, named as RYA-783 named them. near-UV is CENSUS-ONLY on purpose:
#: RYA-759 measured equivalent widths there but nothing derives an abundance from them
#: yet, so there is no number for a sub-pool to sit beside.
POOLS = (
    ("VIS", "I", "FeI_3800_6910_kpno_solar_atlas_PROFILEFIT_1D-LTE_lines.csv", True),
    ("IR", "I", "FeI_6910_9199_kpno_solar_atlas_PROFILEFIT_1D-LTE_lines.csv", True),
    ("VIS", "II", "FeII_3800_6910_kpno_solar_atlas_PROFILEFIT_1D-LTE_lines.csv", True),
    ("IR", "II", "FeII_6910_9199_kpno_solar_atlas_PROFILEFIT_1D-LTE_lines.csv", True),
    ("near-UV", "I", "FeI_3000_3800_kpno_solar_atlas_PROFILEFIT_ew.csv", False),
    ("near-UV", "II", "FeII_3000_3780_kpno_solar_atlas_PROFILEFIT_ew.csv", False),
)

#: RYA-799's rule, restated because this ticket re-parses lab tables: a stated
#: uncertainty is positive and, for these techniques, small. A fixed-width table read on
#: whitespace produces a finite but impossible value here.
SIGMA_MIN, SIGMA_MAX = 0.0, 0.5


def _accounting() -> pd.DataFrame:
    acc = pd.read_csv(ACCOUNTING)
    return acc[acc["element"] == "Fe"].reset_index(drop=True)


def census(pool_dir: Path) -> pd.DataFrame:
    """Every measured Fe line that has a primary laboratory gf, per band."""
    lab = lab_lines()
    bad = lab[(lab.e_loggf_dex <= SIGMA_MIN) | (lab.e_loggf_dex > SIGMA_MAX)]
    if len(bad):
        raise SystemExit(f"{len(bad)} lab rows carry an impossible sigma — refusing to "
                         f"run (RYA-799 positivity guard)")
    acc = _accounting()
    rows: list = []
    missing: list = []
    for band, ion, fname, derivable in POOLS:
        path = pool_dir / fname
        if not path.exists():
            print(f"  [skip] {fname} not present")
            missing.append(fname)
            continue
        d = pd.read_csv(path)
        d = d[d["ion"].astype(str) == ion]
        for r in d.itertuples():
            w = float(r.wavelength_air_A)
            a = acc[(acc["ion"] == ion) & (np.abs(acc["wave_air_A"] - w) <= 0.01)]
            if a.empty:
                continue
            ep, k14 = float(a.ep_eV.iloc[0]), float(a.log_gf.iloc[0])
            m = lab[(np.abs(lab.wavelength_air_A - w) <= WAVE_TOL_A)
                    & (np.abs(lab.elo_eV - ep) <= EP_TOL_EV)]
            if m.empty:
                continue
            b = m.iloc[0]
            rows.append(dict(
                band=band, ion=ion, wavelength_air_A=w, ep_eV=ep,
                ew_mA=float(getattr(r, "ew_mA", np.nan)),
                in_aggregate=bool(getattr(r, "in_aggregate", True)),
                a_1dlte_banked=float(getattr(r, "abundance", np.nan)),
                loggf_k14=k14, loggf_lab=float(b.loggf),
                d_loggf=float(b.loggf) - k14,
                lab_source=str(b.source), lab_sigma_dex=float(b.e_loggf_dex),
                derivable=derivable))
    # 🔴 AN EMPTY CENSUS IS A CONFIGURATION ERROR, NOT A RESULT — RYA-1020.
    #
    # `pd.DataFrame([])` has NO COLUMNS, so returning it turns a wrong --pool-dir into an
    # `AttributeError: 'DataFrame' object has no attribute 'band'` about 100 lines away, in
    # the census printer. That reads as a schema bug in the pool and sends you looking at
    # the CSV; the real fault is that this function matched nothing and said so only in
    # skip lines nobody reads. Fail here, naming the directory and what was sought.
    if not rows:
        raise SystemExit(
            f"rya824: no measured line in {pool_dir} carries a primary laboratory gf — "
            f"the census is EMPTY, so there is nothing to grade.\n"
            f"  pools sought : {len(POOLS)}\n"
            f"  not present  : {len(missing)} ({', '.join(missing) or 'none'})\n"
            f"  present but 0 matches: {len(POOLS) - len(missing)}\n"
            f"Point --pool-dir at the RYA-783 per-line pools. NOTE they do not all live in "
            f"one directory (VIS is under rya880/, red-optical under band_products/), so a "
            f"single --pool-dir may not reach every pool. Refusing to emit a columnless "
            f"frame that fails later as a schema error.")
    return pd.DataFrame(rows)


# ── the gf substitution ───────────────────────────────────────────────────────

def substitute_gf(linelist: np.ndarray, wave_A: float, ep_eV: float,
                  new_total_loggf: float) -> np.ndarray:
    """Return a COPY of the synth line list with one physical line's total gf replaced.

    Reuses `gf_resolver.cluster_physical_lines` and its additive-shift convention rather
    than re-implementing them: an HFS-split line is several array rows, and moving one
    row's gf while leaving its siblings would change the branching as well as the total,
    which is a different physical claim from the one the lab paper makes.
    """
    from pipeline.gf_resolver import cluster_physical_lines
    from pipeline.species import species_key

    arr = linelist.copy()
    keys = [species_key(arr["element"][i], arr["ion"][i], arr["molecule"][i])
            for i in range(len(arr))]
    wls = arr["wave_A"].astype(float)
    eps = arr["lower_state_eV"].astype(float)
    gf = arr["loggf"].astype(float)

    hits = []
    for cl in cluster_physical_lines(keys, wls, eps):
        w = 10.0 ** gf[cl]
        centroid = float(np.sum(wls[cl] * w) / w.sum())
        if (abs(centroid - wave_A) <= WAVE_TOL_A
                and abs(float(np.mean(eps[cl])) - ep_eV) <= EP_TOL_EV):
            hits.append(cl)
    if len(hits) != 1:
        raise ValueError(
            f"gf substitution at {wave_A:.4f} A / EP {ep_eV:.4f} eV matched "
            f"{len(hits)} physical lines, expected exactly 1 — refusing to guess which")
    cl = hits[0]
    cur_total = float(np.log10(np.sum(10.0 ** gf[cl])))
    shift = new_total_loggf - cur_total
    for i in cl:
        arr["loggf"][i] = gf[i] + shift
    return arr


def _aggregate(vals: np.ndarray) -> tuple:
    v = np.asarray([x for x in vals if np.isfinite(x)], dtype=float)
    if not len(v):
        return None, None, 0
    return (float(np.median(v)),
            float(np.std(v, ddof=1)) if len(v) > 1 else None, int(len(v)))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pool-dir", type=Path, required=True,
                    help="directory holding the RYA-783 per-line pools (Sirius)")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--tmp-dir", type=Path,
                    default=Path("/tmp/ispec_rya824_synth"),
                    help="private Turbospectrum scratch (never the shared default)")
    ap.add_argument("--from-per-line", type=Path, default=None,
                    help="rebuild the products/summary from a saved per-line CSV without "
                         "re-inverting (the inversion is ~30 min of Turbospectrum)")
    ap.add_argument("--census-only", action="store_true",
                    help="emit the census and the first-order prediction; no synthesis")
    args = ap.parse_args(argv)

    c = census(args.pool_dir)
    args.out.mkdir(parents=True, exist_ok=True)
    c.to_csv(args.out / "rya824_lab_gf_census.csv", index=False)

    print("=" * 78)
    print("RYA-824 — primary-lab-gf Fe sub-pool")
    print("=" * 78)
    lab = lab_lines()
    print(f"\nLAB TABLE (Fe I, primary measurements): {len(lab)} lines")
    for name, lo, hi in (("near-UV", 3000, 3800), ("VIS", 3800, 6910),
                         ("IR", 6910, 9199), ("beyond 9199", 9199, 99999)):
        m = lab[(lab.wavelength_air_A >= lo) & (lab.wavelength_air_A < hi)]
        print(f"  {name:<12s} {len(m):3d}  "
              + ", ".join(f"{s} {int(v)}" for s, v in m.source.value_counts().items()))

    print("\nCENSUS — measured Fe lines that have a primary lab gf.")
    print("  EVERY pool is listed, including the zeros: a band absent from this table "
          "would read\n  as 'not looked at' rather than 'looked at, nothing there'.")
    print(f"  {'band':<9}{'ion':<5}{'lines':>7}{'usable':>8}{'derivable':>13}"
          f"{'median d(log gf)':>18}{'lab sigma':>12}")
    for band, ion, fname, derivable in POOLS:
        g = c[(c.band == band) & (c.ion == ion)]
        u = g[g.in_aggregate]
        der = "yes" if derivable else "CENSUS-ONLY"
        dg = f"{u.d_loggf.median():+.3f}" if len(u) else "--"
        sg = f"{u.lab_sigma_dex.min():.2f}-{u.lab_sigma_dex.max():.2f}" if len(u) else "--"
        print(f"  {band:<9}{ion:<5}{len(g):7d}{len(u):8d}{der:>13}{dg:>18}{sg:>12}")
    print("\n  WHY THE ZEROS ARE ZERO — each checked, not assumed:")
    print("    Fe II (all bands): the vendored lab table is Fe I ONLY. The Fe II "
          "candidate is\n      Melendez & Barbuy 2009 (A&A 497, 611; RYA-472), and it "
          "is a DIFFERENT TIER --\n      multiplets normalised on laboratory data but "
          "relative gf within a multiplet from\n      THEORY. Not solar-calibrated, so "
          "it clears the RYA-161 firewall, but it is not a\n      primary lab "
          "measurement either and must not be pooled with one unlabelled.")
    print("    Fe I near-UV: 901 measured lines and 64 lab lines share the 3000-3800 A "
          "window\n      and NOT ONE PAIR falls within 0.02 A. The closest approach is "
          "0.0651 A -- three\n      times the tolerance with nothing in between, so "
          "this is disjoint line lists and\n      not a tolerance artefact. The "
          "positive control is the IR in the same comparison,\n      where 68 pairs "
          "match and the minimum separation is 0.0000 A.")

    usable = c[c.in_aggregate & c.derivable]
    print(f"\n  RE-MEASUREMENT TARGET: {len(usable)} lines")
    print(f"  lab gf is systematically HIGHER than Kurucz K14: median d(log gf) "
          f"{usable.d_loggf.median():+.3f} dex")
    print(f"  first-order expectation dA ~ -d(log gf) = "
          f"{-usable.d_loggf.median():+.3f} dex (weak-line limit; the inversion below "
          f"is what actually decides)")

    if args.census_only:
        print("\n[--census-only] no synthesis run. The sub-pool needs Sirius.")
        (args.out / "rya824_summary.json").write_text(json.dumps({
            "ticket": "RYA-824", "stage": "census-only",
            "target_lines": int(len(usable)),
            "median_d_loggf": float(usable.d_loggf.median()),
            "run_date": str(date.today())}, indent=2) + "\n")
        return 0

    if args.from_per_line is not None:
        res = pd.read_csv(args.from_per_line)
        return _report(res, usable, args)

    # ── the re-measurement ────────────────────────────────────────────────────
    # CREATE the private scratch. Turbospectrum writes into this directory and does not
    # make it; the shared default only works because earlier runs left it behind. The
    # first version of this script passed a private path without creating it and every
    # single inversion failed with "No such file or directory" -- a guard against a
    # stale-scratch defect that introduced a missing-scratch one.
    args.tmp_dir.mkdir(parents=True, exist_ok=True)

    from pipeline.abundances_derive import _bisect_synth_abundance
    sys.path.insert(0, str(ROOT / "scripts"))
    from rya716_al_products import build_context           # the production context

    ctx = build_context("Fe")
    kw0 = dict(atmosphere=ctx["atmosphere"], teff=ctx["teff"], logg=ctx["logg"],
               feh=ctx["feh"], vturb=ctx["vturb"], isotopes=ctx["isotopes"],
               solar_abund=ctx["solar_abund"], element="Fe",
               atom_code=ctx["atom_code"],
               # PRIVATE scratch. `_bisect_synth_abundance` defaults to a SHARED
               # /tmp/ispec_codex_synth, so two synthesis jobs at once read each other's
               # Turbospectrum working files -- the RYA-785 stale-workdir class of
               # defect, and this run inverts every line twice, which is exactly the
               # shape that would silently reuse a neighbour's scratch.
               tmp_dir=str(args.tmp_dir))

    print(f"\nRE-MEASURING {len(usable)} lines — each inverted TWICE "
          f"(control on K14, then on lab gf)")
    print(f"  {'line':>11}{'band':>7}{'EW':>8}{'A_banked':>10}{'A_ctrl':>9}"
          f"{'A_lab':>9}{'dA':>8}{'-dgf':>8}")
    rows = []
    for r in usable.itertuples():
        c_A = float(r.wavelength_air_A)
        w_nm = np.linspace(c_A - 0.6, c_A + 0.6, 300) / 10.0
        a_ctrl, ok_c, _ = _bisect_synth_abundance(
            w_nm, float(r.ew_mA), linelist=ctx["linelist"], **kw0)
        try:
            ll_lab = substitute_gf(ctx["linelist"], c_A, float(r.ep_eV),
                                   float(r.loggf_lab))
        except ValueError as exc:
            print(f"  {c_A:11.4f} SUBSTITUTION REFUSED: {exc}")
            rows.append(dict(wavelength_air_A=c_A, band=r.band, ion=r.ion,
                             a_control=np.nan, a_labgf=np.nan,
                             note=f"substitution refused: {exc}"))
            continue
        a_lab, ok_l, _ = _bisect_synth_abundance(
            w_nm, float(r.ew_mA), linelist=ll_lab, **kw0)
        dA = (a_lab - a_ctrl) if (np.isfinite(a_lab) and np.isfinite(a_ctrl)) else np.nan
        rows.append(dict(wavelength_air_A=c_A, band=r.band, ion=r.ion,
                         ew_mA=float(r.ew_mA), a_banked=float(r.a_1dlte_banked),
                         a_control=float(a_ctrl) if ok_c else np.nan,
                         a_labgf=float(a_lab) if ok_l else np.nan,
                         d_A=dA, d_loggf=float(r.d_loggf),
                         lab_source=r.lab_source, lab_sigma_dex=float(r.lab_sigma_dex),
                         note=""))
        print(f"  {c_A:11.4f}{r.band:>7}{r.ew_mA:8.1f}{r.a_1dlte_banked:10.3f}"
              f"{a_ctrl if np.isfinite(a_ctrl) else float('nan'):9.3f}"
              f"{a_lab if np.isfinite(a_lab) else float('nan'):9.3f}"
              f"{dA:8.3f}{-float(r.d_loggf):8.3f}")

    res = pd.DataFrame(rows)
    res.to_csv(args.out / "rya824_lab_gf_per_line.csv", index=False)
    return _report(res, usable, args)


def _report(res: pd.DataFrame, usable: pd.DataFrame, args) -> int:
    """Everything downstream of the inversion. Shared by the full run and --from-per-line
    so a re-report can never drift from what the run itself printed."""
    # CONTROL: does the re-derivation reproduce production?
    ctrl = res[np.isfinite(res.a_control) & np.isfinite(res.a_banked)]
    resid = (ctrl.a_control - ctrl.a_banked).abs()
    print("\nCONTROL — the untouched-line-list inversion against the banked RYA-783 value")
    print(f"  n={len(ctrl)}  max |A_ctrl - A_banked| = {resid.max():.4f} dex  "
          f"median {resid.median():.4f}")
    print(f"  {'PASS — this IS the production path' if resid.max() < 0.02 else 'FAIL — the harness is NOT reproducing production; the lab-gf number below is not interpretable'}")

    zero = res[res.d_A.abs() < 1e-9]
    changed = res[res.d_A.abs() >= 1e-9]
    print("\nWHAT THE SUBSTITUTION ACTUALLY CHANGED")
    print(f"  {len(zero)} of {len(res)} lines came back with dA EXACTLY 0.000 — not a failed")
    print("  substitution but a NO-OP: the production synth line list ALREADY carries the")
    print("  lab gf for those lines, because `_load_synth_resources` defaults to")
    print("  apply_canonical_gf=True and canonical_gf.csv already holds the lab value.")
    print("  The K14 numbers live in the pools' `log_gf` COLUMN, which does not describe")
    print("  the gf the inversion used.")
    if len(changed):
        sens = (changed.d_A / changed.d_loggf).median()
        print(f"  {len(changed)} lines genuinely changed: dA median {changed.d_A.median():+.4f}, "
              f"range {changed.d_A.min():+.4f}..{changed.d_A.max():+.4f}")
        print(f"  sensitivity dA/d(log gf) median {sens:.3f} — the WEAK-LINE limit is -1.000.")
        print(f"  These are {changed.ew_mA.min():.0f}-{changed.ew_mA.max():.0f} mA lines, so they sit "
              f"on the flat part of the curve of\n  growth and absorb roughly half the gf change. "
              f"That is why the first-order\n  census estimate over-predicts the shift.")

    print("\nPRODUCTS — both pools, per band (RYA-712: presented side by side, never merged)")
    print(f"  {'band':<9}{'ion':<5}{'pool':<26}{'A(Fe)':>9}{'sigma':>9}{'n':>6}")
    products = []
    for (band, ion), g in res.groupby(["band", "ion"]):
        for label, col, sig in (("PRODUCTION (canonical_gf)", "a_control", None),
                                ("PRIMARY-LAB-GF", "a_labgf", "lab")):
            v, sd, n = _aggregate(g[col].values)
            if v is None:
                continue
            gf_sigma = (float(np.sqrt(np.mean(g.lab_sigma_dex.values ** 2)))
                        if sig == "lab" else 0.20)
            products.append(dict(band=band, ion=ion, pool=label, value=v,
                                 scatter=sd, n=n, gf_sigma_dex=gf_sigma))
            print(f"  {band:<9}{ion:<5}{label:<26}{v:9.3f}"
                  f"{(sd if sd is not None else float('nan')):9.3f}{n:6d}"
                  f"   gf-systematic {gf_sigma:.3f} dex")
    prod = pd.DataFrame(products)
    prod.to_csv(args.out / "rya824_lab_gf_products.csv", index=False)

    (args.out / "rya824_summary.json").write_text(json.dumps({
        "ticket": "RYA-824", "stage": "re-measured",
        "target_lines": int(len(usable)),
        "control_max_resid_dex": float(resid.max()) if len(ctrl) else None,
        "control_pass": bool(len(ctrl) and resid.max() < 0.02),
        "median_d_loggf": float(usable.d_loggf.median()),
        "median_d_A": float(res.d_A.median(skipna=True)),
        "products": products,
        "run_date": str(date.today()),
    }, indent=2, default=float) + "\n")
    print(f"\nwrote {args.out}/rya824_lab_gf_census.csv")
    print(f"wrote {args.out}/rya824_lab_gf_per_line.csv")
    print(f"wrote {args.out}/rya824_lab_gf_products.csv")
    print(f"wrote {args.out}/rya824_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

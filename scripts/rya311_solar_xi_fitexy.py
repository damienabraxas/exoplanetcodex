#!/usr/bin/env python3
"""
scripts/rya311_solar_xi_fitexy.py
=================================
RYA-311 — MEASURE the solar microturbulence and its uncertainty from our own Fe I data,
by the established reduced-equivalent-width slope method with a both-axis (FITEXY) fit.

    python3 scripts/rya311_solar_xi_fitexy.py --star solar --element Fe --ion I

🔴 WHY THIS EXISTS: BOTH SOLAR xi NUMBERS WERE BORROWED.

The central value is PINNED at 1.0 (RYA-196) and the uncertainty that multiplies dA/dxi in
the Type B budget has never been ours either. It was an uncited 0.05 traced to a RYA-158
spec line ("vturb = 0.9 +/- 0.05" -- a central value we do not even use), replaced by
RYA-1089 with Jofre+2014 Table 2's sigma_vmic = 0.18. The 0.18 is defensible as Jofre's OWN
method, but it measures how SEVEN GBS PIPELINES disagree with each other. That ensemble
spread is not the formal error of a single-pipeline analysis, and borrowing it made
sigma_params -- and therefore the whole solar reported bar -- a property of somebody else's
node scatter. Ryan's ruling (2026-08-27): where other scientists have established the
method, we follow the METHOD, not the number.

THE METHOD (all three legs already have a home in this project)
---------------------------------------------------------------
1. The solve. Our own Science Architecture doc: "Microturbulence xi verified by requiring
   zero slope in A(Fe) vs. log(EW/lambda)". Jofre+2014 determines vmic the same way. So xi
   is the value that zeroes d A(Fe I) / d log(EW/lambda).

2. The fit. Mucciarelli 2011 (A&A 528 A44): the abundance-vs-reduced-EW relation must be
   fitted with uncertainties in BOTH coordinates (Press et al. 1992 `fitexy`), because
   sigma_EW enters the abscissa directly and the ordinate through the curve of growth. A
   naive OLS slope is biased towards zero and its error is too small -- Mucciarelli's whole
   point, and the reason `pipeline/fitexy.py` exists rather than a `np.polyfit` call.

3. The error. delta_xi is the formal error on the slope propagated through the solve:
   delta_xi = sigma_slope / |d slope / d xi|, both evaluated at the solution.

4. The line selection. Mucciarelli: discard lines whose EW error exceeds 10% (they bias the
   solve), and sample the full EW range -- xi is constrained by the STRONG end, where the
   curve of growth turns over. A pool of weak lines cannot determine it at all.

WHAT IS AND IS NOT HELD FIXED
-----------------------------
🔴 THE LINE POOL IS FIXED ONCE, BEFORE THE SWEEP. Every trial xi and every EW perturbation
inverts THE SAME lines. `abundances_derive.invert_linemasks` (RYA-311 split it out of
`_ew_to_abundance` for exactly this) inverts an already-selected linemask array, so the
production admission filters run once at the nominal parameters and never again. A pool
that moved with xi would report a changing line set as a changing microturbulence -- the
RYA-875 shape, where a "comparison" let something other than the compared axis vary.

The atmosphere is likewise interpolated once: ATLAS9 keys on Teff/logg/[M/H]/alpha and
carries no microturbulence, so xi enters ONLY as `microturbulence_vel` in the inversion.

⚠️ THE POOL IS PART OF THE MEASUREMENT AND BOTH POOLS ARE PRE-REGISTERED. `--pool graded`
is the headline (Ryan's standing rule: model and parameter work starts from laboratory-gf
lines) and `--pool all` is the declared cross-check on the same EWs. Both are run and both
are reported; whichever is nicer is NOT then chosen (RYA-161).

WHAT THIS SCRIPT DOES NOT DO
----------------------------
It does not write the measured xi back into `config/stars.yaml`. If the solve disagrees
with the pinned 1.0 that is a science finding which MOVES A(Fe) -- it goes to RYA-196 and
the register as a flagged reconciliation, never as a silent central-value change.
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config.constants import PATHS, get_star_params           # noqa: E402
from pipeline.fitexy import fitexy                            # noqa: E402
from pipeline.line_match import match                         # noqa: E402

#: Mucciarelli 2011: lines whose EW is uncertain by more than this fraction bias the solve.
EW_ERR_FRACTION_MAX = 0.10

#: Default xi grid, km/s. Wide enough to bracket any plausible solar value and to measure
#: d slope / d xi far from the root, which is what the error propagation needs.
XI_MIN, XI_MAX, XI_STEP = 0.60, 1.80, 0.10

OUT_DIR = ROOT / "data" / "results" / "rya311"


def _species_rows(ew: pd.DataFrame, element: str, ion: str) -> pd.DataFrame:
    return ew[(ew.element.astype(str).str.strip() == element)
              & (ew.ion.astype(str).str.strip() == ion)].copy()


def _lab_graded_mask(rows: pd.DataFrame, species: str) -> np.ndarray:
    """Which measured lines carry a LABORATORY oscillator strength (canonical_gf LAB tier).

    Matching goes through `pipeline.line_match`, never a rounded wavelength: `canonical_gf`
    and the measured pools store the same line to different precision and a 2-dp key splits
    matched pairs outright (RYA-1033).
    """
    cg = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    lab = cg[(cg.species == species)
             & cg.gf_tier.astype(str).str.contains("LAB", na=False)].reset_index(drop=True)
    if lab.empty:
        raise SystemExit(f"no LAB-tier {species} rows in canonical_gf -- refusing to call "
                         f"a pool 'graded' when nothing in it is graded.")
    res = match(rows.wavelength_air_A.to_numpy(float),
                lab.wavelength_air_A.to_numpy(float))
    if res.ambiguous:
        raise SystemExit(
            f"{len(res.ambiguous)} measured {species} line(s) match more than one LAB-tier "
            f"canonical_gf row and no EW-pool excitation potential is available to separate "
            f"them. Wavelength alone does not identify an Fe line (RYA-1033); adjudicate "
            f"these rather than letting the solve pick one:\n  "
            + "\n  ".join(f"{w:.4f} -> {c}" for w, c in res.ambiguous[:8]))
    return res.resolved


def solve(star: str, element: str, ion: str, pool: str,
          ew_ceiling_mA: float | None,
          xi_grid: np.ndarray) -> tuple[dict, pd.DataFrame]:
    import pipeline.abundances_derive as ad

    species = f"{element} {ion}"
    sp = get_star_params(star)
    params0 = {"teff_K": float(sp["teff"]), "logg": float(sp["logg"]),
               "feh": float(sp.get("feh_ref", 0.0)),
               "vturb_kms": float(sp.get("xi", 1.0))}
    print(f"\n=== RYA-311  {species}  {star}  pool={pool} ===")
    print(f"  nominal params: Teff {params0['teff_K']} / logg {params0['logg']} / "
          f"[Fe/H] {params0['feh']} / xi_pinned {params0['vturb_kms']} (RYA-196)")

    # ── The EW pool: whatever the production path uses for THIS species ──────────
    ew_all = ad._load_solar_ews()
    rows = _species_rows(ew_all, element, ion)
    if rows.empty:
        raise SystemExit(f"no {species} rows in the production solar EW pool")
    n_start = len(rows)
    cuts: list[dict] = []

    def _cut(mask_drop: np.ndarray, reason: str) -> None:
        for _, r in rows[mask_drop].iterrows():
            cuts.append({"wavelength_air_A": float(r.wavelength_air_A),
                         "ew_mA": float(r.ew_mA),
                         "ew_err_mA": float(r.ew_err_mA),
                         "reason": reason})

    if pool == "graded":
        keep = _lab_graded_mask(rows, species)
        _cut(~keep, "not LAB-tier in canonical_gf (--pool graded)")
        rows = rows[keep]
    elif pool != "all":
        raise SystemExit(f"unknown --pool {pool!r}")

    # Mucciarelli's cut, applied to whatever pool survived above.
    frac = rows.ew_err_mA.to_numpy(float) / rows.ew_mA.to_numpy(float)
    keep = frac <= EW_ERR_FRACTION_MAX
    _cut(~keep, f"EW error > {EW_ERR_FRACTION_MAX:.0%} (Mucciarelli 2011)")
    rows = rows[keep]

    if ew_ceiling_mA is not None:
        keep = rows.ew_mA.to_numpy(float) <= ew_ceiling_mA
        _cut(~keep, f"EW > {ew_ceiling_mA:.0f} mA (production vmic COG ceiling, RYA-330)")
        rows = rows[keep]

    print(f"  pool: {n_start} measured {species} lines -> {len(rows)} after "
          f"{len(cuts)} cut(s)")
    if len(rows) < 5:
        raise SystemExit(f"{len(rows)} lines left -- a slope and its error cannot be "
                         f"measured from this. Refusing rather than quoting one.")

    # ── Fix the pool ONCE at the nominal parameters ──────────────────────────────
    atm = ad._load_atmosphere(params0["teff_K"], params0["logg"],
                              params0["feh"], params0["vturb_kms"])
    linemasks, _, _, _ = ad._ew_to_abundance(rows, params0, atm, star_id=star)
    ew0 = np.asarray(linemasks["ew"], dtype=float).copy()
    err0 = np.asarray(linemasks["ew_err"], dtype=float).copy()
    wl = np.asarray(linemasks["wave_A"], dtype=float).copy()
    if np.any(err0 <= 0):
        raise SystemExit("a line reached the fit with a zero EW error -- FITEXY would give "
                         "it infinite weight. Fix the EW pool, do not floor it here.")
    print(f"  production admission filters kept {len(linemasks)} of {len(rows)}")

    rew = np.log10(ew0 / 1000.0 / wl)
    sigma_rew = err0 / (ew0 * np.log(10.0))

    def _abund(vturb: float, ew: np.ndarray) -> np.ndarray:
        lm = linemasks.copy()
        lm["ew"] = ew
        p = dict(params0, vturb_kms=float(vturb))
        _, a, _, _ = ad.invert_linemasks(lm, p, atm)
        return np.asarray(a, dtype=float)

    # ── Sweep ────────────────────────────────────────────────────────────────────
    sweep = []
    for xi in xi_grid:
        a0 = _abund(xi, ew0)
        ap = _abund(xi, ew0 + err0)
        am = _abund(xi, ew0 - err0)
        ok = np.isfinite(a0) & np.isfinite(ap) & np.isfinite(am)
        # sigma_A from varying the EW by +/-1 sigma_EW -- Mucciarelli 2011's prescription.
        sig_a = np.abs(ap - am) / 2.0
        ok &= sig_a > 0
        if ok.sum() < 5:
            raise SystemExit(f"xi={xi}: only {int(ok.sum())} lines inverted -- refusing")
        fit = fitexy(rew[ok], a0[ok], sigma_rew[ok], sig_a[ok])
        # The unweighted OLS slope on the SAME points, at every grid node. This is what
        # the production vmic bisection fits (`abundances_derive._compute_rew_slope` is a
        # bare `np.polyfit`), so the two roots below are a like-for-like statement of what
        # the both-axis fit changes -- not an argument that it should.
        ols = float(np.polyfit(rew[ok], a0[ok], 1)[0])
        sweep.append({"xi": float(xi), "n": int(ok.sum()),
                      "slope": fit.slope, "sigma_slope": fit.sigma_slope,
                      "sigma_slope_scaled": fit.sigma_slope_scaled,
                      "ols_slope": ols,
                      "chi2_red": fit.chi2_red,
                      "A_mean": float(np.mean(a0[ok])),
                      "A_std": float(np.std(a0[ok], ddof=1))})
        print(f"  xi={xi:.2f}  n={ok.sum():3d}  slope={fit.slope:+.4f} "
              f"+/- {fit.sigma_slope:.4f} (formal)  OLS={ols:+.4f}  "
              f"chi2_red={fit.chi2_red:8.2f}  A={np.mean(a0[ok]):.4f}")

    s = pd.DataFrame(sweep)
    # slope(xi) is smooth and near-linear over this range; a quadratic through the grid
    # gives both the root and the derivative without pretending the relation is a line.
    coef = np.polyfit(s.xi, s.slope, 2)
    roots = np.roots(coef)
    real = [float(r.real) for r in roots if abs(r.imag) < 1e-9
            and xi_grid.min() <= r.real <= xi_grid.max()]
    if len(real) != 1:
        # 🔴 NOT AN ERROR -- A RESULT. A pool whose reduced-EW slope never reaches zero
        # inside a physically sane xi range has told us something: microturbulence alone
        # cannot flatten it. Refusing to publish a xi is right; refusing to publish the
        # measurement that shows why would be RYA-833 (an absence recorded as nothing).
        outside = [float(r.real) for r in roots if abs(r.imag) < 1e-9]
        print(f"\n  NO SOLVE: slope(xi) has {len(real)} root(s) inside "
              f"{xi_grid.min()}-{xi_grid.max()} km/s "
              f"(quadratic roots anywhere: {[round(r, 3) for r in outside]}).")
        print(f"  The slope runs {s.slope.iloc[0]:+.4f} -> {s.slope.iloc[-1]:+.4f} across "
              f"the grid without crossing zero. This pool does not determine xi.")
        return {"ticket": "RYA-311", "star": star, "species": species, "pool": pool,
                "ew_ceiling_mA": ew_ceiling_mA,
                "ew_err_fraction_max": EW_ERR_FRACTION_MAX,
                "n_measured": int(n_start), "n_after_cuts": int(len(rows)),
                "n_fitted": int(sweep[0]["n"]),
                "xi_grid": [float(x) for x in xi_grid],
                "root_found": False,
                "roots_anywhere": outside,
                "slope_first": float(s.slope.iloc[0]),
                "slope_last": float(s.slope.iloc[-1]),
                "verdict": "NO_SOLVE -- the reduced-EW slope never reaches zero in the "
                           "swept xi range; this pool does not constrain microturbulence",
                "sweep": sweep, "cuts": cuts}, pd.DataFrame()
    xi_hat = real[0]
    dslope_dxi = float(np.polyval(np.polyder(coef), xi_hat))

    # Where the SAME sweep would have landed on an unweighted OLS slope -- the production
    # bisection's criterion. Reported as the comparison, never substituted for the answer.
    ols_coef = np.polyfit(s.xi, s.ols_slope, 2)
    ols_real = [float(r.real) for r in np.roots(ols_coef) if abs(r.imag) < 1e-9
                and xi_grid.min() <= r.real <= xi_grid.max()]
    xi_ols = ols_real[0] if len(ols_real) == 1 else None

    # ── Verify at the solution rather than trusting the interpolation ────────────
    a_hat = _abund(xi_hat, ew0)
    ap = _abund(xi_hat, ew0 + err0)
    am = _abund(xi_hat, ew0 - err0)
    sig_a = np.abs(ap - am) / 2.0
    ok = np.isfinite(a_hat) & np.isfinite(sig_a) & (sig_a > 0)
    fit_hat = fitexy(rew[ok], a_hat[ok], sigma_rew[ok], sig_a[ok])
    sigma_xi = abs(fit_hat.sigma_slope / dslope_dxi)
    sigma_xi_scaled = abs(fit_hat.sigma_slope_scaled / dslope_dxi)
    # The Delta-chi2 = 1 interval is not symmetric (chi2(b) is not a parabola once
    # b^2 sigma_x^2 is in the denominator). Both sides are recorded; the headline
    # delta_xi is their mean, and the record says how far apart they were.
    sigma_xi_lo = abs(fit_hat.sigma_slope_lo / dslope_dxi)
    sigma_xi_hi = abs(fit_hat.sigma_slope_hi / dslope_dxi)

    # OLS on the same points, reported ONLY to show what the naive fit would have said.
    ols_slope, _ = np.polyfit(rew[ok], a_hat[ok], 1)
    a_pin = _abund(params0["vturb_kms"], ew0)
    ok_pin = np.isfinite(a_pin)

    res = {
        "ticket": "RYA-311",
        "star": star, "species": species, "pool": pool,
        "ew_source": "production solar EW pool via abundances_derive._load_solar_ews()",
        "ew_err_fraction_max": EW_ERR_FRACTION_MAX,
        "ew_ceiling_mA": ew_ceiling_mA,
        "n_measured": int(n_start),
        "n_after_cuts": int(len(rows)),
        "n_fitted": int(ok.sum()),
        "rew_min": float(rew[ok].min()), "rew_max": float(rew[ok].max()),
        "ew_min_mA": float(ew0[ok].min()), "ew_max_mA": float(ew0[ok].max()),
        "root_found": True,
        "xi_grid": [float(x) for x in xi_grid],
        "xi_measured": float(xi_hat),
        "delta_xi_formal": float(sigma_xi),
        "delta_xi_lo": float(sigma_xi_lo),
        "delta_xi_hi": float(sigma_xi_hi),
        "delta_xi_chi2_scaled": float(sigma_xi_scaled),
        "slope_at_solution": float(fit_hat.slope),
        "sigma_slope_formal": float(fit_hat.sigma_slope),
        "sigma_slope_chi2_scaled": float(fit_hat.sigma_slope_scaled),
        "chi2_red_at_solution": float(fit_hat.chi2_red),
        "dslope_dxi": dslope_dxi,
        "ols_slope_at_solution": float(ols_slope),
        "xi_ols_root": (float(xi_ols) if xi_ols is not None else None),
        "xi_pinned": float(params0["vturb_kms"]),
        "A_at_measured_xi": float(np.mean(a_hat[ok])),
        "A_at_pinned_xi": float(np.mean(a_pin[ok_pin])),
        "A_shift_measured_minus_pinned": float(np.mean(a_hat[ok]) - np.mean(a_pin[ok_pin])),
        "sweep": sweep,
        "cuts": cuts,
    }
    print(f"\n  xi = {xi_hat:.4f} +/- {sigma_xi:.4f} km/s  (formal FITEXY error; "
          f"chi2-scaled alternative {sigma_xi_scaled:.4f})")
    print(f"  residual slope at the solution: {fit_hat.slope:+.5f} "
          f"(OLS on the same points: {ols_slope:+.5f})")
    print(f"  the SAME sweep solved on an unweighted OLS slope instead: xi_OLS = "
          f"{xi_ols if xi_ols is None else round(xi_ols, 4)}  "
          f"(what the production bisection's criterion would give)")
    print(f"  A({species}) {np.mean(a_pin[ok_pin]):.4f} at the pinned xi -> "
          f"{np.mean(a_hat[ok]):.4f} at the measured xi  "
          f"({np.mean(a_hat[ok]) - np.mean(a_pin[ok_pin]):+.4f} dex)")

    per_line = pd.DataFrame({
        "wavelength_air_A": wl[ok], "ew_mA": ew0[ok], "ew_err_mA": err0[ok],
        "rew": rew[ok], "sigma_rew": sigma_rew[ok],
        "A_at_measured_xi": a_hat[ok], "sigma_A_from_ew": sig_a[ok],
    })
    return res, per_line


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--star", default="solar")
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--ion", default="I")
    ap.add_argument("--pool", default="graded", choices=("graded", "all"),
                    help="graded = LAB-tier gf only (headline); all = every measured line")
    ap.add_argument("--ew-ceiling-mA", type=float, default=None,
                    help="optional saturation ceiling; the production vmic bisection uses "
                         "100 mA (RYA-330), which Mucciarelli's 'sample the strong end' "
                         "argues against. Declared sensitivity, not the default.")
    ap.add_argument("--xi-min", type=float, default=XI_MIN)
    ap.add_argument("--xi-max", type=float, default=XI_MAX)
    ap.add_argument("--xi-step", type=float, default=XI_STEP)
    ap.add_argument("--out", default=None, help="JSON path (default data/results/rya311/)")
    a = ap.parse_args()

    grid = np.round(np.arange(a.xi_min, a.xi_max + 0.5 * a.xi_step, a.xi_step), 3)
    res, per_line = solve(a.star, a.element, a.ion, a.pool, a.ew_ceiling_mA, grid)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tag = f"{a.star}_{a.element}{a.ion}_{a.pool}"
    if a.ew_ceiling_mA is not None:
        tag += f"_ceil{int(a.ew_ceiling_mA)}"
    if (a.xi_min, a.xi_max, a.xi_step) != (XI_MIN, XI_MAX, XI_STEP):
        tag += f"_xi{a.xi_min:g}-{a.xi_max:g}"
    out = Path(a.out) if a.out else OUT_DIR / f"rya311_xi_{tag}.json"
    out.write_text(json.dumps(res, indent=2))
    print(f"\n  -> {out.relative_to(ROOT)}")
    if not per_line.empty:
        csv = out.with_name(out.stem + "_per_line.csv")
        per_line.to_csv(csv, index=False)
        print(f"  -> {csv.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

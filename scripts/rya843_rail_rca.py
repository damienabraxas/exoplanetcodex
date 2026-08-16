#!/usr/bin/env python3
"""Why does a synthesis flux-fit run to its abundance bound? — RYA-843.

    python3 scripts/rya843_rail_rca.py --step atlas     # cheap, no synthesis
    python3 scripts/rya843_rail_rca.py --step curves    # ~25 s per line (Turbospectrum)
    python3 scripts/rya843_rail_rca.py --step report    # rebuilds from the CSVs

THE QUESTION
------------
RYA-837's first NIR run reported A(Fe I) = 7.588 +/- 1.977 because several lines fitted
to ~12.49 — `a_hi = A_solar + 5`, the search rail. A fit parked on its bound is not a
measurement, and `_fit_synth_flux`'s `edge_pinned` detector let three of them through
because its tolerance is 1e-2 and they sat 0.017 away.

RYA-843 forbids setting that tolerance by whatever makes this product look clean
(RYA-161) and asks for the ROOT CAUSE first: a line the fitter runs to the bound is
more likely a line-selection failure than a fitter failure (RYA-842).

THE THREE STEPS, AND WHY THEY ARE SEPARATE
------------------------------------------
`atlas` needs no synthesis at all and answers most of it. `curves` is the expensive
half. Splitting them means a bug in the reporting cannot cost the Turbospectrum time —
the RYA-836 lesson, where `~series.fillna(False)` on an object column crashed the report
AFTER all 122 fits and cost nothing because the per-line CSV was already written.

WHAT THE STEPS FOUND (all measured; see docs/science/rya843_rail_rca.md)
-----------------------------------------------------------------------
1. TELLURIC, 5 of the original 6. Observed core depths 0.987-1.001 inside H2O
   11120-11560 A — flux at or below ZERO, no stellar profile to fit. The EW route
   already excluded these via `telluric_reason`; the synthesis route inherited
   `select_lines`, which knows nothing about tellurics. Fixed at SELECTION, which makes
   the detector tolerance moot for them rather than tuning it to hide them.
2. NOT the gf. Linelist `loggf` == `canonical_gf.log_gf` to 0.000 for railers and
   controls alike, all VALD3/single_source. RYA-834's redward extension substituted
   nothing here.
3. NOT saturation, which was the standing hypothesis and is REFUTED by `curves`: the
   synthetic core reaches the observed depth at a sane abundance (11593.588: observed
   0.588, synthetic 0.595 at A=8.996) and chi2 keeps falling well past it. The core is
   not what the fitter is looking at.
4. THE WINDOW BASELINE IS. At +/-1.4 A the window holds 150-270 pixels and the core a
   handful, so chi2 is dominated by how well the BASELINE is reproduced. Every railing
   line sits in a window whose Kitt Peak median flux is 0.73-0.95 against a synthesis
   normalised to unity, and A(Fe) is the only free parameter available to close that
   gap. The rail is a normalisation failure wearing a fitter's clothes.
5. 🔴 AND THE DEFECT IS BIGGER THAN THE RAIL. Two lines that fit at plausible values and
   entered the aggregate — 11119.795 -> 7.833 and 11572.523 -> 7.979 — have chi2 that
   moves 2.2% and 1.4% across EIGHT DEX of iron. They measured nothing and landed
   somewhere believable by luck. The railed fits are the loud half of one defect:
   UNCONSTRAINED FITS ARE ACCEPTED, and only the ones that land far away are visible.

WHY THIS SCRIPT DOES NOT SET A THRESHOLD
----------------------------------------
It reports three candidate metrics per line and picks none. `frac` is the fractional
chi2 rise to the ends of the caller's bracket; it separates cleanly here, but the
bracket is [A_solar-3, A_solar+5], so a minimum sitting near A=10 is close to the +5 end
and has little room to rise — the metric is part sharpness-of-minimum and part
proximity-to-bracket-edge, and the second is a plausibility prior tied to A_solar coming
in through the back door. `sigma_A` is bracket-free and in DEX, but measures PRECISION
not ACCURACY: it accepts 11973.046 at A=10.056, whose minimum is genuinely sharp and
genuinely in the wrong place. `red_chi2` is what the production handler already gates on
and its constant does NOT transfer — the good line 10145.561 fits at red_chi2 = 72.

Nine lines cannot settle that, and RYA-161 forbids setting a cut from one product.
RYA-847 owns the criterion and sets it from a sweep across every synthesis band.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

OUT = ROOT / "data" / "results" / "rya843"
#: The per-line output of the NIR run this ticket is diagnosing, committed as EVIDENCE
#: and deliberately not as a product: RYA-837 and RYA-843 both publish nothing, because
#: the acceptance defect below is unfixed until RYA-847. It is kept out of
#: `data/results/band_products/` so nothing downstream can mistake it for one.
LINES_CSV = OUT / "rya843_nir_run_per_line.csv"
LINELIST = ROOT / "data" / "linelists" / "ispec_ir_9200_13000" / "atomic_lines.tsv"
CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"

#: The NIR entry of `derive_band_products.SYNTH_BANDS`. Named, not re-chosen.
HALF_WIDTH_A = 1.40
INSTRUMENT = "kpno_solar_atlas"
BAND_LO_A, BAND_HI_A = 10000.0, 12935.0

#: Probed abundances, as OFFSETS from the run's own `solar_A` so the grid follows the
#: bracket rather than pinning a literal. The bracket is [solar-3, solar+5]; the grid
#: starts at -1 because nothing in this band fits below it and each point costs a
#: synthesis. Dense near solar where a real minimum lives, sparse out at the rail.
OFFSETS = (-1.0, -0.5, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0)

#: The five that survived the telluric fix, and four controls. The controls are NOT all
#: clean: 11119.795 and 11572.523 sit in depressed windows and fit to plausible values,
#: and they are here precisely because they refuted the tidy story twice — first that
#: continuum depression separates railers from good lines, then that a constraint metric
#: only ever flags lines that look wrong.
RAILERS = (11593.588, 11689.972, 11973.046, 12638.703, 12648.741)
CONTROLS = (10145.561, 12342.916, 11119.795, 11572.523)


def _ctx():
    """RYA-759's solar context, built exactly as the band-product route builds it."""
    from pipeline.nearuv_synth import build_solar_context, gf_provenance
    cat = pd.read_csv(ROOT / "data" / "catalog" / "instrument_catalog.csv")
    row = cat[cat.iloc[:, 0].astype(str) == INSTRUMENT]
    if row.empty:
        raise SystemExit(f"{INSTRUMENT!r} absent from instrument_catalog.csv")
    R = float(row.iloc[0]["resolving_power_max"])
    ll = pd.read_csv(LINELIST, sep="\t", usecols=["wave_A"], low_memory=False)
    w = ll.wave_A[(ll.wave_A >= BAND_LO_A) & (ll.wave_A <= BAND_HI_A)]
    prov = gf_provenance(float(w.min()), float(w.max()))
    return build_solar_context("Fe", R, linelist_file=str(LINELIST),
                               apply_canonical_gf=prov["apply_canonical_gf"])


# ── step: atlas ───────────────────────────────────────────────────────────────

def step_atlas() -> pd.DataFrame:
    """Everything answerable from the line list and the atlas. No synthesis."""
    from rya759_nearuv_synth import _kp_segments
    from rya759_nearuv_fe_product import _load_kp_window
    from measure_band_ew import telluric_reason

    ll = pd.read_csv(LINELIST, sep="\t", low_memory=False)
    cg = pd.read_csv(CANONICAL, low_memory=False)
    allw = ll.wave_A.to_numpy(float)
    segs = _kp_segments()

    rows = []
    for w in sorted(RAILERS + CONTROLS):
        m = ll[(ll.wave_A - w).abs() < 0.005]
        c = cg[((cg.wavelength_air_A.astype(float) - w).abs() < 0.05) &
               (cg.species.astype(str).str.strip() == "Fe I")]
        d = np.abs(allw - w)
        d = d[d > 1e-4]
        ow, of = _load_kp_window(segs, w, HALF_WIDTH_A + 0.4)
        core = np.abs(ow - w) <= 0.30
        imin = int(np.argmin(np.where(core, of, np.inf)))
        rows.append(dict(
            wave_A=w, klass="RAIL" if w in RAILERS else "ctrl",
            loggf_list=float(m.loggf.iloc[0]) if len(m) else np.nan,
            loggf_canonical=float(c.log_gf.iloc[0]) if len(c) else np.nan,
            gf_reference=str(c.loggf_reference.iloc[0]) if len(c) else "NO ROW",
            theo_depth=float(m.theoretical_depth.iloc[0]) if len(m) else np.nan,
            ep_eV=float(m.lower_state_eV.iloc[0]) if len(m) else np.nan,
            nearest_neighbour_A=float(d.min()) if d.size else np.nan,
            obs_core_depth=1.0 - float(of[imin]),
            core_offset_A=float(ow[imin] - w),
            kp_window_median=float(np.median(of)),
            kp_window_p90=float(np.percentile(of, 90)),
            telluric=telluric_reason(w, INSTRUMENT) or "",
        ))
    df = pd.DataFrame(rows)
    df["gf_delta"] = df.loggf_canonical - df.loggf_list
    return df


# ── step: curves ──────────────────────────────────────────────────────────────

def step_curves() -> pd.DataFrame:
    """chi2 as a function of A(Fe), using THE FITTER'S OWN OBJECTIVE.

    `_fit_synth_flux` builds its `chi2` closure locally, so it is captured by spying on
    `minimize_scalar` in that module: the spy records the closure and returns a stub, so
    the bounded search never runs and we drive the objective ourselves. Rebuilding chi2
    here instead would be a copy that agrees with itself while the product does something
    else — the RYA-845 failure, one layer out.

    Spying on `_synth_flux_at_abund` too gets the synthetic core depth out of the SAME
    evaluations, so "the model saturates" is measured rather than argued, at no extra
    synthesis cost.
    """
    import pipeline.abundances_derive as ad
    from rya759_nearuv_synth import _kp_segments
    from rya759_nearuv_fe_product import _load_kp_window

    ctx = _ctx()
    segs = _kp_segments()
    a_solar = float(ctx["solar_A"])
    a_lo, a_hi = max(a_solar - 3.0, 1.0), a_solar + 5.0
    print(f"[bounds] solar_A={a_solar:.4f}  a_lo={a_lo:.4f}  a_hi={a_hi:.4f}", flush=True)

    tmp = "/tmp/ispec_synth_rya843"
    Path(tmp).mkdir(parents=True, exist_ok=True)

    cap: dict = {}
    orig_min, orig_syn = ad.minimize_scalar, ad._synth_flux_at_abund

    class _Stub:
        def __init__(self, x):
            self.x = x

    def spy_min(fun, bounds=None, **kw):
        cap["fun"] = fun
        return _Stub(a_solar)

    def spy_syn(sw, trial_A=None, **kw):
        sf = orig_syn(sw, trial_A=trial_A, **kw)
        cap.setdefault("synth", []).append((np.asarray(sw, float), np.asarray(sf, float)))
        return sf

    ad.minimize_scalar, ad._synth_flux_at_abund = spy_min, spy_syn
    rows = []
    try:
        for w in sorted(RAILERS + CONTROLS):
            cap.pop("fun", None)
            cap["synth"] = []
            ow_A, of = _load_kp_window(segs, w, HALF_WIDTH_A + 0.4)
            ad._fit_synth_flux(
                ow_A / 10.0, np.asarray(of, float), ctx["atmosphere"], ctx["teff"],
                ctx["logg"], ctx["feh"], ctx["vturb"], ctx["linelist"], ctx["isotopes"],
                ctx["solar_abund"], "Fe", int(ctx["atom_code"]),
                (w - HALF_WIDTH_A) / 10.0, (w + HALF_WIDTH_A) / 10.0, a_lo, a_hi,
                float(ctx["resolving_power"]), float(ctx["macroturbulence"]),
                float(ctx["vsini"]), tmp_dir=tmp)
            if "fun" not in cap:
                raise RuntimeError(f"objective never built for {w} — the fitter bailed "
                                   f"before minimize_scalar; nothing to probe")
            fun = cap["fun"]
            n_pix = int(((ow_A >= w - HALF_WIDTH_A) & (ow_A <= w + HALF_WIDTH_A)).sum())
            core = np.abs(ow_A - w) <= 0.05
            obs_depth = 1.0 - float(np.min(of[core])) if core.any() else np.nan
            klass = "RAIL" if w in RAILERS else "ctrl"
            print(f"\n{klass} {w:10.3f}  n_pix {n_pix}  observed core {obs_depth:.3f}",
                  flush=True)
            for off in OFFSETS:
                a = a_solar + off
                cap["synth"] = []
                c = float(fun(a))
                sd = np.nan
                if cap["synth"]:
                    sw, sf = cap["synth"][-1]
                    sel = np.abs(sw * 10.0 - w) <= 0.05
                    if sel.any():
                        sd = 1.0 - float(np.min(sf[sel]))
                rows.append(dict(wave_A=w, klass=klass, A=a, chi2=c, synth_depth=sd,
                                 obs_depth=obs_depth, n_pix=n_pix,
                                 a_lo=a_lo, a_hi=a_hi))
                print(f"    A={a:7.3f}  chi2={c:12.1f}  synth_depth={sd:.3f}", flush=True)
    finally:
        ad.minimize_scalar, ad._synth_flux_at_abund = orig_min, orig_syn
    return pd.DataFrame(rows)


# ── step: report ──────────────────────────────────────────────────────────────

def step_report(curves: pd.DataFrame) -> pd.DataFrame:
    """Three candidate constraint metrics per line. Sets no threshold — see RYA-847."""
    prod = pd.read_csv(LINES_CSV) if LINES_CSV.exists() else None
    rows = []
    for w, g in curves.groupby("wave_A"):
        g = g.sort_values("A").reset_index(drop=True)
        A, c = g.A.to_numpy(), g.chi2.to_numpy()
        i = int(np.argmin(c))
        a0, c0 = float(A[i]), float(c[i])
        dof = max(int(g.n_pix.iloc[0]) - 1, 1)
        red = c0 / dof
        # sigma_A: 1 sigma on A(Fe) from the curvature, with chi2 RESCALED so red_chi2=1.
        # The rescale is required because `sig_scalar` is a 0.01 placeholder, not a real
        # per-pixel error, so a raw Delta-chi2 = 1 means nothing. Delta-chi2 is a
        # DIFFERENCE, so the depressed-baseline offset that makes absolute red_chi2
        # untransferable cancels exactly — that is why this one crosses bands and the
        # absolute threshold does not.
        sig = []
        for step in (-0.5, 0.5):
            a = a0 + step
            if a < A.min() or a > A.max():
                continue
            d = max((float(np.interp(a, A, c)) - c0) / red, 1e-9)
            sig.append(abs(step) / np.sqrt(d))
        pr = prod[(prod.wavelength_air_A - w).abs() < 0.01] if prod is not None else None
        rows.append(dict(
            wave_A=w, klass=g.klass.iloc[0], n_pix=int(g.n_pix.iloc[0]),
            chi2_min=c0, red_chi2=red, at_A=a0,
            interior_minimum=bool(0 < i < len(c) - 1),
            frac_rise_lo=(float(c[0]) - c0) / c0, frac_rise_hi=(float(c[-1]) - c0) / c0,
            frac_rise_weaker=min((float(c[0]) - c0) / c0, (float(c[-1]) - c0) / c0),
            sigma_A=max(sig) if sig else np.nan,
            product_A=float(pr.abundance.iloc[0]) if pr is not None and len(pr) else np.nan,
            in_aggregate=(bool(pr.in_aggregate.iloc[0])
                          if pr is not None and len(pr) else None)))
    df = pd.DataFrame(rows).sort_values("sigma_A", ascending=False)

    print("\n  %-11s %-5s %6s %9s %8s %8s %9s %8s  %s" %
          ("lambda", "class", "n_pix", "red_chi2", "at A", "sigma_A", "frac_weak",
           "prod A", "interior?"))
    for r in df.itertuples():
        print("  %10.3f %-5s %6d %9.1f %8.3f %8.3f %8.1f%% %8.3f  %s%s"
              % (r.wave_A, r.klass, r.n_pix, r.red_chi2, r.at_A, r.sigma_A,
                 100 * r.frac_rise_weaker, r.product_A,
                 "yes" if r.interior_minimum else "NO — min ON THE BOUND",
                 "   <== IN THE AGGREGATE" if r.in_aggregate else ""))
    print("\n  No threshold is set here. The three metrics disagree about which lines "
          "are\n  acceptable, nine lines cannot settle it, and RYA-161 forbids cutting "
          "from one\n  product. RYA-847 sets it from a sweep across every synthesis band.")
    return df


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--step", choices=("atlas", "curves", "report", "all"),
                    default="report")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    if a.step in ("atlas", "all"):
        d = step_atlas()
        d.to_csv(out / "rya843_line_context.csv", index=False)
        with pd.option_context("display.width", 200, "display.max_columns", 30):
            print(d.to_string(index=False))
        print(f"\nwrote {out / 'rya843_line_context.csv'}")

    if a.step in ("curves", "all"):
        d = step_curves()
        d.to_csv(out / "rya843_chi2_curves.csv", index=False)
        print(f"\nwrote {out / 'rya843_chi2_curves.csv'}")

    if a.step in ("report", "all"):
        p = out / "rya843_chi2_curves.csv"
        if not p.exists():
            raise SystemExit(f"{p} missing — run `--step curves` first (it needs "
                             f"Turbospectrum; ~25 s per line).")
        step_report(pd.read_csv(p)).to_csv(out / "rya843_constraint.csv", index=False)
        print(f"wrote {out / 'rya843_constraint.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

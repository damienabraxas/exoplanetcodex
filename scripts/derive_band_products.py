#!/usr/bin/env python3
"""Derive per-(band × treatment) abundance products from measured EWs — RYA-712/713.

    python3 scripts/derive_band_products.py --element Fe --ion I --lo 6910 --hi 9199

Ryan: *"run the IR products now, and the UV"* · *"the two engines should never be added
together and presented, they are separate data products"* · *"diff error bars for diff
bands."*

WHAT THIS PRODUCES
------------------
One row per **treatment**, never merged:

  * **1D-LTE**   — the measured EW inverted through a curve of growth.
  * **ENGINE-A** — 1D-LTE plus the Bergemann/MPIA per-line departure correction.
  * **ENGINE-B** — the spectrum fitted directly by `pipeline.measure.SynthesisHandler`.
                   NO equivalent width is inverted anywhere on this path: 1D-LTE and
                   ENGINE-A invert a measured EW, ENGINE-B fits flux. The EW table supplies
                   only the LINE SET and a per-line strength hint that sets each fitting
                   window; it never supplies the answer.

                   WHICH synthesis is selected by `--engine-b-deck` (RYA-784):
                     `ts-lte`       Turbospectrum LTE flux-fit. What the register calls
                                    Engine B in production for all 27 species. Default.
                     `gerber-nlte`  the TS-native Gerber departure deck (RYA-533). NOT yet
                                    provisioned for Fe — RYA-785 pulls atom.fe607a. Selecting
                                    it today fails LOUDLY rather than silently returning the
                                    LTE number under an NLTE label, which is the one
                                    confusion this flag exists to prevent.

Each carries its own value, line count, and a full `ErrorBudget` with statistical and
systematic kept apart. There is no combined product: `pipeline.band_products` has no
`combine()` and `ErrorBudget.total()` returns a pair.

WHERE THE ABUNDANCES COME FROM
------------------------------
The EWs are whatever the band's controlled handler produced (`*_PROFILEFIT_ew.csv`). The
EW→abundance inversion is the project's own `_bisect_synth_abundance` — Turbospectrum
curve-of-growth, the same engine the EW path uses — not a formula written here.

WHAT IT REFUSES
---------------
* A band whose handler has not passed its optical control (control/frontier rule).
* A line that the measurement quarantined: those are carried into the output with their
  reason and excluded from the aggregate, never dropped (RYA-711).
* An Engine-A product for lines MPIA does not serve: reported as reduced coverage, with
  the count, rather than silently averaging fewer lines.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.band_policy import resolve as resolve_band  # noqa: E402
from pipeline.band_products import build_product, LineMeasurement, products_frame  # noqa: E402
from pipeline.error_budget import build as build_budget  # noqa: E402

EW_DIR = ROOT / "data" / "measured" / "band_ew"
OUT = ROOT / "data" / "results" / "band_products"

# Measured optical residual of the profile fitter against the banked HARPS pool (RYA-713).
# Carried as the harness systematic in every frontier band rather than assumed zero.
PROFILE_FIT_RESIDUAL_DEX = 0.0129


def mpia_delta(element: str, ion: str, waves: np.ndarray,
               cache: str | None = None) -> dict[float, float]:
    """Per-line Engine-A departure corrections. Missing = not served.

    The MPIA scraper needs `bs4`, which is absent from the Sirius interpreter that has
    iSpec. Rather than installing into a pinned environment mid-run, the query can be run
    where bs4 lives and handed over as a cache -- the numbers are identical, and the file
    records which lines MPIA actually served.
    """
    if cache:
        d = json.loads(Path(cache).read_text())
        return {float(k): float(v) for k, v in d.items()}
    sys.path.insert(0, str(ROOT))
    from scripts.build_nlte_grids_mpia import _submit_batch
    from config.constants import STAR_PARAMS
    p = STAR_PARAMS["solar"]
    node = [dict(name="sun", teff_K=int(round(float(p["teff"]))),
                 logg=float(p["logg"]), feh=0.0)]
    code = {"I": "26.01", "II": "26.02"}.get(ion) if element == "Fe" else None
    if code is None:
        return {}
    out: dict[float, float] = {}
    for i in range(0, len(waves), 40):          # MPIA dislikes very long batches
        chunk = [float(w) for w in waves[i:i + 40]]
        try:
            r = _submit_batch(chunk, code, node).get("sun", {})
        except Exception as e:
            print(f"    MPIA batch {i//40} failed: {type(e).__name__}: {str(e)[:80]}")
            continue
        for w in chunk:
            if not r:
                continue
            k = min(r, key=lambda x: abs(x - w))
            v = r[k]
            if abs(k - w) < 0.6 and v is not None and np.isfinite(v) and v != 0.0:
                out[w] = float(v)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--element", required=True)
    ap.add_argument("--ion", default="I")
    ap.add_argument("--lo", type=float, required=True)
    ap.add_argument("--hi", type=float, required=True)
    ap.add_argument("--instrument", default="kpno_solar_atlas")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--mpia-cache", help="pre-computed per-line delta_nlte JSON")
    ap.add_argument("--engine-b-deck", choices=["ts-lte", "gerber-nlte"], default="ts-lte",
                    help="which synthesis is Engine B. 'ts-lte' is the production "
                         "Turbospectrum LTE flux-fit; 'gerber-nlte' is the TS-native "
                         "departure deck and is not provisioned for Fe until RYA-785.")
    ap.add_argument("--skip-engine-b", action="store_true",
                    help="derive 1D-LTE and ENGINE-A only. Engine B refits the spectrum "
                         "per line and is much slower than the EW inversion.")
    a = ap.parse_args()

    stem = f"{a.element}{a.ion}_{int(a.lo)}_{int(a.hi)}_{a.instrument}_PROFILEFIT"
    src = EW_DIR / f"{stem}_ew.csv"
    if not src.exists():
        raise SystemExit(f"no measured EWs at {src}. Run measure_band_profilefit.py first.")
    ew = pd.read_csv(src)
    pol = resolve_band(0.5 * (a.lo + a.hi))
    ok = ew[ew.in_aggregate].copy()
    print(f"{a.element} {a.ion}  {a.lo:.0f}-{a.hi:.0f} A  band={pol.name}  "
          f"{len(ok)} in-aggregate of {len(ew)} measured")
    if not len(ok):
        raise SystemExit("no lines survived measurement -- nothing to derive")

    # ── 1D-LTE ────────────────────────────────────────────────────────────────
    print("\n[1] 1D-LTE via the project's Turbospectrum curve-of-growth...")
    from scripts.control_synthesis_handler import build_context
    ctx = build_context(a.element, a.ion, 500000.0)
    from pipeline.abundances_derive import _bisect_synth_abundance
    acc = pd.read_csv(ROOT / "data" / "audit" / "line_accounting" / "per_line.csv")

    rows: list[LineMeasurement] = []
    for _, r in ok.iterrows():
        c = float(r.wavelength_air_A)
        w_nm = np.linspace(c - 0.6, c + 0.6, 300) / 10.0
        try:
            A, conv, _ = _bisect_synth_abundance(
                w_nm, float(r.ew_mA), ctx["atmosphere"], ctx["teff"], ctx["logg"],
                ctx["feh"], ctx["vturb"], ctx["linelist"], ctx["isotopes"],
                ctx["solar_abund"], a.element, ctx["atom_code"])
        except Exception as e:
            A, conv = float("nan"), False
        lm = LineMeasurement(element=a.element, ion=a.ion, wavelength_air_A=c,
                             instrument=a.instrument, ew_mA=float(r.ew_mA),
                             ew_method=str(r.ew_method), abundance=(A if conv else None),
                             treatment="1D-LTE")
        if not conv or not np.isfinite(A):
            lm.in_aggregate = False
            lm.excluded_reason = ("COG-INVERSION: bisection did not converge, so this EW "
                                  "does not map to an abundance here")
        rows.append(lm)

    p_lte = build_product(a.element, a.ion, a.instrument, pol.name, "1D-LTE", rows,
                          provenance=f"EWs from {src.name}; COG via _bisect_synth_abundance")
    print(f"    1D-LTE: A={p_lte.value}  n={p_lte.n_lines}  excluded={p_lte.n_excluded}")

    # ── ENGINE-A ──────────────────────────────────────────────────────────────
    print("\n[2] ENGINE-A — live MPIA per-line departure corrections...")
    used = [l for l in rows if l.in_aggregate and l.abundance is not None]
    deltas = mpia_delta(a.element, a.ion,
                        np.array([l.wavelength_air_A for l in used]),
                        cache=a.mpia_cache)
    rows_a: list[LineMeasurement] = []
    for l in used:
        # Tolerance match -- the cache keys are rounded, and two catalogues quoting the
        # same transition differ by up to ~0.03 A (RYA-704).
        d = None
        if deltas:
            k = min(deltas, key=lambda x: abs(x - l.wavelength_air_A))
            if abs(k - l.wavelength_air_A) < 0.06:
                d = deltas[k]
        la = LineMeasurement(element=l.element, ion=l.ion,
                             wavelength_air_A=l.wavelength_air_A,
                             instrument=l.instrument, ew_mA=l.ew_mA,
                             ew_method=l.ew_method, treatment="ENGINE-A",
                             abundance=(l.abundance + d) if d is not None else None)
        if d is None:
            la.in_aggregate = False
            la.excluded_reason = ("ENGINE-A-NOT-SERVED: MPIA returns no usable delta_nlte "
                                 "for this line (absent, nan, or a placeholder zero). "
                                 "Reduced coverage, not a failed correction.")
        rows_a.append(la)
    p_a = build_product(a.element, a.ion, a.instrument, pol.name, "ENGINE-A", rows_a,
                        provenance="Bergemann MPIA per-line delta_nlte, live query, solar node")
    print(f"    ENGINE-A: A={p_a.value}  n={p_a.n_lines}  not-served={p_a.n_excluded}")

    # ── ENGINE-B ──────────────────────────────────────────────────────────────
    #
    # RYA-784. This block did not exist: the driver reserved Engine B and produced no
    # value, so RYA-783's Engine-B column could not resolve. The flux-fit engine itself
    # has been Done since RYA-287/338 and its adapter was fixed by RYA-770 (24/30 lines,
    # -0.026 dex against the banked optical answer); nothing had connected the two.
    #
    # It is a DIFFERENT MEASUREMENT, not a correction applied to the previous one:
    # 1D-LTE inverts an EW and ENGINE-A adds a departure term to that same inversion,
    # while this fits the observed flux. So it re-derives from the spectrum and shares
    # only the line set — which is exactly what makes the per-engine spread meaningful
    # rather than tautological (RYA-712).
    #
    # RYA-342: the gate here is FIT QUALITY, not the EW saturation ceiling. That ceiling
    # is an EW-path concept and the handler must never inherit it — flux-space synthesis
    # exists precisely to measure the lines EW filtering kills.
    # RYA-288: broadening comes from `build_context` (per-star, refuses a vmac='fit' rail),
    # so no solar default leaks in.
    p_b = None
    if a.skip_engine_b:
        print("\n[3] ENGINE-B — skipped (--skip-engine-b)")
    elif a.engine_b_deck == "gerber-nlte":
        # RYA-785. The reason this refuses has CHANGED, and the old message is now false:
        # the Fe deck IS provisioned and, as of the second pass, VALIDATED — it reproduces
        # Gerber+2023's own published solar anchor (+0.06) to 0.0011 dex on 8 isolated
        # lines. So "not provisioned" would be a wrong statement, and a wrong reason is
        # worth as much as a wrong number.
        #
        # What is still missing is PLUMBING, named precisely so nobody re-derives it:
        #   * iSpec's Turbospectrum wrapper DOES accept NLTE — `generate_spectrum(...,
        #     nlte_departure_coefficients=...)` writes the departure + nlteinfo files;
        #   * but `_fit_synth_flux`, which this driver's handler calls, has no NLTE
        #     parameter at all, so there is nowhere to hand them in;
        #   * and iSpec's own interpolator reads `input/dep-grid/{El}_nlte_grid_data.h5`,
        #     which is EMPTY on Sirius — our deck is TS-native (atom.* + auxData + .bin),
        #     a different format, currently consumed only by the gate's direct
        #     interpol_modeles_nlte + bsyn calls.
        #
        # ⚠️ And iSpec fails SOFT here: an element whose linelist rows carry no NLTE label
        # lands in `nlte_ignored` and the synthesis proceeds in LTE without raising. That is
        # the RYA-764 trap (2,644 in-band Fe I lines with nlte_label_up='none' running
        # silently in LTE). So the wiring owes a hard check on `nlte_available`, not just a
        # successful call — which is exactly why this stays a refusal rather than a
        # best-effort attempt. Returning the LTE synthesis under an NLTE label is still the
        # single worst outcome available here.
        raise SystemExit(
            "ENGINE-B deck 'gerber-nlte' is VALIDATED but NOT WIRED. The Gerber deck "
            "passed its RYA-785 gate (median +0.0589 vs the published Gerber+2023 solar "
            "anchor +0.06, n=8 isolated lines), so the physics is cleared — but "
            "`_fit_synth_flux` has no NLTE parameter and iSpec's dep-grid store is empty, "
            "so there is no path from this driver to a real NLTE synthesis. Refusing: "
            "returning the LTE synthesis under an NLTE label is the single worst outcome "
            "available here. Use --engine-b-deck ts-lte. See RYA-785 for the wiring ticket.")
    else:
        print(f"\n[3] ENGINE-B — Turbospectrum LTE flux-fit "
              f"(deck={a.engine_b_deck}; NOT the Gerber NLTE deck)...")
        from pipeline.measure import resolve_handler
        from scripts.measure_band_ew import kp_segments, load_kp_window

        # TELLURIC — RYA-786. An earlier version of this block DECLARED
        # `telluric_corrected: True` from the instrument catalog to get past the handler's
        # gate. That was wrong in kind and is removed: `telluric_required = no` for a
        # reference atlas means the tellurics are handled by per-line clean-line SELECTION,
        # NOT that the atlas has been telluric-divided. The KPNO atlas HAS tellurics in it,
        # so declaring it corrected asserted something that never happened.
        #
        # The determination now lives in `pipeline.telluric_policy` and the handler
        # consumes it: the instrument decides whether a CORRECTION STAGE is required, and
        # the enumerated O2/H2O bands quarantine individual lines inside them. Passing the
        # instrument through the context is all this driver has to do.
        handler = resolve_handler(3400.0)      # the synthesis handler
        handler.prepare(pol, {**ctx, "instrument": a.instrument})
        segs = kp_segments()
        rows_b: list[LineMeasurement] = []
        for _, r in ok.iterrows():
            c = float(r.wavelength_air_A)
            try:
                w_obs, f_obs, _ = load_kp_window(segs, c, pad=1.4)
            except Exception as e:
                lb = LineMeasurement(element=a.element, ion=a.ion, wavelength_air_A=c,
                                     instrument=a.instrument, ew_mA=float("nan"),
                                     ew_method="synthesis flux-fit", treatment="ENGINE-B")
                lb.in_aggregate = False
                lb.excluded_reason = f"WINDOW-LOAD: {type(e).__name__}: {str(e)[:60]}"
                rows_b.append(lb)
                continue
            lb = handler.measure_line(
                w_obs, f_obs, element=a.element, ion=a.ion, wavelength_A=c,
                instrument=a.instrument, policy=pol,
                # Kitt Peak is atlas residual flux; HARPS arrives normalised by our own
                # pipeline. Both are already on a continuum and neither wants a second.
                pre_normalised=True,
                context={**ctx, "ew_hint_mA": float(r.ew_mA)})
            lb.treatment = "ENGINE-B"
            rows_b.append(lb)

        # RYA-768: deterministic row order, so the artifact byte-diffs clean.
        rows_b.sort(key=lambda l: (l.wavelength_air_A, l.element, l.ion))
        p_b = build_product(
            a.element, a.ion, a.instrument, pol.name, "ENGINE-B", rows_b,
            provenance=(f"{handler.name} flux-fit against the {a.instrument} spectrum; "
                        f"line set + window hint from {src.name}; deck={a.engine_b_deck} "
                        f"(Turbospectrum LTE, NOT the Gerber TS-native NLTE deck)"))
        print(f"    ENGINE-B: A={p_b.value}  n={p_b.n_lines}  excluded={p_b.n_excluded}")
        if p_b.n_excluded:
            why = pd.Series([l.excluded_reason.split(":")[0]
                             for l in rows_b if not l.in_aggregate]).value_counts()
            for k, v in why.items():
                print(f"      {k:<30} {v}")

    # ── budgets, per product ──────────────────────────────────────────────────
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    summary, budgets = [], {}
    for prod in (p_lte, p_a, p_b):
        if prod is None or prod.value is None:
            continue
        # The harness residual is a property of the HANDLER that produced the number.
        # Engine B never touches the profile fitter, so charging it the profile fitter's
        # measured residual would be attributing someone else's systematic to it.
        is_b = prod.treatment == "ENGINE-B"
        b = build_budget(a.element, 0.5 * (a.lo + a.hi), prod.n_lines,
                         scatter_dex=(prod.sigma or 0.0), gf_graded=False,
                         harness_residual_dex=(0.0 if is_b else PROFILE_FIT_RESIDUAL_DEX),
                         handler=("SynthesisHandler" if is_b else "ProfileFitHandler"))
        stat, syst = b.total()
        summary.append(dict(element=a.element, ion=a.ion, band=pol.name,
                            instrument=a.instrument, treatment=prod.treatment,
                            A=round(prod.value, 3), n_lines=prod.n_lines,
                            n_excluded=prod.n_excluded,
                            stat_dex=round(stat, 4), syst_dex=round(syst, 4),
                            dominant=b.dominant().name if b.dominant() else ""))
        budgets[prod.treatment] = b.describe()
        prod.to_frame().to_csv(out / f"{stem}_{prod.treatment}_lines.csv", index=False)

    df = pd.DataFrame(summary)
    df.to_csv(out / f"{stem}_products.csv", index=False)
    (out / f"{stem}_budgets.txt").write_text("\n\n".join(budgets.values()))

    print("\n=== PRODUCTS (separate, never combined — RYA-712) ===")
    print(df.to_string(index=False))
    print()
    for t, d in budgets.items():
        print(d); print()
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()

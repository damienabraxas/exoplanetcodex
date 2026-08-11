#!/usr/bin/env python3
"""RYA-716 — Al solar products: 1D-LTE and Engine A, kept separate (RYA-712).

    ISPEC_DIR=/srv/codex/engines/ispec_src CODEX_INTAKE_TRACE=1 \
      /mnt/codex-data/venv312/bin/python scripts/rya716_al_products.py

WHY THIS DRIVER AND NOT `derive_band_products.py`
--------------------------------------------------
It was written when `derive_band_products.py` was NOT on main — it rode the unmerged
RYA-759 branch, and using it would have coupled Al to a ticket being worked in parallel.
**RYA-774 has since landed the shared band harness on main (PR #219), so that is no
longer true** and the constraint that produced this driver is gone.

It is kept for what it adds over that driver, not to avoid it. The FINISH-Al brief pins
Al's number to *"the validated core (`_fit_synth_flux` / bisection), NOT the
SynthesisHandler adapter (that is RYA-770's problem)"* — which is what this calls
(`abundances_derive._bisect_synth_abundance`) — and it carries the two behaviours the Al
result turned on: the per-line Engine-A disposition via
`nlte_corrections._mpia_element_delta`, so no line is ever carried as silent LTE, and the
gf-provenance split that isolates the 7836.134 outlier. It invents no physics.

Folding those two behaviours back into `derive_band_products.py` and retiring this driver
is the right follow-up now that the harness is shared. It is deliberately NOT done at
merge time.

WHAT IT EMITS — two products, never combined (RYA-712)
------------------------------------------------------
  * **Al 1D-LTE**   — measured EW inverted through the Turbospectrum curve of growth.
  * **Al ENGINE-A** — 1D-LTE + the registered Amarsi-2020/PySME per-line departure
                      correction. A line the grid does not serve is NOT quietly carried
                      at its LTE value into this product; it is dispositioned and
                      excluded, and the count is reported.
  * **Al ENGINE-B** — DEFERRED to RYA-710 (in progress). Not derived, not written.

Both are read against `SOLAR_ASPLUND2021['Al']`, pulled from config.constants.

CAVEATS CARRIED, NOT HIDDEN
---------------------------
* gf is UNGRADED for these lines -> the systematic is large; the value is not disqualified.
  The pool is NOT homogeneous in gf source, and that turns out to matter: see below.
* RYA-771: the pinned engine quantises the trial abundance to 0.01 dex, so NO line-to-line
  scatter below 0.01 dex is quotable off this run. Reported sigmas are floored at 0.01.

THE gf-SOURCE SPLIT (RYA-716 step 3 finding)
--------------------------------------------
The four in-aggregate lines do not share a gf provenance. Three carry the 1995 J.Phys.B
experimental oscillator strengths; 7836.134 alone carries **K75** (Kurucz 1975
semi-empirical), the ungraded class RYA-161 assigns a 0.1-0.3 dex systematic. That is
exactly where the one outlier sits (+0.12 dex). So the primary aggregate is taken over
the gf-homogeneous experimental subset, and the K75 line is reported separately rather
than averaged in. The split is READ FROM THE LINE LIST, not hardcoded to a wavelength.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

ISPEC_SRC = os.environ.get("ISPEC_DIR", "/srv/codex/engines/ispec_src")
if ISPEC_SRC not in sys.path:
    sys.path.insert(0, ISPEC_SRC)

from config.constants import STAR_PARAMS, SOLAR_ASPLUND2021, NLTE_CORRECTION_ELEMENTS  # noqa: E402
from pipeline import intake_debug as dbg                                               # noqa: E402

EW_DIR = ROOT / "data" / "measured" / "band_ew"
OUT_DIR = ROOT / "data" / "results" / "band_products"
BANDS = [("3800", "6910", "VIS"), ("6910", "9199", "red-optical")]

# RYA-771: the pinned iSpec writes trial abundances with %.2f (twice), so EW(A) is a
# staircase with 0.01 dex treads. Nothing below this is resolvable on the pinned engine.
QUANTISER_FLOOR_DEX = 0.01

# gf reference codes that are semi-empirical/theoretical rather than measured. RYA-161
# puts this class at 0.1-0.3 dex. A line carrying one of these is not averaged into the
# primary product; it is reported on its own so the gf term stays visible.
SEMI_EMPIRICAL_GF_REFS = ("K75", "K88", "KURUCZ")


def gf_provenance(linelist, element: str, wave_A: float, tol: float = 0.05) -> tuple:
    """(reference_code, loggf, is_semi_empirical) for `element` nearest `wave_A`.

    Read from the synthesis line list itself so the disposition is data-driven — no
    wavelength is hardcoded. Where a line has several J-components at one wavelength
    (7836.134 does), the STRONGEST is the one that sets the gf provenance.
    """
    w = np.asarray(linelist["wave_A"], dtype=float)
    el = np.asarray([str(x) for x in linelist["element"]])
    gf = np.asarray(linelist["loggf"], dtype=float)
    idx = [i for i in np.where(np.abs(w - wave_A) < tol)[0]
           if el[i].strip().upper().startswith(element.upper())]
    if not idx:
        return ("NOT-IN-LINELIST", np.nan, False)
    i = max(idx, key=lambda k: gf[k])
    ref = str(linelist["reference_code"][i]).strip()
    semi = any(ref.upper().startswith(p) for p in SEMI_EMPIRICAL_GF_REFS)
    return (ref, float(gf[i]), semi)


def _guard_numpy() -> None:
    """RYA-682: numpy >= 2.3 breaks ispec/abundances.py:132 SILENTLY (0-row artifact,
    exit 0). Refuse to run rather than emit a confident empty product."""
    import numpy
    if tuple(int(x) for x in numpy.__version__.split(".")[:2]) >= (2, 3):
        raise SystemExit(f"numpy {numpy.__version__} is at/above the RYA-682 ceiling 2.3")


def load_measured() -> pd.DataFrame:
    frames = []
    for lo, hi, band in BANDS:
        p = EW_DIR / f"AlI_{lo}_{hi}_kpno_solar_atlas_PROFILEFIT_ew.csv"
        if not p.exists():
            raise SystemExit(f"missing measured EWs: {p}")
        d = pd.read_csv(p)
        d["band"] = band
        frames.append(d)
        dbg.trace_asset(f"measured_ew[{band}]", True, path=str(p),
                        detail=f"{int(d.in_aggregate.sum())} in-aggregate of {len(d)}")
    return pd.concat(frames, ignore_index=True)


def build_context(element: str):
    """Synth context from main's own loaders — nothing invented here."""
    import ispec
    from pipeline.abundances_derive import (_load_atmosphere, _load_synth_resources,
                                            _ISPEC_SOLAR_ABUND_FILE)
    from pipeline.cno_synthesis import _atom_codes

    p = STAR_PARAMS["solar"]
    teff, logg = float(p["teff"]), float(p["logg"])
    if "xi" not in p:
        raise SystemExit("no microturbulence (`xi`) for the Sun — refusing to default it")
    vturb = float(p["xi"])
    feh = float(p.get("feh", p.get("feh_ref", 0.0)))
    linelist, isotopes, chem = _load_synth_resources()
    solar_abund = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    atm = _load_atmosphere(teff, logg, feh, vturb)
    codes = _atom_codes([element], chem, solar_abund)
    dbg.trace_asset("synth_linelist", True, detail=f"{len(linelist)} transitions")
    dbg.trace_asset("atmosphere", True,
                    detail=f"MARCS teff={teff} logg={logg} feh={feh} vturb={vturb}")
    return dict(atmosphere=atm, teff=teff, logg=logg, feh=feh, vturb=vturb,
                linelist=linelist, isotopes=isotopes, solar_abund=solar_abund,
                atom_code=int(codes[element]))


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--aggregate-only", action="store_true",
                    help="re-aggregate from the saved per-line CSV instead of "
                         "re-running the (10-minute) synthesis. The per-line "
                         "inversions are deterministic, so this changes no number — "
                         "it only re-derives the products from them.")
    ap.add_argument("--no-trace", action="store_true",
                    help="do not start an intake trace. Used to CHECK the tracer's "
                         "inertness invariant: this run's numbers must be identical to "
                         "the traced run's. Note CODEX_INTAKE_TRACE=1 is a SEPARATE "
                         "activation path — start_trace() here is unconditional, so "
                         "without this flag a run is traced whatever the env says.")
    args = ap.parse_args()
    element, ion = "Al", "I"

    if args.aggregate_only:
        per_line_p = OUT_DIR / "AlI_kpno_solar_atlas_per_line.csv"
        if not per_line_p.exists():
            raise SystemExit(f"--aggregate-only needs {per_line_p}; run the full pass first")
        df = pd.read_csv(per_line_p)
        print(f"\n[aggregate-only] re-aggregating {len(df)} saved per-line inversions "
              f"from {per_line_p.name} — no synthesis re-run")
        report(df, element, NLTE_CORRECTION_ELEMENTS[element])
        return

    _guard_numpy()
    traced = not args.no_trace
    if traced:
        dbg.start_trace(element=element,
                        product="solar band products (1D-LTE + Engine A)", star="solar")
    else:
        print("\n[--no-trace] tracer NOT started — inertness check run")

    ew = load_measured()
    ok = ew[ew.in_aggregate].copy().sort_values("wavelength_air_A")
    print(f"\n{element} {ion}: {len(ok)} in-aggregate of {len(ew)} measured "
          f"(Kitt Peak solar atlas)")

    ctx = build_context(element)
    from pipeline.abundances_derive import _bisect_synth_abundance
    from pipeline.nlte_corrections import _mpia_element_delta, element_grid_in_bounds

    spec = NLTE_CORRECTION_ELEMENTS[element]
    dbg.trace_asset("nlte_grid", True, path=spec["grid"], detail=spec["ref"][:70])
    in_bounds = element_grid_in_bounds(element, ctx["teff"], ctx["logg"], ctx["feh"])
    dbg.trace_check("nlte_grid_in_bounds", in_bounds,
                    detail=f"solar ({ctx['teff']:.0f}, {ctx['logg']:.2f}, {ctx['feh']:+.2f})")

    kw = dict(atmosphere=ctx["atmosphere"], teff=ctx["teff"], logg=ctx["logg"],
              feh=ctx["feh"], vturb=ctx["vturb"], linelist=ctx["linelist"],
              isotopes=ctx["isotopes"], solar_abund=ctx["solar_abund"],
              element=element, atom_code=ctx["atom_code"])

    rows = []
    print(f"\n{'line':>11} {'band':>12} {'EW_obs':>8} {'A_LTE':>9} "
          f"{'delta_A':>9} {'A_EngA':>9}  disposition")
    for _, r in ok.iterrows():
        c = float(r.wavelength_air_A)
        w_nm = np.linspace(c - 0.6, c + 0.6, 300) / 10.0
        A, conv, n_it = _bisect_synth_abundance(w_nm, float(r.ew_mA), **kw)
        if not conv or not np.isfinite(A):
            dbg.trace_fallback("cog_inversion_failed",
                               f"{element} {c:.3f}: bisection did not converge -> no A",
                               severity="ERROR", species=element, wavelength=c)
            rows.append(dict(wave=c, band=r.band, ew=float(r.ew_mA), a_lte=np.nan,
                             delta=np.nan, a_a=np.nan,
                             disp="COG-INVERSION: did not converge"))
            print(f"{c:11.4f} {r.band:>12} {r.ew_mA:8.2f} {'--':>9} {'--':>9} {'--':>9}"
                  f"  COG-INVERSION FAILED")
            continue

        # Engine A — the registered per-line departure correction. A NaN here is a line
        # the grid does not serve; it is dispositioned, never carried silently as LTE.
        d = _mpia_element_delta(element, c, ctx["teff"], ctx["logg"], ctx["feh"])
        if np.isfinite(d):
            a_a, disp = A + float(d), "ENGINE-A: grid-served"
        else:
            a_a, disp = (np.nan,
                         "ENGINE-A UNCOVERED: no Amarsi-2020 node within 0.15 A "
                         "-> LTE-only, EXCLUDED from the Engine-A product")
            dbg.trace_decision("engine_a_line_disposition", "EXCLUDED",
                               detail=f"{element} {c:.3f}: uncovered by "
                                      f"{spec['grid']} -> not carried as silent LTE",
                               species=element, wavelength=c)
        gf_ref, loggf, semi = gf_provenance(ctx["linelist"], element, c)
        if semi:
            disp += f" | gf={gf_ref} SEMI-EMPIRICAL (RYA-161 0.1-0.3 dex)"
            dbg.trace_decision("gf_provenance_disposition", "EXCLUDED-FROM-PRIMARY",
                               detail=f"{element} {c:.3f}: loggf from {gf_ref} "
                                      f"(semi-empirical) — not averaged into the "
                                      f"gf-homogeneous product",
                               species=element, wavelength=c, gf_ref=gf_ref)
        rows.append(dict(wave=c, band=r.band, ew=float(r.ew_mA), a_lte=float(A),
                         delta=(float(d) if np.isfinite(d) else np.nan),
                         a_a=a_a, gf_ref=gf_ref, loggf=loggf, gf_semi_empirical=semi,
                         disp=disp))
        print(f"{c:11.4f} {r.band:>12} {r.ew_mA:8.2f} {A:9.4f} "
              f"{(f'{d:+9.4f}' if np.isfinite(d) else '       --')} "
              f"{(f'{a_a:9.4f}' if np.isfinite(a_a) else '       --')}  "
              f"gf={gf_ref:<11s}{'SEMI-EMP' if semi else ''}")

    df = pd.DataFrame(rows)
    report(df, element, spec)
    if traced:
        dbg.end_trace()


def report(df: pd.DataFrame, element: str, spec: dict) -> None:
    """Aggregate the per-line inversions into the products and print the comparison.

    Shared by the full pass and by --aggregate-only so there is exactly ONE aggregation
    rule; a second copy is how the two paths would silently drift apart.
    """
    ref = float(SOLAR_ASPLUND2021[element])

    def product(name: str, vals: np.ndarray, note: str) -> dict:
        vals = np.asarray([v for v in vals if np.isfinite(v)], dtype=float)
        n = len(vals)
        if n == 0:
            return dict(product=name, A=np.nan, sigma=np.nan, n_lines=0, note=note)
        sd = float(np.std(vals, ddof=1)) if n > 1 else np.nan
        # RYA-771: the engine cannot resolve below one 0.01 dex tread, so the RAW
        # line-to-line sd is not a measurement of anything finer. Three lines landing on
        # the same tread give sd = 0 exactly — quoting sigma = 0 would claim infinite
        # precision from a quantiser artifact. Floor the scatter at the tread, and derive
        # the quoted statistical error from the FLOORED scatter, never from the raw 0.
        sd_q = (max(sd, QUANTISER_FLOOR_DEX) if n > 1 else QUANTISER_FLOOR_DEX)
        sem_q = sd_q / np.sqrt(n)
        return dict(product=name, A=float(np.mean(vals)), sigma_stat=float(sem_q),
                    sd_lines=sd, sd_quoted=float(sd_q), n_lines=n,
                    sigma_is_floor=bool(not np.isfinite(sd) or sd <= QUANTISER_FLOOR_DEX),
                    note=note)

    # PRIMARY = the gf-homogeneous (experimental-gf) subset. The semi-empirical line is
    # reported beside it, never averaged in — the gf term is the dominant systematic and
    # mixing sources inside one mean hides it.
    graded = df[~df.gf_semi_empirical]
    semi = df[df.gf_semi_empirical]
    p_lte = product("Al 1D-LTE", graded.a_lte.values,
                    "measured EW -> Turbospectrum COG (_bisect_synth_abundance); "
                    "experimental-gf lines only")
    p_a = product("Al ENGINE-A", graded.a_a.values,
                  f"1D-LTE + {spec['grid']} per-line delta ({spec['flag']})")
    p_all = product("Al 1D-LTE (all)", df.a_lte.values,
                    "all in-aggregate lines incl. semi-empirical gf — DIAGNOSTIC ONLY")

    print("\n" + "=" * 78)
    print("PRODUCTS — separate, never combined (RYA-712)")
    print("=" * 78)
    for p in (p_lte, p_a):
        if p["n_lines"] == 0:
            print(f"  {p['product']:<14} NO VALUE — {p['note']}")
            continue
        sd = p["sd_quoted"]
        floored = " (FLOOR: RYA-771 quantiser, not a measured scatter)" \
            if p.get("sigma_is_floor") else ""
        print(f"  {p['product']:<14} A = {p['A']:.4f}   "
              f"sigma_stat = {p['sigma_stat']:.4f}{floored}   "
              f"sd(lines) = {'n/a' if not np.isfinite(sd) else f'{sd:.4f}'}   "
              f"n_lines = {p['n_lines']}")
        print(f"  {'':14} delta vs SOLAR_ASPLUND2021['{element}'] = {ref:.3f}  ->  "
              f"{p['A'] - ref:+.4f} dex")
        print(f"  {'':14} {p['note']}")
    print("  Al ENGINE-B    DEFERRED to RYA-710 (in progress) — not derived, "
          "no artifact written")

    if len(semi):
        print(f"\n  reported beside, NOT averaged in — semi-empirical gf "
              f"(RYA-161 0.1-0.3 dex systematic):")
        for _, s in semi.iterrows():
            print(f"    Al I {s.wave:.4f}  A_LTE = {s.a_lte:.4f}  "
                  f"({s.a_lte - p_lte['A']:+.4f} vs the experimental-gf mean)  "
                  f"gf={s.gf_ref} loggf={s.loggf:+.3f}")
        print(f"    diagnostic all-line mean (do NOT quote): "
              f"A = {p_all['A']:.4f}, n = {p_all['n_lines']}")

    if np.isfinite(p_lte["A"]) and np.isfinite(p_a["A"]):
        print(f"\n  LTE -> Engine-A gap (the NLTE correction, visible): "
              f"{p_a['A'] - p_lte['A']:+.4f} dex")
        served = df[df.delta.notna()]
        if len(served):
            print(f"    per-line deltas applied: "
                  f"{', '.join(f'{w:.3f}: {d:+.4f}' for w, d in zip(served.wave, served.delta))}")

    n_unc = int(df.a_lte.notna().sum() - df.a_a.notna().sum())
    print(f"\n  Engine-A coverage: {int(df.a_a.notna().sum())} served / "
          f"{int(df.a_lte.notna().sum())} with an LTE value  "
          f"({n_unc} uncovered, dispositioned, NOT carried as silent LTE)")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    per_line = OUT_DIR / "AlI_kpno_solar_atlas_per_line.csv"
    df.to_csv(per_line, index=False)
    prod = OUT_DIR / "AlI_kpno_solar_atlas_products.csv"
    pd.DataFrame([p_lte, p_a, p_all]).to_csv(prod, index=False)
    print(f"\n  wrote {per_line.relative_to(ROOT)}")
    print(f"  wrote {prod.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

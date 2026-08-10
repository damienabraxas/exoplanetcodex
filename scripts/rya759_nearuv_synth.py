#!/usr/bin/env python3
"""Near-UV Fe LTE synthesis on the production atmosphere path — RYA-759 Move 2.

    # MAC (spectra live here, RYA-567): cut the observed windows and ship them
    python3 scripts/rya759_nearuv_synth.py --step windows --windows-out nearuv_win.npz

    # SIRIUS (compute lives here): everything else
    venv312/bin/python scripts/rya759_nearuv_synth.py --step all --windows nearuv_win.npz

STEPS
-----
  linelist  build the iSpec near-UV atomic line list from our VALD holdings
  reach     does iSpec + Turbospectrum synthesise at 3000-3800 A at all?  (the one
            genuine unknown — the 4200 A wall was always the GES list, never iSpec)
  nearuv    the FIRST REAL near-UV run: non-zero spectrum, Fe features at the right
            wavelengths, a sane count of lines actually entering the synthesis
  physics   synthetic vs OBSERVED (Kitt Peak) in a near-UV window — position + depth
  loudfail  point it at an empty line list and a bogus model; both must RAISE

The optical control (ATLAS9.Castelli, solar params, banked Fe) is NOT reimplemented
here — it is `scripts/control_synthesis_handler.py`, which already does it and is what
gates the SynthesisHandler. Run that; this script prints its recorded verdict.

WHAT THIS DOES NOT DO
---------------------
It does not quote a near-UV Fe abundance. Move 2's job is the capability and its
controls. And when the abundance does come, expect it LOW: UV Fe I is heavily
over-ionized, so the missing (positive) Engine-A NLTE correction is large. A low
near-UV 1D-LTE value is physically expected — the RCA is "the NLTE correction has not
been applied yet", not "hunt a pipeline bug".
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

LO_A, HI_A = 3000.0, 3780.0
OUT = ROOT / "data" / "audit" / "nearuv_synth"

# Windows for the first real run and the physics check. Chosen as near-UV Fe I lines
# with laboratory-grade VALD gf sources and the largest neighbour gaps available in a
# band whose MEDIAN gap is 0.146 A -- "isolated" here is relative, and that is the
# point: the pseudo-continuum systematic (0.10 dex, error_budget.py) exists because
# nothing in this band is truly isolated.
PROBE_WINDOWS = [
    (3705.57, "Fe I 3705.566"),
    (3737.13, "Fe I 3737.131"),
    (3749.49, "Fe I 3749.485"),
]
PROBE_HALFWIDTH_A = 1.0


def _ispec_guard():
    """RYA-682: numpy >= 2.3 makes iSpec synthesis write a zero-row artifact and exit 0.
    A near-UV run that 'passes' on zero rows is the exact failure this ticket exists to
    stop, so this is checked first and loudly."""
    for p in ("/srv/codex/engines/ispec_src", "/mnt/codex-data/engines/ispec_src"):
        if Path(p).exists() and p not in sys.path:
            sys.path.insert(0, p)
    import numpy
    if tuple(int(x) for x in numpy.__version__.split(".")[:2]) >= (2, 3):
        raise SystemExit(
            f"numpy {numpy.__version__} >= 2.3 breaks ispec/abundances.py:132 and the "
            f"failure is SILENT (RYA-682). Use venv312 or venv_pysme on Sirius.")
    import ispec  # noqa: F401


def step_linelist(args) -> dict:
    from pipeline import nearuv_linelist as nl
    _ispec_guard()
    gf_csv = ROOT / "data" / "linelists" / f"vald_gf_sources_Fe_{int(LO_A)}_{int(HI_A)}.csv"
    out = Path(args.linelist) if args.linelist else None
    path, rep = nl.build(lo_A=LO_A, hi_A=HI_A, out_path=out,
                         gf_sources_csv=gf_csv if gf_csv.exists() else None)
    print(f"\n[linelist] {path}")
    print(f"  {rep['n_lines']:,} lines / {rep['n_species']} species  "
          f"{rep['lo_A']:.2f}-{rep['hi_A']:.2f} A")
    print(f"  parser: {rep['n_parsed']:,} parsed, {rep['n_failures']} malformed "
          f"(reported, never dropped silently)")
    print(f"  vdW=0 (VALD supplied none, TS default applies): {rep['n_vdw_zero']:,}")
    print(f"  vdW>0 (not ABO-packed, carried verbatim):       {rep['n_vdw_positive']:,}")
    print(f"  gamma_rad=0 (VALD supplied none):               {rep['n_rad_zero']:,}")
    print(f"  gf source tags matched from VALD refs:          {rep['gf_sources_matched']:,}")
    print("  top species: " + ", ".join(f"{e} {n}" for e, n in rep["top_species"]))
    return rep


def _context(args, *, linelist_file: str):
    """The SAME context builder the optical control uses (RYA-759 parameterised it)."""
    from scripts.control_synthesis_handler import build_context
    from pipeline.nearuv_synth import gf_provenance
    prov = gf_provenance(LO_A, HI_A)
    print(f"[gf] {prov['detail']}")
    ctx = build_context("Fe", "I", args.resolving_power,
                        linelist_file=linelist_file,
                        apply_canonical_gf=prov["apply_canonical_gf"])
    ctx["model_grid"] = "ATLAS9.Castelli"
    ctx["gf_provenance"] = prov["detail"]
    print(f"[atm] ATLAS9.Castelli at Teff={ctx['teff']:.0f} logg={ctx['logg']:.3f} "
          f"[Fe/H]={ctx['feh']:.2f} xi={ctx['vturb']:.2f}  "
          f"({len(np.atleast_1d(ctx['atmosphere']))} layers)")
    print(f"[brd] R={ctx['resolving_power']:.0f} vmac={ctx['macroturbulence']} "
          f"vsini={ctx['vsini']}   a_start(iSpec internal)={ctx['solar_A']:.3f}")
    return ctx


def step_reach(args) -> dict:
    """Probe the blue limit of iSpec's Turbospectrum wrapper, 3000 -> 3800 A."""
    from pipeline.nearuv_synth import synthesize_band, NearUVSynthesisError
    _ispec_guard()
    ctx = _context(args, linelist_file=args.linelist)
    probes = [3000.0, 3100.0, 3200.0, 3400.0, 3600.0, 3750.0, 3900.0, 4500.0]
    rows = []
    for lo in probes:
        hi = lo + 1.0
        try:
            r = synthesize_band(ctx, lo, hi, element="Fe",
                                trial_A=float(ctx["solar_A"]), step_A=0.02,
                                check_sensitivity=False)
            f = r["flux"]
            rows.append(dict(lo_A=lo, ok=True, n_pts=int(f.size),
                             flux_min=float(np.nanmin(f)),
                             flux_max=float(np.nanmax(f)),
                             n_lines=r["n_lines_in_band"]))
            print(f"  {lo:7.1f} A  OK   {f.size:5d} pts  flux {np.nanmin(f):.4f}-"
                  f"{np.nanmax(f):.4f}  ({r['n_lines_in_band']} lines in band)")
        except (NearUVSynthesisError, Exception) as e:          # noqa: B014
            rows.append(dict(lo_A=lo, ok=False, error=f"{type(e).__name__}: {e}"))
            print(f"  {lo:7.1f} A  FAIL {type(e).__name__}: {str(e)[:150]}")
    ok = [r["lo_A"] for r in rows if r.get("ok")]
    verdict = (f"clean from {min(ok):.0f} A" if ok else "NO probe succeeded")
    print(f"\n[reach] {verdict}")
    return {"probes": rows, "verdict": verdict}


def step_nearuv(args) -> dict:
    """The first real near-UV run. Every prior reading on this list is void."""
    from pipeline.nearuv_synth import synthesize_band
    _ispec_guard()
    ctx = _context(args, linelist_file=args.linelist)
    out = []
    for centre, label in PROBE_WINDOWS:
        lo, hi = centre - PROBE_HALFWIDTH_A, centre + PROBE_HALFWIDTH_A
        r = synthesize_band(ctx, lo, hi, element="Fe",
                            trial_A=float(ctx["solar_A"]), step_A=0.005)
        f, w = r["flux"], r["wave_A"]
        i = int(np.nanargmin(f))
        out.append(dict(label=label, centre_A=centre, n_pts=int(f.size),
                        n_lines_in_band=r["n_lines_in_band"],
                        sensitivity=r["sensitivity"],
                        flux_min=float(np.nanmin(f)), flux_median=float(np.nanmedian(f)),
                        deepest_A=float(w[i]), offset_mA=float((w[i] - centre) * 1000)))
        print(f"  {label:>16s}  {f.size} pts  {r['n_lines_in_band']} lines in +/-1 A  "
              f"flux min {np.nanmin(f):.4f} med {np.nanmedian(f):.4f}  "
              f"deepest at {w[i]:.3f} A ({(w[i]-centre)*1000:+.0f} mA)  "
              f"sens {r['sensitivity']:.4f}")
    return {"windows": out}


def step_windows(args) -> dict:
    """MAC SIDE: cut observed Kitt Peak windows and ship them to Sirius."""
    from scripts.measure_band_ew import kp_segments, load_kp_window
    segs = kp_segments()
    store, kept = {}, []
    for centre, label in PROBE_WINDOWS:
        w, f, _ = load_kp_window(segs, centre, pad=PROBE_HALFWIDTH_A + 0.4)
        store[f"w_{centre:.4f}"] = np.asarray(w, dtype="float32")
        store[f"f_{centre:.4f}"] = np.asarray(f, dtype="float32")
        kept.append(centre)
        print(f"  {label}: {len(w)} px  {w[0]:.2f}-{w[-1]:.2f} A")
    store["index"] = np.array(kept, dtype="float64")
    dest = Path(args.windows_out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dest, **store)
    print(f"\n[windows] wrote {dest} ({dest.stat().st_size/1e3:.0f} kB)")
    return {"path": str(dest), "centres": kept}


def step_physics(args) -> dict:
    """Synthetic vs OBSERVED in a near-UV window: is the line where it should be, and
    is it as deep as it should be? Validation against the sky, not against ourselves."""
    from pipeline.nearuv_synth import synthesize_band
    _ispec_guard()
    if not args.windows:
        raise SystemExit("--windows NPZ is required (make it on the Mac: --step windows)")
    z = np.load(args.windows)
    ctx = _context(args, linelist_file=args.linelist)
    rows = []
    for centre, label in PROBE_WINDOWS:
        kw, kf = f"w_{centre:.4f}", f"f_{centre:.4f}"
        if kw not in z:
            print(f"  {label}: no observed window in the NPZ — skipped")
            continue
        ow, of = z[kw].astype(float), z[kf].astype(float)
        lo, hi = centre - PROBE_HALFWIDTH_A, centre + PROBE_HALFWIDTH_A
        r = synthesize_band(ctx, lo, hi, element="Fe",
                            trial_A=float(ctx["solar_A"]), step_A=0.005)
        sw, sf = r["wave_A"], r["flux"]
        m = (ow >= lo) & (ow <= hi)
        ow, of = ow[m], of[m]
        if ow.size < 5:
            print(f"  {label}: {ow.size} observed px in band — skipped")
            continue
        # Position: nearest observed minimum to the synthetic minimum.
        si = int(np.nanargmin(sf))
        oi = int(np.nanargmin(of))
        # Depth: 1 - flux at the line core, both sides.
        sd, od = 1.0 - float(sf[si]), 1.0 - float(of[oi])
        rows.append(dict(label=label, synth_min_A=float(sw[si]), obs_min_A=float(ow[oi]),
                         dlambda_mA=float((sw[si] - ow[oi]) * 1000),
                         synth_depth=sd, obs_depth=od,
                         depth_ratio=float(sd / od) if od else float("nan"),
                         obs_median_flux=float(np.nanmedian(of)),
                         synth_median_flux=float(np.nanmedian(sf))))
        print(f"  {label:>16s}  synth min {sw[si]:.3f} A vs obs {ow[oi]:.3f} A "
              f"(d={((sw[si]-ow[oi])*1000):+.0f} mA)   depth synth {sd:.3f} / obs {od:.3f} "
              f"= {sd/od if od else float('nan'):.2f}   "
              f"median flux synth {np.nanmedian(sf):.3f} / obs {np.nanmedian(of):.3f}")
    return {"windows": rows}


def step_loudfail(args) -> dict:
    """Proof that the three guards RAISE. A guard nobody has seen fire is a guess."""
    from pipeline import nearuv_synth as ns
    _ispec_guard()
    results = []

    def _expect(name, fn):
        try:
            fn()
        except ns.NearUVSynthesisError as e:
            print(f"  {name:<24s} RAISED  {str(e)[:110]}")
            results.append(dict(case=name, raised=True, error=str(e)[:400]))
            return
        except Exception as e:                                   # noqa: BLE001
            print(f"  {name:<24s} raised {type(e).__name__} (not the guard): "
                  f"{str(e)[:100]}")
            results.append(dict(case=name, raised=True, unexpected=type(e).__name__,
                                error=str(e)[:400]))
            return
        print(f"  {name:<24s} *** DID NOT RAISE — the guard is not wired ***")
        results.append(dict(case=name, raised=False))

    empty = np.zeros(0, dtype=[("wave_A", "<f8"), ("element", "|U4")])
    _expect("empty linelist", lambda: ns.assert_linelist_covers(empty, LO_A, HI_A))
    _expect("linelist misses band", lambda: ns.assert_linelist_covers(
        np.array([(5000.0, "Fe 1")], dtype=[("wave_A", "<f8"), ("element", "|U4")]),
        LO_A, HI_A))
    _expect("bogus model (None)", lambda: ns.assert_atmosphere(
        None, teff=5772.0, logg=4.438, feh=0.0, vturb=1.0, model_grid="BOGUS"))
    _expect("all-zero flux", lambda: ns.assert_usable_flux(np.zeros(500),
                                                           where="loud-fail probe"))
    _expect("empty flux", lambda: ns.assert_usable_flux(np.zeros(0),
                                                        where="loud-fail probe"))
    _expect("flat (insensitive)", lambda: ns.assert_sensitive(
        np.ones(100), np.ones(100), where="loud-fail probe"))

    # A real out-of-grid interpolation, not a mock: ATLAS9 cannot serve Teff 2000 K.
    def _real_atm():
        from pipeline.abundances_derive import _load_atmosphere
        _load_atmosphere(2000.0, 4.438, 0.0, 1.0)
    _expect("real model out-of-grid", _real_atm)
    return {"cases": results}


def step_control(args) -> dict:
    """Report the recorded optical-control verdict — the gate on any near-UV number."""
    p = ROOT / "data" / "audit" / "synthesis_control" / "control_FeI.json"
    if not p.exists():
        print("  no recorded optical control — run scripts/control_synthesis_handler.py")
        return {"present": False}
    d = json.loads(p.read_text())
    print(f"  handler={d['handler']} n={d['n_lines']} dex_offset={d['dex_offset']:+.4f} "
          f"tol={d['tolerance_dex']}  PASSED={d['passed']}")
    print(f"  {d['evidence']}")
    return d


STEPS = {"linelist": step_linelist, "reach": step_reach, "nearuv": step_nearuv,
         "physics": step_physics, "loudfail": step_loudfail, "windows": step_windows,
         "control": step_control}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step", default="all",
                    choices=sorted(STEPS) + ["all"])
    ap.add_argument("--linelist",
                    default=str(ROOT / "data" / "linelists" /
                                f"ispec_nearuv_{int(LO_A)}_{int(HI_A)}" /
                                "atomic_lines.tsv"))
    ap.add_argument("--windows", help="observed Kitt Peak windows NPZ (from --step windows)")
    ap.add_argument("--windows-out", default=str(OUT / "nearuv_windows.npz"))
    ap.add_argument("--resolving-power", type=float, default=None,
                    help="default: the instrument catalog's value for kpno_solar_atlas")
    ap.add_argument("--no-trace", action="store_true",
                    help="skip the RYA-765 intake trace (it is numerically inert)")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()

    if a.resolving_power is None:
        import pandas as pd
        cat = pd.read_csv(ROOT / "data" / "catalog" / "instrument_catalog.csv")
        row = cat[cat.iloc[:, 0].astype(str) == "kpno_solar_atlas"]
        if row.empty:
            raise SystemExit("kpno_solar_atlas absent from data/catalog/instrument_catalog.csv")
        a.resolving_power = float(row.iloc[0]["resolving_power_max"])

    order = ["linelist", "reach", "nearuv", "physics", "loudfail", "control"] \
        if a.step == "all" else [a.step]

    trace = None
    if not a.no_trace and a.step != "windows":
        from pipeline import intake_debug
        trace = intake_debug.start_trace("Fe", "1D-LTE-nearUV", "solar",
                                         extra={"band_A": [LO_A, HI_A],
                                                "model_grid": "ATLAS9.Castelli"})

    results, failed = {}, []
    try:
        for name in order:
            print(f"\n{'='*72}\n== {name}\n{'='*72}")
            try:
                results[name] = STEPS[name](a)
            except Exception as e:                               # noqa: BLE001
                # Report and continue so one blocked step does not hide the rest; the
                # exit code still reflects the failure.
                print(f"  STEP FAILED: {type(e).__name__}: {e}")
                results[name] = {"error": f"{type(e).__name__}: {e}"}
                failed.append(name)
    finally:
        if trace is not None:
            from pipeline import intake_debug
            intake_debug.end_trace()

    if a.step != "windows":
        dest = Path(a.out); dest.mkdir(parents=True, exist_ok=True)
        f = dest / "rya759_nearuv_synth.json"

        def _j(o):
            if isinstance(o, (np.floating, np.integer)):
                return o.item()
            if isinstance(o, np.ndarray):
                return o.tolist()
            return str(o)
        f.write_text(json.dumps(results, indent=2, default=_j))
        print(f"\n[out] {f}")
    if failed:
        raise SystemExit(f"steps failed: {', '.join(failed)}")


if __name__ == "__main__":
    main()

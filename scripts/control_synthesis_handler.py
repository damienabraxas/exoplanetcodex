#!/usr/bin/env python3
"""Control the SynthesisHandler in the optical, BOTH WAYS — RYA-713.

    python3 scripts/control_synthesis_handler.py --element Fe --ion I --n 20

Ryan, 2026-08-09: *"Lets look at it from both angles. I want to see both results."*

The synthesis handler natively produces an **abundance**; the profile fitter produces an
**equivalent width**. Controlling only one of them would leave a real failure mode
invisible, so this runs both against the same lines and reports them side by side.

  ANGLE 1 — EW LEVEL
      synthetic EW at the best-fit abundance   vs   the banked HARPS EW
      Tests: does the synthesis reproduce the observed line STRENGTH?
      Comparable to the profile fitter's control (median ratio 0.971, −0.0129 dex),
      so the two methods can be judged on the same axis.

  ANGLE 2 — ABUNDANCE LEVEL
      fitted A(X) per line                     vs   the banked value (A(Fe I) = 7.466)
      Tests: does the synthesis recover the known ANSWER?
      This exercises the whole chain — atmosphere, line list, oscillator strengths,
      broadening — not just the spectrum match.

WHY BOTH, AND WHAT THEIR DISAGREEMENT MEANS
-------------------------------------------
They can diverge, and the direction is diagnostic:

  EW agrees, abundance does not   -> the spectrum is being matched but converted wrongly:
                                     suspect the gf scale, the atmosphere, or the
                                     curve-of-growth regime.
  abundance agrees, EW does not   -> compensating errors. Treat with MORE suspicion than
                                     a clean failure, because it looks like success.
  both agree                      -> the handler is controlled.
  neither agrees                  -> a plain failure, and the easiest to diagnose.

RUNS ON SIRIUS. Needs iSpec, the MARCS atmospheres and the GES line list.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.band_policy import resolve as resolve_band  # noqa: E402
from pipeline.measure import resolve_handler  # noqa: E402
from pipeline import intake_debug as dbg  # noqa: E402  RYA-770/765 per-line visibility
from pipeline.measure.base import ControlResult, CONTROL_TOLERANCE_DEX  # noqa: E402
from scripts.measure_band_ew import kp_segments, load_kp_window  # noqa: E402
from config.constants import codex_path  # RYA-810 path register

OUT = ROOT / "data" / "audit" / "synthesis_control"
# The control band is the overlap of THREE things, not just the optical:
#   the VIS policy band, the banked HARPS pool, and the SYNTHESIS LINE LIST.
# GES v6 spans 4200-9200 A. Running the control from 3924 A included lines the
# synthesiser cannot see at all -- 4065.381 A synthesised to EW 0.00 against a pool
# value of 73.59 and counted as a failure of the METHOD when it was a coverage gap.
GES_LO, GES_HI = 4200.0, 9200.0
CONTROL_LO, CONTROL_HI = max(3924.0, GES_LO), min(6905.0, GES_HI)


# iSpec is a source tree on Sirius, not a pip install.
ISPEC_SRC = str(codex_path('engines.ispec'))
# RYA-682: numpy >= 2.3 breaks ispec/abundances.py:132, and the failure is SILENT --
# synthesis-v2 writes a 0-usable-row artifact and exits 0. A control that "passes" on
# zero rows is worse than one that fails, so this is checked loudly and first.
NUMPY_CEILING = (2, 3)


def _guard_environment() -> None:
    if ISPEC_SRC not in sys.path:
        sys.path.insert(0, ISPEC_SRC)
    import numpy
    ver = tuple(int(x) for x in numpy.__version__.split(".")[:2])
    if ver >= NUMPY_CEILING:
        raise SystemExit(
            f"numpy {numpy.__version__} is at or above the {NUMPY_CEILING[0]}."
            f"{NUMPY_CEILING[1]} ceiling ratified in RYA-682. Above it, "
            f"ispec/abundances.py:132 fails and synthesis writes a ZERO-ROW artifact "
            f"while exiting 0 — a control run would report success on no data. "
            f"Use an interpreter with numpy < 2.3 (venv312 or venv_pysme on Sirius).")
    try:
        import ispec  # noqa: F401
    except ImportError as e:
        raise SystemExit(f"iSpec not importable even with {ISPEC_SRC} on sys.path: {e}")


def build_context(element: str, ion: str, resolving_power: float) -> dict:
    """Assemble the synth context from the pipeline's own loaders — never invented."""
    _guard_environment()
    from pipeline.abundances_derive import (
        _load_atmosphere, _load_synth_resources, _ISPEC_SOLAR_ABUND_FILE)
    # _atom_codes lives in cno_synthesis, not abundances_derive -- import it from where
    # it actually is rather than duplicating six lines that call iSpec.
    from pipeline.cno_synthesis import _atom_codes, _solar_A as _ispec_solar_A_map


    def _ispec_solar_A(el, chem, sab):
        """A(X) on iSpec's internal scale -- independent of our banked reference."""
        return _ispec_solar_A_map([el], chem, sab)[el]
    import ispec
    from config.constants import STAR_PARAMS

    p = STAR_PARAMS["solar"]
    teff, logg = float(p["teff"]), float(p["logg"])
    # Microturbulence is `xi` in stars.yaml, not `vmic`. My first pass used
    # p.get("vmic", p.get("microturbulence", 1.0)) and fell through to a HARDCODED 1.0.
    # It happens to equal the true solar xi, so this run was not biased -- but the next
    # star would have silently inherited the Sun's value, which is the exact failure
    # RYA-288 forbids for broadening and which applies identically here: microturbulence
    # changes the curve of growth and therefore biases A(X), not just its error bar.
    if "xi" not in p:
        raise SystemExit(
            "no microturbulence (`xi`) in STAR_PARAMS for the Sun. Refusing to default "
            "it -- xi sets the saturation regime and biases the fitted abundance.")
    vturb = float(p["xi"])
    feh = float(p.get("feh", p.get("feh_ref", 0.0)))

    # Broadening from the project's own resolver -- never invented, never fitted here.
    # RYA-288: "Wrong broadening biases the fitted A(X) itself, not just its error bar,
    # so there is NO default." The resolver returns HARPS_R for R, which we DISCARD:
    # resolution is an INSTRUMENT constant and this control runs on Kitt Peak.
    from pipeline.abundances_derive import _resolve_broadening
    _r_unused, vmac, vsini, fit_vmac = _resolve_broadening("solar")
    if fit_vmac:
        raise SystemExit(
            "STAR_PARAMS asks for vmac='fit' on the Sun. Fitting vmac by full-profile "
            "chi2 rails to the upper bound at high SNR (RYA-309), so this control will "
            "not do it -- set a fixed solar vmac or adopt a literature RT value.")

    linelist, isotopes, chem = _load_synth_resources()
    solar_abund = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    atm = _load_atmosphere(teff, logg, feh, vturb)
    codes = _atom_codes([element], chem, solar_abund)

    return dict(atmosphere=atm, teff=teff, logg=logg, feh=feh, vturb=vturb,
                linelist=linelist, isotopes=isotopes, solar_abund=solar_abund,
                atom_code=int(codes[element]), resolving_power=float(resolving_power),
                macroturbulence=float(vmac), vsini=float(vsini),
                # CIRCULARITY GUARD (RYA-713). The starting guess must NOT be the value
                # the control is testing. It was: a_start defaulted to solar_A, which I
                # had set to SOLAR_REFERENCE -- so any line whose chi2 was flat returned
                # the reference EXACTLY and counted as a perfect result. Three of eight
                # lines did precisely that (dA = -3.6e-15), silently pulling the median
                # toward a pass.
                #
                # Seed from iSpec's own internal solar composition instead: an
                # independent starting point that knows nothing about our banked answer.
                solar_A=float(_ispec_solar_A(element, chem, solar_abund)),
                span_dex=1.5)


# Banked reference values. These are the ANSWER the control must recover -- they are NOT
# a target to tune toward, and nothing in the handler can see them.
SOLAR_REFERENCE = {"Fe": 7.466}   # RYA-553, 1D-NLTE on the 3D scale
REFERENCE_SOURCE = {"Fe": "RYA-553 banked A(Fe I) = 7.466 (Asplund 2021 gives 7.46)"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--ion", default="I")
    ap.add_argument("--n", type=int, default=20, help="lines to control on")
    ap.add_argument("--instrument", default="kpno_solar_atlas",
                    choices=["kpno_solar_atlas", "harps"],
                    help="harps uses the pipeline's normalised solar spectrum -- the SAME "
                         "spectrum synth-v2 measured A(Fe I)=7.520 on, so the comparison "
                         "isolates the handler rather than the instrument.")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--extract-windows", metavar="NPZ",
                    help="MAC SIDE: cut the spectral windows and write them here, then exit. "
                         "Spectra live on the Mac (Ryan, RYA-567) and compute lives on "
                         "Sirius, so we ship a few hundred kB of windows rather than "
                         "staging a 2.3 GB atlas.")
    ap.add_argument("--windows", metavar="NPZ",
                    help="SIRIUS SIDE: read windows from this file instead of the atlas.")
    ap.add_argument("--lines-from-banked", action="store_true",
                    help="draw the control lines from the banked synth-v2 per-line "
                         "table instead of by EW range, so the handler is compared "
                         "against the reference on the SAME lines (RYA-770)")
    ap.add_argument("--star", default="solar",
                    help="star namespace for the banked per-line table")
    ap.add_argument("--trace", action="store_true",
                    help="write a local, git-ignored intake trace of every per-line "
                         "accept/reject to debug/intake/ (same as CODEX_INTAKE_TRACE=1)")
    a = ap.parse_args()

    # RYA-770: the handler emits a tracer event for every accepted and rejected line,
    # but those are no-ops unless a trace is actually started -- the env var alone does
    # nothing here, because the gate lives in each driver. Without this the control
    # could run with CODEX_INTAKE_TRACE=1 set and still record nothing, which is the
    # silent-instrumentation trap the tracer exists to prevent.
    tracing = a.trace or os.environ.get("CODEX_INTAKE_TRACE") == "1"
    if tracing:
        dbg.start_trace(a.element, f"synthesis_control_{a.instrument}", "solar",
                        extra={"ticket": "RYA-770", "n_lines": a.n,
                               "driver": "scripts/control_synthesis_handler.py"})
    try:
        _run(a)
    finally:
        if tracing:
            dbg.end_trace()


def _run(a) -> None:
    pool = pd.read_csv(ROOT / "data" / "measured" / "sol_ew_results_v1.csv")
    if a.lines_from_banked:
        # RYA-770: control on the lines the VALIDATED path actually fitted.
        #
        # Selecting by EW range instead (the default below) draws from the measured
        # pool, and those lines share NOTHING with the banked synth-v2 set -- measured:
        # 0 of 12 in common. So "the adapter got 7.306 where synth-v2 got 7.520" was
        # comparing two different line samples and attributing the difference to the
        # handler. Same engine and same spectrum are not enough; the control has to
        # hold the LINES fixed too, or it is not isolating the handler.
        #
        # This takes every line synth-v2 ATTEMPTED, not the ones it accepted, so no
        # selection on the answer enters: the adapter applies its own gate and the
        # comparison is of what survives, not of a set pre-filtered by the reference.
        bank = pd.read_csv(ROOT / "data" / "outputs" / a.star / f"{a.star}_per_line_synth_v2.csv")
        bank = bank[(bank.element == a.element) & (bank.ion == a.ion)].copy()
        if not len(bank):
            raise SystemExit(
                f"--lines-from-banked: no {a.element} {a.ion} rows in the banked "
                f"synth-v2 per-line table. It is a GENERATED artifact (gitignored); "
                f"regenerate it on Sirius before controlling against it.")
        bank = bank.sort_values("wavelength_air_A").reset_index(drop=True)
        idx = np.unique(np.linspace(0, len(bank) - 1, a.n).astype(int))
        bank = bank.iloc[idx]
        # The banked table carries the EW production used for the window, so take it
        # from there rather than re-joining to the measured pool: the join lost 18 of
        # 24 lines to a 0.05 A tolerance, and any line it did match could have been
        # matched to a DIFFERENT pool row than the one production used -- which would
        # silently change the window and therefore what is being compared.
        ref = pd.DataFrame({
            "wavelength_air_A": bank.wavelength_air_A.astype(float).to_numpy(),
            "ew_mA": bank.ew_mA.astype(float).to_numpy(),
            "banked_a_synth": bank.a_synth.astype(float).to_numpy(),
            "banked_red_chi2": bank.red_chi2.astype(float).to_numpy(),
        }).reset_index(drop=True)
        print(f"  controlling on {len(ref)} lines drawn from the BANKED synth-v2 set "
              f"({len(bank)} sampled of {len(bank)}); banked median A = "
              f"{ref.banked_a_synth.median():.3f}")
    else:
        ref = pool[(pool.element == a.element) & (pool.ion == a.ion) &
                   pool.wavelength_air_A.between(CONTROL_LO, CONTROL_HI) &
                   pool.ew_mA.between(10, 90)].copy()      # unsaturated, well measured
        ref = ref.sort_values("wavelength_air_A").reset_index(drop=True)
        idx = np.unique(np.linspace(0, len(ref) - 1, a.n).astype(int))
        ref = ref.iloc[idx].reset_index(drop=True)

    cat = pd.read_csv(ROOT / "data" / "catalog" / "instrument_catalog.csv")
    R = float(cat[cat.iloc[:, 0].astype(str) == a.instrument].iloc[0]["resolving_power_max"])

    # Which reference are we testing against? They are different quantities and must not
    # be conflated: 7.466 is the banked EW-path value (1D-NLTE on the 3D scale); 7.520 is
    # what synth-v2 -- THIS engine, on THIS spectrum -- produced. Controlling the handler
    # against synth-v2 asks "does my driver reproduce the project's validated synthesis?",
    # which is the question. Comparing to 7.466 would fold in NLTE and 3D corrections the
    # synthesis never applied.
    if a.instrument == "harps":
        SOLAR_REFERENCE["Fe"] = 7.520
        REFERENCE_SOURCE["Fe"] = ("banked synth-v2 A(Fe I) = 7.520 (23 lines, HARPS, "
                                  "1D-LTE, nlte_applied=False) -- same engine, same "
                                  "spectrum, so this isolates the handler")

    if a.extract_windows:
        segs = kp_segments()
        store, kept = {}, []
        for _, r in ref.iterrows():
            c = float(r.wavelength_air_A)
            try:
                w, f, _ = load_kp_window(segs, c, pad=1.4)
            except Exception as e:
                print(f"  skip {c:.3f}: {e}"); continue
            store[f"w_{c:.4f}"] = w.astype("float32")
            store[f"f_{c:.4f}"] = f.astype("float32")
            kept.append((c, float(r.ew_mA)))
        store["index"] = np.array(kept, dtype="float64")
        out = Path(a.extract_windows); out.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out, **store)
        mb = out.stat().st_size / 1e6
        print(f"  extracted {len(kept)} windows -> {out} ({mb:.2f} MB)")
        return

    policy = resolve_band(0.5 * (CONTROL_LO + CONTROL_HI))
    handler = resolve_handler(3400.0)          # the synthesis handler
    print(f"controlling {handler.name} on {a.element} {a.ion}, {len(ref)} lines, "
          f"band={policy.name}, R={R:.0f}")
    print(f"reference: {REFERENCE_SOURCE.get(a.element, '(none)')}\n")

    ctx = build_context(a.element, a.ion, R)
    handler.prepare(policy, ctx)

    if a.instrument == "harps":
        from pipeline.abundances_derive import _load_observed_spectrum
        ow_nm, oflux = _load_observed_spectrum("solar")
        ow = ow_nm * 10.0
        print(f"  HARPS normalised solar: {ow.size} pts, {ow.min():.1f}-{ow.max():.1f} A")

        def get_window(c):
            m = np.abs(ow - c) <= 1.4
            if m.sum() < 12:
                raise LookupError(f"only {int(m.sum())} HARPS pixels near {c:.3f} A")
            return ow[m], oflux[m]
    elif a.windows:
        z = np.load(a.windows)
        index = z["index"]
        ref = pd.DataFrame({"wavelength_air_A": index[:, 0], "ew_mA": index[:, 1]})
        def get_window(c):
            return z[f"w_{c:.4f}"].astype(float), z[f"f_{c:.4f}"].astype(float)
        print(f"  reading {len(ref)} pre-extracted windows from {a.windows}")
    else:
        segs = kp_segments()
        def get_window(c):
            w, f, _ = load_kp_window(segs, c, pad=1.4)
            return w, f

    rows = []
    for _, r in ref.iterrows():
        c = float(r.wavelength_air_A)
        try:
            w, f = get_window(c)
        except Exception as e:
            rows.append(dict(wave=c, ew_ref=float(r.ew_mA), ew_synth=float("nan"),
                             abundance=float("nan"), status=f"load: {e}")); continue
        # The window is wing-wide from the line's own strength, so the handler needs to
        # know roughly how strong it is. The banked EW is the honest hint -- it sets the
        # WINDOW, never the answer.
        lm = handler.measure_line(w, f, element=a.element, ion=a.ion, wavelength_A=c,
                                  instrument=a.instrument, policy=policy,
                                  # HARPS arrives globally normalised by our pipeline;
                                  # Kitt Peak arrives as atlas residual flux. Both are
                                  # already on a continuum, so neither wants a second one.
                                  pre_normalised=True,
                                  context={**ctx, "ew_hint_mA": float(r.ew_mA)})
        rows.append(dict(wave=c, ew_ref=float(r.ew_mA), ew_synth=float(lm.ew_mA),
                         abundance=lm.abundance,
                         status="ok" if lm.in_aggregate else lm.excluded_reason[:70]))

    d = pd.DataFrame(rows)
    for col in ("ew_synth", "abundance"):
        if col not in d.columns:
            d[col] = float("nan")
    ok = d[(d.status == "ok") & d.abundance.notna()].copy()
    if not len(ok):
        # Say WHY every line failed instead of dying on a missing column -- a control
        # that cannot report its own failure mode is useless.
        print("  no line produced an abundance. statuses:")
        for st, n in d.status.value_counts().items():
            print(f"    {n:3d}  {st}")
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out / f"synthesis_control_{a.element}{a.ion}.csv", index=False)

    print(f"measured {len(ok)} of {len(d)}\n")
    if not len(ok):
        print("  no lines measured — cannot form either angle")
        return

    # ── ANGLE 1: EW level ────────────────────────────────────────────────────
    ok["ratio"] = ok.ew_synth / ok.ew_ref
    med = float(ok.ratio.median()); mad = float(np.median(np.abs(ok.ratio - med)))
    dex_ew = float(np.log10(med)) if med > 0 else float("nan")

    # ── ANGLE 2: abundance level ─────────────────────────────────────────────
    A_ref = SOLAR_REFERENCE[a.element]
    a_med = float(ok.abundance.median())
    a_sc = float(np.std(ok.abundance, ddof=1)) if len(ok) > 1 else float("nan")
    dex_ab = a_med - A_ref

    print(f"{'':30s}{'result':>12s}{'vs reference':>15s}")
    print(f"{'ANGLE 1 — EW level':30s}")
    print(f"{'  median synth/banked EW':30s}{med:12.3f}{dex_ew:+15.4f} dex")
    print(f"{'  MAD of ratio':30s}{mad:12.3f}")
    print(f"{'  within 10%':30s}{100*(ok.ratio.sub(1).abs()<=0.10).mean():11.0f}%")
    print(f"{'ANGLE 2 — abundance level':30s}")
    print(f"{'  median fitted A(X)':30s}{a_med:12.3f}{dex_ab:+15.4f} dex")
    print(f"{'  line-to-line scatter':30s}{a_sc:12.3f}")
    print(f"{'  reference A(X)':30s}{A_ref:12.3f}")
    print()
    print(f"  profile fitter, for comparison: EW ratio 0.971, -0.0129 dex, MAD 0.060")

    pass_ew = abs(dex_ew) <= CONTROL_TOLERANCE_DEX
    pass_ab = abs(dex_ab) <= CONTROL_TOLERANCE_DEX
    verdict = {
        (True, True):   "BOTH AGREE — handler is controlled",
        (True, False):  "EW AGREES, ABUNDANCE DOES NOT — the spectrum is matched but "
                        "converted wrongly: suspect the gf scale, the atmosphere, or the "
                        "curve-of-growth regime",
        (False, True):  "ABUNDANCE AGREES, EW DOES NOT — compensating errors. Treat with "
                        "MORE suspicion than a clean failure, because it looks like success",
        (False, False): "NEITHER AGREES — a plain failure, and the easiest to diagnose",
    }[(pass_ew, pass_ab)]
    print(f"\n  {verdict}")

    # The handler is controlled on the ABUNDANCE angle: that is what it produces, and
    # what it will be asked for in the near-UV. The EW angle is recorded as evidence.
    result = ControlResult(
        handler=handler.name, element=a.element, n_lines=len(ok),
        median_ratio=med, mad_ratio=mad, dex_offset=dex_ab, passed=pass_ab,
        tolerance_dex=CONTROL_TOLERANCE_DEX,
        evidence=(f"ANGLE 1 EW: median synth/banked {med:.3f} ({dex_ew:+.4f} dex), "
                  f"MAD {mad:.3f}. ANGLE 2 abundance: median A={a_med:.3f} vs "
                  f"{A_ref:.3f} ({dex_ab:+.4f} dex), scatter {a_sc:.3f}. {verdict}"))
    (out / f"control_{a.element}{a.ion}.json").write_text(
        json.dumps(result.__dict__, indent=2))
    print(f"  {result.summary()}")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()

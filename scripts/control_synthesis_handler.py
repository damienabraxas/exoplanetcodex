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
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.band_policy import resolve as resolve_band  # noqa: E402
from pipeline.measure import resolve_handler  # noqa: E402
from pipeline.measure.base import ControlResult, CONTROL_TOLERANCE_DEX  # noqa: E402
from scripts.measure_band_ew import kp_segments, load_kp_window  # noqa: E402

OUT = ROOT / "data" / "audit" / "synthesis_control"
CONTROL_LO, CONTROL_HI = 3924.0, 6905.0


# iSpec is a source tree on Sirius, not a pip install.
ISPEC_SRC = "/srv/codex/engines/ispec_src"
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
    from pipeline.cno_synthesis import _atom_codes
    import ispec
    from config.constants import STAR_PARAMS

    p = STAR_PARAMS["solar"]
    teff, logg = float(p["teff"]), float(p["logg"])
    feh = float(p.get("feh", 0.0))
    vturb = float(p.get("vmic", p.get("microturbulence", 1.0)))

    linelist, isotopes, chem = _load_synth_resources()
    solar_abund = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    atm = _load_atmosphere(teff, logg, feh, vturb)
    codes = _atom_codes([element], chem, solar_abund)

    return dict(atmosphere=atm, teff=teff, logg=logg, feh=feh, vturb=vturb,
                linelist=linelist, isotopes=isotopes, solar_abund=solar_abund,
                atom_code=int(codes[element]), resolving_power=float(resolving_power),
                solar_A=float(SOLAR_REFERENCE[element]), span_dex=1.0)


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
    ap.add_argument("--instrument", default="kpno_solar_atlas")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--extract-windows", metavar="NPZ",
                    help="MAC SIDE: cut the spectral windows and write them here, then exit. "
                         "Spectra live on the Mac (Ryan, RYA-567) and compute lives on "
                         "Sirius, so we ship a few hundred kB of windows rather than "
                         "staging a 2.3 GB atlas.")
    ap.add_argument("--windows", metavar="NPZ",
                    help="SIRIUS SIDE: read windows from this file instead of the atlas.")
    a = ap.parse_args()

    pool = pd.read_csv(ROOT / "data" / "measured" / "sol_ew_results_v1.csv")
    ref = pool[(pool.element == a.element) & (pool.ion == a.ion) &
               pool.wavelength_air_A.between(CONTROL_LO, CONTROL_HI) &
               pool.ew_mA.between(10, 90)].copy()          # unsaturated, well measured
    ref = ref.sort_values("wavelength_air_A").reset_index(drop=True)
    idx = np.unique(np.linspace(0, len(ref) - 1, a.n).astype(int))
    ref = ref.iloc[idx].reset_index(drop=True)

    cat = pd.read_csv(ROOT / "data" / "catalog" / "instrument_catalog.csv")
    R = float(cat[cat.iloc[:, 0].astype(str) == a.instrument].iloc[0]["resolving_power_max"])

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

    if a.windows:
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
        lm = handler.measure_line(w, f, element=a.element, ion=a.ion, wavelength_A=c,
                                  instrument=a.instrument, policy=policy,
                                  pre_normalised=True, context=ctx)
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

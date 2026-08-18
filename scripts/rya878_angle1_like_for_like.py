#!/usr/bin/env python3
"""
RYA-878 — make ANGLE 1 a MEASURED angle, then read what it says.
================================================================
    ISPEC_DIR=... python3 scripts/rya878_angle1_like_for_like.py

The SynthesisHandler control's ANGLE 1 reports a synthetic-vs-observed EW ratio of ~1.433
(+0.1562 dex). RYA-875 settled ANGLE 2 — the abundance residual is 0.0000, paired per
line — and left this one with a specific defect that is NOT "we do not know the cause":

🔴 ANGLE 1 IS AN INFERENCE WEARING AN ANGLE. The claim that the production engine shares
the offset (i.e. that it is not the adapter's fault) could not be measured, because the
production synth-v2 path BANKS NO SYNTHETIC EW — it records `A_X`, `red_chi2`, `status`
and nothing else. There was no engine-side number to compare against, so RYA-875 could
only argue it: the fitted abundance matches the banked one to <=0.009 dex and EW is
monotonic in A, therefore the engine at its own answer would produce nearly the same EW.
Sound, and still an inference.

This measures it instead, and it does so WITHOUT re-fitting. The synthetic EW is a pure
function of (abundance, line data, window, broadening) — no chi2 loop — so the engine-side
EW can be evaluated directly at the BANKED `a_synth` using the same `pipeline.synth_ew`
definition the handler now calls. Two syntheses per line rather than a full refit.

THE THREE COMPARISONS, AND WHAT EACH SEPARATES
----------------------------------------------
    (a) handler EW   vs  engine EW      -> does the ADAPTER reproduce the ENGINE?
                                           like-for-like: both synthetic, same definition.
                                           This is the number RYA-875 could only infer.
    (b) engine EW    vs  observed EW    -> is the offset a property of the MODEL?
                                           if (a) ~ 1 and (b) ~ 1.43, the engine carries
                                           it too and the adapter is exonerated by
                                           measurement rather than by argument.
    (c) handler EW   vs  observed EW    -> ANGLE 1 as it stands today, for continuity.

ALREADY REFUTED — DO NOT RE-TEST (RYA-873/875)
-----------------------------------------------
1. same-species neighbour contamination surviving the difference-of-syntheses;
2. a uniform continuum/normalisation offset (log10 std 0.247, a factor-15 spread);
3. curve-of-growth saturation (corr -0.047 vs reduced EW);
4. line strength (-0.100);   5. red_chi2 (+0.059);   6. sigma_A (-0.169).

⚠️ This script changes no abundance, no bar and no harness residual. ANGLE 2 is closed at
0.0 and RYA-878 does not touch it.
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
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import synth_ew                                      # noqa: E402


OUT = ROOT / "data" / "results" / "rya878"
CONTROL = ROOT / "data" / "audit" / "synthesis_control" / "synthesis_control_FeI.csv"


def _ratio_stats(num, den, label):
    r = (np.asarray(num, float) / np.asarray(den, float))
    m = np.isfinite(r) & (r > 0)
    r = r[m]
    lg = np.log10(r)
    return {"comparison": label, "n": int(r.size),
            "median_ratio": round(float(np.median(r)), 4),
            "dex": round(float(np.log10(np.median(r))), 4),
            "mad_ratio": round(float(np.median(np.abs(r - np.median(r)))), 4),
            "log10_std": round(float(lg.std(ddof=1)), 4) if r.size > 1 else None,
            "min": round(float(r.min()), 4), "max": round(float(r.max()), 4),
            "within_10pct": int(np.sum(np.abs(r - 1.0) <= 0.10))}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--control", type=Path, default=CONTROL)
    ap.add_argument("--star", default="solar")
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--ion", default="I")
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()

    ctl = pd.read_csv(a.control)
    ok = ctl[(ctl.status == "ok") & ctl.ew_synth.notna()].copy()
    if not len(ok):
        raise SystemExit("no accepted control lines")

    banked_p = ROOT / "data" / "outputs" / a.star / f"{a.star}_per_line_synth_v2.csv"
    if not banked_p.exists():
        raise SystemExit(
            f"banked per-line synth-v2 absent: {banked_p}\n"
            f"  It is a GENERATED, gitignored artifact. Link or regenerate it on Sirius.")
    bank = pd.read_csv(banked_p)
    bank = bank[(bank.element == a.element) & (bank.ion == a.ion)]
    m = ok.merge(bank[["wavelength_air_A", "ew_mA", "a_synth"]].rename(
        columns={"wavelength_air_A": "wave", "ew_mA": "ew_obs",
                 "a_synth": "a_engine"}), on="wave", how="left")
    assert m.a_engine.notna().all(), "a control line is absent from the banked set"

    # If the banked table ALREADY carries a synthetic EW, the engine banked it itself and
    # nothing needs recomputing — that is the end state RYA-878 wires up for future runs.
    if "ew_synth_mA" in bank.columns:
        print("[engine EW] the banked table already carries ew_synth_mA — using it")
        m = m.merge(bank[["wavelength_air_A", "ew_synth_mA"]].rename(
            columns={"wavelength_air_A": "wave", "ew_synth_mA": "ew_engine"}),
            on="wave", how="left")
    else:
        print("[engine EW] the banked table predates RYA-878 and carries no synthetic EW")
        print("            — evaluating it at the BANKED a_synth with the same")
        print("              pipeline.synth_ew definition the handler uses. No refit:")
        print("              EW is a pure function of the abundance, not of the search.")
        from control_synthesis_handler import build_context
        from pipeline.measure import resolve_handler
        from pipeline.band_policy import resolve as resolve_band
        cat = pd.read_csv(ROOT / "data" / "catalog" / "instrument_catalog.csv")
        R = float(cat[cat.iloc[:, 0].astype(str) == "harps"].iloc[0]["resolving_power_max"])
        ctx = build_context(a.element, a.ion, R)
        handler = resolve_handler(3400.0)
        handler.prepare(resolve_band(5000.0), ctx)
        kw = handler._synth_kw(ctx, a.element)
        rows = []
        for r in m.itertuples():
            c = float(r.wave)
            wbase, wtop = handler._window(c / 10.0, float(r.ew_obs))
            hw_A = (wtop - wbase) * 10.0 / 2.0
            rows.append(synth_ew.synthetic_ew_mA(
                handler._synth, centre_A=c, abundance=float(r.a_engine),
                fit_half_width_A=hw_A, synth_kwargs=kw))
            print(f"  {c:10.3f}  A_engine={float(r.a_engine):.3f}  "
                  f"EW_engine={rows[-1]:9.2f}  EW_handler={float(r.ew_synth):9.2f}")
        m["ew_engine"] = rows

    stats = [
        _ratio_stats(m.ew_synth, m.ew_engine,
                     "(a) handler EW / engine EW — does the ADAPTER reproduce the ENGINE?"),
        _ratio_stats(m.ew_engine, m.ew_obs,
                     "(b) engine EW / observed EW — is the offset the MODEL's?"),
        _ratio_stats(m.ew_synth, m.ew_obs,
                     "(c) handler EW / observed EW — ANGLE 1 as it stands"),
    ]
    print(f"\n=== ANGLE 1, decomposed ===")
    for s in stats:
        print(f"  {s['comparison']}")
        print(f"      n={s['n']}  median {s['median_ratio']:.3f} ({s['dex']:+.4f} dex)  "
              f"MAD {s['mad_ratio']:.3f}  range {s['min']:.2f}-{s['max']:.2f}  "
              f"within 10%: {s['within_10pct']}/{s['n']}")

    a_stat, b_stat, c_stat = stats
    adapter_faithful = abs(a_stat["dex"]) <= 0.02
    print(f"\n=== verdict ===")
    if adapter_faithful and abs(b_stat["dex"]) > 0.05:
        print(f"  MEASURED, not inferred: the adapter reproduces the engine's synthetic EW")
        print(f"  to {a_stat['dex']:+.4f} dex, and the ENGINE ITSELF sits {b_stat['dex']:+.4f}")
        print(f"  dex above the observed EW. So ANGLE 1's offset is a property of the")
        print(f"  MODEL-vs-OBSERVATION comparison, shared by the production engine — not")
        print(f"  an adapter defect. RYA-875 argued this; it is now a number.")
    elif not adapter_faithful:
        print(f"  🔴 THE ADAPTER DOES NOT REPRODUCE THE ENGINE: {a_stat['dex']:+.4f} dex.")
        print(f"  That REFUTES RYA-875's inference and makes ANGLE 1 an adapter finding")
        print(f"  after all. This is the outcome that would matter — chase it.")
    else:
        print(f"  the engine matches the observed EW ({b_stat['dex']:+.4f} dex) while the")
        print(f"  handler does not — the offset is the ADAPTER's. Chase it.")

    a.out.mkdir(parents=True, exist_ok=True)
    m.to_csv(a.out / "rya878_angle1_by_line.csv", index=False)
    (a.out / "rya878_angle1_summary.json").write_text(json.dumps({
        "ticket": "RYA-878",
        "control": str(a.control), "banked": str(banked_p),
        "n_lines": int(len(m)),
        "engine_ew_source": ("banked column" if "ew_synth_mA" in bank.columns
                             else "evaluated at the banked a_synth (no refit)"),
        "comparisons": stats,
        "adapter_reproduces_engine": bool(adapter_faithful),
        "already_refuted_do_not_retest": [
            "same-species neighbour contamination (RYA-873)",
            "uniform continuum/normalisation offset (log10 std 0.247)",
            "curve-of-growth saturation (corr -0.047 vs reduced EW)",
            "line strength (-0.100)", "red_chi2 (+0.059)", "sigma_A (-0.169)"],
    }, indent=2) + "\n")
    print(f"\n  wrote {a.out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

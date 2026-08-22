#!/usr/bin/env python3
"""Did synthesis recover the lines EW dropped, and which gate had killed them? — RYA-967.

    python3 scripts/rya967_survivors.py --arm kpno_solar_atlas

🔴 A RAW RECOVERY COUNT IS NOT THE RESULT. "Synthesis recovers 300 lines" says nothing
about whether the route is doing what the physics predicts. The claim RYA-958/959 rest on
is specific: EW dies on SATURATED lines because the curve of growth goes flat, while
synthesis fits the profile and inverts no equivalent width, so the saturation ceiling
never applies (`pipeline/measure/synthesis.py:277`). That predicts a SHAPE — recovery
should be high for REW-SATURATED and for OVER-PHYSICAL-WIDTH, and low for
FEATURE-VERIFICATION, because a blend or a ghost is a property of the SPECTRUM that
changes method cannot fix.

So the output is a contingency table: EW loss reason x synth outcome. If saturation
recovers and blends do not, the mechanism is confirmed. If everything recovers equally,
something is wrong with the synth leg — most likely it is fitting whatever is in the
window rather than the named line, which is the RYA-843 failure that returned ~12.49 in
the NIR.

⚠️ WHAT THIS CANNOT SAY. A recovered line is a line the synthesis route ACCEPTED, not a
line whose abundance is right. The two legs' A(Fe) are compared per line, and a large
spread among "recovered" lines is itself a finding, not a rounding error.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

EW_DIR = ROOT / "data" / "measured" / "band_ew"
PROD_DIR = ROOT / "data" / "results" / "band_products"
OUT = ROOT / "data" / "audit" / "rya967_ew_vs_synth"

#: Same key the selector uses, so "this line" means the same thing on both sides.
TOL_A = 0.005


def _verdict(reason: object) -> str:
    """Collapse an EW exclusion sentence to its coded gate (RYA-959's census rule)."""
    if not isinstance(reason, str) or not reason:
        return "IN-AGGREGATE"
    if reason.startswith("REW "):
        return "REW-SATURATED"
    return reason.split(":")[0]


def _match(a: np.ndarray, b: np.ndarray, tol: float = TOL_A) -> np.ndarray:
    """Index into `b` for each element of `a`, or -1. Nearest-within-tolerance, never a
    rounded key (RYA-703/704 are two tickets about exactly that mistake)."""
    order = np.argsort(b)
    bs = b[order]
    out = np.full(a.shape, -1, dtype=int)
    for n, w in enumerate(a):
        i = int(np.searchsorted(bs, w))
        best, bi = np.inf, -1
        for k in (i - 1, i, i + 1):
            if 0 <= k < len(bs) and abs(bs[k] - w) < best:
                best, bi = abs(bs[k] - w), k
        if bi >= 0 and best <= tol:
            out[n] = order[bi]
    return out


def run(arm: str, lo: float, hi: float, ew_name: str, synth_glob: str) -> dict:
    ew_path = EW_DIR / ew_name
    if not ew_path.exists():
        raise SystemExit(f"no EW artifact at {ew_path}")
    ew = pd.read_csv(ew_path)
    ew = ew[ew.wavelength_air_A.between(lo, hi)].copy()
    ew["gate"] = ew.excluded_reason.map(_verdict)

    synth_files = sorted(PROD_DIR.glob(synth_glob))
    if not synth_files:
        raise SystemExit(
            f"no synthesis product matching {synth_glob} in {PROD_DIR} — run "
            f"derive_band_products.py --force-synthesis --lines-from-ew ... first. "
            f"Refusing to report a comparison with one side missing.")
    sy = pd.concat([pd.read_csv(f) for f in synth_files], ignore_index=True)
    sy = sy.drop_duplicates(subset=["wavelength_air_A"], keep="first")

    idx = _match(ew.wavelength_air_A.values.astype(float),
                 sy.wavelength_air_A.values.astype(float))
    ew["in_synth_run"] = idx >= 0
    ew["synth_ok"] = False
    ew["synth_A"] = np.nan
    hit = idx >= 0
    ew.loc[hit, "synth_ok"] = sy.in_aggregate.values[idx[hit]]
    ew.loc[hit, "synth_A"] = sy.abundance.values[idx[hit]]

    # Only lines the synth leg actually ATTEMPTED can be scored. A line absent from the
    # synthesis list was never a candidate, and counting it as "not recovered" would blame
    # the method for a line-list gap (RYA-833: an absence is a hypothesis).
    scored = ew[ew.in_synth_run].copy()
    tab = (scored.groupby("gate")
           .agg(attempted=("synth_ok", "size"),
                recovered=("synth_ok", "sum"))
           .assign(rate=lambda d: (d.recovered / d.attempted).round(3))
           .sort_values("attempted", ascending=False))
    return dict(arm=arm, ew=ew, scored=scored, table=tab,
                n_ew=len(ew), n_not_in_synth=int((~ew.in_synth_run).sum()),
                synth_files=[f.name for f in synth_files])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--arm", default="kpno_solar_atlas")
    ap.add_argument("--lo", type=float, default=4200.0)
    ap.add_argument("--hi", type=float, default=6910.0)
    ap.add_argument("--ew-name", default=None)
    ap.add_argument("--synth-glob", default=None)
    a = ap.parse_args(argv)

    ew_name = a.ew_name or f"FeI_3780_6910_{a.arm}_PROFILEFIT_ew.csv"
    glob = a.synth_glob or f"FeI_{int(a.lo)}_{int(a.hi)}_{a.arm}*1D-LTE_lines.csv"
    r = run(a.arm, a.lo, a.hi, ew_name, glob)

    print(f"\nRYA-967 — EW vs SYNTH on the same lines, {a.arm}, {a.lo:.0f}-{a.hi:.0f} A")
    print(f"  synth product(s): {r['synth_files']}")
    print(f"  EW lines in band: {r['n_ew']};  not attempted by synth "
          f"(absent from the synthesis list): {r['n_not_in_synth']}")
    print(f"\n  EW LOSS REASON            attempted  recovered   rate")
    for gate, row in r["table"].iterrows():
        print(f"    {gate:24s} {int(row.attempted):9d} {int(row.recovered):10d} "
              f"{row.rate:6.3f}")

    both = r["scored"]
    both = both[both.in_aggregate & both.synth_ok].dropna(subset=["synth_A"])
    if len(both):
        print(f"\n  Lines BOTH routes accepted: n={len(both)} — this is the agreement "
              f"check, and it is the one that says whether 'recovered' means 'right'.")
    OUT.mkdir(parents=True, exist_ok=True)
    r["table"].to_csv(OUT / f"survivors_by_reason_{a.arm}.csv")
    r["scored"].to_csv(OUT / f"per_line_{a.arm}.csv", index=False)
    print(f"\n  wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
scripts/line_accounting_rya709.py — RYA-709
===========================================
EVERY ELEMENT, EVERY LINE, AGAINST EVERY INSTRUMENT WE HOLD.

Ryan, 2026-08-09: *"I want to make it clear, we go through every element, and make
sure we are accounting for all the spectral data for those elements. We need to do it
for everything, is that clear?"*

This is that accounting, and it is a SCRIPT rather than a table because the answer
changes whenever an instrument is registered, a line list is updated, or a line is
measured — and a number that was true once and is quoted forever is how the project
got a HARPS-only pool nobody noticed.

What it answers, per element
----------------------------
  usable        distinct features with predicted central depth in [0.05, 0.60] — deep
                enough to carry signal, shallow enough to stay off the flat curve of
                growth. Calibrated on Al, whose four recovered lines sat at 0.21-0.37.
  reachable     of those, how many fall inside an instrument we actually hold, asked
                of `pipeline.coverage` and never of a loaded array.
  measured      how many appear in the committed EW pool.
  UNMEASURED    reachable minus measured — work we can do today with data on disk.
  no_instrument reachable by nothing we hold — an ACQUISITION target, not a task.

The caveat that must travel with every number here
--------------------------------------------------
`measured` counts the **EW pool only**. An element measured by SYNTHESIS has zero pool
lines and is not unmeasured: K I is gold with no pool line at all, via the Kitt Peak
synthesis channel, and C/N/O run multi-arm synthesis. Reading a zero here as "never
measured" would repeat, in the opposite direction, exactly the single-channel mistake
this audit exists to expose.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline import state_surfaces
from pipeline.gf_resolver import annotate_used_gf                       # noqa: E402
from pipeline.coverage import coverage_at                 # noqa: E402

#: Usable-depth window. A plotting/triage cutoff, never a science threshold — nothing
#: downstream reads it, and it is stated here so a reader can re-run with another.
DEPTH_LO, DEPTH_HI = 0.05, 0.60
#: Two linelist entries closer than this are one feature (HFS, isotopic, unresolved).
GROUP_A = 0.05
#: Pool-match tolerance. Deliberately nearest-within-tolerance rather than a rounded
#: key — RYA-703/704 are two tickets about exactly that mistake.
POOL_TOL_A = 0.10

OUT = ROOT / "data" / "audit" / "line_accounting"


def features(ll: pd.DataFrame, element: str) -> pd.DataFrame:
    a = ll[ll.element == element].sort_values("wavelength_air_A").copy()
    if a.empty:
        return a
    a["_k"] = (a.wavelength_air_A.diff().fillna(9e9) > GROUP_A).cumsum()
    return a.groupby("_k").agg(
        w=("wavelength_air_A", "mean"), ion=("ion", lambda s: s.mode().iat[0]),
        gf=("log_gf", "max"), ep=("excitation_potential_eV", "min"),
        d=("central_depth", "max"), n=("wavelength_air_A", "size")).reset_index(drop=True)


def audit(star: str = "solar") -> dict:
    ll = pd.read_csv(ROOT / "data" / "linelists" / "linelist_solar.csv", low_memory=False)
    pool = pd.read_csv(ROOT / "data" / "measured" / "sol_ew_results_v1.csv", comment="#")
    trk = pd.read_csv(ROOT / state_surfaces.TRACKER, comment="#")
    cache: dict[float, tuple[str, ...]] = {}

    def reach(w: float) -> tuple[str, ...]:
        k = round(float(w), 1)
        if k not in cache:
            cache[k] = tuple(i.instrument_id for i in coverage_at(w, star).covering)
        return cache[k]

    rows, per_line = [], []
    for el in sorted(set(trk.element)):
        g = features(ll, el)
        if g.empty:
            continue
        use = g[(g.d >= DEPTH_LO) & (g.d <= DEPTH_HI)].copy()
        pw = pool.loc[pool.element == el, "wavelength_air_A"].astype(float).values
        for _, r in use.iterrows():
            ins = reach(r.w)
            measured = bool(len(pw) and np.min(np.abs(pw - r.w)) <= POOL_TOL_A)
            per_line.append(dict(element=el, ion=r.ion, wave_air_A=round(float(r.w), 4),
                                 log_gf=round(float(r.gf), 3), ep_eV=round(float(r.ep), 4),
                                 predicted_depth=round(float(r.d), 4),
                                 instruments="|".join(ins), measured_in_pool=measured))
            # `log_gf` here is still the INTAKE value (max over the cluster). It is
            # corrected to the gf the inversion actually uses, once, below — doing it
            # per row would re-read the canonical index 21,828 times.
        n_reach = sum(1 for _, r in use.iterrows() if reach(r.w))
        n_meas = sum(1 for _, r in use.iterrows()
                     if len(pw) and np.min(np.abs(pw - r.w)) <= POOL_TOL_A)
        rows.append(dict(element=el, tier=str(trk.loc[trk.element == el, "tier"].iloc[0]),
                         pool_lines=int((pool.element == el).sum()), usable=len(use),
                         reachable=n_reach, measured=n_meas,
                         unmeasured_reachable=n_reach - n_meas,
                         no_instrument=len(use) - n_reach))
    # RYA-825 — make `log_gf` report the gf the INVERSION uses, not the intake value.
    #
    # `_load_synth_resources` defaults `apply_canonical_gf=True`, so the production
    # synthesis path resolves each line's gf out of canonical_gf.csv (RYA-353) BEFORE
    # inverting. This table was reporting the linelist's intake value instead, and the
    # two differ on 3,197 of 21,828 rows across twelve elements. RYA-799 read this column
    # and concluded the Fe pool was "100 % Kurucz K14"; RYA-824 re-inverted and found 18
    # of 29 lab-covered lines were already on the canonical value. The intake value is
    # preserved in `log_gf_intake` — nothing is lost, the primary column just stops lying.
    pl = pd.DataFrame(per_line)
    if len(pl):
        pl, gf_stats = annotate_used_gf(
            pl, wl_col="wave_air_A", ep_col="ep_eV", gf_col="log_gf",
            intake_source="linelist_solar.csv (VALD intake, max over the HFS cluster)")
        print(f"  gf column: {gf_stats['n_resolved']}/{gf_stats['n']} rows resolved to "
              f"canonical, {gf_stats['n_changed']} corrected, "
              f"{gf_stats['n_unresolved']} outside canonical (intake retained)")
    return dict(summary=pd.DataFrame(rows), per_line=pl)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--star", default="solar")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--element", action="append")
    args = ap.parse_args(argv)

    res = audit(args.star)
    s = res["summary"]
    if args.element:
        s = s[s.element.isin(args.element)]
    s = s.sort_values("unmeasured_reachable", ascending=False)

    print(f"{'el':4s}{'tier':12s}{'pool':>5s}{'usable':>7s}{'reach':>7s}{'meas':>6s}"
          f"{'UNMEASURED':>11s}{'no instr':>9s}")
    for _, r in s.iterrows():
        print(f"{r.element:4s}{r.tier:12s}{r.pool_lines:5d}{r.usable:7d}{r.reachable:7d}"
              f"{r.measured:6d}{r.unmeasured_reachable:11d}{r.no_instrument:9d}")
    print(f"\n{'TOTAL':4s}{'':12s}{s.pool_lines.sum():5d}{s.usable.sum():7d}"
          f"{s.reachable.sum():7d}{s.measured.sum():6d}"
          f"{s.unmeasured_reachable.sum():11d}{s.no_instrument.sum():9d}")
    print("\nCAVEAT: `measured` counts the EW POOL ONLY. A synthesis-measured element has "
          "zero pool lines and is NOT unmeasured — K I is gold with no pool line at all.")

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    res["summary"].to_csv(out / "summary.csv", index=False)
    res["per_line"].to_csv(out / "per_line.csv", index=False)
    (out / "thresholds.json").write_text(json.dumps(
        dict(ticket="RYA-709", depth_window=[DEPTH_LO, DEPTH_HI], group_A=GROUP_A,
             pool_tol_A=POOL_TOL_A,
             depth_window_basis="calibrated on Al's four recovered lines (0.21-0.37)",
             measured_means="present in the committed EW pool; synthesis channels are NOT counted"),
        indent=2))
    print(f"\nwrote {out}/summary.csv, per_line.csv, thresholds.json "
          f"({len(res['per_line'])} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

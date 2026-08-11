#!/usr/bin/env python3
"""RYA-760 — check the suspect-source log gf against NIST ASD, verbatim by curl.

RYA-760 says the FMW offset is "checkable against NIST ASD directly". It is, with one
caveat that has to be stated or the check is worthless:

    **FMW is Fuhr, Martin & Wiese — a NIST compilation.**

So if ASD still serves the FMW-era values, agreement proves nothing: it is the same
number twice. What agreement WOULD rule out is a transcription error in our linelist.
And DISagreement is decisive in the useful direction: it means ASD has been revised since
(Ruffoni-2014 FTS and friends were folded into ASD) and our GES entry is simply stale, in
which case the fix is a linelist update and not a tier decision at all.

Pulled verbatim by curl per the project's recipe: unit=0 (Angstrom -- unit=1 silently
returns the wrong range), format=3 (tab-delimited). bibrefs=1 returns HTTP 500 on this
endpoint, so the reference column is not available here and the circularity above cannot
be resolved from the response itself.

log gf = log10(g_i * f_ik) from the returned g_i and fik columns.
"""
from __future__ import annotations

import argparse
import csv
import io
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, os.environ.get("ISPEC_DIR", "/srv/codex/engines/ispec_src"))

ASD = "https://physics.nist.gov/cgi-bin/ASD/lines1.pl"


def asd_query(spectrum: str, lo: float, hi: float) -> pd.DataFrame:
    args = ["curl", "-sS", "-G", ASD,
            "--data-urlencode", f"spectra={spectrum}",
            "--data-urlencode", f"low_w={lo:.3f}",
            "--data-urlencode", f"upp_w={hi:.3f}",
            "--data-urlencode", "unit=0", "--data-urlencode", "de=0",
            "--data-urlencode", "format=3", "--data-urlencode", "line_out=0",
            "--data-urlencode", "remove_js=on", "--data-urlencode", "en_unit=1",
            "--data-urlencode", "output=0", "--data-urlencode", "page_size=15",
            "--data-urlencode", "show_obs_wl=1", "--data-urlencode", "A_out=1",
            "--data-urlencode", "f_out=on", "--data-urlencode", "allowed_out=1",
            "--data-urlencode", "g_out=on", "--data-urlencode", "enrg_out=on",
            "--data-urlencode", "J_out=on",
            "--data-urlencode", "submit=Retrieve Data"]
    out = subprocess.run(args, capture_output=True, text=True, timeout=120).stdout
    if not out.strip():
        return pd.DataFrame()
    rows = list(csv.DictReader(io.StringIO(out), delimiter="\t"))
    return pd.DataFrame(rows)


def _num(x):
    try:
        return float(str(x).strip().strip('"').replace("+", ""))
    except Exception:
        return np.nan


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sources", default="FMW,MRW,GESB82c")
    ap.add_argument("--lo", type=float, default=6910.0)
    ap.add_argument("--hi", type=float, default=9199.9)
    ap.add_argument("--max-per-source", type=int, default=14)
    ap.add_argument("--window", type=float, default=0.15,
                    help="+/- A around each line to ask ASD for")
    a = ap.parse_args()

    from pipeline.abundances_derive import _load_synth_resources
    ll, _, _ = _load_synth_resources()
    w = np.asarray(ll["wave_A"], float)
    el = np.asarray([str(x).strip() for x in ll["element"]])
    ref = np.asarray([str(x).strip() for x in ll["reference_code"]])
    gf = np.asarray(ll["loggf"], float)
    fe1 = np.array([e.upper().startswith("FE 1") for e in el])
    band = fe1 & (w >= a.lo) & (w <= a.hi)

    print("NOTE: FMW = Fuhr, Martin & Wiese, itself a NIST compilation. Agreement with")
    print("      ASD is therefore NOT independent confirmation; disagreement is the")
    print("      informative outcome (ASD revised, our GES entry stale).\n")

    rows = []
    for src in a.sources.split(","):
        m = band & np.array([str(r).split("+")[0] == src for r in ref])
        idx = np.where(m)[0][: a.max_per_source]
        print(f"  {src}: {int(m.sum())} in band, querying {len(idx)}")
        for i in idx:
            c = float(w[i])
            try:
                t = asd_query("Fe I", c - a.window, c + a.window)
            except Exception as e:
                print(f"    {c:.3f}: query failed {type(e).__name__}")
                continue
            if t.empty:
                rows.append(dict(source=src, wave=c, loggf_ges=gf[i],
                                 loggf_nist=np.nan, note="no ASD row"))
                continue
            cols = {k.strip().strip('"'): k for k in t.columns}
            wcol = next((cols[k] for k in cols if k.startswith("obs_wl")), None)
            fcol = next((cols[k] for k in cols if k == "fik"), None)
            gcol = next((cols[k] for k in cols if k in ("g_i", "gi")), None)
            if not (wcol and fcol and gcol):
                rows.append(dict(source=src, wave=c, loggf_ges=gf[i],
                                 loggf_nist=np.nan, note="columns missing"))
                continue
            t["_w"] = t[wcol].map(_num)
            t["_f"] = t[fcol].map(_num)
            t["_g"] = t[gcol].map(_num)
            t = t[t._w.notna() & t._f.notna() & t._g.notna() & (t._f > 0)]
            if t.empty:
                rows.append(dict(source=src, wave=c, loggf_ges=gf[i],
                                 loggf_nist=np.nan, note="no usable fik"))
                continue
            j = (t._w - c).abs().idxmin()
            lg = float(np.log10(t.loc[j, "_g"] * t.loc[j, "_f"]))
            rows.append(dict(source=src, wave=c, loggf_ges=float(gf[i]),
                             loggf_nist=lg, nist_wave=float(t.loc[j, "_w"]),
                             note=""))

    d = pd.DataFrame(rows)
    ok = d[d.loggf_nist.notna()].copy()
    print("\n" + "=" * 78)
    print("GES log gf vs NIST ASD, per source")
    print("=" * 78)
    if ok.empty:
        print("  no ASD matches returned — cannot adjudicate this way")
    else:
        ok["delta"] = ok.loggf_ges - ok.loggf_nist
        print(f"{'source':<10}{'matched':>8}{'median d(GES-NIST)':>20}{'max|d|':>9}")
        for s, g in ok.groupby("source"):
            print(f"{s:<10}{len(g):>8}{g.delta.median():>20.3f}{g.delta.abs().max():>9.3f}")
        print("\n  per line:")
        print(f"{'source':<10}{'wave':>10}{'GES':>9}{'NIST':>9}{'delta':>9}")
        for _, r in ok.sort_values(["source", "wave"]).iterrows():
            print(f"{r.source:<10}{r.wave:>10.3f}{r.loggf_ges:>9.3f}"
                  f"{r.loggf_nist:>9.3f}{r.delta:>9.3f}")
        print("\n  A gf too LOW by d makes the derived abundance too HIGH by ~d.")
        print("  So a NEGATIVE delta here would explain the positive A offset.")
    miss = d[d.loggf_nist.isna()]
    if len(miss):
        print(f"\n  {len(miss)} line(s) unmatched: "
              f"{miss.note.value_counts().to_dict()}")


if __name__ == "__main__":
    main()

"""The element's Asplund line set measured the way Asplund measures it (RYA-1350).

Asplund, Amarsi & Grevesse 2021 fit the DISK-CENTRE intensity (mu=1); every spectrum we fit
is the integrated flux, the geometry of a point-source target. RYA-1311 measured the
difference on Ca once, by hand (`scripts/diag_flux_vs_mu1.py`); this runs it for any element,
as an orchestrator step (`run_pipeline.py --element El --asplund-mu1`).

For each ion with a solar-literature `<EL>_AGSS21` set and each band the Liege atlas covers
(3000-10000 A), the production 1D-LTE synthesis cell (`derive_band_products`, pool
set:<EL>_AGSS21) runs on the Delbouille/Liege BASS2000 disk-centre atlas with Turbospectrum
intensity at mu=1 (RYA-441 injection, wrapping the fit and the local-continuum synthesis),
vsini 0, R 500000; the local continuum is referenced to our own synthesis (the IAG atlas is
flux). Writes data/results/orchestrator/budget/solar_<El>/asplund_mu1.json: per line our
mu=1 value, A21's 1D MARCS and 3D values where the set records them, and the mean.

    CODEX_SYNTH_CACHE=0 ISPEC_DIR=... python scripts/asplund_mu1.py Si
    (one cell, used by the above:)  ... python scripts/asplund_mu1.py --cell Si I VIS <out>
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]

LIEGE = "/mnt/codex-data/solar_reference/delbouille_bass2000/delbouille_bass2000_visible.csv.gz"
LIEGE_SPAN_A = (3000.0, 10000.0)
#: which production cell each band borrows (its line list, window, pool); only the spectrum
#: and the geometry change
CELL = {"near-UV": ("kpno_solar_atlas", "solar_kpno_molecfit_corrected"),
        "VIS": ("harps", "solar_harps_molecfit_corrected"),
        "red-optical": ("kpno_solar_atlas", "solar_kpno_molecfit_corrected"),
        "NIR": ("kpno_solar_atlas", "solar_kpno_molecfit_corrected")}


def _set(symbol: str, ion: str) -> dict | None:
    from pipeline import run_matrix as rm
    name = f"{symbol.upper()}_AGSS21"
    return next((r for r in rm.line_sets_for(symbol, ion) if r["set_name"] == name), None)


def plan(symbol: str) -> list[tuple[str, str, float, float]]:
    """(ion, band, lo, hi) for every band x ion holding an Asplund set line on Liege."""
    from pipeline import band_policy
    out = []
    for ion in ("I", "II"):
        s = _set(symbol, ion)
        if s is None:
            continue
        waves = s.get("coverage_waves") or s["waves"]
        for p in band_policy.POLICIES:
            if p.name not in CELL:
                continue
            lo, hi = max(p.lo_A, LIEGE_SPAN_A[0]), min(p.hi_A, LIEGE_SPAN_A[1])
            if lo < hi and any(lo <= w < hi for w in waves):
                out.append((ion, p.name, lo, hi))
    return out


def run_cell(symbol: str, ion: str, band: str, out: Path) -> None:
    import numpy as np
    import pandas as pd
    from pipeline import band_policy
    from pipeline.run_descriptor import RunDescriptor, resolve

    assert os.environ.get("CODEX_SYNTH_CACHE") == "0", \
        "the synth cache key ignores flux vs intensity: run with CODEX_SYNTH_CACHE=0"
    pol = next(p for p in band_policy.POLICIES if p.name == band)
    instrument, holding = CELL[band]
    d = RunDescriptor(element=symbol, ion=ion, instrument=instrument, holding=holding,
                      lo_A=max(pol.lo_A, LIEGE_SPAN_A[0]), hi_A=min(pol.hi_A, LIEGE_SPAN_A[1]),
                      engine_deck="ts-lte", method="synthesis", pool=f"set:{symbol.upper()}_AGSS21")
    step = [s for s in resolve(d, interpreter=sys.executable, ispec_dir=os.environ["ISPEC_DIR"]).steps
            if s["name"] == "derive_products"][0]
    args = list(step["args"])
    assert "--lines-from-set" in args, args
    args[args.index("--out") + 1] = str(out)

    import pipeline.nearuv_synth as ns
    _bsc = ns.build_solar_context

    def bsc(*a, **k):
        ctx = _bsc(*a, **k)
        ctx["vsini"] = 0.0                             # one point on the disk: no rotation
        return ctx
    ns.build_solar_context = bsc

    import pipeline.local_continuum as lc
    lc.solar_atlas_reference = lambda *a, **k: None    # the IAG atlas is flux, not mu=1

    import derive_band_products as dbp
    import measure_band_ew as mbe
    with gzip.open(LIEGE, "rt") as fh:
        L = pd.read_csv(fh).to_numpy(float)
    LW, LF = L[:, 0], L[:, 1]
    spec = mbe.HoldingSpec(holding_id="solar_delbouille_liege", reader="liege_mu1",
                           pre_normalised=True, span_A=LIEGE_SPAN_A)

    def lwx(instrument, centre, pad, segs=None, allow_uncorrected=False, holding=None):
        m = (LW >= centre - pad) & (LW <= centre + pad)
        if not m.any():
            raise LookupError(f"Liege serves nothing at {centre:.1f} A")
        return mbe.Window(LW[m], LF[m], "holding=solar_delbouille_liege BASS2000 mu=1", spec)
    mbe.select_holding = lambda *a, **k: spec
    mbe.load_window_ex = lwx

    _rc = pd.read_csv

    def rc(p, *a, **k):
        df = _rc(p, *a, **k)
        if str(p).endswith("instrument_catalog.csv"):
            df.loc[df.iloc[:, 0].astype(str) == instrument, "resolving_power_max"] = 500000
        return df
    pd.read_csv = rc

    from ace_co_mu1_bracket import turbospectrum_intensity_at_mu
    sys.argv = ["derive_band_products.py"] + args
    with turbospectrum_intensity_at_mu(1.0) as st:
        dbp.main()
    print(f"[asplund-mu1] intensity injections fired={st['fired']} remapped={st['remapped']}")
    assert st["remapped"] > 0, "no intensity remap -> this would be flux"
    del np


#: C/N/O are measured by pipeline.cno_synthesis (molecular bands + atomic cross-checks), not
#: by a line set; their mu=1 run is the VIS CNO region on Liege. A21 Table 3 gives one 1D
#: MARCS LTE value per INDICATOR; only the indicators A21 also uses on this spectrum pair
#: (CN red, and the near-UV OH/NH A-X, have no A21 counterpart: its CN/OH/NH are other systems).
CNO = ("C", "N", "O")
CNO_A21 = {"CH_Gband": "CH (A-X)", "CI_5052": "C I", "CI_5380": "C I",
           "C2_Swan": "C2 Swan", "OI_6300": "[O I]"}


def run_cno_cell(out: Path) -> None:
    """The VIS C/N/O region fitted to the Liege disk-centre atlas at mu=1 (one run serves
    C, N and O: the fit is coupled through the molecular equilibrium)."""
    import dataclasses
    import pipeline.cno_synthesis as cs
    import pipeline.local_continuum as lc
    from ace_co_mu1_bracket import turbospectrum_intensity_at_mu

    base = cs.REGIONS["vis"]
    reg = dataclasses.replace(base, name="vis_liege_mu1", instrument="delbouille_liege",
                              R=500000.0, holding="solar_delbouille_liege",
                              notes="RYA-1350 Asplund geometry: Liege disk centre, mu=1")
    cs.REGIONS[reg.name] = reg
    cs.REGION_DIAGNOSTICS[reg.name] = cs.REGION_DIAGNOSTICS["vis"]
    with gzip.open(LIEGE, "rt") as fh:
        import pandas as pd
        L = pd.read_csv(fh).to_numpy(float)
    m = (L[:, 0] >= reg.wave_min_A) & (L[:, 0] <= reg.wave_max_A)
    # the atlas is observed-frame (disk-centre gravitational redshift less convective
    # blueshift, ~+0.26 km/s on Fe I cores) and the CNO fit has no shift term: to rest first
    # (the RYA-478 rule the UVES arm follows)
    from pipeline.abundances_derive import _measure_rv_kms
    rv = _measure_rv_kms(L[m, 0] / 10.0, L[m, 1])
    print(f"[asplund-mu1] Liege VIS frame {rv:+.3f} km/s -> rest", flush=True)
    w_nm = L[m, 0] / 10.0 / (1.0 + rv / 299792.458)
    _load = cs._load_region_spectrum
    cs._load_region_spectrum = (lambda star, region: (w_nm, L[m, 1])
                                if region.name == reg.name else _load(star, region))
    _rb = cs._resolve_broadening

    def rb(star):
        R, vmac, _vsini, fit_vmac = _rb(star)
        return R, vmac, 0.0, fit_vmac                   # one point on the disk: no rotation
    cs._resolve_broadening = rb
    lc.solar_atlas_reference = lambda *a, **k: None     # the IAG atlas is flux, not mu=1

    out.mkdir(parents=True, exist_ok=True)
    with turbospectrum_intensity_at_mu(1.0) as st:
        res = cs.run_cno("solar", reg.name, out_dir=out, with_systematics=False,
                         tmp_dir=str(out / "ts_tmp"))
    print(f"[asplund-mu1] intensity injections fired={st['fired']} remapped={st['remapped']}")
    assert st["remapped"] > 0, "no intensity remap -> this would be flux"
    rows = [{k: (v if isinstance(v, (str, int, float, bool)) or v is None else str(v))
             for k, v in r.items()} for r in res.per_band]
    (out / "cno_mu1_per_band.json").write_text(json.dumps(
        {"abundances": res.abundances, "converged": res.converged, "flags": res.flags,
         "measured_elements": list(res.measured_elements), "per_band": rows},
        indent=1, default=str) + "\n")


def collect_cno(symbol: str, work: Path) -> dict:
    a21 = {}
    with (ROOT / "data/reference/agss21_cited_sets/agss21_table3_cno.csv").open(newline="") as fh:
        for r in csv.DictReader(fh):
            a21[r["indicator"]] = float(r["lte_MARCS"]) if r["lte_MARCS"] else None
    d = json.loads((work / "cno_VIS" / "cno_mu1_per_band.json").read_text())
    lines = []
    for r in d["per_band"]:
        if r.get("element") != symbol:
            continue
        a = r.get("A_X")
        ok = isinstance(a, (int, float)) and a == a and r.get("constrained", True) is not False
        ind = CNO_A21.get(r.get("key"))
        lines.append({"indicator": r.get("key"), "role": r.get("role"),
                      "A_mu1": a if ok else None, "red_chi2": r.get("red_chi2"),
                      "sigma_fit": r.get("sigma_fit"),
                      "status": r.get("status"), "A21_indicator": ind,
                      "A21_marcs_lte": a21.get(ind) if ind else None})
    used = [x["A_mu1"] for x in lines if x["A_mu1"] is not None]
    pl = [x for x in lines if x["A_mu1"] is not None and x["A21_marcs_lte"] is not None]
    paired = [x["A_mu1"] - x["A21_marcs_lte"] for x in pl]
    # one weak indicator can pass the constraint gate with a fit sigma of dex: the offset is
    # carried WITH that sigma, so it reads as what it is rather than as a disagreement
    sig = [x["sigma_fit"] for x in pl if isinstance(x["sigma_fit"], (int, float))]
    return {"element": symbol, "geometry": "disk-centre intensity mu=1 (Delbouille/Liege BASS2000)",
            "treatment": "1D-LTE", "route": "pipeline.cno_synthesis VIS region (coupled C/N/O)",
            "mean_minus_A21_sigma_fit": (round(sum(v * v for v in sig) ** 0.5 / len(sig), 3)
                                         if pl and len(sig) == len(pl) else None),
            "a21_basis": "AGSS21 Table 3 lte_MARCS, per indicator",
            "n_lines": len(used),
            "mean_A_mu1": round(sum(used) / len(used), 3) if used else None,
            "mean_minus_A21_marcs_lte": round(sum(paired) / len(paired), 3) if paired else None,
            "n_paired": len(paired), "lines": lines}


def main_cno(symbol: str) -> int:
    """C/N/O share one coupled run; it is reused if a sibling element already made it."""
    budget = ROOT / "data/results/orchestrator/budget"
    work = budget / "solar_CNO_asplund_mu1"
    rc = 0
    if not (work / "cno_VIS" / "cno_mu1_per_band.json").exists():
        env = {**os.environ, "CODEX_SYNTH_CACHE": "0"}
        print(f"[asplund-mu1] C/N/O VIS region on Liege (shared by C, N, O)", flush=True)
        rc = subprocess.run([sys.executable, __file__, "--cno-cell", str(work / "cno_VIS")],
                            env=env).returncode
        if rc:
            return rc
    res = collect_cno(symbol, work)
    res["cells"] = [{"ion": "C/N/O", "band": "VIS"}]
    out = budget / f"solar_{symbol}" / "asplund_mu1.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1) + "\n")
    print(f"[asplund-mu1] {symbol}: mu=1 mean {res['mean_A_mu1']} on {res['n_lines']} indicators; "
          f"vs A21 1D MARCS {res['mean_minus_A21_marcs_lte']} +/- {res['mean_minus_A21_sigma_fit']} "
          f"(fit sigma, n={res['n_paired']})")
    return rc


def collect(symbol: str, work: Path) -> dict:
    """Per-line mu=1 values from the cells' 1D-LTE line files, joined to A21's per-line
    values where the set records them."""
    lit = {}
    for f in (ROOT / "data/reference/agss21_cited_sets/agss21_cited_lines.csv",):
        with f.open(newline="") as fh:
            for r in csv.DictReader(fh):
                if r["species"].split()[0] == symbol:
                    lit[(r["species"], round(float(r["wavelength_air_A"]), 1))] = r
    lines = []
    for f in sorted(work.rglob("*_1D-LTE_lines.csv")):
        with f.open(newline="") as fh:
            for r in csv.DictReader(fh):
                sp = f"{r['element']} {r['ion']}"
                w = float(r["wavelength_air_A"])
                a = lit.get((sp, round(w, 1)), {})
                lines.append({"species": sp, "wavelength_air_A": w,
                              "A_mu1": float(r["abundance"]), "in_aggregate": r["in_aggregate"] == "True",
                              "red_chi2": float(r["red_chi2"]) if r.get("red_chi2") else None,
                              "A21_marcs_lte": float(a["A_marcs_lte"]) if a.get("A_marcs_lte") else None,
                              "A21_3d": float(a["A_3d_nlte_or_lte"]) if a.get("A_3d_nlte_or_lte") else None})
    used = [x["A_mu1"] for x in lines if x["in_aggregate"]]
    paired = [x["A_mu1"] - x["A21_marcs_lte"] for x in lines
              if x["in_aggregate"] and x["A21_marcs_lte"] is not None]
    return {"element": symbol, "geometry": "disk-centre intensity mu=1 (Delbouille/Liege BASS2000)",
            "treatment": "1D-LTE", "n_lines": len(used),
            "mean_A_mu1": round(sum(used) / len(used), 3) if used else None,
            "mean_minus_A21_marcs_lte": round(sum(paired) / len(paired), 3) if paired else None,
            "n_paired": len(paired), "lines": lines}


def main(argv: list[str]) -> int:
    if argv and argv[0] == "--cell":
        _, symbol, ion, band, out = argv
        run_cell(symbol, ion, band, Path(out))
        return 0
    if argv and argv[0] == "--cno-cell":
        run_cno_cell(Path(argv[1]))
        return 0
    symbol = argv[0]
    if symbol in CNO:
        return main_cno(symbol)
    work = ROOT / f"data/results/orchestrator/budget/solar_{symbol}/asplund_mu1"
    cells = plan(symbol)
    if not cells:
        print(f"[asplund-mu1] {symbol}: no Asplund set line on the Liege atlas; nothing to run")
        return 0
    env = {**os.environ, "CODEX_SYNTH_CACHE": "0"}
    rc = 0
    for ion, band, lo, hi in cells:
        out = work / f"{ion}_{band}"
        out.mkdir(parents=True, exist_ok=True)
        print(f"[asplund-mu1] {symbol} {ion} {band} {lo:.0f}-{hi:.0f} A", flush=True)
        r = subprocess.run([sys.executable, __file__, "--cell", symbol, ion, band, str(out)], env=env)
        rc = rc or r.returncode
    res = collect(symbol, work)
    res["cells"] = [{"ion": i, "band": b} for i, b, _, _ in cells]
    (work.parent / "asplund_mu1.json").write_text(json.dumps(res, indent=1) + "\n")
    print(f"[asplund-mu1] {symbol}: mu=1 mean {res['mean_A_mu1']} on {res['n_lines']} lines; "
          f"vs A21 1D MARCS per line {res['mean_minus_A21_marcs_lte']} (n={res['n_paired']})")
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

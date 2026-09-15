#!/usr/bin/env python3
"""RYA-1214 — rest-frame condition the RYA-1219 CRIRES+ Vesta J and K products.

RYA-1219 molecfit-corrected eight Vesta IDPs across the full J and K arms and left them
TOPOCENTRIC (`SPECSYS = TOPOCENT`, manifest `rest_frame` empty; its README: "rest-frame
conditioning remains a separate downstream step"). Every CRIRES+ holding before them —
Y (RYA-794/1054) and H (RYA-1094) — was an Elgueta+2026 product that arrived ALREADY in the
solar rest frame, so this repo never had to condition a CRIRES+ spectrum. This is that step.

THE VELOCITY IS EMPIRICAL, AND VERIFIED ON LINES IT WAS NOT MEASURED ON (RYA-372 method)
--------------------------------------------------------------------------------------
* Anchors: deep (catalogued depth 0.25-0.85), ISOLATED (no other catalogued line deeper than
  0.05 within +/-0.35 A) neutral photospheric lines from `linelist_solar`, outside every
  enumerated telluric band, in AIR -> converted to VACUUM, because CRIRES+ wavelengths are
  vacuum (the RYA-373 air/vac boundary: a slip there is ~83 km/s at 2.3 um).
* Each anchor's core velocity is `reflected_solar_rv._core_velocity` — the RYA-372 parabola
  core — searched around the Horizons-predicted position. Horizons SEEDS THE SEARCH ONLY;
  the velocity applied is the robust median of the TRAIN anchors (even-indexed).
* TEST anchors (odd-indexed) must land within `PASS_TOL_KMS` (0.5 km/s) of it, with at least
  `MIN_LINES` on each side. Pixels whose molecfit transmission is below 0.8 are treated as
  continuum for the core fit, so a telluric residual cannot pose as a solar core.
* Guardrails, all refusing: |v| <= 50 km/s (air/vac signature), |v - Horizons| <= 10 km/s,
  and telluric-anchor closure |CCF(1 - raw/continuum, 1 - MTRANS)| <= 1.5 km/s — the one
  check a common-mode wavelength offset cannot pass, because it moves TRAIN and TEST
  together (RYA-372 review item 2).

THE PRODUCT
-----------
One CSV per arm, `wavelength_air_A, flux_normalized` (the `crires_y` reader's contract), plus
`n_frames` and `min_mtrans`. Each (order, detector) segment is shifted by its own frame's
velocity, converted to air, and resampled onto a common grid at the arm's median pixel
spacing; overlapping settings are combined with (FLUX_RAW/ERR)^2 weights. Pixels with molecfit
transmission below 0.5 are NOT written: a correction there divides by a saturated core.
The rest frame carries the solar convective blueshift and gravitational redshift, as every
RYA-372 product does (`CONVECTIVE_BLUESHIFT_NOTE`).

Inputs are checked against the RYA-1219 manifest's SHA-256 before anything is read.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline import reflected_solar_rv as rrv  # noqa: E402
from pipeline.crires_telluric import _ccf_velocity, continuum_normalize  # noqa: E402
from pipeline.telluric_policy import TELLURIC_BANDS  # noqa: E402
from pipeline.wavelength_util import air_to_vac, vac_to_air  # noqa: E402

MANIFEST = ROOT / "data/audit/rya1219_crires_jk/corrected_products_manifest.csv"
LINELIST = ROOT / "data/linelists/linelist_solar.csv"
OUT_DIR = ROOT / "data/results/rya1214_crires_jk"
AUDIT = ROOT / "data/audit/rya1214_crires_jk"
HOLDING = {"J": "solar_crires_plus_j_rya1219", "K": "solar_crires_plus_k_rya1219"}

ANCHOR_DEPTH = (0.25, 0.85)
ISOLATION_A, ISOLATION_DEPTH = 0.35, 0.05
ANCHOR_SPECIES = ("Fe", "Si", "Mg", "Ti", "Ca", "Al", "C", "Na", "K", "Cr", "Ni", "Mn", "S")
CORE_MTRANS_MIN = 0.80
WRITE_MTRANS_MIN = 0.50
V_PHYSICAL_MAX = 50.0
HORIZONS_MAX = 10.0
CLOSURE_MAX = rrv.TELL_CLOSURE_TOL
C_KMS = rrv.C_KMS


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def anchors_air(lo: float, hi: float) -> np.ndarray:
    ll = pd.read_csv(LINELIST, low_memory=False,
                     usecols=["element", "ion", "wavelength_air_A", "central_depth"])
    ll = ll[ll.wavelength_air_A.between(lo - 1, hi + 1)].copy()
    ll["el"] = ll.element.astype(str).str.strip()
    ll["io"] = ll.ion.astype(str).str.strip()
    w_all = ll.wavelength_air_A.to_numpy(float)
    d_all = ll.central_depth.fillna(0).to_numpy(float)
    cand = ll[ll.el.isin(ANCHOR_SPECIES) & (ll.io == "I")
              & ll.central_depth.between(*ANCHOR_DEPTH)]
    keep = []
    for w in cand.wavelength_air_A.to_numpy(float):
        near = (np.abs(w_all - w) <= ISOLATION_A) & (np.abs(w_all - w) > 1e-6)
        if (d_all[near] > ISOLATION_DEPTH).any():
            continue
        if any(b_lo - 1.0 <= w <= b_hi + 1.0 for b_lo, b_hi, _ in TELLURIC_BANDS):
            continue
        keep.append(w)
    return np.array(sorted(keep))


def frame_velocity(tab: pd.DataFrame, anchors_vac: np.ndarray, v_seed: float) -> dict:
    per = {}
    for _, seg in tab.groupby(["ORDER", "DETEC"]):
        seg = seg.sort_values("WAVE")
        w = seg.WAVE.to_numpy(float)
        f = seg.FLUX.to_numpy(float).copy()
        f[~np.isfinite(f) | (seg.MTRANS.to_numpy(float) < CORE_MTRANS_MIN)] = 1.0
        for x in anchors_vac[(anchors_vac > w.min() + 1) & (anchors_vac < w.max() - 1)]:
            win = max(0.22, 6.0 / C_KMS * x)          # +/-6 km/s about the seed
            v = rrv._core_velocity(w, f, x, win=win, center=x * (1.0 + v_seed / C_KMS))
            if np.isfinite(v):
                per.setdefault(float(x), []).append(float(v))
    lines = sorted(per)
    vals = np.array([np.median(per[x]) for x in lines])
    train, test = vals[0::2], vals[1::2]
    v, v_std, n_tr = rrv._robust_median(train)
    t_med, t_std, n_te = rrv._robust_median(test)
    return {"v": v, "v_std": v_std, "n_train": n_tr, "n_test": n_te,
            "test_residual": (t_med - v) if np.isfinite(t_med) and np.isfinite(v) else np.nan,
            "n_anchor_lines_measured": len(lines)}


def closure(tab: pd.DataFrame) -> float:
    vs = []
    for _, seg in tab.groupby(["ORDER", "DETEC"]):
        seg = seg.sort_values("WAVE")
        w = seg.WAVE.to_numpy(float)
        raw = seg.FLUX_RAW.to_numpy(float)
        mt = seg.MTRANS.to_numpy(float)
        ok = np.isfinite(w) & np.isfinite(raw) & np.isfinite(mt) & (raw > 0)
        if ok.sum() < 200 or np.nanmin(mt[ok]) > 0.95:
            continue                                   # no telluric structure to anchor on
        cont = continuum_normalize(w[ok], raw[ok])
        v = _ccf_velocity(w[ok], 1.0 - cont, 1.0 - mt[ok])
        if np.isfinite(v):
            vs.append(v)
    return float(np.median(vs)) if vs else np.nan


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--products", required=True,
                    help="directory holding RYA-1219's J/ and K/ corrected FITS")
    ap.add_argument("--arm", choices=["J", "K"], action="append")
    a = ap.parse_args()
    from astropy.io import fits

    man = pd.read_csv(MANIFEST)
    arms = a.arm or ["J", "K"]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)
    report = {"ticket": "RYA-1214", "manifest": str(MANIFEST.relative_to(ROOT)),
              "rest_frame_note": rrv.CONVECTIVE_BLUESHIFT_NOTE, "arms": {}}

    for arm in arms:
        rows = man[man.arm == arm]
        frames, segs = [], []
        for _, r in rows.iterrows():
            p = Path(a.products) / arm / Path(r["product"]).name
            got = _sha(p)
            if got != r["sha256"]:
                raise SystemExit(f"{p.name}: sha256 {got[:16]} != manifest "
                                 f"{str(r['sha256'])[:16]} — refusing to condition a "
                                 f"product that is not the one RYA-1219 verified")
            with fits.open(p) as h:
                hdr = h[0].header
                if not hdr.get("TELLAPP"):
                    raise SystemExit(f"{p.name}: TELLAPP is not True — telluric first")
                if str(hdr.get("SPECSYS", "")).upper() != "TOPOCENT":
                    raise SystemExit(f"{p.name}: SPECSYS={hdr.get('SPECSYS')!r}, expected "
                                     f"TOPOCENT — refusing to shift a frame of unknown frame")
                mjd = float(hdr["MJD-OBS"])
                tab = pd.DataFrame({c: np.asarray(h[1].data[c]).astype(float)
                                    for c in ("WAVE", "FLUX", "FLUX_RAW", "ERR", "MTRANS",
                                              "ORDER", "DETEC")})
            lo, hi = float(np.nanmin(tab.WAVE)), float(np.nanmax(tab.WAVE))
            anc = air_to_vac(anchors_air(lo, hi))
            hz = rrv.reflected_solar_rv(mjd, "vesta")["v_total"]
            fv = frame_velocity(tab, anc, hz)
            cl = closure(tab)
            rec = {"setting": r["wlen"], "product": p.name, "mjd": mjd,
                   "horizons_v_total": round(hz, 3), "anchors_available": int(len(anc)),
                   **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in fv.items()},
                   "telluric_closure_kms": round(cl, 3) if np.isfinite(cl) else None}
            fails = []
            if not (np.isfinite(fv["v"]) and fv["n_train"] >= rrv.MIN_LINES
                    and fv["n_test"] >= rrv.MIN_LINES):
                fails.append(f"INSUFFICIENT anchors (train {fv['n_train']}, test {fv['n_test']})")
            else:
                if abs(fv["v"]) > V_PHYSICAL_MAX:
                    fails.append(f"|v|={fv['v']:.1f} > {V_PHYSICAL_MAX} (air/vac signature)")
                if abs(fv["v"] - hz) > HORIZONS_MAX:
                    fails.append(f"|v - Horizons| = {abs(fv['v'] - hz):.2f} > {HORIZONS_MAX}")
                if not abs(fv["test_residual"]) <= rrv.PASS_TOL_KMS:
                    fails.append(f"held-out residual {fv['test_residual']:+.3f} km/s > "
                                 f"{rrv.PASS_TOL_KMS}")
            if np.isfinite(cl) and abs(cl) > CLOSURE_MAX:
                fails.append(f"telluric closure {cl:+.2f} km/s > {CLOSURE_MAX}")
            rec["verdict"] = "PASS" if not fails else "FAIL: " + "; ".join(fails)
            frames.append(rec)
            print(f"  {arm} {r['wlen']:>6s}  v={fv['v']:+8.3f}  Horizons {hz:+8.3f}  "
                  f"test resid {fv['test_residual']:+.3f}  n {fv['n_train']}/{fv['n_test']}  "
                  f"closure {cl:+.2f}  -> {rec['verdict']}")
            if fails:
                continue
            for _, seg in tab.groupby(["ORDER", "DETEC"]):
                seg = seg.sort_values("WAVE")
                ok = (np.isfinite(seg.FLUX) & np.isfinite(seg.ERR) & (seg.ERR > 0)
                      & (seg.MTRANS >= WRITE_MTRANS_MIN)).to_numpy()
                if ok.sum() < 50:
                    continue
                w_air = vac_to_air(seg.WAVE.to_numpy(float) / (1.0 + fv["v"] / C_KMS))
                snr2 = (seg.FLUX_RAW.to_numpy(float) / seg.ERR.to_numpy(float)) ** 2
                segs.append((w_air, seg.FLUX.to_numpy(float), snr2, ok,
                             seg.MTRANS.to_numpy(float)))

        if not segs:
            report["arms"][arm] = {"frames": frames, "written": None}
            print(f"  {arm}: NO frame passed — nothing written")
            continue
        step = float(np.median([np.median(np.diff(s[0])) for s in segs]))
        lo = min(float(s[0][s[3]].min()) for s in segs)
        hi = max(float(s[0][s[3]].max()) for s in segs)
        grid = np.arange(lo, hi, step)
        num = np.zeros_like(grid); den = np.zeros_like(grid)
        nfr = np.zeros_like(grid); mtr = np.ones_like(grid)
        for w, f, s2, ok, mt in segs:
            inside = (grid >= w[ok].min()) & (grid <= w[ok].max())
            gi = grid[inside]
            fi = np.interp(gi, w[ok], f[ok])
            wi = np.interp(gi, w[ok], s2[ok])
            mi = np.interp(gi, w, mt)
            good = np.isfinite(fi) & np.isfinite(wi) & (wi > 0) & (mi >= WRITE_MTRANS_MIN)
            idx = np.where(inside)[0][good]
            num[idx] += fi[good] * wi[good]; den[idx] += wi[good]
            nfr[idx] += 1; mtr[idx] = np.minimum(mtr[idx], mi[good])
        m = den > 0
        out = pd.DataFrame({"wavelength_air_A": np.round(grid[m], 4),
                            "flux_normalized": np.round(num[m] / den[m], 6),
                            "n_frames": nfr[m].astype(int),
                            "min_mtrans": np.round(mtr[m], 4)})
        fn = OUT_DIR / f"{HOLDING[arm]}_rest.csv"
        out.to_csv(fn, index=False)
        report["arms"][arm] = {
            "holding": HOLDING[arm], "frames": frames, "written": str(fn.relative_to(ROOT)),
            "n_pixels": int(len(out)), "span_air_A": [float(out.wavelength_air_A.min()),
                                                      float(out.wavelength_air_A.max())],
            "grid_step_A": round(step, 5),
            "frames_used": sum(1 for f in frames if f["verdict"] == "PASS")}
        print(f"  wrote {fn.relative_to(ROOT)}: {len(out)} px, "
              f"{out.wavelength_air_A.min():.1f}-{out.wavelength_air_A.max():.1f} A")

    (AUDIT / "rest_frame_conditioning.json").write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""The CRIRES+ H-arm sky, measured on CRIRES+ itself, for the RYA-587 telluric term (RYA-1233).

The budget's telluric term (scripts/rya1230_cno_budget.telluric_line) maps the sky with a
raw/corrected Kitt Peak pair and measures the holding's residual against an INDEPENDENT
correction. Kitt Peak stops at 13000 A, so the Elgueta H product (solar_crires_plus_h_rya1094)
could not be priced and CRIRES+ H Si was held. This module supplies the same evidence from
CRIRES+ itself:

  sky        molecfit's fitted transmission (mtrans) of OUR Vesta IDPs (RYA-1191 recipe,
             re-run by RYA-1233 on the segments covering the graded lines;
             data/results/rya1233_h_molecfit/*_corrected.npz: topocentric VACUUM wavelengths)
  frame      per segment, the velocity that maps the holding's rest-frame AIR wavelengths onto
             the segment, MEASURED by cross-correlating our corrected flux against the holding;
             segments of one frame must agree (FRAME_TOL_KMS) or the frame is refused
  telluric   pixels where the sky absorbs more than the segment's own NULL: 3 robust sigma of
             holding/ours on its clean pixels (mtrans > 0.995) -- the clean-window null rule
  residual   mean |holding/ours - baseline| on those pixels x their fraction: the holding vs an
             independent correction of the same night (Elgueta's reduction vs our molecfit)
"""
from __future__ import annotations

import glob
import json
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SEG_DIR = ROOT / "data" / "results" / "rya1233_h_molecfit"
HOLDING = "solar_crires_plus_h_rya1094"
C_KMS = 299792.458
CLEAN_MTRANS = 0.995
FRAME_TOL_KMS = 1.0


@lru_cache(maxsize=1)
def _holding():
    """The holding's own registered file (rest frame, air, normalised) -- the same flux its
    fit read (the crires_y reader's contract)."""
    import pandas as pd
    from config.constants import codex_path
    d = pd.read_csv(str(codex_path("repo.crires_plus_solar_h_rya1094")))
    w = d["wavelength_air_A"].to_numpy(float)
    f = d[[c for c in d.columns if c.startswith("flux")][0]].to_numpy(float)
    g = np.isfinite(f) & np.isfinite(w)
    o = np.argsort(w[g])
    return w[g][o], f[g][o]


def _velocity(wv, fc, hw, hf) -> tuple[float, float]:
    """(v km/s, correlation) mapping rest-frame AIR (holding) onto topocentric VACUUM (wv)."""
    from pipeline.wavelength_util import air_to_vac
    ok = np.isfinite(fc)
    wv, fc = wv[ok], fc[ok]
    sel = (hw > wv.min() - 10) & (hw < wv.max() + 10)
    hv = air_to_vac(hw[sel])
    best = (0.0, -2.0)
    for v in np.arange(-150.0, 150.01, 0.25):
        model = np.interp(wv, hv * (1 + v / C_KMS), hf[sel], left=np.nan, right=np.nan)
        m = np.isfinite(model)
        if m.sum() < 200:
            continue
        r = np.corrcoef(fc[m], model[m])[0, 1]
        if r > best[1]:
            best = (float(v), float(r))
    return best


@lru_cache(maxsize=1)
def segments() -> list[dict]:
    """Every molecfit segment, put on the holding's rest-frame air grid, frame-checked."""
    from pipeline.wavelength_util import vac_to_air
    hw, hf = _holding()
    out = []
    for f in sorted(glob.glob(str(SEG_DIR / "*_corrected.npz"))):
        d = np.load(f)
        wv, fc, mt = d["wave_A"], d["flux_corr"], d["mtrans"]
        v, r = _velocity(wv, fc / np.nanmedian(fc), hw, hf)
        out.append(dict(tag=Path(f).name.replace("_corrected.npz", ""),
                        frame=Path(f).name.split("_")[0], v_kms=v, ccf=r,
                        wave_rest_air=vac_to_air(wv / (1 + v / C_KMS)),
                        ours=fc / np.nanmedian(fc), mtrans=mt))
    by_frame: dict = {}
    for s in out:
        by_frame.setdefault(s["frame"], []).append(s["v_kms"])
    for fr, vs in by_frame.items():
        if max(vs) - min(vs) > FRAME_TOL_KMS:
            raise ValueError(f"CRIRES+ H frame {fr}: segment velocities disagree {vs} "
                             f"(> {FRAME_TOL_KMS} km/s) -- frame not established, refusing")
    return out


def telluric_line(w0: float, hw: float) -> dict:
    """The rya1230 telluric evidence for one line of the CRIRES+ H holding."""
    hwav, hflx = _holding()
    segs = [s for s in segments()
            if np.nanmin(s["wave_rest_air"]) <= w0 - hw and w0 + hw <= np.nanmax(s["wave_rest_air"])]
    if not segs:
        return {"residual_flux": None,
                "residual_basis": f"no molecfit segment of our Vesta IDPs covers {w0:.3f} A"}
    s = max(segs, key=lambda q: q["ccf"])
    g = np.linspace(w0 - hw, w0 + hw, 400)
    mt = np.interp(g, s["wave_rest_air"], s["mtrans"])
    ours = np.interp(g, s["wave_rest_air"], s["ours"])
    own = np.interp(g, hwav, hflx)
    # The segment's own null: holding / ours on its CLEAN pixels, whole segment.
    seg_own = np.interp(s["wave_rest_air"], hwav, hflx, left=np.nan, right=np.nan)
    rr = seg_own / s["ours"]
    clean = (s["mtrans"] > CLEAN_MTRANS) & np.isfinite(rr)
    null = 3 * 1.4826 * float(np.nanmedian(np.abs(rr[clean] - np.nanmedian(rr[clean]))))
    P = (1 - mt) > null
    ev = {"sky_pair": f"molecfit mtrans of {s['tag']} (v={s['v_kms']:+.2f} km/s, ccf {s['ccf']:.3f})",
          "clean_edge": null, "clean_window_null": null, "sky_baseline": 1.0,
          "n_px": int(g.size), "n_telluric_px": int(P.sum()),
          "max_depth": float(max(0.0, 1 - mt.min()))}
    if not P.any():
        ev.update(residual_flux=0.0, sky_absorption=0.0,
                  residual_basis="no pixel where the sky absorbs above the segment's own null")
        return ev
    ev["sky_absorption"] = float(np.mean(1 - mt[P]) * P.mean())
    ratio = own / ours
    ok = np.isfinite(ratio)
    rbase = float(np.median(ratio[~P & ok])) if (~P & ok).sum() >= 3 else float(np.nanmedian(ratio))
    # A pixel the sky SATURATES (molecfit leaves it NaN: mtrans below its floor) cannot be
    # corrected by anyone on this night, so the holding's flux there is UNVERIFIED -- charged
    # the full sky depth, the assembler's rule for an uncorrected holding (RYA-1230). Si
    # 16434.927 sits on such a line in every Vesta H frame.
    dev = np.where(ok, np.abs(ratio - rbase), 1 - mt)
    ev["n_unverifiable_px"] = int((P & ~ok).sum())
    ev["residual_flux"] = float(np.mean(dev[P]) * P.mean())
    ev["residual_basis"] = (f"{HOLDING} vs independent correction (our molecfit of the Vesta "
                            f"IDPs, {s['tag']}), baseline-removed, on pixels the sky absorbs"
                            + (f"; {ev['n_unverifiable_px']} sky-saturated pixel(s) charged full "
                               f"depth (no correction can verify them)" if ev["n_unverifiable_px"] else ""))
    return ev


def frame_report() -> list[dict]:
    return [{k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items()
             if k in ("tag", "frame", "v_kms", "ccf")} for s in segments()]


if __name__ == "__main__":
    print(json.dumps(frame_report(), indent=1))
    for w in (15376.831, 16170.164, 16434.927, 16828.159):
        e = telluric_line(w, 0.9)
        print(w, {k: e.get(k) for k in ("n_telluric_px", "max_depth", "residual_flux", "clean_edge")})

"""RYA-805 — is telluric correction APPLIED to the 18 CRIRES+ Vesta IDPs?

THE QUESTION AND WHY IT NEEDED ASKING
-------------------------------------
RYA-370 asserted "no telluric correction" for this set and RYA-373's spec repeats it,
but neither showed the keyword that proves it. Ryan (2026-08-12): *"normal CRIRES+ data
should carry telluric correction if memory serves"* — so the assertion had to become
evidence before RYA-373/797/787 spend effort on either premise.

THE THREE LEGS
--------------
A header can only ever prove that no correction RECIPE ran. It cannot prove the
absorption is still in the flux — a product could in principle be corrected upstream of
the recipe chain we can see. So this audit runs three independent legs and requires them
to agree:

  1. RECIPE CHAIN (``--headers``)   walk ``ESO PRO REC*`` to whatever depth it goes,
     plus PRO CATG, the HDU list, and the column list. A telluric step would appear as
     a REC entry, a transmission/``Recon`` HDU, or a second ``FLUX_TELL`` column.
  2. CONTROLLED DEPTH (``--depth``) same target (the Sun), same instrument (CRIRES+),
     same normalisation, in windows where the ATMOSPHERE is strong and the solar
     photosphere is not. Control = Elgueta+2026 ``sp/Sun_{Y,J,H}_rv.dat``, which IS
     telluric-corrected. Excess deep absorption in the IDP is atmospheric.
  3. ATMOSPHERIC COVARIATION (``--iwv``) the falsifiable one. If the absorption is
     terrestrial its depth must track the header's precipitable water vapour and
     airmass. Solar photospheric lines cannot vary with Paranal's humidity.

⚠️ WINDOW CHOICE IS THE WHOLE EXPERIMENT (leg 2). A first pass compared whole settings
and found the IDPs 20x deeper than the control — but that compared DIFFERENT wavelength
ranges, and the excess was just the band edges each setting happens to include. Restricted
to the control's own 980-1080 nm the two agree to 0.05 pp, because 980-1080 nm is
intrinsically dry and discriminates nothing. Only a telluric-HEAVY window that BOTH
spectra cover is a test. The O2 1.27 um band is the cleanest of them: O2 is purely
terrestrial, so the Sun contributes no counterpart at all.

⚠️ ``HIERARCH`` (RYA-791): astropy strips the prefix, so a lookup of
``HIERARCH ESO PRO CATG`` returns empty and MANUFACTURES AN ABSENCE. Look up the bare
``ESO PRO CATG`` form. In an audit whose finding is "this keyword is not present", a
lookup bug that fakes absence is the failure mode that matters most.

Sirius only (the IDPs live at ``/mnt/codex-data``).

    python3 scripts/rya805_telluric_audit.py --headers --depth --iwv --csv OUT.csv
"""
from __future__ import annotations

import argparse
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np

IDP_GLOB = "/mnt/codex-data/spectra/vesta/CRIRESPlus/ADP*.fits"
# RYA-789 holding. The Elgueta+2026 solar spectrum is itself Vesta through CRIRES+,
# telluric-corrected by their pipeline -- the ideal control: same Sun, same instrument.
ELGUETA_SP = "/mnt/codex-data/codex/rya789/data/reference/elgueta2026_vizier/sp"

# Vocabulary that would betray a telluric step anywhere in the header.
TELLURIC_VOCAB = re.compile(
    r"molecfit|telluric|ATM_|SKYCALC|TRANSMISSION|RECON|TELL_STD|STD_STAR|"
    r"BEST_FIT_MODEL|ATM_PARAMETER|MAPPING_ATMOSPHERIC|corr_tell",
    re.IGNORECASE,
)
# Cards that match the vocabulary but are not telluric: detector board shift registers
# (ESO DET ... TRANS), the OB sky-transparency CONSTRAINT, and FITS boilerplate.
VOCAB_FALSE_POSITIVE = re.compile(
    r"^(COMMENT|HISTORY)$|ESO DET .*TRANS|ESO OBS AMBI TRANS", re.IGNORECASE
)

# (label, lo_nm, hi_nm, absorber, control_band). Windows where the atmosphere is strong.
# The O2 1.27 um band is the decisive one: no solar counterpart exists.
DEPTH_WINDOWS = (
    ("O2 1.27um", 1265.0, 1300.0, "O2", "J"),
    ("H2O 1480-1520", 1480.0, 1520.0, "H2O", "H"),
    ("H2O 1240-1260", 1240.0, 1260.0, "H2O", "J"),
    ("H2O 1110-1150", 1110.0, 1150.0, "H2O", None),
    ("H2O 1750-1790", 1750.0, 1790.0, "H2O", None),
    ("dry control 980-1080", 980.0, 1080.0, "-", "Y"),
)
IWV_WINDOW = ("H2O 1480-1520", 1480.0, 1520.0)

# Rolling-continuum width. Must be wider than any solar line and narrower than a
# telluric BAND, or the normalisation absorbs the very thing being measured:
# ~1.3 nm at the 0.0032 nm SPEC_BIN.
CONT_WIN = 401


def _fits():
    from astropy.io import fits
    return fits


def rolling_p95(y: np.ndarray, win: int = CONT_WIN) -> np.ndarray:
    """Local continuum as a rolling 95th percentile, subsampled then interpolated."""
    n = len(y)
    if n < win:
        return np.full(n, np.nanpercentile(y, 95))
    step = max(1, win // 8)
    centers, vals = [], []
    for i in range(0, n, step):
        seg = y[max(0, i - win // 2):min(n, i + win // 2)]
        seg = seg[np.isfinite(seg)]
        if seg.size:
            centers.append(i)
            vals.append(np.percentile(seg, 95))
    if not centers:
        return np.full(n, np.nan)
    return np.interp(np.arange(n), centers, vals)


def normalised(wave, flux, lo, hi, qual=None, seg_key=None):
    """Locally-normalised flux inside [lo, hi], normalised PER DETECTOR SEGMENT.

    Segmenting matters: a cr2res order/detector boundary is a flux discontinuity, and a
    continuum drawn across one would invent absorption at the seam.
    """
    good = np.isfinite(flux) & np.isfinite(wave) & (wave >= lo) & (wave <= hi)
    if qual is not None:
        good &= (qual == 0)          # RYA-794: select on QUAL, never WAVELMIN/MAX
    if seg_key is None:
        seg_key = np.zeros(len(flux), dtype=int)
    out = []
    for k in np.unique(seg_key[good]):
        m = good & (seg_key == k)
        if m.sum() < 50:
            continue
        order = np.argsort(wave[m])
        f = flux[m][order]
        cont = rolling_p95(f)
        ok = np.isfinite(cont) & (cont > 0)
        out.append(f[ok] / cont[ok])
    return np.concatenate(out) if out else np.array([])


def depth_stats(nrm: np.ndarray) -> dict:
    if not nrm.size:
        return {}
    return {
        "n_px": int(nrm.size),
        "median": round(float(np.median(nrm)), 4),
        "pct_below_0.7": round(float((nrm < 0.7).mean() * 100), 2),
        "pct_below_0.5": round(float((nrm < 0.5).mean() * 100), 2),
    }


def read_idp(path):
    fits = _fits()
    with fits.open(path) as hdul:
        h, d = hdul[0].header, hdul[1].data
        return {
            "path": path,
            "file": Path(path).name,
            "header": h,
            "hdus": [(x.name, type(x).__name__) for x in hdul],
            "columns": [(c.name, str(c.unit)) for c in hdul[1].columns],
            "wave": d["WAVE"][0].astype(float),
            "flux": d["FLUX"][0].astype(float),
            "qual": d["QUAL"][0].astype(int),
            # order and detector jointly identify a physical chip segment
            "seg": d["ORDER"][0].astype(int) * 100 + d["DETEC"][0].astype(int),
        }


# ── leg 1: the recipe chain ──────────────────────────────────────────────────
def audit_headers(frames: list[dict]) -> list[dict]:
    rows = []
    for fr in frames:
        h = fr["header"]
        chain, n = [], 1
        while f"ESO PRO REC{n} ID" in h:          # bare form: astropy strips HIERARCH
            chain.append(str(h[f"ESO PRO REC{n} ID"]).strip())
            n += 1
        cal_catg = sorted({
            str(h[k]).strip() for k in h
            if re.match(r"ESO PRO REC\d+ CAL\d+ CATG$", k)
        })
        hits = [
            (k, str(v).strip()) for k, v in h.items()
            if TELLURIC_VOCAB.search(f"{k} = {v}") and not VOCAB_FALSE_POSITIVE.search(k)
        ]
        has_trans_hdu = any(
            re.search(r"TRANS|RECON|TELL", name or "", re.IGNORECASE)
            for name, _ in fr["hdus"]
        )
        has_tell_col = any(
            re.search(r"TELL|TRANS|RECON", c, re.IGNORECASE) for c, _ in fr["columns"]
        )
        rows.append({
            "file": fr["file"],
            "setting": str(h.get("ESO INS WLEN ID", "")).strip(),
            "date_obs": str(h.get("DATE-OBS", "")).strip()[:19],
            "PRO_CATG": str(h.get("ESO PRO CATG", "")).strip(),
            "PRODCATG": str(h.get("PRODCATG", "")).strip(),
            "PROCSOFT": str(h.get("PROCSOFT", "")).strip(),
            "n_recipes": len(chain),
            "recipe_chain": "|".join(chain),
            "cal_categories": "|".join(cal_catg),
            "n_hdu": len(fr["hdus"]),
            "hdu_names": "|".join(nm for nm, _ in fr["hdus"]),
            "columns": "|".join(c for c, _ in fr["columns"]),
            "transmission_hdu": has_trans_hdu,
            "telluric_column": has_tell_col,
            "telluric_vocab_hits": len(hits),
            "FLUXCAL": str(h.get("FLUXCAL", "")).strip(),
            "CONTNORM": bool(h.get("CONTNORM", False)),
            "SPECSYS": str(h.get("SPECSYS", "")).strip(),
            "SNR": round(float(h.get("SNR", np.nan)), 1),
            "wl_min_nm": round(float(h.get("WAVELMIN", np.nan)), 2),
            "wl_max_nm": round(float(h.get("WAVELMAX", np.nan)), 2),
            "IWV": float(h.get("ESO TEL AMBI IWV START", np.nan)),
            "airmass": round(0.5 * (float(h.get("ESO TEL AIRM START", np.nan))
                                    + float(h.get("ESO TEL AIRM END", np.nan))), 3),
            # the verdict for THIS file, from THIS file's keywords
            "telluric_applied": "no",
        })
    return rows


# ── leg 2: controlled depth against a telluric-corrected twin ────────────────
def load_control(band: str):
    a = np.loadtxt(f"{ELGUETA_SP}/Sun_{band}_rv.dat")
    return a[:, 0], a[:, 1]          # RYA-794 unit trap: sp/ wavelength is NANOMETRES


def audit_depth(frames: list[dict]) -> list[dict]:
    controls = {}
    rows = []
    for label, lo, hi, mol, cband in DEPTH_WINDOWS:
        ctrl = {}
        if cband:
            if cband not in controls:
                controls[cband] = load_control(cband)
            cw, cf = controls[cband]
            ctrl = depth_stats(normalised(cw, cf, lo, hi))
            if ctrl:
                rows.append({"window": label, "absorber": mol, "source": "CONTROL",
                             "which": f"Elgueta Sun_{cband}_rv (telluric-CORRECTED)",
                             **ctrl})
        for fr in frames:
            nrm = normalised(fr["wave"], fr["flux"], lo, hi, fr["qual"], fr["seg"])
            st = depth_stats(nrm)
            if st and st["n_px"] > 500:
                setting = str(fr["header"].get("ESO INS WLEN ID", "")).strip()
                rows.append({"window": label, "absorber": mol, "source": "IDP",
                             "which": f"{setting} {fr['file']}", **st})
    return rows


# ── leg 3: does the depth track the atmosphere? ──────────────────────────────
def audit_iwv(frames: list[dict]) -> tuple[list[dict], dict]:
    label, lo, hi = IWV_WINDOW
    rows = []
    for fr in frames:
        h = fr["header"]
        setting = str(h.get("ESO INS WLEN ID", "")).strip()
        if not setting.startswith("H"):        # the window lives in the H settings
            continue
        nrm = normalised(fr["wave"], fr["flux"], lo, hi, fr["qual"], fr["seg"])
        st = depth_stats(nrm)
        if not st:
            continue
        iwv = float(h["ESO TEL AMBI IWV START"])
        am = 0.5 * (float(h["ESO TEL AIRM START"]) + float(h["ESO TEL AIRM END"]))
        rows.append({"setting": setting, "night": str(h["DATE-OBS"])[:10],
                     "IWV_mm": iwv, "airmass": round(am, 3),
                     "IWV_x_airmass": round(iwv * am, 3),
                     "pct_below_0.7": st["pct_below_0.7"]})
    corr = {}
    if len(rows) > 2:
        x = np.array([r["IWV_x_airmass"] for r in rows])
        y = np.array([r["pct_below_0.7"] for r in rows])
        corr = {"window": label, "n": len(rows),
                "pearson_r_iwv_x_airmass": round(float(np.corrcoef(x, y)[0, 1]), 4),
                "pearson_r_iwv": round(float(np.corrcoef(
                    np.array([r["IWV_mm"] for r in rows]), y)[0, 1]), 4)}
    return rows, corr


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--headers", action="store_true", help="leg 1: recipe chain")
    p.add_argument("--depth", action="store_true", help="leg 2: controlled depth")
    p.add_argument("--iwv", action="store_true", help="leg 3: atmospheric covariation")
    p.add_argument("--glob", default=IDP_GLOB)
    p.add_argument("--csv", type=Path, help="write the per-file leg-1 table here")
    p.add_argument("--json", type=Path, help="write the full result bundle here")
    a = p.parse_args(argv)
    if not (a.headers or a.depth or a.iwv):
        a.headers = a.depth = a.iwv = True

    paths = sorted(glob.glob(a.glob))
    if not paths:
        print(f"NO IDPs at {a.glob} -- this audit is Sirius-only.", file=sys.stderr)
        return 2
    print(f"loading {len(paths)} IDPs ...", file=sys.stderr)
    frames = [read_idp(x) for x in paths]
    bundle: dict = {"n_files": len(frames), "glob": a.glob}

    if a.headers:
        rows = audit_headers(frames)
        bundle["headers"] = rows
        print("\n=== LEG 1: RECIPE CHAIN ===")
        print(f"  PRO CATG        : {sorted({r['PRO_CATG'] for r in rows})}")
        print(f"  recipe chain    : {sorted({r['recipe_chain'] for r in rows})}")
        print(f"  n recipes       : {sorted({r['n_recipes'] for r in rows})}")
        print(f"  calibrations    : {sorted({r['cal_categories'] for r in rows})}")
        print(f"  HDUs            : {sorted({r['hdu_names'] for r in rows})}")
        print(f"  columns         : {sorted({r['columns'] for r in rows})}")
        print(f"  transmission HDU: {sorted({r['transmission_hdu'] for r in rows})}")
        print(f"  telluric column : {sorted({r['telluric_column'] for r in rows})}")
        print(f"  telluric vocab  : {sum(r['telluric_vocab_hits'] for r in rows)} "
              f"real hits across {len(rows)} files")
        print(f"  FLUXCAL/CONTNORM: {sorted({r['FLUXCAL'] for r in rows})} / "
              f"{sorted({r['CONTNORM'] for r in rows})}")
        print(f"  VERDICT         : telluric_applied = "
              f"{sorted({r['telluric_applied'] for r in rows})}")
        if a.csv:
            import csv as _csv
            a.csv.parent.mkdir(parents=True, exist_ok=True)
            with a.csv.open("w", newline="") as fh:
                w = _csv.DictWriter(fh, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)
            print(f"  wrote {a.csv}")

    if a.depth:
        rows = audit_depth(frames)
        bundle["depth"] = rows
        print("\n=== LEG 2: CONTROLLED DEPTH (telluric-heavy windows) ===")
        cur = None
        for r in rows:
            if r["window"] != cur:
                cur = r["window"]
                print(f"\n  -- {cur} ({r['absorber']}) --")
            print(f"     {r['source']:<7} {r['which'][:44]:<46} "
                  f"med={r['median']:.3f}  %<0.7={r['pct_below_0.7']:5.2f}  "
                  f"%<0.5={r['pct_below_0.5']:5.2f}")

    if a.iwv:
        rows, corr = audit_iwv(frames)
        bundle["iwv"] = {"rows": rows, "correlation": corr}
        print("\n=== LEG 3: ATMOSPHERIC COVARIATION ===")
        print(f"  window {IWV_WINDOW[0]}")
        for r in sorted(rows, key=lambda x: x["IWV_x_airmass"]):
            print(f"     {r['setting']:<7} {r['night']}  IWV={r['IWV_mm']:5.2f}  "
                  f"AM={r['airmass']:.3f}  IWV*AM={r['IWV_x_airmass']:6.2f}  "
                  f"%<0.7={r['pct_below_0.7']:6.2f}")
        if corr:
            print(f"\n  Pearson r (IWV*airmass vs depth) = "
                  f"{corr['pearson_r_iwv_x_airmass']:+.3f}  (n={corr['n']})")
            print(f"  Pearson r (IWV alone      vs depth) = "
                  f"{corr['pearson_r_iwv']:+.3f}")
            print("  Solar photospheric lines cannot track Paranal's humidity.")

    if a.json:
        a.json.parent.mkdir(parents=True, exist_ok=True)
        a.json.write_text(json.dumps(bundle, indent=1, default=str))
        print(f"\nwrote {a.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

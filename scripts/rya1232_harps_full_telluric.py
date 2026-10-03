#!/usr/bin/env python3
"""RYA-1232 -- full-range telluric correction of the HARPS direct-Sun exposures.

    python scripts/rya1232_harps_full_telluric.py <exposures_dir> <gdas> <iag_atlas> \
        <workroot> <corrected_dir> [--only NAME]

🔴 WHY THIS EXISTS. RYA-931 corrected HARPS for O2 B ONLY. Against the RYA-1230 full-coverage
Kitt Peak transmission, `solar_harps_molecfit_corrected` was byte-identical to raw at
H2O 5850-6030 (sky depth 0.26), O2 gamma + H2O 6250-6630 (0.52) and the weaker H2O bands
-- every HARPS VIS product (Fe, [O I] 6300, C I) was measured through uncorrected sky.

🔴 WHY NOT ONE BIG molecfit FIT. Adding the H2O windows to the molecfit_model fit failed on
every LSF start (RYA-1232, exposure .735): the windows share ONE Gaussian kernel, the
solar residuals in the H2O windows (rms 6-9%) outweigh their weak telluric lines, and the
kernel railed to its 100 px bound, taking the O2 column to 0.25-0.45 and wrecking the
O2 B fit that works on its own. FIT_RES_GAUSS=FALSE cannot hold the kernel (it disables
convolution, RYA-931).

THE STANDARD ESO TWO-STEP INSTEAD:
  1. molecfit_model on the proven O2 B window ONLY (RYA-931's fit, its LSF ladder and
     admissibility gates unchanged), with H2O listed at the night's GDAS column (held);
  2. molecfit_calctrans evaluates that solution -- kernel, wavelength solution, O2 column
     -- over the FULL exposure, 3782-6913 A;
  3. the H2O column is then MEASURED, not assumed: calctrans is run again with H2O
     removed, giving T_H2O = T_all / T_O2, and the scale s in T_O2 * T_H2O**s is fitted
     against the IAG telluric-free solar atlas (Baker+2020, broadened to the fitted kernel)
     over the H2O windows, with a free linear continuum per window. The final
     transmission is re-evaluated by calctrans at H2O x s.

Output is RYA-931's format exactly (FLUX divided, MTRANS + QUARMASK extensions, TELL*
header keys), so `rya931_normalize_corrected.py` builds the holding unchanged.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import rya931_correct_exposures as C  # noqa: E402
import rya931_molecfit_model as M  # noqa: E402

O2B_WINDOW = "6857-6911"
#: H2O-dominated windows where the column scale is fitted (KP-measured sky depth 0.06-0.28),
#: placed between the Na D lines and clear of the H-alpha core.
H2O_FIT_WINDOWS_A = ((5900.0, 5985.0), (6465.0, 6535.0), (5695.0, 5790.0))
#: a telluric-free window for the HARPS-vs-atlas registration (KP sky depth < 0.5%)
REGISTRATION_WINDOW_A = (6100.0, 6200.0)
#: the atlas is only trusted where it is near continuum: deep solar cores carry the
#: atlas-vs-HARPS line-shape mismatch, not the atmosphere
ATLAS_MIN_FLUX = 0.90
#: transmission this close to 1 is left untouched (no division, no count)
UNITY = 1.0 - 1.0e-5


def _esorex_env():
    esorex = M.resolve_esorex()
    bindir = str(Path(esorex).parent)
    env = dict(os.environ, PATH=f"{bindir}:{os.environ.get('PATH', '')}")
    return esorex, env


def stage_a(source: Path, gdas: Path, atlas: Path, workdir: Path, starts,
            molecules: str = "O2:fit,H2O:1.0", prefix: str = "stageA"):
    """RYA-931's O2 B fit, H2O held at the GDAS column. Same ladder, same gates.
    Returns (dir, params, attempts); dir is None when no start is admissible."""
    env = dict(os.environ, RYA931_INCLUDE_A=O2B_WINDOW, RYA931_MOLECULES=molecules)
    attempts = []
    for start in starts:
        d = workdir / f"{prefix}_{start:.3f}".replace(".", "p")
        if not (d / "out" / "BEST_FIT_MODEL.fits").exists():
            r = subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "rya931_molecfit_model.py"),
                 str(source), str(gdas), str(d), "--solar-atlas", str(atlas),
                 "--frame", "AIR_RV", "--rv-sign", "-1", "--lsf-init-fwhm-pix", f"{start}"],
                text=True, capture_output=True, env=env)
            if r.returncode:
                raise SystemExit(f"{source.name}: stage A molecfit_model failed\n{r.stdout}\n{r.stderr}")
        p = {str(x["parameter"]): float(x["value"])
             for x in fits.open(d / "out" / "BEST_FIT_PARAMETERS.fits")[1].data}
        ok = (p["gaussfwhm"] >= C.ACCEPT_MIN_LSF_FWHM_PIX
              and p["reduced_chi2"] <= C.ACCEPT_MAX_REDUCED_CHI2
              and abs(p["rel_mol_col_O2"] - 1.0) <= C.ACCEPT_MAX_O2_COLUMN_DEVIATION)
        attempts.append({"lsf_start_pix": start, "gaussfwhm": p["gaussfwhm"],
                         "reduced_chi2": p["reduced_chi2"],
                         "rel_mol_col_O2": p["rel_mol_col_O2"], "admissible": ok})
        if ok:
            return d, p, attempts
    return None, None, attempts


def full_science(source: Path, dest: Path) -> tuple[np.ndarray, np.ndarray]:
    """The whole exposure as a molecfit SCIENCE table (RYA-931's writer, full extent)."""
    old = os.environ.get("RYA931_INCLUDE_A")
    os.environ["RYA931_INCLUDE_A"] = "3000-8000"
    try:
        M.write_inputs(source, None, dest)
    finally:
        if old is None:
            os.environ.pop("RYA931_INCLUDE_A", None)
        else:
            os.environ["RYA931_INCLUDE_A"] = old
    t = fits.open(dest / "science.fits")[1].data
    return np.asarray(t["lambda"], float) * 1.0e4, np.asarray(t["flux"], float)


def _with_donor_h2o(stage_out: Path, donor_out: Path, run_dir: Path) -> None:
    """O2-only fallback: add the H2O column of a same-night, same-GDAS exposure's profile
    (interpolated on height) and an H2O row to MODEL_MOLECULES. The H2O AMOUNT is still
    fitted against the atlas afterwards; only the GDAS profile SHAPE is borrowed."""
    with fits.open(stage_out / "ATM_PARAMETERS.fits") as h, fits.open(donor_out / "ATM_PARAMETERS.fits") as dn:
        t = h[1].data
        dd = dn[1].data
        h2o = np.interp(np.asarray(t["HGT"], float), np.asarray(dd["HGT"], float), np.asarray(dd["H2O"], float))
        cols = list(h[1].columns) + [fits.Column(name="H2O", format="D", array=h2o)]
        hdus = [h[0].copy(), fits.BinTableHDU.from_columns(cols, header=h[1].header, name=h[1].name)]
        hdus += [x.copy() for x in h[2:]]
        fits.HDUList(hdus).writeto(run_dir / "ATM_PARAMETERS.fits", overwrite=True)
    with fits.open(stage_out / "MODEL_MOLECULES.fits") as h:
        m = h[1].data
        names = [str(x).strip() for x in m["LIST_MOLEC"]] + ["H2O"]
        fit = [int(x) for x in m["FIT_MOLEC"]] + [0]
        rel = [float(x) for x in m["REL_COL"]] + [1.0]
        fits.HDUList([h[0].copy(), fits.BinTableHDU.from_columns([
            fits.Column(name="LIST_MOLEC", format="4A", array=names),
            fits.Column(name="FIT_MOLEC", format="1J", array=np.array(fit, dtype=np.int32)),
            fits.Column(name="REL_COL", format="1D", array=rel)], header=h[1].header)]
        ).writeto(run_dir / "MODEL_MOLECULES.fits", overwrite=True)


def calctrans(stage_out: Path, science_dir: Path, run_dir: Path, h2o_scale: float,
              donor_out: Path | None = None) -> np.ndarray:
    """molecfit_calctrans over the full SCIENCE table with the stage-A solution and the
    H2O profile scaled by `h2o_scale` (0 => O2 only). Returns mtrans row-for-row."""
    if (run_dir / "out" / "TELLURIC_DATA.fits").exists():
        return np.asarray(fits.open(run_dir / "out" / "TELLURIC_DATA.fits")[1].data["mtrans"], float)
    (run_dir / "out").mkdir(parents=True, exist_ok=True)
    shutil.copy(science_dir / "science.fits", run_dir / "science.fits")
    if donor_out is None:
        shutil.copy(stage_out / "MODEL_MOLECULES.fits", run_dir / "MODEL_MOLECULES.fits")
    # 🔴 The wavelength correction fitted on the 54 A O2 B chip is a polynomial in THAT
    # chip's normalised coordinate; applied to the 3100 A full-exposure chip it shifted the
    # model ~0.95 A at 5284 A and decorrelated every H2O line (fitted scale 0.005). HARPS
    # wavelengths are good to m/s and the barycentric shift is carried by the AIR_RV frame,
    # so calctrans uses the IDENTITY correction.
    with fits.open(stage_out / "BEST_FIT_PARAMETERS.fits") as h:
        h = fits.HDUList([x.copy() for x in h])
        t = h[1].data
        for i, name in enumerate(t["parameter"]):
            if str(name).strip() == "chip 1, coef 0":
                t["value"][i] = 0.0
            elif str(name).strip() == "chip 1, coef 1":
                t["value"][i] = 1.0
        h.writeto(run_dir / "BEST_FIT_PARAMETERS.fits", overwrite=True)
    if donor_out is not None:
        _with_donor_h2o(stage_out, donor_out, run_dir)
    with fits.open(run_dir / "ATM_PARAMETERS.fits" if donor_out is not None
                   else stage_out / "ATM_PARAMETERS.fits") as h:
        h = fits.HDUList([x.copy() for x in h])
        for ext in h[1:]:
            if ext.data is not None and "H2O" in ext.columns.names:
                ext.data["H2O"] = ext.data["H2O"] * max(h2o_scale, 1.0e-6)
        h.writeto(run_dir / "ATM_PARAMETERS.fits", overwrite=True)
    (run_dir / "calctrans.sof").write_text(
        "science.fits SCIENCE\nMODEL_MOLECULES.fits MODEL_MOLECULES\n"
        "ATM_PARAMETERS.fits ATM_PARAMETERS\nBEST_FIT_PARAMETERS.fits BEST_FIT_PARAMETERS\n")
    esorex, env = _esorex_env()
    cmd = [esorex, f"--output-dir={(run_dir / 'out').resolve()}", "--suppress-prefix=TRUE",
           "molecfit_calctrans", "--USE_ONLY_INPUT_PRIMARY_DATA=FALSE",
           "--MAPPING_ATMOSPHERIC=0,1", "--MAPPING_CONVOLVE=0,1", "--CHIP_EXTENSIONS=FALSE",
           "calctrans.sof"]
    r = subprocess.run(cmd, cwd=run_dir, env=env, text=True, capture_output=True)
    (run_dir / "esorex.stdout.txt").write_text(r.stdout)
    (run_dir / "esorex.stderr.txt").write_text(r.stderr)
    (run_dir / "command.json").write_text(json.dumps({"cmd": cmd, "h2o_scale": h2o_scale}) + "\n")
    td = run_dir / "out" / "TELLURIC_DATA.fits"
    if r.returncode or not td.exists():
        raise SystemExit(f"calctrans failed in {run_dir} (rc={r.returncode}); see esorex.stdout.txt")
    return np.asarray(fits.open(td)[1].data["mtrans"], float)


def atlas_on(wave_a: np.ndarray, atlas: Path, fwhm_a: float) -> np.ndarray:
    with fits.open(atlas) as h:
        d = h[1].data
        aw = M.vac_to_air(1.0e8 / np.asarray(d["v"], float))
        af = np.asarray(d["s"], float)
    o = np.argsort(aw)
    aw, af = aw[o], af[o]
    step = float(np.median(np.diff(aw)))
    af = gaussian_filter1d(af, fwhm_a / 2.3548 / step)
    return np.interp(wave_a, aw, af, left=np.nan, right=np.nan)


def register(wave, flux, atlas_fn) -> float:
    """HARPS-vs-atlas shift (A) on a telluric-free window, by chi2 over a shift grid."""
    lo, hi = REGISTRATION_WINDOW_A
    k = (wave >= lo) & (wave <= hi) & np.isfinite(flux)
    best = None
    for sh in np.arange(-0.10, 0.1001, 0.002):
        a = atlas_fn(wave[k] - sh)
        g = np.isfinite(a)
        A = np.vstack([a[g], a[g] * (wave[k][g] - lo)]).T
        c, *_ = np.linalg.lstsq(A, flux[k][g], rcond=None)
        r = float(np.sum((flux[k][g] - A @ c) ** 2))
        if best is None or r < best[0]:
            best = (r, sh)
    return float(best[1])


def fit_h2o_scale(wave, flux, atl, t_o2, t_h2o) -> tuple[float, float, int]:
    """ln(flux/atlas) - ln T_O2 = per-window linear continuum + s * ln T_H2O."""
    rows, ys, cols = [], [], []
    nwin = len(H2O_FIT_WINDOWS_A)
    for i, (lo, hi) in enumerate(H2O_FIT_WINDOWS_A):
        k = ((wave >= lo) & (wave <= hi) & np.isfinite(flux) & (flux > 0) & np.isfinite(atl)
             & (atl > ATLAS_MIN_FLUX) & (t_o2 > 0.5) & (t_h2o > 0.3))
        if k.sum() < 50:
            continue
        x = np.log(t_h2o[k])
        y = np.log(flux[k]) - np.log(atl[k]) - np.log(t_o2[k])
        X = np.zeros((k.sum(), 2 * nwin + 1))
        X[:, 2 * i] = 1.0
        X[:, 2 * i + 1] = (wave[k] - lo) / (hi - lo)
        X[:, -1] = x
        rows.append(X)
        ys.append(y)
    X, y = np.vstack(rows), np.concatenate(ys)
    used = np.any(X != 0, axis=0)
    used[-1] = True
    X = X[:, used]
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ c
    cov = np.linalg.inv(X.T @ X) * float(np.var(res, ddof=X.shape[1]))
    return float(c[-1]), float(np.sqrt(cov[-1, -1])), int(len(y))


def correct(source: Path, gdas: Path, atlas: Path, workdir: Path, corrected_dir: Path,
            accepted_start: float | None) -> dict:
    # RYA-931's ladder, then finer steps: holding H2O in the molecule list moves the
    # collapsed-LSF basin slightly, and .745 landed in it from all seven RYA-931 starts.
    extra = (6.65, 6.85, 6.9, 6.7, 6.6, 7.1, 6.4, 6.95, 6.8)
    starts = ([accepted_start] if accepted_start else []) + [
        s for s in (*C.LSF_START_LADDER_PIX, *extra) if s != accepted_start]
    sa_dir, params, attempts = stage_a(source, gdas, atlas, workdir, starts)
    donor = None
    if sa_dir is None:
        # RYA-931's exact O2-only fit (admissible for every exposure in RYA-931), with the
        # H2O profile borrowed from an exposure of the same night and GDAS profile.
        sa_dir, params, more = stage_a(source, gdas, atlas, workdir, starts,
                                       molecules="O2:fit", prefix="stageAo2")
        attempts += more
        if sa_dir is None:
            raise SystemExit(f"{source.name}: no admissible O2 B solution; attempts={json.dumps(attempts)}")
        cands = sorted(p for p in workdir.parent.glob("*/stageA_*/out/ATM_PARAMETERS.fits")
                       if "H2O" in fits.open(p)[1].columns.names)
        if not cands:
            raise SystemExit(f"{source.name}: O2-only fallback needs a same-night H2O profile donor")
        donor = cands[0].parent
    sa_out = sa_dir / "out"
    sci = workdir / "science_full"
    wave, flux_sci = full_science(source, sci)

    t_all = calctrans(sa_out, sci, workdir / "ct_h2o_1", 1.0, donor)
    t_o2 = calctrans(sa_out, sci, workdir / "ct_h2o_0", 0.0, donor)
    t_h2o = np.where(t_o2 > 0, t_all / np.where(t_o2 > 0, t_o2, 1.0), 1.0)

    pix_a = float(np.median(np.diff(wave)))
    fwhm_a = params["gaussfwhm"] * pix_a

    def atlas_fn(w):
        return atlas_on(w, atlas, fwhm_a)
    shift = register(wave, flux_sci, atlas_fn)
    atl = atlas_fn(wave - shift)
    s, s_err, n_fit = fit_h2o_scale(wave, flux_sci, atl, t_o2, t_h2o)
    if not (0.2 <= s <= 5.0):
        raise SystemExit(f"{source.name}: fitted H2O scale {s:.3f} is unphysical; refusing")
    trans_full = calctrans(sa_out, sci, workdir / f"ct_h2o_{s:.4f}".replace(".", "p"), s, donor)

    # dT and T_min exactly as RYA-931 measures them, on the stage-A O2 B model
    model = fits.open(sa_out / "BEST_FIT_MODEL.fits")[1].data
    with fits.open(source) as hdul:
        header = hdul[0].header.copy()
        table = hdul[1].data
        w0 = np.asarray(table["WAVE"][0], float)
        flux = np.asarray(table["FLUX"][0], float)
        o2_keep = ((w0 >= 6857.0 - M.EXTRACT_MARGIN_A) & (w0 <= 6911.0 + M.EXTRACT_MARGIN_A)
                   & np.isfinite(flux))
        if int(o2_keep.sum()) != len(model):
            raise SystemExit(f"{source.name}: stage-A model rows do not map onto the exposure")
        dt = C.measure_dt(model, w0[o2_keep], C.solar_mask_from(sa_out, w0[o2_keep]))
        t_min = dt / C.CORRECTED_RELATIVE_ERROR_BUDGET

        keep = (w0 >= 3000.0 - M.EXTRACT_MARGIN_A) & (w0 <= 8000.0 + M.EXTRACT_MARGIN_A) & np.isfinite(flux)
        if int(keep.sum()) != trans_full.size or not np.allclose(w0[keep], wave):
            raise SystemExit(f"{source.name}: calctrans rows do not map onto the exposure; refusing")
        transmission = np.ones_like(w0)
        transmission[keep] = np.where(trans_full > 0, trans_full, 1.0)
        # Inside the O2 B fit window keep RYA-931's own fitted transmission (its fitted
        # wavelength solution is right THERE); calctrans with the identity correction left
        # a -0.109 mean residual on O2 B sky pixels against -0.014 for the fitted model.
        m_t = np.asarray(model["mtrans"], float)
        o2_idx = np.where(o2_keep)[0]
        transmission[o2_idx[m_t > 0]] = m_t[m_t > 0]
        active = transmission < UNITY
        applied = active & (transmission >= t_min)
        quarantined = active & (transmission < t_min)
        corrected = flux.copy()
        corrected[applied] = flux[applied] / transmission[applied]
        corrected[quarantined] = np.nan

        header["HISTORY"] = "RYA-1232 full-range telluric correction: molecfit_model O2 B + molecfit_calctrans"
        header["HISTORY"] = f"RYA-1232 H2O column scale {s:.4f} +/- {s_err:.4f} (fit vs IAG atlas); T_min={t_min:.4f}"
        header["HISTORY"] = f"RYA-1232 GDAS {gdas.name} md5 {C.md5(gdas)}"
        header["HISTORY"] = f"RYA-1232 source {source.name} md5 {C.md5(source)}"
        header["TELLCORR"] = (True, "RYA-1232 full-range telluric correction applied")
        header["TELLENG"] = ("molecfit-4.4.4", "telluric correction engine")
        header["TELLTMIN"] = (t_min, "min transmission corrected; below = NaN")
        header["TELLNPX"] = (int(applied.sum()), "pixels divided by transmission")
        header["TELLNQAR"] = (int(quarantined.sum()), "pixels quarantined (NaN)")
        header["TELLH2OS"] = (s, "H2O column scale vs GDAS, fitted")
        out = corrected_dir / source.name
        table = table.copy()
        table["FLUX"][0] = corrected
        mt = fits.ImageHDU(data=transmission.astype(np.float64), name="MTRANS")
        mt.header["TELLENG"] = ("molecfit-4.4.4", "telluric correction engine")
        mt.header["TELLTMIN"] = (t_min, "min transmission corrected")
        mt.header["COMMENT"] = "transmission divided out of FLUX; 1.0 = untouched"
        qh = fits.ImageHDU(data=quarantined.astype(np.uint8), name="QUARMASK")
        qh.header["COMMENT"] = "1 = transmission below TELLTMIN; FLUX set NaN"
        fits.HDUList([fits.PrimaryHDU(header=header),
                      fits.BinTableHDU(data=table, header=hdul[1].header.copy()),
                      mt, qh]).writeto(out, overwrite=True)

    def depth(lo, hi):
        k = (w0 >= lo) & (w0 <= hi)
        return float(1.0 - np.nanmin(transmission[k])) if k.any() else float("nan")
    return {
        "exposure": source.name, "source_md5": C.md5(source), "corrected_md5": C.md5(out),
        "mjd_obs": float(header["MJD-OBS"]),
        "airmass_start": float(header.get("HIERARCH ESO TEL AIRM START", np.nan)),
        "rel_mol_col_O2": params["rel_mol_col_O2"], "reduced_chi2_o2b": params["reduced_chi2"],
        "lsf_gauss_fwhm_pix": params["gaussfwhm"], "lsf_start_attempts": attempts,
        "h2o_profile_donor": str(donor) if donor else None,
        "atlas_shift_A": shift, "h2o_scale": s, "h2o_scale_err": s_err, "h2o_fit_pixels": n_fit,
        "dt_transmission_sigma": dt, "t_min": t_min,
        "n_corrected": int(applied.sum()), "n_quarantined": int(quarantined.sum()),
        "max_depth": {"H2O_5900": depth(5850, 6030), "O2g_6280": depth(6270, 6300),
                      "H2O_6500": depth(6440, 6560), "O2B_6870": depth(6866, 6912),
                      "blue_4400_5000": depth(4400, 5000)},
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("exposures", type=Path)
    ap.add_argument("gdas", type=Path)
    ap.add_argument("atlas", type=Path)
    ap.add_argument("workroot", type=Path)
    ap.add_argument("corrected_dir", type=Path)
    ap.add_argument("--only", default=None)
    ap.add_argument("--rya931-report", type=Path,
                    default=ROOT / "data/audit/rya931_molecfit_runtime/exposure_correction.json")
    a = ap.parse_args()
    a.corrected_dir.mkdir(parents=True, exist_ok=True)
    accepted = {}
    if a.rya931_report.exists():
        accepted = {r["exposure"]: r.get("lsf_start_accepted_pix")
                    for r in json.loads(a.rya931_report.read_text())}
    srcs = sorted(a.exposures.glob("ADP*.fits"))
    if a.only:
        srcs = [s for s in srcs if a.only in s.name]
    report = []
    for s in srcs:
        r = correct(s, a.gdas.resolve(), a.atlas.resolve(),
                    a.workroot / s.stem.replace(".", "_"), a.corrected_dir, accepted.get(s.name))
        print(json.dumps({k: v for k, v in r.items() if k != "lsf_start_attempts"}), flush=True)
        report.append(r)
    out = a.workroot / ("rya1232_exposure_correction.json" if not a.only else f"rya1232_{a.only}.json")
    out.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()

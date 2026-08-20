#!/usr/bin/env python3
"""Telluric-correct every HARPS direct-solar exposure at the exposure level (RYA-931).

Runs `molecfit_model` per exposure via scripts/rya931_molecfit_model.py, then
divides the O2 B band out of each ADP with that exposure's own fitted
transmission.  Originals are opened read-only and never rewritten.

Two properties of the molecfit product drive the arithmetic here:

* `mtrans` is 0 outside the fitted WAVE_INCLUDE region, not 1.  Dividing the
  whole spectrum by it would zero everything blueward of 6857 A.  The
  correction is therefore applied only inside the fitted region, which is also
  the only place a registered telluric band falls for HARPS (O2 B 6867-6884 is
  the sole entry of pipeline.telluric_policy.TELLURIC_BANDS inside 3782-6912 A).
* Deep cores cannot be corrected.  Dividing by a transmission T amplifies the
  model's own transmission error dT into a relative flux error dT/T, so pixels
  below a measured threshold are quarantined (set NaN), never forced into a
  corrected state.  The threshold is derived per exposure from that exposure's
  measured dT, not hard-coded.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

O2_B_LO_A, O2_B_HI_A = 6867.0, 6884.0
EXTRACT_LO_A, EXTRACT_HI_A = 6840.0, 6912.0
CLEAN_CONTROL_A = (6800.0, 6860.0)
# Budget: a corrected pixel may carry at most this relative flux error from the
# transmission model.  Quarantine threshold T_min = dT_measured / BUDGET.
CORRECTED_RELATIVE_ERROR_BUDGET = 0.05

# The molecfit fit for this band is bistable.  Alongside the physical solution
# (Gaussian LSF ~6.77 px, reduced chi2 ~0.81, O2 column ~1.00x the reference) sits
# a degenerate one that converges just as happily with a ZERO-width LSF, reduced
# chi2 ~11.5 and the O2 column 3.7% off.  Which basin an exposure lands in is
# chaotic in the starting LSF width -- 6.767 px and 6.773 px send the same
# exposure to different answers -- so a single starting value cannot be trusted.
#
# These are physical admissibility criteria, declared before the fits are run,
# not targets to tune toward.  A line-spread function cannot have zero width, and
# O2 is a well-mixed gas whose column is known to far better than 10%: a fit that
# violates either is refuted by physics regardless of how cleanly it converged.
# Only the STARTING POINT is varied on retry; no fit parameter, threshold or
# result is adjusted.
LSF_START_LADDER_PIX = (6.773, 6.500, 7.000, 6.250, 7.250, 6.000, 7.500)
ACCEPT_MIN_LSF_FWHM_PIX = 1.0
ACCEPT_MAX_REDUCED_CHI2 = 2.0
ACCEPT_MAX_O2_COLUMN_DEVIATION = 0.10


def md5(path: Path) -> str:
    d = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            d.update(chunk)
    return d.hexdigest()


def measure_dt(model, air_lam: np.ndarray, solar_mask: np.ndarray) -> float:
    """Robust transmission-domain accuracy on solar-clean, in-band pixels.

    `air_lam` is passed in because BEST_FIT_MODEL's own `lambda` column is in
    vacuum; selecting the O2 B window on it would shift the window by ~1.9 A.
    """
    band = ((air_lam >= O2_B_LO_A) & (air_lam <= O2_B_HI_A)
            & (model["mscal"] > 0) & ~solar_mask)
    residual = (model["flux"][band] - model["mflux"][band]) / model["mscal"][band]
    return float(1.4826 * np.median(np.abs(residual - np.median(residual))))


def solar_mask_from(products: Path, lam_a: np.ndarray) -> np.ndarray:
    path = products / "WAVE_EXCLUDE.fits"
    mask = np.zeros(lam_a.size, dtype=bool)
    if not path.exists():
        return mask
    for row in fits.open(path)[1].data:
        mask |= (lam_a >= row[0] * 1.0e4) & (lam_a <= row[1] * 1.0e4)
    return mask


def correct_one(source: Path, gdas: Path, atlas: Path, workdir: Path,
                corrected_dir: Path, lsf_fwhm_pix: float | None = None) -> dict:
    ladder = [lsf_fwhm_pix] if lsf_fwhm_pix else list(LSF_START_LADDER_PIX)
    attempts, products = [], None
    for start in ladder:
        attempt_dir = workdir / f"start_{start:.3f}".replace(".", "p")
        if not (attempt_dir / "out" / "BEST_FIT_MODEL.fits").exists():
            run = subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "rya931_molecfit_model.py"),
                 str(source), str(gdas), str(attempt_dir), "--solar-atlas", str(atlas),
                 "--frame", "AIR_RV", "--rv-sign", "-1",
                 "--lsf-init-fwhm-pix", f"{start}"],
                text=True, capture_output=True)
            if run.returncode:
                raise SystemExit(
                    f"{source.name}: molecfit_model failed\n{run.stdout}\n{run.stderr}")
        p = {str(r["parameter"]): float(r["value"])
             for r in fits.open(attempt_dir / "out" / "BEST_FIT_PARAMETERS.fits")[1].data}
        admissible = (p["gaussfwhm"] >= ACCEPT_MIN_LSF_FWHM_PIX
                      and p["reduced_chi2"] <= ACCEPT_MAX_REDUCED_CHI2
                      and abs(p["rel_mol_col_O2"] - 1.0) <= ACCEPT_MAX_O2_COLUMN_DEVIATION)
        attempts.append({"lsf_start_pix": start, "gaussfwhm": p["gaussfwhm"],
                         "reduced_chi2": p["reduced_chi2"],
                         "rel_mol_col_O2": p["rel_mol_col_O2"], "admissible": admissible})
        if admissible:
            products = attempt_dir / "out"
            break
    if products is None:
        raise SystemExit(
            f"{source.name}: no physically admissible molecfit solution from any of "
            f"{len(ladder)} starting points; attempts={json.dumps(attempts)}")
    model = fits.open(products / "BEST_FIT_MODEL.fits")[1].data
    params = {str(r["parameter"]): float(r["value"])
              for r in fits.open(products / "BEST_FIT_PARAMETERS.fits")[1].data}
    trans = np.asarray(model["mtrans"], dtype=float)
    fitted = trans > 0.0            # mtrans == 0 marks "outside the fitted region"

    with fits.open(source) as hdul:
        header = hdul[0].header.copy()
        table = hdul[1].data
        wave = np.asarray(table["WAVE"][0], dtype=float)
        flux = np.asarray(table["FLUX"][0], dtype=float)
        err = np.asarray(table["ERR"][0], dtype=float)

        # Map the transmission back by ROW INDEX, not by wavelength.  molecfit
        # reports BEST_FIT_MODEL's `lambda` in VACUUM (the column is our declared
        # air grid times the refractive index, ~1.000277 here), so matching it
        # against the air ADP grid silently mis-registers the correction by
        # ~1.9 A -- about 190 HARPS pixels.  The output table is row-for-row the
        # SCIENCE table we wrote, which the flux column proves exactly, so the
        # `keep` mask that built it is the correct and exact mapping.
        keep = (wave >= EXTRACT_LO_A) & (wave <= EXTRACT_HI_A) & np.isfinite(flux)
        air_lam = wave[keep]
        dt = measure_dt(model, air_lam, solar_mask_from(products, air_lam))
        t_min = dt / CORRECTED_RELATIVE_ERROR_BUDGET
        if int(keep.sum()) != trans.size:
            raise SystemExit(f"{source.name}: model has {trans.size} rows but the "
                             f"extract mask selects {int(keep.sum())}; refusing to guess")
        if not np.array_equal(np.asarray(model["flux"]), flux[keep]):
            raise SystemExit(f"{source.name}: BEST_FIT_MODEL flux is not row-identical "
                             "to the extracted flux; the row mapping is not safe")
        transmission = np.ones_like(wave)
        transmission[keep] = np.where(fitted, trans, 1.0)
        in_fit = np.zeros(wave.size, dtype=bool)
        in_fit[keep] = fitted
        applied = in_fit & (transmission >= t_min)
        quarantined = in_fit & (transmission < t_min)

        corrected = flux.copy()
        corrected[applied] = flux[applied] / transmission[applied]
        corrected[quarantined] = np.nan

        header["HISTORY"] = "RYA-931 exposure-level telluric correction (ESO molecfit 4.4.4)"
        header["HISTORY"] = f"RYA-931 transmission from BEST_FIT_MODEL mtrans; T_min={t_min:.4f}"
        header["HISTORY"] = f"RYA-931 GDAS {gdas.name} md5 {md5(gdas)}"
        header["HISTORY"] = f"RYA-931 source {source.name} md5 {md5(source)}"
        header["TELLCORR"] = (True, "RYA-931 telluric correction applied")
        header["TELLENG"] = ("molecfit-4.4.4", "telluric correction engine")
        header["TELLTMIN"] = (t_min, "min transmission corrected; below = NaN")
        header["TELLNPX"] = (int(applied.sum()), "pixels divided by transmission")
        header["TELLNQAR"] = (int(quarantined.sum()), "pixels quarantined (NaN)")

        out = corrected_dir / source.name
        table = table.copy()
        table["FLUX"][0] = corrected
        # Persist the transmission that was divided out, as its own extension.
        # This is what makes the product self-describing: pipeline.telluric_intake
        # reads a non-unity transmission extension as positive evidence of
        # APPLIED, so the holding's telluric state is derived from the product
        # rather than asserted in a registry column (RYA-806/927 C4).  It also
        # lets a downstream line-level check query the local transmission, which
        # is what the ratified two-level policy requires.
        # An IMAGE extension, not a bintable: pipeline.telluric_intake can read
        # an image's values and so reach its strongest verdict ("a transmission
        # extension is present and is NOT all-unity").  From a bintable it can
        # only establish presence, and the module is explicit that presence is
        # NOT proof -- values decide.  Row order matches extension 1 exactly.
        mtrans_hdu = fits.ImageHDU(data=transmission.astype(np.float64), name="MTRANS")
        mtrans_hdu.header["TELLENG"] = ("molecfit-4.4.4", "telluric correction engine")
        mtrans_hdu.header["TELLTMIN"] = (t_min, "min transmission corrected")
        mtrans_hdu.header["COMMENT"] = "transmission divided out of FLUX; 1.0 = untouched"
        quar_hdu = fits.ImageHDU(data=quarantined.astype(np.uint8), name="QUARMASK")
        quar_hdu.header["COMMENT"] = "1 = transmission below TELLTMIN; FLUX set NaN"
        fits.HDUList([fits.PrimaryHDU(header=header),
                      fits.BinTableHDU(data=table, header=hdul[1].header.copy()),
                      mtrans_hdu, quar_hdu]).writeto(out, overwrite=True)

    def stats(values, lo, hi, continuum):
        k = (wave >= lo) & (wave <= hi) & np.isfinite(values)
        v = values[k] / continuum
        return {"n": int(k.sum()), "min": float(v.min()), "median": float(np.median(v)),
                "pct_below_0.5": float(100 * np.mean(v < 0.5)),
                "pct_below_0.95": float(100 * np.mean(v < 0.95))}

    cont = float(np.percentile(flux[(wave >= 6840) & (wave <= 6866)], 95))
    return {
        "exposure": source.name, "source_md5": md5(source), "corrected_md5": md5(out),
        "mjd_obs": float(header["MJD-OBS"]), "airmass_start": float(header.get("HIERARCH ESO TEL AIRM START", np.nan)),
        "rel_mol_col_O2": params["rel_mol_col_O2"], "ppmv_O2": params["ppmv_O2"],
        "reduced_chi2": params["reduced_chi2"], "rms_rel_to_mean": params["rms_rel_to_mean"],
        "lsf_gauss_fwhm_pix": params["gaussfwhm"],
        "lsf_start_attempts": attempts,
        "lsf_start_accepted_pix": attempts[-1]["lsf_start_pix"],
        "dt_transmission_sigma": dt, "t_min": t_min,
        "n_corrected": int(applied.sum()), "n_quarantined": int(quarantined.sum()),
        "o2b_before": stats(flux, O2_B_LO_A, O2_B_HI_A, cont),
        "o2b_after": stats(corrected, O2_B_LO_A, O2_B_HI_A, cont),
        "control_before": stats(flux, *CLEAN_CONTROL_A, cont),
        "control_after": stats(corrected, *CLEAN_CONTROL_A, cont),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("exposures", type=Path)
    ap.add_argument("gdas", type=Path)
    ap.add_argument("atlas", type=Path)
    ap.add_argument("workroot", type=Path)
    ap.add_argument("corrected_dir", type=Path)
    ap.add_argument("--lsf-init-fwhm-pix", type=float, default=None)
    args = ap.parse_args()

    args.corrected_dir.mkdir(parents=True, exist_ok=True)
    sources = sorted(args.exposures.glob("ADP*.fits"))
    if not sources:
        raise SystemExit(f"no ADP*.fits under {args.exposures}")
    report = [correct_one(s, args.gdas.resolve(), args.atlas.resolve(),
                          args.workroot / s.stem.replace(".", "_"), args.corrected_dir,
                          args.lsf_init_fwhm_pix)
              for s in sources]
    out = args.workroot / "rya931_exposure_correction.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

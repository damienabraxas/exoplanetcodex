#!/usr/bin/env python3
"""Telluric-correct one band of the 1984 Kitt Peak solar flux atlas (RYA-940).

The 1984 atlas differs from the HARPS case (RYA-931) in one way that governs the
whole design: **it carries no observation metadata**. No airmass, no date, no
site keywords — the only non-numeric text in all 251 segments was a saved HTTP
500 page. So the ratified observation-night-GDAS rule cannot be met, and Ryan
ratified the alternative on 2026-08-21: fit the molecular columns freely against
a standard atmosphere and referee the result against a telluric-free reference
rather than against provenance.

That is sound for correction specifically, because airmass and column density are
degenerate in the transmission: molecfit only has to reproduce the absorption we
can see, and `rel_mol_col` absorbs whatever path length produced it. It does mean
the fitted column is NOT an atmospheric measurement here, and it is therefore NOT
usable as the physical referee that RYA-931 leaned on for HARPS. What replaces it
is external agreement, per region, named explicitly by `--referee`.

Three atlas properties a naive reader gets wrong, all from the recovered README:

* wavelength is AIR and column 2 is a PSEUDO-residual flux — unity IS the
  continuum, so nothing here should be re-normalised;
* segments carry 4 nm PLUS 0.1 nm of overlap at the red end, so concatenating
  them naively double-counts every seam;
* the grid is quasi-even (step spread 7-15% inside one segment). molecfit's
  kernel convolution assumes regular sampling, so the fit runs on a resampled
  even grid — and the transmission is then applied back on the ORIGINAL grid, so
  the delivered product keeps the atlas's own sampling. Resampling is a fitting
  device, never a change to the product.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.rya931_molecfit_model import resolve_esorex  # noqa: E402
from scripts.rya938_kp_crosscheck import read_kurucz2005, read_iag  # noqa: E402

#: McMath/Pierce Solar Telescope, Kitt Peak, Arizona (recovered 1984 README).
SITE = {"latitude_deg": 31.9584, "longitude_deg": -111.5967, "elevation_m": 2096.0}

#: Physical admissibility, declared before any fit runs (RYA-931's discipline).
#: NOTE the difference from HARPS: the O2-column-equals-one referee is NOT
#: available here, because with airmass unknown the column absorbs it. What
#: remains physical is that a line-spread function cannot have zero width and a
#: model must actually describe the data.
ACCEPT_MIN_LSF_FWHM_PIX = 0.8
#: NOT reduced chi2. This atlas ships no uncertainties, so chi2 is computed against
#: molecfit's ASSUMED DEFAULT_ERROR -- gating on it gates on our own assumption, and
#: the same fit can be "rejected" or "accepted" purely by changing that constant.
#: `rms_rel_to_mean` is the fractional residual: error-independent, directly
#: interpretable, and it separates the cases cleanly in practice (a runaway fit gives
#: 0.43, the physical one 0.05).
ACCEPT_MAX_RMS_REL_TO_MEAN = 0.10
#: With airmass unknown it is absorbed into the fitted column, so RYA-931's "a
#: well-mixed gas must return 1.0" referee is NOT available here. What survives is
#: weaker but still physical: the column must correspond to a POSSIBLE airmass.
#: Below 1 sec z is impossible; above 5 is implausible for a solar atlas.
ACCEPT_MIN_REL_COL = 0.5
ACCEPT_MAX_REL_COL = 5.0
#: The atlas resolving power, MEASURED: the O2 A fit returns a Gaussian FWHM of
#: 6.040 px at 0.00850 A/px and 7640 A, i.e. R = 7640 / (6.040 x 0.00850) = 148,000.
#: R is an INSTRUMENT property; the FWHM in PIXELS is not, because the atlas
#: sampling runs from 0.0067 to 0.0183 A across 2960-13000 A. A fixed pixel ladder
#: is therefore wrong everywhere except where it was tuned -- and it failed exactly
#: that way on H2O 7160-7340, where every start (max 4.0 px) sat below the 5.26 px
#: the instrument implies there, and the LSF collapsed to zero in all six.
ATLAS_RESOLVING_POWER = 148_000.0
LSF_START_MULTIPLIERS = (1.0, 0.8, 1.25, 0.65, 1.6, 0.5)

#: A corrected pixel may carry at most this relative flux error from the
#: transmission model; T_min is derived per band from the measured residual.
CORRECTED_RELATIVE_ERROR_BUDGET = 0.05

SEGMENT_MAIN_NM = 4.0          # each lmNNNN owns 4 nm; the rest is overlap

#: Intrinsic-solar masking. molecfit fits continuum + telluric ONLY, so any solar
#: photospheric line left inside the fit region is charged to the atmosphere. On
#: the first RYA-940 pilot, with no mask, all six starting points railed the
#: Gaussian LSF at its 100-pixel bound (reduced chi2 ~1096, O2 column 0.25) --
#: the fitter smoothing its model into mush trying to cover solar structure. The
#: HARPS leg only worked because RYA-931 masked solar lines with the IAG atlas.
SOLAR_MASK_DEPTH = 0.97        # reference flux below this is intrinsic solar
SOLAR_MASK_PAD_A = 0.08        # widen each masked interval by this much


def md5(path: Path) -> str:
    d = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            d.update(chunk)
    return d.hexdigest()


def load_segments(atlas_dir: Path, lo_a: float, hi_a: float) -> tuple[np.ndarray, np.ndarray, list]:
    """Atlas flux over [lo, hi] in AIR Angstrom, overlap de-duplicated.

    Each file owns exactly its 4 nm; the trailing 0.1 nm is overlap with the next
    segment and is dropped rather than averaged, so a seam cannot be counted twice
    or silently smoothed.
    """
    waves, fluxes, used = [], [], []
    for path in sorted(atlas_dir.glob("lm[0-9]*")):
        start_nm = float(path.name[2:])
        own_lo, own_hi = start_nm * 10.0, (start_nm + SEGMENT_MAIN_NM) * 10.0
        if own_hi < lo_a or own_lo > hi_a:
            continue
        data = np.loadtxt(path)
        if data.ndim != 2 or data.shape[1] < 2:
            raise SystemExit(f"{path.name}: not an N x >=2 numeric table (corrupt segment?)")
        w = data[:, 0] * 10.0
        f = data[:, 1]
        own = (w >= own_lo) & (w < own_hi)        # drop the 0.1 nm overlap tail
        waves.append(w[own]); fluxes.append(f[own])
        used.append({"segment": path.name, "md5": md5(path),
                     "own_range_A": [own_lo, own_hi], "n_kept": int(own.sum())})
    if not waves:
        raise SystemExit(f"no atlas segment covers {lo_a}-{hi_a} A under {atlas_dir}")
    w = np.concatenate(waves); f = np.concatenate(fluxes)
    order = np.argsort(w); w, f = w[order], f[order]
    keep = (w >= lo_a) & (w <= hi_a)
    if np.any(np.diff(w[keep]) <= 0):
        raise SystemExit("atlas wavelengths are not strictly increasing after de-overlap")
    return w[keep], f[keep], used


def solar_mask_intervals(reference: Path, kind: str, lo: float, hi: float,
                         wave: np.ndarray) -> list[tuple[float, float]]:
    """Intervals of intrinsic SOLAR absorption, in Angstrom, over [lo, hi].

    The reference is telluric-free, so what it shows IS the Sun. It supplies a
    MASK only -- never flux, and never a wavelength scale for the product.
    """
    if kind == "kurucz2005":
        rw, rf = read_kurucz2005(reference, lo, hi, True)   # vacuum grid (RYA-938)
    elif kind == "iag":
        rw, rf = read_iag(reference, lo, hi)
    else:
        raise SystemExit(f"unknown solar reference kind {kind!r}")
    if rw.size < 10:
        raise SystemExit(f"solar reference {reference.name} covers none of {lo}-{hi} A. "
                         "Above 10000 A neither Kurucz 2005 nor IAG reaches; that band "
                         "has no external referee and must be handled explicitly.")
    order = np.argsort(rw); rw, rf = rw[order], rf[order]

    # Kurucz 2005 is absolute irradiance, IAG is normalised. Put whichever it is on
    # a local pseudo-continuum so one depth threshold means the same thing for both.
    window = max(3, int(round(2.0 / np.median(np.diff(rw)))))
    pad = window // 2
    padded = np.pad(rf, pad, mode="edge")
    strides = np.lib.stride_tricks.sliding_window_view(padded, window)[: rf.size]
    continuum = np.percentile(strides, 95, axis=1)
    normalised = rf / np.where(continuum > 0, continuum, np.nan)

    on_atlas = np.interp(wave, rw, normalised, left=np.nan, right=np.nan)
    step = float(np.median(np.diff(wave)))
    bad = np.isfinite(on_atlas) & (on_atlas < SOLAR_MASK_DEPTH)
    width = max(1, int(round(SOLAR_MASK_PAD_A / step)))
    from scipy.ndimage import maximum_filter1d
    bad = maximum_filter1d(bad.astype(np.uint8), size=width) > 0
    edges = np.diff(np.r_[False, bad, False].astype(np.int8))
    starts, stops = np.where(edges == 1)[0], np.where(edges == -1)[0] - 1
    return [(float(wave[a]), float(wave[b])) for a, b in zip(starts, stops)]


def even_grid(w: np.ndarray, f: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Resample onto the median step. A fitting device, not a product change."""
    step = float(np.median(np.diff(w)))
    grid = np.arange(w[0], w[-1] + 0.5 * step, step)
    return grid, np.interp(grid, w, f), step


def write_inputs(w: np.ndarray, f: np.ndarray, molecules: list[str], include: list[tuple],
                 exclude: list[tuple], destination: Path) -> Path:
    destination.mkdir(parents=True, exist_ok=True)
    finite = np.isfinite(f)
    if not finite.any():
        raise SystemExit("no finite flux in the requested band; refusing an empty SCIENCE table")
    primary = fits.PrimaryHDU()
    science = fits.BinTableHDU.from_columns([
        fits.Column(name="lambda", format="1D", unit="um", array=w[finite] / 1.0e4),
        fits.Column(name="flux", format="1D", array=f[finite]),
    ], name="SCIENCE")
    fits.HDUList([primary, science]).writeto(destination / "science.fits", overwrite=True)

    fits.BinTableHDU.from_columns([
        fits.Column(name="LOWER_LIMIT", format="1D", array=[a / 1.0e4 for a, _ in include]),
        fits.Column(name="UPPER_LIMIT", format="1D", array=[b / 1.0e4 for _, b in include]),
    ]).writeto(destination / "wave_include.fits", overwrite=True)
    fits.BinTableHDU.from_columns([
        fits.Column(name="LIST_MOLEC", format="8A", array=molecules),
        fits.Column(name="FIT_MOLEC", format="1J",
                    array=np.ones(len(molecules), dtype=np.int32)),
        fits.Column(name="REL_COL", format="1D", array=np.ones(len(molecules))),
    ]).writeto(destination / "molecules.fits", overwrite=True)

    sof_lines = ["science.fits SCIENCE", "wave_include.fits WAVE_INCLUDE",
                 "molecules.fits MOLECULES"]
    if exclude:
        fits.BinTableHDU.from_columns([
            fits.Column(name="LOWER_LIMIT", format="1D", array=[a / 1.0e4 for a, _ in exclude]),
            fits.Column(name="UPPER_LIMIT", format="1D", array=[b / 1.0e4 for _, b in exclude]),
        ]).writeto(destination / "wave_exclude.fits", overwrite=True)
        sof_lines.insert(2, "wave_exclude.fits WAVE_EXCLUDE")
    sof = destination / "model.sof"
    sof.write_text("\n".join(sof_lines) + "\n")
    return sof


def run_molecfit(sof: Path, products: Path, lsf_pix: float,
                 ref_atmos: str) -> subprocess.CompletedProcess:
    esorex = resolve_esorex()
    products.mkdir(parents=True, exist_ok=True)
    cmd = [
        esorex, f"--output-dir={products.resolve()}", "--suppress-prefix=TRUE",
        "molecfit_model",
        "--COLUMN_LAMBDA=lambda", "--COLUMN_FLUX=flux", "--COLUMN_DFLUX=NULL",
        "--DEFAULT_ERROR=0.01", "--WLG_TO_MICRON=1.0",
        # The atlas is AIR and its own frame is the observatory's: no RV term.
        "--WAVELENGTH_FRAME=AIR",
        # DO NOT FIT A CONTINUUM. Column 2 of this atlas is a PSEUDO-RESIDUAL FLUX:
        # unity IS the continuum and the product ships it. Fitting a second one is
        # the RYA-911 double-continuum defect, and it fails the same way here --
        # with a degree-3 polynomial free over a 141 A window that is one-third
        # saturated O2 A band, the "continuum" followed the band down (mscal median
        # 0.832, min 0.657 against an atlas median of 0.93), so the band was
        # explained by continuum rather than by absorption and the fitted O2 column
        # collapsed to 0.24 while the LSF railed at its 100-pixel bound.
        "--FIT_CONTINUUM=0", "--CONTINUUM_N=0", "--CONTINUUM_CONST=1.0",
        # DO NOT FIT THE WAVELENGTH SOLUTION. An FTS atlas's wavelength scale is its
        # single greatest strength, and RYA-938 verified this one independently
        # against the IAG atlas at a lag of -0.000 A. Left free it is a degree of
        # freedom with nothing to constrain it once solar lines are masked, and it
        # RAILED: `chip 1, coef 0` pinned at exactly -0.050000 (its bound, quoted to
        # +/-1e-6), imposing a spurious -3.3 A displacement on the model. That is
        # what a wide cross-correlation found -- a sharp peak at -3.316 A with
        # correlation 0.735, against 0.40 at zero lag. The fit then gave up after 8
        # iterations having improved chi2 by 0.02%.
        "--FIT_WLC=0", "--WLC_N=0", "--WLC_CONST=0.0",
        # An FTS has no slit and no boxcar; leaving those free is what let the
        # HARPS fit drive a Lorentzian to its bound and eat the column (RYA-931).
        "--FIT_RES_BOX=FALSE", "--RES_BOX=0.0",
        "--FIT_RES_LORENTZ=FALSE", "--RES_LORENTZ=0.0",
        "--FIT_RES_GAUSS=TRUE", f"--RES_GAUSS={lsf_pix:.4f}",
        "--KERNMODE=FALSE", "--KERNFAC=5.0",
        # Every observing condition is supplied as a VALUE: there are no headers
        # to read, and TELESCOPE_ANGLE is 90 deg because the airmass is unknown
        # and is absorbed by the fitted column.
        "--OBSERVING_DATE_KEYWORD=NONE", "--OBSERVING_DATE_VALUE=1983.0",
        "--UTC_KEYWORD=NONE", "--UTC_VALUE=43200.0",
        "--TELESCOPE_ANGLE_KEYWORD=NONE", "--TELESCOPE_ANGLE_VALUE=90.0",
        "--RELATIVE_HUMIDITY_KEYWORD=NONE", "--RELATIVE_HUMIDITY_VALUE=20.0",
        "--PRESSURE_KEYWORD=NONE", "--PRESSURE_VALUE=790.0",
        "--TEMPERATURE_KEYWORD=NONE", "--TEMPERATURE_VALUE=10.0",
        "--MIRROR_TEMPERATURE_KEYWORD=NONE", "--MIRROR_TEMPERATURE_VALUE=10.0",
        "--SLIT_WIDTH_KEYWORD=NONE", "--SLIT_WIDTH_VALUE=1.0",
        "--ELEVATION_KEYWORD=NONE", f"--ELEVATION_VALUE={SITE['elevation_m']}",
        "--LONGITUDE_KEYWORD=NONE", f"--LONGITUDE_VALUE={SITE['longitude_deg']}",
        "--LATITUDE_KEYWORD=NONE", f"--LATITUDE_VALUE={SITE['latitude_deg']}",
        f"--REFERENCE_ATMOSPHERIC={ref_atmos}",
        # 'none' = the bundled monthly-average profile. This IS the standard
        # atmosphere the HARPS leg was forbidden to use, taken knowingly here
        # because there is no observation night to retrieve.
        "--GDAS_PROFILE=none",
        sof.name,
    ]
    bindir = str(Path(esorex).parent)
    env = dict(os.environ, PATH=f"{bindir}:{os.environ.get('PATH', '')}")
    libdir = Path(esorex).parent.parent / "lib"
    if libdir.is_dir():
        env["LD_LIBRARY_PATH"] = f"{libdir}:{os.environ.get('LD_LIBRARY_PATH', '')}"
    return subprocess.run(cmd, cwd=sof.parent, env=env, text=True, capture_output=True)


def apply_and_verify(chosen: Path, grid: np.ndarray, resampled: np.ndarray,
                     wave: np.ndarray, flux: np.ndarray, step: float,
                     lo: float, hi: float, exclude: list,
                     reference: Path, ref_kind: str, output: Path) -> dict:
    """Divide the fitted transmission out, ON THE ORIGINAL ATLAS GRID.

    The fit ran on an even grid; the product must not inherit that resampling. The
    transmission is mapped back to the fit grid BY ROW -- BEST_FIT_MODEL's own
    `lambda` is in VACUUM and matching on it mis-registers by ~2 A (RYA-931) -- and
    is then interpolated onto the atlas's own uneven wavelengths.
    """
    model = fits.open(chosen / "out" / "BEST_FIT_MODEL.fits")[1].data
    trans_fit = np.asarray(model["mtrans"], dtype=float)
    finite = np.isfinite(resampled)
    if trans_fit.size != int(finite.sum()):
        raise SystemExit(f"model has {trans_fit.size} rows but {int(finite.sum())} "
                         "pixels were written; the row mapping is not safe")
    if not np.allclose(np.asarray(model["flux"], float), resampled[finite], rtol=0, atol=1e-9):
        raise SystemExit("BEST_FIT_MODEL flux is not row-identical to the input")
    fit_wave = grid[finite]
    inside = trans_fit > 0.0            # mtrans is 0 OUTSIDE the fitted region, not 1

    transmission = np.ones_like(wave)
    covered = (wave >= fit_wave[inside].min()) & (wave <= fit_wave[inside].max())
    transmission[covered] = np.interp(wave[covered], fit_wave[inside], trans_fit[inside])

    # dT measured on solar-CLEAN in-band pixels: a solar line left in the residual
    # would inflate it and quarantine pixels that are perfectly correctable.
    solar = np.zeros(wave.size, dtype=bool)
    for a, b in exclude:
        solar |= (wave >= a) & (wave <= b)
    band = (wave >= lo) & (wave <= hi) & covered & ~solar
    model_on_atlas = np.interp(wave, fit_wave[inside],
                               np.asarray(model["mflux"], float)[inside])
    residual = flux[band] - model_on_atlas[band]
    dt = float(1.4826 * np.median(np.abs(residual - np.median(residual))))
    t_min = dt / CORRECTED_RELATIVE_ERROR_BUDGET

    applied = covered & (transmission >= t_min)
    quarantined = covered & (transmission < t_min)
    corrected = flux.copy()
    corrected[applied] = flux[applied] / transmission[applied]
    corrected[quarantined] = np.nan

    ref = np.full(wave.size, np.nan)
    rw = np.array([])
    if ref_kind == "kurucz2005":
        rw, rf = read_kurucz2005(reference, lo - 5, hi + 5, True)
    elif ref_kind == "iag":
        rw, rf = read_iag(reference, lo - 5, hi + 5)
    if rw.size > 10:
        order = np.argsort(rw); rw, rf = rw[order], rf[order]
        win = max(3, int(round(2.0 / np.median(np.diff(rw)))))
        pad = win // 2
        cont = np.percentile(
            np.lib.stride_tricks.sliding_window_view(np.pad(rf, pad, mode="edge"), win)[: rf.size],
            95, axis=1)
        ref = np.interp(wave, rw, rf / np.where(cont > 0, cont, np.nan),
                        left=np.nan, right=np.nan)

    def stats(values):
        k = (wave >= lo) & (wave <= hi) & np.isfinite(values)
        v = values[k]
        return {"n": int(k.sum()), "min": float(v.min()), "median": float(np.median(v)),
                "pct_below_0.5": float(100 * np.mean(v < 0.5)),
                "pct_below_0.95": float(100 * np.mean(v < 0.95))}

    def vs_ref(values):
        k = (wave >= lo) & (wave <= hi) & np.isfinite(values) & np.isfinite(ref)
        if k.sum() < 10:
            return None
        return {"n": int(k.sum()),
                "rms": float(np.sqrt(np.mean((values[k] - ref[k]) ** 2))),
                "median_offset": float(np.median(values[k] - ref[k]))}

    out = output / "corrected"
    out.mkdir(parents=True, exist_ok=True)
    np.savetxt(out / f"kp1984_corrected_{int(lo)}_{int(hi)}.txt",
               np.column_stack([wave, corrected, transmission]),
               fmt="%12.5f %14.6f %10.6f",
               header=("air_wavelength_A  corrected_residual_flux  transmission" + chr(10) +
                      "RYA-940 telluric-corrected 1984 Kitt Peak atlas. NaN = quarantined "
                      f"(transmission < T_min={t_min:.4f}), never divided."))
    return {
        "dt_transmission_sigma": dt, "t_min": t_min,
        "n_corrected": int(applied.sum()), "n_quarantined": int(quarantined.sum()),
        "before": stats(flux), "after": stats(corrected),
        "vs_reference_before": vs_ref(flux), "vs_reference_after": vs_ref(corrected),
        "reference_kind": ref_kind,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("atlas_dir", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--lo", type=float, required=True, help="band low edge, air Angstrom")
    ap.add_argument("--hi", type=float, required=True, help="band high edge, air Angstrom")
    ap.add_argument("--pad", type=float, default=25.0,
                    help="clean continuum shoulder on each side, Angstrom")
    ap.add_argument("--molecules", default="O2", help="comma-separated, e.g. O2,H2O")
    ap.add_argument("--ref-atmos", default="std.fits",
                    choices=["std.fits", "equ.fits", "tro.fits"],
                    help="MIPAS reference atmosphere; Kitt Peak at 32N is mid-latitude, "
                         "so std is the default rather than the Paranal-tuned equ")
    ap.add_argument("--solar-reference", type=Path, default=None,
                    help="telluric-free solar reference used ONLY to mask intrinsic "
                         "solar lines out of the atmospheric fit")
    ap.add_argument("--solar-reference-kind", default="kurucz2005",
                    choices=["kurucz2005", "iag", "none"],
                    help="'none' skips solar masking AND external verification. Required "
                         "beyond 10000 A, where neither Kurucz 2005 nor IAG reaches, and "
                         "the resulting correction is therefore UNVALIDATED against any "
                         "external product -- which the manifest records explicitly.")
    args = ap.parse_args()

    lo, hi = args.lo - args.pad, args.hi + args.pad
    wave, flux, segments = load_segments(args.atlas_dir, lo, hi)
    grid, resampled, step = even_grid(wave, flux)

    if args.solar_reference_kind == "none":
        exclude = []
    else:
        exclude = solar_mask_intervals(args.solar_reference, args.solar_reference_kind,
                                       lo, hi, grid)
    masked_px = sum(int(round((b - a) / step)) + 1 for a, b in exclude)
    expected_fwhm_pix = (0.5 * (args.lo + args.hi) / ATLAS_RESOLVING_POWER) / step
    ladder = [round(m * expected_fwhm_pix, 3) for m in LSF_START_MULTIPLIERS]
    attempts, chosen = [], None
    for start in ladder:
        attempt = args.output / f"start_{start:.2f}".replace(".", "p")
        sof = write_inputs(grid, resampled, args.molecules.split(","),
                           [(lo, hi)], exclude, attempt / "in")
        result = run_molecfit(sof, attempt / "out", start, args.ref_atmos)
        (attempt / "esorex.stdout.txt").write_text(result.stdout)
        (attempt / "esorex.stderr.txt").write_text(result.stderr)
        best = attempt / "out" / "BEST_FIT_MODEL.fits"
        if result.returncode or not best.exists():
            attempts.append({"lsf_start_pix": start, "returncode": result.returncode,
                             "product": best.exists(), "admissible": False})
            continue
        p = {str(r["parameter"]): float(r["value"])
             for r in fits.open(attempt / "out" / "BEST_FIT_PARAMETERS.fits")[1].data}
        columns = [v for k, v in p.items() if k.startswith("rel_mol_col_")]
        ok = (p["gaussfwhm"] >= ACCEPT_MIN_LSF_FWHM_PIX
              and p["rms_rel_to_mean"] <= ACCEPT_MAX_RMS_REL_TO_MEAN
              and all(ACCEPT_MIN_REL_COL <= c <= ACCEPT_MAX_REL_COL for c in columns))
        attempts.append({"lsf_start_pix": start, "gaussfwhm": p["gaussfwhm"],
                         "reduced_chi2": p["reduced_chi2"],          # recorded, NOT gated
                         "rms_rel_to_mean": p["rms_rel_to_mean"],
                         # NOT an airmass unless the gas is WELL MIXED. rel_mol_col
                         # scales airmass x (abundance / reference abundance). For O2
                         # the abundance term is 1, so the column IS an airmass proxy
                         # and three separate O2 bands agreeing is a real cross-check.
                         # For H2O it is airmass x (PWV / reference PWV), and calling
                         # it an airmass would invent an agreement that does not exist
                         # -- the H2O bands return 1.9-2.6 against O2's 1.30-1.35, and
                         # that is water vapour, not a contradiction.
                         "column_scale": round(max(columns), 4) if columns else None,
                         "column_is_airmass_proxy": all(
                             k.split("rel_mol_col_")[-1] == "O2"
                             for k in p if k.startswith("rel_mol_col_")),
                         "returncode": 0, "admissible": ok,
                         **{k: p[k] for k in p if k.startswith("rel_mol_col_")}})
        if ok:
            chosen = attempt
            break

    correction = None
    if chosen is not None:
        correction = apply_and_verify(chosen, grid, resampled, wave, flux, step,
                                      args.lo, args.hi, exclude,
                                      args.solar_reference, args.solar_reference_kind,
                                      args.output)

    manifest = {
        "ticket": "RYA-940", "band_A": [args.lo, args.hi], "fit_window_A": [lo, hi],
        "molecules": args.molecules.split(","), "reference_atmosphere": args.ref_atmos,
        "gdas": "none (bundled monthly average) - NO observation-night profile exists "
                "for this atlas; recorded deviation ratified 2026-08-21",
        "site": SITE, "segments": segments,
        "n_atlas_pixels": int(wave.size), "resample_step_A": step,
        "n_fit_pixels": int(grid.size),
        "resampling_note": "fit runs on the even grid; transmission is applied back on "
                           "the ORIGINAL atlas grid, which the product keeps",
        "solar_reference": str(args.solar_reference) if args.solar_reference else None,
        "solar_reference_kind": args.solar_reference_kind,
        "solar_reference_md5": md5(args.solar_reference) if args.solar_reference else None,
        "externally_validated": args.solar_reference_kind != "none",
        "solar_mask_intervals": len(exclude),
        "solar_mask_pixels": masked_px,
        "solar_mask_fraction": round(masked_px / max(1, grid.size), 4),
        "atlas_resolving_power": ATLAS_RESOLVING_POWER,
        "expected_lsf_fwhm_pix": expected_fwhm_pix,
        "lsf_start_ladder_pix": ladder,
        "lsf_attempts": attempts,
        "accepted": str(chosen) if chosen else None,
        "correction": correction,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "fit_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"band": [args.lo, args.hi], "segments": len(segments),
                      "atlas_px": int(wave.size), "fit_px": int(grid.size),
                      "solar_mask_intervals": len(exclude),
                      "solar_mask_fraction": round(masked_px / max(1, grid.size), 4),
                      "attempts": attempts, "accepted": bool(chosen)}, indent=2))
    if chosen is None:
        raise SystemExit("no physically admissible fit from any starting point")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run ESO Molecfit `molecfit_model` on one HARPS direct-solar exposure (RYA-931).

RYA-927 left the runtime blocked at rc=13 with molecfit reporting

    mf_config_chk_ref_atmos: Unable to open rat data file .../profiles/mipas/equ.fits!
    molecfit_model: Bad Specifications of Molecules and Reference Atmosphere

which reads as a broken MIPAS reference atmosphere.  It is not.  The HARPS
Phase-3 S1D `ERR` column is entirely NaN in all ten exposures, so a
finite-error filter selects zero pixels; molecfit then fails in
`molecfit_data_extract_spectrum_from_cpl_table` ("Access beyond boundaries",
"Parse 1D spectrum from BINTABLE failed") and only *afterwards* prints the
reference-atmosphere line off an already-dirty CPL error state.  The MIPAS
message is a downstream artefact of an empty SCIENCE table.

This runner therefore refuses to emit an empty SCIENCE table, and reports the
real reason when a source column is unusable.  It never substitutes a standard
atmosphere: GDAS_PROFILE is required and checksummed.

Wavelength frame: HARPS S1D is air, SPECSYS=BARYCENT, so the telluric absorber
is *not* at rest in the file's frame.  `--frame` selects how that is declared to
molecfit; `AIR_RV` carries the observatory radial velocity in
`OBS_ERF_RV_VALUE`, whose sign convention (positive when the target-Earth
distance increases) is the negative of the ESO DRS BERV.  The sign is not
asserted here -- scripts/rya931_frame_sign_test.py measures it.
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
from scipy.ndimage import gaussian_filter1d, maximum_filter1d

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.wavelength_util import vac_to_air  # noqa: E402

#: RYA-939 — resolve esorex, never hard-code it.  molecfit lives in different
#: places on the two machines (Homebrew on the Mac, an ESO source-kit prefix on
#: Sirius), and a Mac path baked into the runner made this whole leg
#: Mac-only for no reason.  Order: $ESOREX, then PATH, then the two known
#: install prefixes.  Loud on absence, naming what was tried -- a missing engine
#: must not read as a failed correction.
def _esorex_candidates() -> tuple[str, ...]:
    from config.constants import codex_root
    return (str(codex_root("eso_pipelines") / "bin" / "esorex"),   # Sirius kit
            str(Path("/opt") / "homebrew" / "bin" / "esorex"))     # Mac brew tap


def resolve_esorex() -> str:
    from shutil import which
    explicit = os.environ.get("ESOREX")
    if explicit:
        if not Path(explicit).is_file():
            raise SystemExit(f"$ESOREX={explicit!r} is not a file")
        return explicit
    found = which("esorex")
    if found:
        return found
    for candidate in _esorex_candidates():
        if Path(candidate).is_file():
            return candidate
    raise SystemExit(
        "esorex not found. Looked at $ESOREX, PATH, and "
        + ", ".join(_esorex_candidates())
        + ". Install the ESO molecfit kit or set $ESOREX. Failing here rather "
          "than reporting an uncorrected product -- a missing engine is not a "
          "correction result.")
# 🔴 RYA-1191 — THE FIT WINDOW IS DERIVED FROM THE BAND REGISTRY, NOT TYPED HERE.
#
# It used to be one hardcoded pair, `INCLUDE_LO_A, INCLUDE_HI_A = 6857.0, 6911.0`, chosen
# when the intake manifest could truthfully record
# `registered_telluric_bands_in_instrument_range: ['O2 B 6867-6884']` -- one band, so one
# window. RYA-1193 then added **o2gamma 6270-6300** to `telluric_policy.TELLURIC_BANDS`,
# and nothing re-scoped this fit. Measured consequence (RYA-1191): HARPS carries the O2
# gamma lines at **+12.0 sigma** against a displaced null, IDENTICALLY in its raw and its
# "corrected" product -- molecfit changed the flux there by 1e-5, because 6270-6300 was
# never inside the window it was given.
#
# That is the same defect this ticket exists for, one level up: a BAND-LIMITED correction
# consumed as if it were global. A hardcoded window cannot notice that the registry grew,
# so the window is now computed from the registry and the fit follows it automatically.
# The molecule list needs no change -- the O2 gamma band IS O2.
#
# ⚠️ The SHOULDER is the reason a band is padded rather than used raw: molecfit fits a
# continuum polynomial across the include region, and it must be anchored on telluric-free
# pixels either side or the continuum absorbs the band.
_SHOULDER_A = 10.0
_FITTED_MOLECULES = frozenset({"O2"})

#: 🔴 REGIONS THE REGISTRY UNDERSTATES, KEPT AS A FLOOR SO THIS CAN ONLY EVER WIDEN.
#: `TELLURIC_BANDS` is a LINE-SELECTION registry: its edges say where a line should be
#: refused, not where a fit should be constrained. For O2 B it lists 6867-6884, but
#: RYA-931 fitted 6857-6911 because the band measurably runs from 6866 A to beyond the
#: S1D red cutoff at 6912.8 A. Deriving the window from the registry alone would have
#: SHRUNK a working correction by 17 A -- a silent regression dressed as a cleanup. The
#: derived regions are unioned with these, never substituted for them.
_EMPIRICAL_REGIONS = (
    (6857.0, 6911.0, "O2 B — RYA-931 measured extent; the registered 6867-6884 "
                     "understates it and the band runs past the S1D cutoff"),
)


def include_regions_A(lo_cut: float = 0.0, hi_cut: float = 1.0e9) -> list[tuple[float, float]]:
    """Fit regions (Å) = every registered telluric band of a FITTED molecule, widened by a
    continuum shoulder, UNIONED with the empirically established regions above.

    Single-sourced on `telluric_policy.TELLURIC_BANDS` so the fit cannot silently fall
    behind the registry again, and floored on `_EMPIRICAL_REGIONS` so it cannot silently
    fall behind the data either. Overlapping regions are merged — molecfit takes
    WAVE_INCLUDE as disjoint intervals.

    ⚠️ THE MOLECULE TEST IS ON THE NAME'S FIRST TOKEN, NOT A SUBSTRING. `"O2" in "CO2"`
    is True, and the substring version duly pulled the CO2 15700-16100 band into a HARPS
    O2 fit — a band 9000 A outside the instrument.
    """
    from pipeline.telluric_policy import TELLURIC_BANDS
    want = []
    for lo, hi, name in TELLURIC_BANDS:
        tok = str(name).split()[0].upper() if str(name).strip() else ""
        if tok not in _FITTED_MOLECULES:
            continue
        want.append((float(lo) - _SHOULDER_A, float(hi) + _SHOULDER_A))
    want.extend((a, b) for a, b, _ in _EMPIRICAL_REGIONS)
    clipped = []
    for a, b in want:
        a, b = max(a, lo_cut), min(b, hi_cut)
        if b - a > 1.0:
            clipped.append((a, b))
    clipped.sort()
    merged: list[list[float]] = []
    for a, b in clipped:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return [(round(a, 3), round(b, 3)) for a, b in merged]


#: The widest span the science table must carry: every fit region plus a margin.
EXTRACT_MARGIN_A = 17.0
SOLAR_ATLAS_DEPTH = 0.985      # atlas flux below this is intrinsic solar absorption
SOLAR_ATLAS_PAD_A = 0.12       # widen each solar interval by this much


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def solar_exclusions(atlas: Path, wave_a: np.ndarray) -> list[tuple[float, float]]:
    """Intrinsic-solar absorption intervals, in microns, over the include region.

    Molecfit fits continuum + telluric only; solar photospheric lines inside the
    fit region are pulled into the telluric column density.  The Baker+2020 IAG
    telluric-free atlas marks where those lines are, without importing its flux.
    """
    with fits.open(atlas) as hdul:
        data = hdul[1].data
        atlas_wave = vac_to_air(1.0e8 / np.asarray(data["v"], dtype=float))
        atlas_flux = np.asarray(data["s"], dtype=float)
    order = np.argsort(atlas_wave)
    atlas_wave, atlas_flux = atlas_wave[order], atlas_flux[order]
    atlas_flux = gaussian_filter1d(atlas_flux, 0.03 / float(np.median(np.diff(atlas_wave))))
    solar = np.interp(wave_a, atlas_wave, atlas_flux, left=np.nan, right=np.nan)
    pixel = float(np.median(np.diff(wave_a)))
    bad = np.isfinite(solar) & (solar < SOLAR_ATLAS_DEPTH)
    bad = maximum_filter1d(bad.astype(np.uint8),
                           size=max(1, int(round(SOLAR_ATLAS_PAD_A / pixel)))) > 0
    _inc = include_regions_A()
    _in_any = np.zeros_like(wave_a, dtype=bool)
    for _a, _b in _inc:
        _in_any |= (wave_a >= _a) & (wave_a <= _b)
    bad &= _in_any
    edges = np.diff(np.r_[False, bad, False].astype(np.int8))
    starts, stops = np.where(edges == 1)[0], np.where(edges == -1)[0] - 1
    return [(float(wave_a[a] / 1.0e4), float(wave_a[b] / 1.0e4))
            for a, b in zip(starts, stops)]


def write_inputs(source: Path, atlas: Path | None, destination: Path) -> dict:
    destination.mkdir(parents=True, exist_ok=True)
    with fits.open(source) as hdul:
        header = hdul[0].header
        if header.get("INSTRUME") != "HARPS":
            raise SystemExit(f"{source.name}: INSTRUME={header.get('INSTRUME')!r}, expected HARPS")
        if header.get("SPECSYS") != "BARYCENT":
            raise SystemExit(f"{source.name}: SPECSYS={header.get('SPECSYS')!r}, expected BARYCENT")
        row = hdul[1].data[0]
        wave = np.asarray(row["WAVE"], dtype=float)
        flux = np.asarray(row["FLUX"], dtype=float)
        error = np.asarray(row["ERR"], dtype=float)
        berv = float(header["HIERARCH ESO DRS BERV"])
        _regs = include_regions_A()
        _elo = min(a for a, _ in _regs) - EXTRACT_MARGIN_A
        _ehi = max(b for _, b in _regs) + EXTRACT_MARGIN_A
        keep = (wave >= _elo) & (wave <= _ehi) & np.isfinite(flux)

        # RYA-931 root cause: ERR is 100% NaN in every HARPS solar exposure.
        # Selecting on finite ERR yields a zero-row table, which molecfit
        # misreports as a reference-atmosphere failure.  Weight uniformly and
        # say so, rather than silently dropping every pixel.
        finite_error = int(np.count_nonzero(np.isfinite(error[keep]) & (error[keep] > 0)))
        error_usable = finite_error == int(keep.sum()) and finite_error > 0
        if not keep.any():
            raise SystemExit(
                f"{source.name}: no finite FLUX pixel in {_elo:.1f}-{_ehi:.1f} A; "
                "refusing to write an empty SCIENCE table")

        primary = fits.PrimaryHDU()
        wanted = {"MJD-OBS", "UTC", "INSTRUME", "SPECSYS", "ESO TEL ALT",
                  "ESO TEL AMBI RHUM", "ESO TEL AMBI PRES START", "ESO TEL AMBI TEMP",
                  "ESO TEL GEOELEV", "ESO TEL GEOLON", "ESO TEL GEOLAT"}
        for key in header:
            if key in wanted:
                primary.header[key] = header[key]
        columns = [
            fits.Column(name="lambda", format="1D", unit="um", array=wave[keep] / 1.0e4),
            fits.Column(name="flux", format="1D", unit="count", array=flux[keep]),
        ]
        if error_usable:
            columns.append(fits.Column(name="dflux", format="1D", unit="count",
                                       array=error[keep]))
        science = fits.BinTableHDU.from_columns(columns, name="SCIENCE")
        fits.HDUList([primary, science]).writeto(destination / "science.fits", overwrite=True)

    n_rows = int(keep.sum())
    exclusions = solar_exclusions(atlas, wave[keep]) if atlas else []

    # 🔴 MULTI-REGION. molecfit's WAVE_INCLUDE is a TABLE of intervals and always was;
    # writing one row was a choice made when the registry listed one band, not a limit of
    # the recipe. Regions are clipped to what this exposure actually covers, so a band
    # outside the S1D range is dropped here rather than handed to molecfit as empty.
    _cov_lo, _cov_hi = float(np.nanmin(wave[keep])), float(np.nanmax(wave[keep]))
    _regions = [(a, b) for a, b in include_regions_A(_cov_lo, _cov_hi) if b > a + 1.0]
    if not _regions:
        raise SystemExit(
            f"no registered telluric band of {_FITTED_MOLECULES} lies inside this "
            f"exposure's {_cov_lo:.1f}-{_cov_hi:.1f} A coverage — refusing to run a fit "
            f"with an empty WAVE_INCLUDE rather than emit a no-op 'correction'.")
    fits.BinTableHDU.from_columns([
        fits.Column(name="LOWER_LIMIT", format="1D",
                    array=[a / 1.0e4 for a, _ in _regions]),
        fits.Column(name="UPPER_LIMIT", format="1D",
                    array=[b / 1.0e4 for _, b in _regions]),
    ]).writeto(destination / "wave_include.fits", overwrite=True)
    fits.BinTableHDU.from_columns([
        fits.Column(name="LIST_MOLEC", format="4A", array=["O2"]),
        fits.Column(name="FIT_MOLEC", format="1J", array=np.array([1], dtype=np.int32)),
        fits.Column(name="REL_COL", format="1D", array=[1.0]),
    ]).writeto(destination / "molecules.fits", overwrite=True)
    # esorex splits SOF lines on whitespace, so an absolute path under
    # "~/Documents/Exoplanet Codex/" is torn in two ("Could not open the input
    # file '.../Exoplanet' with tag 'Codex/...'").  Name the frames relatively
    # and run the recipe with cwd=destination.
    sof_lines = ["science.fits SCIENCE",
                 "wave_include.fits WAVE_INCLUDE",
                 "molecules.fits MOLECULES"]
    if exclusions:
        fits.BinTableHDU.from_columns([
            fits.Column(name="LOWER_LIMIT", format="1D", array=[x[0] for x in exclusions]),
            fits.Column(name="UPPER_LIMIT", format="1D", array=[x[1] for x in exclusions]),
        ]).writeto(destination / "wave_exclude.fits", overwrite=True)
        sof_lines.insert(2, "wave_exclude.fits WAVE_EXCLUDE")
    sof = destination / "model.sof"
    sof.write_text("\n".join(sof_lines) + "\n")
    return {"sof": sof, "berv_kms": berv, "science_rows": n_rows,
            "pixel_step_a": float(np.median(np.diff(wave[keep]))),
            "source_error_column_usable": error_usable,
            "source_error_finite_pixels": finite_error,
            "solar_exclusion_count": len(exclusions),
            "solar_exclusion_intervals_um": exclusions}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("gdas", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--solar-atlas", type=Path, default=None)
    parser.add_argument("--frame", default="AIR_RV", choices=["AIR", "AIR_RV"])
    parser.add_argument("--rv-sign", type=float, default=-1.0,
                        help="OBS_ERF_RV_VALUE = rv_sign * BERV")
    parser.add_argument("--resolving-power", type=float, default=115000.0,
                        help="HARPS R; sets the initial Gaussian LSF FWHM in pixels")
    parser.add_argument("--lsf-init-fwhm-pix", type=float, default=None,
                        help="start the Gaussian LSF fit from this measured FWHM in "
                             "pixels instead of from the nominal resolving power. The "
                             "LSF is always FITTED: FIT_RES_GAUSS=FALSE does not hold "
                             "the kernel at the given width, it disables convolution "
                             "(same 6.767 px gives reduced chi2 0.81 fitted vs 894 "
                             "fixed, and collapses the O2 column to 0.004).")
    args = parser.parse_args()

    if not args.gdas.is_file():
        raise SystemExit(f"GDAS profile {args.gdas} missing; RYA-931 forbids a "
                         "standard-atmosphere fallback")
    inputs, products = args.output / "in", args.output / "out"
    products.mkdir(parents=True, exist_ok=True)
    built = write_inputs(args.source.resolve(),
                         args.solar_atlas.resolve() if args.solar_atlas else None,
                         inputs.resolve())

    esorex = resolve_esorex()
    obs_rv = args.rv_sign * built["berv_kms"] if args.frame == "AIR_RV" else 0.0
    # LSF init only: the widest fit region's centre. The kernel is FITTED (RYA-931), so
    # this sets a starting width, not the answer.
    _regs = include_regions_A()
    _widest = max(_regs, key=lambda r: r[1] - r[0])
    band_centre_a = 0.5 * (_widest[0] + _widest[1])
    gauss_fwhm_pix = (args.lsf_init_fwhm_pix if args.lsf_init_fwhm_pix else
                      (band_centre_a / args.resolving_power) / built["pixel_step_a"])
    cmd = [
        esorex, f"--output-dir={products.resolve()}",
        # RYA-939 — name products by their PRO CATG, explicitly.  esorex takes this
        # from ~/.esorex/esorex.rc when not given, and the two machines disagree: the
        # Mac has no rc at all (compiled default), while a fresh ESO source-kit install
        # writes suppress-prefix=FALSE.  Under FALSE the recipe still succeeds and
        # still writes all ten products -- as out_0000.fits ... out_0009.fits -- so
        # BEST_FIT_MODEL.fits is simply absent and the run reads as a failure with
        # rc=0.  Passing it here makes the product names a property of this recipe
        # invocation rather than of whoever's home directory it ran in.
        "--suppress-prefix=TRUE",
        "molecfit_model",
        "--COLUMN_LAMBDA=lambda", "--COLUMN_FLUX=flux",
        f"--COLUMN_DFLUX={'dflux' if built['source_error_column_usable'] else 'NULL'}",
        "--DEFAULT_ERROR=0.01", "--WLG_TO_MICRON=1.0",
        f"--WAVELENGTH_FRAME={args.frame}", f"--OBS_ERF_RV_VALUE={obs_rv}",
        "--FIT_CONTINUUM=1", "--CONTINUUM_N=2", "--FIT_WLC=1", "--WLC_N=1",
        # HARPS is fibre-fed with a near-Gaussian LSF and no slit, so the
        # boxcar and Lorentzian components have no physical counterpart here.
        # Left at their defaults they are free parameters with nothing to
        # constrain them: the first RYA-931 run drove the Lorentzian FWHM to
        # its 100-pixel bound and paid for the over-smoothed model by
        # collapsing the O2 column to 3.6% of atmospheric.
        "--FIT_RES_BOX=FALSE", "--RES_BOX=0.0",
        "--FIT_RES_LORENTZ=FALSE", "--RES_LORENTZ=0.0",
        "--FIT_RES_GAUSS=TRUE", f"--RES_GAUSS={gauss_fwhm_pix:.4f}",
        "--KERNMODE=FALSE", "--KERNFAC=5.0",
        "--MIRROR_TEMPERATURE_KEYWORD=NONE", "--MIRROR_TEMPERATURE_VALUE=20.9",
        "--SLIT_WIDTH_KEYWORD=NONE", "--SLIT_WIDTH_VALUE=1.0",
        "--ELEVATION_VALUE=2400.0", "--LONGITUDE_VALUE=-70.7345",
        "--LATITUDE_VALUE=-29.2584", f"--GDAS_PROFILE={args.gdas.resolve()}",
        built["sof"].name,
    ]
    bindir = str(Path(esorex).parent)
    env = dict(os.environ, PATH=f"{bindir}:{os.environ.get('PATH', '')}")
    libdir = Path(esorex).parent.parent / "lib"
    if libdir.is_dir():                     # ESO source-kit installs need this
        env["LD_LIBRARY_PATH"] = f"{libdir}:{os.environ.get('LD_LIBRARY_PATH', '')}"
    result = subprocess.run(cmd, cwd=inputs, env=env, text=True, capture_output=True)
    (args.output / "esorex.stdout.txt").write_text(result.stdout)
    (args.output / "esorex.stderr.txt").write_text(result.stderr)

    manifest = {
        "ticket": "RYA-931",
        "source": str(args.source.resolve()), "source_md5": md5(args.source),
        "gdas_profile": str(args.gdas.resolve()), "gdas_md5": md5(args.gdas),
        "standard_atmosphere_fallback": False,
        "solar_exclusion_atlas": str(args.solar_atlas.resolve()) if args.solar_atlas else None,
        "solar_exclusion_atlas_md5": md5(args.solar_atlas) if args.solar_atlas else None,
        "esorex": esorex,
        "wavelength_frame": args.frame, "berv_kms": built["berv_kms"],
        "resolving_power": args.resolving_power,
        "pixel_step_a": built["pixel_step_a"],
        "lsf_gauss_fwhm_pix_initial": gauss_fwhm_pix,
        "lsf_box_fitted": False, "lsf_lorentz_fitted": False,
        "lsf_gauss_fitted": True,
        "obs_erf_rv_value_kms": obs_rv,
        "science_rows": built["science_rows"],
        "source_error_column_usable": built["source_error_column_usable"],
        "source_error_finite_pixels": built["source_error_finite_pixels"],
        "default_error_relative": None if built["source_error_column_usable"] else 0.01,
        "default_error_role": "fit weighting only; never an EW uncertainty",
        "solar_exclusion_count": built["solar_exclusion_count"],
        "solar_exclusion_intervals_um": built["solar_exclusion_intervals_um"],
        "command": cmd, "returncode": result.returncode,
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    best = products / "BEST_FIT_MODEL.fits"
    if result.returncode:
        raise SystemExit(f"molecfit_model FAILED (rc={result.returncode}); see {args.output}")
    if not best.exists():
        # Distinct from a failed recipe, and it used to be reported as one
        # ("failed (rc=0)"), which sent RYA-939 hunting a fit problem that did not
        # exist.  Name what is actually wrong: the recipe ran, the products are there,
        # they are not called what we expected.
        produced = sorted(f.name for f in products.glob("*.fits"))
        raise SystemExit(
            f"molecfit_model SUCCEEDED (rc=0) but {best.name} is missing from "
            f"{products}. Produced: {produced or 'nothing'}. Generic out_NNNN.fits "
            f"names mean --suppress-prefix did not take effect.")
    print(best)


if __name__ == "__main__":
    main()

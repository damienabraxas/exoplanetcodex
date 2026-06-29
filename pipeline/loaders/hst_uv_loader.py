"""
pipeline/loaders/hst_uv_loader.py
=================================
Loader for HST STIS/COS ultraviolet 1D spectra (Procyon FUV/NUV; RYA-471).
The UV analogue of the RYA-272 UVES loader: it turns a staged, audited HASP
``_cspec.fits`` into a pipeline ``SpectrumData`` in the pipeline wavelength frame,
conditioned through the standing RYA-426 UV stage.

What this loader is (and is NOT)
--------------------------------
  * It LOADS + CONDITIONS the UV arm: read the native-vacuum HASP product, take it
    to the pipeline frame (NUV vac->air, FUV stays vacuum — RYA-426 keystone),
    mask the chromospheric emission cores, and stamp the RYA-426 conditioning
    manifest (gates + analysis_ready) as provenance. provenance = MEASURED.
  * It does NOT measure C/N/O. The FUV C I / O I / S I lines are synthesis-only
    (no continuum in the FUV — RYA-426 gate 5) and carry a LARGE negative NLTE
    correction supplied by the Amarsi C/O grid (RYA-359), which is not yet on
    disk: ``cno_synthesis.amarsi_grid_backend`` loud-fails by design. So the UV
    arm stays DEFERRED in the RYA-464 registry on those two downstream blockers
    (NLTE grid + FUV pseudo-continuum) even though the loader is built and
    smoke-proven here. Flipping the arm to ready is a one-line change once the
    grid lands — the loader, conditioning, registry dispatch and diagnostics are
    all wired below.

Data discipline (RYA-222 / RYA-262)
-----------------------------------
  * Files come ONLY from the RYA-222 science-ready whitelist
    (``data/audit/procyon_hst/procyon_hst_science_ready.csv``). The "Procyon HST"
    MAST tree ALSO holds 55 Cnc (srho01cnc / rho-cnc) frames — NEVER glob the
    directory; the whitelist + the per-file target guard are the safeguard.
  * Target guard: TARGNAME must confirm Procyon (HD 61421 / alpha CMi). A
    55 Cnc / unknown target raises loudly — a UV frame is never routed by header.

References
----------
  HST STIS: STIS Instrument Handbook; HASP coadd products.
  Vac<->air: Birch & Downs 1994 (via pipeline.uv_conditioning, RYA-303/426).
  Linear: RYA-222 (data audit) / RYA-262 (UV coverage gate) / RYA-471 (this loader).
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from astropy.io import fits

import config.constants as const
from . import base_loader
from .base_loader import SpectrumData
from pipeline import uv_conditioning as uvc
from pipeline import uv_line_selection as uvls

ROOT = Path(str(const.ROOT))

# RYA-222 science-ready whitelist (the ONLY source of file paths — never glob).
SCIENCE_READY_CSV = ROOT / 'data' / 'audit' / 'procyon_hst' / 'procyon_hst_science_ready.csv'

# Accepted target identities for Procyon (the per-file guard against 55 Cnc / strays).
_PROCYON_TARGETS = ('HD61421', 'HD 61421', 'HD061421', 'PROCYON', 'ALF CMI', 'ALPHA CMI',
                    'ALF-CMI', 'ALPHACMI', 'HR2943')
_FORBIDDEN_TARGETS = ('CNC', '55CNC', 'RHO', 'SRHO')   # 55 Cnc lives in the same tree


class HstUvProductError(RuntimeError):
    """Raised when a file is not a loadable Procyon STIS/COS UV product (grating /
    target / structure guard failure). Names the file and the specific failed guard."""


# ── low-level HASP cspec reader ───────────────────────────────────────────────
def _read_cspec(path: Path):
    """Read one HASP/STIS ``_cspec.fits`` -> (wave_vac_A, flux, err|None, targname).
    Wavelengths are native VACUUM Angstrom (HST UV convention); rows are sorted by wl."""
    with fits.open(path) as h:
        names = [e.name for e in h]
        sci = h['SCI'].data if 'SCI' in names else h[1].data
        w = np.asarray(sci['WAVELENGTH']).ravel().astype(float)
        f = np.asarray(sci['FLUX']).ravel().astype(float)
        cols = sci.columns.names
        e = np.asarray(sci['ERROR']).ravel().astype(float) if 'ERROR' in cols else None
        targ = str(h[0].header.get('TARGNAME', '') or h[0].header.get('OBJECT', ''))
    order = np.argsort(w)
    return w[order], f[order], (e[order] if e is not None else None), targ


def _norm_target(s: str) -> str:
    return ''.join(ch for ch in str(s).upper() if ch.isalnum())


def _assert_procyon(targname: str, filename: str) -> None:
    """Loud target guard: confirm Procyon, reject 55 Cnc / strays (never route by glob)."""
    t = _norm_target(targname)
    if any(_norm_target(b) in t for b in _FORBIDDEN_TARGETS):
        raise HstUvProductError(
            f"{filename}: TARGNAME={targname!r} is a non-Procyon target (55 Cnc / rho Cnc "
            f"share the 'Procyon HST' tree) — refusing to load it as Procyon (RYA-222).")
    if not any(_norm_target(a) in t or t in _norm_target(a) for a in _PROCYON_TARGETS):
        raise HstUvProductError(
            f"{filename}: TARGNAME={targname!r} does not confirm Procyon (HD 61421 / alpha CMi). "
            f"A UV frame is never routed by header — verify it is on the RYA-222 whitelist.")


# ── the loader ────────────────────────────────────────────────────────────────
class HstUvLoader:
    """Load a single HST STIS/COS UV ``_cspec.fits`` into a conditioned SpectrumData.

        spec = HstUvLoader(path, grating='E140M').load()
        # spec.wave_A — pipeline frame (NUV air, FUV vacuum), Angstrom
        # spec.flux   — flux-calibrated; chromospheric cores -> NaN (RYA-426 gate 4)
        # spec.meta   — provenance incl. conditioning manifest + analysis_ready
    """

    def __init__(self, filepath, grating: str, *, star: str = 'procyon',
                 distance_pc: float = 3.51):
        self.filepath = Path(filepath)
        self.grating = str(grating).upper()
        self.star = star
        self.distance_pc = float(distance_pc)

    def load(self) -> SpectrumData:
        fn = self.filepath.name
        if self.grating in uvc.EXCLUDED_GRATINGS:
            raise HstUvProductError(f"{fn}: grating {self.grating} is permanently excluded (RYA-426 §2).")
        if self.grating not in uvc.STIS_UV_GRATINGS:
            raise HstUvProductError(
                f"{fn}: unknown UV grating {self.grating!r} — add it to STIS_UV_GRATINGS "
                f"(uv_conditioning, single source of truth) before loading.")
        if not self.filepath.exists():
            raise HstUvProductError(f"{self.filepath}: file not present on disk.")

        wave_vac, flux, err, targ = _read_cspec(self.filepath)
        _assert_procyon(targ, fn)

        # Pipeline frame (NUV vac->air, FUV stays vacuum — the RYA-426 keystone).
        wave_pf = uvc.to_pipeline_frame(wave_vac)

        # RYA-426 conditioning manifest (gates + analysis_ready) on the RAW vacuum grid.
        cond = uvc.condition_uv_spectrum(
            wave_vac, flux, err, star=self.star, grating=self.grating, instrument='STIS',
            distance_pc=self.distance_pc, is_binary=False, identity_confirmed=True, write=True)

        # Mask chromospheric emission cores -> NaN (a filled core is a silent overestimate).
        core = uvc.chromospheric_core_mask(wave_pf)
        flux_out = flux.copy()
        flux_out[core] = np.nan
        err_out = None
        if err is not None:
            err_out = err.copy()
            err_out[core] = np.nan

        gcfg = uvc.STIS_UV_GRATINGS[self.grating]
        meta = {
            'instrument'         : 'STIS',
            'object'             : targ or 'Procyon',
            'date_obs'           : 'UNKNOWN',          # HASP coadd: per-exposure dates in the audit
            'exptime_s'          : 0.0,
            'snr_summary'        : float(cond.gates.get('scattered_light', {}).get('neg_flux_fraction', -1.0)),
            'apero_version'      : 'hasp',             # base contract's reduction-software slot
            'berv_kms'           : 0.0,                # HASP products are already rest-frame coadds
            'wave_units_raw'     : 'angstrom (vacuum)',
            'telluric_corrected' : False,              # space UV — no telluric (always recorded)
            'wapiti_applied'     : False,
            'filepath'           : str(self.filepath),
            # UV-specific provenance
            'provenance'         : 'measured',
            'grating'            : self.grating,
            'regime'             : gcfg['regime'],
            'detector'           : gcfg['detector'],
            'resolution_R'       : float(gcfg['R']),
            'frame_out'          : cond.frame_out,     # 'air>=2000/vac<2000'
            'native_frame'       : 'vacuum',
            'masked_fraction'    : cond.masked_fraction,
            'analysis_ready'     : cond.analysis_ready,
            'ew_allowed'         : cond.gates.get('synthesis_gate', {}).get('ew_allowed', False),
            'conditioning_schema': cond.schema,
            'conditioning_gates' : {k: v.get('passed', None) if isinstance(v, dict) else None
                                    for k, v in cond.gates.items()},
            'distance_pc'        : self.distance_pc,
        }
        return SpectrumData(wave_A=wave_pf, flux=flux_out, err=err_out, meta=meta)


# ── whitelist access (NEVER glob the tree) ────────────────────────────────────
def science_ready_rows(grating: str = None):
    """RYA-222 whitelisted, science-grating, target-confirmed, on-disk rows. Optional
    grating filter. The ONLY way file paths enter the loader — the tree holds 55 Cnc frames."""
    if not SCIENCE_READY_CSV.exists():
        raise HstUvProductError(
            f"RYA-222 science-ready whitelist absent ({SCIENCE_READY_CSV}). Stage + audit the "
            f"Procyon HST frames (RYA-222) before loading — never glob the MAST tree.")
    rows = []
    with open(SCIENCE_READY_CSV) as fh:
        for r in csv.DictReader(fh):
            if str(r.get('is_science_grating', '')).strip() not in ('True', 'true', '1'):
                continue
            if str(r.get('target_confirmed', '')).strip() not in ('True', 'true', '1'):
                continue
            if not Path(r['filepath']).exists():
                continue
            if grating is not None and str(r.get('opt_elem', '')).upper() != grating.upper():
                continue
            rows.append(r)
    return rows


def _pick_arm_frame(rows):
    """Pick the arm spectrum from whitelist rows: PREFER the merged HASP coadd
    (``_cspec.fits``, a monotonic 1D product) over a raveled echelle ``x1d`` (overlapping
    orders), then the widest coverage. Falls back to widest if no cspec is present."""
    cspec = [r for r in rows if str(r.get('filename', r['filepath'])).lower().endswith('_cspec.fits')]
    pool = cspec or rows
    return max(pool, key=lambda r: float(r['wl_max_A']) - float(r['wl_min_A']))


def load_procyon_uv_arm(grating: str = 'E140M', *, return_spectrum: bool = False):
    """Load Procyon's UV arm for a grating from the RYA-222 whitelist -> (wave_nm, flux)
    for the RYA-464 resolver, or the full SpectrumData with ``return_spectrum=True``.

    The widest-coverage whitelisted frame for the grating is used (the merged coadd:
    E140M 1140-1730 covers C I 1657 + O I 1355). Flux is normalised to a working
    PSEUDO-continuum (robust high percentile) so the resolver has a unit-ish baseline;
    in the FUV there is NO true continuum (RYA-426 gate 5), so this is explicitly a
    pseudo-continuum and a science-grade C I synthesis is OWED a proper FUV
    pseudo-continuum + the Amarsi NLTE grid (RYA-359). Flagged, never silently 'normalised'.
    """
    rows = science_ready_rows(grating=grating)
    if not rows:
        raise HstUvProductError(
            f"no whitelisted, on-disk, target-confirmed Procyon frame for grating {grating!r}.")
    row = _pick_arm_frame(rows)
    spec = HstUvLoader(row['filepath'], grating=grating).load()
    if return_spectrum:
        return spec
    wave_nm = spec.wave_A / 10.0
    flux = np.asarray(spec.flux, float)
    finite = flux[np.isfinite(flux)]
    pseudo_cont = np.nanpercentile(finite, 95) if finite.size else 1.0
    norm = flux / pseudo_cont if pseudo_cont else flux
    return wave_nm, norm


# ── UV diagnostic set for the RYA-464 registry (from RYA-190, synthesis FUV lines) ──
def uv_arm_diagnostics():
    """Build the UV arm's Diagnostic tuple from RYA-190's USABLE FUV synthesis lines.
    Single source of truth = uv_line_selection (RYA-190); this only maps verdict=USE
    synthesis lines into cno_synthesis Diagnostic objects. Imported lazily to avoid a
    loader<->synthesis import cycle. C I 1657 is the primary; O I 1355 / S I 1474 ride along.
    NLTE is GRID_OWED (Amarsi RYA-359) — flagged loud, never applied as a scalar here."""
    from pipeline.cno_synthesis import Diagnostic
    role_primary = {'C I': 1657.38}        # the Procyon C advantage over the solar cited composite
    out = []
    for d in uvls.usable_diagnostics():
        if d['regime'] != 'FUV' or d['method'] != 'synthesis':
            continue                       # FUV synthesis lines only (NUV NH 3360 is a coverage gap)
        lam = float(d['wavelength_A'])
        role = 'primary' if role_primary.get(d['species']) == lam else 'cross_check'
        out.append(Diagnostic(
            key=f"{d['species'].replace(' ', '')}_{int(round(lam))}_fuv",
            element=d['element'], kind='atomic',
            windows_A=((lam - 1.0, lam + 1.0),),     # vacuum window (FUV stays vacuum)
            use_molecules=False, role=role,
            nlte_flag='cI_fuv_grid_owed' if d['element'] == 'C' else 'fuv_grid_owed',
            nlte_ref=(f"{d['nlte_ref']}; GRID_OWED ~{d['nlte_expected_dex']:+.2f} dex — Amarsi "
                      f"C/O grid (RYA-359) not on disk, NOT applied (flagged, never silent LTE)"
                      if d.get('nlte_expected_dex') is not None else
                      f"{d['nlte_ref'] or 'RYA-190'}; GRID_OWED — Amarsi grid (RYA-359) owed"),
            reference=d.get('reference', 'RYA-190')))
    return tuple(out)


# ── smoke harness (real Procyon STIS; verifies vac->air on a known line) ──────
def _smoke(grating: str = 'E140M'):
    print(f"[hst_uv_loader] RYA-471 smoke — grating {grating}")
    # vac->air sanity on a KNOWN NUV line (Mg II k) + FUV-stays-vacuum on C I 1657
    mgk_air = float(uvc.to_pipeline_frame(np.array([2796.3543]))[0])
    ci_pf = float(uvc.to_pipeline_frame(np.array([1657.38]))[0])
    ok_air = abs(mgk_air - 2795.528) < 0.01
    ok_vac = abs(ci_pf - 1657.38) < 1e-6
    print(f"  [vac->air] Mg II k 2796.3543(vac) -> {mgk_air:.3f}(air) [known 2795.528]  {'OK' if ok_air else 'FAIL'}")
    print(f"  [boundary] C I 1657.38 stays vacuum -> {ci_pf:.3f}  {'OK' if ok_vac else 'FAIL'}")
    try:
        rows = science_ready_rows(grating=grating)
    except HstUvProductError as exc:
        print(f"  [data] {exc}\n  (vac->air verified above; data step skipped)")
        return ok_air and ok_vac
    if not rows:
        print(f"  [data] no whitelisted on-disk frame for {grating}; data step skipped")
        return ok_air and ok_vac
    spec = load_procyon_uv_arm(grating, return_spectrum=True)
    lo, hi = spec.wave_range_A
    print(f"  [load] {Path(spec.meta['filepath']).name}  targ={spec.meta['object']!r}")
    print(f"         range={lo:.2f}-{hi:.2f} A  npix={spec.n_pixels}  regime={spec.meta['regime']}  "
          f"frame={spec.meta['frame_out']}")
    print(f"         masked_frac={spec.meta['masked_fraction']}  analysis_ready={spec.meta['analysis_ready']}  "
          f"ew_allowed={spec.meta['ew_allowed']}  provenance={spec.meta['provenance']}")
    ci_covered = lo <= 1657.38 <= hi
    print(f"  [coverage] C I 1657.38 in arm: {ci_covered}  (the measured-C Procyon advantage)")
    diags = uv_arm_diagnostics()
    print(f"  [diagnostics] {len(diags)} FUV-synthesis lines wired from RYA-190: "
          f"{[d.key for d in diags]}")
    return ok_air and ok_vac and ci_covered


if __name__ == '__main__':
    import argparse
    import sys
    ap = argparse.ArgumentParser(description='RYA-471 HST UV (STIS/COS) loader')
    ap.add_argument('--smoke', action='store_true', help='end-to-end on real Procyon STIS')
    ap.add_argument('--grating', default='E140M')
    args = ap.parse_args()
    if args.smoke:
        sys.exit(0 if _smoke(args.grating) else 1)
    ap.print_help()

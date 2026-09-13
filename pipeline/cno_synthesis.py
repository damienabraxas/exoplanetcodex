"""
pipeline/cno_synthesis.py
=========================
Region-aware C / N / O synthesis engine — Turbospectrum flux-fit (RYA-237).

C/N/O is the flagship science (the C/O ratio drives the rocky-planet-composition
thesis and the 55 Cnc controversy). C, N and O are coupled through molecular
equilibrium (CO / CN / CH), so we **synthesize** rather than invert EWs:
Turbospectrum's equation of state solves the molecular partial pressures
internally and self-consistently at each set of A(C)/A(N)/A(O). There is no
hand-coded "CO correction" — the iteration is just re-fit + re-synthesize until
A(C)/A(N)/A(O) are self-consistent (the EOS recomputes equilibrium each call).

This module builds the **region-aware engine** and validates it on the
**HARPS-VIS arm** (the data we have now: solar Dumusque HARPS + Procyon HARPS).
The engine is built region-aware — per-instrument LSF, a telluric-correction
gate for IR arms, and per-region product output — so the UV (STIS) and IR
(CO-band) arms plug in without rework as their data clears the campaign gates
(RYA-351 / RYA-162 / RYA-119). Those arms are sequenced by data readiness, not
deferred.

Plugs into (already on main):
  * synth-v2 flux-fitting core — `abundances_derive._synth_flux_at_abund` /
    `_fit_synth_flux` (RYA-285/287); we reuse its resources + broadening.
  * single-source gf via `gf_resolver` (RYA-353) — the GES linelist loaded by
    `_load_synth_resources()` is already rescaled to canonical gf.
  * per-star broadening (RYA-288) — `_resolve_broadening`, fail-loud (no solar
    default leak into Procyon).
  * the sign-corrected NLTE module (RYA-339) — wired as a pluggable per-arm hook.

Molecular bands: iSpec's Turbospectrum wrapper auto-includes
`input/linelists/turbospectrum/molecules/*.bsyn` whenever `use_molecules=True`.
Coverage is band-specific (RYA-360 measured spans; molecules-dir README = Gerber
et al. 2023 / Masseron VALD compilation, 420–920 nm):
  * CH, CN, C2, OH, NH — held as **400–950 nm ELECTRONIC bands only** (measured
    spans ≈4200–9200 Å). Their mid-IR ro-vibrational fundamentals are NOT held
    (0 rows in the OH/NH/CH fundamental windows — confirmed by measurement,
    RYA-360/499); acquiring those is RYA-503.
  * CO — the mid-IR exception: `CO_IR_Li2015.dat` (ExoMol Li+2015, RYA-236) is a
    genuine mid-IR ro-vibrational list (spans NIR overtones through the ~4.6 µm
    fundamental).
All are vendored + guarded (RYA-360: data/linelists/molecular/turbospectrum/ +
the [molecular] stewardship invariant). Setting C/N/O as fixed abundances feeds
babsma's molecular equilibrium, so the band depths respond to A(C)/A(N)/A(O) —
that coupling is the whole point of synthesizing CNO.

NLTE — VIS arm is LTE-correct by physics (run scope, Ryan 2026-06-19; this is
correct treatment, NOT a silent fallback — see `VIS_NLTE_POLICY`):
  * [O I] 6300 — forbidden, LTE-insensitive.
  * CH G-band / C2 Swan / CN red — molecular bands (no molecular NLTE grid).
  * C I 5052/5380 — optical C I, ΔNLTE < 0.03 dex (Alexeeva & Mashonkina 2015);
    stamp `cI_vis_lte_assumed` so it is revisited when RYA-359 lands and the
    red-arm O I 777 cross-check runs.
The 339+359 post-hoc / in-synthesis NLTE path is exercised when the
red-optical / IR / UV arms land (O I 777 triplet, IR/FUV C I) — that is where
corrections are large and the grid is mandatory. The NLTE application is kept
pluggable per arm (`nlte_backend`) so those arms slot in without rework.

Linear issue: RYA-237   (Procyon-VIS shakedown: RYA-348)
"""

import os
import sys
import json
import time
import argparse
import warnings
from dataclasses import dataclass, field
from pathlib import Path

from pipeline import _runtime as _rt   # RYA-514: force-fork + single-thread BLAS (before numpy)
import numpy as np
from pipeline._numcompat import trapezoid as _trapezoid  # numpy>=2 removed np.trapz (RYA-313)
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.interpolate import interp1d

from config.constants import (
    PATHS, SOLAR_ASPLUND2021, HARPS_R, ISPEC_DIR, get_star_params, ROOT,
)
# Reuse the synth-v2 core (RYA-285/287/288/353) — same atmosphere interpolation,
# linelist (canonical-gf rescaled), isotopes, observed-spectrum loader and
# per-star broadening resolver the Fe synthesis path uses. No parallel machinery.
from pipeline.abundances_derive import (
    ispec,
    _load_atmosphere,
    _load_synth_resources,
    _load_observed_spectrum,
    _resolve_broadening,
    _ISPEC_SOLAR_ABUND_FILE,
)
from pipeline.gf_resolver import resolve as resolve_gf   # RYA-365: canonical Ni gf assert
from pipeline import nlte_cno   # RYA-359: vendored Amarsi 2019 C I / O I 3D-NLTE grid (Phase-A C I correction)
from pipeline.spectra_normalize import fit_continuum   # RYA-371: arm co-add continuum (same machinery as HARPS)
from pipeline.loaders import hst_uv_loader as _hst_uv   # RYA-471: HST STIS/COS UV arm loader + diagnostics

# Molecular line lists (RYA-236) — iSpec globs this dir when use_molecules=True.
_MOLECULES_DIR = ISPEC_DIR / 'input' / 'linelists' / 'turbospectrum' / 'molecules'

# χ²ᵣ fit-quality gate per band (the established synth convention, RYA-342).
from config.constants import SYNTH_CHI2_GATE


# ── Reflected-solar arm loader (RYA-371 Phase A multi-arm) ────────────────────
# Co-add + continuum-normalize the RYA-372 rest-frame conditioned ESPRESSO/UVES
# frames into ONE normalized spectrum per arm for synthesis. The 372 --write emits
# ONLY PASS frames (held-out rest-frame residual < 0.5 km/s); we re-assert rest-frame
# readiness from the manifest and LOUD-FAIL otherwise — the Phase-A "wiring check
# FIRST" rule (never synthesize against raw velocity-shifted spectra). FLUXCAL flux is
# rescaled to ~unity before fit_continuum (which floors the continuum at 1e-6; the raw
# ESPRESSO FLUXCAL flux is ~1e-12, so un-rescaled it normalizes to ~0).
_REFLECTED_DIR = Path(__file__).resolve().parents[1] / 'data' / 'processed' / 'reflected_solar'


def _load_reflected_solar_arm(inst: str, closure_tol_kms: float = 1.5) -> tuple:
    """RYA-372 rest-frame ESPRESSO/UVES → (wave_nm, flux_norm), co-added + normalized.
    Loud-fails if the rest-frame product is absent or no PASS frame is rest-frame
    verified (no silent use of unconditioned data)."""
    man = _REFLECTED_DIR / f'vesta_{inst}_manifest.csv'
    if not man.exists():
        raise FileNotFoundError(
            f"RYA-372 rest-frame manifest absent: {man}. Run `python -m "
            f"pipeline.reflected_solar_rv --set vesta_{inst} --write` first — Phase A "
            f"consumes the 372 rest-frame product, NOT raw velocity-shifted spectra.")
    mf = pd.read_csv(man)
    passed = mf[mf['status'].astype(str) == 'PASS']
    if passed.empty:
        raise RuntimeError(f"No PASS (rest-frame-verified) {inst} frame in {man.name}; "
                           f"refusing to synthesize against unconditioned data.")
    cl = pd.to_numeric(passed['closure_resid'], errors='coerce').abs()
    if cl.notna().any() and cl.max() > closure_tol_kms:
        raise RuntimeError(f"{inst}: max telluric-closure residual {cl.max():.2f} km/s "
                           f"> {closure_tol_kms} — rest frame NOT verified, refusing.")
    frames = []
    for _, r in passed.iterrows():
        p = _REFLECTED_DIR / str(r['output'])
        if not p.exists():
            continue
        d = pd.read_csv(p, comment='#')
        w = d['wavelength_air_A'].to_numpy(float)
        fl = d['flux'].to_numpy(float)
        m = np.isfinite(fl) & (fl > 0)
        if m.sum() > 1000:
            frames.append((w[m], fl[m]))
    if not frames:
        raise RuntimeError(f"{inst}: PASS manifest rows present but no readable rest-frame files.")
    lo = min(w.min() for w, _ in frames)
    hi = max(w.max() for w, _ in frames)
    grid = np.arange(lo, hi, 0.02)
    stack = np.array([np.interp(grid, w, fl, left=np.nan, right=np.nan) for w, fl in frames])
    coadd = np.nanmedian(stack, axis=0)
    ok = np.isfinite(coadd) & (coadd > 0)
    g, c = grid[ok], coadd[ok]
    c = c / np.nanmedian(c)              # FLUXCAL → ~unity (fit_continuum 1e-6 floor)
    cont = fit_continuum(g, c)
    print(f"  [arm-load] {inst}: co-added {len(frames)} PASS frame(s) over "
          f"{g.min():.0f}-{g.max():.0f} A; continuum-normalized (median "
          f"{np.median(c / cont):.3f}); max closure {cl.max():.2f} km/s")
    return g / 10.0, (c / cont)


def _load_procyon_uves_arm() -> tuple:
    """RYA-348 Phase 2 — Procyon UVES red arm → (wave_nm, flux_norm), the O I 777 PRIMARY-O
    arm. Resolves the RYA-272 registry's oi_anchor epoch (the single CLEAN telluric-verdict
    frame), loads it through the RYA-272 UVESLoader (BERV applied → barycentric rest frame,
    air Angstrom, GES-quarantine guards), then continuum-normalizes the same way the
    reflected-solar arm does. NEVER reflected-solar Vesta (RYA-464): a real Procyon UVES
    spectrum or a loud failure. Telluric: the anchor is the CLEAN epoch (RYA-271/272 audit);
    EXCLUDE/CORRECTABLE epochs are never selected here."""
    from pipeline.loaders.uves_loader import UVESLoader
    from pipeline.abundances_derive import _measure_rv_kms
    reg_path = ROOT / 'data' / 'spectra' / 'procyon' / 'uves_registry.csv'
    if not reg_path.exists():
        raise FileNotFoundError(
            f"RYA-272 UVES registry absent: {reg_path}. Build it (pipeline UVES intake) "
            f"before resolving the Procyon UVES arm — no silent Vesta fallback (RYA-464).")
    reg = pd.read_csv(reg_path)
    anchor = reg[(reg['oi_anchor'].astype(str).str.lower() == 'true') &
                 (reg['oi_telluric_verdict'].astype(str).str.upper() == 'CLEAN')]
    if anchor.empty:
        raise RuntimeError(
            f"No CLEAN oi_anchor epoch in {reg_path.name} — the O I 777 primary-O arm needs a "
            f"telluric-CLEAN frame (RYA-271). EXCLUDE/CORRECTABLE epochs are not auto-used.")
    row = anchor.iloc[0]
    uves_dir = ROOT.parent / 'data' / 'spectra' / 'exoplanetcodex-data' / 'Procyon' / 'Procyon UVES'
    fpath = uves_dir / str(row['filename'])
    if not fpath.exists():
        raise FileNotFoundError(
            f"Procyon UVES anchor staged in registry but absent on disk: {fpath}. Stage the "
            f"RYA-272 UVES frames before the run (data store, gitignored).")
    spec = UVESLoader(fpath).load()
    obj = str(spec.meta.get('object', '')).lower()
    if 'procyon' not in obj and 'vesta' in obj:        # anti-silent-Vesta belt-and-suspenders
        raise RuntimeError(f"UVES anchor object={spec.meta.get('object')!r} is reflected-solar — "
                           f"refusing to synthesize Procyon against Vesta (RYA-464).")
    w = np.asarray(spec.wave_A, float)
    fl = np.asarray(spec.flux, float)
    m = np.isfinite(w) & np.isfinite(fl) & (fl > 0)
    w, fl = w[m], fl[m]
    fl = fl / np.nanmedian(fl)                          # FLUXCAL → ~unity (fit_continuum floor)
    cont = fit_continuum(w, fl)
    fln = fl / cont
    # Bring to the STELLAR REST frame (RYA-478): UVESLoader applies BERV (-> barycentric) but
    # NOT the stellar systemic RV, so the lines sit blueshifted ~0.12 A at 7773 (Procyon RV).
    # The synthesis is rest-frame; an uncorrected RV misaligns line cores and biases the
    # flux-space fit (the χ²ᵣ-147 / inflated-A(O) bug). Measure from Fe I centroids and shift
    # (the same correction _load_observed_spectrum applies to the HARPS arm, RYA-309).
    rv = _measure_rv_kms(w / 10.0, fln)
    w_rest = w / (1.0 + rv / _C_KMS) if abs(rv) > 0.05 else w
    print(f"  [arm-load] UVES Procyon: anchor {row['filename']} ({row['date_obs']} "
          f"{row['setting']}, SNR {float(row['snr']):.0f}, telluric {row['oi_telluric_verdict']}); "
          f"{w.min():.0f}-{w.max():.0f} A, BERV {spec.meta['berv_kms']:+.2f} km/s applied, "
          f"stellar RV {rv:+.2f} km/s -> rest; continuum-normalized (median {np.nanmedian(fln):.3f})")
    return w_rest / 10.0, fln


_C_KMS = 299792.458


# ── Region awareness ──────────────────────────────────────────────────────────
# Built now, lit per arm as data clears (RYA-351/162/119). Each arm carries its
# instrument LSF, a telluric-clearance requirement, and the NLTE backend the arm
# routes through. Products are emitted PER REGION (never coadded across
# resolutions — combined at the abundance layer, RYA-282 presentation decision).

@dataclass(frozen=True)
class RegionConfig:
    name: str                       # 'vis' | 'red' | 'ir' | 'uv'
    instrument: str                 # 'HARPS' | 'UVES' | 'CRIRES+' | 'STIS' ...
    R: float                        # instrumental resolving power (LSF)
    wave_min_A: float
    wave_max_A: float
    telluric_correction_required: bool
    nlte_backend: str               # 'lte_by_design' (VIS) | 'amarsi_grid' | 'in_synthesis'
    notes: str = ''


HARPS_VIS = RegionConfig(
    name='vis', instrument='HARPS', R=float(HARPS_R),
    wave_min_A=3780.0, wave_max_A=6910.0,
    telluric_correction_required=False,   # optical HARPS; the chosen VIS windows
                                          # are free of significant telluric bands
    nlte_backend='lte_by_design',         # VIS CNO lines are LTE-correct (see module docstring)
    notes='HARPS 3780-6910 Angstrom; solar Dumusque + Procyon. RYA-237 deliverable.',
)

#: 🔴 RYA-1214 — THE NEAR-UV MOLECULAR REGION. Ryan, 2026-09-12: *"PRIORITY for CNO
#: molecular: near-UV and VIS are where the molecular bands are strongest ... OH A-X
#: 306-330 nm."*
#:
#: Of the four bands that comment names, TWO are inside 3000-3780 A: OH A-X (3060-3300)
#: and NH A-X (~3360). CN violet (~3883) and CH A-X (3870-4320) are redward of the edge —
#: and CH A-X IS the VIS G-band. So this region is specifically an OXYGEN and NITROGEN
#: opportunity, and the only route to a molecular A(O) or A(N) that needs no region above
#: 1 um.
#:
#: ⚠️ THE PRODUCTS OF THIS REGION ARE UPPER BOUNDS, NOT ABUNDANCES, AND THE REASON IS
#: MEASURED. The near-UV synthesis carries LESS opacity than the observation: RYA-1204/1207
#: measured the molecular lever at -0.050 dex (paired) on iron and concluded it STILL
#: under-corrects afterwards ("synth/obs above 1 on 24 of 40 windows"; "a floor on the
#: missing opacity, not a correction stopped at a target"), and RYA-1189 measured the band
#: blend-dominated at 59 lines/A with 0 of 10 clean side-bands. A fit compensates missing
#: opacity by RAISING the abundance, so the bias direction is known IN ADVANCE. Carried as
#: `opacity_deficit_upper_bound` on every result this region emits.
NEARUV_KP = RegionConfig(
    name='nearuv', instrument='kpno_solar_atlas', R=500000.0,
    wave_min_A=3000.0, wave_max_A=3780.0,
    telluric_correction_required=False,   # Kurucz 2005 is telluric-corrected at source
                                          # (RYA-933); nothing terrestrial reaches 3000-3780
    nlte_backend='lte_by_design',         # molecular bands: no molecular NLTE grid exists
    notes='Kitt Peak / Kurucz 2005 residual atlas 3000-3780 A. OH A-X + NH A-X molecular '
          'bands. RYA-1214. Products are UPPER BOUNDS — see the opacity-deficit note.',
)

#: Half-width of an IR CN fit window, in Angstrom. Wider than the near-UV's 1.5 because the
#: regime is 300x less crowded (0.2 atomic lines/A against 59) and narrower than the
#: red-optical's 1.1 because these lines are 1-20 mA and a wide window buys blend, not signal.
_CN_IR_PAD_A = 1.5
#: Compute cap on how many windows one diagnostic fits. `_fit_element` synthesises EVERY
#: window on EVERY chi2 evaluation, so cost is linear in total window width: the near-UV's
#: 2x3 A diagnostic took ~120 s, so 48x3 A would be ~48 min per fit and ~4 h over 5
#: equilibrium iterations. Capped, and the cap is REPORTED with the count it dropped —
#: a silently truncated line set is a different measurement (RYA-842).
_CN_IR_MAX_WINDOWS = 12


def _cn_ir_windows():
    """AGSS21's OWN CN lines, clustered into fit windows. Returns (kp_windows, iag_windows).

    🔴 THE WINDOWS ARE AGSS21's LINE POSITIONS, NOT POSITIONS I CHOSE. Clustering their 59
    published CN wavelengths and padding by `_CN_IR_PAD_A` guarantees the fit measures the
    features they measured — the same principle as `--lines-from-set` for atomic lines, and
    the reason this is a replication rather than a new selection. `linelist_solar` has no rows
    at 1.1-1.3 um to dominance-test against, so the source's own selection IS the rule here.

    ⚠️ THE H2O BAND IS EXCLUDED, NOT CORRECTED. `telluric_policy.TELLURIC_BANDS` puts H2O at
    11120-11560 A, which swallows 5 of the 59 lines. Both holdings are telluric-corrected and
    could be argued into fitting there; 440 A of the band's 2333 is not worth leaning on a
    correction for when 1893 A is enumerated-clean.

    Selection among what is left is by AGSS21's OWN PUBLISHED EQUIVALENT WIDTH, descending —
    a source-published strength ranking, not a quantity of ours. Capped at
    `_CN_IR_MAX_WINDOWS`.
    """
    import csv as _csv
    from pathlib import Path as _P
    src = (_P(__file__).resolve().parents[1] / 'data' / 'reference' / 'amarsi2021_cno'
           / 'derived' / 'amarsi2021_cno_molecular_lines.csv')
    rows = []
    with src.open(newline='') as fh:
        for r in _csv.DictReader(fh):
            if r['element_parameter'] != 'logepsN' or r['species'] != 'CN':
                continue
            lv = float(r['wavelength_vac_nm']) * 10.0
            s2 = (1e4 / lv) ** 2
            air = lv / (1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2))
            rows.append((air, float(r['equivalent_width_pm']) * 10.0))
    if not rows:
        raise RuntimeError(
            "no CN rows in the Amarsi 2021 holding — the IR CN windows ARE its line "
            "positions, so an empty read must not silently become an empty region.")
    rows.sort()
    H2O_LO, H2O_HI = 11120.0, 11560.0
    clusters, cur = [], [rows[0]]
    for w, ew in rows[1:]:
        if w - cur[-1][0] <= 2 * _CN_IR_PAD_A:
            cur.append((w, ew))
        else:
            clusters.append(cur); cur = [(w, ew)]
    clusters.append(cur)
    cand = []
    for c in clusters:
        lo = min(x[0] for x in c) - _CN_IR_PAD_A
        hi = max(x[0] for x in c) + _CN_IR_PAD_A
        if hi > H2O_LO and lo < H2O_HI:
            continue                               # inside H2O — excluded, not corrected
        cand.append((round(lo, 2), round(hi, 2), max(x[1] for x in c)))
    kp = sorted(sorted(cand, key=lambda t: -t[2])[:_CN_IR_MAX_WINDOWS])
    iag_all = [c for c in cand if c[1] <= 11083.0]
    iag = sorted(sorted(iag_all, key=lambda t: -t[2])[:_CN_IR_MAX_WINDOWS])
    return (tuple((a, b) for a, b, _ in kp), tuple((a, b) for a, b, _ in iag))


#: 🔴 RYA-1214 — THE INFRARED CN REGIONS. Ryan: *"visible N is going to be hard, hence the
#: IR."* This is where AGSS21's nitrogen actually lives.
#:
#: AGSS21's CN indicator is the A-X band at 10872-13204 A (Amarsi et al. 2021 Table 2, all 59
#: CN lines in band "(0-0)"). The optical `CN_red` diagnostic above fits 6125/6195 A — the CN
#: A-X RED system — and ZERO of AGSS21's 59 lines fall in it. Same molecule, different band.
#:
#: Two regions because two holdings reach different amounts of it, MEASURED by probing each
#: holding rather than reading declared spans (three of them declare none):
#:
#:   NIR_CN_KP    solar_kpno_molecfit_corrected   CLEAN   the WHOLE band   54 clean lines
#:   NIR_CN_IAG   solar_iag                       CLEAN   to 11083 A       16 clean lines
#:
#: ⚠️ CRIRES+ CANNOT SERVE THIS BAND AT ALL, and that is a measurement. Its Y arm ends at
#: 10796 A — 76 A blueward of AGSS21's first CN line — and its H arm starts at 15007 A, 1803 A
#: redward of the last. The raw Vesta IDPs reach further but are BLOCKED
#: (RestFrameNotConditioned), not merely uncorrected. The instrument catalogue says CRIRES+
#: spans 9500-53000 A, which is exactly why this was probed per HOLDING.
#:
#: WHY THE IR IS THE RIGHT PLACE, in one number: the NIR atomic list carries 757 lines over
#: 3800 A (0.2 lines/A) against the near-UV's 59 lines/A. The blend contamination that makes
#: the optical N I lines unmeasurable — N I carrying 0.9-5.9% of its own fit window — is a
#: property of a crowded band, and this band is not crowded.
_CN_IR_WINDOWS_KP, _CN_IR_WINDOWS_IAG = _cn_ir_windows()

NIR_CN_KP = RegionConfig(
    name='nir_cn_kp', instrument='kpno_solar_atlas', R=500000.0,
    wave_min_A=10872.0, wave_max_A=13205.0,
    telluric_correction_required=True,    # H2O 11120-11560 sits inside the band; the fit
                                          # windows AVOID it rather than lean on a correction
    nlte_backend='lte_by_design',         # molecular band: no molecular NLTE grid exists
    notes='Kitt Peak 1984 composite, molecfit-corrected. AGSS21 CN A-X 10872-13204 A. '
          'Fit windows are clustered around AGSS21 OWN 59 CN lines, H2O 11120-11560 excluded.',
)
NIR_CN_IAG = RegionConfig(
    name='nir_cn_iag', instrument='iag_fts_solar_atlas', R=700000.0,
    wave_min_A=10872.0, wave_max_A=11083.0,
    telluric_correction_required=True,
    nlte_backend='lte_by_design',
    notes='IAG FTS to its own red edge 11083.46 A. 16 of AGSS21 59 CN lines, all telluric-'
          'clean (nearest band H2O starts at 11120 A, 37 A redward of the edge).',
)

REGIONS = {'vis': HARPS_VIS, 'nearuv': NEARUV_KP,
           'nir_cn_kp': NIR_CN_KP, 'nir_cn_iag': NIR_CN_IAG}


# ── Diagnostic registry (HARPS-VIS, wavelength-correct) ───────────────────────
# Windows are AIR wavelengths in Angstrom, the fit sub-windows for each band.
# `depends_on` lists elements whose abundance must be fixed before this diagnostic
# is fit (the equilibrium coupling — CN needs A(C); [O I] is tied to A(C) via CO).

@dataclass(frozen=True)
class Diagnostic:
    key: str
    element: str                    # 'C' | 'N' | 'O'
    kind: str                       # 'molecular_band' | 'atomic' | 'forbidden_blend'
    windows_A: tuple                # tuple of (lo, hi) Angstrom fit sub-windows
    use_molecules: bool
    role: str                       # 'primary' | 'cross_check'
    nlte_flag: str
    nlte_ref: str
    depends_on: tuple = ()
    pinned_blends: tuple = ()        # elements pinned from EW/canonical (e.g. Ni for [O I])
    reference: str = ''


VIS_DIAGNOSTICS = (
    # ── Carbon ────────────────────────────────────────────────────────────────
    Diagnostic(
        key='CH_Gband', element='C', kind='molecular_band',
        windows_A=((4303.5, 4306.5), (4310.0, 4313.0)),
        use_molecules=True, role='primary',
        nlte_flag='lte_molecular_band',
        nlte_ref='molecular band — no NLTE grid (LTE-by-design)',
        reference='CH A-X G-band 4290-4315 (Masseron+2014/2022); primary solar/Procyon A(C)',
    ),
    Diagnostic(
        key='CI_5052', element='C', kind='atomic',
        windows_A=((5051.3, 5053.0),),
        use_molecules=False, role='cross_check',
        nlte_flag='cI_vis_lte_assumed',
        nlte_ref='Alexeeva & Mashonkina 2015: optical C I ΔNLTE < 0.03 dex',
        reference='C I 5052.17 atomic; LTE cross-check vs CH (revisit when RYA-359 lands)',
    ),
    Diagnostic(
        key='CI_5380', element='C', kind='atomic',
        windows_A=((5379.3, 5381.3),),
        use_molecules=False, role='cross_check',
        nlte_flag='cI_vis_lte_assumed',
        nlte_ref='Alexeeva & Mashonkina 2015: optical C I ΔNLTE < 0.03 dex',
        reference='C I 5380.34 atomic; LTE cross-check vs CH',
    ),
    Diagnostic(
        key='C2_Swan', element='C', kind='molecular_band',
        windows_A=((5160.0, 5166.0),),
        use_molecules=True, role='cross_check',
        nlte_flag='lte_molecular_band',
        nlte_ref='molecular band — no NLTE grid (LTE-by-design)',
        reference='C2 Swan (0,0) bandhead 5165 (metal-rich cross-check, 55 Cnc)',
    ),
    # ── Nitrogen ──────────────────────────────────────────────────────────────
    Diagnostic(
        key='CN_red', element='N', kind='molecular_band',
        windows_A=((6125.0, 6130.0), (6195.0, 6200.0)),
        use_molecules=True, role='primary', depends_on=('C',),
        nlte_flag='lte_molecular_band',
        nlte_ref='molecular band — no NLTE grid (LTE-by-design)',
        reference='CN red A-X 6000-6200 (Brooke+Sneden 2014); primary N, given A(C). '
                  'Also deblends Li 6707 (RYA-103).',
    ),
    # ── Oxygen ────────────────────────────────────────────────────────────────
    Diagnostic(
        key='OI_6300', element='O', kind='forbidden_blend',
        windows_A=((6299.5, 6301.0),),
        use_molecules=True, role='primary', depends_on=('C',),
        pinned_blends=('Ni',),
        nlte_flag='lte_forbidden_insensitive',
        nlte_ref='[O I] forbidden — LTE-insensitive',
        reference='[O I] 6300.30 + Ni I 6300.34 joint synthesis (A(Ni) pinned). '
                  'Ni I 6300.34 gf resolves via gf_resolver = Johansson+2003 '
                  'log gf -2.11 (RYA-365 adjudication; the [O I]-blend lab authority).',
    ),
)


# ── Near-UV molecular diagnostics (Kitt Peak 3000-3780 A) — RYA-1214 ──────────
#
# 🔴 THE WINDOWS ARE CHOSEN BY THE TARGET MOLECULE'S DOMINANCE, MEASURED, NOT BY EYE.
# A band fit is a measurement of X only if X's lines carry the window. Counted on RYA-1207's
# own near-UV .bsyn lists, per 3 A window (the width VIS_DIAGNOSTICS uses), with every
# candidate scored whether it won or lost -- `scripts/rya1214_nearuv_molecular_sizing.py`:
#
#     ADOPTED  OH 3063-3066   OH  30  NH 10  CH  0  CN  0   dominance 3.0
#     ADOPTED  OH 3122-3125   OH  26  NH  7  CH  6  CN  0   dominance 2.0
#     rejected OH 3144-3147   OH  19  NH  3  CH 48  CN 14   dominance 0.29  CH-dominated
#     rejected OH 3180-3183   OH  10  NH  2  CH 26  CN  0   dominance 0.36  CH-dominated
#     ADOPTED  NH 3358-3361   NH 109  OH 11  CH  0  CN  6   dominance 6.41
#     ADOPTED  NH 3370-3373   NH  75  OH  8  CH  0  CN  0   dominance 9.38
#
# Two OH candidates were REFUSED on that measurement, and they are the two a wider-is-better
# instinct would have taken: 3144 and 3180 hold more total molecular lines than either
# adopted OH window and are dominated by CH, so a fit there would move A(O) to fit carbon.
#
# ⚠️ OH's dominance (2-3x) is far weaker than NH's (6-9x). The near-UV OH A-X region is
# genuinely contaminated by NH and CH, and A(C)/A(N) are pinned before O is fit for exactly
# that reason -- but it is the reason the OH bound is the looser of the two.
NEARUV_DIAGNOSTICS = (
    Diagnostic(
        key='NH_AX', element='N', kind='molecular_band',
        windows_A=((3358.0, 3361.0), (3370.0, 3373.0)),
        use_molecules=True, role='primary', depends_on=('C',),
        nlte_flag='lte_molecular_band',
        nlte_ref='molecular band — no NLTE grid (LTE-by-design)',
        reference='NH A-X (0,0) band head 3360 A; near-UV molecular N. UPPER BOUND — the '
                  'near-UV synthesis under-corrects opacity (RYA-1204/1207/1189).',
    ),
    Diagnostic(
        key='OH_AX', element='O', kind='molecular_band',
        windows_A=((3063.0, 3066.0), (3122.0, 3125.0)),
        use_molecules=True, role='primary', depends_on=('C', 'N'),
        nlte_flag='lte_molecular_band',
        nlte_ref='molecular band — no NLTE grid (LTE-by-design)',
        reference='OH A-X (0,0) 3064 A + (1,1) 3123 A; near-UV molecular O. UPPER BOUND, '
                  'and the looser of the two: OH dominance here is 2-3x against NH 6-9x.',
    ),
)


# ── ESPRESSO + UVES arms (RYA-371 Phase A multi-arm; reflected-solar, RYA-372) ─
# These consume the RYA-372 rest-frame conditioned co-add (_load_reflected_solar_arm),
# NOT the HARPS normalized spectrum. Per the RYA-455 O-handling amendment, O I 777
# (ESPRESSO) is the PRIMARY O; [O I] 6300 is a continuum-limited cross-check only.
ESPRESSO_OPT = RegionConfig(
    name='espresso', instrument='ESPRESSO', R=140000.0,
    wave_min_A=3780.0, wave_max_A=7890.0, telluric_correction_required=False,
    nlte_backend='amarsi_grid',
    notes='ESPRESSO reflected-solar Vesta (RYA-370/372); O I 777 primary O + [O I] 6300 '
          'cross-check + C I. O I 777 sits in an H2O region (RYA-373 optical telluric owed).',
)
UVES_OPT = RegionConfig(
    name='uves', instrument='UVES', R=70000.0,
    wave_min_A=3760.0, wave_max_A=9470.0, telluric_correction_required=False,
    nlte_backend='amarsi_grid',
    notes='UVES reflected-solar Vesta (RYA-370/372); N from N I 8216 + CN red. NH 3360 is '
          'a DATA-GAP (UVES-346 blue under SNR floor, did not condition — RYA-369). '
          'R is a representative value across the co-added dichroic settings.',
)

ESPRESSO_DIAGNOSTICS = (
    Diagnostic(
        key='OI_777', element='O', kind='atomic',
        windows_A=((7771.0, 7772.9), (7773.2, 7775.0), (7774.6, 7776.3)),
        use_molecules=False, role='primary',
        nlte_flag='oI_777_amarsi2019',
        nlte_ref='O I 777 triplet 3D-NLTE — Amarsi 2019 (RYA-359); large negative',
        reference='O I 7771.94/7774.17/7775.39 (ESPRESSO); PRIMARY solar A(O) (RYA-455).',
    ),
    Diagnostic(
        key='OI_6300', element='O', kind='forbidden_blend',
        windows_A=((6299.5, 6301.0),), use_molecules=True, role='cross_check',
        pinned_blends=('Ni',),
        nlte_flag='lte_forbidden_continuum_limited',
        nlte_ref='[O I] forbidden; continuum-limited (RYA-447->455); adopt Caffau 2015 8.73',
        reference='[O I] 6300.30 + Ni I 6300.34 (ESPRESSO); CONTINUUM-LIMITED cross-check.',
    ),
    Diagnostic(
        key='CI_5052', element='C', kind='atomic', windows_A=((5051.3, 5053.0),),
        use_molecules=False, role='cross_check',
        nlte_flag='cI_amarsi2019', nlte_ref='C I 3D-NLTE — Amarsi 2019 (RYA-359)',
        reference='C I 5052.17 (ESPRESSO); C cross-check vs HARPS.',
    ),
    Diagnostic(
        key='CI_5380', element='C', kind='atomic', windows_A=((5379.3, 5381.3),),
        use_molecules=False, role='cross_check',
        nlte_flag='cI_amarsi2019', nlte_ref='C I 3D-NLTE — Amarsi 2019 (RYA-359)',
        reference='C I 5380.34 (ESPRESSO); C cross-check vs HARPS.',
    ),
)

UVES_DIAGNOSTICS = (
    Diagnostic(
        key='NI_8216', element='N', kind='atomic', windows_A=((8216.0, 8216.7),),
        use_molecules=False, role='primary',
        nlte_flag='nI_lte_flagged',
        nlte_ref='N I 8216; no NLTE grid wired — LTE-flagged (NLTE owed, RYA-369)',
        reference='N I 8216.34 (UVES red); best available N (NH 3360 DATA-GAP, RYA-369).',
    ),
    Diagnostic(
        key='CN_red', element='N', kind='molecular_band',
        windows_A=((6125.0, 6130.0), (6195.0, 6200.0)),
        use_molecules=True, role='cross_check', depends_on=('C',),
        nlte_flag='lte_molecular_band',
        nlte_ref='molecular band — no NLTE grid (LTE)',
        reference='CN red A-X (UVES); N cross-check given A(C).',
    ),
)


# ── UV (HST STIS) + IR (CRIRES+/SPIRou) arms — DECLARED, wired as loaders+data land ──
# RYA-464: declared so the per-star registry can carry them; ready=False until their
# loader (RYA-184 loaders/ package) + staged+audited data exist. Wavelength spans are the
# factual instrument ranges (RYA-351 coverage); diagnostics are populated by the audits
# (UV RYA-262, IR RYA-425) — NOT fabricated here (no invented UV/IR line data).
HST_UV = RegionConfig(
    name='uv', instrument='STIS', R=114000.0,
    wave_min_A=1150.0, wave_max_A=3200.0, telluric_correction_required=False,
    nlte_backend='amarsi_grid',
    notes='HST STIS/COS UV (RYA-222 data in hand; RYA-262 audit). FUV C I / O I / S I '
          '= MEASURED C/O/S (Procyon advantage over the solar cited composite). Loader BUILT '
          '(RYA-471, pipeline.loaders.hst_uv_loader); synthesis gated on the Amarsi NLTE grid '
          '(RYA-359) + FUV pseudo-continuum (RYA-426 gate 5).',
)
CRIRES_IR = RegionConfig(
    name='ir', instrument='CRIRES+', R=86000.0,
    wave_min_A=15000.0, wave_max_A=24000.0, telluric_correction_required=True,
    nlte_backend='lte_by_design',
    notes='IR CO first-overtone 2.3um + OH/CN (C cross-check + 12C/13C). TELLURIC-GATED '
          '(cr2res+molecfit / APERO+Wapiti, RYA-373). APOGEE H-band = weak-CO only (RYA-351).',
)

# Procyon UVES red optical diagnostic set: O I 777 (PRIMARY O) + [O I] 6300 cross-check +
# C I 5052/5380 + N I 8216 + CN red — composed by REUSING the existing Diagnostic objects
# (ESPRESSO red-optical C/O set + the UVES N set), no new line data (RYA-464 reuse rule).
#: Region -> its own diagnostic set. `run_cno` reads THIS rather than VIS_DIAGNOSTICS,
#: which is what makes a second region expressible at all (RYA-1214).
NIR_CN_KP_DIAGNOSTICS = (
    Diagnostic(
        key='CN_AX_IR', element='N', kind='molecular_band',
        windows_A=_CN_IR_WINDOWS_KP,
        use_molecules=True, role='primary', depends_on=('C',),
        nlte_flag='lte_molecular_band',
        nlte_ref='molecular band — no NLTE grid (LTE-by-design)',
        reference='CN A-X 10872-13204 A (Brooke+2014 gf, validated against Amarsi 2021 '
                  'Table 2 to a median 0.0064 dex). THE band AGSS21 nitrogen rests on; '
                  'windows are AGSS21 own line positions, H2O 11120-11560 excluded.',
    ),
)
NIR_CN_IAG_DIAGNOSTICS = (
    Diagnostic(
        key='CN_AX_IR', element='N', kind='molecular_band',
        windows_A=_CN_IR_WINDOWS_IAG,
        use_molecules=True, role='primary', depends_on=('C',),
        nlte_flag='lte_molecular_band',
        nlte_ref='molecular band — no NLTE grid (LTE-by-design)',
        reference='CN A-X 10872-11083 A on the IAG FTS, to its own red edge. 16 of AGSS21 '
                  '59 CN lines, all telluric-clean.',
    ),
)

REGION_DIAGNOSTICS = {'vis': VIS_DIAGNOSTICS, 'nearuv': NEARUV_DIAGNOSTICS,
                      'nir_cn_kp': NIR_CN_KP_DIAGNOSTICS,
                      'nir_cn_iag': NIR_CN_IAG_DIAGNOSTICS}


def primary_by_element(diagnostics) -> dict:
    """{element: the one Diagnostic whose role is 'primary'} — RYA-1214.

    REFUSES an ambiguous set rather than picking. Two primaries for one element in one
    region means the region does not state which measurement IS the product, and taking
    the first would make that choice silently — the RYA-780 rule applied to diagnostics
    instead of to lines. An element with NO primary is simply absent from the mapping;
    the caller pins it and records that it did.
    """
    out: dict = {}
    for d in diagnostics:
        if d.role != 'primary':
            continue
        if d.element in out:
            raise ValueError(
                f"two primary diagnostics for {d.element} in one region "
                f"({out[d.element].key!r} and {d.key!r}). A region must state which "
                f"measurement IS the product; refusing to pick by order.")
        out[d.element] = d
    return out


def holding_for_region(region: RegionConfig) -> str:
    """Which HOLDING serves this region — ONE decision, read by the loader AND the gate.

    🔴 WHICH KITT PEAK HOLDING, DECIDED BY WHAT THE BAND NEEDS, AND THEY ARE NOT
    INTERCHANGEABLE. `solar_kpno_kurucz2005_corrected` is the residual atlas and stops at
    10010 A; the AGSS21 CN band starts at 10872. The 1984 composite
    (`solar_kpno_molecfit_corrected`) reaches 13204 and is the only CLEAN Kitt Peak holding
    that does — `solar_kpno` reaches it too but is CONTROL_ONLY, not a science basis
    (RYA-1026). Decided on the band's own edges, never defaulted: serving the near-UV from
    the composite, or the IR from the residual atlas, would each be a product naming a
    spectrum it was not measured on (RYA-904).

    Factored out because the telluric gate has to ask about the SAME holding the loader will
    read. Two copies of this choice is how a gate clears one spectrum and a fit measures
    another (RYA-845).
    """
    if region.instrument == 'kpno_solar_atlas':
        return ('solar_kpno_kurucz2005_corrected' if region.wave_max_A <= 10010.0
                else 'solar_kpno_molecfit_corrected')
    if region.instrument == 'iag_fts_solar_atlas':
        return 'solar_iag'
    if region.instrument.lower().startswith('harps'):
        return 'solar_harps_molecfit_corrected'
    raise ArmNotWired(
        f"no holding declared for region {region.name!r} on {region.instrument!r}")


def _load_region_spectrum(star_id: str, region: RegionConfig):
    """The observed spectrum for (star, region) -> (wave_nm, flux).

    🔴 RYA-1214 — `run_cno` called `_load_observed_spectrum(star_id)`, the HARPS loader,
    for every region. So `--region nearuv` would have fitted 3000-3780 A windows against a
    3780-6910 A spectrum: `_fit_element` slices per window, finds fewer than 5 pixels in
    each, and returns status='failed' with 'no observed pixels in windows'. Loud, but for
    the wrong reason — it would have read as "the near-UV has no data" when the near-UV
    atlas is right there.
    """
    if region.instrument.lower().startswith('harps'):
        return _load_observed_spectrum(star_id)
    if region.instrument == 'kpno_solar_atlas':
        return _load_kpno_atlas_arm(star_id, region)
    if region.instrument == 'iag_fts_solar_atlas':
        return _load_generic_atlas_arm(star_id, region, holding_for_region(region))
    raise ArmNotWired(
        f"no spectrum loader for region {region.name!r} on instrument "
        f"{region.instrument!r}. Wire one in `_load_region_spectrum` rather than letting "
        f"it fall through to the HARPS loader (RYA-1214).")


def _load_generic_atlas_arm(star_id: str, region: RegionConfig, holding: str):
    """A named SOLAR atlas holding over a region's band -> (wave_nm, flux) — RYA-1214.

    One reader for every atlas region, so the near-UV OH/NH run, the IR CN runs and the
    atomic band-product route all see byte-identical flux from a given holding. A second
    reader for one holding is how two products of one spectrum drift apart (RYA-845).

    ⚠️ `load_window_ex(instrument, CENTRE, PAD)` takes a centre and a HALF-WIDTH, not
    (lo, hi). Passing the band edges asks for `10872 +/- 13205 A` and the coverage check
    refuses it — loudly, which is the only reason a units slip in an argument pair was a
    two-minute fix rather than a wrong spectrum.

    ⚠️ NO CONTINUUM IS FITTED. These holdings ship their own, and `prenormalised_guard`
    refuses a second one: placing one on Kitt Peak tilted a band 4% blue-to-red and cost
    0.0218 dex (RYA-933/1026).
    """
    if 'solar' not in star_id.lower() and 'sun' not in star_id.lower():
        raise ArmNotWired(
            f"{star_id}: {holding} is a SOLAR atlas. Refusing to synthesize {star_id} "
            f"against it (RYA-464's no-silent-substitution rule).")
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
    from measure_band_ew import load_window_ex                      # noqa: E402
    centre = 0.5 * (region.wave_min_A + region.wave_max_A)
    pad = 0.5 * (region.wave_max_A - region.wave_min_A)
    win = load_window_ex(region.instrument, centre, pad, holding=holding)
    w_A = np.asarray(win.wave, float)
    flux = np.asarray(win.flux, float)
    m = np.isfinite(w_A) & np.isfinite(flux) & (flux > 0)
    if int(m.sum()) < 1000:
        raise ArmNotWired(
            f"{holding} returned only {int(m.sum())} usable px over "
            f"{region.wave_min_A}-{region.wave_max_A} A — refusing to fit a band on a "
            f"spectrum that is mostly absent.")
    print(f"  [arm-load] {region.instrument} / {win.holding.holding_id}: "
          f"{int(m.sum())} px over {w_A[m].min():.1f}-{w_A[m].max():.1f} A "
          f"(pre-normalised at source — no continuum fitted, RYA-1026)")
    return w_A[m] / 10.0, flux[m]


def _load_kpno_atlas_arm(star_id: str, region: RegionConfig):
    """The Kitt Peak / Kurucz-2005 residual atlas over a region's band -> (wave_nm, flux).

    Reads through `measure_band_ew.load_window_ex`, the SAME loader the band-product route
    uses, so the near-UV molecular run and the near-UV atomic products see byte-identical
    flux. A second reader for one holding is how two products of one spectrum drift
    (RYA-845), and this holding in particular has a history: it was served from the wrong
    file for months because `0irrad.readme` does not list the residual atlas (RYA-933).

    ⚠️ NO CONTINUUM IS FITTED HERE. `solar_kpno_kurucz2005_corrected` ships its own
    (`irradrelwl.dat` is the residual atlas), and `pipeline.prenormalised_guard` REFUSES a
    second continuum on this class — placing one tilted the band 4% blue-to-red and cost
    0.0218 dex (RYA-933/1026).
    """
    if 'solar' not in star_id.lower() and 'sun' not in star_id.lower():
        raise ArmNotWired(
            f"{star_id}: the Kitt Peak atlas is the SOLAR flux atlas. Refusing to "
            f"synthesize {star_id} against it (RYA-464's no-silent-substitution rule).")
    return _load_generic_atlas_arm(star_id, region, holding_for_region(region))


PROCYON_UVES_DIAGNOSTICS = ESPRESSO_DIAGNOSTICS + UVES_DIAGNOSTICS


# ── Per-star multi-instrument arm registry (RYA-464 — the unlock) ──────────────
# THE FIX: run_cno/run_phase_a registered a HARDCODED single region ('vis') and the
# espresso/uves arms loaded reflected-solar Vesta UNCONDITIONALLY. So a non-solar star
# would silently synthesize against SUNLIGHT. This replaces that with a PER-STAR arm
# registry: each arm declares its RegionConfig, diagnostics, spectrum loader, and a
# readiness flag. The solar/Vesta path is one case (all arms ready via the existing
# loaders → bit-identical). A non-solar arm whose loader+data+audit are not yet present is
# DECLARED but ready=False → run_phase_a DEFERS it with the reason and NEVER falls back to
# reflected-solar geometry (the loud-fail that closes the latent bug).

class ArmNotWired(RuntimeError):
    """Raised when an arm is requested for a star but its loader/data/audit isn't ready —
    instead of silently substituting reflected-solar (Vesta) data (RYA-464)."""


@dataclass(frozen=True)
class ArmWiring:
    name: str                       # 'harps' | 'uves' | 'uv' | 'ir' | 'espresso'
    region: RegionConfig
    diagnostics: tuple
    loader: str                     # 'harps_normalized' | 'reflected_solar' | 'uves_rya272'
                                    #   | 'hst_stis' | 'ir_crires'
    ready: bool
    defer_reason: str = ''
    provenance: str = ''            # 'measured' | 'cited' | per-arm note


# Reflected-solar (Vesta) arms — the existing RYA-371/372 solar path, untouched.
_SOLAR_ARMS = {
    'harps':    ArmWiring('harps', HARPS_VIS, VIS_DIAGNOSTICS, 'harps_normalized', True,
                          provenance='measured (HARPS Dumusque)'),
    'espresso': ArmWiring('espresso', ESPRESSO_OPT, ESPRESSO_DIAGNOSTICS, 'reflected_solar', True,
                          provenance='measured (Vesta reflected-solar, RYA-370/372)'),
    'uves':     ArmWiring('uves', UVES_OPT, UVES_DIAGNOSTICS, 'reflected_solar', True,
                          provenance='measured (Vesta reflected-solar, RYA-370/372)'),
}

# Procyon — HARPS VIS is runnable today; UVES/UV/IR are DECLARED but gated on their
# loaders + staged+audited data (RYA-272 / HST loader / IR telluric), so ready=False.
_PROCYON_ARMS = {
    'harps': ArmWiring('harps', HARPS_VIS, VIS_DIAGNOSTICS, 'harps_normalized', True,
                       provenance='measured (HARPS ADP, RYA-273)'),
    'uves':  ArmWiring('uves', UVES_OPT, PROCYON_UVES_DIAGNOSTICS, 'uves_rya272', True,
                       provenance='measured (UVES RED760 oi_anchor 2013-10-08, CLEAN telluric, '
                                  'BERV-applied; RYA-272 loader + RYA-271 audit). O I 777 = primary O.'),
    'uv':    ArmWiring('uv', HST_UV, _hst_uv.uv_arm_diagnostics(), 'hst_stis', False,
                       defer_reason='loader BUILT + Procyon STIS staged/conditioned (RYA-471, '
                                    'smoke-proven on E140M: C I 1657 covered, vac->air verified). '
                                    'Synthesis still gated on (1) the Amarsi C/O NLTE grid (RYA-359) '
                                    '— FUV C I carries a large negative correction; amarsi_grid_backend '
                                    'loud-fails by design, and (2) the FUV pseudo-continuum (RYA-426 '
                                    'gate 5, synthesis-not-EW). Flip to ready=True once RYA-359 lands.',
                       provenance='measured (loader built; synthesis gated on RYA-359)'),
    'ir':    ArmWiring('ir', CRIRES_IR, (), 'ir_crires', False,
                       defer_reason='telluric-gated (RYA-373); no 2.3um CO overtone staged '
                                    '(RYA-351: APOGEE weak-CO only); IR conditioning RYA-425.',
                       provenance='measured (pending)'),
}

STAR_ARMS = {'solar': _SOLAR_ARMS, 'procyon': _PROCYON_ARMS}


def star_arm_registry(star_id: str) -> dict:
    """Per-star arm registry (RYA-464). Substring-matched like the other star resolvers.
    Fails loud for an undeclared star — no silent solar-geometry default."""
    s = star_id.strip().lower()
    if s in STAR_ARMS:
        return STAR_ARMS[s]
    for k, v in STAR_ARMS.items():
        if k in s or s in k:
            return v
    raise KeyError(
        f"No multi-instrument arm registry for star_id={star_id!r}. Declare its arms in "
        f"STAR_ARMS (RYA-464) before a multi-arm run — refusing to assume solar geometry.")


def available_arms(star_id: str) -> tuple:
    """(ready_arm_names, {deferred_name: reason}) for `star_id` — the runtime arm map."""
    reg = star_arm_registry(star_id)
    ready = tuple(n for n, a in reg.items() if a.ready)
    deferred = {n: a.defer_reason for n, a in reg.items() if not a.ready}
    return ready, deferred


def resolve_arm_spectrum(star_id: str, arm: ArmWiring):
    """Load the observed spectrum for (star, arm) → (wave_nm, flux). Dispatch by loader
    kind; LOUD-FAIL (ArmNotWired) for an unready arm or a non-solar star routed at the
    reflected-solar loader — the bug RYA-464 closes (no silent Vesta substitution)."""
    if not arm.ready:
        raise ArmNotWired(f"{star_id}/{arm.name}: {arm.defer_reason}")
    if arm.loader == 'harps_normalized':
        return _load_observed_spectrum(star_id)
    if arm.loader == 'reflected_solar':
        if 'solar' not in star_id.lower() and 'sun' not in star_id.lower():
            raise ArmNotWired(
                f"{star_id}/{arm.name}: reflected-solar (Vesta) loader is solar-only; a "
                f"non-solar star must use its own instrument loader (RYA-464). Refusing to "
                f"synthesize {star_id} against reflected sunlight.")
        return _load_reflected_solar_arm(arm.region.instrument.lower())
    if arm.loader == 'hst_stis':                 # RYA-471: HST STIS/COS UV arm
        if 'procyon' not in star_id.lower():
            raise ArmNotWired(
                f"{star_id}/{arm.name}: the HST STIS UV loader is Procyon-only today "
                f"(RYA-222 whitelist). Stage + audit {star_id}'s UV frames before wiring (RYA-464).")
        return _hst_uv.load_procyon_uv_arm('E140M')   # FUV grating covering C I 1657 + O I 1355
    if arm.loader == 'uves_rya272':              # RYA-348 Phase 2: Procyon UVES O I 777 primary-O
        if 'procyon' not in star_id.lower():
            raise ArmNotWired(
                f"{star_id}/{arm.name}: the RYA-272 UVES loader is wired for Procyon's staged "
                f"anchor only. Stage + audit {star_id}'s UVES frames before wiring (RYA-464).")
        return _load_procyon_uves_arm()
    raise ArmNotWired(
        f"{star_id}/{arm.name}: loader {arm.loader!r} not implemented yet "
        f"(IR CRIRES build pending).")


# ── VIS NLTE policy — LTE-by-design, pluggable per arm ────────────────────────
# This is the `lte_by_design` backend. It applies ZERO correction and stamps the
# physics-justified flag per diagnostic. The red/IR/UV arms pass a different
# backend (Amarsi post-hoc grid RYA-359, or in-synthesis departures RYA-361/363)
# with the SAME signature: (diagnostic, A_lte, params) -> (A_nlte, delta, flag, ref).

def vis_lte_backend(diag: 'Diagnostic', a_lte: float, params: dict) -> tuple:
    """LTE-by-design VIS NLTE backend. delta = 0; flag justified per diagnostic.

    Correct treatment, NOT a silent fallback (Ryan handoff 2026-06-19): the VIS
    CNO lines are LTE or near-LTE by their physics. Does NOT loud-fail on an
    absent C/O grid — the grid (RYA-359) is mandatory only for the red/IR/UV arms.
    """
    return float(a_lte), 0.0, diag.nlte_flag, diag.nlte_ref


def amarsi_grid_backend(diag: 'Diagnostic', a_lte: float, params: dict) -> tuple:
    """Placeholder for the red/IR-arm Amarsi C I / O I post-hoc grid (RYA-359).

    Loud-fails until the grid lands — the red-arm O I 777 / FUV C I corrections
    are large and negative; routing them through LTE silently is exactly the
    failure this guards. Wired here so the per-arm NLTE path is real and pluggable;
    NOT used by the VIS deliverable.
    """
    raise NotImplementedError(
        f"Amarsi C/O NLTE grid backend (RYA-359) not yet available for diagnostic "
        f"'{diag.key}'. The red-optical/IR/UV arms require it (O I 777, FUV C I — "
        f"large negative corrections). The VIS arm is LTE-by-design (use "
        f"vis_lte_backend); do not route VIS through here."
    )


NLTE_BACKENDS = {'lte_by_design': vis_lte_backend, 'amarsi_grid': amarsi_grid_backend}


# ── Phase-A cited correction layer (RYA-371) ──────────────────────────────────
# The VIS synthesis MEASURES in 1D-LTE; Phase A applies the CITED 3D/NLTE
# correction on top, per diagnostic, each tagged with its source. VALIDATE-DON'T-
# TUNE: every value here is published / vendored, NEVER fitted to the Asplund
# anchors. Three kinds, by what the literature actually provides for the line:
#   * vendored grid delta — C I 5052/5380: Amarsi 2019 3D-NLTE (pipeline.nlte_cno,
#     RYA-359). A real interpolated Delta = A(3D-NLTE) - A(1D-LTE).
#   * cited 3D anchor     — [O I] 6300: Caffau et al. 2015 (A&A 579 A88, POSP III)
#     full-3D solar A(O)=8.73 with OUR EXACT atomic data (Ni I -2.11 Johansson
#     2003 + [O I] -9.717 Storey & Zeippen 2000; RYA-367). NO hardcoded 3D-1D
#     offset — Caffau 2015 publishes the absolute, not a grid node (RYA-367 rule).
#   * 3D-offset-owed      — CH/CN/C2 molecular bands: no vendored solar 3D grid →
#     reported 1D-LTE, flagged owed (honest, not silently called LTE).

# Cited per-line full-3D solar anchors: {key: (A_3d, unc, flag, source)}. The
# published absolute for the exact line + our atomic data, surfaced as the
# reconciled value — a citation, not a fit.
CITED_3D_ANCHORS = {
    'OI_6300': (
        8.73, 0.05, '3d_lte_caffau2015',
        'Caffau et al. 2015, A&A 579 A88 (POSP III): [O I] 630 nm CO5BOLD full-3D, '
        'Ni I -2.11 (Johansson 2003) + [O I] -9.717 (Storey & Zeippen 2000) = our '
        'atomic data; in gate 8.69+/-0.05. RYA-367 (no hardcoded 3D-1D offset).'),
}
# Atomic C I / O I lines that carry a vendored Amarsi-2019 3D-NLTE grid delta:
# {diagnostic key: (species, representative air wavelength A)}. O I 777 is the
# ESPRESSO PRIMARY O (RYA-455 amendment) — large negative NLTE (~-0.17 solar).
_ATOMIC_GRID_KEYS = {'CI_5052': ('CI', 5052.17), 'CI_5380': ('CI', 5380.34),
                     'OI_777': ('OI', 7773.0)}


def apply_cited_corrections(per_band, params, region) -> list:
    """Attach the cited Phase-A correction to each VIS diagnostic. Returns a list of
    records {key, element, role, a_lte, kind, a_corr, delta, flag, source}. Cited /
    vendored values only — never fitted (RYA-371 validate-don't-tune)."""
    teff = float(params['teff_K']); logg = float(params['logg'])
    feh = float(params['feh']); vmic = float(params['vturb_kms'])
    out = []
    for r in per_band:
        key = r.get('key'); a_lte = r.get('A_X')
        rec = {'key': key, 'element': r.get('element'), 'role': r.get('role'),
               'a_lte': a_lte, 'kind': None, 'a_corr': a_lte, 'delta': 0.0,
               'flag': r.get('nlte_flag'), 'source': r.get('nlte_ref')}
        if not (a_lte is not None and np.isfinite(a_lte)):
            out.append(rec); continue
        if key in CITED_3D_ANCHORS:                       # cited full-3D anchor
            a3d, unc, flag, src = CITED_3D_ANCHORS[key]
            rec.update(kind='cited_3d_anchor', a_corr=a3d, delta=round(a3d - a_lte, 3),
                       flag=flag, source=src, unc=unc)
        elif key in _ATOMIC_GRID_KEYS:                    # vendored Amarsi-2019 grid delta (C I / O I)
            species, wave = _ATOMIC_GRID_KEYS[key]
            try:
                label = nlte_cno.resolve_line(species, wave)
                delta = nlte_cno.cno_nlte_delta(species, label, teff, logg, feh, vmic, a_lte)
                if np.isfinite(delta):
                    nlte_cno.assert_cno_sign(species, label, delta)
                    rec.update(kind='amarsi2019_grid', a_corr=round(a_lte + delta, 3),
                               delta=round(delta, 3), flag='3d_nlte_amarsi2019',
                               source=f'Amarsi, Nissen & Skuladottir 2019 A&A 630 A104, '
                                      f'{species} {label} 3D-NLTE leg ({nlte_cno.select_leg(teff)})')
                else:                                     # outside 4D hull → flag, no silent LTE
                    rec.update(kind='grid_out_of_hull',
                               flag='amarsi2019_out_of_hull', source='Amarsi 2019 grid: query outside 4D hull')
            except Exception as exc:                       # noqa: BLE001 — surface, never fake
                rec.update(kind='grid_error', flag='amarsi2019_error', source=f'grid error: {exc}')
        elif r.get('nlte_flag') == 'lte_molecular_band':  # molecular band, no vendored 3D grid
            rec.update(kind='3d_offset_owed', flag='lte_molecular_band_3d_offset_owed',
                       source='molecular band; no vendored solar 3D-LTE offset grid → '
                              'reported 1D-LTE, 3D offset OWED (Asplund 2005b CH/C2 3D-1D '
                              '0.00..-0.15; not applied — would be uncited)')
        out.append(rec)
    return out


# ── Abundance state + low-level synthesis ─────────────────────────────────────

_ISPEC_SCALE_OFFSET = 12.036   # A(X) = log(N/Ntot)_iSpec + 12.036 (matches abundances_derive)


def _atom_codes(elements, chem_elements, solar_abund) -> dict:
    """{element: iSpec atom code} via create_free_abundances_structure."""
    codes = {}
    for el in elements:
        s = ispec.create_free_abundances_structure([el], chem_elements, solar_abund)
        codes[el] = int(s['code'][0])
    return codes


def _solar_A(elements, chem_elements, solar_abund) -> dict:
    """{element: A(X) on the iSpec internal (Asplund 2009) scale}."""
    out = {}
    for el in elements:
        s = ispec.create_free_abundances_structure([el], chem_elements, solar_abund)
        out[el] = float(s['Abund'][0]) + _ISPEC_SCALE_OFFSET
    return out


def _fixed_ab(state: dict, codes: dict) -> np.recarray:
    """Build a multi-element fixed-abundance recarray on the iSpec SPECTRUM scale.

    state: {element: A(X)} for every element to pin (C, N, O, Ni). All listed
    elements are passed to Turbospectrum as fixed abundances; the abundance under
    fit is just the one we vary between synth calls. babsma uses these to solve the
    molecular equilibrium, so CH/CN/CO band strengths track the set C/N/O.
    """
    elems = list(state.keys())
    fa = np.recarray(len(elems), dtype=[('code', int), ('Abund', float),
                                        ('element', '|U30')])
    for i, el in enumerate(elems):
        fa['code'][i]    = codes[el]
        fa['Abund'][i]   = float(state[el]) - _ISPEC_SCALE_OFFSET
        fa['element'][i] = el
    return fa


def _synth_window(sw_nm, atm, params, ll, iso, sab, fixed_ab,
                  broadening, use_molecules, tmp_dir) -> np.ndarray:
    """Normalized Turbospectrum flux over sw_nm at the given fixed composition.

    Asserts code='turbospectrum' (no SPECTRUM fallback) and threads the
    region-aware instrumental R + per-star vmac/vsini broadening through.
    """
    R, vmac, vsini = broadening
    return ispec.generate_spectrum(
        sw_nm, atm,
        float(params['teff_K']), float(params['logg']), float(params['feh']), 0.0,
        ll, iso, sab, fixed_ab,
        microturbulence_vel=float(params['vturb_kms']),
        macroturbulence=vmac, vsini=vsini, R=R,
        verbose=0, code='turbospectrum',
        use_molecules=use_molecules, tmp_dir=tmp_dir,
    )


# ── Window flux fit (single free element) ─────────────────────────────────────
# Fit A(free_el) by minimizing reduced χ² between observed normalized flux and
# the broadened synthetic flux over the diagnostic's sub-windows. All other CNO
# elements (and pinned blends) are held fixed at `state` — the EOS recomputes
# molecular equilibrium each eval. Single free parameter (the abundance), per the
# RYA-287 convention. σ is the constant 0.01-flux model-adequacy floor (RYA-287):
# its scale shifts χ²ᵣ magnitude but NOT the χ² minimum location (the fitted A).

_SIGMA_FLUX = 0.01
_WSTEP_NM = 0.0002          # fine synthesis grid (0.002 Angstrom)

#: RYA-847 — the curvature sigma and the probe step now live in ONE place,
#: `pipeline.fit_constraint`, and are imported rather than defined here.
#:
#: This module used to own them. RYA-843 then found a THIRD copy of the same
#: accept/reject arithmetic (this one, plus `abundances_derive._fit_synth_flux` and the
#: translation in `measure/synthesis.py`), and copies of a rule are how the three drifted
#: into disagreeing about what an acceptable fit is. The re-export keeps this module's
#: public surface unchanged for its existing callers while there is exactly one
#: definition (RYA-845: declare it once).
from pipeline.fit_constraint import (                              # noqa: E402
    CURVATURE_PROBE_STEP_DEX, measure_constraint,
    # RYA-1214 — RE-EXPORTED, NOT CALLED HERE. `_fit_element` calls `measure_constraint`,
    # which computes this internally; importing it too would NOT create a second sigma.
    # It stays because `tests/test_curvature_sigma_rya848.py` imports it from THIS module
    # (it was a re-export before the call moved), and dropping it broke CI on a test about
    # a function whose behaviour never changed. A removed re-export is an API change.
    curvature_sigma)  # noqa: F401
# RYA-1214 — the SAME decider the band-product route and the Engine-B handler call.
# `STALE.md` names "RYA-847's synthesis constraint gate" as what retires this path's
# sigma clip; importing it rather than re-implementing the check is what makes the three
# routes unable to disagree (RYA-845/847).
from pipeline.constraint_gate import verdict as constraint_verdict   # noqa: E402


def _fit_element(obs_w_nm, obs_f, atm, params, free_el, state, codes,
                 windows_A, use_molecules, broadening, a_lo, a_hi,
                 ll, iso, sab, tmp_dir) -> dict:
    """Minimize χ²ᵣ over A(free_el) across the diagnostic windows.

    Returns {A_X, red_chi2, sigma_fit, n_pix, n_eval, status}. No silent fallback:
    a synthesis error → status='failed', A_X=nan (never substituted).
    """
    # Pre-slice observed pixels per window (rest-frame), build synthesis grids.
    segs = []
    for (lo_A, hi_A) in windows_A:
        lo_nm, hi_nm = lo_A / 10.0, hi_A / 10.0
        m = (obs_w_nm >= lo_nm) & (obs_w_nm <= hi_nm)
        if m.sum() < 5:
            continue
        sw = np.arange(lo_nm, hi_nm + _WSTEP_NM * 0.5, _WSTEP_NM)
        segs.append((obs_w_nm[m], obs_f[m], sw))
    if not segs:
        return {'A_X': np.nan, 'red_chi2': np.nan, 'sigma_fit': np.nan,
                'n_pix': 0, 'n_eval': 0, 'status': 'failed',
                'reason': 'no observed pixels in windows'}

    n_pix = int(sum(len(ow) for ow, _, _ in segs))
    fail = [None]
    n_eval = [0]
    local = dict(state)

    def chi2(a_x):
        n_eval[0] += 1
        local[free_el] = float(a_x)
        fa = _fixed_ab(local, codes)
        tot = 0.0
        for ow, of, sw in segs:
            try:
                sf = _synth_window(sw, atm, params, ll, iso, sab, fa,
                                   broadening, use_molecules, tmp_dir)
            except Exception as exc:          # noqa: BLE001 — surface, never fake
                fail[0] = str(exc)
                return 1e30
            sf_i = interp1d(sw, sf, bounds_error=False, fill_value=1.0)(ow)
            r = (of - sf_i) / _SIGMA_FLUX
            tot += float(np.nansum(r * r))
        return tot

    res = minimize_scalar(chi2, bounds=(a_lo, a_hi), method='bounded',
                          options={'xatol': 1e-3})
    if fail[0] is not None:
        return {'A_X': np.nan, 'red_chi2': np.nan, 'sigma_fit': np.nan,
                'n_pix': n_pix, 'n_eval': int(n_eval[0]), 'status': 'failed',
                'reason': f'synthesis error: {fail[0]}'}

    a_best = float(res.x)
    dof = max(n_pix - 1, 1)
    chi2_min = chi2(a_best)
    red_chi2 = chi2_min / dof
    # σ from the χ² curvature. The rationale lives in `curvature_sigma`'s docstring and
    # is deliberately NOT restated here: this value is the published σ_stat for N and O,
    # and a second copy of the reasoning beside the call is how a number drifts from its
    # own justification (RYA-845).
    #
    # 🔴 RYA-1214 — THE np.clip IS GONE, AND ITS OWN STATED CONDITION IS WHAT RETIRES IT.
    #
    # The clip read `np.clip(sigma_fit, 0.0, 1.0)` and was left in place deliberately,
    # annotated "RYA-847 owns removing this" because that ticket "needs its sweep's
    # constraint metric to define unconstrained; fixing it here would back-fit that
    # criterion from CNO". RYA-847 HAS LANDED. Its sweep ran over 9 cells and 581
    # synthesis lines, and its answer was that NO transferable threshold exists —
    # `constraint_gate.SYNTH_CONSTRAINT` is None PERMANENTLY, not pending — with the
    # NON-MINIMUM check standing as the one criterion that survived, because zero is the
    # boundary between "chi2 rose away from the answer" and "it did not" and so cannot be
    # tuned. So the thing the clip was waiting for is decided, and it is decided in a form
    # that does not need the clip at all.
    #
    # WHAT THE CLIP WAS DOING, MEASURED: `solar_vis_cno_product.csv` publishes
    # sigma_stat = 1.000 for nitrogen (CN_red) and sigma_fit = 1.000 for CI_5380 — a
    # SENTINEL WEARING THE SHAPE OF A MEASURED 1 DEX. `data/audit/cno_synthesis/STALE.md`
    # names exactly this as the last known defect in this path and names the constraint
    # gate as what retires it. An unconstrained fit now says so through the gate instead
    # of through a number that looks like an uncertainty.
    #
    # ⚠️ `measure_constraint` REPLACES the bare `curvature_sigma` call rather than being
    # added beside it. Both compute σ from the same curvature, and calling each would
    # make two σ for one fit — the RYA-845 shape. It also returns the frac_rise the gate
    # needs, which the bare call never produced, which is why this path had no way to ask
    # the question before.
    m = measure_constraint(chi2, a_best=a_best, chi2_min=chi2_min, red_chi2=red_chi2,
                           a_lo=a_lo, a_hi=a_hi)
    cv = constraint_verdict(m)
    edge = m.edge_distance_dex < 1e-2
    status = 'edge_pinned' if edge else ('unconstrained' if not cv.ok else 'ok')
    return {'A_X': round(a_best, 3), 'red_chi2': round(float(red_chi2), 3),
            'sigma_fit': (round(float(m.sigma_A), 3) if np.isfinite(m.sigma_A)
                          else float('nan')),
            'n_pix': n_pix, 'n_eval': int(n_eval[0]),
            'status': status,
            'frac_rise_weaker': float(m.frac_rise_weaker),
            'edge_distance_dex': float(m.edge_distance_dex),
            'constrained': bool(cv.ok),
            'reason': '' if cv.ok else cv.reason}


# ── Preflight assertions (no silent fallback) ─────────────────────────────────

def _molecules_cover(windows_A) -> bool:
    """True iff at least one .bsyn molecular file spans each requested window.

    Filenames encode nm ranges, e.g. 12CH_400-450.bsyn. iSpec globs them at synth
    time; we fail loud here if a molecular band is requested but uncovered.
    """
    import re
    files = list(_MOLECULES_DIR.glob('*.bsyn'))
    spans = []
    for f in files:
        m = re.match(r'.*_(\d+)-(\d+)\.bsyn', f.name)
        if m:
            spans.append((float(m.group(1)), float(m.group(2))))
    for (lo_A, hi_A) in windows_A:
        lo_nm, hi_nm = lo_A / 10.0, hi_A / 10.0
        if not any(s_lo <= lo_nm and hi_nm <= s_hi for s_lo, s_hi in spans):
            return False
    return True


def preflight(region: RegionConfig, star_id: str, diagnostics) -> dict:
    """Assert the CRITICAL no-silent-fallback invariants; return broadening tuple.

    - Turbospectrum only (RT code is hardcoded to 'turbospectrum' in _synth_window).
    - gf via gf_resolver (canonical table present; the GES linelist is rescaled in
      _load_synth_resources — RYA-353).
    - per-star broadening (RYA-288, fail-loud — never a solar default for Procyon).
    - molecular lists cover every molecular band (RYA-236, at the iSpec path).
    - telluric clearance for arms that require it (IR; VIS does not).
    """
    print(f"  [preflight] region={region.name} ({region.instrument}, R={region.R:.0f})")

    # gf_resolver / canonical gf table
    canon = PATHS.get('canonical_gf') if hasattr(PATHS, 'get') else None
    canon = canon or (Path(PATHS['linelist_solar']).parent / 'canonical_gf.csv')
    if not Path(str(canon)).exists():
        raise FileNotFoundError(
            f"Canonical gf table not found ({canon}); the synth path must resolve gf "
            f"via gf_resolver (RYA-353). Refusing to run with unresolved/duplicated gf.")
    print(f"  [preflight] gf_resolver canonical table present: {Path(str(canon)).name}")

    # per-star broadening (fail-loud)
    R_star, vmac, vsini, fit_vmac = _resolve_broadening(star_id)
    if fit_vmac:
        raise NotImplementedError(
            f"[{star_id}] vmac='fit' (RYA-309 §3.3) — RT vmac fit not wired in CNO "
            f"engine; refusing to use the init guess as a fixed vmac (RYA-288).")
    broadening = (float(region.R), float(vmac), float(vsini))   # region LSF + per-star
    print(f"  [preflight] broadening (per-star, RYA-288): R={region.R:.0f}, "
          f"vmac={vmac} km/s, vsini={vsini} km/s")

    # molecular coverage for molecular bands
    for d in diagnostics:
        if d.use_molecules and not _molecules_cover(d.windows_A):
            raise FileNotFoundError(
                f"No molecular .bsyn list covers {d.key} windows {d.windows_A} at "
                f"{_MOLECULES_DIR}. RYA-236 lists must be present (verify at iSpec "
                f"tool path; secured by RYA-360). Refusing to synthesize the band "
                f"without its molecular opacity.")
    print(f"  [preflight] molecular lists cover all molecular bands "
          f"({_MOLECULES_DIR.name}/*.bsyn)")

    # ── telluric clearance gate ───────────────────────────────────────────────
    # 🔴 RYA-1214 — IT ASKS THE HOLDING NOW. This was an UNCONDITIONAL raise: any region
    # with `telluric_correction_required=True` could not run, whatever spectrum it was
    # pointed at, because when it was written no such region existed and the flag meant
    # "not wired". A "clearance flag" is also the wrong mechanism — a flag is an assertion
    # by the caller, and RYA-1026 is explicit that this fact is read "through
    # telluric_policy.applied_state, never inferred" and that VERIFIED state outranks a
    # declaration (RYA-1194/1196).
    #
    # So the gate resolves the holding the loader will actually read and asks TWO questions:
    # is that holding a corrected science basis, and do any FIT WINDOWS sit inside an
    # enumerated telluric band. The second matters even on a corrected holding, because
    # telluric reachability is a per-OBSERVATION property and a correction is not a
    # guarantee (RYA-1193).
    if region.telluric_correction_required:
        from pipeline import telluric_display_policy as _tdp
        from pipeline import telluric_policy as _tp
        _hold = holding_for_region(region)
        _state = _tdp.display_state(_hold)
        if _state not in ('CLEAN', 'CLEAN_WITH_ANOMALY'):
            raise RuntimeError(
                f"Region {region.name} requires telluric-corrected input and its holding "
                f"{_hold!r} is display_state={_state!r}. Refusing to fit CNO over a "
                f"spectrum that is not a corrected science basis (RYA-1026). This is the "
                f"holding's VERIFIED state, not a flag the caller passed.")
        _inside = [(lo, hi, d.key) for d in diagnostics for (lo, hi) in d.windows_A
                   if _tp.in_telluric_band(0.5 * (lo + hi))]
        if _inside:
            raise RuntimeError(
                f"Region {region.name}: {len(_inside)} fit window(s) fall inside an "
                f"enumerated telluric band — {_inside[:4]}. The holding is corrected, but "
                f"telluric reachability is per-OBSERVATION and a correction is not a "
                f"guarantee (RYA-1193). Move the window or state the exception; refusing "
                f"to fit there by default.")
        print(f"  [preflight] telluric gate: {_hold} is {_state}, and 0 of "
              f"{sum(len(d.windows_A) for d in diagnostics)} fit window(s) fall inside an "
              f"enumerated telluric band (telluric_policy.TELLURIC_BANDS)")
    else:
        print(f"  [preflight] telluric gate: not required for {region.name} (optical)")

    # NLTE backend resolves
    if region.nlte_backend not in NLTE_BACKENDS:
        raise KeyError(f"Unknown NLTE backend '{region.nlte_backend}'")
    print(f"  [preflight] NLTE backend = {region.nlte_backend}")
    return broadening


# ── Engine ────────────────────────────────────────────────────────────────────

@dataclass
class CNOResult:
    star_id: str
    region: str
    abundances: dict = field(default_factory=dict)       # element -> A(X)
    per_band: list = field(default_factory=list)         # list of band dicts
    iterations: int = 0
    converged: bool = False
    flags: list = field(default_factory=list)
    provenance: dict = field(default_factory=dict)
    uncertainty: dict = field(default_factory=dict)      # element -> {stat, sys, tot}
    phase_a_corrections: list = field(default_factory=list)  # RYA-371 cited 3D/NLTE per diagnostic
    #: RYA-1214 — which elements this REGION actually measured, and which it only pinned.
    #: A region need not carry a primary for all three (the near-UV has none for carbon),
    #: and `abundances` holds a number either way — so without these two a pinned seed
    #: reads exactly like a measurement (RYA-833).
    measured_elements: tuple = ()
    pinned_elements: tuple = ()


def _seed_abundances(star_id, params, codes, solar_A_ispec, feh) -> dict:
    """Initial A(C/N/O/Ni). C from the solar anchor scaled by [Fe/H]; N, O likewise;
    Ni pinned (canonical/EW). The CH/CN/[O I] fits refine C/N/O from here."""
    seed = {}
    for el in ('C', 'N', 'O', 'Ni'):
        seed[el] = float(SOLAR_ASPLUND2021[el]) + (feh if star_id != 'solar' else 0.0)
    return seed


def run_cno(star_id: str, region_name: str = 'vis', *,
            params_override: dict = None, max_iter: int = 5,
            with_systematics: bool = True, out_dir: Path = None,
            tmp_dir: str = '/tmp/ispec_cno', pins: dict = None) -> CNOResult:
    """Region-aware C/N/O synthesis for `star_id` over `region_name`.

    Flow: seed → fit CH (A(C)) → fit CN given A(C) → fit [O I]+Ni given A(C) →
    re-synthesize (EOS recomputes equilibrium) and re-fit until A(C)/A(N)/A(O)
    converge (Δ < 0.01 dex, max `max_iter`); flag co_equilibrium_not_converged
    otherwise. C I 5052/5380 (LTE) and C2 Swan are fit once as cross-checks.
    """
    region = REGIONS[region_name]
    # 🔴 RYA-1214 — DERIVED FROM THE REGION, NOT PINNED TO VIS. This read
    # `diagnostics = VIS_DIAGNOSTICS` unconditionally, so `--region` chose the RegionConfig
    # (instrument, R, band edges, NLTE backend) while the DIAGNOSTICS stayed HARPS-optical.
    # A second region was therefore not merely unwired: adding one would have fitted VIS
    # windows against a near-UV spectrum and reported it under the new region's name.
    diagnostics = REGION_DIAGNOSTICS[region_name]
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)     # RYA-344: TS tmp_dir must exist
    out_dir = Path(out_dir) if out_dir else (Path(PATHS['solar_ew']).parent.parent /
                                             'audit' / 'cno_synthesis')
    out_dir.mkdir(parents=True, exist_ok=True)

    rec = get_star_params(star_id)
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    if params_override:
        params.update(params_override)
    feh = params['feh']

    print(f"\n{'='*72}\n  C/N/O synthesis — {star_id} / {region.name} "
          f"(Teff={params['teff_K']:.0f} logg={params['logg']:.2f} "
          f"[Fe/H]={feh:+.2f} xi={params['vturb_kms']:.2f})\n{'='*72}")

    broadening = preflight(region, star_id, diagnostics)
    nlte_backend = NLTE_BACKENDS[region.nlte_backend]

    atm = _load_atmosphere(params['teff_K'], params['logg'], feh, params['vturb_kms'])
    ll, iso, chem = _load_synth_resources()
    sab = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    # Same reason as the diagnostics above: `_load_observed_spectrum` is the HARPS loader,
    # so a near-UV region would have been fitted against an optical spectrum that does not
    # even cover its windows. Dispatched on the region's own instrument.
    obs_w, obs_f = _load_region_spectrum(star_id, region)

    codes = _atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)
    solar_A_ispec = _solar_A(('C', 'N', 'O', 'Ni'), chem, sab)
    state = _seed_abundances(star_id, params, codes, solar_A_ispec, feh)
    # RYA-1214 — a pin REPLACES the solar seed for an element this region does not measure.
    # Applied before the first synthesis, because the seed is what feeds babsma's molecular
    # equilibrium: an override applied after the atmosphere is built changes nothing, which
    # is the inert-override failure `rya1120_xi_campaign` documents for xi.
    for _el, _v in (pins or {}).items():
        state[_el] = float(_v)
    print(f"  seed A: " + "  ".join(f"{e}={state[e]:.2f}" for e in ('C', 'N', 'O', 'Ni')))
    if pins:
        print(f"  PINNED (input, not measured): "
              + "  ".join(f"A({e})={v:.3f}" for e, v in sorted(pins.items())))

    result = CNOResult(star_id=star_id, region=region.name)
    by_key = {d.key: d for d in diagnostics}

    def _fit(diag, a_center):
        a_lo, a_hi = a_center - 1.2, a_center + 1.2
        t0 = time.time()
        r = _fit_element(obs_w, obs_f, atm, params, diag.element, state, codes,
                         diag.windows_A, diag.use_molecules, broadening,
                         a_lo, a_hi, ll, iso, sab, tmp_dir)
        r['key'] = diag.key
        r['element'] = diag.element
        r['role'] = diag.role
        r['wall_s'] = round(time.time() - t0, 1)
        a_nlte, delta, flag, ref = nlte_backend(diag, r['A_X'], params)
        r.update({'A_X_nlte': a_nlte, 'nlte_delta': delta,
                  'nlte_flag': flag, 'nlte_ref': ref})
        return r

    # ── CNO equilibrium iteration, in the region's own dependency order ───────
    # 🔴 DERIVED FROM `role`, NOT NAMED. This was
    #     primary = {'C': by_key['CH_Gband'], 'N': by_key['CN_red'], 'O': by_key['OI_6300']}
    # which is three VIS keys, so any other region raised KeyError before it could run —
    # and a region that happened to reuse a key would have silently borrowed VIS's
    # diagnostic. `primary_by_element` reads the declared roles and refuses an ambiguous
    # or absent one rather than picking.
    primary = primary_by_element(diagnostics)
    # An element with no primary IN THIS REGION is not measured here. Its seed value is
    # still needed (it is an input to the molecular equilibrium of the ones that ARE
    # measured), so it is PINNED and recorded as pinned — never reported as a measurement.
    measured_elements = tuple(e for e in ('C', 'N', 'O') if e in primary)
    pinned_elements = tuple(e for e in ('C', 'N', 'O') if e not in primary)
    if pinned_elements:
        result.flags.append(
            'pinned_not_measured:' + ','.join(f'{e}={state[e]:.3f}' for e in pinned_elements))
        print(f"  region {region.name!r} has no primary for "
              f"{', '.join(pinned_elements)} — PINNED at "
              + ", ".join(f"A({e})={state[e]:.3f}" for e in pinned_elements)
              + " and NOT reported as measured")
    converged = False
    last = {}
    for it in range(1, max_iter + 1):
        print(f"\n  ── iteration {it} ──")
        deltas = {}
        for el in measured_elements:
            d = primary[el]
            center = state[el]
            r = _fit(d, center)
            last[d.key] = r
            # 🔴 RYA-1214 — A FIT THE GATE REJECTED MUST NOT PIN THE EQUILIBRIUM.
            # `np.isfinite(A_X)` was the only test, and `minimize_scalar` returns a finite
            # number for a flat objective just as readily as for a real minimum — which is
            # what the retired sigma clip was papering over. In THIS loop the consequence
            # compounds rather than staying local: `state[el]` is pinned into the next
            # element's synthesis (CH -> CN|C -> [O I]|C), so an unconstrained carbon fit
            # would set the molecular equilibrium that nitrogen and oxygen are then
            # measured against. An undetermined abundance must not become a fixed input.
            if np.isfinite(r['A_X']) and r.get('constrained', True):
                deltas[el] = abs(r['A_X'] - state[el])
                state[el] = r['A_X']     # pin for the next element (EOS coupling)
            else:
                deltas[el] = np.nan
                if np.isfinite(r['A_X']) and not r.get('constrained', True):
                    result.flags.append(f'{el}_primary_unconstrained:{d.key}')
                    print(f"    {d.key:10s} REFUSED as the {el} pin — {r['reason'][:88]}")
            print(f"    {d.key:10s} A({el})={r['A_X']}  χ²ᵣ={r['red_chi2']}  "
                  f"σfit={r['sigma_fit']}  [{r['status']}]  ({r['wall_s']}s)")
        result.iterations = it
        finite = [v for v in deltas.values() if np.isfinite(v)]
        if finite and max(finite) < 0.01:
            converged = True
            print(f"  converged after {it} iter (max Δ={max(finite):.4f} dex)")
            break
        print(f"  Δ this iter: " + "  ".join(f"{e}={deltas[e]:.3f}" for e in deltas))
    result.converged = converged
    if not converged:
        result.flags.append('co_equilibrium_not_converged')

    # ── Cross-checks (LTE), derived from `role` rather than named ─────────────
    _cross = tuple(d.key for d in diagnostics if d.role == 'cross_check')
    if _cross:
        print(f"\n  ── cross-checks (LTE): {', '.join(_cross)} ──")
    for key in _cross:
        d = by_key[key]
        r = _fit(d, state['C'])
        last[d.key] = r
        print(f"    {d.key:10s} A(C)={r['A_X']}  χ²ᵣ={r['red_chi2']}  "
              f"[{r['status']}]  flag={r['nlte_flag']}  ({r['wall_s']}s)")

    result.per_band = [last[d.key] for d in diagnostics if d.key in last]
    result.abundances = {'C': state['C'], 'N': state['N'], 'O': state['O']}
    result.measured_elements = measured_elements
    result.pinned_elements = pinned_elements

    # ── Phase-A cited correction layer (RYA-371): 1D-LTE → cited 3D/NLTE ───────
    corrections = apply_cited_corrections(result.per_band, params, region)
    result.phase_a_corrections = corrections
    print(f"\n  ── Phase-A cited corrections ({region.instrument} arm) — "
          f"validate-don't-tune (cited/vendored only) ──")
    print(f"    {'diagnostic':10s} {'el':2s} {'1D-LTE':>7s} {'corr':>7s} "
          f"{'recon':>7s}  kind / source")
    for c in corrections:
        al = f"{c['a_lte']:.3f}" if isinstance(c['a_lte'], float) and np.isfinite(c['a_lte']) else '  N/A '
        ac = f"{c['a_corr']:.3f}" if isinstance(c['a_corr'], float) and np.isfinite(c['a_corr']) else '  N/A '
        dl = f"{c['delta']:+.3f}" if np.isfinite(c.get('delta', np.nan)) else '   -- '
        print(f"    {c['key']:10s} {c['element']:2s} {al:>7s} {dl:>7s} {ac:>7s}  "
              f"[{c['kind']}] {(c['source'] or '')[:64]}")

    # ── Uncertainty budget (Type A statistical + Type B systematic) ───────────
    result.uncertainty = _uncertainty_budget(
        result, last, by_key, star_id, params, rec, state, codes, obs_w, obs_f,
        atm, broadening, ll, iso, sab, tmp_dir, with_systematics,
        diagnostics=diagnostics)

    # ── C/O ratio ─────────────────────────────────────────────────────────────
    aC, aO = state['C'], state['O']
    co = float(10 ** (aC - aO)) if np.isfinite(aC) and np.isfinite(aO) else np.nan
    result.abundances['C/O'] = round(co, 3)

    # ── Provenance ─────────────────────────────────────────────────────────────
    result.provenance = {
        'engine': 'pipeline.cno_synthesis (RYA-237)',
        'rt_code': 'turbospectrum',
        'region': region.name, 'instrument': region.instrument, 'R_LSF': region.R,
        'atomic_linelist': 'GESv6_atom_hfs_iso (canonical-gf via gf_resolver, RYA-353)',
        'molecular_lists': f'{_MOLECULES_DIR.name}/*.bsyn (RYA-236: CH/CN/C2/CO/OH/NH)',
        'broadening': {'R': broadening[0], 'vmac': broadening[1], 'vsini': broadening[2],
                       'source': rec.get('source', ''), 'rule': 'per-star RYA-288'},
        'nlte': {'backend': region.nlte_backend,
                 'policy': 'VIS synthesis is 1D-LTE; Phase-A cited corrections applied '
                           'on top (RYA-371): C I 5052/5380 Amarsi-2019 3D-NLTE grid; '
                           '[O I] 6300 cited Caffau-2015 full-3D anchor 8.73 (RYA-367, '
                           'no hardcoded offset); CH/CN/C2 molecular 3D offset OWED'},
        'phase_a_corrections': result.phase_a_corrections,
        'solar_reference': 'Asplund 2021 (A&A 653, A141) via SOLAR_ASPLUND2021',
        'params': params,
        'caveats': ['Ni I 6300.34 gf in the [O I] blend resolves via gf_resolver '
                    '= Johansson+2003 log gf -2.11 (RYA-365). Broader store-#2 '
                    'per-star-master reroute remains the RYA-353 follow-on umbrella.'],
    }

    _write_product(result, out_dir)
    return result


# ── Uncertainty budget ────────────────────────────────────────────────────────

def _uncertainty_budget(result, last, by_key, star_id, params, rec, state, codes,
                        obs_w, obs_f, atm, broadening, ll, iso, sab, tmp_dir,
                        with_systematics, diagnostics=VIS_DIAGNOSTICS) -> dict:
    """Type A (statistical) + Type B (systematic stellar-param sensitivity).

    Statistical:
      C — std across the carbon diagnostics (CH primary + C I 5052/5380 + C2) / √N
          (multi-diagnostic scatter is the honest random term).
      N, O — single primary band → the χ² fit 1σ (sigma_fit).
    Systematic (opt-in, primary band per element): refit the primary band at
      Teff+e_teff, logg+e_logg, [Fe/H]+e_feh, ξ±0.1; quadrature of the ΔA. A
      perturbation whose parameter error is negligible (e.g. solar e_teff=1,
      e_logg≈0) is skipped — keeps the solar validate cheap, real for Procyon.
    """
    budget = {}
    # Statistical
    # 🔴 RYA-1214 — DERIVED PER ELEMENT, NOT THREE HARDCODED VIS KEYS. This read
    # `('CH_Gband','CI_5052','CI_5380','C2_Swan')` for carbon and `last['CN_red']` /
    # `last['OI_6300']` for N and O, so in any other region every sigma_stat would have
    # come back NaN from a `.get` that could never hit — an UNMEASURED bar on a band that
    # measured perfectly well. The rule itself is unchanged and stated once: where an
    # element has TWO OR MORE diagnostics that measured something, its random term is
    # their scatter / sqrt(N); where it has one, it is that fit's own curvature.
    #
    # And the membership rule is the gate's verdict, in both branches: this pooled every
    # finite `A_X`, so an unconstrained band widened (or narrowed) the scatter of the one
    # element whose sigma_stat is a scatter rather than a curvature.
    prim = primary_by_element(diagnostics)
    stat = {}
    for el in ('C', 'N', 'O'):
        vals = [last[d.key]['A_X'] for d in diagnostics
                if d.element == el and d.key in last
                and np.isfinite(last[d.key]['A_X'])
                and last[d.key].get('constrained', True)]
        if len(vals) >= 2:
            stat[el] = float(np.std(vals, ddof=1) / np.sqrt(len(vals)))
        elif el in prim and prim[el].key in last:
            stat[el] = float(last[prim[el].key].get('sigma_fit', np.nan))
        else:
            stat[el] = float('nan')

    # Systematic
    sys_budget = {e: 0.0 for e in ('C', 'N', 'O')}
    if with_systematics:
        primary = prim                      # RYA-1214: the region's own, same derivation
        e_teff = float(rec.get('e_teff', 0.0))
        e_logg = float(rec.get('e_logg', 0.0))
        e_feh = float(rec.get('e_feh', 0.0))
        perturbs = []
        if e_teff >= 5.0:
            perturbs.append(('teff_K', e_teff))
        if e_logg >= 0.02:
            perturbs.append(('logg', e_logg))
        if e_feh >= 0.02:
            perturbs.append(('feh', e_feh))
        perturbs.append(('vturb_kms', 0.1))   # always include a ξ sensitivity term
        print(f"\n  ── systematics: {len(perturbs)} perturbation(s) × "
              f"{len(primary)} primary band(s) ──")
        for el in sorted(primary):
            d = primary[el]
            base = state[el]
            terms = []
            for pkey, pstep in perturbs:
                pp = dict(params); pp[pkey] = pp[pkey] + pstep
                try:
                    atm_p = _load_atmosphere(pp['teff_K'], pp['logg'], pp['feh'],
                                             pp['vturb_kms'])
                except Exception:
                    continue
                r = _fit_element(obs_w, obs_f, atm_p, pp, el, state, codes,
                                 d.windows_A, d.use_molecules, broadening,
                                 base - 1.0, base + 1.0, ll, iso, sab, tmp_dir)
                if np.isfinite(r['A_X']):
                    terms.append(abs(r['A_X'] - base))
            sys_budget[el] = float(np.sqrt(np.sum(np.square(terms)))) if terms else 0.0
            print(f"    {el}: σ_sys={sys_budget[el]:.3f} dex ({len(terms)} terms)")

    # 🔴 RYA-1214 — AN UNMEASURABLE sigma_stat IS NOT ZERO, AND THIS IS WHERE THE
    # RETIRED CLIP'S SENTINEL WOULD HAVE COME BACK WEARING THE OPPOSITE SIGN.
    #
    # This read `s_a = stat[el] if np.isfinite(stat[el]) else 0.0`. With the
    # `np.clip(..., 0.0, 1.0)` in place an unmeasurable curvature arrived here as 1.000
    # and was reported as a 1 dex bar; with the clip removed it arrives as NaN and this
    # line would have reported 0.000 — the TIGHTEST possible bar for the LEAST constrained
    # possible fit, which is verbatim the defect `fit_constraint.ConstraintMetrics`
    # documents ("the old code returned 0.000 for a railed fit ... and clipped a flat
    # objective to a plausible-looking 1.000"). Removing the clip alone would have made
    # this path strictly worse.
    #
    # So an unmeasured term stays UNMEASURED and says so, and sigma_tot is withheld rather
    # than computed from a stand-in: a total that silently omits its statistical part is
    # not a smaller uncertainty, it is a different quantity (RYA-907).
    for el in ('C', 'N', 'O'):
        s_b = sys_budget[el]
        if not np.isfinite(stat[el]):
            budget[el] = {
                'stat': None, 'sys': round(s_b, 3), 'tot': None,
                'stat_state': 'UNMEASURED',
                'stat_note': ('the chi2 curvature was not measurable for this element\'s '
                              'primary band (railed fit, or an objective with no curvature '
                              'to invert), so there is no statistical term to report and '
                              'sigma_tot is withheld. NOT zero and NOT the 1.000 the '
                              'retired np.clip used to emit (RYA-848/RYA-1214).'),
            }
            continue
        s_a = float(stat[el])
        budget[el] = {'stat': round(s_a, 3), 'sys': round(s_b, 3),
                      'stat_state': 'MEASURED',
                      'tot': round(float(np.sqrt(s_a ** 2 + s_b ** 2)), 3)}
    return budget


# ── Output ────────────────────────────────────────────────────────────────────

def _write_product(result: CNOResult, out_dir: Path) -> None:
    rows = []
    for el in ('C', 'N', 'O'):
        unc = result.uncertainty.get(el, {})
        # 🔴 RYA-1214 — A PINNED SEED IS NOT A MEASUREMENT, AND `abundances` CANNOT TELL
        # THEM APART. The near-UV region has no carbon primary, so A(C) in that run is
        # whatever was pinned into the molecular equilibrium — a real and necessary INPUT,
        # and a number this CSV would otherwise publish in the same column as the two
        # abundances the region actually fitted (RYA-833).
        measured = (el in (result.measured_elements or ('C', 'N', 'O')))
        rows.append({
            'element': el, 'A_X': result.abundances.get(el),
            'value_state': 'MEASURED' if measured else 'PINNED_INPUT_NOT_MEASURED',
            'sigma_stat': unc.get('stat'), 'sigma_sys': unc.get('sys'),
            'sigma_tot': unc.get('tot'),
            # RYA-1214: a blank sigma_stat must say WHY it is blank. "Not measured" and
            # "not written out" look identical in a CSV otherwise (RYA-833).
            'sigma_stat_state': unc.get('stat_state', ''),
        })
    rows.append({'element': 'C/O', 'A_X': result.abundances.get('C/O'),
                 'sigma_stat': None, 'sigma_sys': None, 'sigma_tot': None,
                 'sigma_stat_state': 'NOT_PROPAGATED'})
    prod = pd.DataFrame(rows)
    band = pd.DataFrame(result.per_band)
    base = out_dir / f'{result.star_id}_{result.region}_cno'
    prod.to_csv(f'{base}_product.csv', index=False)
    band.to_csv(f'{base}_per_band.csv', index=False)
    with open(f'{base}_provenance.json', 'w') as fh:
        json.dump({'abundances': result.abundances, 'iterations': result.iterations,
                   'converged': result.converged, 'flags': result.flags,
                   'uncertainty': result.uncertainty,
                   'provenance': result.provenance}, fh, indent=2, default=str)
    print(f"\n  [out] {base}_product.csv / _per_band.csv / _provenance.json")


# ── Validation (solar VIS gates) ──────────────────────────────────────────────

# Acceptance anchors — Asplund 2021 solar VIS reference (RYA-162 target table).
SOLAR_VIS_GATES = {
    'C': (8.46, 0.05), 'N': (7.83, 0.07), 'O': (8.69, 0.05), 'C/O': (0.59, 0.08),
}


def validate_solar(result: CNOResult) -> bool:
    """Print the solar-VIS gate table; CH vs C I agreement. Returns all-pass bool."""
    print(f"\n{'='*72}\n  SOLAR-VIS ACCEPTANCE GATES (Asplund 2021)\n{'='*72}")
    print(f"  {'qty':6s} {'derived':>10s} {'σ_tot':>7s}  {'target':>8s} {'±':>5s}  result")
    all_pass = True
    for q, (tgt, terr) in SOLAR_VIS_GATES.items():
        val = result.abundances.get(q)
        if val is None or not np.isfinite(val):
            print(f"  {q:6s} {'nan':>10s}  — gate INDETERMINATE"); all_pass = False; continue
        sig = result.uncertainty.get(q, {}).get('tot', 0.0) if q != 'C/O' else 0.0
        within = abs(val - tgt) <= (terr + (sig or 0.0) + 1e-9)
        all_pass = all_pass and within
        print(f"  {q:6s} {val:>10.3f} {sig or 0.0:>7.3f}  {tgt:>8.2f} {terr:>5.2f}  "
              f"{'PASS' if within else 'FAIL'} (Δ={val-tgt:+.3f})")

    # CH vs C I agreement
    bands = {b['key']: b for b in result.per_band}
    ch = bands.get('CH_Gband', {}).get('A_X')
    ci = [bands[k]['A_X'] for k in ('CI_5052', 'CI_5380')
          if k in bands and np.isfinite(bands[k].get('A_X', np.nan))]
    if ch is not None and np.isfinite(ch) and ci:
        ci_mean = float(np.mean(ci))
        agree = abs(ch - ci_mean) <= 0.15
        print(f"\n  CH G-band A(C)={ch:.3f}  vs  C I 5052/5380 A(C)={ci_mean:.3f}  "
              f"Δ={ch-ci_mean:+.3f}  {'AGREE' if agree else 'DISAGREE'} (±0.15)")
    print(f"\n  CO-equilibrium iterations: {result.iterations}  "
          f"(converged={result.converged}); flags={result.flags or 'none'}")
    print(f"  Overall: {'ALL GATES PASS' if all_pass else 'GATES NOT ALL MET'}")
    return all_pass


# ── [O I] 6300 blend partition (RYA-365) ─────────────────────────────────────
# The Ni I 6300.34 gf feeding the [O I] 6300 joint fit must resolve through
# gf_resolver (canonical = Johansson+2003 -2.11). This diagnostic "sees the
# blend": it decomposes total / Ni / [O I] contribution + inferred A(O), BEFORE
# (the stale NIST-grade-B gf -2.310 that RYA-353/354 had seeded) vs AFTER (the
# canonical Johansson-2003 -2.11). The before/after is in-memory only — the live
# linelist already carries the canonical value; we transiently override the one
# Ni row to show the decomposition. No hardcoded gf in the live path.

_NI_6300_WL = 6300.342          # Å (air)
_NI_6300_EP = 4.266             # eV
_STALE_NI_GF = -2.310           # pre-RYA-365 canonical (NIST ASD grade B) — diagnostic only
_OI_PART_WIN_A = (6300.00, 6300.55)   # tight window over the [O I]+Ni blend core for EW


def _ni6300_idx(ll) -> int:
    """Row index of Ni I 6300.34 in the synth linelist (raises if absent)."""
    w = ll['wave_A'].astype(float)
    ep = ll['lower_state_eV'].astype(float)
    for i in np.where((np.abs(w - _NI_6300_WL) <= 0.02) & (np.abs(ep - _NI_6300_EP) <= 0.05))[0]:
        if str(ll['element'][i]).strip().startswith('Ni') and int(ll['ion'][i]) == 1:
            return int(i)
    raise RuntimeError(
        "Ni I 6300.34 absent from the synth linelist — cannot run the [O I] blend "
        "partition (RYA-365). The joint fit requires the Ni blend partner.")


def _ew_mA(sw_nm, flux) -> float:
    """Equivalent width in mÅ of (1 - normalized flux) over sw_nm (a nm grid)."""
    return float(_trapezoid(1.0 - np.asarray(flux), np.asarray(sw_nm)) * 1.0e4)


def oi_blend_partition(star_id: str = 'solar', *, tmp_dir: str = '/tmp/ispec_cno',
                       a_C: float = None) -> dict:
    """Decompose the [O I] 6300.30 + Ni I 6300.34 blend and fit A(O), BEFORE
    (stale Ni gf -2.310) vs AFTER (canonical gf_resolver = Johansson+2003 -2.11).

    C/N pinned at the solar anchor (CO coupling to [O I] is negligible at solar);
    A(Ni) pinned from the canonical solar Ni abundance (Asplund 2021, sourced),
    not a literal. Returns {'before': {...}, 'after': {...}, 'canon_gf': ...}.
    """
    region = REGIONS['vis']
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    rec = get_star_params(star_id)
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    feh = params['feh']
    off = (feh if star_id != 'solar' else 0.0)

    print(f"\n{'='*72}\n  [O I] 6300 BLEND PARTITION — {star_id} (RYA-365)\n{'='*72}")
    broadening = preflight(region, star_id, [d for d in VIS_DIAGNOSTICS if d.key == 'OI_6300'])

    atm = _load_atmosphere(params['teff_K'], params['logg'], feh, params['vturb_kms'])
    ll, iso, chem = _load_synth_resources()        # already canonical-gf (RYA-353)
    sab = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    obs_w, obs_f = _load_observed_spectrum(star_id)
    codes = _atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)

    # CRITICAL: the live Ni gf MUST equal the gf_resolver canonical (no hardcoded copy).
    ni_i = _ni6300_idx(ll)
    canon_gf = float(resolve_gf((28, 1), _NI_6300_WL, _NI_6300_EP))
    live_gf = float(ll['loggf'][ni_i])
    if abs(live_gf - canon_gf) > 1e-4:
        raise AssertionError(
            f"Ni I 6300.34 synth gf {live_gf:+.4f} != gf_resolver canonical "
            f"{canon_gf:+.4f} — the [O I] path is NOT resolving via gf_resolver "
            f"(RYA-365 invariant violated).")
    print(f"  Ni I 6300.34 gf via gf_resolver = {canon_gf:+.3f} "
          f"(Johansson+2003); live synth row matches ✓")

    a_Ni = float(SOLAR_ASPLUND2021['Ni']) + off      # pinned, canonical solar Ni
    a_C0 = float(a_C if a_C is not None else SOLAR_ASPLUND2021['C']) + (0.0 if a_C is not None else off)
    a_N0 = float(SOLAR_ASPLUND2021['N']) + off
    print(f"  pinned: A(Ni)={a_Ni:.2f} (Asplund2021 solar Ni)  A(C)={a_C0:.2f}  A(N)={a_N0:.2f}")

    win = (_OI_PART_WIN_A,)
    fit_win = ((6299.5, 6301.0),)        # same fit window run_cno uses for OI_6300
    sw = np.arange(_OI_PART_WIN_A[0] / 10.0, _OI_PART_WIN_A[1] / 10.0 + _WSTEP_NM * 0.5,
                   _WSTEP_NM)

    out = {'canon_gf': canon_gf}
    for tag, ni_gf in (('before', _STALE_NI_GF), ('after', canon_gf)):
        ll_v = ll.copy()
        ll_v['loggf'][ni_i] = ni_gf
        seed_O = float(SOLAR_ASPLUND2021['O']) + off
        rfit = _fit_element(obs_w, obs_f, atm, params, 'O',
                            {'C': a_C0, 'N': a_N0, 'O': seed_O, 'Ni': a_Ni}, codes,
                            fit_win, True, broadening,
                            seed_O - 1.2, seed_O + 1.2, ll_v, iso, sab, tmp_dir)
        a_O = rfit['A_X']
        st_full = {'C': a_C0, 'N': a_N0, 'O': a_O, 'Ni': a_Ni}
        st_noNi = dict(st_full, Ni=a_Ni - 6.0)       # remove Ni (6 dex weaker → ~0)
        st_noO = dict(st_full, O=a_O - 6.0)          # remove [O I]
        f_full = _synth_window(sw, atm, params, ll_v, iso, sab, _fixed_ab(st_full, codes),
                               broadening, True, tmp_dir)
        f_noNi = _synth_window(sw, atm, params, ll_v, iso, sab, _fixed_ab(st_noNi, codes),
                               broadening, True, tmp_dir)
        f_noO = _synth_window(sw, atm, params, ll_v, iso, sab, _fixed_ab(st_noO, codes),
                              broadening, True, tmp_dir)
        ew_tot = _ew_mA(sw, f_full)
        ew_ni = ew_tot - _ew_mA(sw, f_noNi)
        ew_o = ew_tot - _ew_mA(sw, f_noO)
        depth = float(1.0 - np.min(f_full))
        out[tag] = {'ni_gf': ni_gf, 'A_O': a_O, 'red_chi2': rfit['red_chi2'],
                    'ew_total_mA': round(ew_tot, 2), 'ew_Ni_mA': round(ew_ni, 2),
                    'ew_OI_mA': round(ew_o, 2), 'core_depth': round(depth, 3),
                    'ni_frac': round(ew_ni / ew_tot, 3) if ew_tot else float('nan')}
    return out


def print_oi_partition(part: dict) -> bool:
    """Print the before/after blend-partition table; return whether O clears the gate."""
    tgt, terr = SOLAR_VIS_GATES['O']
    print(f"\n  [O I] 6300 blend partition (EW over {_OI_PART_WIN_A[0]}-{_OI_PART_WIN_A[1]} Å):")
    print(f"  {'case':<26}{'Ni gf':>8}{'EW_tot':>9}{'EW_Ni':>8}{'EW_[OI]':>9}"
          f"{'Ni frac':>9}{'A(O)':>8}{'χ²ᵣ':>8}")
    for tag, label in (('before', 'before (NIST B)'), ('after', 'after (Johansson03)')):
        p = part[tag]
        print(f"  {label:<26}{p['ni_gf']:>+8.3f}{p['ew_total_mA']:>9.2f}"
              f"{p['ew_Ni_mA']:>8.2f}{p['ew_OI_mA']:>9.2f}{p['ni_frac']:>9.3f}"
              f"{p['A_O']:>8.3f}{p['red_chi2']:>8.3f}")
    aO = part['after']['A_O']
    within = abs(aO - tgt) <= terr + 1e-9
    print(f"\n  A(O)☉: before={part['before']['A_O']:.3f}  →  after={aO:.3f}  "
          f"(target {tgt:.2f} ± {terr:.2f}, Δ={aO - tgt:+.3f})")
    print(f"  GATE: {'PASS — solar O cleared' if within else 'FAIL — RCA finding, do NOT force (see RYA-354)'}")
    return within


# ── CLI ───────────────────────────────────────────────────────────────────────

# ── Phase-A multi-arm orchestration (RYA-371) ─────────────────────────────────

def _fit_arm(region, diagnostics, obs_w, obs_f, params, fixed_state, atm, ll, iso, sab,
             codes, broadening, tmp_dir):
    """Fit each diagnostic of an arm on its PRE-LOADED co-added rest-frame spectrum
    (obs_w, obs_f — resolved per-star by resolve_arm_spectrum, RYA-464; the free element
    varies, fixed_state holds the rest, e.g. A(C) from HARPS for CN). Returns
    (per_band, corrections) — corrections via the same cited layer as HARPS."""
    per_band = []
    for d in diagnostics:
        st = dict(fixed_state)
        center = float(st.get(d.element, SOLAR_ASPLUND2021[d.element]))
        t0 = time.time()
        r = _fit_element(obs_w, obs_f, atm, params, d.element, st, codes,
                         d.windows_A, d.use_molecules, broadening,
                         center - 1.0, center + 1.0, ll, iso, sab, tmp_dir)
        r.update(key=d.key, element=d.element, role=d.role,
                 nlte_flag=d.nlte_flag, nlte_ref=d.nlte_ref,
                 wall_s=round(time.time() - t0, 1))
        per_band.append(r)
        print(f"    {d.key:10s} A({d.element})={r['A_X']}  χ²ᵣ={r['red_chi2']}  "
              f"[{r['status']}]  ({r['wall_s']}s)")
    return per_band, apply_cited_corrections(per_band, params, region)


def _print_cross_arm_table(per_arm: dict) -> dict:
    """Cross-arm CNO agreement map (RYA-455 O handling). Per element: the reconciled
    A(X) from each arm/indicator (cited corrections applied), primary vs cross-check,
    the spread, and a verdict — disagreement is REPORTED, never averaged. Returns a
    summary dict (the differential-backbone seed)."""
    anchors = {'C': 8.46, 'N': 7.83, 'O': 8.69}     # Asplund 2021 validation targets
    rows = {'C': [], 'N': [], 'O': []}
    for arm, d in per_arm.items():
        for c in d.get('corrections', []):
            el = c.get('element')
            if el in rows and c.get('a_corr') is not None and np.isfinite(c['a_corr']):
                rows[el].append((arm, c['key'], c['role'], float(c['a_corr']), c['kind']))
    print(f"\n{'='*72}\n  CROSS-ARM CNO AGREEMENT — cited corrections applied "
          f"(validate-don't-tune)\n{'='*72}")
    summary = {}
    for el in ('C', 'N', 'O'):
        print(f"\n  {el}  (Asplund {anchors[el]}):")
        if not rows[el]:
            print("    (no reconciled indicator this run)")
            summary[el] = {'indicators': [], 'verdict': 'NO-DATA'}
            continue
        for arm, key, role, a, kind in rows[el]:
            print(f"    {arm:9s} {key:11s} {role:11s} A({el})={a:.3f}  "
                  f"Δ={a - anchors[el]:+.3f}  [{kind}]")
        prim = [a for _, _, role, a, _ in rows[el] if role == 'primary']
        allv = [a for *_, a, _ in rows[el]]
        spread = max(allv) - min(allv)
        verdict = 'AGREE (≤0.10)' if spread <= 0.10 else 'FLAGGED-DISAGREEMENT (reported, not averaged)'
        pnote = (f"primary mean {np.mean(prim):.3f} (Δ {np.mean(prim) - anchors[el]:+.3f})"
                 if prim else "NO primary indicator this arm-set")
        print(f"    → {len(allv)} indicator(s), spread {spread:.3f} dex; {pnote}; {verdict}")
        summary[el] = {'indicators': [{'arm': a, 'key': k, 'role': ro, 'A': v, 'kind': ki}
                                      for a, k, ro, v, ki in rows[el]],
                       'spread': round(spread, 3), 'primary_mean': (round(float(np.mean(prim)), 3) if prim else None),
                       'verdict': verdict}
    return summary


def run_phase_a(star_id: str = 'solar', arms=None,
                tmp_dir: str = '/tmp/ispec_cno', out_dir: Path = None) -> dict:
    """Phase-A multi-arm CNO, PER-STAR (RYA-464 generalization of the RYA-371 solar path).
    Arms are resolved from the star's STAR_ARMS registry, not a hardcoded list: each arm
    declares its region, diagnostics, spectrum loader, and readiness. Ready arms run; arms
    whose loader/data/audit aren't present are DEFERRED with their reason (never silently
    run against reflected-solar geometry). `arms` (names) restricts the run; None = the
    star's full declared set. O handling per RYA-455: O I 777 primary, [O I] 6300 cross-check.

    Solar is one case of this mechanism (harps + espresso + uves all ready via the existing
    loaders) → bit-identical to the pre-RYA-464 solar run."""
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    registry = star_arm_registry(star_id)
    requested = tuple(a.strip().lower() for a in arms) if arms else tuple(registry)

    rec = get_star_params(star_id)
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    atm = _load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    ll, iso, chem = _load_synth_resources()
    sab = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    codes = _atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)
    _, vmac, vsini, _ = _resolve_broadening(star_id)

    # Announce the per-star arm map (anti-silent-assumption, RYA-270/464 discipline).
    _ready, _deferred = available_arms(star_id)
    print(f"\n{'='*72}\n  PHASE-A ARM REGISTRY — {star_id} "
          f"({len(registry)} declared region(s))\n{'='*72}")
    for n, a in registry.items():
        tag = 'READY' if a.ready else f'DEFERRED — {a.defer_reason}'
        print(f"  [{ 'x' if a.ready else ' ' }] {n:9s} {a.region.instrument:9s} "
              f"{a.region.wave_min_A:.0f}-{a.region.wave_max_A:.0f} A  {tag}")

    per_arm = {}
    deferred = {}
    # HARPS first — the molecular CNO-equilibrium engine (only if requested + ready).
    if 'harps' in requested and 'harps' in registry and registry['harps'].ready:
        print(f"\n{'#'*72}\n#  ARM — HARPS (CH/CN/[O I] + C I), molecular CNO equilibrium\n{'#'*72}")
        h = run_cno(star_id, 'vis', tmp_dir=tmp_dir)
        per_arm['harps'] = {'abundances': h.abundances, 'corrections': h.phase_a_corrections}

    fixed = {'Ni': 6.20}
    for el in ('C', 'N', 'O'):                       # CN needs A(C); pin from HARPS, else anchor
        fixed[el] = float(per_arm.get('harps', {}).get('abundances', {}).get(el, SOLAR_ASPLUND2021[el]))

    for name in requested:
        if name == 'harps':
            continue
        arm = registry.get(name)
        if arm is None:
            deferred[name] = 'not declared for this star (add to STAR_ARMS, RYA-464)'
            print(f"\n#  ARM — {name}: SKIPPED — {deferred[name]}")
            continue
        if not arm.ready:                            # loud defer, NEVER fall back to Vesta
            deferred[name] = arm.defer_reason
            print(f"\n#  ARM — {arm.region.instrument} ({name}): DEFERRED — {arm.defer_reason}")
            continue
        print(f"\n{'#'*72}\n#  ARM — {arm.region.instrument} ({name}; pinned "
              f"A(C)={fixed['C']:.2f})\n{'#'*72}")
        obs_w, obs_f = resolve_arm_spectrum(star_id, arm)   # per-star loader dispatch
        pb, corr = _fit_arm(arm.region, arm.diagnostics, obs_w, obs_f, params, fixed,
                            atm, ll, iso, sab, codes, (arm.region.R, vmac, vsini), tmp_dir)
        per_arm[name] = {'per_band': pb, 'corrections': corr}

    summary = _print_cross_arm_table(per_arm)
    if deferred:
        print(f"\n  DEFERRED ARMS ({len(deferred)}): "
              + "; ".join(f"{n} ({r})" for n, r in deferred.items()))

    out_dir = Path(out_dir) if out_dir else (Path(PATHS['solar_ew']).parent.parent /
                                             'audit' / 'cno_synthesis')
    out_dir.mkdir(parents=True, exist_ok=True)
    rep = {'ticket': 'RYA-371 Phase A / RYA-464 per-star arms', 'star': star_id,
           'arms_requested': list(requested),
           'arms_run': [a for a in per_arm], 'arms_deferred': deferred,
           'o_handling': 'RYA-455: O I 777 primary, [O I] 6300 continuum-limited cross-check',
           'cross_arm': summary,
           'per_arm': {a: d.get('corrections', []) for a, d in per_arm.items()}}
    out_name = ('solar_phase_a_cross_arm.json' if 'solar' in star_id.lower()
                else f'{star_id}_phase_a_cross_arm.json')
    (out_dir / out_name).write_text(json.dumps(rep, indent=2, default=str))
    print(f"\n  [out] {out_dir / out_name}")
    return per_arm


def main(argv=None):
    ap = argparse.ArgumentParser(description='Region-aware C/N/O synthesis (RYA-237)')
    ap.add_argument('--star', default='solar')
    ap.add_argument('--region', default='vis', choices=sorted(REGIONS))
    ap.add_argument('--species', default=None,
                    help='restrict to a channel, e.g. "O I" → the [O I] 6300 blend '
                         'partition (Ni vs [O I], before/after canonical gf). '
                         'Omit for the full C/N/O run.')
    ap.add_argument('--validate', action='store_true',
                    help='print the solar-VIS acceptance gate table')
    ap.add_argument('--no-systematics', action='store_true',
                    help='skip the Type-B stellar-parameter sensitivity refits')
    ap.add_argument('--pin', action='append', default=None, metavar='EL=VALUE',
                    help="RYA-1214: pin an element's abundance as an INPUT instead of "
                         "seeding it from the solar table. A region without a primary for "
                         "that element does not measure it, and its value still enters the "
                         "molecular equilibrium of the ones that ARE measured — so the "
                         "near-UV OH/NH run wants A(C) from the VIS CH G-band rather than "
                         "from a solar default. Recorded as PINNED_INPUT_NOT_MEASURED in "
                         "the product, never as an abundance.")
    ap.add_argument('--max-iter', type=int, default=5)
    ap.add_argument('--out', default=None)
    ap.add_argument('--arms', default=None,
                    help='Phase-A multi-arm run (RYA-371/464): comma-separated arm names '
                         'from the STAR\'s registry (e.g. "harps,uves"), or "all" for the '
                         'star\'s full declared set. Per-arm C/N/O + cross-arm agreement; '
                         'arms whose loader/data are not ready are DEFERRED, not run. Omit '
                         'for the single HARPS-VIS run.')
    ap.add_argument('--list-arms', action='store_true',
                    help='print the star\'s per-star arm registry (ready / deferred) and exit')
    args = ap.parse_args(argv)

    # --list-arms → just report the per-star arm registry (RYA-464) and exit.
    if args.list_arms:
        reg = star_arm_registry(args.star)
        ready, deferred = available_arms(args.star)
        print(f"\nArm registry — {args.star} ({len(reg)} declared region(s)):")
        for n, a in reg.items():
            tag = 'READY' if a.ready else f'DEFERRED — {a.defer_reason}'
            print(f"  {n:9s} {a.region.instrument:9s} "
                  f"{a.region.wave_min_A:.0f}-{a.region.wave_max_A:.0f} A  {tag}")
        print(f"\n  ready: {list(ready)}   deferred: {list(deferred)}")
        return {'ready': ready, 'deferred': deferred}

    # --arms → Phase A per-star multi-arm CNO (RYA-371/464), validated vs the star's registry.
    if args.arms:
        reg = star_arm_registry(args.star)
        if args.arms.strip().lower() == 'all':
            arms = tuple(reg)
        else:
            arms = tuple(a.strip().lower() for a in args.arms.split(',') if a.strip())
            unknown = [a for a in arms if a not in reg]
            if unknown:
                raise SystemExit(f"Unknown arm(s) {unknown} for {args.star}; declared arms: "
                                 f"{list(reg)} (STAR_ARMS, RYA-464).")
        return run_phase_a(args.star, arms=arms,
                           out_dir=Path(args.out) if args.out else None)

    # --species "O I" → the focused [O I] 6300 blend-partition diagnostic (RYA-365).
    if args.species and args.species.replace(' ', '').upper() in ('OI', 'O'):
        if args.region != 'vis':
            raise SystemExit("[O I] 6300 blend partition is a VIS diagnostic (--region vis).")
        part = oi_blend_partition(args.star, tmp_dir='/tmp/ispec_cno')
        if args.validate:
            print_oi_partition(part)
        return part

    pins = {}
    for kv in (args.pin or ()):
        el, _, val = str(kv).partition('=')
        el = el.strip()
        if el not in ('C', 'N', 'O', 'Ni') or not val:
            raise SystemExit(f"--pin {kv!r} is not EL=VALUE for one of C / N / O / Ni")
        pins[el] = float(val)
    result = run_cno(args.star, args.region, max_iter=args.max_iter,
                     with_systematics=not args.no_systematics,
                     out_dir=Path(args.out) if args.out else None,
                     pins=pins or None)
    if args.validate and args.star == 'solar':
        validate_solar(result)
    elif args.validate:
        print(f"\n  --validate: gate table is defined for solar only; "
              f"{args.star} compared against GBS downstream (RYA-348).")
    return result


if __name__ == '__main__':
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        main()

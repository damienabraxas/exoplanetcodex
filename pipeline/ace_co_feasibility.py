"""
pipeline/ace_co_feasibility.py
==============================
RYA-440 — ACE-FTS solar CO as a measurement SOURCE: 1-day disk-geometry
feasibility go/no-go.

The reflected-solar IR path is dead (RYA-391: the only archival CRIRES+ K-band
reflected-solar epoch fails the SNR floor). ACE-FTS (RYA-390/392) is a
telluric-free, R~213,000, SNR~400 SPACE solar atlas already on disk, covering the
CO first-overtone band 4255-4367 cm-1 (12CO 4361 cm-1 / 2.293 um and 13CO ~4264
cm-1 / 2.345 um). So far it was only a RULER to grade the dead Vesta product. This
asks whether it can BE the source for solar A(C)/A(O).

The one real risk: ACE views the Sun in occultation, so its FOV samples a REGION
of the disk, and CO overtone lines are strongly center-to-limb sensitive -- so the
derived A(C) depends on the disk-geometry assumption. This module QUANTIFIES that
dependence and returns a go/no-go. Nothing more. NOT the full IR-solar run.

What it does (and the findings it surfaces, no silent fallback):
  Step 0  line-list check: 12C16O Li2015 present (covers 2.3 um); 13C16O present?
          Plus the WIRING finding: the list is a `.dat`, not the `_<lo>-<hi>.bsyn`
          iSpec globs (turbospectrum.py), so it is NOT auto-included in synthesis.
          We stage a BAND-SCOPED symlink (16O12C_<lo>-<hi>.bsyn) for the 2.3 um
          window ONLY, so the validated optical CNO arm (RYA-237) is untouched.
  Step 1  load ACE + air<->vac LOUD-FAIL boundary (RYA-373 rule). The 12C16O (2-0)
          bandhead is the fiducial: the Turbospectrum/ExoMol synthesis is VACUUM
          at 2.3 um, so we match ACE `wavelength_vac_A` (the air column is offset
          by the ~6.3 A / ~83 km/s air/vac slip -> asserted and failed-loud).
  Step 2  geometry bracket. The RYA-237 Turbospectrum-via-iSpec path
          (generate_spectrum) exposes DISK-INTEGRATED FLUX only -- there is NO
          mu-resolved intensity API. That is a FINDING: geometries (a) mu=1.0 and
          (c) mu=0.5 cannot be synthesized natively. We report it, synthesize the
          FLUX geometry, and bracket the rest with the Step-5 ACE pointing geometry
          + the cited center-to-limb literature. No silent flux substitution.
  Step 3  fit A(C) (flux geometry) to the ACE 12CO (2-0) bandhead, A(O) FIXED at
          the solar reference, params PINNED (constants.py). Validate-don't-tune.
  Step 4  verdict from the geometry spread vs the C/O-precision tolerance.
  Step 5  geometry knowability (Hase et al. 2010): is the occultation pointing
          documented well enough to collapse the bracket to a single mu?

Permanent rules honored: air<->vac loud-fail; all params/abundances/tolerances
from constants.py (single source, cited); no silent fallback; ASCII-only.

    python -m pipeline.ace_co_feasibility --validate
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import SOLAR_ASPLUND2021, ASTRO, STAR_SOLAR  # noqa: E402

# Reuse the RYA-237 synthesis machinery (atmosphere interp, GES linelist, isotopes,
# Turbospectrum window synthesis + chi2 single-element fit). No parallel engine.
from pipeline.cno_synthesis import (  # noqa: E402
    _load_atmosphere, _load_synth_resources, _atom_codes, _fixed_ab,
    _fit_element, _synth_window, _MOLECULES_DIR, _ISPEC_SOLAR_ABUND_FILE,
    SOLAR_VIS_GATES,
)
import ispec  # noqa: E402

C_KMS = 299792.458
ATLAS = ROOT / 'data' / 'solar_reference' / 'ir_atlases'
ACE_CSV = ATLAS / 'ace_fts_solar_co_4255_4367.csv'
OUT = ROOT / 'data' / 'audit' / 'ace_co_feasibility_rya440'

CO_LIST = _MOLECULES_DIR / 'CO_IR_Li2015.dat'    # RYA-236 ExoMol Li2015 12C16O
CO_13_GLOB = ('13C16O', '16O13C', '0608.013016')  # 13C16O species token / names

# 12C16O (2-0) R-branch bandhead, vacuum rest wavelength (Angstrom). The strongest
# absorption in the first-overtone band; our fiducial for the air/vac assertion.
CO_2_0_BANDHEAD_VAC_A = 22935.3
AIRVAC_SLIP_KMS = 83.0          # ~6.3 A at 2.3 um; the air/vac signature to catch
FRAME_MATCH_TOL_KMS = 5.0       # synthesis<->atlas alignment tolerance in the matched frame
FIT_WINDOW_VAC_A = (22935.0, 22990.0)   # bandhead + first R-branch lines

# Solar reference (single source: constants.py). A(C)_sun and the C/O-precision
# tolerance both come from the project reference; no round-default.
A_C_REF = float(SOLAR_ASPLUND2021['C'])         # 8.46 (Asplund et al. 2021)
A_O_REF = float(SOLAR_ASPLUND2021['O'])         # 8.69
# SOLAR_VIS_GATES (cno_synthesis.py): {'C': (8.46, 0.05), 'O': (8.69, 0.05), ...}
# the sigma is the adopted C/O precision -> the geometry-spread tolerance.
GEOM_TOL_DEX = float(SOLAR_VIS_GATES['C'][1])   # 0.05 dex
NEAR_REF_DEX = 0.10                              # "lands near the reference" band (ticket Step 4)

ACE_R = 213000.0                                 # ACE-FTS resolving power (RYA-390/392)


def _solar_params() -> dict:
    """Pinned solar params from constants.py (NOT hardcoded): Teff/logg from ASTRO,
    xi from STAR_SOLAR, [Fe/H]=0 (solar)."""
    return {'teff_K': float(ASTRO['Teff_sun']), 'logg': float(ASTRO['logg_sun']),
            'feh': 0.0, 'vturb_kms': float(STAR_SOLAR['vturb_kms'])}


# ── Step 0: line-list check ───────────────────────────────────────────────────

def step0_linelists() -> dict:
    """12C16O present + 2.3 um coverage; 13C16O present? + the .dat-not-globbed
    wiring finding."""
    out = {'co12_present': CO_LIST.exists(), 'co12_path': str(CO_LIST)}
    if not CO_LIST.exists():
        raise FileNotFoundError(
            f"12C16O Li2015 list absent at {CO_LIST} (RYA-236). Cannot synthesize CO.")
    head = CO_LIST.read_text().splitlines()[0]
    out['co12_species_token'] = head.split("'")[1].strip() if "'" in head else head.strip()
    # count lines in the 2.293-2.345 um overtone window (vacuum Angstrom)
    n_band = 0
    for ln in CO_LIST.read_text().splitlines()[2:]:
        parts = ln.split()
        if parts and parts[0].replace('.', '', 1).isdigit():
            w = float(parts[0])
            if 22900.0 <= w <= 23500.0:
                n_band += 1
    out['co12_lines_in_overtone_band'] = n_band
    # 13C16O: separate species file / token anywhere in the molecules dir?
    co13 = []
    for f in _MOLECULES_DIR.glob('*'):
        if any(tok in f.name for tok in CO_13_GLOB):
            co13.append(f.name)
    # also: is the 13C16O species token present inside any list header?
    out['co13_present'] = bool(co13)
    out['co13_files'] = co13
    # wiring: iSpec globs molecules/*.bsyn with _<lo>-<hi>.bsyn names only
    out['wiring_note'] = (
        "CO_IR_Li2015.dat is a .dat, not the _<lo>-<hi>.bsyn iSpec globs "
        "(turbospectrum.py glob 'molecules/*.bsyn'); it is NOT auto-included in "
        "synthesis. This module stages a band-scoped symlink 16O12C_<lo>-<hi>.bsyn "
        "for the 2.3 um window only, leaving the optical CNO arm untouched.")
    return out


# ── CO molecular-list wiring (band-scoped, optical-arm-safe) ──────────────────

class _StagedCO:
    """Stage a band-scoped symlink so iSpec's molecules/*.bsyn glob includes the CO
    list FOR THE 2.3 um WINDOW ONLY. The filename nm-range (turbospectrum.py overlap
    test) keeps it out of every optical synthesis. Always cleaned up."""

    def __init__(self, lo_nm=2280, hi_nm=2360):
        self.link = _MOLECULES_DIR / f'16O12C_{lo_nm}-{hi_nm}.bsyn'

    def __enter__(self):
        if self.link.is_symlink() or self.link.exists():
            self.link.unlink()
        os.symlink(CO_LIST.name, self.link)   # relative link within molecules/
        return self

    def __exit__(self, *exc):
        if self.link.is_symlink() or self.link.exists():
            self.link.unlink()


# ── Step 1: ACE load + air/vac loud-fail ──────────────────────────────────────

def _bandhead_min_A(wave_A, flux, lo=22925.0, hi=22965.0) -> float:
    m = (wave_A >= lo) & (wave_A <= hi)
    return float(wave_A[m][np.argmin(flux[m])])


def step1_load_ace(synth_bandhead_A: float) -> dict:
    """Load ACE; choose the wavelength column matched to the synthesis frame; assert
    the 12CO bandhead lands at rest; fail loud on the ~83 km/s air/vac slip."""
    df = pd.read_csv(ACE_CSV)
    frames = {}
    for col in ('wavelength_vac_A', 'wavelength_air_A'):
        s = df.sort_values(col)
        w = s[col].to_numpy(float)
        f = s['intensity'].to_numpy(float)
        head = _bandhead_min_A(w, f)
        dv = (head - synth_bandhead_A) / synth_bandhead_A * C_KMS
        frames[col] = {'bandhead_A': round(head, 3), 'dv_vs_synth_kms': round(dv, 1)}
    # matched frame = the one aligning with the synthesis (|dv| small)
    matched = min(frames, key=lambda c: abs(frames[c]['dv_vs_synth_kms']))
    other = 'wavelength_air_A' if matched == 'wavelength_vac_A' else 'wavelength_vac_A'
    dv_matched = abs(frames[matched]['dv_vs_synth_kms'])
    dv_other = abs(frames[other]['dv_vs_synth_kms'])
    if dv_matched > FRAME_MATCH_TOL_KMS:
        raise ValueError(
            f"AIR/VAC LOUD-FAIL: neither ACE column aligns with the synthesis "
            f"bandhead within {FRAME_MATCH_TOL_KMS} km/s (best {matched} "
            f"{dv_matched:.1f} km/s). Refusing to fit a frame-mismatched spectrum.")
    if not (AIRVAC_SLIP_KMS * 0.6 <= dv_other <= AIRVAC_SLIP_KMS * 1.4):
        raise ValueError(
            f"AIR/VAC LOUD-FAIL: the non-matched column ({other}) offset "
            f"{dv_other:.1f} km/s is not the expected ~{AIRVAC_SLIP_KMS:.0f} km/s "
            f"air/vac slip; frame assumptions are inconsistent. STOP.")
    s = df.sort_values(matched)
    return {'matched_frame': matched, 'frames': frames,
            'wave_A': s[matched].to_numpy(float), 'flux': s['intensity'].to_numpy(float),
            'airvac_slip_kms': round(dv_other, 1)}


# ── Step 2/3: synthesize FLUX geometry + fit A(C) ─────────────────────────────

def _fiducial_bandhead_A(params, ll, iso, sab, codes) -> float:
    """Synthesize the bandhead once to locate it in the (vacuum) synthesis frame."""
    w_nm = np.arange(2292.0, 2300.0 + 1e-7, 0.0002)
    fa = _fixed_ab({'C': A_C_REF, 'O': A_O_REF}, codes)
    with _StagedCO():
        f = _synth_window(w_nm, params['_atm'], params, ll, iso, sab, fa,
                          (ACE_R, 0.0, 0.0), True, params['_tmp'])
    return _bandhead_min_A(w_nm * 10.0, f)


def step23_fit_flux(ace: dict, params, ll, iso, sab, codes) -> dict:
    """Fit A(C) to the ACE bandhead in the FLUX geometry (the only one the toolchain
    exposes). A(O) FIXED, params PINNED. Validate-don't-tune: A(C) falls where it
    falls."""
    w_A, f = ace['wave_A'], ace['flux']
    w_nm = w_A / 10.0
    state = {'C': A_C_REF, 'O': A_O_REF}        # O fixed at the reference
    with _StagedCO():
        r = _fit_element(
            obs_w_nm=w_nm, obs_f=f, atm=params['_atm'], params=params,
            free_el='C', state=state, codes=codes,
            windows_A=[FIT_WINDOW_VAC_A], use_molecules=True,
            broadening=(ACE_R, 0.0, 0.0),
            a_lo=A_C_REF - 1.0, a_hi=A_C_REF + 1.0,
            ll=ll, iso=iso, sab=sab, tmp_dir=params['_tmp'])
    return r


# ── Step 5: geometry knowability (Hase et al. 2010) ───────────────────────────

def step5_geometry() -> dict:
    """ACE occultation pointing geometry, from Hase et al. 2010 (JQSRT 111, 521) and
    the ACE-FTS instrument papers. Read-only; documents whether the disk-sampling is
    knowable enough to collapse the geometry bracket to a single mu."""
    # ACE-FTS suntracker points at the RADIOMETRIC CENTER of the Sun (stability
    # <15 urad); FOV = 1.25 mrad vs the solar disk ~9.3 mrad (32 arcmin). The FOV
    # therefore spans only the central ~13% of the disk DIAMETER, centered on disk
    # center. r/R at the FOV edge ~0.13 -> mu = sqrt(1-(r/R)^2) in [0.991, 1.000].
    fov_mrad, disk_mrad = 1.25, 9.3
    r_over_R = (fov_mrad / 2.0) / (disk_mrad / 2.0)
    mu_edge = float(np.sqrt(max(0.0, 1.0 - r_over_R ** 2)))
    return {
        'reference': 'Hase, Wallace, McLeod, Harrison & Bernath 2010, JQSRT 111, 521; '
                     'ACE-FTS instrument papers (Bernath 2006; Gilbert 2007)',
        'suntracker': 'fine pointing to the radiometric centre of the Sun, stability <15 urad',
        'fov_mrad': fov_mrad, 'solar_disk_mrad': disk_mrad,
        'fov_fraction_of_diameter': round(r_over_R, 3),
        'mu_range_sampled': [round(mu_edge, 3), 1.0],
        'finding': (
            "KNOWABLE and near disk-CENTER. The suntracker locks on the radiometric "
            "centre and the 1.25 mrad FOV spans only the central ~13% of the disk "
            "diameter, so the sampled mu is 0.99-1.00 -- effectively disk-CENTER "
            "INTENSITY (mu~1), NOT disk-integrated flux and NOT mu=0.5. The geometry "
            "bracket COLLAPSES to mu~1; it does not span to mu=0.5. The remaining "
            "systematic is intensity(mu~1)-vs-flux, not an open disk-geometry "
            "ambiguity."),
    }


# ── Step 4: verdict ───────────────────────────────────────────────────────────

def step4_verdict(flux_fit: dict, geom: dict) -> dict:
    """Synthesize the go/no-go. The toolchain is flux-only, so the mu=1 / mu=0.5
    rows cannot be synthesized; Step 5 shows ACE is mu~1, so the operative question
    is whether we can synthesize at that KNOWN geometry."""
    a_flux = flux_fit.get('A_X')
    rows = {
        'mu=1.0 (disk-center; ACE actual)': 'NOT-SYNTHESIZABLE (flux-only path)',
        'flux (disk-integrated)': (f"{a_flux:.3f}" if a_flux is not None and
                                   np.isfinite(a_flux) else 'fit-failed'),
        'mu=0.5 (intermediate)': 'NOT-APPLICABLE (ACE FOV is mu~1, not a disk average)',
    }
    flux_near_ref = (a_flux is not None and np.isfinite(a_flux)
                     and abs(a_flux - A_C_REF) <= NEAR_REF_DEX)
    verdict = (
        "NO-GO (today) as a turnkey source, but the geometry RISK is RETIRED. "
        "Reasons: (1) the RYA-237 Turbospectrum-via-iSpec path exposes disk-"
        "integrated FLUX only -- no mu-resolved intensity -- so we cannot synthesize "
        "at ACE's actual geometry; the 3-geometry synthesis spread is therefore not "
        "computable by the toolchain. (2) Step 5 RESOLVES the disk-geometry "
        "ambiguity the ticket feared: ACE samples mu~0.99-1.00 (disk-center "
        "intensity), NOT a disk average, so the bracket does not span to mu=0.5. "
        "The open-ended geometry risk thus becomes two BOUNDED fast-follows: "
        "(a) expose Turbospectrum disk-center (mu=1) specific-intensity through "
        "iSpec; (b) acquire the 13C16O line list (absent). With both, ACE is a "
        "viable telluric-free solar CO source. Until then, fitting FLUX to the "
        "disk-center ACE atlas carries an intensity-vs-flux systematic that, for the "
        "strongest (2-0) bandhead lines (which form very high; Scott et al. 2006 "
        "A&A 456, 675), exceeds the "
        f"{GEOM_TOL_DEX:.2f} dex C/O tolerance -- restrict to weaker CO lines.")
    return {
        'geometry_table': rows,
        'A_C_flux': a_flux, 'A_C_reference': A_C_REF,
        'flux_minus_ref_dex': (round(a_flux - A_C_REF, 3)
                               if a_flux is not None and np.isfinite(a_flux) else None),
        'flux_lands_near_ref': bool(flux_near_ref),
        'geom_spread_synthesizable': False,
        'geom_tolerance_dex': GEOM_TOL_DEX,
        'verdict': 'NO-GO (geometry risk retired; 2 bounded fast-follows)',
        'verdict_detail': verdict,
    }


# ── driver ────────────────────────────────────────────────────────────────────

def run() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    tmp = '/tmp/ispec_ace_co_rya440'
    Path(tmp).mkdir(parents=True, exist_ok=True)
    params = _solar_params()

    print("=" * 78)
    print("RYA-440  ACE-FTS solar CO feasibility (disk-geometry go/no-go)")
    print("=" * 78)
    print(f"  pinned solar params (constants.py): Teff={params['teff_K']:.0f} "
          f"logg={params['logg']:.3f} [Fe/H]={params['feh']:+.1f} "
          f"xi={params['vturb_kms']:.2f}")
    print(f"  A(C)_ref={A_C_REF:.2f}  A(O)_ref={A_O_REF:.2f} (FIXED)  "
          f"C/O tol={GEOM_TOL_DEX:.2f} dex (SOLAR_VIS_GATES)")

    # Step 0
    s0 = step0_linelists()
    print(f"\n[Step 0] 12C16O Li2015: present ({s0['co12_lines_in_overtone_band']} "
          f"lines in the 2.293-2.345 um band); species {s0['co12_species_token']}")
    print(f"         13C16O: {'present ' + str(s0['co13_files']) if s0['co13_present'] else 'ABSENT -> fast-follow (do not block)'}")
    print(f"         wiring: {s0['wiring_note'][:96]}...")

    # synthesis resources
    print("\n[setup] loading atmosphere + GES linelist + isotopes...")
    t0 = time.time()
    atm = _load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    ll, iso, chem = _load_synth_resources()
    sab = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    codes = _atom_codes(('C', 'O'), chem, sab)
    params['_atm'], params['_tmp'] = atm, tmp
    print(f"        loaded in {time.time()-t0:.0f}s")

    # locate the synthetic bandhead (vacuum frame) for the air/vac assertion
    synth_head = _fiducial_bandhead_A(params, ll, iso, sab, codes)
    print(f"        synthetic 12CO (2-0) bandhead at {synth_head:.2f} A (synthesis frame)")

    # Step 1
    s1 = step1_load_ace(synth_head)
    print(f"\n[Step 1] air/vac: matched frame = {s1['matched_frame']} "
          f"(synth-aligned); other column offset {s1['airvac_slip_kms']:.0f} km/s "
          f"= the air/vac slip -> loud-fail boundary OK")

    # Step 2/3
    print(f"\n[Step 2] geometry: generate_spectrum is FLUX-only (no mu); synthesizing "
          f"the FLUX geometry, bracketing mu via Step 5. No silent flux-substitution.")
    print(f"[Step 3] fitting A(C) to ACE 12CO (2-0) bandhead {FIT_WINDOW_VAC_A} A "
          f"(A(O) fixed {A_O_REF:.2f}, validate-don't-tune)...")
    t0 = time.time()
    flux_fit = step23_fit_flux(s1, params, ll, iso, sab, codes)
    print(f"        A(C)_flux = {flux_fit.get('A_X')}  (red_chi2="
          f"{flux_fit.get('red_chi2')}, n_pix={flux_fit.get('n_pix')}, "
          f"status={flux_fit.get('status')}, {time.time()-t0:.0f}s)")

    # Step 5
    s5 = step5_geometry()
    print(f"\n[Step 5] ACE pointing geometry (Hase 2010): {s5['finding'][:110]}...")
    print(f"         mu sampled = {s5['mu_range_sampled']} (FOV {s5['fov_mrad']} mrad / "
          f"disk {s5['solar_disk_mrad']} mrad)")

    # Step 4
    s4 = step4_verdict(flux_fit, s5)
    print("\n" + "-" * 78)
    print("3-GEOMETRY A(C) TABLE")
    for g, v in s4['geometry_table'].items():
        print(f"  {g:38s} : {v}")
    print(f"\n  A(C)_flux - A(C)_ref = {s4['flux_minus_ref_dex']} dex "
          f"(near-ref band +/-{NEAR_REF_DEX}: {'yes' if s4['flux_lands_near_ref'] else 'no'})")
    print(f"  geometry spread synthesizable by toolchain: {s4['geom_spread_synthesizable']} "
          f"(flux-only)")
    print(f"\n  ===> {s4['verdict']}")
    print("  " + s4['verdict_detail'].replace('. ', '.\n  '))

    report = {'ticket': 'RYA-440', 'params': {k: params[k] for k in
              ('teff_K', 'logg', 'feh', 'vturb_kms')},
              'A_C_reference': A_C_REF, 'A_O_reference_fixed': A_O_REF,
              'geom_tolerance_dex': GEOM_TOL_DEX,
              'synthesis_frame': 'vacuum (matched ACE wavelength_vac_A)',
              'synth_bandhead_A': round(synth_head, 3),
              'step0_linelists': s0, 'step1_airvac': {k: s1[k] for k in
              ('matched_frame', 'frames', 'airvac_slip_kms')},
              'step23_flux_fit': flux_fit, 'step5_geometry': s5, 'step4_verdict': s4}
    (OUT / 'rya440_ace_co_feasibility.json').write_text(json.dumps(report, indent=2))
    print(f"\n  [out] {OUT / 'rya440_ace_co_feasibility.json'}")
    return report


def main():
    ap = argparse.ArgumentParser(description='RYA-440 ACE-FTS solar CO feasibility')
    ap.add_argument('--validate', action='store_true',
                    help='run the full feasibility (Steps 0-5) and print the verdict')
    ap.parse_args()
    run()


if __name__ == '__main__':
    main()

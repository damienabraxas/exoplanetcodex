"""
scripts/ace_co_mu1_bracket.py
=============================
RYA-441 — ACE-FTS disk-center (mu=1) intensity CO fit: the REAL go/no-go.

RYA-440 returned A(C)_flux = 8.497 but fit it in the FLUX geometry, while 440's
Step 5 proved ACE samples mu ~ 0.99-1.00 (disk-center INTENSITY): the encouraging
number is in the WRONG geometry. CO overtone lines strengthen toward the limb, so
disk-integrated flux carries STRONGER CO than mu=1 intensity at fixed A(C); fitting
the weak-CO disk-center ACE spectrum with a strong-CO flux synthesis biases A(C)
DOWN. This script refits ACE in the correct mu=1 intensity geometry. That IS the
verdict.

DISPOSABLE DIAGNOSTIC (throwaway; NOT wired into RYA-237 production).

Toolchain path (Step 0): direct Turbospectrum specific-intensity at mu=1, reached
WITHOUT touching iSpec's flux wrapper or the flux production path. iSpec drives
bsyn_lu with 'INTENSITY/FLUX:' 'Flux'; bsyn natively also computes specific
intensity. We intercept the bsyn stdin and (1) flip 'Flux' -> 'Intensity', (2)
inject a MU-POINTS file with nangles=1, mu=1.0. bsyn then writes the mu=1 intensity
as a trailing column of its RESULTFILE (bsynb.f: in intensity mode col 2 stays the
FLUX -- "output spectrum is flux spectrum!" -- and the per-mu normalized intensity
is appended). We remap that intensity column into col 2 in-place before iSpec reads
it, so the entire RYA-440 harness (continuum, convolution, chi2 A(C) fit) then runs
on mu=1 intensity. No edit to iSpec source; no change to the flux path.

Sign tripwire (Step 2): mu=1 intensity has WEAKER CO than flux, so A(C)_mu1 MUST be
HIGHER than 8.497. A lower value contradicts CO center-to-limb physics -> STOP/RCA.

Mandatory caveat (Step 3): a 1D MARCS mu=1-vs-flux bracket is a FLOOR on the real
geometry sensitivity. CO center-to-limb is largely a 3D/granulation effect that 1D
underpredicts; a small 1D bracket does NOT close the 3D uncertainty (residual flag).

    python scripts/ace_co_mu1_bracket.py --validate
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from contextlib import contextmanager
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# Reuse the RYA-440 harness (ACE load, air<->vac loud-fail, fit window, constants,
# fiducial-bandhead locator, the CO-list band-scoped symlink) verbatim.
from pipeline import ace_co_feasibility as ace  # noqa: E402
from pipeline.cno_synthesis import (  # noqa: E402
    _load_atmosphere, _load_synth_resources, _atom_codes, _fit_element,
    _ISPEC_SOLAR_ABUND_FILE,
)
import ispec  # noqa: E402
import ispec.synth.turbospectrum as ts  # noqa: E402

OUT = ROOT / 'data' / 'audit' / 'ace_co_mu1_rya441'
A_C_FLUX_440 = 8.497          # RYA-440 flux-geometry fit (single source: 440 report JSON)

_FLUX_TOKEN = b"'INTENSITY/FLUX:' 'Flux'"
_MU_POINTS_FILE = Path('/tmp/ispec_mu1_rya441/mupoints.dat')


@contextmanager
def turbospectrum_intensity_at_mu(mu: float):
    """Patch the bsyn stdin to compute specific intensity at `mu`, and remap the
    resulting per-mu intensity column into col 2 of the RESULTFILE so iSpec's
    (flux-shaped) reader transparently consumes the mu intensity. Restores on exit.
    Babsma calls (no INTENSITY/FLUX token) pass through untouched."""
    _MU_POINTS_FILE.parent.mkdir(parents=True, exist_ok=True)
    _MU_POINTS_FILE.write_text(f"1\n{mu:.4f}\n")
    inject = (b"'INTENSITY/FLUX:' 'Intensity'\n'MU-POINTS:' '"
              + str(_MU_POINTS_FILE).encode() + b"'")
    real_popen = ts.subprocess.Popen
    stats = {'fired': 0, 'remapped': 0}

    def patched(*a, **k):
        proc = real_popen(*a, **k)
        orig = proc.communicate

        def communicate(input=None, timeout=None):
            resultfile = None
            if input is not None and _FLUX_TOKEN in input:
                input = input.replace(_FLUX_TOKEN, inject)
                stats['fired'] += 1
                m = re.search(rb"'RESULTFILE :' '([^']+)'", input)
                resultfile = m.group(1).decode() if m else None
            out = orig(input=input, timeout=timeout)
            if resultfile:
                _remap_intensity_column(resultfile, stats)
            return out

        proc.communicate = communicate
        return proc

    ts.subprocess.Popen = patched
    try:
        yield stats
    finally:
        ts.subprocess.Popen = real_popen


def _remap_intensity_column(resultfile: str, stats: dict) -> None:
    """bsyn intensity RESULTFILE: col0 wave, col1 normalized FLUX, col2/3 absolute,
    col4 normalized intensity at the requested mu. Rewrite as (wave, mu-intensity)
    so iSpec's data[:,1] reader picks up the intensity. No-op (loud) if the file is
    not in intensity format -> never silently substitute flux."""
    try:
        d = np.loadtxt(resultfile)        # np.loadtxt skips the '# mu-points' header
    except Exception:
        return
    if d.ndim != 2 or d.shape[1] < 5:
        raise RuntimeError(
            f"mu=1 intensity remap FAILED: RESULTFILE {resultfile} has shape "
            f"{getattr(d, 'shape', None)}, expected >=5 cols (wave, flux, 2x abs, "
            f"mu-intensity). bsyn did not emit intensity -> refusing to fit flux as mu=1.")
    wave, mu_int = d[:, 0], d[:, 4]
    np.savetxt(resultfile, np.column_stack([wave, mu_int]), fmt='%.4f %.6f')
    stats['remapped'] += 1


def run(mu: float = 1.0) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    tmp = '/tmp/ispec_ace_co_mu1_rya441'
    Path(tmp).mkdir(parents=True, exist_ok=True)
    params = ace._solar_params()

    print("=" * 78)
    print("RYA-441  ACE-FTS disk-center (mu=1) intensity CO fit -- the real go/no-go")
    print("=" * 78)
    print(f"  pinned solar params (constants.py): Teff={params['teff_K']:.0f} "
          f"logg={params['logg']:.3f} [Fe/H]={params['feh']:+.1f} xi={params['vturb_kms']:.2f}")
    print(f"  A(C)_ref={ace.A_C_REF:.2f}  A(O)_ref={ace.A_O_REF:.2f} (FIXED)  "
          f"flux baseline A(C)={A_C_FLUX_440:.3f} (RYA-440)")

    # Step 0 -- toolchain access
    print(f"\n[Step 0] toolchain: DIRECT mu={mu:.1f} specific-intensity via bsyn, by "
          f"flipping INTENSITY/FLUX->Intensity + injecting MU-POINTS (nangles=1, "
          f"mu={mu:.2f}) and remapping the intensity column. iSpec flux wrapper + "
          f"flux production path untouched.")

    # synthesis resources (same as 440)
    print("\n[setup] loading atmosphere + GES linelist + isotopes...")
    atm = _load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    ll, iso, chem = _load_synth_resources()
    sab = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    codes = _atom_codes(('C', 'O'), chem, sab)
    params['_atm'], params['_tmp'] = atm, tmp

    # locate the synthetic bandhead (vacuum frame) -> air/vac loud-fail (440 Step 1)
    synth_head = ace._fiducial_bandhead_A(params, ll, iso, sab, codes)
    s1 = ace.step1_load_ace(synth_head)
    print(f"[Step 1] air/vac: matched frame = {s1['matched_frame']}; other column "
          f"offset {s1['airvac_slip_kms']:.0f} km/s (air/vac slip) -> loud-fail OK")

    # fit A(C) in mu=1 intensity geometry (validate-don't-tune: A(O) fixed, params pinned)
    print(f"[Step 1'] fitting A(C) to ACE 12CO (2-0) bandhead {ace.FIT_WINDOW_VAC_A} A "
          f"in mu={mu:.1f} INTENSITY geometry...")
    state = {'C': ace.A_C_REF, 'O': ace.A_O_REF}
    with ace._StagedCO():
        with turbospectrum_intensity_at_mu(mu) as stats:
            fit = _fit_element(
                obs_w_nm=s1['wave_A'] / 10.0, obs_f=s1['flux'], atm=atm, params=params,
                free_el='C', state=state, codes=codes,
                windows_A=[ace.FIT_WINDOW_VAC_A], use_molecules=True,
                broadening=(ace.ACE_R, 0.0, 0.0),
                a_lo=ace.A_C_REF - 1.0, a_hi=ace.A_C_REF + 1.0,
                ll=ll, iso=iso, sab=sab, tmp_dir=tmp)
    a_mu1 = fit.get('A_X')
    print(f"        bsyn intensity injections fired={stats['fired']}, "
          f"intensity-column remaps={stats['remapped']}")
    if stats['remapped'] == 0:
        raise RuntimeError("mu=1 intensity remap never ran -> would be a silent flux "
                           "substitution. STOP.")
    print(f"        A(C)_mu1 = {a_mu1}  (red_chi2={fit.get('red_chi2')}, "
          f"n_pix={fit.get('n_pix')}, status={fit.get('status')})")

    # Step 2 -- sign tripwire
    delta_vs_flux = (a_mu1 - A_C_FLUX_440) if a_mu1 is not None and np.isfinite(a_mu1) else None
    sign_ok = delta_vs_flux is not None and delta_vs_flux >= -0.005   # allow tiny 1D float noise
    print(f"\n[Step 2] sign tripwire: A(C)_mu1 - A(C)_flux = {delta_vs_flux:+.3f} dex "
          f"(physics: mu=1 weaker CO -> MUST be >= 0)")
    if not sign_ok:
        raise RuntimeError(
            f"SIGN-WRONG: A(C)_mu1 ({a_mu1}) < A(C)_flux ({A_C_FLUX_440}). Contradicts "
            f"CO center-to-limb physics. STOP and RCA (mu handling / column remap / "
            f"continuum) -- not reporting a sign-wrong result as the verdict.")

    # Step 3 -- verdict
    delta_vs_ref = a_mu1 - ace.A_C_REF if a_mu1 is not None and np.isfinite(a_mu1) else None
    go = delta_vs_ref is not None and abs(delta_vs_ref) <= ace.NEAR_REF_DEX
    verdict = 'GO' if go else 'NO-GO'
    caveat = (
        "MANDATORY 1D-vs-3D caveat: this is a 1D MARCS mu=1-vs-flux bracket, a FLOOR "
        "on the real geometry sensitivity. CO center-to-limb is largely a "
        "3D/granulation effect that 1D underpredicts, so this small 1D bracket does "
        "NOT close the 3D uncertainty -- carried as a residual. The 3D-true mu=1 "
        "A(C) could sit further off-reference; a GO here is a 1D GO pending a "
        "3D-CO check, never a turnkey clearance.")

    print("\n" + "-" * 78)
    print(f"  A(C) mu=1 intensity  = {a_mu1}")
    print(f"  delta vs flux-8.497  = {delta_vs_flux:+.3f} dex  (sign check: "
          f"{'PASS' if sign_ok else 'FAIL'})")
    print(f"  delta vs Asplund ref = {delta_vs_ref:+.3f} dex  (band +/-{ace.NEAR_REF_DEX})")
    print(f"\n  ===> {verdict} (1D)")
    print("  " + caveat.replace('. ', '.\n  '))

    report = {'ticket': 'RYA-441', 'mu': mu, 'toolchain_path': 'direct mu=1 bsyn intensity',
              'params': {k: params[k] for k in ('teff_K', 'logg', 'feh', 'vturb_kms')},
              'A_C_reference': ace.A_C_REF, 'A_O_reference_fixed': ace.A_O_REF,
              'A_C_flux_440': A_C_FLUX_440,
              'A_C_mu1': a_mu1, 'fit': fit,
              'delta_vs_flux': delta_vs_flux, 'sign_check_pass': bool(sign_ok),
              'delta_vs_reference': delta_vs_ref, 'near_ref_band': ace.NEAR_REF_DEX,
              'verdict': verdict, 'verdict_is_1D': True, 'caveat_1d_vs_3d': caveat,
              'airvac': {k: s1[k] for k in ('matched_frame', 'airvac_slip_kms')}}
    (OUT / 'rya441_ace_co_mu1.json').write_text(json.dumps(report, indent=2))
    print(f"\n  [out] {OUT / 'rya441_ace_co_mu1.json'}")
    return report


def main():
    ap = argparse.ArgumentParser(description='RYA-441 ACE-FTS mu=1 intensity CO fit')
    ap.add_argument('--validate', action='store_true',
                    help='run the mu=1 fit + sign check + verdict')
    ap.add_argument('--mu', type=float, default=1.0, help='intensity angle (default 1.0)')
    args = ap.parse_args()
    run(mu=args.mu)


if __name__ == '__main__':
    main()

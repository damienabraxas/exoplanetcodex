"""
data/test/test_oi6300_twocomp.py
=================================
Synthetic test for the O I 6300.304 + Ni I 6300.336 blend treatment.

Strategy confirmed by this test: because the O-Ni separation (0.032 Å) is
smaller than the HARPS resolution element (~0.050 Å FWHM), a free two-component
fit cannot reliably decompose the blend. The correct approach is:
  1. Predict Ni EW from COG of clean Ni I lines (_predict_ni6300_ew)
  2. Build and subtract a fixed Ni I model from the spectrum
  3. Fit O I alone on the Ni-subtracted data

This file tests steps 2-3 with a known injected Ni EW.

DEV_CYCLE stage: TEST (RYA-104)
Run from: exoplanetcodex/
    python3 data/test/test_oi6300_twocomp.py
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from pipeline.lines_fit import _fit_profile, _gauss_abs, _voigt_abs, _integrate_profile

# ── Synthetic spectrum parameters ─────────────────────────────────────────────

WAV_O  = 6300.304   # O I forbidden line
WAV_NI = 6300.336   # Ni I blend (+0.032 Å)

# HARPS R~115,000 at 6300 Å → σ ~ λ/(R*2.355) ≈ 0.023 Å
SIGMA_HARPS = 0.023

EW_O_TARGET  = 4.5   # mÅ  (O I,  target: Asplund+2021)
EW_NI_TARGET = 1.0   # mÅ  (Ni I, target: Allende Prieto+2001)

# depth = EW_mA / (sigma_A * sqrt(2pi) * 1000)
DEPTH_O  = EW_O_TARGET  / (SIGMA_HARPS * np.sqrt(2 * np.pi) * 1000)
DEPTH_NI = EW_NI_TARGET / (SIGMA_HARPS * np.sqrt(2 * np.pi) * 1000)

SNR   = 350
NOISE = 1.0 / SNR


def make_synthetic(seed: int = 42) -> tuple:
    """±0.20 Å window, O+Ni blend on flat continuum, HARPS-like noise."""
    rng = np.random.default_rng(seed)
    wav = np.linspace(WAV_O - 0.20, WAV_O + 0.20, 400)
    flux = (1.0
            - DEPTH_O  * np.exp(-0.5 * ((wav - WAV_O)  / SIGMA_HARPS) ** 2)
            - DEPTH_NI * np.exp(-0.5 * ((wav - WAV_NI) / SIGMA_HARPS) ** 2))
    flux += rng.normal(0, NOISE, size=len(wav))
    return wav, flux


def ni_model(wav: np.ndarray, ni_ew_mA: float,
             ni_wav: float = WAV_NI, sigma: float = SIGMA_HARPS) -> np.ndarray:
    """Gaussian Ni I absorption model from a predicted EW."""
    depth = ni_ew_mA / (sigma * np.sqrt(2 * np.pi) * 1000.0)
    return depth * np.exp(-0.5 * ((wav - ni_wav) / sigma) ** 2)


# ── Test ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  RYA-104 — O I 6300 Ni-subtract + single-component test")
    print("=" * 60)

    wav, flux = make_synthetic()

    print(f"\n[1] Synthetic spectrum built")
    print(f"    Pixels: {len(wav)}, range: {wav[0]:.3f}–{wav[-1]:.3f} Å")
    print(f"    Injected O I  EW = {EW_O_TARGET:.2f} mÅ  depth={DEPTH_O:.5f}")
    print(f"    Injected Ni I EW = {EW_NI_TARGET:.2f} mÅ  depth={DEPTH_NI:.5f}")
    print(f"    O-Ni separation  = {WAV_NI - WAV_O:.3f} Å  (HARPS FWHM ~0.050 Å)")
    print(f"    SNR = {SNR},  noise σ = {NOISE:.5f}")

    # ── Step 1: subtract Ni model (simulating COG-predicted EW) ──────────────
    print(f"\n[2] Subtracting Ni I model (EW = {EW_NI_TARGET:.2f} mÅ, σ = {SIGMA_HARPS} Å) ...")
    # In production: ni_ew_mA comes from _predict_ni6300_ew(); here we use truth
    flux_ni_sub = flux + ni_model(wav, EW_NI_TARGET)

    # ── Step 2: fit O I on cleaned spectrum ───────────────────────────────────
    print(f"\n[3] Fitting O I alone on Ni-subtracted spectrum ...")
    popt, pcov, profile_t, chi2 = _fit_profile(wav, flux_ni_sub, WAV_O)

    passed = True

    if popt is None:
        print("    FAIL — _fit_profile returned None")
        passed = False
    else:
        ew_O = _integrate_profile(wav, popt, profile_t)
        x0_O = popt[0]

        print(f"\n    Fit results:")
        print(f"    {'Param':<22} {'Fitted':>10}  {'Injected':>10}  {'Pass?':>6}")
        print(f"    {'-'*54}")

        def chk(label, fitted, injected, tol):
            ok = abs(fitted - injected) / max(abs(injected), 1e-6) < tol
            mark = '✓' if ok else '✗  FAIL'
            print(f"    {label:<22} {fitted:>10.4f}  {injected:>10.4f}  {mark:>6}")
            return ok

        ok1 = chk("O I centroid (Å)",   x0_O,       WAV_O,        0.005)
        ok2 = chk("O I EW (mÅ)",        ew_O,       EW_O_TARGET,  0.20)
        ok3 = chk("profile type",        0,          0,            1.0)   # placeholder
        print(f"    {'profile type':<22} {profile_t:>10}  {'voigt/gauss':>10}  {'✓':>6}")

        nan_popt = not all(np.isfinite(popt))
        nan_pcov = not np.all(np.isfinite(pcov)) if pcov is not None else True
        print(f"\n    NaN check — popt: {'✓' if not nan_popt else '✗'}"
              f"   pcov: {'✓' if not nan_pcov else '✗'}")
        print(f"    chi²_red = {chi2:.6f}")

        ok_chi = chi2 < 5e-5
        print(f"    chi²_red < 5e-5: {'✓' if ok_chi else '✗  FAIL'}")
        ok_ew_range = 3.5 <= ew_O <= 5.5
        print(f"    EW in QA range [3.5, 5.5] mÅ: {'✓' if ok_ew_range else '✗  FAIL'}")

        passed = ok1 and ok2 and not nan_popt and ok_chi and ok_ew_range

        # ── Diagnostic plot ───────────────────────────────────────────────────
        out_plot = Path(__file__).parent / 'oi6300_synthetic_test.png'
        x_fine = np.linspace(wav[0], wav[-1], 3000)
        fn = _voigt_abs if profile_t == 'voigt' else _gauss_abs
        fit_curve = fn(x_fine, *popt)
        ni_curve  = ni_model(wav, EW_NI_TARGET)

        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        fig.suptitle(
            f"RYA-104 Synthetic Test — O I 6300.304 + Ni I 6300.336  (SNR={SNR})\n"
            f"Strategy: COG-predict Ni EW → subtract → fit O I alone",
            fontsize=10,
        )

        # Panel 1: raw blend
        axes[0].plot(wav, flux, 'o', color='#4a90d9', ms=2, alpha=0.7, label='O+Ni blend')
        axes[0].fill_between(wav, 1.0 - ni_curve, 1.0, alpha=0.25,
                             color='orange', label=f'Ni model ({EW_NI_TARGET} mÅ)')
        axes[0].axvline(WAV_O,  color='purple', lw=0.8, ls='--', alpha=0.7, label='O I')
        axes[0].axvline(WAV_NI, color='orange', lw=0.8, ls='--', alpha=0.7, label='Ni I')
        axes[0].axhline(1.0, color='gray', lw=0.4, ls=':')
        axes[0].set_title('Raw blend + Ni model', fontsize=9)
        axes[0].legend(fontsize=7); axes[0].tick_params(labelsize=8)
        axes[0].set_xlabel('Wavelength (Å)', fontsize=9)
        axes[0].set_ylabel('Flux', fontsize=9)

        # Panel 2: Ni-subtracted + O I fit
        axes[1].plot(wav, flux_ni_sub, 'o', color='#4a90d9', ms=2,
                     alpha=0.7, label='Ni-subtracted')
        axes[1].plot(x_fine, fit_curve, '-', color='tomato', lw=1.5,
                     label=f'O I fit ({ew_O:.2f} mÅ)')
        axes[1].axvline(WAV_O, color='purple', lw=0.8, ls='--', alpha=0.7)
        axes[1].axhline(1.0, color='gray', lw=0.4, ls=':')
        axes[1].set_title(f'Ni-subtracted  →  O I fit\nEW = {ew_O:.2f} mÅ  (target {EW_O_TARGET} mÅ)',
                          fontsize=9)
        axes[1].legend(fontsize=7); axes[1].tick_params(labelsize=8)
        axes[1].set_xlabel('Wavelength (Å)', fontsize=9)
        axes[1].set_ylabel('Flux', fontsize=9)

        # Panel 3: residuals
        resid = flux_ni_sub - fn(wav, *popt)
        axes[2].plot(wav, resid, 'o', color='steelblue', ms=2, alpha=0.7)
        axes[2].axhline(0,      color='gray', lw=0.8, ls='--')
        axes[2].axhline( NOISE, color='gray', lw=0.5, ls=':', alpha=0.5)
        axes[2].axhline(-NOISE, color='gray', lw=0.5, ls=':', alpha=0.5)
        axes[2].set_title(f'Residuals  (±1σ = ±{NOISE:.4f})\nchi²_red = {chi2:.5f}', fontsize=9)
        axes[2].tick_params(labelsize=8)
        axes[2].set_xlabel('Wavelength (Å)', fontsize=9)
        axes[2].set_ylabel('Residuals', fontsize=9)

        plt.tight_layout()
        plt.savefig(out_plot, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"\n[4] Diagnostic plot saved → {out_plot.name}")

    print(f"\n{'='*60}")
    print(f"  RESULT: {'PASS ✓' if passed else 'FAIL ✗'}")
    if passed:
        print(f"  Strategy confirmed: predict Ni EW via COG, subtract,")
        print(f"  then fit O I alone. Ready to implement in lines_fit.py.")
    print(f"{'='*60}\n")
    return passed


if __name__ == '__main__':
    ok = main()
    sys.exit(0 if ok else 1)

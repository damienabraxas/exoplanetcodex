"""Single-call SPECTRUM test for solar vs Procyon — RYA-263 diagnostic."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.constants import ISPEC_DIR
sys.path.insert(0, str(ISPEC_DIR))
import ispec
import numpy as np
import pandas as pd
from pipeline.abundances_derive import _load_atmosphere, _build_ispec_line_regions, _ew_to_abundance

SOLAR_ABUND_FILE = str(ISPEC_DIR / 'input' / 'abundances' / 'Asplund.2009' / 'stdatom.dat')


def run(label, teff, logg, feh, vmic, ew_csv, target_wave):
    ew_df = pd.read_csv(ew_csv)
    if 'blend_flag' in ew_df.columns:
        ew_df = ew_df[ew_df['blend_flag'] == False]
    ew_df = ew_df[(ew_df['ew_mA'] >= 5) & (ew_df['ew_mA'] <= 300)]

    row = ew_df[(ew_df['wavelength_air_A'] - target_wave).abs() < 0.1]
    if len(row) == 0:
        # fall back to any Fe I line
        row = ew_df[(ew_df['element'] == 'Fe') & (ew_df['ion'] == 'I')]
    row = row.iloc[[0]]

    lm = _build_ispec_line_regions(row)
    if len(lm) == 0:
        print(f"[{label}] No linelist match for {target_wave:.3f} Å"); return

    atm = _load_atmosphere(teff, logg, feh, vmic)
    solar_abund = ispec.read_solar_abundances(SOLAR_ABUND_FILE)

    print(f"\n[{label}]  Teff={teff} logg={logg} feh={feh} vmic={vmic}")
    print(f"  Line: {lm['wave_A'][0]:.4f} Å  EW={lm['ew'][0]:.2f} mÅ  note={lm['note'][0]}")

    s, n, xh, xfe = ispec.determine_abundances(
        atm, teff, logg, feh, 0.0, lm, solar_abund,
        microturbulence_vel=vmic, verbose=0, code='spectrum',
        tmp_dir='/tmp/ispec_codex',
    )
    print(f"  spec_abund   = {s}  (SPECTRUM native, log N/N_H≈0)")
    print(f"  normal_abund = {n}  (iSpec 'normal scale', H=12)")
    print(f"  x_over_h     = {xh}  ([X/H])")
    print(f"  x_over_fe    = {xfe}  ([X/Fe])")

    # Also run full _ew_to_abundance on ALL elements (as the pipeline does)
    params = {'teff_K': teff, 'logg': logg, 'feh': feh, 'vturb_kms': vmic}
    lm_all, norm_all, xh_all, _ = _ew_to_abundance(ew_df, params, atm)
    notes = np.array([str(x) for x in lm_all['note']])
    fe1 = np.where(notes == 'Fe 1')[0]
    print(f"\n  Full _ew_to_abundance ({len(fe1)} Fe I lines):")
    print(f"  normal_abund Fe I: min={norm_all[fe1].min():.4f}  max={norm_all[fe1].max():.4f}  mean={norm_all[fe1].mean():.4f}")
    print(f"  x_over_h Fe I:     min={xh_all[fe1].min():.4f}  max={xh_all[fe1].max():.4f}  mean={xh_all[fe1].mean():.4f}")
    n_zero = int((norm_all[fe1] == 0.0).sum())
    print(f"  n_zero in normal_abund: {n_zero}/{len(fe1)}")


if __name__ == '__main__':
    print("=== PROCYON (post-fix: Fe II 5376.466 Å flagged) ===")
    run('procyon', 6530, 3.96, -0.04, 1.66,
        'data/processed/procyon_ew.csv', 5379.574)

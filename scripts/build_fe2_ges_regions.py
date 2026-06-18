"""
scripts/build_fe2_ges_regions.py
=================================
Vets orphaned Fe II lines from solar_ew_results and appends passing lines
as GES region entries to the Codex-extended line regions file.

Vetting criteria for Fe II (science decision RYA-228 follow-up, June 2026):
  - EW in [10, 120] mÅ  (COG reliability window)
  - ew_err < ew_mA       (fit must have converged)
  - vald_proximity_flag criterion: DROPPED for Fe II
      Rationale: existing GES Fe II regions contain prox up to 0.635.
      EW and error gates are the meaningful quality filters.
      Blend contamination is handled by the COG ceiling.
  - Minimum passing lines: 5 (ionisation balance requires adequate Fe II coverage)

Output: appends to turbospectrum_synth_good_for_params_codex_extended.txt
        (created from the all_extended base file if it does not exist)
"""

import math
import shutil
import pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config.constants import PATHS, ISPEC_DIR
from pipeline.species import species_key   # RYA-345 canonical species matcher

_FE_II = (26, 2)


def _is_fe2(element, ion=None) -> bool:
    """True iff (element, ion) / a region note / a species code is Fe II,
    routed through the canonical normalizer (RYA-345) — encoding-agnostic."""
    try:
        return species_key(element, ion) == _FE_II
    except ValueError:
        return False

# ── Paths (derived from PATHS / ISPEC_DIR in constants.py) ──────────────────
EW_FILE  = PATHS['solar_ew']
LINELIST = PATHS['linelist_solar']

_GES_BASE   = (ISPEC_DIR / 'input' / 'regions' / '42000_VALD' /
               'turbospectrum_synth_good_for_params_all_extended.txt')
GES_FILE    = (ISPEC_DIR / 'input' / 'regions' / '42000_VALD' /
               'turbospectrum_synth_good_for_params_codex_extended.txt')

MATCH_TOL_A = 0.15   # Å — matching tolerance between EW file and GES / linelist
EW_MIN      = 10.0   # mÅ
EW_MAX      = 120.0  # mÅ
MIN_PASSING = 5      # CRITICAL STOP threshold

FIT_HALF_WINDOW_A = 0.40   # ± Å around line centre for fitting window

# eV → cm⁻¹ conversion (exact: 1 eV = 8065.544 cm⁻¹)
EV_TO_CM1 = 8065.544
# hc in eV·Å (used to compute upper level from lower level + photon energy)
HC_EV_A   = 12398.42

# ── Create codex_extended GES file if it doesn't exist ──────────────────────
if not GES_FILE.exists():
    shutil.copy2(_GES_BASE, GES_FILE)
    print(f"Created {GES_FILE.name} from {_GES_BASE.name}")
else:
    print(f"Using existing {GES_FILE.name}")

# ── Load data ─────────────────────────────────────────────────────────────────
ew_df = pd.read_csv(EW_FILE)
ll_df = pd.read_csv(LINELIST)

# Fe II only (RYA-345: canonical species match, encoding-agnostic)
fe2_ew = ew_df[[_is_fe2(e, i) for e, i in zip(ew_df['element'], ew_df['ion'])]].copy()
fe2_ll = ll_df[[_is_fe2(e, i) for e, i in zip(ll_df['element'], ll_df['ion'])]].copy()

# ── Identify already-matched lines (present in existing GES file) ─────────────
ges_wavs = []
with open(GES_FILE) as f:
    next(f)  # skip header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2 and _is_fe2(parts[0]):  # element col identifies Fe II
            try:
                ges_wavs.append(float(parts[1]))  # wave_A col
            except ValueError:
                pass

ges_wav_series = pd.Series(ges_wavs)


def nearest_match(wav, candidates, tol):
    if candidates.empty:
        return None
    diffs = (candidates - wav).abs()
    idx = diffs.idxmin()
    return idx if diffs[idx] <= tol else None


already_matched = []
orphaned = []
for _, row in fe2_ew.iterrows():
    idx = nearest_match(row['wavelength_air_A'], ges_wav_series, MATCH_TOL_A)
    if idx is not None:
        already_matched.append(row)
    else:
        orphaned.append(row)

orphaned_df = pd.DataFrame(orphaned)
print(f"\nFe II EW lines: {len(fe2_ew)}  |  Already in GES: {len(already_matched)}  |  Orphaned: {len(orphaned_df)}")

# ── Vet orphaned lines ────────────────────────────────────────────────────────
# NOTE: vald_proximity_flag criterion intentionally omitted for Fe II.
# See module docstring for science rationale.

results = []
for _, row in orphaned_df.iterrows():
    ew      = row['ew_mA']
    ew_err  = row['ew_err_mA']
    wav     = row['wavelength_air_A']

    reasons = []
    if not (EW_MIN <= ew <= EW_MAX):
        reasons.append(f'EW={ew:.1f} outside [{EW_MIN:.0f},{EW_MAX:.0f}] mA')
    if ew_err >= ew:
        reasons.append(f'ew_err={ew_err:.2f} >= EW={ew:.2f}')

    status = 'PASS' if not reasons else 'FAIL'
    results.append({'wavelength_air_A': wav, 'ew_mA': ew, 'ew_err_mA': ew_err,
                    'status': status, 'reasons': '; '.join(reasons)})

results_df = pd.DataFrame(results)
passing = results_df[results_df['status'] == 'PASS']
failing = results_df[results_df['status'] == 'FAIL']

print(f"\nVetting results (prox criterion dropped for Fe II):")
print(results_df[['wavelength_air_A', 'ew_mA', 'ew_err_mA', 'status', 'reasons']].to_string(index=False))
print(f"\nPASS: {len(passing)}  FAIL: {len(failing)}")

if len(passing) < MIN_PASSING:
    print(f"\nCRITICAL STOP: {len(passing)} passing lines < threshold of {MIN_PASSING}.")
    print("Do not append to GES file. Post results to RYA-243 and stop.")
    raise SystemExit(1)

# ── Match passing lines to linelist for atomic data ───────────────────────────
new_entries = []
for _, row in passing.iterrows():
    wav = row['wavelength_air_A']
    diffs = (fe2_ll['wavelength_air_A'] - wav).abs()
    best_idx = diffs.idxmin()
    if diffs[best_idx] > MATCH_TOL_A:
        print(f"WARNING: No linelist match within {MATCH_TOL_A} A for Fe II {wav:.5f} A — skipping")
        continue
    ll_row = fe2_ll.loc[best_idx]
    new_entries.append({
        'wavelength_air_A'        : wav,
        'ew_mA'                   : row['ew_mA'],
        'ew_err_mA'               : row['ew_err_mA'],
        'excitation_potential_eV' : ll_row['excitation_potential_eV'],
        'log_gf'                  : ll_row['log_gf'],
        'damping_rad'             : float(ll_row['damping_rad']),
        'damping_stark'           : float(ll_row['damping_stark']),
        'damping_vdW'             : float(ll_row['damping_vdW']),
        'central_depth'           : float(ll_row.get('central_depth', 0.1)),
        'loggf_source'            : str(ll_row.get('loggf_source', 'VALD3')),
    })

print(f"\n{len(new_entries)} lines ready to append to GES regions file.")

# ── Print GES file format for verification ────────────────────────────────────
print("\nExisting GES file format (first 3 lines):")
with open(GES_FILE) as f:
    for _ in range(3):
        l = f.readline()
        print(repr(l))

# ── Count Fe II entries before append ─────────────────────────────────────────
fe2_before = len(ges_wavs)
print(f"\nFe II entries in GES before append: {fe2_before}")

# ── Build GES-format rows for new entries ─────────────────────────────────────

def _fmt_sci(x):
    """Format a float in scientific notation matching GES style (e.g. 3.47E+08)."""
    return f"{x:.2E}"


def _build_ges_row(entry):
    """
    Construct an 83-column tab-separated GES row for a new Fe II entry.

    Atomic columns filled from linelist_solar.csv.  Observational columns
    (fitting window, EW, profile parameters) filled from our solar_ew.csv
    measurements or sensible defaults.  Quantum numbers (j, Lande g) not
    available in linelist_solar — set to 0.0.
    """
    wav_A  = entry['wavelength_air_A']
    wav_nm = wav_A / 10.0

    lo_eV  = entry['excitation_potential_eV']
    lo_cm1 = lo_eV * EV_TO_CM1
    hi_eV  = lo_eV + HC_EV_A / wav_A
    hi_cm1 = hi_eV * EV_TO_CM1

    loggf   = entry['log_gf']
    d_rad   = entry['damping_rad']
    d_stark = entry['damping_stark']
    d_vdW   = entry['damping_vdW']
    depth   = entry['central_depth']
    ew      = entry['ew_mA']
    ew_err  = entry['ew_err_mA']

    turbo_rad = _fmt_sci(10.0 ** d_rad)

    # Fitting window
    base_nm = (wav_A - FIT_HALF_WINDOW_A) / 10.0
    top_nm  = (wav_A + FIT_HALF_WINDOW_A) / 10.0

    # Profile parameters (Gaussian approximation for solar Fe II)
    # FWHM at HARPS-like R~115000 plus thermal+macro broadening → ~0.020 nm
    fwhm_nm  = 0.020
    sig_nm   = fwhm_nm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    fwhm_kms = fwhm_nm * 299792.458 / wav_nm
    R_val    = round(wav_nm / fwhm_nm)

    # EWR = log10(EW_Å / lambda_Å)
    ew_A  = ew / 1000.0
    ewr   = math.log10(ew_A / wav_A) if ew_A > 0 else -6.0
    integ = ew_A / 10.0   # EW in nm

    snr = ew / ew_err if ew_err > 0 else 0.0

    ref_code = f"VALD3  E 0"

    cols = [
        # 1–13: element and quantum numbers
        'Fe 2',         # element
        f'{wav_A:.3f}', # wave_A
        f'{wav_nm:.4f}',# wave_nm
        f'{loggf:.3f}', # loggf
        f'{lo_eV:.4f}', # lower_state_eV
        f'{lo_cm1:.3f}',# lower_state_cm1
        '0.0',          # lower_j
        f'{hi_eV:.4f}', # upper_state_eV
        f'{hi_cm1:.3f}',# upper_state_cm1
        '0.0',          # upper_j
        '0.0',          # upper_g
        '0.0',          # lande_lower
        '0.0',          # lande_upper
        # 14: transition type
        'AO',
        # 15–20: damping
        turbo_rad,       # turbospectrum_rad
        f'{d_rad:.3f}',  # rad
        f'{d_stark:.3f}',# stark
        f'{d_vdW:.3f}',  # waals
        f'{d_vdW:.3f}',  # waals_single_gamma_format
        f'{d_vdW:.3f}',  # turbospectrum_fdamp
        # 21–23: fudge, theoretical
        '1.000',
        f'{depth:.3f}',  # theoretical_depth
        '0.000',         # theoretical_ew
        # 24–27: orbital types, molecule, isotope
        's', 'p', 'F', '0',
        # 28–31: ion and species codes
        '2',
        '26.1',   # spectrum_moog_species
        '26.1',   # turbospectrum_species
        '26.01',  # width_species
        # 32: reference
        ref_code,
        # 33–37: NLTE
        'T', '0', '0', 'none', 'none',
        # 38–43: support flags
        'T', 'T', 'T', 'T', 'T', 'T',
        # 44–46: fitting window
        f'{wav_nm:.4f}',  # wave_peak
        f'{base_nm:.4f}', # wave_base
        f'{top_nm:.4f}',  # wave_top
        # 47–50: note and pixel indices
        'Fe 2', '0', '0', '0',
        # 51–52: depth
        f'{depth:.3f}',  # depth
        f'{depth:.3f}',  # relative_depth
        # 53–56: fit window + pixel indices (repeat)
        f'{base_nm:.4f}', f'{top_nm:.4f}', '0', '0',
        # 57–65: profile parameters
        f'{wav_nm:.4f}',   # mu
        f'{sig_nm:.4f}',   # sig
        f'{-depth:.4f}',   # A
        '1.000',           # baseline
        '9999.0000',       # gamma
        f'{sig_nm:.4f}',   # mu_err (reuse sig as conservative estimate)
        f'{fwhm_nm:.4f}',  # fwhm
        f'{fwhm_kms:.4f}', # fwhm_kms
        f'{R_val}',        # R
        # 66–68: fit depth + flux
        f'{depth:.3f}',    # depth_fit
        f'{depth:.3f}',    # relative_depth_fit
        f'{integ:.4f}',    # integrated_flux
        # 69–76: EW and flux stats
        f'{ewr:.2f}',      # ewr
        f'{ew:.1f}',       # ew
        f'{ew_err:.1f}',   # ew_err
        f'{snr:.1f}',      # snr
        '1.000',           # mean_flux
        '1.000',           # mean_flux_continuum
        '0.0000',          # diff_wavelength
        '0.005',           # rms
        # 77–80: telluric
        '0.0000', '0.0000', '0', '0.000',
        # 81–83: grouped, reference_for_group, discarded
        'False', 'False', 'False',
    ]

    assert len(cols) == 83, f"Expected 83 columns, got {len(cols)}"
    return '\t'.join(str(c) for c in cols)


# ── Append to GES file ────────────────────────────────────────────────────────
with open(GES_FILE, 'a') as f:
    for entry in new_entries:
        row = _build_ges_row(entry)
        f.write(row + '\n')

# ── Count Fe II entries after append ─────────────────────────────────────────
fe2_after = 0
with open(GES_FILE) as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2 and parts[0].strip() == 'Fe 2':
            fe2_after += 1

print(f"\nFe II entries in GES after append:  {fe2_after}  (+{fe2_after - fe2_before} new)")
print("\nDone. Post results to Linear issue RYA-243 per end-of-session requirements.")

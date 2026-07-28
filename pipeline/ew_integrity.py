"""
pipeline/ew_integrity.py
========================
RYA-458 — the EW-VERIFICATION layer (EW half of the RYA-451 epic; the continuum
half RYA-452..454 concluded [O I] is continuum-limited).

WHAT THIS IS
------------
A per-line EW-INTEGRITY QA pass over the solar run. It FLAGS measurement anomalies
and assigns dispositions; it NEVER adjusts a measured EW and NEVER shrinks a value
toward a literature EW. The literature reference is a CROSS-CHECK, not a tuning
target. `assert_no_ew_mutation` proves the EWs are untouched.

THE CARDINAL RULE — validate, don't tune (RYA-451/395)
------------------------------------------------------
Every threshold below is fixed on PHYSICAL / fit-quality grounds, set BEFORE any
abundance is read, documented, and applied UNIFORMLY across all elements. None is a
function of A(X) or of the literature EW it is checked against. The layer's only
outputs are flags, dispositions, and exclusions with named reasons.

TWO MECHANISMS (internal first, literature second)
--------------------------------------------------
(A) Abundance-internal per-line consistency (universal, needs no external table):
      BAD_FIT       — fit reduced chi-square above a fixed ceiling (EW profile fit for
                      EW-path lines; cited synthesis flux-fit chi2r for the C/N/O anchors,
                      e.g. C I 5380's ESPRESSO chi2r~103).
      ABUND_OUTLIER — the line's implied A(X) sits >N robust-sigma (median/MAD) from its
                      element-ion's robust mean (internal self-consistency; fixed N).
      COG_FLAG      — reduced EW above the linear curve-of-growth knee (saturated), or
                      EW above the absolute saturation ceiling.
(B) Literature-EW cross-check where a CITED reference EW exists (seed table
    data/curation/ew_reference/solar_ew_reference.csv): flag LIT_DEVIATION beyond a
    fixed band. A missing reference is silent, not a failure.

THE THREE CHARTER CASES (the reason this exists)
------------------------------------------------
  C I 5380  -> BAD_FIT (cited chi2r~103 + RYA-454 2.8% continuum strike; EW-pool
               profile 149.5 mA is saturated) -> excluded from the C verdict.
  Li 6707   -> CN-blend contaminated -> UPPER_LIMIT disposition, never a point value.
  Eu 6645   -> after the RYA-102 HFS-summing fix, RECOVERED if the EW lands in the
               documented ~6-10 mA window, else FITTER_INCOMPLETE.

Linear issue: RYA-458
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from config.constants import LINE_SCORE_PARAMS, PIPELINE  # noqa: E402

# ─────────────────────────────────────────────────────────────────────────────
# FIXED, DOCUMENTED EW-INTEGRITY THRESHOLDS
# Set on fit-quality / curve-of-growth grounds, BEFORE any abundance is read,
# uniform across elements. Reuse the pipeline-wide constants where the physics is
# already encoded (single source) rather than inventing a parallel number.
# ─────────────────────────────────────────────────────────────────────────────
# BAD_FIT — EW PROFILE fit reduced chi-square ceiling. Reuse the RYA-220 line-score
# "chi2 scores 0" floor: a profile fit worse than this is unusable.
EW_FIT_CHI2_MAX = float(LINE_SCORE_PARAMS['chi2_floor'])        # 10.0
# BAD_FIT — SYNTHESIS flux-fit reduced chi-square ceiling for the C/N/O anchors. The
# CNO flux fits run over 150-1000 px and are not per-DOF-normalised, so the ceiling is
# set high enough to flag only a genuinely broken fit (C I 5380 ESPRESSO chi2r~103),
# not the merely-large fits of valid anchors (CI_5052 ~56, OI_6300 ~66).
SYNTH_FIT_CHI2R_MAX = 100.0
# ABUND_OUTLIER — robust (median/MAD) per-line deviation, in robust sigma. Fixed N.
ABUND_OUTLIER_NSIGMA = 3.0
ABUND_OUTLIER_MIN_LINES = 4                                     # need a real distribution
_MAD_TO_SIGMA = 1.4826
_ABUND_MAD_FLOOR = 0.03                                         # dex — don't divide by ~0 scatter
# COG_FLAG — linear curve-of-growth knee (Gray, Stellar Photospheres) + absolute
# saturation ceiling. Both reused from the curation / pipeline single source.
REW_LINEAR_CEILING = -4.90
SAT_EW_CEILING_MA = float(PIPELINE['vmic_ew_ceiling_mA'])       # 100 mA
# LIT_DEVIATION — fractional band on |measured - reference| / reference for a cited
# point-value reference. Wide on purpose: a cross-check that catches a gross mismatch,
# never a tuning target.
LIT_DEVIATION_FRAC = 0.50
# Eu HFS recovery window (RYA-102) — the disposition gate for the Eu charter case.
EU_RECOVER_LO_MA = 6.0
EU_RECOVER_HI_MA = 10.0

_REF_PATH = _REPO / 'data' / 'curation' / 'ew_reference' / 'solar_ew_reference.csv'
_REF_WAV_TOL = 0.05                                            # A — match a reference to a line

# The named flags / dispositions (closed vocabulary — every flag is one of these).
FLAGS = ('BAD_FIT', 'ABUND_OUTLIER', 'COG_FLAG', 'LIT_DEVIATION')
DISPOSITIONS = ('UPPER_LIMIT', 'BAD_FIT', 'RECOVERED', 'FITTER_INCOMPLETE')
# Flags that EXCLUDE a line from its element verdict (quality strikes).
EXCLUDING_FLAGS = frozenset({'BAD_FIT', 'ABUND_OUTLIER', 'COG_FLAG'})


def load_reference_table(path=None) -> pd.DataFrame:
    """The cited cross-check reference table (provenance header skipped)."""
    p = Path(path) if path is not None else _REF_PATH
    df = pd.read_csv(p, comment='#')
    df['wavelength_air_A'] = df['wavelength_air_A'].astype(float)
    return df


def _match_reference(ref: pd.DataFrame, element, ion, wav) -> pd.Series | None:
    sub = ref[(ref['element'] == element) & (ref['ion'] == ion)]
    if sub.empty:
        return None
    d = (sub['wavelength_air_A'] - float(wav)).abs()
    j = int(d.values.argmin())
    if float(d.iloc[j]) > _REF_WAV_TOL:
        return None
    return sub.iloc[j]


def _rew(ew_mA, wav) -> float:
    if not (np.isfinite(ew_mA) and ew_mA > 0 and np.isfinite(wav) and wav > 0):
        return np.nan
    return float(np.log10(ew_mA / 1000.0 / wav))


def _abund_outlier_mask(df: pd.DataFrame) -> np.ndarray:
    """Per (element, ion) robust (median/MAD) outlier mask on the implied A(X).
    Internal self-consistency only — compares each line to its OWN species' robust
    mean, never to an external anchor."""
    out = np.zeros(len(df), dtype=bool)
    if 'a_lte' not in df.columns:
        return out
    a = df['a_lte'].to_numpy(dtype=float)
    for (_el, _ion), idx in df.groupby(['element', 'ion']).groups.items():
        ii = np.asarray(list(idx))
        vals = a[ii]
        fin = np.isfinite(vals)
        if fin.sum() < ABUND_OUTLIER_MIN_LINES:
            continue
        med = np.median(vals[fin])
        mad = np.median(np.abs(vals[fin] - med)) * _MAD_TO_SIGMA
        scale = max(mad, _ABUND_MAD_FLOOR)
        dev = np.abs(vals - med) / scale
        out[ii] = np.where(np.isfinite(dev), dev > ABUND_OUTLIER_NSIGMA, False)
    return out


def flag_ew_integrity(per_line_df: pd.DataFrame, reference=None) -> pd.DataFrame:
    """Add the EW-integrity flags + disposition to a per-line frame. Returns a COPY
    with new columns; the input EWs are never modified.

    Expected columns: element, ion, wavelength_air_A, ew_mA. Optional and used when
    present: ew_err_mA, chi2 (EW profile fit), synth_chi2r (synthesis flux fit),
    a_lte (implied A(X) for ABUND_OUTLIER), blend_flag, notes.

    New columns:
      ew_integrity   — comma-joined flags from FLAGS (empty string = clean)
      ew_disposition — one of DISPOSITIONS, or '' (none)
      ew_excluded    — bool: carries an EXCLUDING_FLAG (BAD_FIT/ABUND_OUTLIER/COG_FLAG)
      ew_reference_mA, ew_lit_delta_mA, ew_reason — cross-check provenance
    """
    df = per_line_df.copy().reset_index(drop=True)
    n = len(df)
    if n == 0:
        for c in ('ew_integrity', 'ew_disposition', 'ew_excluded',
                  'ew_reference_mA', 'ew_lit_delta_mA', 'ew_reason'):
            df[c] = pd.Series(dtype='object')
        return df
    ref = reference if reference is not None else load_reference_table()

    outlier = _abund_outlier_mask(df)

    flags_col, disp_col, refmA_col, litd_col, reason_col = [], [], [], [], []
    for i, row in df.iterrows():
        el, ion = row['element'], row['ion']
        wav = float(row['wavelength_air_A'])
        ew = float(row['ew_mA']) if np.isfinite(row.get('ew_mA', np.nan)) else np.nan
        flags, reasons = [], []
        disposition = ''
        ref_mA, lit_delta = np.nan, np.nan

        # ── (A) internal per-line consistency ────────────────────────────────
        chi2 = float(row.get('chi2', np.nan)) if pd.notna(row.get('chi2', np.nan)) else np.nan
        synth_chi2r = (float(row.get('synth_chi2r', np.nan))
                       if pd.notna(row.get('synth_chi2r', np.nan)) else np.nan)
        if np.isfinite(chi2) and chi2 > EW_FIT_CHI2_MAX:
            flags.append('BAD_FIT'); reasons.append(f'profile chi2 {chi2:.2f} > {EW_FIT_CHI2_MAX:.0f}')
        if np.isfinite(synth_chi2r) and synth_chi2r > SYNTH_FIT_CHI2R_MAX:
            flags.append('BAD_FIT'); reasons.append(f'synth chi2r {synth_chi2r:.0f} > {SYNTH_FIT_CHI2R_MAX:.0f}')

        rew = _rew(ew, wav)
        if (np.isfinite(rew) and rew > REW_LINEAR_CEILING) or (np.isfinite(ew) and ew > SAT_EW_CEILING_MA):
            flags.append('COG_FLAG')
            reasons.append(f'REW {rew:.2f} > {REW_LINEAR_CEILING} / EW {ew:.1f} > {SAT_EW_CEILING_MA:.0f} (saturated)')

        if bool(outlier[i]):
            flags.append('ABUND_OUTLIER'); reasons.append(f'>{ABUND_OUTLIER_NSIGMA:.0f} robust-sigma in A(X)')

        # ── (B) cited literature cross-check + charter disposition ───────────
        r = _match_reference(ref, el, ion, wav)
        if r is not None:
            ref_mA = float(r['reference_ew_mA']) if pd.notna(r.get('reference_ew_mA', np.nan)) else np.nan
            kind = str(r.get('ew_kind', '') or '')
            # cited synthesis fit anomaly (e.g. C I 5380 chi2r=103) -> BAD_FIT
            cited_chi2r = float(r['fit_chi2r']) if pd.notna(r.get('fit_chi2r', np.nan)) else np.nan
            if np.isfinite(cited_chi2r) and cited_chi2r > SYNTH_FIT_CHI2R_MAX and 'BAD_FIT' not in flags:
                flags.append('BAD_FIT'); reasons.append(f'cited synth chi2r {cited_chi2r:.0f} (RYA-371/454)')
            # point-value cross-check (never for an upper limit)
            if np.isfinite(ref_mA) and np.isfinite(ew) and kind == 'value' and ref_mA > 0:
                lit_delta = ew - ref_mA
                if abs(lit_delta) / ref_mA > LIT_DEVIATION_FRAC:
                    flags.append('LIT_DEVIATION')
                    reasons.append(f'measured {ew:.1f} vs cited {ref_mA:.1f} mA '
                                   f'({lit_delta:+.1f}, > {LIT_DEVIATION_FRAC:.0%})')
            # fixed charter disposition from the table
            tdisp = str(r.get('disposition', '') or '')
            if tdisp in DISPOSITIONS:
                disposition = tdisp
            # Eu HFS recovery (RYA-102): computed disposition for the charter line
            if bool(r.get('charter_case', False)) and el == 'Eu' and not disposition:
                if np.isfinite(ew) and EU_RECOVER_LO_MA <= ew <= EU_RECOVER_HI_MA:
                    disposition = 'RECOVERED'
                    reasons.append(f'HFS-summed EW {ew:.1f} mA in [{EU_RECOVER_LO_MA:.0f},'
                                   f'{EU_RECOVER_HI_MA:.0f}] (RYA-102 recovered)')
                else:
                    disposition = 'FITTER_INCOMPLETE'
                    reasons.append(f'HFS EW {ew:.1f} mA outside [{EU_RECOVER_LO_MA:.0f},'
                                   f'{EU_RECOVER_HI_MA:.0f}] — HFS not summed across components')

        flags_col.append(','.join(dict.fromkeys(flags)))     # de-dup, preserve order
        disp_col.append(disposition)
        refmA_col.append(round(ref_mA, 2) if np.isfinite(ref_mA) else np.nan)
        litd_col.append(round(lit_delta, 2) if np.isfinite(lit_delta) else np.nan)
        reason_col.append('; '.join(reasons))

    df['ew_integrity'] = flags_col
    df['ew_disposition'] = disp_col
    df['ew_excluded'] = [bool(set(f.split(',')) & EXCLUDING_FLAGS) if f else False for f in flags_col]
    df['ew_reference_mA'] = refmA_col
    df['ew_lit_delta_mA'] = litd_col
    df['ew_reason'] = reason_col

    assert_no_ew_mutation(per_line_df, df)
    return df


def assert_no_ew_mutation(before: pd.DataFrame, after: pd.DataFrame) -> None:
    """The cardinal guard: the layer flags, it never edits an EW. Raise if any
    measured ew_mA changed between input and output (row-aligned)."""
    b = before.reset_index(drop=True)['ew_mA'].to_numpy(dtype=float)
    a = after.reset_index(drop=True)['ew_mA'].to_numpy(dtype=float)
    if len(b) != len(a):
        raise AssertionError(f"EW-integrity layer changed the row count "
                             f"({len(b)} -> {len(a)}) — it must only ADD flag columns.")
    bad = ~((np.isnan(b) & np.isnan(a)) | (b == a))
    if bad.any():
        k = int(bad.argmax())
        raise AssertionError(
            f"EW-integrity layer MUTATED a measured EW (row {k}: {b[k]} -> {a[k]}). "
            f"The layer must flag/exclude only, never adjust an EW (RYA-451/458 cardinal rule).")


def charter_summary(flagged: pd.DataFrame) -> dict:
    """Pull the three charter-case dispositions out of a flagged frame for reporting."""
    out = {}
    cases = [('C', 'I', 5380.337, 'C_I_5380'),
             ('Li', 'I', 6707.840, 'Li_6707'),
             ('Eu', 'II', 6645.127, 'Eu_6645')]
    for el, ion, wav, key in cases:
        sub = flagged[(flagged['element'] == el) & (flagged['ion'] == ion)]
        if sub.empty:
            out[key] = {'present': False}
            continue
        d = (sub['wavelength_air_A'] - wav).abs()
        r = sub.iloc[int(d.values.argmin())]
        out[key] = {
            'present': True,
            'wavelength_air_A': float(r['wavelength_air_A']),
            'ew_mA': float(r['ew_mA']),
            'ew_integrity': str(r['ew_integrity']),
            'ew_disposition': str(r['ew_disposition']),
            'ew_excluded': bool(r['ew_excluded']),
            'ew_reason': str(r['ew_reason']),
        }
    return out

"""
pipeline/rejection_ledger.py
============================
Per-line rejection ledger for the EW fitter (RYA-429).

Guarantee (the "no-silent-drop" invariant): no line in the input linelist
exits the fitter unaccounted-for. Every in-coverage line that does NOT
produce a measured EW is recorded here with an explicit, controlled-vocabulary
reason and the numeric value that tripped the gate.

    count(in-coverage input lines) == count(measured) + count(rejected-with-reason)

Motivated by the Sr II silent-drop RCA: all four canonical Sr II lines
(4077 / 4161 / 4215 / 4305) sat inside coverage yet none reached the EW pool,
and the two MODERATE subordinate lines (4161 / 4305) vanished with no record.
Same failure family as RYA-329 (mq_grade computed-but-never-applied leak).

Scope of "input linelist": the priority>0 set that lines_fit.run() constructs
and literally calls "the line list" (05.._runtime prints "N priority>0 lines").
Lines with priority<=0 are graded OUT upstream at linelist-build time
(RYA-354) — a documented pre-fitter gate, not a silent fitter drop — and are
outside this guarantee. Accounting granularity is the fitter's own dedup unit:
(element, ion, wavelength rounded to 0.1 A), because two same-species lines
within 0.1 A are a single unresolved feature at HARPS resolution.

Related: RYA-430 (Sr II line-selection science — owns WHY 4161/4305 are hard),
RYA-422 (Sr I vs Sr II ionization-balance science diagnostic). This module is
INFRASTRUCTURE only: it records reasons and raises the stage-presence audit
flag; it does not adjudicate the Sr science.

Linear issue: RYA-429
"""

import warnings

import numpy as np
import pandas as pd


# ── Controlled rejection-reason vocabulary ────────────────────────────────────
# A rejected line MUST carry exactly one of these. 'other' additionally REQUIRES
# a non-empty freetext note (a blank/placeholder reason defeats the ledger).
REJECTION_REASONS = (
    'below_priority_threshold',  # priority<=0, or central_depth<min_fit_depth: pre-fit selection gate
    'blend_overlap',             # blend_flag=True (not force-measured), or modeled/subtracted as a blend partner
    'failed_chi2',               # profile fit did not converge / produced a non-finite EW
    'saturated_core',            # measured EW above the max reliable EW (strong / saturated line)
    'low_snr',                   # measured EW below the min reliable EW floor
    'continuum_fail',            # local continuum re-normalisation could not be established
    'out_of_coverage',           # line outside the spectral range, or too few in-window pixels (gap)
    'other',                     # anything else — REQUIRES a freetext note
)

_UNCLASSIFIED_NOTE = 'UNCLASSIFIED'  # marker in note for a genuine silent drop the instrumentation missed


def line_key(element, ion, wavelength_air_A) -> str:
    """Accounting unit: (element, ion, wavelength rounded to 0.1 A). Uses numpy
    rounding to match _build_worklist's `.round(1)` dedup EXACTLY — printf-style
    `%.1f` and numpy round disagree on exact-half values (e.g. 5039.95 -> 5039.9
    vs 5040.0), which would split one unresolved feature into two buckets."""
    return f"{np.round(float(wavelength_air_A), 1)}|{str(element)}|{str(ion)}"


def _clean_value(value):
    if value is None:
        return np.nan
    try:
        v = float(value)
    except (TypeError, ValueError):
        return np.nan
    return v if np.isfinite(v) else np.nan


class RejectionLedger:
    """Accumulates per-line rejections for one star. First reason logged for a
    given line_key wins (the earliest gate that tripped)."""

    def __init__(self, star_key: str):
        self.star_key = star_key
        self._rows = []
        self._keys = set()

    def __contains__(self, key: str) -> bool:
        return key in self._keys

    @property
    def keys(self) -> set:
        return set(self._keys)

    def reject(self, element, ion, wavelength_air_A, reason,
               value=None, note='') -> None:
        if reason not in REJECTION_REASONS:
            raise ValueError(
                f"unknown rejection reason {reason!r}; "
                f"allowed: {REJECTION_REASONS}")
        note = '' if note is None else str(note)
        if reason == 'other' and not note.strip():
            raise ValueError("reason 'other' requires a non-empty freetext note")
        key = line_key(element, ion, wavelength_air_A)
        if key in self._keys:
            return                       # earliest gate already recorded
        self._keys.add(key)
        self._rows.append(dict(
            element=str(element),
            ion=str(ion),
            wavelength_air_A=round(float(wavelength_air_A), 5),
            reason=reason,
            value=_clean_value(value),
            note=note,
        ))

    def drop_keys(self, keys) -> None:
        """Remove already-logged rejections for these keys (measured wins)."""
        keys = set(keys)
        if not keys:
            return
        self._keys -= keys
        self._rows = [r for r in self._rows
                      if line_key(r['element'], r['ion'],
                                  r['wavelength_air_A']) not in keys]

    def to_frame(self) -> pd.DataFrame:
        cols = ['element', 'ion', 'wavelength_air_A', 'reason', 'value', 'note']
        if not self._rows:
            return pd.DataFrame(columns=cols)
        return (pd.DataFrame(self._rows)[cols]
                .sort_values(['wavelength_air_A', 'element', 'ion'])
                .reset_index(drop=True))


# ── Reconciliation + no-silent-drop invariant ─────────────────────────────────

def _measured_keys(measured_df: pd.DataFrame) -> set:
    """line_keys of rows that produced a FINITE EW (a real measurement)."""
    if measured_df is None or len(measured_df) == 0:
        return set()
    df = measured_df
    ew = pd.to_numeric(df['ew_mA'], errors='coerce')
    ok = df[np.isfinite(ew)]
    return {line_key(r['element'], r['ion'], r['wavelength_air_A'])
            for _, r in ok.iterrows()}


def reconcile(lines_priority: pd.DataFrame,
              measured_df: pd.DataFrame,
              ledger: RejectionLedger,
              spec_wav: np.ndarray,
              min_fit_depth: float) -> dict:
    """
    Close the ledger: assign a reason to every in-coverage input line that was
    not measured, then verify the no-silent-drop invariant.

    Fitter-stage drops (fit failure, EW out of range, in-window pixel gap) are
    already in `ledger`. Here we sweep the lines that never reached the fitter
    (excluded by the worklist's blend / depth pre-filter) and derive their
    reason from the linelist attributes, mirroring _build_worklist exactly. Any
    line that passed the worklist filter yet is neither measured nor already in
    the ledger is a genuine SILENT DROP: it is recorded as reason='other' with
    an UNCLASSIFIED note and counted so the invariant fails loudly.

    Returns a report dict (also serialised to the committed audit JSON).
    """
    lo, hi = float(np.min(spec_wav)), float(np.max(spec_wav))

    df = lines_priority.copy()
    df['_blend'] = (df['blend_flag'].astype(str).str.strip().str.lower() == 'true')
    df['central_depth'] = pd.to_numeric(df['central_depth'], errors='coerce')
    inc = df[(df['wavelength_air_A'] >= lo) & (df['wavelength_air_A'] <= hi)].copy()
    inc['_key'] = [line_key(e, i, w) for e, i, w in
                   zip(inc['element'], inc['ion'], inc['wavelength_air_A'])]

    measured = _measured_keys(measured_df)
    ledger.drop_keys(measured)           # measured always wins over any earlier drop log

    # Any results row with a NON-finite EW is not a measurement -> record it.
    if measured_df is not None and len(measured_df):
        ew = pd.to_numeric(measured_df['ew_mA'], errors='coerce')
        for _, r in measured_df[~np.isfinite(ew)].iterrows():
            k = line_key(r['element'], r['ion'], r['wavelength_air_A'])
            if k not in measured:
                ledger.reject(r['element'], r['ion'], r['wavelength_air_A'],
                              'failed_chi2', value=r.get('chi2'),
                              note='profile fit produced a non-finite EW')

    # Sweep in-coverage keys that are neither measured nor already rejected.
    for key, grp in inc.groupby('_key', sort=False):
        if key in measured or key in ledger:
            continue
        rep = grp.loc[grp['central_depth'].idxmax()]
        passed_worklist = ((~grp['_blend']) &
                           (grp['central_depth'] >= min_fit_depth)).any()
        if passed_worklist:
            # Entered the fitter but no drop reason was recorded -> real silent drop.
            ledger.reject(rep['element'], rep['ion'], rep['wavelength_air_A'],
                          'other', value=rep['central_depth'],
                          note=(f'{_UNCLASSIFIED_NOTE}: passed the worklist filter '
                                'but was neither measured nor rejected'))
        elif bool(rep['_blend']):
            ledger.reject(rep['element'], rep['ion'], rep['wavelength_air_A'],
                          'blend_overlap', value=np.nan,
                          note='blend_flag=True; excluded from worklist (not force-measured)')
        else:
            ledger.reject(rep['element'], rep['ion'], rep['wavelength_air_A'],
                          'below_priority_threshold', value=rep['central_depth'],
                          note=f'central_depth<{min_fit_depth} (min_fit_depth pre-fit selection gate)')

    rej = ledger.to_frame()
    rej_keys = set(ledger.keys)

    inc_keys = set(inc['_key'])
    n_inc = len(inc_keys)
    n_measured = len(measured & inc_keys)
    n_rejected = len(rej_keys & inc_keys)

    unaccounted = sorted(inc_keys - measured - rej_keys)      # must be empty
    n_unclassified = int(rej['note'].str.startswith(_UNCLASSIFIED_NOTE).sum()) \
        if len(rej) else 0
    invariant_holds = (n_measured + n_rejected == n_inc) and not unaccounted

    if len(rej):
        _in = rej.apply(lambda r: line_key(r['element'], r['ion'],
                                           r['wavelength_air_A']) in inc_keys, axis=1)
        reason_counts = {str(k): int(v) for k, v in
                         rej[_in]['reason'].value_counts().items()}
    else:
        reason_counts = {}

    return dict(
        star_key=ledger.star_key,
        coverage_A=[round(lo, 3), round(hi, 3)],
        min_fit_depth=float(min_fit_depth),
        n_in_coverage=n_inc,
        n_measured=n_measured,
        n_rejected=n_rejected,
        invariant_holds=bool(invariant_holds),
        n_unaccounted=len(unaccounted),
        unaccounted_sample=unaccounted[:20],
        n_unclassified=n_unclassified,
        reason_counts=reason_counts,
    )


# ── Ionization-stage-presence audit flag (Step 3) ─────────────────────────────

def ionization_stage_presence_check(lines_priority: pd.DataFrame,
                                    measured_df: pd.DataFrame,
                                    spec_wav: np.ndarray,
                                    emit: bool = True) -> list:
    """
    For any element whose in-coverage linelist carries BOTH a neutral (ion 'I')
    and an ionized (ion != 'I') stage, if one stage returns ZERO measured lines
    raise a LOUD audit flag (warn, do not crash). Infrastructure guard only;
    the Sr I vs Sr II ionization-balance SCIENCE stays with RYA-422.

    Returns a list of flag dicts (one per element that tripped it).
    """
    lo, hi = float(np.min(spec_wav)), float(np.max(spec_wav))
    ll = lines_priority[(lines_priority['wavelength_air_A'] >= lo) &
                        (lines_priority['wavelength_air_A'] <= hi)]

    def _neutral(ion):
        return str(ion).strip() == 'I'

    meas = measured_df if measured_df is not None else pd.DataFrame(
        columns=['element', 'ion', 'ew_mA'])
    if len(meas):
        meas = meas[np.isfinite(pd.to_numeric(meas['ew_mA'], errors='coerce'))]

    flags = []
    for elem, grp in ll.groupby('element'):
        has_neutral  = grp['ion'].map(_neutral).any()
        has_ionized  = (~grp['ion'].map(_neutral)).any()
        if not (has_neutral and has_ionized):
            continue
        m = meas[meas['element'] == elem]
        n_neutral = int(m['ion'].map(_neutral).sum())
        n_ionized = int((~m['ion'].map(_neutral)).sum())
        if n_neutral == 0 or n_ionized == 0:
            empty = 'neutral' if n_neutral == 0 else 'ionized'
            flag = dict(element=str(elem),
                        neutral_measured=n_neutral,
                        ionized_measured=n_ionized,
                        empty_stage=empty,
                        linelist_neutral=int(grp['ion'].map(_neutral).sum()),
                        linelist_ionized=int((~grp['ion'].map(_neutral)).sum()))
            flags.append(flag)
            if emit:
                msg = (f"IONIZATION-STAGE-PRESENCE flag [{elem}]: "
                       f"{empty} stage has ZERO measured lines "
                       f"(neutral measured={n_neutral}, ionized measured={n_ionized}; "
                       f"linelist has both stages in coverage). "
                       f"See RYA-422 for the ionization-balance science diagnostic.")
                warnings.warn(msg, stacklevel=2)
                print(f"  !! {msg}")
    return flags

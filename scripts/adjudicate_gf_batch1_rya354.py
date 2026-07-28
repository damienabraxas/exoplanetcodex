#!/usr/bin/env python3
"""
scripts/adjudicate_gf_batch1_rya354.py
======================================
RYA-354 Batch 1 (Na / Mg / Al) — per-line gf adjudication to best-available
NIST-graded values on the single-source gf table (data/linelists/canonical_gf.csv,
RYA-353). The RE-PRIORITIZED scope: the 10 gf-data-limited solar EW-pool elements;
Batch 1 is the high-confidence light metals.

CARDINAL RULES (RYA-354):
  - adopt CITED graded values only; the underlying reference is NAMED (never
    "VALD3"/"GES v6" as if a source). For Na/Mg/Al the underlying graded source is
    the NIST critical compilation Kelleher & Podobedova 2008 (JPCRD 37) behind NIST ASD.
  - validate-don't-tune: the differential barely moves (gf cancels); no value pulled
    toward Asplund. A line that adjudicates and still reads high is a finding.
  - uncitable -> FLAG, never guess (keep the existing sourced value, mark the flag).

Per-line NIST ASD source data (v5.11, transcribed verbatim; log gf = log10(gi*fik)):
  Al I 6696.015  fik 1.35e-2  gi 2  Acc C+  -> log gf -1.569   [JPCRD 37,709;  L4737]
  Mg I 6319.236  fik 1.58e-3  gi 3  Acc D   -> log gf -2.324   [JPCRD 37,1285; L6095]
  Mg I 6319.493  fik 5.25e-4  gi 3  Acc E   -> log gf -2.803   [JPCRD 37,1285; L6095]
  Mg I 6083, 6630, Al I 6631  -> NO citable graded NIST gf  -> FLAG.

Out: data/linelists/canonical_gf.csv (in place, Batch-1 rows only) +
     data/audit/gf_adjudication/batch1_namgal_rya354.csv (the before/after table).
"""
import csv
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
CANON = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
OUT_AUDIT = ROOT / 'data' / 'audit' / 'gf_adjudication' / 'batch1_namgal_rya354.csv'

# References are kept COMMA-FREE (the table is unquoted CSV; underlying source named).
_AL = 'NIST ASD v5.11; Kelleher & Podobedova 2008 JPCRD 37 p709 (Al/Si)'
_MG = 'NIST ASD v5.11; Kelleher & Podobedova 2008 JPCRD 37 p1285 (Na/Mg)'

# 0-based column indices in canonical_gf.csv.
C_SPECIES, C_WAVE, C_LOGGF, C_REF, C_GRADE, C_ADJ = 3, 4, 7, 8, 9, 11

# (species, wavelength) -> decision dict.
#   ADJUDICATE: set log_gf/nist_grade/loggf_reference/adjudication_status to the cited NIST value.
#   FLAG:       keep value, set adjudication_status='flagged_uncitable_rya354', note why.
#   CONFIRM:    already NIST-graded (NIST A) — no change, recorded for the audit table.
DECISIONS = {
    # --- Na I: both already NIST grade A; the blocker is saturation, NOT gf ---
    ('Na I', 5682.633): dict(action='CONFIRM', note='already NIST A; blocker=SAT+BLEND (EW 159 mA), not gf'),
    ('Na I', 5688.205): dict(action='CONFIRM', note='already NIST A; blocker=SAT+BLEND (EW 163 mA), not gf'),
    # --- Mg I ---
    ('Mg I', 5711.088): dict(action='CONFIRM', note='already NIST A; blocker=SAT (EW 91 mA, REW -4.80), not gf'),
    ('Mg I', 6083.522): dict(action='FLAG', note='NO Mg I line near 6083 in NIST ASD -> uncitable (genuine void); '
                             'keeps Kurucz KP value; also HIERR-culled'),
    ('Mg I', 6319.237): dict(action='ADJUDICATE', log_gf=-2.324, grade='D', ref=f'{_MG} grade D L6095',
                             note='NIST grade D (fik 1.58e-3, gi 3); value == prior VALD3 -2.324 (no move); '
                             'D is below the graded-reliable bar + line is SAT-culled'),
    ('Mg I', 6319.493): dict(action='ADJUDICATE', log_gf=-2.803, grade='E', ref=f'{_MG} grade E L6095',
                             note='NIST grade E (fik 5.25e-4, gi 3); value == prior -2.803 (no move); SAT-culled'),
    ('Mg I', 6630.833): dict(action='FLAG', note='NIST lists Mg I 6630.834 but with NO oscillator strength / '
                             'no graded gf -> uncitable; keeps Kurucz KP value; also HIERR-culled'),
    # --- Al I ---
    ('Al I', 6631.218): dict(action='FLAG', note='NO Al I line near 6631 in NIST ASD -> uncitable (likely mis-ID / '
                             'Kurucz-only line); keeps K75 value; flag for manual line-ID verification'),
    ('Al I', 6696.185): dict(action='ADJUDICATE', log_gf=-1.569, grade='C+', ref=f'{_AL} grade C+ L4737',
                             note='NIST Al I 6696.015 grade C+ (fik 1.35e-2, gi 2 -> log gf -1.569); prior K75 '
                             '-1.576 ~ identical (no move). NB pool wavelength 6696.185 is +0.17 A vs NIST 6696.015 '
                             '(same transition) -> wavelength-ID flag, gf adjudicated.'),
}

# columns the adjudication is ALLOWED to write — nothing else may change.
_WRITABLE = ('log_gf', 'nist_grade', 'loggf_reference', 'adjudication_status')


def _find_target(parsed, species, wave, tol=0.02):
    """Return the data-line index (0-based, excluding header) of the matching row,
    via a robust CSV parse — or None."""
    best, bi = tol + 1, None
    for i, f in enumerate(parsed):
        if f[C_SPECIES] != species:
            continue
        try:
            d = abs(float(f[C_WAVE]) - wave)
        except (ValueError, IndexError):
            continue
        if d < best:
            best, bi = d, i
    return bi if best <= tol else None


def main():
    # Read raw lines AND a robust CSV parse (handles the few comma-containing fields).
    raw = CANON.read_text().splitlines(keepends=True)
    header, body_raw = raw[0], raw[1:]
    parsed = list(csv.reader(body_raw))

    rows, changed = [], 0
    for (species, wave), dec in DECISIONS.items():
        idx = _find_target(parsed, species, wave)
        if idx is None:
            raise SystemExit(f"FATAL: {species} {wave} not found in canonical_gf — refusing to proceed")
        f = parsed[idx]
        rec = dict(species=species, wavelength_air_A=float(f[C_WAVE]), action=dec['action'],
                   loggf_before=float(f[C_LOGGF]), grade_before=f[C_GRADE], ref_before=f[C_REF],
                   adj_before=f[C_ADJ])
        if dec['action'] in ('ADJUDICATE', 'FLAG'):
            # the target rows are simple (no embedded commas) — reconstruct textually so
            # every UNCHANGED field stays byte-identical (no pandas type/format drift).
            assert len(f) == len(body_raw[idx].rstrip('\n').split(',')), \
                f"row {species} {wave} has an embedded comma — refusing textual edit"
            nf = list(f)
            if dec['action'] == 'ADJUDICATE':
                nf[C_LOGGF] = f"{round(float(dec['log_gf']), 3)}"
                nf[C_GRADE] = dec['grade']
                nf[C_REF] = dec['ref']
                nf[C_ADJ] = 'adjudicated_rya354'
            else:
                nf[C_ADJ] = 'flagged_uncitable_rya354'
            assert not any(',' in x for x in nf), "new field contains a comma"
            newline = ','.join(nf) + ('\n' if body_raw[idx].endswith('\n') else '')
            body_raw[idx] = newline
            changed += 1
        rec.update(loggf_after=round(float(nf[C_LOGGF]) if dec['action'] == 'ADJUDICATE'
                                      else float(f[C_LOGGF]), 3),
                   grade_after=(dec.get('grade') if dec['action'] == 'ADJUDICATE' else rec['grade_before']),
                   ref_after=(dec.get('ref') if dec['action'] == 'ADJUDICATE' else rec['ref_before']),
                   adj_after=('adjudicated_rya354' if dec['action'] == 'ADJUDICATE'
                              else 'flagged_uncitable_rya354' if dec['action'] == 'FLAG' else rec['adj_before']),
                   note=dec['note'])
        rec['delta_loggf'] = round(rec['loggf_after'] - rec['loggf_before'], 4)
        rows.append(rec)

    CANON.write_text(header + ''.join(body_raw))
    OUT_AUDIT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT_AUDIT, index=False)
    changed_idx = [r for r in rows if r['action'] in ('ADJUDICATE', 'FLAG')]

    n_adj = sum(r['action'] == 'ADJUDICATE' for r in rows)
    n_flag = sum(r['action'] == 'FLAG' for r in rows)
    n_conf = sum(r['action'] == 'CONFIRM' for r in rows)
    print(f"RYA-354 Batch 1 (Na/Mg/Al): {len(rows)} lines  ->  {n_adj} adjudicated, "
          f"{n_flag} flagged-uncitable, {n_conf} confirmed-already-NIST-A")
    print(f"  canonical_gf.csv: {len(changed_idx)} rows changed textually (writable cols only); "
          f"all other rows byte-identical")
    for r in rows:
        print(f"  [{r['action']:10s}] {r['species']} {r['wavelength_air_A']:.3f}  "
              f"gf {r['loggf_before']:+.3f}->{r['loggf_after']:+.3f} (d{r['delta_loggf']:+.3f})  "
              f"grade {r['grade_before']}->{r['grade_after']}")
    print(f"  audit -> {OUT_AUDIT.relative_to(ROOT)}")


if __name__ == '__main__':
    main()

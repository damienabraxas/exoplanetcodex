#!/usr/bin/env python3
"""
scripts/adjudicate_gf_batch2_mncuv_rya354.py
============================================
RYA-354 Batch 2 (Mn / Cu / V — iron-peak, HFS) — per-line gf adjudication on the
single-source gf table (data/linelists/canonical_gf.csv, RYA-353).

THE BATCH-2 FINDING (settled line-by-line from IN-REPO evidence — the RYA-395 pre-grade
and RYA-398 graded non-Fe curation diagnostics, no external transcription):

  • Mn I — GENUINELY gf-GRADE-limited, but only 3 specific lines. Of the 9 measured solar
    Mn I lines, 3 are CLEAN pre-grade (RYA-395) yet culled by the RYA-398 independent-gf
    firewall purely because their gf is Kurucz-2007 (gf_tier LOW / GF_UNVERIFIED), not an
    independently graded value: Mn I 5495.008, 5720.179, 6306.365 (weak, REW -5.1..-5.4,
    low blend, hfs_n=1 → no HFS-sum complication). The other 6 are SAT/BLEND/HIERR-culled
    (NOT gf-fixable — the Na/Mg saturation pattern). So Mn unlocks (0 -> 3 clean ->
    measured value, the Batch-1 Al precedent) IF those 3 get an independent graded gf.
    Their value EXISTS (Kurucz; GES == Kurucz), so this is GRADE-owed, not a value void —
    it needs the verbatim independent graded pull (NIST ASD / Den Hartog, Lawler, Sneden,
    Wood, Nave & Cowan 2011 ApJS 194 35, Mn I lab gf). Per the cardinal rule that pull is
    Ryan's manual action; we FLAG (flagged_manual_pull_rya354), we do NOT invent a value.

  • Cu, V — NOT gf-grade-limited. Zero measured Cu/V lines reach the curation/firewall
    stage (absent from BOTH RYA-395 and RYA-398 diagnostics) → the blocker is UPSTREAM
    (line selection / EW measurement / region matching), so gf adjudication cannot unlock
    them. (Cu additionally shows a genuine VALD-vs-GES value DISAGREEMENT, median |Δ|=0.275
    dex, all 5 lines HFS-split — a real value conflict, but moot until lines are measured.)
    This corrects the RYA-456 framing for Cu/V: not a gf void, an upstream measurement gap.

CARDINAL RULES (RYA-354): adopt CITED graded values only (underlying reference named —
never "VALD3"/"GES v6" as a source); validate-don't-tune; uncitable/unpullable -> FLAG,
never guess. NB WebFetch answers via a small summarising model that can hallucinate table
values — unacceptable for canonical atomic data, so per-line graded values are NOT pulled
in-session; they are flagged for the verbatim manual pull.

Out: data/linelists/canonical_gf.csv (adjudication_status only, the 3 Mn unlock lines) +
     data/audit/gf_adjudication/batch2_mncuv_rya354.csv (the per-line analysis table).
"""
import csv
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
CANON = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
CURATION = ROOT / 'data' / 'curation' / 'nonfe_pools'
OUT_AUDIT = ROOT / 'data' / 'audit' / 'gf_adjudication' / 'batch2_mncuv_rya354.csv'

# 0-based column indices in canonical_gf.csv (same layout as Batch 1).
C_SPECIES, C_WAVE, C_LOGGF, C_REF, C_GRADE, C_ADJ = 3, 4, 7, 8, 9, 11

# Named authoritative graded sources for the manual pull (Ryan's action) — comma-free.
_MN_SRC = 'NIST ASD; Den Hartog Lawler Sneden Wood Nave & Cowan 2011 ApJS 194 35 (Mn I lab gf)'

# The 3 Mn I lines that are CLEAN pre-grade but culled ONLY by the RYA-398 grade firewall
# (K07 Kurucz gf, GF_UNVERIFIED, LOW tier). Flag for the independent graded pull. Value is
# NOT changed (no verified independent value in-session) — adjudication_status only.
MN_UNLOCK_LINES = (5495.008, 5720.179, 6306.365)

_WRITABLE_STATUS_ONLY = 'adjudication_status'


def _find_target(parsed, species, wave, tol=0.02):
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


def _mn_per_line_audit():
    """The 9 measured solar Mn I lines + their cull reason, straight from the RYA-395/398
    in-repo curation tables — the evidence for the grade-vs-not-grade split."""
    pre = pd.read_csv(CURATION / 'Mn_cull_rya395.csv', comment='#')
    rows = []
    for _, r in pre.iterrows():
        wl = float(r['wavelength_air_A'])
        # clean pre-grade (kept) but lost to the RYA-398 grade firewall == unlock candidate
        grade_culled = any(abs(wl - u) <= 0.02 for u in MN_UNLOCK_LINES)
        rows.append(dict(species='Mn I', wavelength_air_A=round(wl, 3),
                         ew_mA=r.get('ew_mA'), rew=round(float(r['rew']), 3) if pd.notna(r.get('rew')) else None,
                         log_gf=r.get('log_gf'), gf_tier=r.get('gf_tier'),
                         loggf_reference=r.get('loggf_reference'),
                         pre_grade_kept=bool(r.get('kept')), cull_reason=r.get('cull_reason'),
                         blocker=('GRADE_ONLY -> unlock candidate' if grade_culled
                                  else 'SAT/BLEND/HIERR (not gf-fixable)'),
                         action=('FLAG_MANUAL_PULL' if grade_culled else 'NO_ACTION_not_gf_limited'),
                         manual_pull_source=(_MN_SRC if grade_culled else ''),
                         note=('clean pre-grade (RYA-395) but culled by the RYA-398 graded '
                               'firewall (K07 Kurucz, GF_UNVERIFIED) -> needs independent graded gf'
                               if grade_culled else 'culled upstream of the grade firewall')))
    return rows


def main():
    raw = CANON.read_text().splitlines(keepends=True)
    header, body_raw = raw[0], raw[1:]
    parsed = list(csv.reader(body_raw))

    # --- canonical_gf: flag the 3 Mn unlock lines (adjudication_status ONLY) ---
    changed = 0
    for wl in MN_UNLOCK_LINES:
        idx = _find_target(parsed, 'Mn I', wl)
        if idx is None:
            raise SystemExit(f"FATAL: Mn I {wl} not found in canonical_gf — refusing to proceed")
        f = parsed[idx]
        assert len(f) == len(body_raw[idx].rstrip('\n').split(',')), \
            f"Mn I {wl}: embedded comma — refusing textual edit"
        nf = list(f)
        nf[C_ADJ] = 'flagged_manual_pull_rya354'
        assert not any(',' in x for x in nf), "new field contains a comma"
        # value/grade/reference UNTOUCHED — validate-don't-tune, no invented graded value
        assert nf[C_LOGGF] == f[C_LOGGF] and nf[C_GRADE] == f[C_GRADE] and nf[C_REF] == f[C_REF]
        body_raw[idx] = ','.join(nf) + ('\n' if body_raw[idx].endswith('\n') else '')
        changed += 1
    CANON.write_text(header + ''.join(body_raw))

    # --- audit table: Mn per-line + the Cu/V upstream-blocked finding ---
    rows = _mn_per_line_audit()
    for el, nreg, dmed, hfsmax, note in (
        ('Cu I', 5, 0.275, 10, 'NOT gf-grade-limited: 0 measured lines reach curation '
         '(absent RYA-395/398) -> upstream measurement/selection blocker. ALSO genuine '
         'VALD-vs-GES value disagreement median 0.275 dex; all HFS-split. Source for a '
         'future pull: NIST ASD / Kock & Richter 1968 (Cu I).'),
        ('V I', 33, 0.003, 9, 'NOT gf-grade-limited: 0 measured lines reach curation -> '
         'upstream blocker. VALD/GES agree (median 0.003); 21/38 single-component. Source '
         'for a future pull: Lawler Wood Den Hartog Nave Sneden & Cowan 2014 ApJS 215 20 (V I).'),
        ('V II', 13, 0.003, 1, 'NOT gf-grade-limited: 0 measured lines reach curation -> '
         'upstream blocker. No HFS (single-component). Source for a future pull: Wood Lawler '
         'Sneden & Cowan 2014 ApJS 214 18 (V II).')):
        rows.append(dict(species=el, wavelength_air_A=None, ew_mA=None, rew=None, log_gf=None,
                         gf_tier=None, loggf_reference='VALD3/GES', pre_grade_kept=False,
                         cull_reason='not_measured', blocker='UPSTREAM (not gf)',
                         action='NO_ACTION_not_gf_limited', manual_pull_source='',
                         note=note))

    OUT_AUDIT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT_AUDIT, index=False)

    n_flag = sum(r['action'] == 'FLAG_MANUAL_PULL' for r in rows)
    print(f"RYA-354 Batch 2 (Mn/Cu/V): canonical_gf adjudication_status flagged on "
          f"{changed} Mn I unlock-candidate lines (value/grade UNTOUCHED).")
    print(f"  Mn I: {n_flag} GRADE-ONLY culled (unlock candidates, flagged_manual_pull) + "
          f"{sum(r['species']=='Mn I' for r in rows) - n_flag} SAT/BLEND (not gf-fixable)")
    print(f"  Cu / V: NOT gf-grade-limited — 0 measured lines reach curation (upstream blocker)")
    for r in rows:
        if r['action'] == 'FLAG_MANUAL_PULL':
            print(f"  [FLAG manual-pull] {r['species']} {r['wavelength_air_A']:.3f}  "
                  f"gf {r['log_gf']} ({r['loggf_reference']}, tier {r['gf_tier']}) -> {_MN_SRC}")
    print(f"  audit -> {OUT_AUDIT.relative_to(ROOT)}")


if __name__ == '__main__':
    main()

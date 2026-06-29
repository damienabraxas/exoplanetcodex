#!/usr/bin/env python3
"""
scripts/adjudicate_gf_batch3_zrbayeu_rya354.py
==============================================
RYA-354 Batch 3 (Zr / Ba / Y / Eu — the neutron-capture / s-process elements) —
per-line gf adjudication on the single-source gf table (data/linelists/canonical_gf.csv).

THE BATCH-3 FINDING (settled line-by-line from IN-REPO evidence — the solar EW pool
(sol_ew_results_v1.csv) + the canonical_gf provenance/grade + REW saturation, no external
transcription):

  Unlike Cu/V (Batch 2, ZERO measured lines), Zr/Ba/Y/Eu ARE measured — 12 solar lines
  (Zr 6, Y 3, Ba 2, Eu 1). But ALL four elements still produce NO abundance (CURATION-OWED,
  no value). The blocker is NOT gf grade; gf adjudication cannot unlock them. Two binding
  blockers, both upstream of gf:

    1. SATURATION — 7 of the 12 lines sit on the flat COG (REW > -4.85; Zr I 4542/4815,
       Zr II 4629/5350/5372, Y I 6191, Ba I 6498 at EW 120-266 mA). On a saturated line the
       EW is gf-INSENSITIVE, so a graded gf does not improve the abundance — the Batch-1
       Na/Mg precedent, extended to the heavy elements.

    2. REGIONS-NOT-WIRED (in_regions=False) — every line EXCEPT Ba II 5853 is absent from
       the non-Fe regions file the EW abundance path matches against, so it never reaches
       the derivation regardless of grade (the RYA-456 "regions file 0 for Eu/Ba/Sr"
       finding). These s-process lines are HFS/isotope-split (Ba II hfs_n=6, Eu II hfs_n=11)
       → the single-profile EW path is the wrong tool, the same root cause as Cu/V (the
       unlock route is HFS-resolved synthesis, RYA-466, or the RYA-161/162 differential
       survey — NOT gf grade).

  So Batch 3 ADJUDICATES the gf QUALITY (the umbrella's deliverable) without claiming an
  unlock (there is none until the measurement path is fixed — honest, the B2 pattern):

   • Ba II 5853.668 — already NIST ASD grade A (adj=nist). Best-available; CONFIRMED, no
     change. (Strong + HFS → blend/saturation-limited, not gf-limited.)
   • Eu II 6645.091 — carries LWHS = Lawler, Wickliffe, den Hartog & Sneden 2001 (ApJS 137,
     351), THE authoritative Eu II HFS lab gf (the source NIST ASD cites for Eu II). The
     value (log gf 0.4208, 11 HFS components) is already in place — STAMP the adjudication
     (pending -> adjudicated_lwhs2001). The cleanest line in the set (REW -5.99) and the
     RYA-102/458 HFS recovery target — best-available graded gf now recorded.
   • Zr II 5372.466 — carries LNAJ = Ljung, Nilsson, Asplund & Johansson 2006 (A&A 456,
     1181), the modern Zr II lab gf — STAMP adjudicated (lnaj2006). (Saturated, so moot for
     the abundance, but the gf grade is recorded best-available.)
   • Y I 4348/6191/6435 — gf is K06 (Kurucz 2006, semi-empirical / UNGRADED). No graded lab
     value in-repo -> FLAG for the verbatim manual pull (flagged_manual_pull_rya354). Named
     candidate source: Hannaford Lowe Grevesse Biemont & Whaling 1982 (ApJ 261, 736; Y I/II
     lab gf — verify coverage, Y I graded lab data is thin). Per the cardinal rule the pull
     is Ryan's manual action; we FLAG, never invent.
   • Zr I 4028/4542/4815 (VALD3 / BGHL=Biemont+1981) and Zr II 4629/5350 (CC/CCout=older
     Corliss-Cowley) — sourced values, mostly saturated; documented saturation-limited, no
     gf action (the value is cited, the blocker is not gf).

VALIDATE-DON'T-TUNE: value/grade/reference FROZEN on every row (the graded values are
already in place; we never invent or transcribe a gf). adjudication_status only, via the
RYA-354 byte-safe surgical textual edit (csv.writer would rewrite all 145k line endings; the
mixed-type columns also corrupt under pandas to_csv). Solar Fe anchor untouched (no value
changes -> bit-identical science).

Out: data/linelists/canonical_gf.csv (adjudication_status on 5 lines) +
     data/audit/gf_adjudication/batch3_zrbayeu_rya354.csv (the per-line analysis table).
"""
import csv
import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
CANON = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
EW_POOL = ROOT / 'data' / 'measured' / 'sol_ew_results_v1.csv'
OUT_AUDIT = ROOT / 'data' / 'audit' / 'gf_adjudication' / 'batch3_zrbayeu_rya354.csv'

# Saturation threshold on the reduced EW (log10 EW/lambda). Above the knee the curve of
# growth is flat -> EW is gf-insensitive (the Batch-1 Na/Mg saturation finding).
REW_SATURATED = -4.85
REW_CLEAN = -5.20

# The surgical adjudication edits (adjudication_status ONLY; value/grade/ref frozen).
#   (species, wavelength, new_status, cited_source_for_audit)
_LWHS = 'Lawler Wickliffe den Hartog & Sneden 2001 ApJS 137 351 (Eu II HFS lab gf; the NIST ASD authority)'
_LNAJ = 'Ljung Nilsson Asplund & Johansson 2006 A&A 456 1181 (Zr II lab gf)'
_YI_PULL = 'NIST ASD -> Hannaford Lowe Grevesse Biemont & Whaling 1982 ApJ 261 736 (Y I/II lab gf; verify Y I coverage - graded Y I lab data is thin)'

ADJUDICATIONS = [
    ('Eu II', 6645.0905, 'adjudicated_lwhs2001_rya354', _LWHS),
    ('Zr II', 5372.466,  'adjudicated_lnaj2006_rya354', _LNAJ),
    ('Y I',   4348.786,  'flagged_manual_pull_rya354',  _YI_PULL),
    ('Y I',   6191.718,  'flagged_manual_pull_rya354',  _YI_PULL),
    ('Y I',   6435.004,  'flagged_manual_pull_rya354',  _YI_PULL),
]

# The 12 measured solar lines (species, wavelength, EW_mA) — from sol_ew_results_v1.csv.
MEASURED = [
    ('Zr I', 4028.95, 12.51), ('Zr I', 4542.22, 202.95), ('Zr I', 4815.04, 122.1),
    ('Zr II', 4629.079, 127.24), ('Zr II', 5350.089, 237.47), ('Zr II', 5372.466, 266.13),
    ('Y I', 4348.786, 37.02), ('Y I', 6191.718, 189.94), ('Y I', 6435.004, 9.04),
    ('Ba II', 5853.668, 74.62), ('Ba I', 6498.757, 166.7), ('Eu II', 6645.127, 6.8),
]


def _canon_rows():
    rows = list(csv.reader(open(CANON)))
    return rows[0], {n: i for i, n in enumerate(rows[0])}, rows


def _rew_class(wl, ew):
    rew = math.log10(ew * 1e-3 / wl)
    return rew, ('SATURATED' if rew > REW_SATURATED else
                 ('clean' if rew < REW_CLEAN else 'strong'))


def _find(rows, c, species, wl, tol=0.05):
    for i, r in enumerate(rows[1:], start=1):
        if r[c['species']] == species:
            try:
                if abs(float(r[c['wavelength_air_A']]) - wl) < tol:
                    return i
            except (ValueError, IndexError):
                continue
    return None


def adjudicate(apply: bool = True) -> list:
    """Apply the adjudication_status edits (byte-safe textual; value/grade/ref frozen) and
    return the full per-line audit for all 12 measured lines."""
    hdr, c, rows = _canon_rows()
    raw = CANON.read_text().splitlines(keepends=True)

    edit_idx = {}
    for species, wl, status, source in ADJUDICATIONS:
        i = _find(rows, c, species, wl)
        if i is None:
            raise SystemExit(f"FATAL: {species} {wl} not in canonical_gf — refusing to proceed")
        if len(raw[i].rstrip('\n').split(',')) != len(hdr) or len(rows[i]) != len(hdr):
            raise SystemExit(f"{species} {wl}: embedded comma — refusing textual edit")
        edit_idx[(species, round(wl, 3))] = (i, status, source)

    if apply:
        for (species, wl), (i, status, _src) in edit_idx.items():
            fields = raw[i].rstrip('\n').split(',')
            assert len(fields) == len(hdr), "field-count drift"
            fields[c['adjudication_status']] = status
            assert not any(',' in f for f in fields), "new field contains a comma"
            eol = '\n' if raw[i].endswith('\n') else ''
            raw[i] = ','.join(fields) + eol
        CANON.write_text(''.join(raw))

    # Build the per-line audit (all 12 measured lines). The EW-pool wavelength can differ
    # slightly from the canonical air wavelength (e.g. Eu II 6645.127 vs 6645.0905), so the
    # adjudication is matched by the canonical ROW INDEX (_find), not an exact wavelength key.
    edits_by_index = {i: (st, src) for (i, st, src) in edit_idx.values()}
    audit = []
    for species, wl, ew in MEASURED:
        i = _find(rows, c, species, wl)
        r = rows[i] if i is not None else None
        rew, sat = _rew_class(wl, ew)
        st, src = edits_by_index.get(i, (None, ''))
        ref = r[c['loggf_reference']] if r else '?'
        grade = r[c['nist_grade']] if r else ''
        in_reg = r[c['in_regions']] if r else '?'
        gf = r[c['log_gf']] if r else '?'
        if st and st.startswith('adjudicated'):
            disp = f'STAMP graded ({ref}) — best-available cited, value frozen'
        elif st and st.startswith('flagged'):
            disp = 'FLAG manual pull — gf is ungraded Kurucz (K06)'
        elif grade == 'A':
            disp = 'CONFIRMED NIST grade A — best-available, no change'
        elif sat == 'SATURATED':
            disp = f'saturation-limited (REW {rew:.2f}) — gf-moot, no action'
        else:
            disp = f'sourced ({ref}); blocker is in_regions={in_reg}/HFS, not gf'
        audit.append(dict(species=species, wavelength_air_A=round(wl, 3), ew_mA=ew,
                          rew=round(rew, 2), saturation=sat, log_gf=gf,
                          loggf_reference=ref, nist_grade=grade, in_regions=in_reg,
                          new_adjudication_status=(st or '(unchanged)'),
                          cited_source=src, disposition=disp))
    return audit


def main():
    audit = adjudicate(apply=True)
    OUT_AUDIT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(audit).to_csv(OUT_AUDIT, index=False)

    n_stamp = sum(a['new_adjudication_status'].startswith('adjudicated') for a in audit)
    n_flag = sum(a['new_adjudication_status'].startswith('flagged') for a in audit)
    n_sat = sum(a['saturation'] == 'SATURATED' for a in audit)
    print("RYA-354 Batch 3 (Zr/Ba/Y/Eu — s-process):")
    print(f"  {len(audit)} measured solar lines; {n_sat} SATURATED (gf-moot); "
          f"{n_stamp} stamped graded; {n_flag} flagged manual-pull. value/grade FROZEN.")
    for a in audit:
        print(f"  {a['species']:6} {a['wavelength_air_A']:9.3f}  REW {a['rew']:+.2f} "
              f"{a['saturation']:9}  gf={a['log_gf']:>8} ({a['loggf_reference']})  "
              f"in_reg={a['in_regions']}  -> {a['disposition']}")
    print(f"  audit -> {OUT_AUDIT.relative_to(ROOT)}")
    print("  FINDING: gf grade is NOT the binding blocker — saturation (7/12) + "
          "in_regions=False (HFS s-process) are. No unlock; gf quality adjudicated. "
          "Eu II 6645 best-available = LWHS (Lawler+2001); Ba II 5853 = NIST-A; "
          "Zr II 5372 = LNAJ (Ljung+2006); Y I = ungraded K06 (flagged).")


if __name__ == '__main__':
    main()

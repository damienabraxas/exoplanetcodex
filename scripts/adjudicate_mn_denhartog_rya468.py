#!/usr/bin/env python3
"""
scripts/adjudicate_mn_denhartog_rya468.py
=========================================
RYA-468 — Mn I gf disposition on the single-source canonical_gf table.

Two findings (Ryan pulled the Den Hartog+2011 PDF):
  CONFIRMED VOID — Mn I 5495.008 / 5720.179 / 6306.365 are NOT in NIST ASD nor in
    Den Hartog+2011 Table 3 (which spans 2576-6021 A; 6306 is past the ceiling).
    Genuinely Kurucz-only → record the void (no value to adopt; never guessed).
  ADJUDICATE — Mn I 6013.51 / 6016.67 / 6021.82 (the e6S->z6P multiplet, EP~3 eV,
    "ideal for Mn in high-metallicity stars") to the Den Hartog+2011 lab values:
        6013.51  log gf -0.354   (was VALD3 -0.352)
        6016.67  log gf -0.181   (was VALD3 -0.183)
        6021.82  log gf -0.054   (was VALD3 -0.054)
    cited: Den Hartog, Lawler, Sobeck, Sneden & Cowan 2011, ApJS 194, 35 (unc 0.016
    dex). The value barely moves (validate-don't-tune); the REF changes VALD3 -> the
    graded lab source, so the line clears the RYA-398 independent-gf firewall.

Surgical TEXTUAL edit (RYA-354 discipline): only the matched rows change, every other
row byte-identical (no pandas round-trip), refs kept COMMA-FREE.
"""
import csv
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
CANON = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
OUT_AUDIT = ROOT / 'data' / 'audit' / 'gf_adjudication' / 'mn_denhartog_rya468.csv'

_DH = 'Den Hartog Lawler Sobeck Sneden Cowan 2011 ApJS 194 35 (unc 0.016 dex)'
C_SPECIES, C_WAVE, C_LOGGF, C_REF, C_GRADE, C_ADJ = 3, 4, 7, 8, 9, 11

DECISIONS = {
    # confirmed lab-gf void (Kurucz-only; NIST + Den Hartog out-of-coverage)
    ('Mn I', 5495.008): dict(action='VOID', note='no NIST/DH line (NIST void, DH Table 3 ends 6021 A); Kurucz-only'),
    ('Mn I', 5720.179): dict(action='VOID', note='no NIST/DH line; Kurucz-only'),
    ('Mn I', 6306.365): dict(action='VOID', note='past Den Hartog 6021 A ceiling + NIST void; Kurucz-only'),
    # adjudicate to Den Hartog+2011 (e6S->z6P multiplet)
    ('Mn I', 6013.51): dict(action='ADJUDICATE', log_gf=-0.354, note='DH+2011 e6S->z6P J2.5->1.5'),
    ('Mn I', 6016.67): dict(action='ADJUDICATE', log_gf=-0.181, note='DH+2011 e6S->z6P J2.5->2.5'),
    ('Mn I', 6021.82): dict(action='ADJUDICATE', log_gf=-0.054, note='DH+2011 e6S->z6P J2.5->3.5'),
}


def _find(parsed, species, wave, tol=0.02):
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
    raw = CANON.read_text().splitlines(keepends=True)
    header, body = raw[0], raw[1:]
    parsed = list(csv.reader(body))
    rows = []
    for (species, wave), dec in DECISIONS.items():
        idx = _find(parsed, species, wave)
        if idx is None:
            raise SystemExit(f"FATAL: {species} {wave} not in canonical_gf")
        f = parsed[idx]
        rec = dict(species=species, wavelength_air_A=float(f[C_WAVE]), action=dec['action'],
                   loggf_before=float(f[C_LOGGF]), ref_before=f[C_REF], adj_before=f[C_ADJ], note=dec['note'])
        # target rows are simple (no embedded comma) -> reconstruct textually, byte-exact elsewhere
        assert len(f) == len(body[idx].rstrip('\n').split(',')), f"{species} {wave} has an embedded comma"
        nf = list(f)
        if dec['action'] == 'ADJUDICATE':
            nf[C_LOGGF] = f"{round(float(dec['log_gf']), 3)}"
            nf[C_REF] = _DH
            nf[C_ADJ] = 'adjudicated_rya468'
        else:  # VOID
            nf[C_ADJ] = 'void_lab_gf_confirmed_rya468'
        assert not any(',' in x for x in nf), "new field has a comma"
        body[idx] = ','.join(nf) + ('\n' if body[idx].endswith('\n') else '')
        rec.update(loggf_after=float(nf[C_LOGGF]), ref_after=nf[C_REF], adj_after=nf[C_ADJ],
                   delta_loggf=round(float(nf[C_LOGGF]) - rec['loggf_before'], 4))
        rows.append(rec)

    CANON.write_text(header + ''.join(body))
    OUT_AUDIT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT_AUDIT, index=False)
    n_adj = sum(r['action'] == 'ADJUDICATE' for r in rows)
    n_void = sum(r['action'] == 'VOID' for r in rows)
    print(f"RYA-468 Mn: {n_adj} adjudicated to Den Hartog+2011, {n_void} recorded confirmed-void "
          f"(only these {len(rows)} rows changed; all else byte-identical)")
    for r in rows:
        print(f"  [{r['action']:10s}] Mn I {r['wavelength_air_A']:.3f}  gf {r['loggf_before']:+.3f}->"
              f"{r['loggf_after']:+.3f} (d{r['delta_loggf']:+.3f})  ref->{r['ref_after'][:32]}")
    print(f"  audit -> {OUT_AUDIT.relative_to(ROOT)}")


if __name__ == '__main__':
    main()

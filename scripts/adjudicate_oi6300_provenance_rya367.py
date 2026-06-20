#!/usr/bin/env python3
"""
scripts/adjudicate_oi6300_provenance_rya367.py
==============================================
RYA-367 Step 1 — make the [O I] 6300.304 gf provenance cite its PRIMARY source.

The canonical value (log gf = -9.717) is already correct and already a cited
primary (NIST ASD v5.11 grade A) — NOT the generic 'VALD3 / -9.776' that the
stale store-#2 (linelist_solar.csv) carries. NIST's transition probability for
the [O I] 6300.30 forbidden line is the Storey & Zeippen 2000 (MNRAS 312, 813)
value, log gf = -9.717. This makes that primary reference explicit inline (the
VALUE is unchanged — this is a provenance pin, not a value adjudication).

Surgical single-line edit (no pandas.to_csv round-trip → clean 1-line diff).
Idempotent. `--check` verifies without writing.
"""
import argparse
import csv
import io
import sys
from pathlib import Path

_CANON = Path(__file__).resolve().parents[1] / 'data' / 'linelists' / 'canonical_gf.csv'
OI_LINE_ID = 'gf_108773'        # [O I] 6300.304, EP 0.0 (verified below)

NEW_REF = 'NIST ASD v5.11 grade A (Storey & Zeippen 2000, MNRAS 312, 813)'
NEW_STATUS = 'adjudicated_rya367'
EXPECTED_LOG_GF = '-9.717'      # S&Z 2000 — value MUST already match (no change)


def _header_and_row():
    lines = _CANON.read_text().split('\n')
    header = next(csv.reader([lines[0]]))
    for ln in lines[1:]:
        if ln and next(csv.reader([ln]))[0] == OI_LINE_ID:
            return header, dict(zip(header, next(csv.reader([ln])))), ln
    raise SystemExit(f"FAIL: line_id {OI_LINE_ID} not found — refusing to edit.")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    args = ap.parse_args(argv)

    header, rec, old_line = _header_and_row()
    # identity + value guard (no value change, no silent mis-edit)
    assert rec['species'].strip() == 'O I', rec['species']
    assert abs(float(rec['wavelength_air_A']) - 6300.304) <= 0.01, rec['wavelength_air_A']
    assert rec['log_gf'] == EXPECTED_LOG_GF, (
        f"[O I] log_gf is {rec['log_gf']}, expected {EXPECTED_LOG_GF} (S&Z 2000). "
        f"This script only pins provenance; a value mismatch is a separate finding.")

    print(f"[O I] 6300.304 (line_id={OI_LINE_ID}):")
    print(f"  log_gf={rec['log_gf']} (unchanged)  ref={rec['loggf_reference']!r} "
          f"status={rec['adjudication_status']!r}")
    print(f"  target ref={NEW_REF!r}")

    is_current = (rec['loggf_reference'] == NEW_REF
                  and rec['adjudication_status'] == NEW_STATUS)
    if args.check:
        print("CHECK: OK." if is_current else "CHECK: provenance not yet pinned to S&Z 2000.")
        return 0 if is_current else 1
    if is_current:
        print("Already pinned (idempotent).")
        return 0

    new_rec = dict(rec)
    new_rec['loggf_reference'] = NEW_REF
    new_rec['adjudication_status'] = NEW_STATUS
    buf = io.StringIO()
    csv.writer(buf, lineterminator='').writerow([new_rec[c] for c in header])
    text = _CANON.read_text()
    assert text.count(old_line) == 1, "old [O I] row not uniquely located — aborting"
    _CANON.write_text(text.replace(old_line, buf.getvalue(), 1))
    print(f"  → pinned provenance to S&Z 2000 (value -9.717 unchanged); "
          f"surgical 1-line edit to {_CANON.name}")
    return 0


if __name__ == '__main__':
    sys.exit(main())

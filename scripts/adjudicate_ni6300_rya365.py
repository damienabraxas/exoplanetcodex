#!/usr/bin/env python3
"""
scripts/adjudicate_ni6300_rya365.py
===================================
RYA-365 — adjudicate the canonical gf of Ni I 6300.34 to the dedicated lab
measurement for this [O I]-blended line.

Finding (RYA-365, RCA of the RYA-237 solar-O FAIL):
  The single-source migration (RYA-353) seeded Ni I 6300.342 from a generic
  "NIST ASD v5.11 grade B" value, log gf = -2.310. But for THIS specific line —
  the Ni I blend under the [O I] 6300 forbidden line used for the solar O
  abundance — the dedicated laboratory authority is:

      Johansson, Litzén, Lundberg & Zhang 2003, ApJ 584, L107
      (arXiv:astro-ph/0301382) — experimental f-value + isotopic structure for
      the Ni I line blended with [O I] at 6300 Å: log gf = -2.11.

  The GES synth seed already carried -2.110 (i.e. Johansson 2003); the migration
  overwrote it with the 0.2-dex-weaker NIST value, which under-strengthens the Ni
  component, so the [O I]+Ni joint fit over-attributes the blend to O → A(O)☉
  reads ~0.17 dex high. This is an RYA-354-class per-line adjudication to the
  best-available value, NOT a tune-to-anchor: -2.11 is chosen because it is the
  ticket-named dedicated lab reference (and matches the GES seed), independent of
  the O outcome; solar O landing in its gate is the confirmation.

Idempotent. Asserts exactly one matching row. Prints before/after.

Usage:
    python3 scripts/adjudicate_ni6300_rya365.py            # apply
    python3 scripts/adjudicate_ni6300_rya365.py --check    # verify only (no write)
"""
import argparse
import csv
import io
import sys
from pathlib import Path

_CANON = Path(__file__).resolve().parents[1] / 'data' / 'linelists' / 'canonical_gf.csv'

# Target line (by stable line_id) + adjudicated value ---------------------------
NI_LINE_ID = 'gf_101075'        # Ni I 6300.342, EP 4.266 (verified by the asserts below)

NEW_VALUES = {
    'log_gf': '-2.11',
    'loggf_reference': 'Johansson, Litzen, Lundberg & Zhang 2003, ApJ 584, L107 '
                       '(arXiv:astro-ph/0301382)',
    'nist_grade': '',                      # lab measurement, not a NIST grade
    'seed_source': 'Johansson2003',
    'adjudication_status': 'adjudicated_rya365',
}
STALE_LOG_GF = '-2.31'

# Surgical, single-line text edit (the canonical table is ~137k rows; a pandas
# round-trip reformats every float and produces a useless whole-file diff). We
# rewrite exactly the one Ni I 6300.34 row, leaving every other byte untouched.


def _read_lines():
    text = _CANON.read_text()
    return text.split('\n'), text.endswith('\n')


def _header_and_row():
    lines, _ = _read_lines()
    header = next(csv.reader([lines[0]]))
    for ln in lines[1:]:
        if not ln:
            continue
        row = next(csv.reader([ln]))
        if row and row[0] == NI_LINE_ID:
            return header, dict(zip(header, row)), ln
    raise SystemExit(f"FAIL: line_id {NI_LINE_ID} not found — refusing to edit.")


def _assert_target(rec: dict):
    # Verify identity before touching it (no silent mis-edit / fabrication).
    assert rec['species'].strip() == 'Ni I', rec['species']
    assert abs(float(rec['wavelength_air_A']) - 6300.342) <= 0.02, rec['wavelength_air_A']
    assert abs(float(rec['excitation_potential_eV']) - 4.266) <= 0.02, rec['excitation_potential_eV']


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true',
                    help='verify the canonical value only; non-zero exit if stale')
    args = ap.parse_args(argv)

    if not _CANON.exists():
        raise SystemExit(f"canonical gf table not found: {_CANON}")
    header, rec, old_line = _header_and_row()
    _assert_target(rec)

    print(f"Ni I 6300.34 (line_id={NI_LINE_ID}):")
    print(f"  current  log_gf={rec['log_gf']}  ref={rec['loggf_reference']!r}  "
          f"status={rec['adjudication_status']!r}")
    print(f"  GES seed gf_synth_ges={rec['gf_synth_ges']}  "
          f"VALD gf_linelist_vald={rec['gf_linelist_vald']}")
    print(f"  target   log_gf={NEW_VALUES['log_gf']}  ref={NEW_VALUES['loggf_reference']!r}")

    is_current = (abs(float(rec['log_gf']) - float(NEW_VALUES['log_gf'])) < 1e-6
                  and rec['seed_source'] == NEW_VALUES['seed_source'])

    if args.check:
        if is_current:
            print("CHECK: canonical Ni I 6300.34 = Johansson 2003 (-2.11). OK.")
            return 0
        print("CHECK: canonical Ni I 6300.34 is NOT the Johansson 2003 value (stale).")
        return 1

    if is_current:
        print("Already adjudicated to Johansson 2003 — no change (idempotent).")
        return 0

    delta = float(NEW_VALUES['log_gf']) - float(rec['log_gf'])
    new_rec = dict(rec)
    new_rec.update(NEW_VALUES)
    # Re-serialize the single row with the SAME csv dialect (quote only when needed).
    buf = io.StringIO()
    csv.writer(buf, lineterminator='').writerow([new_rec[c] for c in header])
    new_line = buf.getvalue()

    text = _CANON.read_text()
    assert text.count(old_line) == 1, "old row not uniquely located — aborting"
    _CANON.write_text(text.replace(old_line, new_line, 1))
    print(f"  → wrote log_gf={NEW_VALUES['log_gf']} (Δ {delta:+.3f} dex, stronger Ni); "
          f"surgical single-line edit to {_CANON.name}")
    return 0


if __name__ == '__main__':
    sys.exit(main())

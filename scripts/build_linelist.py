#!/usr/bin/env python3
"""
scripts/build_linelist.py
=========================
Build linelist_master.csv from a VALD3 long-format extract.

Usage
-----
    python scripts/build_linelist.py \\
        --star 55cnc \\
        --vald data/linelists/vald_55cnc_raw.txt

    # With explicit overrides
    python scripts/build_linelist.py \\
        --star 55cnc \\
        --vald data/linelists/vald_55cnc_raw.txt \\
        --nist data/linelists/nist_crosscheck.csv \\
        --elements data/config/elements_master.json \\
        --out data/linelists/linelist_master.csv \\
        --min-depth 0.0

Output schema matches data/linelists/README.md column spec.
NIST-only science lines (depth < 1% threshold) are injected regardless of
--min-depth, since they are required for oxygen/carbon/phosphorus science.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import PIPELINE

IONMAP = {'1': 'I', '2': 'II', '3': 'III'}

OUTPUT_COLS = [
    'element', 'ion', 'wavelength_air_A', 'excitation_potential_eV',
    'log_gf', 'loggf_source', 'nist_grade', 'damping_rad', 'damping_stark',
    'damping_vdW', 'central_depth', 'blend_flag', 'priority', 'notes',
]

# Science-critical lines absent from VALD (predicted depth < 1% at 55 Cnc params).
# These are injected unconditionally from NIST ASD values.
NIST_INJECTIONS = [
    {
        'element': 'O', 'ion': 'I', 'wavelength_air_A': 6300.304,
        'excitation_potential_eV': 0.000, 'log_gf': -9.717,
        'loggf_source': 'NIST', 'nist_grade': 'A',
        'damping_rad': 0.000, 'damping_stark': 0.000, 'damping_vdW': 0.000,
        'central_depth': 0.005, 'blend_flag': True, 'priority': 1,
        'notes': 'Forbidden line; primary oxygen indicator; subtract Ni I 6300.336 contribution before EW measurement',
    },
    {
        'element': 'Ni', 'ion': 'I', 'wavelength_air_A': 6300.336,
        'excitation_potential_eV': 4.266, 'log_gf': -2.110,
        'loggf_source': 'NIST', 'nist_grade': 'B',
        'damping_rad': 7.800, 'damping_stark': -5.600, 'damping_vdW': -7.800,
        'central_depth': 0.004, 'blend_flag': True, 'priority': 2,
        'notes': 'O I 6300 blend partner; model jointly with O I 6300.304',
    },
    {
        'element': 'C', 'ion': 'I', 'wavelength_air_A': 5380.337,
        'excitation_potential_eV': 7.685, 'log_gf': -1.616,
        'loggf_source': 'NIST', 'nist_grade': 'B',
        'damping_rad': 7.950, 'damping_stark': -5.700, 'damping_vdW': -7.700,
        'central_depth': 0.007, 'blend_flag': False, 'priority': 1,
        'notes': 'High excitation (7.685 eV); sensitive to Teff errors ±50 K',
    },
    {
        'element': 'P', 'ion': 'I', 'wavelength_air_A': 6034.040,
        'excitation_potential_eV': 1.873, 'log_gf': -0.310,
        'loggf_source': 'NIST', 'nist_grade': 'C',
        'damping_rad': 7.900, 'damping_stark': -5.400, 'damping_vdW': -7.600,
        'central_depth': 0.006, 'blend_flag': False, 'priority': 1,
        'notes': 'Phosphorus tracer; requires S/N > 300; NIST grade C only',
    },
    {
        'element': 'P', 'ion': 'I', 'wavelength_air_A': 6043.120,
        'excitation_potential_eV': 1.873, 'log_gf': -0.460,
        'loggf_source': 'NIST', 'nist_grade': 'C',
        'damping_rad': 7.900, 'damping_stark': -5.400, 'damping_vdW': -7.600,
        'central_depth': 0.005, 'blend_flag': False, 'priority': 1,
        'notes': 'Phosphorus tracer; requires S/N > 300; NIST grade C only',
    },
]


def parse_vald(path: Path) -> list:
    """
    Parse VALD3 long-format output into a list of line dicts.

    VALD long format: 4 lines per transition.
    Data line starts with 'El Num', and has 13 numeric fields:
        wavelength, log_gf, E_low, J_low, E_up, J_up,
        lande_lower, lande_upper, lande_mean,
        rad_damp, stark_damp, waals_damp, central_depth
    """
    lines = []
    n_skipped = 0

    with open(path, encoding='utf-8', errors='replace') as f:
        for raw in f:
            line = raw.strip()
            if not line.startswith("'"):
                continue

            end_q = line.find("',")
            if end_q == -1:
                continue

            spec_parts = line[1:end_q].split()
            if len(spec_parts) != 2:
                continue
            elem, ion_num = spec_parts[0], spec_parts[1]

            rest = line[end_q + 2:]
            fields = [x.strip() for x in rest.split(',') if x.strip()]
            if len(fields) < 13:
                n_skipped += 1
                continue

            try:
                wl    = float(fields[0])
                loggf = float(fields[1])
                ep    = float(fields[2])
                # fields[3] = J_low, fields[4] = E_up, fields[5] = J_up  (skip)
                # fields[6,7,8] = lande factors                           (skip)
                rad   = float(fields[9])
                stark = float(fields[10])
                waals = float(fields[11])
                depth = float(fields[12])
            except (ValueError, IndexError):
                n_skipped += 1
                continue

            lines.append({
                'element': elem,
                'ion': IONMAP.get(ion_num, ion_num),
                'wavelength_air_A': wl,
                'excitation_potential_eV': ep,
                'log_gf': loggf,
                'loggf_source': 'VALD3',
                'nist_grade': '',
                'damping_rad': rad,
                'damping_stark': stark,
                'damping_vdW': waals,
                'central_depth': depth,
                'blend_flag': False,
                'priority': 0,
                'notes': '',
            })

    if n_skipped:
        print(f"  (skipped {n_skipped} malformed data lines)")
    return lines


def load_elements(path: Path) -> dict:
    """Load elements_master.json. Returns {symbol: element_dict}."""
    with open(path) as f:
        data = json.load(f)
    return {e['symbol']: e for e in data['elements']}


def filter_and_assign(lines: list, elements: dict) -> list:
    """
    Assign priorities from elements_master.json.
    Non-target elements are kept (priority=0) for blend coverage.
    Element filtering happens at load time via loader.py.
    """
    for row in lines:
        meta = elements.get(row['element'])
        row['priority'] = meta['priority'] if meta else 0
    return lines


def flag_blends(lines: list, tol_A: float = 0.10) -> list:
    """Flag any line within tol_A Å of another line as a blend candidate."""
    if not lines:
        return lines
    sorted_lines = sorted(lines, key=lambda r: r['wavelength_air_A'])
    for i in range(1, len(sorted_lines)):
        sep = sorted_lines[i]['wavelength_air_A'] - sorted_lines[i - 1]['wavelength_air_A']
        if sep < tol_A:
            sorted_lines[i]['blend_flag'] = True
            sorted_lines[i - 1]['blend_flag'] = True
    return sorted_lines


def crosscheck_nist(lines: list, nist_path: Path) -> list:
    """
    Apply NIST grades from nist_crosscheck.csv to matching lines.
    Matches on element, ion, and wavelength within 0.010 Å.
    Skips lines beginning with '#' (comment lines in the CSV header).
    """
    import io
    nist_rows = []
    with open(nist_path) as f:
        non_comment = [l for l in f if not l.lstrip().startswith('#')]
    reader = csv.DictReader(io.StringIO(''.join(non_comment)))
    for row in reader:
        try:
            nist_rows.append({
                'element': row['element'].strip(),
                'ion': row['ion'].strip(),
                'wl': float(row['wavelength_air_A']),
                'grade': row['nist_grade'].strip(),
            })
        except (KeyError, ValueError):
            continue

    for row in lines:
        for nr in nist_rows:
            if (row['element'] == nr['element']
                    and row['ion'] == nr['ion']
                    and abs(row['wavelength_air_A'] - nr['wl']) < 0.010):
                row['nist_grade'] = nr['grade']
                break
    return lines


def inject_nist_lines(lines: list, tol_A: float = 0.005) -> tuple:
    """
    Add NIST-only science lines not present in VALD extract.
    Returns (updated_lines, n_injected).
    Existing lines within tol_A of an injection wavelength are skipped.
    """
    existing_wls = [r['wavelength_air_A'] for r in lines]
    n_injected = 0
    for inj in NIST_INJECTIONS:
        wl = inj['wavelength_air_A']
        already_present = any(abs(wl - ew) < tol_A for ew in existing_wls)
        if not already_present:
            lines.append(dict(inj))
            existing_wls.append(wl)
            n_injected += 1
    return lines, n_injected


def write_master(lines: list, out_path: Path) -> None:
    """Write master CSV sorted by wavelength."""
    lines.sort(key=lambda r: r['wavelength_air_A'])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLS, extrasaction='ignore')
        writer.writeheader()
        for row in lines:
            writer.writerow(row)
    print(f"Wrote {len(lines):,} lines → {out_path}")


def run_qa(lines: list) -> None:
    """Print a QA summary of the built line list."""
    from collections import Counter
    elements = Counter(r['element'] for r in lines)
    sources = Counter(r['loggf_source'] for r in lines)
    grades = Counter(r['nist_grade'] or '(none)' for r in lines)
    priorities = Counter(r['priority'] for r in lines)
    n_blends = sum(1 for r in lines if r['blend_flag'])

    print(f"\n── QA Summary ──────────────────────────────────────────")
    print(f"  Total lines   : {len(lines):,}")
    print(f"  Source        : {dict(sources)}")
    print(f"  NIST grades   : {dict(grades)}")
    print(f"  Blend flags   : {n_blends}")
    print(f"  By priority   : {dict(sorted(priorities.items()))}")
    print(f"  Elements ({len(elements)}): {sorted(elements)}")
    wls = [r['wavelength_air_A'] for r in lines]
    print(f"  Wav range     : {min(wls):.3f} – {max(wls):.3f} Å")
    print(f"────────────────────────────────────────────────────────\n")


def main():
    parser = argparse.ArgumentParser(
        description='Build linelist_master.csv from VALD3 long-format extract.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--star',      default='55cnc',
                        help='Star identifier (default: 55cnc)')
    parser.add_argument('--vald',      default='data/linelists/vald_55cnc_raw.txt',
                        help='VALD3 long-format input file')
    parser.add_argument('--nist',      default='data/linelists/nist_crosscheck.csv',
                        help='NIST crosscheck CSV for grade assignment')
    parser.add_argument('--elements',  default='data/config/elements_master.json',
                        help='Elements master JSON (target element list + priorities)')
    parser.add_argument('--out',       default='data/linelists/linelist_master.csv',
                        help='Output master CSV')
    parser.add_argument('--min-depth', type=float, default=0.0,
                        help='Drop VALD lines with central_depth below this threshold '
                             '(NIST-injected lines are always kept; default: 0.0 = keep all)')
    parser.add_argument('--qa',        action='store_true',
                        help='Print QA summary after build')
    args = parser.parse_args()

    vald_path     = ROOT / args.vald
    nist_path     = ROOT / args.nist
    elements_path = ROOT / args.elements
    out_path      = ROOT / args.out

    # ── Parse ────────────────────────────────────────────────────────────────
    print(f"[1/6] Parsing VALD3: {vald_path.name}")
    raw = parse_vald(vald_path)
    print(f"      {len(raw):,} transitions parsed")

    # ── Load element targets ─────────────────────────────────────────────────
    print(f"[2/6] Loading elements: {elements_path.name}")
    elements = load_elements(elements_path)
    print(f"      {len(elements)} target elements")

    # ── Filter to target elements ────────────────────────────────────────────
    print(f"[3/6] Filtering to target elements")
    lines = filter_and_assign(raw, elements)
    n_target = sum(1 for r in lines if r['priority'] > 0)
    n_nontarget = len(lines) - n_target
    print(f"      {n_target:,} target-element lines (priority assigned); {n_nontarget} non-target kept for blend coverage")

    # ── Flag blends ──────────────────────────────────────────────────────────
    print(f"[4/6] Flagging blends (threshold: 0.10 Å)")
    lines = flag_blends(lines)
    n_blends = sum(1 for r in lines if r['blend_flag'])
    print(f"      {n_blends:,} blend flags set")

    # ── Apply NIST grades ────────────────────────────────────────────────────
    print(f"[5/6] Applying NIST grades from: {nist_path.name}")
    lines = crosscheck_nist(lines, nist_path)
    n_graded = sum(1 for r in lines if r['nist_grade'])
    print(f"      {n_graded} lines graded")

    # ── Inject NIST-only science lines ───────────────────────────────────────
    print(f"[5/6] Injecting NIST-only science lines")
    lines, n_injected = inject_nist_lines(lines)
    print(f"      {n_injected} lines injected (O I, Ni I, C I, P I x2)")

    # ── Depth filter (VALD lines only) ───────────────────────────────────────
    if args.min_depth > 0:
        n_before = len(lines)
        lines = [r for r in lines
                 if r['central_depth'] >= args.min_depth or r['loggf_source'] == 'NIST']
        print(f"[6/6] Depth filter (>= {args.min_depth}): "
              f"{n_before - len(lines):,} lines removed, {len(lines):,} kept")
    else:
        print(f"[6/6] No depth filter applied")

    # ── Write ────────────────────────────────────────────────────────────────
    write_master(lines, out_path)

    if args.qa:
        run_qa(lines)


if __name__ == '__main__':
    main()

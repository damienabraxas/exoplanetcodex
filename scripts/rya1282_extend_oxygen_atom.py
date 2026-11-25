#!/usr/bin/env python3
"""Build a radiative-only O41f extension for experimental wavelength tests."""
import argparse
import csv
import hashlib
import json
import math
import re
import sys

import numpy as np
from scipy.constants import Rydberg
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.multi_atom import read_multi_atom

TERMS = {('5p', '3P'), ('5d', '5D'), ('6s', '5S'), ('5d', '3D')}
TRIPLET_6S = ('6s', '3S')
# Configuration, spectroscopic designation and LS-series key in NORAD O I.
NORAD_TERMS = {
    ('5p', '3P'): ('d3P', (3, 1, 0, 4)),
    ('5d', '5D'): ('x5Do', (5, 2, 1, 3)),
    ('6s', '5S'): ('w5So', (5, 0, 1, 4)),
    ('5d', '3D'): ('w3Do', (3, 2, 1, 4)),
    TRIPLET_6S: ('w3So', (3, 0, 1, 4)),
}
DOI = '10.1051/0004-6361/202245645'


def _table(path):
    rows = (line for line in path.read_text().splitlines() if line and not line.startswith('#'))
    reader = csv.DictReader(rows, delimiter='\t')
    return [{k: v.strip() for k, v in r.items()} for r in reader if r['recno'].strip().isdigit()]


def _key(config, term, j, parity):
    return config.split('.')[-1], term, int(j), parity


def _norad_cross_sections(directory, terms):
    energies = (directory / 'o1.en.ls.txt').read_text()
    rows = (directory / 'o1.ptpx.txt').read_text().splitlines()
    headers = [(i, tuple(map(int, line.split()))) for i, line in enumerate(rows)
               if re.fullmatch(r'\s*(?:\d+\s+){3}\d+\s*', line)]
    by_key = {key: (start, end) for (start, key), (end, _) in zip(headers, headers[1:])}
    sections = {}
    for term in terms:
        designation, key = NORAD_TERMS[term]
        match = re.search(rf'^\s+2s22p3 \(4So\)\s+{term[0]}\s+{designation}\s+(-[\d.]+)\s*$',
                          energies, re.MULTILINE)
        if match is None:
            raise ValueError(f'Unidentified NORAD term {term}')
        start, end = by_key[key]
        binding = float(rows[start + 2].split()[0])
        count = int(rows[start + 1].split()[1])  # second count includes appended high-energy points
        table = np.array([list(map(float, row.split())) for row in rows[start + 3:end] if row.strip()])
        if (table.shape != (count, 2) or not np.isfinite(table).all()
                or np.any(np.diff(table[:, 0]) <= 0) or np.any(table[:, 1] < 0)
                or abs(binding + float(match[1])) > 5.1e-6
                or abs(binding - table[0, 0]) > 5.1e-8):
            raise ValueError(f'Invalid NORAD ground-target continuum for {term}')
        sections[term] = (key, table)
    return sections


def build(base, levels_path, lines_path, output, norad=None, *, add_6s_triplet=False):
    atom = read_multi_atom(base)
    if add_6s_triplet:
        if atom.symbol != 'O' or len(atom.levels) != 53 or norad is None:
            raise ValueError('The 6s triplet trial needs an O53 atom and NORAD cross sections')
        terms, expected = {TRIPLET_6S}, 1
    elif (atom.symbol, len(atom.levels), len(atom.lines), len(atom.continua)) != ('O', 41, 256, 22):
        raise ValueError('Expected the O41f baseline')
    else:
        terms, expected = TERMS, 12
    first_ion = next(level.index for level in atom.levels if level.stage == 2)
    levels = {_key(r['State'], r['Term'], r['J'], r['Parity']): r for r in _table(levels_path)}
    added = sorted(
        (key for key in levels if key[:2] in terms),
        key=lambda key: (float(levels[key]['ENIST']), key),
    )
    if len(added) != expected:
        raise ValueError(f'Expected {expected} resolved levels for {sorted(terms)}')
    index = {}
    for key, row in levels.items():
        if not row['ENIST']:
            continue
        energy = float(row['ENIST'])
        matches = [level.index for level in atom.levels if level.stage == 1
                   and abs(level.E_cm - energy) < 0.002 and level.g == 2 * key[2] + 1]
        if len(matches) > 1:
            raise ValueError(f'Ambiguous O41f level for {key}')
        if matches:
            index[key] = matches[0]

    if add_6s_triplet and not TERMS.issubset({key[:2] for key in index}):
        raise ValueError('Base atom is missing one or more O53 terms')
    added_set = set(added)
    new_levels = []
    for i, key in enumerate(added, first_ion):
        if key in index:
            raise ValueError(f'Already present in O41f: {key}')
        index[key] = i
        energy = float(levels[key]['ENIST'])
        new_levels.append(f"{energy:14.3f} {2 * key[2] + 1:7.1f} 'Li2023 {' '.join(map(str, key))}' 1")

    new_lines, links, omitted = [], [], []
    seen = set()
    for r in _table(lines_path):
        upper = _key(r['Upper'], r['Termu'], r['Ju'], r['Pu'])
        lower = _key(r['Lower'], r['Terml'], r['Jl'], r['Pl'])
        if not (upper in added_set or lower in added_set):
            continue
        if upper not in index or lower not in index:
            omitted.append(int(r['recno']))
            continue
        j, i = index[upper], index[lower]
        if (j, i) in seen:
            raise ValueError(f'Duplicate E1 endpoints: {upper}, {lower}')
        seen.add((j, i))
        wavelength_nm = 1e7 / (float(levels[upper]['ENIST']) - float(levels[lower]['ENIST']))
        if not math.isclose(wavelength_nm, float(r['lambda']), abs_tol=0.00006):
            raise ValueError(f'E1 wavelength mismatch at row {r["recno"]}')
        f = 10 ** float(r['loggf-B']) / (2 * lower[2] + 1)
        natural = sum(1 / float(levels[k]['taul']) for k in (upper, lower) if levels[k]['taul'])
        new_lines.append(f'{j:4} {i:4} {f:.8E}  31  15.0  0.0  0  {natural:.8E}  0  0  VOIGT')
        links.append(dict(li_row=int(r['recno']), upper=j, lower=i, wavelength_vac_nm=wavelength_nm,
                          loggf_b=float(r['loggf-B']), li_class=r['Acc']))

    # Keep stage ordering and remap the two original ionic endpoints in all rate blocks.
    def shifted(i):
        i = int(i)
        return i + len(added) if i >= first_ion else i

    def remap(line):
        return re.sub(r'^(\s*)(\d+)(\s+)(\d+)',
                      lambda m: f'{m[1]}{shifted(m[2])}{m[3]}{shifted(m[4])}', line)

    bf_rows, bf_links = [], []
    if norad is not None:
        ion = shifted(first_ion)
        ion_energy = atom.level(first_ion).E_cm
        sections = _norad_cross_sections(norad, terms)
        wavelength_floor = min(float(c.wavelength[-1]) for c in atom.continua)
        for key in added:
            source_key, table = sections[key[:2]]
            edge_cm = ion_energy - float(levels[key]['ENIST'])
            # Use the observed threshold and preserve the tabulated excess photon
            # energy. LS approximation: each resolved J gets the same cross section.
            photon_cm = edge_cm + (table[:, 0] - table[0, 0]) * (Rydberg / 100)
            wavelength = 1e8 / photon_cm
            keep = wavelength >= wavelength_floor
            wavelength = wavelength[keep]
            sigma = table[keep, 1] * 1e-18  # Mb to cm²
            i = index[key]
            bf_rows.append(f'{ion:4} {i:4} {sigma[0]:.8E} {len(wavelength):5} -1.0')
            bf_rows.extend(f'{w:.10f} {x:.8E}' for w, x in zip(wavelength, sigma))
            bf_links.append(dict(upper=ion, lower=i, norad_key=source_key,
                                 edge_vac_nm=wavelength[0] * .1, points=len(wavelength),
                                 omitted_shortwave_points=len(table) - len(wavelength)))
    rows = base.read_text().splitlines()
    data = [i for i, row in enumerate(rows) if row.strip() and not row.lstrip().startswith('*')]
    rows[data[2]] = f'    {len(atom.levels) + len(added)} {len(atom.lines) + len(links)} {len(atom.continua) + len(bf_links)} 0'
    continuum = 3 + len(atom.levels) + len(atom.lines)
    for _ in atom.continua:
        k = data[continuum]
        rows[k] = remap(rows[k])
        continuum += 1 + int(rows[k].split()[3])
    gencol = next(n for n in range(continuum, len(data)) if rows[data[n]].strip() == 'GENCOL')
    collision_kinds = {collision.kind for collision in atom.collisions}
    for n in range(gencol + 1, len(data) - 1):
        if rows[data[n]].split()[0] in collision_kinds:
            k = data[n + 1]
            rows[k] = remap(rows[k])
    if bf_rows:
        bf_end = data[gencol - 1] + 1
        rows[bf_end:bf_end] = ['* NORAD O II ground-target photoionization'] + bf_rows
    bb_end = data[3 + len(atom.levels) + len(atom.lines) - 1] + 1
    rows[bb_end:bb_end] = ['* Li 2023 E1 links for added levels'] + new_lines
    ion_start = data[3 + first_ion - 1]
    rows[ion_start:ion_start] = new_levels
    output.mkdir(parents=True, exist_ok=False)
    path = output / ('atom.o54-radiative-sensitivity' if add_6s_triplet else 'atom.o53-radiative')
    path.write_text('* EXPERIMENTAL RADIATIVE EXTENSION: NO NEW COLLISIONS\n'
                    + '\n'.join(rows) + '\n')
    inputs = [base, levels_path, lines_path]
    if norad is not None:
        inputs.extend(norad / name for name in ('o1.en.ls.txt', 'o1.ptpx.txt'))
    receipt = dict(status='RADIATIVE_ONLY_HOLD', li_doi=DOI,
                   inputs={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in inputs},
                   added_levels=[dict(index=index[k], state=k, energy_cm=float(levels[k]['ENIST']))
                                 for k in added],
                   added_links=links, omitted_li_rows=omitted, added_ground_target_continua=bf_links,
                   norad_level_doi='10.1103/PhysRevA.58.3766' if bf_links else None,
                   norad_partial_source='1999ApJS..120..131N' if bf_links else None,
                   missing=['excited-target photoionization and recombination' if bf_links
                            else 'photoionization and recombination', 'electron and hydrogen collisions',
                            'pressure broadening', 'radiative links to levels absent or collapsed in O41f'])
    (output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', required=True, type=Path)
    parser.add_argument('--li-dir', required=True, type=Path)
    parser.add_argument('--norad-dir', type=Path)
    parser.add_argument('--add-6s-triplet', action='store_true',
                        help='add the 6s triplet term radiatively to an experimental O53 atom')
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    print(build(args.base, args.li_dir / 'tablea1.tsv', args.li_dir / 'tablea2.tsv', args.out,
                norad=args.norad_dir, add_6s_triplet=args.add_6s_triplet))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Attach an experimental LS O+H excitation table to the radiative O53 atom."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
from scipy.constants import Boltzmann, c, h

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.multi_atom import read_multi_atom

NEW_TERMS = {'5p.3P', '6s.5So', '5d.5Do', '5d.3Do'}
CM_TO_K = 100 * h * c / Boltzmann


def _terms(atom):
    terms = defaultdict(list)
    for level in atom.levels:
        if level.stage != 1:
            continue
        if level.label.startswith('Li2023 '):
            _, orbital, term, _, parity = level.label.split()
            term = f'{orbital}.{term}' + ('o' if parity == '-' else '')
        else:
            match = re.fullmatch(r'Level\s+\d+ =\s+(\d[spd])([35][SPD])\d?(\*)?', level.label)
            if match is None:
                continue
            orbital, spin_term, odd = match.groups()
            term = f'{orbital}.{spin_term}' + ('o' if odd else '')
        terms[term].append(level)
    return terms


def _read_rates(path):
    lines = path.read_text().splitlines()
    _, _, *columns = lines[0].split()
    temperatures = np.array(columns, dtype=float)
    if (len(temperatures) < 2 or not np.isfinite(temperatures).all()
            or temperatures[0] <= 0 or np.any(np.diff(temperatures) <= 0)):
        raise ValueError('Expected an increasing positive temperature axis')
    pairs = {}
    for line in lines[1:]:
        lower, upper, *rates = line.split()
        rate = np.array(rates, dtype=float)
        if ((lower, upper) in pairs or rate.shape != temperatures.shape
                or not np.isfinite(rate).all() or np.any(rate < 0)):
            raise ValueError(f'Invalid O+H excitation row {lower} -> {upper}')
        pairs[lower, upper] = rate
    return temperatures, pairs


def add_rates(atom_path, rates_path, out, expected_pairs=None):
    atom = read_multi_atom(atom_path)
    terms = _terms(atom)
    if len(atom.levels) != 53 or not NEW_TERMS.issubset(terms):
        raise ValueError('Expected the radiative O53 candidate')
    temps, pairs = _read_rates(rates_path)
    new_pairs = {pair: rate for pair, rate in pairs.items() if NEW_TERMS.intersection(pair)}
    if expected_pairs is not None and len(new_pairs) != expected_pairs:
        raise ValueError(f'Expected {expected_pairs} new LS pairs, found {len(new_pairs)}')
    if not new_pairs:
        raise ValueError('No new O+H pairs')
    old = {(min(row.i, row.j), max(row.i, row.j)) for row in atom.collisions if row.kind == 'CH'}
    rows, links = [], []
    for (lower, upper), up in new_pairs.items():
        if lower not in terms or upper not in terms:
            raise ValueError(f'Unresolved O53 LS pair {lower} -> {upper}')
        g_upper = sum(level.g for level in terms[upper])
        for low in terms[lower]:
            for high in terms[upper]:
                if (low.index, high.index) in old or high.E_cm <= low.E_cm:
                    raise ValueError(f'Overlapping or unordered CH pair {low.index} -> {high.index}')
                # Populate final J by its fraction of the LS statistical weight.
                # MULTI CH then stores de-excitation. Use measured J splittings
                # so the adapter recovers detailed balance at the component level.
                down = up * (low.g / g_upper) * np.exp((high.E_cm - low.E_cm) * CM_TO_K / temps)
                if not np.isfinite(down).all():
                    raise ValueError(f'Nonfinite downward rates for {high.index} -> {low.index}')
                rows.extend(['CH', f'{high.index:4} {low.index:4} ' +
                             ' '.join(f'{k:.9e}' for k in down)])
                links.append([high.index, low.index])
                old.add((low.index, high.index))

    source_lines = atom_path.read_text().splitlines()
    end = next(i for i, line in enumerate(source_lines) if line.strip() == 'END')
    temp_row = f'{len(temps)} ' + ' '.join(map(str, temps))
    source_lines[end:end] = ['* Experimental O+H impulse rates for new levels', 'TEMP', temp_row, *rows]
    out.mkdir(parents=True, exist_ok=False)
    atom_out = out / 'atom.o53-hydrogen'
    atom_out.write_text('* DIAGNOSTIC O53: EXPERIMENTAL HYDROGEN COLLISIONS\n' +
                        '\n'.join(source_lines) + '\n')
    receipt = dict(status='H_IMPULSE_DIAGNOSTIC_HOLD', ls_pairs=len(new_pairs), ch_components=links,
                   temperatures_K=temps.tolist(), fine_structure='statistical final-J excitation weights',
                   inputs={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in (atom_path, rates_path)})
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return atom_out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--atom', type=Path, required=True)
    parser.add_argument('--rates', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--expected-pairs', type=int)
    args = parser.parse_args()
    print(add_rates(args.atom, args.rates, args.out, args.expected_pairs))


if __name__ == '__main__':
    main()

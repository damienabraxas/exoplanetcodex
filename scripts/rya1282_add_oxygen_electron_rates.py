#!/usr/bin/env python3
"""Build an allowed electron excitation sensitivity variant of experimental O53."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.special import exp1

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.multi_atom import read_multi_atom
from scripts.rya1282_add_oxygen_hydrogen_rates import CM_TO_K, NEW_TERMS, _terms

# Van Regemorter 1962, Table 2 (neutral targets), with its limiting forms.
Y = np.array([.005, .01, .02, .04, .1, .2, .4, 1., 2., 4., 10.])
SMALL_Y_FACTOR = np.sqrt(3) / (2 * np.pi)
P = np.array([SMALL_Y_FACTOR * exp1(.005), 1.160, .956, .758, .493, .331,
              .209, .100, .063, .040, .023])
TEMP_K = np.array([2000., 3000., 4000., 5000., 6000., 7000., 8000.,
                   10000., 12000., 15000., 20000.])


def neutral_gaunt(y):
    y = np.asarray(y)
    if not np.isfinite(y).all() or np.any(y <= 0):
        raise ValueError('Expected a positive excitation threshold over kT')
    middle = np.exp(np.interp(np.log(y), np.log(Y), np.log(P)))
    return np.where(y < Y[0], SMALL_Y_FACTOR * exp1(y),
                    np.where(y > Y[-1], .066 / np.sqrt(y), middle))


def multi_ce(delta_cm, f_abs, temperatures):
    # Equation 18 gives q_down = 20.60 lambda_cm^3 A P(y) / sqrt(T).
    # Substitute A = 6.670e15 (g_low/g_high) f / lambda_A^2
    # and divide out MULTI's (g_low/g_high) sqrt(T), in cgs.
    wavelength_a = 1e8 / delta_cm
    return (20.60 * 6.670e-9 * wavelength_a * f_abs
            * neutral_gaunt(delta_cm * CM_TO_K / temperatures) / temperatures)


def add_electron_rates(path, out):
    atom = read_multi_atom(path)
    terms = _terms(atom)
    term_by_level = {level.index: name for name, levels in terms.items() for level in levels}
    if len(atom.levels) != 53 or not NEW_TERMS.issubset(terms):
        raise ValueError('Expected an experimental O53 candidate')
    old = {(row.j, row.i) for row in atom.collisions if row.kind == 'CE'}
    rows, links = [], []
    for line in atom.lines:
        low_term = term_by_level.get(line.i)
        high_term = term_by_level.get(line.j)
        if not NEW_TERMS.intersection((low_term, high_term)) or not (low_term and high_term):
            continue
        if low_term.split('.')[1][0] != high_term.split('.')[1][0]:
            continue  # The LS dipole recipe does not supply spin-changing collisions.
        if (line.j, line.i) in old:
            raise ValueError(f'Existing CE pair {line.j}->{line.i}')
        delta = atom.level(line.j).E_cm - atom.level(line.i).E_cm
        rates = multi_ce(delta, line.f, TEMP_K)
        rows.extend(['CE', f'{line.j:4} {line.i:4} ' + ' '.join(f'{k:.9e}' for k in rates)])
        links.append([line.j, line.i])
        old.add((line.j, line.i))
    if not links:
        raise ValueError('No supported new E1 pairs')
    data = path.read_text().splitlines()
    end = next(i for i, line in enumerate(data) if line.strip() == 'END')
    data[end:end] = ['* Experimental Van Regemorter neutral O I electron excitation',
                     'TEMP', f'{len(TEMP_K)} ' + ' '.join(map(str, TEMP_K)), *rows]
    out.mkdir(parents=True, exist_ok=False)
    candidate = out / 'atom.o53-electron-sensitivity'
    candidate.write_text('* DIAGNOSTIC O53: EMPIRICAL LS-ALLOWED ELECTRON EXCITATION\n'
                         + '\n'.join(data) + '\n')
    report = dict(status='ELECTRON_EXCITATION_SENSITIVITY_HOLD', ce_components=links,
                  temperatures_K=TEMP_K.tolist(), source='1962ApJ...136..906V Eq. 18 and Table 2',
                  missing=['spin-changing and forbidden electron excitation', 'electron ionization'],
                  original_atom_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    (out / 'receipt.json').write_text(json.dumps(report, indent=2) + '\n')
    return candidate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--atom', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    print(add_electron_rates(args.atom, args.out))


if __name__ == '__main__':
    main()

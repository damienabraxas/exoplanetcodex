#!/usr/bin/env python3
"""Add three NIST [O I] transitions between levels already represented by Bride."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.multi_atom import read_multi_atom

# Ritz air wavelength in the held NIST export, upper and lower O41f/O53 indices.
BRANCHES = [(5577.339, 5, 4), (6363.776, 4, 2), (6391.733, 4, 3)]
DECAY_REFS = {'T4539,T5081'}


def build(source, nist, out):
    atom = read_multi_atom(source)
    rows = []
    for table in nist:
        with table.open() as file:
            rows.extend(csv.DictReader(file, delimiter='\t'))
    branches = {}
    for wavelength, j, i in BRANCHES:
        target = [r for r in rows if r['wavelength_ritz_A'] and
                  abs(float(r['wavelength_ritz_A']) - wavelength) < .00001 and
                  r['ref_transition_probability'] in DECAY_REFS]
        if not target or any(line.j == j and line.i == i for line in atom.lines):
            raise ValueError(f'Missing NIST data or duplicate atom transition for {wavelength}')
        for r in target:
            lower, upper = map(float, r['ei_ek_raw'].split('  -'))
            if (abs(atom.level(i).E_eV-lower)>1e-5 or abs(atom.level(j).E_eV-upper)>1e-5
                    or atom.level(i).g!=float(r['gi'])):
                raise ValueError(f'Level identity mismatch at {wavelength}')
        branches[wavelength] = (j, i, target)

    def decay_sum(index):
        energy = atom.level(index).E_eV
        return sum(float(r['aki_s-1']) for r in rows
                   if r['aki_s-1'] and r['ei_ek_raw'] and
                   r['ref_transition_probability'] in DECAY_REFS and
                   abs(float(r['ei_ek_raw'].split('  -')[1])-energy)<1e-5)

    added, receipt = [], []
    for wavelength, (j, i, ritz) in branches.items():
        f = sum(float(row['fik']) for row in ritz)
        width = decay_sum(j)+decay_sum(i)
        if f<=0 or width<=0:
            raise ValueError(f'Invalid NIST rates at {wavelength}')
        added.append(f'{j:4} {i:4} {f:.9e}  251  100.0  3.0  0  {width:.9e}  0  0  VOIGT')
        receipt.append(dict(air_A=wavelength, upper=j, lower=i, f=f, nominal_radiative_width_s=width,
                            A_components_s=[float(row['aki_s-1']) for row in ritz]))
    lines = source.read_text().splitlines()
    data = [i for i,line in enumerate(lines) if line.strip() and not line.lstrip().startswith('*')]
    nlev, nlin, ncont, nfix = map(int, lines[data[2]].split())
    lines[data[2]] = f'{nlev} {nlin+len(added)} {ncont} {nfix}'
    end = data[3+nlev+nlin-1]+1
    lines[end:end] = ['* NIST forbidden [O I] lines']+added
    out.mkdir(exist_ok=False)
    path = out/'atom.o53-forbidden'
    path.write_text('\n'.join(lines)+'\n')
    (out/'receipt.json').write_text(json.dumps(dict(status='EXPERIMENTAL_HOLD', lines=receipt,
        inputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (source, *nist)}),indent=2)+'\n')
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--atom', type=Path, required=True)
    p.add_argument('--nist', type=Path, action='append', required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    print(build(a.atom,a.nist,a.out))


if __name__ == '__main__':
    main()

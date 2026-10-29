"""Compare all held O41 CH0 entries with pinned author forward/reverse tables."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_charge_contract import read_author_rates, TERM_MAP
    import lightweaver.constants as c

    directory = require_codex_path('experiments.bride_barklem_oxygen')
    grid, tables, receipt = read_author_rates(directory)
    atom_path = require_codex_path('experiments.bride_oxygen_atom')
    if fingerprint(atom_path)['sha256'] != 'd058b75bdec80c5317523a6d7f185a64ca59144d3d494a17eeabee5f10d04b57':
        raise ValueError('Changed source atom')
    atom = read_multi_atom(atom_path)
    rows = []
    for levels, state, term_g in TERM_MAP:
        if sum(atom.levels[j-1].g for j in levels) != term_g:
            raise ValueError('Incomplete statistical-weight group')
        for j in levels:
            raw = next(r for r in atom.collisions if r.kind == 'CH0' and r.j == j)
            if raw.i != 40 or not np.array_equal(raw.temperature, grid):
                raise ValueError('Changed raw endpoints/grid')
            lower, ion = atom.levels[j-1], atom.levels[39]
            for index, temperature in enumerate(grid):
                forward = tables[index, 18, state-1]
                mn = tables[index, state-1, 18]
                # Candidate encoding, independently sourced reverse and raw O energies.
                predicted = mn*ion.g/lower.g*np.exp(
                    -(ion.E_cm-lower.E_cm)*c.HC_CM/(c.KBoltzmann*temperature))
                if min(forward, mn, predicted) <= 0:
                    raise ValueError('Missing positive author source coefficient')
                value = raw.rates[index]
                rows.append(dict(level_1based=j, level_label=lower.label,
                    author_state_1based=state, temperature_K=temperature,
                    raw_CH0_cm3_s=value, author_ion_pair_cm3_s=forward,
                    author_MN_term_cm3_s=mn,
                    reverse_Boltzmann_candidate_cm3_s=predicted,
                    raw_over_candidate=value/predicted,
                    raw_over_author_ion_pair=value/forward,
                    fine_structure_g=lower.g, term_g=term_g))
    if len(rows) != 240:
        raise ValueError('Incomplete source comparison')
    ratio = np.array([r['raw_over_candidate'] for r in rows])
    direct = np.array([r['raw_over_author_ion_pair'] for r in rows])
    report = dict(schema='codex.bride_charge_contract/1',scientific_admission='HOLD',
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        generator=fingerprint(__file__),adapter=fingerprint(ROOT/'pipeline/bride_charge_contract.py'),
        atom=fingerprint(atom_path),author_acquisition=receipt,entries=len(rows),
        candidate='raw_CH0 ~= q_MN_term * g_OII/g_OI_J * exp[-(E_OII-E_OI_J)/kT]',
        candidate_ratio_range=[float(ratio.min()),float(ratio.max())],
        candidate_max_absolute_relative_difference=float(np.max(abs(ratio-1))),
        raw_over_published_ion_pair_range=[float(direct.min()),float(direct.max())],
        conclusion='Strong numerical identification of an encoded reverse/Boltzmann relationship, not exact producer identity. No coefficient changes or rate adoption.',
        source='https://github.com/barklem/public-data/tree/cd2f54ed46687659017b343763e95478ed93653e/inelastic-O%2BH',
        assumptions='Manual LS/J term mapping explicit in TERM_MAP; source tables rounded. No best-fit factor, energy offset or abundance fit.')
    args.output.mkdir(parents=True,exist_ok=False)
    with (args.output/'coefficient_comparison.csv').open('x') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    report['comparison']=fingerprint(args.output/'coefficient_comparison.csv')
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()

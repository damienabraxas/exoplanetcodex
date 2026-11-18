"""Per-line OpenMULTI3D 1D NLTE corrections at the ORCH 1D-LTE abundance, any element.

Delta = A_NLTE(EW*) - A_LTE, with EW* = EW_LTE(A_LTE) on the engine's own LTE curve of growth
(differential within one engine, as in Amarsi+2019 / Engine A). Inputs: the cog.csv of
run_openmulti3d_cog.py, its spec (line windows), and ORCH *_1D-LTE_lines.csv files. Each ORCH
line maps to the spec line whose integration window contains it.

    python scripts/analyze_openmulti3d_nlte_corrections.py --cog <dir> --spec <spec.json> \
        --lte <file> [<file> ...] --output <dir>
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def delta(cog, atom, line, a_lte):
    d = cog[(cog.atom == atom) & (cog.line == line)].sort_values('A')
    if not d.A.min() <= a_lte <= d.A.max():
        raise ValueError(f'{line} A_LTE {a_lte} outside grid')
    log_ew = np.interp(a_lte, d.A, np.log10(d.ew_lte_mA))
    if not np.all(np.diff(np.log10(d.ew_nlte_mA)) > 0):
        raise ValueError(f'{atom} {line} NLTE COG not monotonic')
    a_nlte = float(np.interp(log_ew, np.log10(d.ew_nlte_mA), d.A))
    if not d.A.min() < a_nlte < d.A.max():
        raise ValueError(f'{line} A_NLTE {a_nlte} at grid edge')
    return a_nlte - a_lte, 10 ** log_ew


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cog', type=Path, required=True)
    p.add_argument('--spec', type=Path, required=True)
    p.add_argument('--lte', type=Path, nargs='+', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    cog = pd.read_csv(a.cog / 'cog.csv')
    spec = json.loads(a.spec.read_text())
    atoms = list(spec['atoms'])

    def spec_line(wl):
        hits = [l['name'] for l in spec['lines'] if l['integrate_air_A'][0] <= wl <= l['integrate_air_A'][1]]
        if len(hits) != 1:
            raise ValueError(f'{wl}: {len(hits)} spec windows')
        return hits[0]

    rows = []
    for f in a.lte:
        lte = pd.read_csv(f)
        ea_path = Path(str(f).replace('_1D-LTE_lines', '_ENGINE-A_lines'))
        ea = pd.read_csv(ea_path).set_index('wavelength_air_A') if ea_path.exists() else None
        hold = f.name.split('_SYNTH')[0].split('solar_atlas_')[-1]
        for _, r in lte[lte.in_aggregate].iterrows():
            row = dict(holding=hold, wavelength_air_A=r.wavelength_air_A, engine_line=spec_line(r.wavelength_air_A),
                       A_LTE=r.abundance, engineA_delta=None if ea is None else round(float(ea.loc[r.wavelength_air_A, 'nlte_delta_dex']), 4))
            for at in atoms:
                d, ew = delta(cog, at, row['engine_line'], r.abundance)
                row[f'delta_{at}'] = round(d, 4)
                row[f'A_NLTE_{at}'] = round(r.abundance + d, 4)
            row['engine_lte_ew_mA'] = round(ew, 2)
            rows.append(row)
    t = pd.DataFrame(rows)
    a.output.mkdir(parents=True, exist_ok=False)
    t.to_csv(a.output / 'nlte_corrections.csv', index=False)
    summary = {}
    for (hold, is_triplet), g in t.groupby(['holding', t.wavelength_air_A < 7800]):
        key = f"{hold}:{'777' if is_triplet else '8446'}"
        summary[key] = dict(n=len(g), A_LTE_mean=round(g.A_LTE.mean(), 4),
                            A_engineA_mean=None if g.engineA_delta.isna().any() else round((g.A_LTE + g.engineA_delta).mean(), 4),
                            **{f'A_NLTE_{at}_mean': round(g[f'A_NLTE_{at}'].mean(), 4) for at in atoms})
    (a.output / 'report.json').write_text(json.dumps(dict(
        schema='codex.openmulti3d_nlte_corrections/1', generator_sha256=sha(__file__), cog_report=json.loads((a.cog / 'report.json').read_text()),
        cog_sha256=sha(a.cog / 'cog.csv'), lte_inputs={str(f): sha(f) for f in a.lte}, summary=summary,
        definition='delta = A_NLTE(EW_LTE(A_LTE)) - A_LTE on OpenMULTI3D flux COGs (log-EW interpolation)'), indent=2) + '\n')
    print(t.to_string(index=False))
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()

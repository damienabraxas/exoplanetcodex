"""Read-only wider Teff/gravity preflight against actual held collision domains."""
import argparse,json,sys,platform
from pathlib import Path
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Mac-only candidate preflight')
    from config.constants import require_codex_path
    from pipeline.bride_reference import read_marcs_structure,fingerprint
    from pipeline.bride_parameter_checks import interpolate_marcs_structure,marcs_filename,assert_collision_coverage
    from pipeline.multi_atom import read_multi_atom
    plan=json.loads(args.plan.read_text());atom_path=require_codex_path('experiments.bride_oxygen_atom');atom=read_multi_atom(atom_path)
    records=[];sources={}
    for case in plan['cases']:
        nodes=[]
        for t in case['bracket_teff']:
            for g in case['bracket_logg']:
                path=require_codex_path('grids.marcs_standard')/marcs_filename(t,g);sources[path.name]=fingerprint(path)
                nodes.append(((t,g),read_marcs_structure(path)))
        s,interpolation=interpolate_marcs_structure(nodes,(case['teff'],case['logg']),thermodynamic_scale='log_positive_layers')
        T=s['temperature_K'];bad=[]
        for c in atom.collisions:
            grid=atom.temperature if c.temperature is None else c.temperature
            mask=(T<grid[0])|(T>grid[-1])
            if mask.any():bad.append(dict(kind=c.kind,i=c.i,j=c.j,domain_K=[float(grid[0]),float(grid[-1])],depth_indices=np.flatnonzero(mask).tolist()))
        try:assert_collision_coverage(atom,T)
        except ValueError as exc:
            if not bad:raise
            disposition='REFUSED_COLLISION_DOMAIN';reason=str(exc)
        else:
            if bad:raise ValueError('Independent domain enumeration disagrees')
            disposition='DOMAIN_COVERED_ONLY';reason='Temperature coverage only; no provider, RT/SE or physical atmosphere admission'
        records.append(dict(case=case,disposition=disposition,reason=reason,ndep=s['ndep'],temperature_range_K=[float(min(T)),float(max(T))],interpolation=interpolation,uncovered_table_count=len(bad),uncovered_by_kind=dict(Counter(c['kind']for c in bad)),uncovered=bad))
    args.output.mkdir(parents=True,exist_ok=False)
    report=dict(schema='codex.bride_parameter_domain/1',scientific_admission='HOLD',generator=fingerprint(__file__),plan=fingerprint(args.plan),atom=fingerprint(atom_path),
        interpolator=fingerprint(ROOT/'pipeline/bride_parameter_checks.py'),parser=fingerprint(ROOT/'pipeline/multi_atom.py'),source_models=list(sources.values()),cases=records,
        caveats='Read-only candidate preflight. Covered temperature domain is necessary, not sufficient for valid opacity/EOS/binding/atmosphere or oxygen solution. No extrapolation, truncation or new runs admitted.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    for r in records:print(r['case']['name'],r['temperature_range_K'],r['disposition'],r['uncovered_table_count'])

if __name__=='__main__':main()

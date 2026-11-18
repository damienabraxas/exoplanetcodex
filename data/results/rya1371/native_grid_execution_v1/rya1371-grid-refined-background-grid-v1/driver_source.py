"""Union the inherited full O grid with native background-opacity table knots."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint,load_plan_engine
    from pipeline.hop_background_knots import read_background_knots
    from pipeline.hop_packet import MAX_PACKET_WAVELENGTHS
    from config.constants import require_codex_path
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['replay','plan','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();plan=json.loads(a.plan.read_text());rd=json.loads((a.replay/'report.json').read_text())
    if not rd['computational_completed']:raise ValueError('Incomplete inherited replay')
    native=require_codex_path(plan['native_analysis_resource']);checks={str(Path(__file__)):fingerprint(__file__),str(a.plan.resolve()):fingerprint(a.plan),str(a.replay/'report.json'):fingerprint(a.replay/'report.json'),str(a.replay/'rates.npz'):fingerprint(a.replay/'rates.npz'),str(ROOT/'pipeline/hop_background_knots.py'):fingerprint(ROOT/'pipeline/hop_background_knots.py')}
    for name in ['arrays.npz','report.json']:
        path=native/name;checks[str(path)]=fingerprint(path)
        if rd['inputs'].get(str(path))!=fingerprint(path):raise ValueError('Wrong background-grid state')
    meta=rd['context_audit']['background_records'][0]['packet']['inputs']
    candidates=[Path(x)for x in meta if x.endswith('/input/absdat/TS_absdat.txt')]
    if len(candidates)!=1:raise ValueError('Ambiguous qualified opacity table')
    source=candidates[0];checks[str(source)]=fingerprint(source)
    if checks[str(source)]!=meta[str(source)]:raise ValueError('Changed qualified opacity table')
    _,app,_,_=load_plan_engine(dict(engine='openmulti3d',profile_variant='accurate-kernel-control'))
    for name in ['composition_mod.f90','bftrans_mod.f90']:checks[str(app/'HOP'/name)]=fingerprint(app/'HOP'/name)
    tables=read_background_knots(source);knots=np.unique(np.concatenate([t['wave_vac_nm']for t in tables]))
    with np.load(a.replay/'rates.npz')as d:old=d['wave_vac_nm'].copy();pop=d['populations'].copy();star=d['nstar'].copy()
    selected=knots[(knots>=old[0])&(knots<=old[-1])];wave=np.unique(np.concatenate([old,selected]))
    if len(wave)>MAX_PACKET_WAVELENGTHS:raise ValueError('Source-knot union exceeds unchanged20000-point cap')
    np.testing.assert_array_equal(wave[np.searchsorted(wave,old)],old)
    components=[dict(label=t['label'],header_line_1based=t['header_line_1based'],knots_vac_nm=t['wave_vac_nm'].tolist())for t in tables]
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed source-knot input')
    a.output.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(a.output/'background_grid.npz',wave_vac_nm=wave,fine_wave_vac_nm=old,populations=pop,nstar=star)
    report=dict(generator=fingerprint(__file__),inputs=checks,computational_completed=True,scientific_admission='HOLD',grid_control=dict(method='vacuum-nm-background-knots',inherited_points=len(old),table_count=len(tables),unique_table_knots=len(knots),knots_inside_support=len(selected),next_points=len(wave),scope='All source absdat table knots inside inherited support, no selected species/interval or reference-fit. Native get_table_sp clamps/sorts table knots; here every unique knot is retained without close-node merging. KVADL is interpolation order, not AIR conversion. No opacity/state/source change; no removed inherited nodes.'),components=components)
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(report['grid_control']))


if __name__=='__main__':main()

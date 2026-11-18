"""Dry source-preserving native continuum-grid proposal; not a solver plan."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint,load_engine
    from pipeline.hop_background_knots import read_background_knots
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();checks={}
    def pin(path):
        path=Path(path);checks[str(path)]=fingerprint(path);return path
    pin(__file__);pin(ROOT/'pipeline/hop_background_knots.py')
    run=require_codex_path('experiments.openmulti3d_oxygen_shared16_run');doc=json.loads(pin(run/'report.json').read_text())
    if not doc['computational_completed']:raise ValueError('Incomplete native baseline')
    old_path=pin(run/'native/save/all_lam.bin')
    entry=next(r for r in doc['output_files']if r['path']=='native/save/all_lam.bin')
    if fingerprint(old_path)!=entry['receipt']:raise ValueError('Changed native grid')
    old=np.fromfile(old_path,dtype='<f8')
    if len(old)!=39841 or not np.isfinite(old).all()or np.any(np.diff(old)<=0):raise ValueError('Wrong native baseline grid')
    _,app,_,_=load_engine('openmulti3d')
    for path in [app/'HOP/spectrum_mod.f90',app/'HOP/op_math_mod.f90',app/'SRC/atom_mod.f90',app/'SRC/ray_patch_mod.f90']:pin(path)
    text=(app/'HOP/spectrum_mod.f90').read_text()
    import re
    match=re.search(r'\blam_tolerance\s*=\s*([0-9.]+)',text)
    if match is None:raise ValueError('Missing source merge tolerance')
    source_literal=float(match[1])
    params=pin(run/'native/params.nml').read_text()
    consumed=re.search(r'\bLAM_TOLERANCE\s*=\s*([0-9.Ee+\-]+)',params)
    if consumed is None:raise ValueError('Missing consumed native merge tolerance')
    tol=float(consumed[1])
    if source_literal!=.001 or tol!=float(np.float32(.001)):raise ValueError('Changed qualified native merge tolerance')
    table_lines=pin(run/'native/atoms/O.txt').read_text().splitlines();start=next(i for i,s in enumerate(table_lines)if s.split()==['i','lam_air','alpha','j'])+1;tables={}
    for line in table_lines[start:]:
        f=line.split()
        if len(f)!=4:break
        tables.setdefault((int(f[0])+1,int(f[3])+1),[]).append(float(f[1]))
    if len(tables)!=22:raise ValueError('Wrong continuum graph')
    active=np.zeros(len(old)-1,bool)
    for values in tables.values():active|=(old[:-1]>=min(values))&(old[1:]<=max(values))
    midpoints=old[:-1][active]+np.diff(old)[active]*.5
    source=pin(app/'input/absdat/TS_absdat.txt');components=read_background_knots(source)
    knots=np.unique(np.concatenate([r['wave_vac_nm']*10 for r in components]));knots=knots[(knots>=old[0])&(knots<=old[-1])]
    # Reject proposed additions too near an inherited node or accepted addition,
    # so native grid coalescing cannot replace an original frequency sample.
    accepted=[];rejected=0
    for candidate in np.unique(np.r_[knots,midpoints]):
        i=np.searchsorted(old,candidate)
        near=(i<len(old)and abs(old[i]-candidate)<=tol)or(i>0 and abs(old[i-1]-candidate)<=tol)
        if near or(accepted and candidate-accepted[-1]<=tol):rejected+=1;continue
        accepted.append(candidate)
    wave=np.sort(np.r_[old,accepted]);np.testing.assert_array_equal(wave[np.searchsorted(wave,old)],old)
    a.output.mkdir(parents=True,exist_ok=False)
    grid=a.output/'frequency_grid_vac_A.txt';np.savetxt(grid,wave,fmt='%.17e',header='Experimental native grid proposal: vacuum Angstrom; all original nodes retained')
    np.testing.assert_array_equal(np.loadtxt(grid),wave)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed native-grid proposal input')
    report=dict(generator=fingerprint(__file__),inputs=checks,computational_completed=True,scientific_admission='HOLD',original_points=len(old),proposed_points=len(wave),added_points=len(accepted),continuum_panels_with_midpoint_proposals=int(active.sum()),source_background_tables=len(components),in_support_background_knots=len(knots),rejected_close_candidates=rejected,source_merge_tolerance_A=tol,all_original_nodes_exact=True,grid=fingerprint(grid),point_count_ratio=float(len(wave)/len(old)),scope='Dry proposal only, not an executable native plan or RT/SE packet. All original39841nodes/endpoints plus all56table opacity-knot and all22BF midpoint candidates, respecting native0.001A merging. No source atom/cross sections/atmosphere/data edit. Intended source line_grid_file hook replaces line grid; using the full original union avoids silently removing BB nodes, but actual native reproduction/memory/window invariance must still be tested. Native original600s/one-thread/2GB/tolerance controls remain required; point ratio is not measured memory/runtime. Existing HOP20000 packet cap is unchanged; this native grid is not sent to HOP and no engine guard is bypassed. Before execution, add a strictly validated opt-in adapter/control and verify emitted native grid retains originals. No adopted physical/numerical result.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({key:value for key,value in report.items()if key!='inputs'}))


if __name__=='__main__':main()

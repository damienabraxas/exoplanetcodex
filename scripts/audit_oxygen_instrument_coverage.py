"""Map every registered instrument/mode against explicit current Oxygen evidence."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint,output_grid_vacuum_A
    from pipeline.instrument_model_coverage import catalog_interval_vac_nm,interval_intersections
    from pipeline.multi_atom import read_multi_atom
    from config.constants import require_codex_path
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--probe-analysis',type=Path,required=True);p.add_argument('--probe-run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();checks={}
    paths=[ROOT/'data/catalog/instrument_catalog.csv',ROOT/'data/catalog/instrument_modes.csv'];catalogs=[]
    for path in paths:
        with path.open(newline='')as f:catalogs.append(list(csv.DictReader(f)))
        checks[str(path)]=fingerprint(path)
    instruments,modes=catalogs;ids={r['instrument_id']for r in instruments}
    if len(ids)!=len(instruments)or any(r['instrument_id']not in ids for r in modes):raise ValueError('Duplicate/unknown catalog instrument identity')
    atom_path=require_codex_path('experiments.bride_oxygen_atom');atom=read_multi_atom(atom_path);centres=np.array([1e7/(atom.level(l.j).E_cm-atom.level(l.i).E_cm)for l in atom.lines]);checks[str(atom_path)]=fingerprint(atom_path)
    replay=require_codex_path('experiments.hop_native_physical_replay')/'rates.npz'
    with np.load(replay)as d:grid=d['wave_vac_nm'].copy()
    checks[str(replay)]=fingerprint(replay)
    for directory in [a.probe_analysis,a.probe_run]:checks[str(directory/'report.json')]=fingerprint(directory/'report.json')
    producer=json.loads((a.probe_run/'report.json').read_text());analysis=json.loads((a.probe_analysis/'report.json').read_text())
    if not producer['computational_completed']or not analysis['spectra']or analysis['native_report']!=fingerprint(a.probe_run/'report.json'):raise ValueError('Unqualified numerical output probe')
    for row in producer['output_files']:
        q=a.probe_run/row['path']
        if fingerprint(q)!=row['receipt']:raise ValueError('Changed multiband producer output')
        checks[str(q)]=row['receipt']
    plan_path=ROOT/'config/experiments/rya1371_openmulti3d_oxygen_accurate_profile.json';plan=json.loads(plan_path.read_text());original_windows=plan['output_grid']['windows'];checks[str(plan_path)]=fingerprint(plan_path)
    old_windows=[catalog_interval_vac_nm(lo/10.,hi/10.,'air').tolist()for lo,hi in original_windows]
    name=producer['plan']['name']
    if Path(name).name!=name:raise ValueError('Unsafe probe plan name')
    probe_plan=ROOT/'config/experiments'/name;pr=json.loads(probe_plan.read_text())
    if fingerprint(probe_plan)!=producer['plan']:raise ValueError('Probe plan identity mismatch')
    checks[str(probe_plan)]=fingerprint(probe_plan);probe_windows=pr['output_grid']['windows']
    with np.load(a.probe_analysis/'spectra.npz')as d:
        waves=d['wave_vac_A'].copy()/10.;intensity=d['Inu_cgs'].copy();mask=d['air_frame_supported'].copy()
        if not np.isfinite(waves).all()or not np.isfinite(intensity).all()or np.any(intensity<=0)or not np.array_equal(mask,waves>200.):raise ValueError('Invalid multi-band output/coordinate mask')
    expected,_=output_grid_vacuum_A(pr['output_grid']);np.testing.assert_allclose(waves,expected/10.,rtol=0,atol=2e-12);checks[str(a.probe_analysis/'spectra.npz')]=fingerprint(a.probe_analysis/'spectra.npz')
    rows=[]
    for kind,entries in [('instrument',instruments),('mode',modes)]:
        for r in entries:
            mode=kind=='mode';frame=r['wavelength_frame']if mode else r['wavelength_frame_default'];bounds=catalog_interval_vac_nm(float(r['wavelength_min_nm']),float(r['wavelength_max_nm']),frame)
            intersections=interval_intersections(bounds,old_windows);smoke=interval_intersections(bounds,probe_windows)
            rows.append(dict(kind=kind,id=r['mode_id']if mode else r['instrument_id'],instrument_id=r['instrument_id'],catalog_status=r['mode_status']if mode else r['codex_status'],catalog_min_nm=r['wavelength_min_nm'],catalog_max_nm=r['wavelength_max_nm'],catalog_frame=frame,vacuum_min_nm=float(bounds[0]),vacuum_max_nm=float(bounds[1]),bands=r['bands_supported'],resolving_power=r['resolving_power']if mode else r['resolving_power_min']+'..'+r['resolving_power_max'],source_issue_ids=r['source_issue_ids'],source_documentation=r['source_documentation'],source_atomic_O_BB_centres_in_envelope=int(np.count_nonzero((centres>=bounds[0])&(centres<=bounds[1]))),internal_SE_grid_envelope_contains_bounds=bool(bounds[0]>=grid[0]and bounds[1]<=grid[-1]),old_optical_diagnostic_window_intersections=json.dumps(intersections),new_numerical_smoke_window_intersections=json.dumps(smoke),dense_full_band_synthesis='NOT_COMPLETED',independent_full_band_validation='NOT_ESTABLISHED',scope='Nominal catalog envelope only, not actual setup/order/holding coverage; preserve catalog status including rejected instruments. O-only atomic model lacks general molecular/background line completeness.'))
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'instrument_mode_coverage.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    report=dict(generator=fingerprint(__file__),inputs=checks,instruments=len(instruments),modes=len(modes),rows=len(rows),catalog_target_nm=[min(float(r['wavelength_min_nm'])for r in instruments),max(float(r['wavelength_max_nm'])for r in instruments)],internal_SE_grid_nm=[float(grid[0]),float(grid[-1])],old_diagnostic_windows_air_A=original_windows,new_numerical_probe_windows_vac_nm=probe_windows,probe_points=len(waves),probe_air_unsupported_samples=int(np.count_nonzero(~mask)),source_BB_lines=len(atom.lines),source_BF_continua=len(atom.continua),scientific_admission='HOLD',computational_completed=True,requirements=['Target all registered instruments and modes; regenerate as catalog expands.','Keep sparse SE support, numerical forward output, full-band synthesis and independent validation distinct.','Per-product gaps, wavelength/velocity frame, intensity/flux, LSF/R, continuum/tellurics and model/atom/molecular applicability must qualify separately.','Catalog metadata dates/statuses remain inherited, not freshly verified manufacturer capabilities. Rejected entries are not re-enabled.','Current O41 photospheric diagnostic is not universal UV/chromosphere/molecular/3D physics; no silent fallback or production admission.'])
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed coverage input')
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['instruments','modes','catalog_target_nm','internal_SE_grid_nm','probe_points']}))
if __name__=='__main__':main()

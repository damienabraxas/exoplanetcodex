"""Reviewable catalog gap table and sparse UV/IR numerical-output figure."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['coverage','analysis','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();doc=json.loads((a.coverage/'report.json').read_text());checks={str(a.coverage/'report.json'):fingerprint(a.coverage/'report.json'),str(a.coverage/'instrument_mode_coverage.csv'):fingerprint(a.coverage/'instrument_mode_coverage.csv'),str(a.analysis/'spectra.npz'):fingerprint(a.analysis/'spectra.npz'),str(a.analysis/'report.json'):fingerprint(a.analysis/'report.json')}
    if not doc['computational_completed']or checks[str(a.analysis/'spectra.npz')]!=doc['inputs'].get(str(a.analysis/'spectra.npz')):raise ValueError('Coverage/probe identity mismatch')
    with(a.coverage/'instrument_mode_coverage.csv').open(newline='')as f:rows=list(csv.DictReader(f))
    if len(rows)!=doc['rows']:raise ValueError('Incomplete catalog coverage matrix')
    with np.load(a.analysis/'spectra.npz')as d:wave=d['wave_vac_A'].copy()/10.;normal=d['Inu_cgs']/d['continuum_Inu_cgs']
    a.output.mkdir(parents=True,exist_ok=False);text=['# Oxygen engine: instrument coverage target and current gaps','',f"Catalog target: {doc['instruments']} instruments/reference atlases and {doc['modes']} modes; nominal {doc['catalog_target_nm'][0]:g}--{doc['catalog_target_nm'][1]:g} nm.",'','Sparse atomic SE support and numerical output probes do not certify dense full-band spectra. Current model is O-only photospheric disk-centre intensity; molecular/background line completeness, full flux, instrumental response and observed validation remain owed. Catalog status is inherited; rejected entries are not re-enabled. Actual setup/order/holding gaps remain authoritative.','']
    for kind in ['instrument','mode']:
        text+=['## '+kind.capitalize()+'s','','| ID | Catalog status | Nominal range (nm), catalog frame | Atomic O centres | Old diagnostic overlaps | New numerical probe overlaps | Full-band independently validated |','|---|---|---|---:|---:|---:|---|']
        for r in rows:
            if r['kind']!=kind:continue
            text.append(f"| {r['id']} | {r['catalog_status']} | {r['catalog_min_nm']}--{r['catalog_max_nm']}, {r['catalog_frame']} | {r['source_atomic_O_BB_centres_in_envelope']} | {len(json.loads(r['old_optical_diagnostic_window_intersections']))} | {len(json.loads(r['new_numerical_smoke_window_intersections']))} | Not established |")
        text+=['']
    (a.output/'coverage.md').write_text('\n'.join(text)+'\n')
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    windows=doc['new_numerical_probe_windows_vac_nm'];fig,axes=plt.subplots(len(windows),1,figsize=(10,12),constrained_layout=True)
    for ax,(lo,hi)in zip(axes,windows):
        mask=(wave>=lo-1e-12)&(wave<=hi+1e-12)
        if np.count_nonzero(mask)<2 or not np.isfinite(normal[mask]).all():raise ValueError('Missing numerical probe window')
        ax.plot(wave[mask],normal[mask]);ax.set_xlabel('VACUUM wavelength / nm');ax.set_ylabel('I / continuum');ax.set_title(f'{lo:g}--{hi:g} nm; numerical output only');ax.ticklabel_format(axis='x',style='plain',useOffset=False)
    fig.suptitle('O41/MARCS5750g4.5A8.70; full-BF/BB theoretical disk-centre output\nSparse UV/NIR/IR and catalog-edge probes; no observed/full-band/instrument validation');fig.savefig(a.output/'uv_ir_numerical_probes.png',dpi=150);plt.close(fig)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed coverage-view input')
    (a.output/'report.json').write_text(json.dumps(dict(generator=fingerprint(__file__),inputs=checks,scientific_admission='HOLD',computational_completed=True),indent=2)+'\n')
    print('Created complete instrument/mode gap table and sparse UV/IR probe figure')
if __name__=='__main__':main()

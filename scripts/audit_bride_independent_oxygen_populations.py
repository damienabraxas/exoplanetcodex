"""All-level audit of qualified Bride columns against an exact held MULTI record.

Departure comparisons are diagnostic: source producer LTE/opacity/rate arrays are
not supplied by the record. Never infer a physical pass from a small b difference.
"""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def compare_departures(populations,nstar,reference,tau):
    """Validate distinct arrays and preserve absolute/relative and normalization effects."""
    n,s,b=[np.asarray(x,float)for x in [populations,nstar,reference]];t=np.asarray(tau,float)
    if n.shape!=s.shape or n.shape!=b.shape or n.ndim!=2 or t.shape!=(n.shape[1],):raise ValueError('Wrong matched level/depth layout')
    if any(not np.isfinite(x).all()or np.any(x<=0)for x in [n,s,b,t]):raise ValueError('Invalid positive population/reference array')
    current=n/s;foreign=b*s;scale=n.sum(axis=0)/foreign.sum(axis=0);mapped=foreign*scale
    return dict(local_b=current,reference_b=b,delta_b=current-b,relative_b=current/b-1,reference_mapped_populations=mapped,reference_mapping_scale=scale,logtau500=np.log10(t))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--source',type=Path,default=ROOT/'data/results/rya1281/solar_full_cap_v1/rya1281-solar-full-cap-v1/baseline');p.add_argument('--reference',type=Path);a=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_rate_audit import process_flow
    refpath=a.reference or require_codex_path('experiments.bride_oxygen_reference');ref=json.loads(refpath.read_text());r=json.loads((a.source/'report.json').read_text());checks=[]
    def check(path,receipt):
        if fingerprint(path)!=receipt:raise ValueError('Changed source '+str(path))
        checks.append((path,receipt))
    check(refpath,fingerprint(refpath));check(a.source/'report.json',fingerprint(a.source/'report.json'));check(a.source/'populations.npz',r['coupled']['populations']);check(a.source/'independent_rates.npz',r['independent']['rates']);atompath=require_codex_path('experiments.bride_oxygen_atom');check(atompath,r['atom']);atom=read_multi_atom(atompath)
    if ref['atom']!=r['atom']or ref['reference']['abundance']!=atom.abundance:raise ValueError('Different source atom/abundance')
    if r['case_parameters']!=dict(name='baseline',teff=5750.,logg=4.5,native=True)or ref['atmosphere_file']not in r['source_nodes']:raise ValueError('Not the exact native benchmark node')
    for key in ['height_m','temperature_K','electron_pressure_Pa']:
        np.testing.assert_array_equal(r['state'][key],ref['atmosphere'][key])
    if max(abs(np.log10(ref['reference']['tau'])-np.asarray(ref['atmosphere']['logtau_500'])))>2e-8:raise ValueError('Wrong optical-depth identity')
    with np.load(a.source/'populations.npz')as d:arrays=compare_departures(d['populations'],d['nstar'],ref['reference']['departures'],ref['reference']['tau']);n=d['populations'].copy();s=d['nstar'].copy()
    with np.load(a.source/'independent_rates.npz')as d:net,gross=process_flow(d['collision_rates_s_inv']+d['radiative_rates_s_inv'],n)
    rate=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)));conservation=float(np.max(abs(n.sum(axis=0)/s.sum(axis=0)-1)))
    if rate>=r['controls']['rtol']or conservation>1e-10:raise ValueError('Unqualified local source solution')
    mask=(arrays['logtau500']>=-3)&(arrays['logtau500']<=0);summary=[];rows=[]
    for i,level in enumerate(atom.levels):
        for k,t in enumerate(arrays['logtau500']):rows.append(dict(level=i+1,label=level.label,depth=k,logtau500=t,local_b=arrays['local_b'][i,k],reference_b=arrays['reference_b'][i,k],delta_b=arrays['delta_b'][i,k],relative_b=arrays['relative_b'][i,k],local_population_fraction=n[i,k]/n[:,k].sum()))
        summary.append(dict(level=i+1,label=level.label,all_depth_max_absolute_b_delta=float(np.max(abs(arrays['delta_b'][i]))),diagnostic_max_absolute_b_delta=float(np.max(abs(arrays['delta_b'][i,mask]))),diagnostic_max_relative_b_delta=float(np.max(abs(arrays['relative_b'][i,mask])))))
    for path,receipt in checks:
        if fingerprint(path)!=receipt:raise ValueError('Source changed during population audit')
    a.output.mkdir(parents=True,exist_ok=False)
    for name,data in [('all_levels.csv',rows),('summary.csv',summary)]:
        with(a.output/name).open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(data[0]),lineterminator='\n');w.writeheader();w.writerows(data)
    np.savez_compressed(a.output/'population_comparison.npz',**arrays)
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for ax,i in zip(axes.flat,[5,7,8,9]):
        ax.plot(arrays['logtau500'],arrays['local_b'][i],label='Bride, fresh physical balance');ax.plot(arrays['logtau500'],arrays['reference_b'][i],'--',label='Held independent MULTI');ax.axvspan(-3,0,color='gray',alpha=.1);ax.set_xlabel('log tau500');ax.set_ylabel('b');ax.set_title(atom.levels[i].label);ax.legend(fontsize=8)
    fig.suptitle('Exact source node/atom; departure agreement does not establish producer physics parity');fig.tight_layout();fig.savefig(a.output/'triplet_departures.png',dpi=150);plt.close(fig)
    doc=dict(schema='codex.bride_independent_oxygen_populations/1',scientific_admission='HOLD',generator=fingerprint(__file__),source_report=fingerprint(a.source/'report.json'),reference=fingerprint(refpath),atom=fingerprint(atompath),source_state_arrays_exact=['height_m','temperature_K','electron_pressure_Pa'],independent_local_rate_residual=rate,conservation=conservation,diagnostic_logtau500=[-3,0],triplet_levels_1based=[6,8,9,10],triplet_diagnostic_max_abs_b_delta=float(np.max(abs(arrays['delta_b'][[5,7,8,9]][:,mask]))),all_levels_all_depths_max_abs_b_delta=float(np.max(abs(arrays['delta_b']))),foreign_mapping_scale_range=[float(arrays['reference_mapping_scale'].min()),float(arrays['reference_mapping_scale'].max())],limitations=['Diagnostic [-3,0] is the original comparison interval, not a measured contribution function.','Foreign departure coefficients do not supply native LTE populations, actual opacity, rates, J or exact executable identity.','Foreign b times local LTE needs explicitly recorded number-normalization before any transfer replay; no SE-balance claim for that mapping.','No tolerance inferred from observed agreement and no physical model admission.'])
    (a.output/'report.json').write_text(json.dumps(doc,indent=2,allow_nan=False)+'\n');print(json.dumps(doc))

if __name__=='__main__':main()

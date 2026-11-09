"""Same-physics Ng-on/off baseline audit; no stellar-parameter inference."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from pipeline.bride_reference import fingerprint


def load(directory):
    r=json.loads((directory/'report.json').read_text())
    if r['case']!='baseline'or r.get('dry_run'):raise ValueError('Require completed baseline')
    for name,receipt in [('flux.csv',r['flux']),('populations.npz',r['coupled']['populations']),('independent_rates.npz',r['independent']['rates']),('history.jsonl',r['coupled']['history'])]:
        if fingerprint(directory/name)!=receipt:raise ValueError('Changed baseline artifact '+name)
    with np.load(directory/'populations.npz')as d:pop={k:d[k].copy()for k in d.files}
    with np.load(directory/'independent_rates.npz')as d:c=d['collision_rates_s_inv'].copy();off=c+d['radiative_rates_s_inv']
    i=np.arange(len(off))
    if np.any(off<0)or np.any(off[i,i]!=0):raise ValueError('Require nonnegative physical offdiagonal rates')
    incoming=np.einsum('ijd,jd->id',off,pop['populations']);outgoing=off.sum(axis=0)*pop['populations'];gross=incoming+outgoing
    if np.any(gross<=0):raise ValueError('Disconnected physical rate equation')
    residual=float(np.max(abs(incoming-outgoing)/gross));conservation=float(np.max(abs(pop['populations'].sum(axis=0)/pop['nstar'].sum(axis=0)-1)))
    if residual>=r['controls']['rtol']or conservation>1e-10:raise ValueError('Physical replay failed')
    np.testing.assert_allclose(residual,r['independent']['physical_rate_residual'],rtol=0,atol=4096*np.finfo(float).eps)
    history=[json.loads(line)for line in(directory/'history.jsonl').read_text().splitlines()]
    flux=np.genfromtxt(directory/'flux.csv',delimiter=',',names=True)
    return r,pop,c,history,flux,dict(residual=residual,conservation=conservation)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--on',type=Path,required=True);p.add_argument('--off',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    on,op,oc,oh,of,orate=load(a.on);off,fp,fc,fh,ff,frate=load(a.off)
    for key in ['state','case_parameters','source_nodes','cap_length_m','atom','context_helper','state_helper','backend','LW_constants','source_collision_contract','iteration_plugin','warm_report','warm_populations','initialization_mode']:
        if on.get(key)!=off.get(key):raise ValueError('Ng control changes physical input '+key)
    expected=dict(on['controls']);expected['ng']=False
    if not on['controls']['ng']or off['controls']!=expected:raise ValueError('Expected Ng-only control difference')
    for key in ['nstar','initial_populations']:np.testing.assert_array_equal(op[key],fp[key])
    np.testing.assert_array_equal(oc,fc)
    np.testing.assert_array_equal(of['wavelength_air_A'],ff['wavelength_air_A'])
    df=float(np.max(abs(ff['normalized_flux']-of['normalized_flux'])));dc=float(np.max(abs(ff['continuum_Fnu_W_m2_Hz']/of['continuum_Fnu_W_m2_Hz']-1)))
    plan=ROOT/'config/experiments/rya1281_solar_matched_analysis.json';screens=json.loads(plan.read_text())
    if df>screens['endpoint_max_normalized_flux_difference']or dc>screens['endpoint_max_relative_continuum_difference']:raise ValueError('Existing endpoint parity screen failed')
    rows=[];summary=[]
    for label,r,h,physical in [('Ng_on',on,oh,orate),('Ng_off',off,fh,frate)]:
        increases=[]
        for j,row in enumerate(h):
            rows.append(dict(control=label,iteration=row['iteration'],population_change=row['population_change'],old_radiation_rate_residual=row['physical_rate_residual'],J_change=row['J_change'],cap_change=row['cap_boundary_change']))
            if j and row['population_change']>h[j-1]['population_change']:increases.append(row['iteration'])
        summary.append(dict(control=label,iterations=len(h),seconds=r['telemetry']['seconds'],fresh_physical=physical,population_change_increase_iterations=increases))
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'history.csv').open('x',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(16,4.5))
    for label,h in [('Ng on',oh),('Ng off',fh)]:
        iterations=[r['iteration']for r in h]
        axes[0].semilogy(iterations,np.ma.masked_less_equal([r['population_change']for r in h],0),label=label)
        axes[1].semilogy(iterations,np.ma.masked_less_equal([max(r['J_change'],r['cap_boundary_change'])for r in h],0),label=label)
    for ax,title in zip(axes[:2],['Population change','Maximum J / cap change']):
        ax.axhline(on['controls']['rtol'],color='gray',ls='--',label='Original joint threshold');ax.set_title(title);ax.set_xlabel('Outer iteration');ax.legend()
    axes[2].plot(ff['wavelength_air_A'],ff['normalized_flux']-of['normalized_flux']);axes[2].set_title('Ng off minus Ng on: normalized flux');axes[2].set_xlabel('AIR wavelength (A)')
    fig.suptitle('Native baseline numerical control; source physics and abundance unchanged\nBoth runs qualified; no parameter or physical-accuracy inference');fig.tight_layout();fig.savefig(a.output/'ng_baseline.png',dpi=140);plt.close(fig)
    result=dict(schema='codex.bride_ng_baseline/1',scientific_admission='HOLD',generator=fingerprint(__file__),source_reports=[fingerprint(a.on/'report.json'),fingerprint(a.off/'report.json')],screen_plan=fingerprint(plan),summaries=summary,max_normalized_flux_difference=df,max_relative_continuum_difference=dc,identities=['Physical inputs, original starting populations, LTE populations, collision arrays and wavelengths exact; only native Ng flag differs.'],limitations=['Single native atmosphere and two bounded runs, not a general performance or convergence proof.','History rate residual uses old radiation; only independently replayed final residuals qualify the endpoints.','Stellar T/g/xi/abundance inference and full four-case factorial analysis remain separate gates.'])
    (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result))

if __name__=='__main__':main()

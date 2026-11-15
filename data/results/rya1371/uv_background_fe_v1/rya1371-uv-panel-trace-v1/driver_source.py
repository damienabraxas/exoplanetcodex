"""Trace ranked continuum quadrature panels against source alpha, J and HOP opacity."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint
    from pipeline.multi_atom import read_multi_atom,to_lightweaver
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.hop_packet import compute_packet
    from pipeline.adaptive_continuum_grid import nested_indices
    import lightweaver.constants as c
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['replay','audit','native-analysis','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--top-panels',type=int,default=5);a=p.parse_args()
    if not 1<=a.top_panels<=20:raise ValueError('Require bounded panel trace')
    checks={str(Path(__file__)):fingerprint(__file__)}
    for directory,names in [(a.replay,['report.json','rates.npz']),(a.audit,['report.json','adaptive_grid.npz']),(a.native_analysis,['report.json','arrays.npz'])]:
        for name in names:checks[str(directory/name)]=fingerprint(directory/name)
    rd=json.loads((a.replay/'report.json').read_text());ad=json.loads((a.audit/'report.json').read_text())
    if not rd['computational_completed']or not ad['computational_completed']or ad['inputs'].get(str(a.replay/'rates.npz'))!=fingerprint(a.replay/'rates.npz'):raise ValueError('Unmatched UV panel audit')
    for name in ['report.json','arrays.npz']:
        if rd['inputs'].get(str(a.native_analysis/name))!=fingerprint(a.native_analysis/name):raise ValueError('Wrong native state')
    with np.load(a.replay/'rates.npz')as d:wave=d['wave_vac_nm'].copy();J=d['Jnu'].copy();star=d['nstar'].copy();pop=d['populations'].copy()
    with np.load(a.audit/'adaptive_grid.npz')as d:
        base=d['base_wave_vac_nm'].copy();scores=d['panel_scores'].copy();np.testing.assert_array_equal(wave,d['fine_wave_vac_nm']);np.testing.assert_array_equal(pop,d['populations']);np.testing.assert_array_equal(star,d['nstar'])
    nested_indices(base,wave);chosen=np.argsort(-scores,kind='stable')[:a.top_panels]
    with np.load(a.native_analysis/'arrays.npz')as d:
        state=d['consumed_state_cgs'][::-1].astype(float);T=state[:,2];ne=state[:,0]*1e6;rho=state[:,1]*1e3;vmic=state[:,7]*1e3
    atom_path=require_codex_path('experiments.bride_oxygen_atom');atom=read_multi_atom(atom_path);model=to_lightweaver(atom);checks[str(atom_path)]=fingerprint(atom_path)
    ref_path=require_codex_path('experiments.bride_oxygen_reference');ref=json.loads(ref_path.read_text());atmos=require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name'];checks[str(ref_path)]=fingerprint(ref_path);checks[str(atmos)]=fingerprint(atmos)
    if checks[str(atmos)]!=ref['atmosphere_file']:raise ValueError('Changed source atmosphere')
    structure=read_marcs_structure(atmos);W=1.0000000018626451
    # Read the ray sum from the qualified replay receipt rather than renormalize.
    W=float(np.sum(rd['context_audit']['angles']['weight']))
    winners={};candidates=[]
    for transition in model.continua:
        mask=(wave>=transition.minLambda)&(wave<=transition.lambdaEdge);local=wave[mask];alpha=transition.alpha(local);field=J[mask];i=transition.i;j=transition.j
        nu=c.CLight/(local[:,None]*1e-9);pref=4*np.pi/c.HPlanck*alpha[:,None]/local[:,None];beta=np.exp(-c.HPlanck*nu/(c.KBoltzmann*T));bare=2*c.HPlanck*nu**3/c.CLight**2
        integrands={'up':pref*field,'down':pref*star[i][None,:]/star[j][None,:]*beta*(bare*W+field)}
        for direction,F in integrands.items():
            total=np.sum((F[1:]+F[:-1])*.5*np.diff(local)[:,None],axis=0)
            for panel in chosen:
                lo,hi=base[panel:panel+2]
                if lo<local[0]or hi>local[-1]:continue
                endpoints=np.searchsorted(local,[lo,hi]);np.testing.assert_array_equal(local[endpoints],[lo,hi]);start,end=endpoints;xx=local[start:end+1];ff=F[start:end+1]
                if (len(xx)-1)%2:raise ValueError('Non-nested trace panel')
                small=(ff[1:]+ff[:-1])*.5*np.diff(xx)[:,None];coarse=(ff[::2][1:]+ff[::2][:-1])*.5*(xx[2::2]-xx[:-2:2])[:,None]
                by_depth=np.sum(abs(small[::2]+small[1::2]-coarse),axis=0)/total;z=int(np.argmax(by_depth));score=float(by_depth[z])
                candidates.append(dict(panel=int(panel),lower_1based=i+1,upper_1based=j+1,direction=direction,depth_top_to_bottom=z,indicator=score))
                if panel not in winners or score>winners[panel]['indicator']:
                    winners[panel]=dict(panel=int(panel),lower_1based=i+1,upper_1based=j+1,direction=direction,depth_top_to_bottom=z,indicator=score,local=xx.copy(),alpha=alpha[start:end+1].copy(),J=field[start:end+1,z].copy(),integrand=ff[:,z].copy(),knots=np.asarray(transition.wavelengthGrid),knot_alpha=np.asarray(transition.alphaGrid))
    a.output.mkdir(parents=True,exist_ok=False);panels=[];sample_rows=[];packet_receipts=[]
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(a.top_panels,3,figsize=(13,2.6*a.top_panels),layout='constrained',squeeze=False)
    for row,panel in enumerate(chosen):
        w=winners[panel];z=w['depth_top_to_bottom'];local=w['local'];lo,hi=base[panel:panel+2]
        if not np.isclose(w['indicator'],scores[panel],rtol=2e-12,atol=0):raise ValueError('Panel score does not reproduce source indicator')
        packet,meta=compute_packet(T,ne,rho,local,structure['abundances_dex'],atom,vmic_m_s=float(vmic[0]));packet_receipts.append(meta)
        for key,fp in meta.get('inputs',{}).items():checks[key]=fp
        knot_mask=(w['knots']>=lo)&(w['knots']<=hi);knots=w['knots'][knot_mask];kalpha=w['knot_alpha'][knot_mask]
        centres=[float(l.lambda0)for l in model.lines if lo<=l.lambda0<=hi]
        r={k:v for k,v in w.items()if k not in ['local','alpha','J','integrand','knots','knot_alpha']}
        r.update(lo_vac_nm=float(lo),hi_vac_nm=float(hi),sample_points=len(local),source_knots=knots.tolist(),O_BB_centres_vac_nm=centres,temperature_K=float(T[z]),J_max_over_min=float(np.max(w['J'])/np.min(w['J'])),alpha_max_over_min=float(np.max(w['alpha'])/np.min(w['alpha'])),background_chi_max_over_min=float(np.max(packet['chi_m_inv'][:,z])/np.min(packet['chi_m_inv'][:,z])))
        panels.append(r)
        for k,lam in enumerate(local):sample_rows.append(dict(panel=int(panel),vacuum_nm=float(lam),alpha_m2=float(w['alpha'][k]),Jnu_SI=float(w['J'][k]),rate_integrand=float(w['integrand'][k]),background_chi_m_inv=float(packet['chi_m_inv'][k,z]),background_eta_SI=float(packet['eta_SI'][k,z]),background_sca_m_inv=float(packet['scat_m_inv'][k,z]),is_source_knot=bool(lam in knots)))
        axes[row,0].plot(local,w['alpha'],'o-',ms=2);axes[row,0].scatter(knots,kalpha,marker='x',c='red',label='source knots');axes[row,0].set_ylabel('alpha [m²]');axes[row,0].legend(fontsize=7)
        axes[row,1].plot(local,w['J'],'o-',ms=2);axes[row,1].set_ylabel('Jnu [SI]')
        axes[row,2].plot(local,packet['chi_m_inv'][:,z],'o-',ms=2);axes[row,2].set_ylabel('Background chi [m⁻¹]')
        for ax in axes[row]:ax.set_xlabel('Vacuum nm');ax.ticklabel_format(axis='x',style='plain',useOffset=False)
        axes[row,0].set_title(f"BF {w['lower_1based']}→{w['upper_1based']}, {w['direction']}, depth {z}")
    fig.savefig(a.output/'uv_panels.png',dpi=140);plt.close(fig)
    for name,rows in [('samples.csv',sample_rows),('candidates.csv',candidates)]:
        with(a.output/name).open('w',newline='')as f:writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed UV trace input')
    report=dict(generator=fingerprint(__file__),inputs=checks,panels=panels,HOP_packet_receipts=packet_receipts,computational_completed=True,scientific_admission='HOLD',scope='Source-table, sampled-field and independently evaluated original HOP background structure at the highest estimated continuum panels. No RT/SE/refinement/physics changes; packet chi includes scattering, not an added extra scattering term. Indicators reproduced before interpretation; correlations are not unique causal attribution. No solar/abundance/Fe validation.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(panels))


if __name__=='__main__':main()

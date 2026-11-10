"""Compare receipted TS component budgets and local LTE at matched spectral samples."""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure,make_reference_context
    from pipeline.bride_ts_budget import read_ts_opacity,read_ts_components
    from pipeline.bride_marcs_opacity import read_marcs_opacity
    from pipeline.bride_continuum_budget import continuum_budget
    from pipeline.multi_atom import read_multi_atom
    from pipeline.wavelength_util import air_to_vac
    from lightweaver.atomic_table import PeriodicTable
    report=json.loads(args.run.read_text())
    if not report['reporting_preserves_serialized_coefficients']:raise ValueError('Unverified reporting build')
    if fingerprint(ROOT/'pipeline/bride_ts_budget.py')!=report['parser']:raise ValueError('Changed output parser')
    folder=require_codex_path('experiments.bride_marcs_opacity')
    intake=json.loads((folder/'manifest.json').read_text());opa_path=folder/intake['source']['name']
    model_path=require_codex_path('grids.marcs_standard')/intake['atmosphere']['name']
    if fingerprint(model_path)!=report['atmosphere'] or fingerprint(opa_path)!=intake['source']:
        raise ValueError('Changed native inputs')
    structure=read_marcs_structure(model_path);opa=read_marcs_opacity(opa_path)
    source=read_multi_atom(require_codex_path('experiments.bride_oxygen_atom'))
    source.abundance=structure['abundances_dex'][7]
    rows=[];budgets=[];cases=[]
    for case in report['cases']:
        if not case['name'].endswith('_reporting'):continue
        directory=args.run.parent/case['name']
        for name,key in [('continuum.opa','output'),('detail.txt','detail')]:
            if fingerprint(directory/name)!=case[key]:raise ValueError('Changed TS output')
        native=read_ts_opacity(directory/'continuum.opa')
        component=read_ts_components(directory/'detail.txt',native['ndep'],native['nwave'])
        # Babsma prints the standard wavelength first, followed by every output
        # sample, once per native depth. Require the entire sequence, not nearest
        # rounded labels which could silently match a different wavelength.
        labels=np.rint(np.r_[5000.,native['wavelength_air_A']])
        np.testing.assert_array_equal(component['records'][:,:,0],np.broadcast_to(labels,(native['ndep'],len(labels))))
        i=1;wave_air=native['wavelength_air_A'][i];wave_vac=float(air_to_vac(np.array([wave_air]))[0])
        ctx=make_reference_context(structure,source,rays=2,qscale=1,lte=True,
            continuum_only=True,extra_wavelengths_vac_nm=np.array([wave_vac*.1]))
        before=np.array(ctx.detailedAtoms[0].n).copy();ctx.depthData.fill=True
        ctx.formal_sol_gamma_matrices(lambdaIterate=True)
        np.testing.assert_array_equal(ctx.detailedAtoms[0].n,before)
        index=int(np.searchsorted(np.asarray(ctx.spect.wavelength),wave_vac*.1))
        assert np.asarray(ctx.spect.wavelength)[index]==wave_vac*.1
        chi=np.asarray(ctx.depthData.chi)[index,0,0]
        local_abs=(chi-np.asarray(ctx.background.sca)[index])/structure['rho_kg_m3']
        bf,_=continuum_budget(ctx);local_fe=bf['Fe(Z=26)'][index]/structure['rho_kg_m3']
        fe=ctx.kwargs['eqPops'].atomicPops[PeriodicTable[26]]
        neutral=np.array([l.stage==0 for l in fe.model.levels])
        local_fe_per_kg=fe.n[neutral].sum(axis=0)/structure['rho_kg_m3']
        state=component['states'];T,Pe,rho,U,fe_pressure,feii_pressure=state.T
        # Same rhokt convention as the inspected native DETABS source, retaining
        # its constant vintage. dyn/cm2 / (erg/K*K) is cm^-3; -> m^-3 / kg/m3.
        ts_fe_per_kg=fe_pressure/(1.38066e-16*T)*1e6/(rho*1000)
        fe_index=component['labels'].index('FeI')
        fractions=component['fractions'][:,i+1,:]
        ts_abs=native['absorption_m2_kg'][i]
        intermediate_abs=component['records'][:,i+1,2]*.1
        ts_fe=intermediate_abs*fractions[:,fe_index]
        if component['post'] is None:raise ValueError('Post-DETABS bookkeeping missing')
        post=component['post'][:,i+1]
        np.testing.assert_allclose(intermediate_abs+post[:,1]*.1+post[:,2]*.1,
            post[:,3]*.1,rtol=2*np.finfo(np.float32).eps,atol=0)
        np.testing.assert_allclose(post[:,3]*.1,ts_abs,rtol=6e-5,atol=0)
        coordinate=wave_vac if case['name'].startswith('OPA_is_vacuum')else wave_air
        marcs=np.array([np.interp(coordinate,opa['wavelength_source_A'],opa['absorption_m2_kg'][:,d]) for d in range(native['ndep'])])
        for d in range(native['ndep']):
            rows.append(dict(case=case['name'],depth_1based=d+1,logtau500=structure['logtau_500'][d],
                wavelength_TS_air_A=wave_air,wavelength_conditional_OPA_A=coordinate,
                ts_absorption_m2_kg=ts_abs[d],marcs_absorption_m2_kg=marcs[d],local_absorption_m2_kg=local_abs[d],
                TS_over_MARCS_absorption=ts_abs[d]/marcs[d],TS_FeI_fraction=ts_fe[d]/ts_abs[d],
                DETABS_FeI_fraction=fractions[d,fe_index],
                post_DETABS_H_bf_m2_kg=post[d,1]*.1,post_DETABS_NH_m2_kg=post[d,2]*.1,
                TS_FeI_m2_kg=ts_fe[d],local_FeI_m2_kg=local_fe[d],
                TS_over_local_FeI_opacity=ts_fe[d]/local_fe[d],
                TS_over_local_neutral_Fe_per_kg=ts_fe_per_kg[d]/local_fe_per_kg[d],
                TS_over_local_net_opacity_per_neutral_Fe=(ts_fe[d]/ts_fe_per_kg[d])/(local_fe[d]/local_fe_per_kg[d])))
            for j,label in enumerate(component['labels']):
                kind='absorption' if j<len(component['labels'])-4 else 'scattering'
                total=intermediate_abs[d] if kind=='absorption'else component['records'][d,i+1,3]*.1
                budgets.append(dict(case=case['name'],depth_1based=d+1,component=label,kind=kind,
                    native_fraction=fractions[d,j],coefficient_m2_kg=total*fractions[d,j]))
        cases.append(dict(name=case['name'],Fe_partition_range=[float(U.min()),float(U.max())],
            local_Fe_model_levels=len(fe.model.levels),local_Fe_continua=len(fe.model.continua),
            top=rows[-native['ndep']]))
    args.output.mkdir(parents=True,exist_ok=False)
    result=dict(schema='codex.bride_ts_budget_analysis/1',scientific_admission='HOLD',
        generator=fingerprint(__file__),parser=fingerprint(ROOT/'pipeline/bride_ts_budget.py'),
        context_builder=fingerprint(ROOT/'pipeline/bride_reference.py'),run=fingerprint(args.run),cases=cases,
        limitations='Conditional OPA frame, distinct EOS/density and constant vintages retained. Net opacity per neutral atom includes level weights and cross-sections, not a measured single-level cross-section. No opacity substitution or O NLTE run.')
    for name,table in [('comparison.csv',rows),('native_components.csv',budgets)]:
        with (args.output/name).open('x')as stream:
            writer=csv.DictWriter(stream,fieldnames=list(table[0]),lineterminator='\n');writer.writeheader();writer.writerows(table)
        result[name]=fingerprint(args.output/name)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(cases,indent=2))


if __name__=='__main__':main()

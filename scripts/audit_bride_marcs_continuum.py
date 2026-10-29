"""Conditional LTE continuum comparison; retains unverified OPA wavelength frame.

Both frame hypotheses are explicit diagnostic interventions. No adoption, SE
iteration, opacity replacement or claim to native MULTI producer parity.
"""
import argparse
import copy
import csv
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure,make_reference_context
    from pipeline.bride_marcs_opacity import read_marcs_opacity
    from pipeline.multi_atom import read_multi_atom
    from pipeline.wavelength_util import air_to_vac
    from pipeline.bride_continuum_budget import continuum_budget
    import lightweaver
    import lightweaver.LwCompiled as backend
    from lightweaver import rh_atoms
    directory=require_codex_path('experiments.bride_marcs_opacity')
    intake=json.loads((directory/'manifest.json').read_text())
    path=directory/intake['source']['name']
    if fingerprint(path)!=intake['source']:raise ValueError('Changed opacity source')
    model_path=require_codex_path('grids.marcs_standard')/intake['atmosphere']['name']
    if fingerprint(model_path)!=intake['atmosphere']:raise ValueError('Changed atmosphere')
    structure=read_marcs_structure(model_path);opa=read_marcs_opacity(path)
    for key in ('temperature_K','electron_pressure_Pa','rho_kg_m3'):
        np.testing.assert_array_equal(opa[key],structure[key])
    atom_path=require_codex_path('experiments.bride_oxygen_atom')
    if fingerprint(atom_path)['sha256']!='d058b75bdec80c5317523a6d7f185a64ca59144d3d494a17eeabee5f10d04b57':
        raise ValueError('Changed held oxygen atom')
    source=read_multi_atom(atom_path)
    mask=(opa['wavelength_source_A']>=1000)&(opa['wavelength_source_A']<=25000)
    source_wave=opa['wavelength_source_A'][mask]
    args.output.mkdir(parents=True,exist_ok=False)
    report=dict(schema='codex.bride_marcs_continuum_audit/1',scientific_admission='HOLD',
        generator=fingerprint(__file__),context_builder=fingerprint(ROOT/'pipeline/bride_reference.py'),
        parser=fingerprint(ROOT/'pipeline/bride_marcs_opacity.py'),intake=fingerprint(directory/'manifest.json'),
        source=intake['source'],atmosphere=intake['atmosphere'],atom=fingerprint(atom_path),
        passive_atom_factories=fingerprint(rh_atoms.__file__),
        backend_library=fingerprint(backend.__file__),
        lightweaver_version=lightweaver.__version__,threads=1,SE_updates=False,
        wavelength_medium='UNVERIFIED; two conditional hypotheses, neither adopted',
        wavelength_source_range_A=[float(source_wave[0]),float(source_wave[-1])],
        ndepth=structure['ndep'],nwave=len(source_wave),cases=[],
        limitations='Independent MARCS LTE continuum, not producer ABSDAT. O bound-free explicitly included. Distinct EOS/partition functions/passive atom data remain. No fitting, abundance adoption, opacity substitution or observed-spectrum change.')
    rows=[];budget_rows=[]
    for frame,wave in [('source_is_vacuum',source_wave*.1),
                       ('source_is_air',air_to_vac(source_wave)*.1)]:
        for composition,abundance in [('matched_MARCS',structure['abundances_dex'][7]),
                                       ('benchmark_O870',source.abundance)]:
            start=time.perf_counter();model=copy.deepcopy(source);model.abundance=abundance
            ctx=make_reference_context(structure,model,rays=3,qscale=1,lte=True,
                continuum_only=True,extra_wavelengths_vac_nm=wave)
            atom=ctx.detailedAtoms[0];before=np.array(atom.n).copy()
            if not report.get('passive_inventory'):
                report['passive_inventory']=[dict(element=str(a.element),nlevels=len(a.levels),
                    ncontinua=len(a.continua),nlines=len(a.lines),
                    continuum_edges_vac_nm=[float(t.lambdaEdge) for t in a.continua])
                    for a in ctx.kwargs['spect'].radSet.passiveAtoms]
            ctx.depthData.fill=True
            ctx.formal_sol_gamma_matrices(lambdaIterate=True)
            np.testing.assert_array_equal(atom.n,before)
            grid=np.asarray(ctx.spect.wavelength);indices=np.searchsorted(grid,wave)
            np.testing.assert_array_equal(grid[indices],wave)
            full=np.asarray(ctx.depthData.chi)[indices]
            total=full[:,0,0,:].copy()
            np.testing.assert_allclose(full,np.broadcast_to(total[:,None,None,:],full.shape),rtol=1e-12,atol=0)
            scattering=np.asarray(ctx.background.sca)[indices].copy()
            background=np.asarray(ctx.background.chi)[indices].copy()
            bf,sca=continuum_budget(ctx)
            bf={key:values[indices] for key,values in bf.items()}
            sca={key:values[indices] for key,values in sca.items()}
            oxygen=total-background;absorption=total-scattering
            if (not np.isfinite(full).all() or np.any(total<=0) or np.any(scattering<0)
                    or np.any(absorption<0) or np.any(oxygen < -1e-12*total)):
                raise ValueError('Invalid continuum coefficient capture')
            rho=structure['rho_kg_m3'][None,:]
            absorption/=rho;scattering/=rho;oxygen/=rho;background/=rho
            bf={key:values/rho for key,values in bf.items()}
            sca={key:values/rho for key,values in sca.items()}
            # Positive residuals are expected: H-/H/H2 free-free/other basic
            # absorption, and molecular Rayleigh are not reconstructed here.
            other_absorption=absorption-oxygen-sum(bf.values())
            other_scattering=scattering-sum(sca.values())
            if (np.any(other_absorption < -1e-10*absorption)
                    or np.any(other_scattering < -1e-10*scattering)):
                raise ValueError('Reconstructed component budget exceeds actual backend')
            a_ratio=absorption/opa['absorption_m2_kg'][mask]
            s_ratio=scattering/opa['scattering_m2_kg'][mask]
            name=frame+'_'+composition
            array_path=args.output/(name+'.npz')
            np.savez_compressed(array_path,wavelength_source_A=source_wave,
                wavelength_conditional_vac_nm=wave,logtau500=structure['logtau_500'],
                local_absorption_m2_kg=absorption,local_scattering_m2_kg=scattering,
                explicit_oxygen_bf_m2_kg=oxygen,local_background_extinction_m2_kg=background,
                marcs_absorption_m2_kg=opa['absorption_m2_kg'][mask],
                marcs_scattering_m2_kg=opa['scattering_m2_kg'][mask],
                other_absorption_m2_kg=other_absorption,other_scattering_m2_kg=other_scattering,
                **{'passive_bf_'+key:values for key,values in bf.items()},
                **{'scattering_'+key:values for key,values in sca.items()})
            case=dict(name=name,oxygen_abundance_dex=abundance,array=fingerprint(array_path),
                elapsed_seconds=time.perf_counter()-start,
                absorption_ratio_range=[float(a_ratio.min()),float(a_ratio.max())],
                scattering_ratio_range=[float(s_ratio.min()),float(s_ratio.max())])
            report['cases'].append(case)
            for w in (1302.,1355.6,2772.477,5000.,7772.,10000.,16000.,22000.):
                wi=int(np.argmin(abs(source_wave-w)))
                for tau in (-4.9,-3.,-1.,0.):
                    d=int(np.argmin(abs(structure['logtau_500']-tau)))
                    rows.append(dict(case=name,wavelength_source_A=source_wave[wi],
                        conditional_vac_nm=wave[wi],logtau500=structure['logtau_500'][d],
                        local_absorption_m2_kg=absorption[wi,d],marcs_absorption_m2_kg=opa['absorption_m2_kg'][mask][wi,d],
                        absorption_local_over_marcs=a_ratio[wi,d],scattering_local_over_marcs=s_ratio[wi,d],
                        explicit_oxygen_bf_fraction_absorption=oxygen[wi,d]/absorption[wi,d]))
                    for kind,components,total_coefficient in [('passive_bound_free',bf,absorption),
                                                               ('scattering',sca,scattering)]:
                        for component,value in components.items():
                            budget_rows.append(dict(case=name,wavelength_source_A=source_wave[wi],
                                logtau500=structure['logtau_500'][d],kind=kind,component=component,
                                coefficient_m2_kg=value[wi,d],fraction_of_local_total=value[wi,d]/total_coefficient[wi,d]))
    with (args.output/'selected_coefficients.csv').open('x')as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    report['selected_coefficients']=fingerprint(args.output/'selected_coefficients.csv')
    with (args.output/'selected_budgets.csv').open('x')as stream:
        writer=csv.DictWriter(stream,fieldnames=list(budget_rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(budget_rows)
    report['selected_budgets']=fingerprint(args.output/'selected_budgets.csv')
    report['budget_reconstruction']=fingerprint(ROOT/'pipeline/bride_continuum_budget.py')
    source_root=require_codex_path('engines.bride_eval')/'lw_src'/'Source'
    report['inspected_backend_sources']=[fingerprint(source_root/name) for name in
        ('FastBackground.cpp','Background.cpp','LwMiddleLayer.pyx')]
    ts_source=require_codex_path('engines.ispec')/'synthesizer'/'turbospectrum'/'source'
    report['inspected_independent_TS_sources']=[fingerprint(ts_source/name) for name in
        ('babsma.f','detabs.f','hbop.f','inabs.f','readmo.f')]
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report['cases'],indent=2))


if __name__=='__main__':main()

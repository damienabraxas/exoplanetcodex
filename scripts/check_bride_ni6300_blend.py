"""Mac trace-LTE Ni6300 opacity prototype, mixed with fixed NLTE oxygen."""
from pathlib import Path
import argparse
import csv
import json
import platform
import subprocess
import sys
import time
from importlib.metadata import version
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);args=p.parse_args()
    if platform.system()!='Darwin' or version('lightweaver')!='0.16.1':raise RuntimeError('Mac/LW0.16.1 required')
    from config.constants import require_codex_path,SOLAR_ASPLUND2021
    from pipeline.bride_reference import fingerprint,read_marcs_structure,make_reference_context
    from pipeline.bride_parameter_checks import controlled_structure,interpolate_marcs_structure,marcs_filename
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_lte_blend import ni_provider
    from pipeline.bride_observables import flux_quadrature,integrate_flux,window_integral
    from pipeline.wavelength_util import air_to_vac
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    import lightweaver as lw
    from lightweaver.utils import get_data_path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    canonical=ROOT/'data/linelists/canonical_gf.csv'
    with canonical.open() as stream:
        matches=[r for r in csv.DictReader(stream) if r['species']=='Ni I' and abs(float(r['wavelength_air_A'])-6300.342)<.001]
    if len(matches)!=1:raise ValueError('Ambiguous canonical Ni6300')
    row=matches[0]
    if int(row['hfs_n_components'])!=1:raise ValueError('Prototype requires canonical single combined line')
    # Explicit installed-iSpec asset, not a fallback for absent data.ges_v6_atomic_lines.
    ges=require_codex_path('engines.ispec')/'input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv'
    with ges.open() as stream:
        gesrows=[r for r in csv.DictReader(stream,delimiter='\t') if r['element']=='Ni 1' and abs(float(r['wave_A'])-float(row['wavelength_air_A']))<.001]
    if len(gesrows)!=1:raise ValueError('Missing unique GES natural damping source')
    line=dict(wavelength_air_A=float(row['wavelength_air_A']),
        lambda0_vac_nm=float(air_to_vac(np.array([float(row['wavelength_air_A'])]))[0]*.1),
        excitation_eV=float(row['excitation_potential_eV']),loggf=float(row['log_gf']),
        gamma_rad_s_inv=float(gesrows[0]['turbospectrum_rad']))
    if not np.isfinite(list(line.values())).all() or line['gamma_rad_s_inv']<=0:raise ValueError('Invalid Ni line data')
    atom_path=require_codex_path('experiments.bride_oxygen_atom')
    if fingerprint(atom_path)['sha256']!='d058b75bdec80c5317523a6d7f185a64ca59144d3d494a17eeabee5f10d04b57':raise ValueError('Changed O atom')
    source=read_multi_atom(atom_path)
    holding=require_codex_path('data.kurucz2005_residual')
    if fingerprint(holding)['sha256']!='75dd93479b5d3fe63c2e7828f5dbb3ce29f4c672f4f1147e0d30ce41bc8f5adb':raise ValueError('Changed observation')
    receipt_path=ROOT/'data/results/rya1281/mac_parameter_checks_v1/analysis.json'
    reviewed=json.loads(receipt_path.read_text())
    args.output.mkdir(parents=True,exist_ok=False)
    report=dict(schema='codex.bride_ni6300_prototype/1',status='RUNNING',scientific_admission='HOLD',
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        generator=fingerprint(__file__),provider=fingerprint(ROOT/'pipeline/bride_lte_blend.py'),
        context_builder=fingerprint(ROOT/'pipeline/bride_reference.py'),canonical=fingerprint(canonical),
        canonical_row=row,installed_ispec_ges=fingerprint(ges),ges_row=gesrows[0],line=line,
        kurucz_partition=fingerprint(Path(get_data_path())/'pf_Kurucz.input'),
        A_Ni=float(SOLAR_ASPLUND2021['Ni']),Ni_abundance_origin='config.constants.SOLAR_ASPLUND2021; trace opacity only',
        O_atom=fingerprint(atom_path),observation=fingerprint(holding),reviewed_populations=fingerprint(receipt_path),
        limitations=['Natural Voigt only: pressure broadening not implemented',
                    'Single combined canonical Ni line: no isotope splitting or other blends',
                    'Kurucz LTE Ni partition/Saha at fixed MARCS ne/nH; no EOS/atmosphere update',
                    'O populations remain prior RH-convention controls; convention discrepancy unresolved',
                    'No fitted abundance, continuum, width or RV; not a complete O6300 analysis'],cases=[])
    def save():(args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    save();wave=np.arange(6299.,6302.+.001, .002);vac=air_to_vac(wave)*.1
    store=require_codex_path('grids.marcs_standard');grid=require_codex_path('experiments.bride_mac_grid')
    for name,relative in [('native5750','rya1281-v1/baseline'),('solar5772','rya1281-mesh-v1/meshfine_T5772')]:
        origin=grid/relative;m=json.loads((origin/'manifest.json').read_text());case=m['case']
        receipt=next(c for c in reviewed['cases'] if c['name']==case['name'])
        if m['status']!='CONVERGED_EXPERIMENT' or fingerprint(origin/'manifest.json')!=receipt['source_manifest'] or fingerprint(origin/'populations.npz')!=receipt['original_population_file']:
            raise ValueError('Changed/unconverged O source')
        if case.get('interpolate'):
            nodes=[]
            for t in case['bracket_teff']:
                for g in case['bracket_logg']:
                    path=store/marcs_filename(t,g)
                    if fingerprint(path) not in m['atmosphere']['sources']:raise ValueError('Changed interpolation source')
                    nodes.append(((t,g),read_marcs_structure(path)))
            structure,_=interpolate_marcs_structure(nodes,(case['teff'],case['logg']),thermodynamic_scale=case['thermodynamic_scale'])
        else:
            path=store/m['atmosphere']['name']
            if fingerprint(path)!=m['atmosphere']:raise ValueError('Changed atmosphere')
            structure=read_marcs_structure(path)
        structure=controlled_structure(structure,xi=case['xi'],drop_bottom=case['drop_bottom'])
        saved=np.load(origin/'populations.npz');pop=saved['populations']
        np.testing.assert_array_equal(saved['temperature_K'],structure['temperature_K'])
        curves={};iterations={};start=time.perf_counter()
        for mode,no_oxygen,use_ni in [('continuum',True,False),('oxygen_only',False,False),('nickel_only',True,True),('oxygen_plus_nickel',False,True)]:
            ctx=make_reference_context(structure,source,rays=3,qscale=1,fixed_populations=pop,continuum_only=no_oxygen)
            if use_ni:
                kwargs=dict(ctx.kwargs);kwargs['backgroundProvider']=ni_provider(line,report['A_Ni']);ctx=lw.Context(**kwargs)
            for j in range(1000):
                if time.perf_counter()-start>300:raise TimeoutError('Blend case exceeded 300 seconds')
                u=ctx.formal_sol_gamma_matrices()
                if u.dJMax<1e-8:break
            else:raise RuntimeError('Blend radiation did not converge')
            np.testing.assert_array_equal(ctx.detailedAtoms[0].n,pop)
            mu,w=flux_quadrature(5);curves[mode]=integrate_flux(ctx.compute_rays(wavelengths=vac,mus=mu),mu,w)
            if not np.isfinite(curves[mode]).all() or np.any(curves[mode]<=0):raise ValueError('Invalid blend flux')
            iterations[mode]=j+1
        normal={mode:curve/curves['continuum'] for mode,curve in curves.items()}
        wo,fo=read_kurucz2005_residual(holding,wave[0],wave[-1]);observed=np.interp(wave,wo,fo)
        metrics={mode:1000*window_integral(wave,curve,6300.004,6300.604) for mode,curve in normal.items()}
        metrics['observed_blended_window']=1000*window_integral(wo,fo,6300.004,6300.604)
        # A shared physical transfer is generally different from multiplying residuals.
        product=normal['oxygen_only']*normal['nickel_only']
        nonfactor=float(np.max(abs(normal['oxygen_plus_nickel']-product)))
        names=list(normal)
        np.savetxt(args.output/f'{name}_profiles.csv',np.column_stack([wave]+[normal[k] for k in names]+[product,observed]),delimiter=',',
            header='air_A,'+','.join(names)+',product_of_separate_profiles,kitt_peak',comments='')
        report['cases'].append(dict(name=name,source_manifest=fingerprint(origin/'manifest.json'),
            source_populations=fingerprint(origin/'populations.npz'),iterations=iterations,
            elapsed_seconds=time.perf_counter()-start,window_area_mA=metrics,
            combined_vs_product_max_absolute_difference=nonfactor));save()
        fig,axes=plt.subplots(2,1,sharex=True,figsize=(9,6),height_ratios=[3,1])
        for mode in names[1:]:axes[0].plot(wave,normal[mode],label=mode)
        axes[0].plot(wo,fo,color='black',lw=.8,label='Kitt Peak corrected')
        focus=(wave>=6300.)&(wave<=6300.65)
        low=min(float(observed[focus].min()),*(float(curve[focus].min()) for curve in normal.values()))
        axes[0].set(xlim=(6300.,6300.65),ylabel='Normalized flux',title=f'{name}: trace LTE Ni + fixed NLTE O; natural Voigt only')
        axes[0].set_ylim(low-.015,1.005)
        axes[1].ticklabel_format(axis='x',useOffset=False,style='plain')
        axes[0].legend(fontsize=7);axes[1].plot(wave,normal['oxygen_plus_nickel']-product)
        axes[1].set(xlabel='Air wavelength (Angstrom)',ylabel='Mixed - product')
        fig.tight_layout();fig.savefig(args.output/f'{name}_blend.png',dpi=150);plt.close(fig)
        print(name,metrics,nonfactor,flush=True)
    report['status']='DIAGNOSTIC_COMPLETE';save()


if __name__=='__main__':main()

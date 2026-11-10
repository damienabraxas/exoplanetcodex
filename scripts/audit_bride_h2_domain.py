"""Qualify held H2 chemistry domain; no out-of-domain zero or RT substitution."""
import argparse,json,sys,inspect,platform
from importlib.metadata import version
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def domain_mask(temperature,limits):
    T=np.asarray(temperature)
    if not np.isfinite(T).all()or np.any(T<=0):raise ValueError('Invalid temperature')
    return (T>=limits[0])&(T<=limits[1])


def require_full_domain(temperature,limits):
    if not domain_mask(temperature,limits).all():raise ValueError('H2 fit does not cover complete column; no zero/extrapolation fallback')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if platform.system()!='Darwin' or version('lightweaver')!='0.16.1':raise RuntimeError('Pinned Mac/LW required')
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_opacity_binding import binding_digest
    import lightweaver as lw
    from lightweaver import rh_atoms
    from lightweaver.atomic_table import PeriodicTable
    from lightweaver.molecule import MolecularTable,Molecule
    from lightweaver.utils import get_data_path
    path=ROOT/'data/results/rya1281/solar_bound_v1/runs/rya1281-solar-bound-v1/solar5772_bound_Fe_sampled/experimental_input.json'
    b=json.loads(path.read_text());binding=b['local_binding']
    if binding_digest(binding)!=b['binding_sha256']:raise ValueError('Changed input binding')
    T=np.array(binding['fields']['temperature_K']);height=np.array(binding['fields']['height_m'])
    ne=np.array(binding['ne_m3']);nh=np.array(binding['nH_m3'])
    source=Path(get_data_path())/'DefaultMolecules/H2.molecule';molecule=Molecule(str(source))
    covered=domain_mask(T,molecule.formationTempRange);index=np.flatnonzero(covered)
    if not len(index):raise ValueError('No covered depths')
    a=lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,depthScale=height[covered].copy(),temperature=T[covered].copy(),
        ne=ne[covered].copy(),nHTot=nh[covered].copy(),vturb=np.full(len(index),binding['scalars']['vturb_m_s']),vlos=np.zeros(len(index)))
    abundance=lw.AtomicAbundance({PeriodicTable[i+1]:v for i,v in enumerate(binding['native_abundances_dex'])})
    radset=lw.RadiativeSet([rh_atoms.H_6_atom(),rh_atoms.He_9_atom()],abundance=abundance)
    baseline=radset.compute_eq_pops(a);with_h2=radset.compute_eq_pops(a,mols=MolecularTable([str(source)]))
    np.testing.assert_array_equal(a.ne,ne[covered]);np.testing.assert_array_equal(a.nHTot,nh[covered])
    before=baseline.atomicPops[PeriodicTable['H']]
    h0before=before.n[[l.stage==0 for l in before.model.levels]].sum(axis=0)
    h2=with_h2['H2'];state=with_h2.atomicPops[PeriodicTable['H']]
    h0=state.n[[l.stage==0 for l in state.model.levels]].sum(axis=0)
    atoms=state.n.sum(axis=0);hminus=with_h2.HminPops
    residual=(atoms+2*h2+hminus)/nh[covered]-1
    if not np.isfinite(h2).all()or np.any(h2<0):raise ValueError('Invalid covered chemistry')
    values=np.full((len(T),8),np.nan);values[:,0]=np.arange(len(T));values[:,1]=T;values[:,2]=covered
    values[covered,3]=h2;values[covered,4]=h0;values[covered,5]=h2/nh[covered];values[covered,6]=residual;values[covered,7]=h0before
    args.output.mkdir(parents=True,exist_ok=False);out=args.output/'covered_chemistry.csv'
    np.savetxt(out,values,delimiter=',',header='depth_index,temperature_K,covered,nH2_m3,nHI_with_H2_m3,nH2_over_total_H,nuclei_residual_including_Hminus,nHI_without_H2_m3',comments='')
    result=dict(schema='codex.bride_h2_domain_audit/1',scientific_admission='HOLD',generator=fingerprint(__file__),binding=fingerprint(path),
        molecule=fingerprint(source),molecule_code=fingerprint(inspect.getsourcefile(Molecule)),chemistry_code=fingerprint(inspect.getsourcefile(lw.RadiativeSet)),
        H_He_atoms=fingerprint(rh_atoms.__file__),domain_K=molecule.formationTempRange,covered_count=int(covered.sum()),
        uncovered_depth_indices=np.flatnonzero(~covered).tolist(),uncovered_temperatures_K=T[~covered].tolist(),full_column_admission=False,
        max_covered_nH2_over_total_H=float(np.max(h2/nh[covered])),max_covered_nuclei_residual=float(np.max(abs(residual))),chemistry_native_relative_step_limit=1e-3,
        max_covered_neutral_H_depletion=float(np.max(1-h0/h0before)),lightweaver_version=version('lightweaver'),
        covered_chemistry=fingerprint(out),interpretation='In-domain H2-only LTE chemical audit at fixed ne/nH, no RT or O update. NaN in CSV denotes explicitly uncovered chemistry, never zero. Native fit returns zero outside range; full-column consumer must refuse. No high-T extension, partial H2 patch or new continuum installed.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print('Covered',result['covered_count'],'uncovered',result['uncovered_depth_indices'],'H2max',result['max_covered_nH2_over_total_H'],'conservation',result['max_covered_nuclei_residual'])


if __name__=='__main__':main()

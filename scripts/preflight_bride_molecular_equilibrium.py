"""Audit molecular fit ranges; CO+H2 LTE conservation on explicitly supported states."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from dataclasses import asdict
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.bride_binding_controls import BindingControl
    from pipeline.bride_H_nlte_control import make_context
    from pipeline.multi_atom import read_multi_atom
    from lightweaver.molecule import Molecule,MolecularTable
    from lightweaver.atomic_table import PeriodicTable
    from lightweaver.utils import get_default_molecule_path
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--state',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    state_receipt=fingerprint(a.state);state=BindingControl(**json.loads(a.state.read_text())['state']);temperature=np.array(state.temperature_K)
    directory=require_codex_path('engines.lightweaver_default_molecules')
    if directory.resolve()!=Path(get_default_molecule_path()).resolve():raise ValueError('Installed/registered molecule data mismatch')
    inventory=[];inputs={}
    for path in sorted(directory.glob('*.molecule')):
        if path.stem=='CO_NLTE':continue
        mol=Molecule(str(path));lo,hi=mol.formationTempRange;inputs[str(path)]=fingerprint(path)
        inventory.append(dict(name=mol.name,source=fingerprint(path),range_K=[lo,hi],unsupported_depths=int(np.sum((temperature<lo)|(temperature>hi)))))
    if len(inventory)!=18:raise ValueError('Unexpected installed molecule inventory')
    co=next(x for x in inventory if x['name']=='CO')
    if co['unsupported_depths']:raise ValueError('CO-only diagnostic outside source fit range')
    h2=next(x for x in inventory if x['name']=='H2');lo,hi=h2['range_K'];supported=(temperature>=lo)&(temperature<=hi)
    excluded=np.flatnonzero(~supported).tolist()
    # LW's molecular solver explicitly requires H among the molecular nuclei.
    # A CO+H2 diagnostic is defined only on the intersection of both fit ranges;
    # these independent EOS states are not a truncated RT atmosphere.
    values=asdict(state)
    for key,value in values.items():
        if isinstance(value,tuple):values[key]=tuple(np.asarray(value)[supported])
    state=BindingControl(**values);temperature=temperature[supported]
    ref=json.loads(require_codex_path('experiments.bride_oxygen_reference').read_text());structure=read_marcs_structure(require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name'])
    ctx,_=make_context(state,structure,read_multi_atom(require_codex_path('experiments.bride_oxygen_atom')),rays=5,qscale=1,extra_vac_nm=np.array([500.,845.]),cap_length_m=0.,ng=False)
    spectrum=ctx.kwargs['spect'];empty=ctx.kwargs['eqPops'];mols=MolecularTable([str(directory/'CO.molecule'),str(directory/'H2.molecule')]);eq=spectrum.radSet.compute_eq_pops(ctx.atmos.pyAtmos,mols=mols)
    O=PeriodicTable['O'];C=PeriodicTable['C'];CO=np.asarray(eq['CO']).copy();before=empty.atomicPops[O].nStar.sum(axis=0);atomic=eq.atomicPops[O].nStar.sum(axis=0)
    if any(not np.isfinite(x).all()or np.any(x<0)for x in [CO,*eq.molecularPops,eq.atomicPops[O].nStar,eq.atomicPops[C].nStar]):raise ValueError('Nonfinite/negative chemical equilibrium output')
    errors={}
    for element in [O,C]:
        total=empty.atomicPops[element].nTotal;after=eq.atomicPops[element].nStar.sum(axis=0)+CO
        errors[str(element)]=float(np.max(abs(after/total-1)))
    if any(not np.isfinite(v)or v>1e-8 for v in errors.values()):raise ValueError('CO-only element conservation failed')
    a.output.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(a.output/'arrays.npz',height_m=np.array(state.height_m),temperature_K=temperature,empty_O_nstar_m3=empty.atomicPops[O].nStar,CO_only_O_nstar_m3=eq.atomicPops[O].nStar,CO_m3=CO,atomic_O_total_m3=atomic,empty_O_total_m3=before)
    for path,receipt in inputs.items():
        if fingerprint(path)!=receipt:raise ValueError('Changed molecular data')
    if fingerprint(a.state)!=state_receipt:raise ValueError('Changed thermodynamic input')
    report=dict(generator=fingerprint(__file__),scientific_admission='HOLD',source_state=state_receipt,inventory=inventory,inputs=inputs,excluded_original_depth_indices=excluded,included_states=int(supported.sum()),full_network_range_admission='HOLD'if any(x['unsupported_depths']for x in inventory)else 'ranges_only_supported',explicit_control='CO+H2 LTE on supported independent EOS states; not fallback full chemistry or truncated RT',element_conservation_max_relative=errors,CO_O_reservoir_fraction_range=[float(np.min(CO/before)),float(np.max(CO/before))],limitations=['No RT/SE update, molecular line opacity or molecular NLTE.','Incomplete CO+H2 network is not a scientifically complete EOS.','Unsupported states explicitly excluded only from the independent EOS diagnostic; no coefficient extrapolation or full-column claim.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['full_network_range_admission','element_conservation_max_relative','CO_O_reservoir_fraction_range']}))


if __name__=='__main__':main()

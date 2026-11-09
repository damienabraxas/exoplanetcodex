"""Prepare a held exact-node comparison manifest on Sirius; never run a solver."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    from config.constants import assert_on_sirius,require_codex_path
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_reference import fingerprint,read_marcs_structure,read_exact_departure_record
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model-stem',required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    assert_on_sirius('Bride independent reference intake')
    if Path(args.model_stem).name!=args.model_stem:raise ValueError('Filename stem only')
    held=require_codex_path('grids.gerber_ts')
    receipt_path=held/'O_gerber2023.prov.json'
    receipt=json.loads(receipt_path.read_text())
    atom_path=held/receipt['files']['model_atom']['name']
    if hashlib.md5(atom_path.read_bytes()).hexdigest()!=receipt['files']['model_atom']['md5']:
        raise ValueError('Held oxygen atom does not match acquisition pin')
    atom=read_multi_atom(atom_path)
    atmosphere_path=require_codex_path('grids.marcs_standard')/(args.model_stem+'.mod')
    atmosphere=read_marcs_structure(atmosphere_path)
    if atmosphere['model_id']!=args.model_stem:raise ValueError('MARCS input identity mismatch')
    aux_path=held/receipt['files']['grid_1d_aux']['name']
    if hashlib.md5(aux_path.read_bytes()).hexdigest()!=receipt['files']['grid_1d_aux']['md5']:
        raise ValueError('O aux does not match acquisition pin')
    binary=held/receipt['files']['grid_1d_bin_zip']['name'].removesuffix('.zip')
    reference=read_exact_departure_record(binary,aux_path.read_text(),
        model_id=args.model_stem,abundance=atom.abundance,nlevels=len(atom.levels),ndep=atmosphere['ndep'])
    def plain(v):return v.tolist()if isinstance(v,np.ndarray)else v
    doc=dict(schema='codex.bride_oxygen_reference/1',science_executed=False,
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        generator=fingerprint(__file__),atom=fingerprint(atom_path),receipt=fingerprint(receipt_path),
        aux=fingerprint(aux_path),atmosphere_file=fingerprint(atmosphere_path),
        atmosphere={k:plain(v)for k,v in atmosphere.items()},
        reference={k:plain(v)for k,v in reference.items()},
        binary=dict(name=binary.name,bytes=binary.stat().st_size,
            header_and_record_identity_checked=True,full_uncompressed_hash_not_recomputed=True),
        tau_axis_diagnostics=dict(max_delta_logtau_ross=float(np.max(np.abs(np.log10(reference['tau'])-atmosphere['logtau_ross']))),max_delta_logtau_500=float(np.max(np.abs(np.log10(reference['tau'])-atmosphere['logtau_500'])))),
        status='Read-only independent MULTI record; exact physics/opacity/EOS parity remains to validate; no solver comparison declared passed')
    with args.output.open('x')as f:json.dump(doc,f,indent=2,allow_nan=False);f.write('\n')

if __name__=='__main__':main()

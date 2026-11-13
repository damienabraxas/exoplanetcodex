"""Replay FAC build/example and O level identities; no rates or accuracy claim."""
import argparse
import hashlib
import json
from pathlib import Path

from rya1354_native_parts import load_parts
from rya1354_fac_bootstrap import PIN,level_inventory
from rya1354_fac_oxygen_levels import input_text,inspect_inventory

ROOT=Path(__file__).resolve().parents[1]
BUILD='evidence/fac_bootstrap_v17/'
OXYGEN='runs/fac_oxygen_levels_v17/'


def summarize(base):
    build,bm=load_parts(base/'fac_complete_v17')
    oxygen,om=load_parts(base/'fac_oxygen_complete_v17')
    receipt=json.loads(build[BUILD+'receipt.json'])
    o=json.loads(oxygen[OXYGEN+'receipt.json'])
    if (receipt['status']!='FAC_BUILD_AND_UPSTREAM_INVENTORY_COMPLETED_NOT_O_VALIDATED'
            or receipt['pin']!=PIN or receipt['tracked_source_diff']
            or o['status']!='FAC_O_LEVEL_INVENTORY_COMPLETED_REVIEW_REQUIRED' or o['pin']!=PIN
            or o['binary_sha256']!=receipt['binary']['sha256']
            or o['bootstrap_receipt_sha256']!=hashlib.sha256(build[BUILD+'receipt.json']).hexdigest()):
        raise ValueError('FAC build/oxygen provenance or scope differs')
    if oxygen[OXYGEN+'input.sf'].decode()!=input_text():
        raise ValueError('FAC oxygen physics/configuration input differs')
    if any(name.endswith(('.tr','.tr.b','.ct')) for name in oxygen):
        raise ValueError('Unexpected transition in level-only scope')
    reference=level_inventory(build[BUILD+'source_documents/demo/structure/ref/ne.lev'].decode())
    actual=level_inventory(build[BUILD+'upstream_structure/ne.lev'].decode())
    if actual!=reference or len(actual)!=receipt['upstream_level_inventory']['levels']:
        raise ValueError('FAC upstream example literal level replay differs')
    levels=inspect_inventory(oxygen[OXYGEN+'o.lev'].decode())
    if levels!=o['levels']:
        raise ValueError('FAC oxygen receipt level inventory differs')
    return dict(schema='rya1354.fac-feasibility-audit.v17',pin=PIN,
        binary_sha256=receipt['binary']['sha256'],upstream_example_levels=len(actual),
        upstream_printed_level_rows_exact=True,oxygen_levels=levels,
        oxygen_physics=o['physics'],build_parts=bm['parts'],oxygen_parts=om['parts'],
        scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        numerical_convergence_admitted=False,physical_LS_terms_admitted=False,
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interpretation='Independent implementation bootstrapped and minimal O configuration/parity/J inventory reproduced. No FAC O rates/lifetime, GRASP state matching, converged comparison, calibration or grade. Native potential stopping diagnostics and physical state mapping remain owed.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.write_text(json.dumps(summarize(args.base),indent=2)+'\n')

"""Freeze native Fe I LTE opacity and neutral populations for a bounded O control."""
import argparse
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
    from pipeline.bride_reference import fingerprint,read_marcs_structure,marcs_number_densities
    from pipeline.bride_ts_budget import read_ts_opacity,read_ts_components
    from pipeline.wavelength_util import air_to_vac
    report=json.loads(args.run.read_text())
    if report.get('support')!='nuv200_350' or not report['reporting_preserves_serialized_coefficients']:
        raise ValueError('Requires verified broad NUV baseline/reporting pair')
    case=next(c for c in report['cases']if c['name']=='NUV_source_air_reporting')
    directory=args.run.parent/case['name']
    for name,key in [('continuum.opa','output'),('detail.txt','detail')]:
        if fingerprint(directory/name)!=case[key]:raise ValueError('Changed native output')
    data=read_ts_opacity(directory/'continuum.opa')
    comp=read_ts_components(directory/'detail.txt',data['ndep'],data['nwave'])
    path=require_codex_path('grids.marcs_standard')/report['atmosphere']['name']
    if fingerprint(path)!=report['atmosphere']:raise ValueError('Changed source atmosphere')
    structure=read_marcs_structure(path);ne,nH=marcs_number_densities(structure)
    fe=comp['labels'].index('FeI');native=comp['records'][:,1:,2].T*.1*comp['fractions'][:,1:,fe].T
    state=comp['states'];T,Pe,rho,U,pFe,pFeII=state.T
    perkg=pFe/(1.38066e-16*T)*1e6/(rho*1000)
    wave=air_to_vac(data['wavelength_air_A'])*.1
    if np.any(native<0) or not np.isfinite(native).all() or np.any(perkg<=0):raise ValueError('Invalid Fe packet')
    np.testing.assert_allclose(comp['post'][:,:,3][:,1:].T*.1,data['absorption_m2_kg'],rtol=6e-5,atol=0)
    args.output.mkdir(parents=True,exist_ok=False);array=args.output/'packet.npz'
    np.savez_compressed(array,wavelength_vac_nm=wave,opacity_m2_kg=native,
        opacity_net_per_neutral_Fe_m2=native/perkg[None,:],native_neutral_Fe_per_kg=perkg,
        expected_temperature_K=structure['temperature_K'],expected_ne_m3=ne,expected_nH_m3=nH,
        expected_rho_kg_m3=structure['rho_kg_m3'],native_temperature_K=T,native_rho_kg_m3=rho*1000)
    result=dict(schema='codex.bride_fe_packet/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        run=fingerprint(args.run),source_detail=case['detail'],source_opacity=case['output'],
        atmosphere=report['atmosphere'],array=fingerprint(array),
        support_vac_nm=[float(wave[0]),float(wave[-1])],ndep=data['ndep'],nwave=data['nwave'],
        interpretation='Native TS FeI bf component only. Raw specific opacity and net opacity per native neutral Fe retained. Exact native atmospheric node only; no depth/T/ne interpolation. Finite support intervention, not full producer parity or source-independent cross-sections.')
    (args.output/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()

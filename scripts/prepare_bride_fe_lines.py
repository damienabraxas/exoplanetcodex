"""Intake exact GES source matches for canonical single-component O777-region Fe I."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def source_packet():
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.wavelength_util import air_to_vac
    from pipeline.bride_trace_fe_lines import validate_lines
    canonical=ROOT/'data/linelists/canonical_gf.csv'
    source=require_codex_path('engines.ispec')/'input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv'
    with canonical.open() as f:rows=[r for r in csv.DictReader(f) if r['species']=='Fe I' and 7763<=float(r['wavelength_air_A'])<=7783]
    with source.open() as f:original=[r for r in csv.DictReader(f,delimiter='\t') if r['element']=='Fe 1' and 7763<=float(r['wave_A'])<=7783]
    selected=[];excluded=[];evidence=[]
    for c in rows:
        if int(c['hfs_n_components'])!=1:
            excluded.append(dict(canonical_row=c,reason='Combined component row; explicit source component treatment owed'));continue
        matches=[r for r in original if abs(float(r['wave_A'])-float(c['wavelength_air_A']))<1e-6
            and abs(float(r['lower_state_eV'])-float(c['excitation_potential_eV']))<1e-6
            and abs(float(r['loggf'])-float(c['log_gf']))<1e-6]
        if len(matches)!=1:raise ValueError('Missing/ambiguous exact GES match: '+c['line_id'])
        r=matches[0];gamma=float(r['turbospectrum_rad'])
        if gamma<=0 or abs(np.log10(gamma)-float(r['rad']))>.005:raise ValueError('Unverified natural damping units')
        selected.append(dict(line_id=c['line_id'],lambda0_vac_nm=float(air_to_vac(float(c['wavelength_air_A']))/10),
            excitation_eV=float(c['excitation_potential_eV']),loggf=float(c['log_gf']),gamma_rad_s_inv=gamma))
        evidence.append(dict(canonical_row=c,source_row=r))
    validate_lines(selected)
    result=dict(schema='codex.bride_trace_fe_lines/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        canonical=fingerprint(canonical),source=fingerprint(source),lines=selected,evidence=evidence,excluded=excluded,
        support_vac_nm=list(map(float,air_to_vac(np.array([7763.,7783.]))/10)),
        interpretation='Exact source/canonical line data only. Natural Voigt; pressure/HFS/isotope/molecule and other species absent. No fitted Fe abundance, damping or opacity scale.')
    return result


def read_packet(path):
    actual=json.loads(Path(path).read_text())
    if actual!=source_packet():raise ValueError("Trace packet differs from fresh source/canonical intake")
    return actual


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",type=Path,required=True);args=p.parse_args()
    result=source_packet()
    args.output.mkdir(parents=True,exist_ok=False);(args.output/'lines.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Source-matched',len(result['lines']),'excluded combined-component',len(result['excluded']))


if __name__=='__main__':main()

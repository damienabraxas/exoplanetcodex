"""Explicit vacuum-to-air census of the separately held CO transition source."""
import argparse,csv,json,shlex,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.molecular_lists import VENDORED_DIR,MANIFEST_PATH,load_manifest
    from pipeline.molecular_canonical import canonical_rows
    from scripts.migrate_molecular_transitions_rya1130 import metadata_for
    from pipeline.wavelength_util import vac_to_air
    plan_path=ROOT/'config/experiments/rya1281_co_source.json';plan=json.loads(plan_path.read_text())
    entry=load_manifest()['molecules']['CO']
    if plan['schema']!='codex.bride_co_source_plan/1' or not entry['conversion']['recipe'].startswith('vacuum'):raise ValueError('CO frame not established by declared conversion')
    if len(entry['files'])!=1:raise ValueError('Changed CO holdings')
    filename=entry['files'][0];path=VENDORED_DIR/entry['vendored_subdir']/filename
    installed=require_codex_path('engines.ispec')/'input/linelists/turbospectrum/molecules'/filename
    if path.read_bytes()!=installed.read_bytes():raise ValueError('Installed CO differs from vendor')
    with path.open()as stream:expected=int(shlex.split(stream.readline())[-1])
    lo,hi=plan['range_air_A'];rows=[];parsed=0;raw_count=0
    for row in canonical_rows(path,'CO',metadata_for(entry)):
        parsed+=1;vac=row['wavelength']['value_A']
        if not np.isfinite(vac)or vac<=0:raise ValueError('Invalid CO vacuum wavelength')
        if lo<=vac<=hi:raw_count+=1
        air=float(vac_to_air(vac))
        if lo<=air<=hi:
            rows.append(dict(transition_id=row['transition_id'],isotopologue=row['isotopologue'],stored_vacuum_A=vac,
                converted_air_A=air,lower_energy_eV=row['lower_energy_eV'],stored_loggf=row['transition_probability']['log_gf'],
                raw_label=row['transition_identity']['raw_label'],intake_status=row['intake_status']))
    if parsed!=expected or parsed!=entry['line_count']:raise ValueError('CO declared row count mismatch')
    result=dict(schema='codex.bride_co_source/1',scientific_admission='HOLD',generator=fingerprint(__file__),plan=fingerprint(plan_path),
        manifest=fingerprint(MANIFEST_PATH),source=fingerprint(path),installed=fingerprint(installed),conversion=entry['conversion'],
        wavelength_kernel=fingerprint(ROOT/'pipeline/wavelength_util.py'),canonical_parser=fingerprint(ROOT/'pipeline/molecular_canonical.py'),
        parsed_rows=parsed,raw_vacuum_window_count=raw_count,converted_air_window_count=len(rows),
        automatic_bsyn_glob_selects_file=path.suffix=='.bsyn',scope=plan['scope'],
        qualification='Canonical source status retained; consumer automatic routing and separate dedicated CO production paths must not be conflated')
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with (args.output/'rows.csv').open('x')as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    print(json.dumps(dict(parsed=parsed,raw_window=raw_count,air_window=len(rows),automatic_bsyn_selected=result['automatic_bsyn_glob_selects_file'])))


if __name__=='__main__':main()

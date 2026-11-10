"""Receipt the complete derived R/C/J archive held outside git's size ceiling."""
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    key='experiments.bride_thermal_refined_rates';path=require_codex_path(key)
    r=json.loads(args.report.read_text())
    if r['schema']!='codex.bride_thermal_flux/1' or not r['numeric_refinement']or fingerprint(path)!=r['coupled']['rates']:raise ValueError('Wrong/changed external refinement rates')
    result=dict(schema='codex.bride_external_derived_artifact/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        resource_key=key,artifact=fingerprint(path),source_report=fingerprint(args.report),source_writer=r['coupled']['helper'],
        scope='Full derived radiative/collision/J arrays, not raw observation or adopted product. Retained in registered Mac storage under repository5MiB ceiling; no renaming or truncating numerical content.')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x')as f:f.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Verified external bytes',result['artifact']['bytes'],'SHA256',result['artifact']['sha256'])

if __name__=='__main__':main()

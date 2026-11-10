"""Receipt immutable native outputs without rewriting source padding."""
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--first',type=Path,required=True);p.add_argument('--final',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_feii_profile import native_recipe_record
    accepted=json.loads((args.final/'report.json').read_text())
    if (args.first/'report.json').exists():raise ValueError('Original v1 was not a rejected source audit')
    first=(args.first/'stdout.txt').read_text();final=(args.final/'stdout.txt').read_text()
    try:native_recipe_record(first)
    except ValueError:pass
    else:raise ValueError('Original malformed recipe unexpectedly accepted')
    recipe=native_recipe_record(final)
    widths=lambda text:[s for s in text.splitlines() if s.split() and s.split()[0]=='FEII_WIDTH']
    if len(widths(first))!=54 or widths(first)!=widths(final):raise ValueError('Native numerical records changed')
    if (args.first/'input.txt').read_bytes()!=(args.final/'input.txt').read_bytes():raise ValueError('Native inputs changed')
    for name,field in [('driver.f','driver'),('stdout.txt','stdout'),('input.txt','input.txt'),('feii_oracle','binary')]:
        key='input' if field=='input.txt' else field
        if fingerprint(args.final/name)!=accepted[key]:raise ValueError('Changed accepted native evidence')
    result=dict(schema='codex.bride_feii_native_receipts/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        first_resource='experiments.bride_feii_width_source',final_resource='experiments.bride_feii_width_final',
        first_status='REJECTED_RECIPE_PARSE',final_status='ACCEPTED_ARITHMETIC_ONLY',accepted_report=fingerprint(args.final/'report.json'),
        first_code='9689ef23',delimiter_fix='842ede414e49b25051710a0ec7e894deda8f26b2',
        first={name:fingerprint(args.first/name) for name in ('driver.f','stdout.txt','input.txt','feii_oracle')},
        final={name:fingerprint(args.final/name) for name in ('driver.f','stdout.txt','input.txt','feii_oracle')},
        final_recipe=list(recipe),native_width_records_equal=54,input_bytes_equal=True,
        scope='Raw source-derived drivers/stdout retain padding externally, never stripped or overwritten. Receipts are not rewritten native values or physical profile admission.')
    args.output.mkdir(parents=True,exist_ok=False);(args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()

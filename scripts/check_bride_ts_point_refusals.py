"""Exercise native point-adapter refusal, retaining failed diagnostic outputs."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from pipeline.bride_reference import fingerprint
    run = json.loads(args.run.read_text())
    build = json.loads((args.build/'manifest.json').read_text())
    if fingerprint(args.build/'manifest.json') != run['build'] or run['status'] != 'DIAGNOSTIC_COMPLETE':
        raise ValueError('Changed build or incomplete point run')
    case = next(c for c in run['cases'] if c['name'] == 'native_identity')
    entry = next(r for r in case['runs'] if r['kind'] == 'baseline')
    directory = args.run.parent/'native_identity_baseline'
    for name, key in (('point_state.txt', 'points'), ('input.txt', 'input')):
        if fingerprint(directory/name) != entry[key]:
            raise ValueError('Changed native probe fixture')
    binary = args.build/'baseline'/'babsma_lu'
    if fingerprint(binary) != next(e for e in build['builds'] if e['name'] == 'baseline')['program']:
        raise ValueError('Changed point binary')
    lines = (directory/'point_state.txt').read_text().splitlines(keepends=True)
    probes = [('wrong_count', str(int(lines[0])-1)+'\n'+''.join(lines[1:]), 72),
              ('truncated', ''.join(lines[:-1]), 73), ('extra_row', ''.join(lines)+'1\n', 76)]
    args.output.mkdir(parents=True, exist_ok=False)
    result = dict(schema='codex.bride_ts_point_refusals/1', scientific_admission='HOLD', generator=fingerprint(__file__),
                  source_run=fingerprint(args.run), build=run['build'], probes=[])
    for name, points, expected in probes:
        out = args.output/name
        out.mkdir()
        (out/'DATA').symlink_to((args.run.parent/'DATA').resolve(), target_is_directory=True)
        (out/'point_state.txt').write_text(points)
        text = (directory/'input.txt').read_text()
        (out/'input.txt').write_text(text)
        proc = subprocess.run([str(binary.resolve())], input=text, cwd=out, capture_output=True, text=True, timeout=60)
        (out/'stdout.log').write_text(proc.stdout)
        (out/'stderr.log').write_text(proc.stderr)
        if proc.returncode != expected:
            raise ValueError(f'{name}: expected native refusal{expected}, got{proc.returncode}')
        result['probes'].append(dict(name=name, returncode=proc.returncode, expected=expected,
            points=fingerprint(out/'point_state.txt'), stderr=fingerprint(out/'stderr.log')))
    (args.output/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    print('All three malformed point-table probes refused with exact native stop codes')


if __name__ == '__main__':
    main()

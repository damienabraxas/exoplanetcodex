"""Build isolated input-only TS derivatives from the verified native/reporting builds."""
import argparse
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if platform.system() != 'Darwin':
        raise RuntimeError('Mac-only isolated diagnostic')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_ts_points import point_driver_source
    base = require_codex_path('experiments.bride_ts_budget_build')
    original = json.loads((base/'manifest.json').read_text())
    if original['source_commit'] != '21269948750e2177fe32c5df38e1d9bb81363374':
        raise ValueError('Changed native source pin')
    args.output.mkdir(parents=True, exist_ok=False)
    result = dict(schema='codex.bride_ts_point_build/1', scientific_admission='HOLD',
        base_manifest=fingerprint(base/'manifest.json'), generator=fingerprint(__file__),
        adapter=fingerprint(ROOT/'pipeline/bride_ts_points.py'), builds=[])
    for entry in original['builds']:
        name = entry['name']
        source = base/name
        for filename, key in (('babsma_lu', 'program'), ('babsma.f', 'babsma_source'), ('detabs.f', 'detabs_source'), ('absko.f', 'absko_source')):
            if fingerprint(source/filename) != entry[key]:
                raise ValueError('Changed verified base build: '+filename)
        originals = {r['name']: r for r in original['source_files']}
        for path in source.iterdir():
            if path.suffix in ('.f', '.f90', '.inc', '.par') and path.name not in ('babsma.f', 'detabs.f', 'absko.f'):
                if path.name not in originals or fingerprint(path) != originals[path.name]:
                    raise ValueError('Changed copied native source/include: '+path.name)
        directory = args.output/name
        shutil.copytree(source, directory)
        path = directory/'babsma.f'
        path.write_text(point_driver_source(path.read_text()))
        # Reuse exact existing native compiler recipe; no new flags or libraries.
        run = subprocess.run(entry['command'], cwd=directory, capture_output=True, text=True, timeout=90)
        (directory/'compiler.log').write_text(run.stdout+run.stderr)
        if run.returncode:
            raise RuntimeError('Point adapter build failed: '+name)
        result['builds'].append(dict(name=name, program=fingerprint(directory/'babsma_lu'),
            driver=fingerprint(path), compiler_log=fingerprint(directory/'compiler.log'), command=entry['command']))
    (args.output/'manifest.json').write_text(json.dumps(result, indent=2)+'\n')
    print('Built isolated baseline/reporting point adapters')


if __name__ == '__main__':
    main()

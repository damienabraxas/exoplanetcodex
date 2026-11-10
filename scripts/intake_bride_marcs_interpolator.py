"""Receipt the official MARCS interpolation archive without executing its contents."""
import argparse
import json
from pathlib import Path, PurePosixPath
import shutil
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
URL = 'https://marcs.astro.uu.se/documents/auxiliary/interpol_marcs.tar.gz'
SHA256 = '8a3d0fbd15a378830ece79c8cdc3c8ee60d74f2fd385f03b1bb3fbecb762f189'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from pipeline.bride_reference import fingerprint
    if fingerprint(args.archive)['sha256'] != SHA256:
        raise ValueError('Changed official archive; new intake review required')
    with tarfile.open(args.archive) as archive:
        members = archive.getmembers()
        for member in members:
            p = PurePosixPath(member.name)
            if p.is_absolute() or '..' in p.parts or not (member.isfile() or member.isdir()):
                raise ValueError('Unsafe archive member')
        args.output.mkdir(parents=True, exist_ok=False)
        shutil.copyfile(args.archive, args.output/args.archive.name)
        files = []
        for member in members:
            path = args.output/member.name
            if member.isdir():
                path.mkdir(parents=True, exist_ok=True)
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                with archive.extractfile(member) as source, path.open('xb') as dest:
                    shutil.copyfileobj(source, dest)
                files.append(dict(relative_path=member.name, **fingerprint(path)))
    result = dict(schema='codex.bride_marcs_interpolator_intake/1', scientific_admission='HOLD',
        generator=fingerprint(__file__), source_url=URL, archive=fingerprint(args.output/args.archive.name), files=files,
        execution=False, redistribution_license='UNVERIFIED; source kept externally; no vendored upstream code or manual',
        purpose='Official source/manual inspection and prospective independent interpolation comparison. Not the native atmosphere-generation source, not physical admission. Upstream shell/control scripts are data only and were not executed.')
    (args.output/'manifest.json').write_text(json.dumps(result, indent=2)+'\n')
    print('Official archive safely receipted; no upstream execution')


if __name__ == '__main__':
    main()

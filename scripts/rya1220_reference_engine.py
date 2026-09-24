"""Make an isolated iSpec view with Brooke CN coverage for N I 10108 A.

Original engine files and canonical gf are untouched. Reuses the RYA-1214
primary-table parser and its independent energy/gf checks. The extra molecular
file is deliberately limited to a previously uncovered interval.
"""
import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--engine-root', type=Path, required=True)
    ap.add_argument('--ispec-root', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    sys.path[:0] = [str(a.engine_root), str(a.engine_root/'scripts')]
    import rya1214_build_cn_ir_bsyn as builder
    frame = builder.parse_brooke()
    energy = builder.assert_e_includes_vibration(frame)
    gf_check = builder.validate_against_amarsi(frame)
    selected = frame[frame.lam_air.between(10080, 10130)].sort_values('lam_air')
    if selected.empty:
        raise ValueError('No primary CN rows at the reference line')

    def view(src, dst, expand):
        dst.mkdir(exist_ok=True)
        for item in src.iterdir():
            target = dst/item.name
            if item.name in expand:
                view(item, target, expand[item.name])
            else:
                target.symlink_to(item.resolve(), target_is_directory=item.is_dir())

    # Physical Python source copy is needed because iSpec resolves __file__ to
    # locate its molecular directory. Large model/binary data remain symlinks.
    root = a.out/'ispec'; root.mkdir()
    for item in a.ispec_root.iterdir():
        if item.name == 'ispec':
            shutil.copytree(item, root/item.name, ignore=shutil.ignore_patterns('__pycache__'))
        elif item.name == 'input':
            view(item, root/'input', {'linelists': {'turbospectrum': {'molecules': {}}}})
        else:
            (root/item.name).symlink_to(item.resolve(), target_is_directory=item.is_dir())
    molecules = root/'input/linelists/turbospectrum/molecules'
    import re
    for p in molecules.glob('*.bsyn'):
        match = re.search(r'_(\d+)-(\d+)\.bsyn$', p.name)
        if match and int(match[1]) <= 1013 and int(match[2]) >= 1008:
            raise ValueError(f'Existing overlapping list requires adjudication: {p}')
    dest = molecules/'12C14N_1008-1013.bsyn'
    with dest.open('x') as f:
        f.write(f"'         {builder.SPECIES_CODE}'    1   {len(selected)}\n")
        f.write("'Brooke2014 Table4 CN; RYA-1220 isolated diagnostic coverage'\n")
        for r in selected.itertuples():
            f.write(f"{r.lam_air:12.3f} {r.Elow_eV:9.5f} {r.loggf:7.3f}    0.000 "
                    f"{r.gu:6.1f}  0.00E+00 'X' 'X'   0.0    0.0\n")
    selected.to_csv(a.out/'cn_primary_rows.csv', index=False)
    manifest = {'primary': str(builder.SRC),
                'primary_sha256': hashlib.sha256(builder.SRC.read_bytes()).hexdigest(),
                'parser_sha256': hashlib.sha256(Path(builder.__file__).read_bytes()).hexdigest(),
                'line_count': len(selected), 'energy_check': energy, 'gf_check': gf_check,
                'molecular_file': str(dest),
                'molecular_sha256': hashlib.sha256(dest.read_bytes()).hexdigest(),
                'shared_engine_modified': False, 'status': 'diagnostic_only',
                'limitations': '12C14N only; isotope sensitivity and other molecules not established'}
    (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))


if __name__ == '__main__':
    main()

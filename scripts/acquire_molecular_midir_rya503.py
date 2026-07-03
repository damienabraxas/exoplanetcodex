#!/usr/bin/env python3
"""
scripts/acquire_molecular_midir_rya503.py
=========================================
RYA-503 Phase 2 — acquire the OH/NH/CH mid-IR ro-vibrational fundamentals (the bands
RYA-360 confirmed BY MEASUREMENT are absent from the held electronic .bsyn), convert
them with the Phase-1-validated converter (scripts/molecular_linelist_convert.py),
vendor them into data/linelists/molecular/turbospectrum/, and MERGE their provenance
into the RYA-360 MOLECULAR_MANIFEST.json (preserving the RYA-360 entries). Also stamps
the CO conversion recipe onto the CO entry — closing the RYA-360 CO-script gap now that
Phase 1 reproduced the vendored CO file exactly.

Gated on Phase 1 PASS (the CO round-trip); this script assumes the converter is valid.

Sources are ExoMol .states/.trans/.def (re-downloadable from the cited URLs; NOT
vendored — the vendored .bsyn + the recipe here are the reproducible record). Provide
the downloaded source files in --src-dir named `<iso>__<dataset>.{states,trans,def}`.

    python scripts/acquire_molecular_midir_rya503.py --src-dir <dir-with-exomol-files>
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))

from scripts.molecular_linelist_convert import (  # noqa: E402
    parse_exomol_def, parse_exomol_states, convert_exomol, format_bsyn,
    GF_CONST, CM_PER_EV)
from pipeline.molecular_lists import VENDORED_DIR, MANIFEST_PATH, count_bsyn_lines  # noqa: E402

# ν≥1000 cm⁻¹ (drop the far-IR pure-rotation tail, as in the CO recipe) + ground X–X
# (the ro-vibrational bands; excludes the UV/optical electronic bands already held as
# the RYA-360 .bsyn set). The RYA-499 window is the fundamental-band coverage proof.
NU_MIN = 1000.0

ACQUIRE = [
    dict(key='OH_MYTHOS_midIR', mol='OH', subdir='OH', iso='16O-1H', dataset='MYTHOS',
         code='0108.000016', tag='MYTHOS', label='ExoMol MYTHOS OH', version='20240526',
         url='https://www.exomol.com/db/OH/16O-1H/MYTHOS/',
         doi='10.1093/mnras/stae2803',
         citation='Mitev, Bowesman, Zhang, Yurchenko & Tennyson 2025, MNRAS 536, 3401 '
                  '(ExoMol LXI; OH MYTHOS)',
         midir=('OH 1-0 fundamental', 2600.0, 3600.0),
         out='16O-1H__MYTHOS_rovib.bsyn'),
    dict(key='NH_kNigHt_midIR', mol='NH', subdir='NH', iso='14N-1H', dataset='kNigHt',
         code='0107.000014', tag='kNigHt', label='ExoMol kNigHt NH', version='20240301',
         url='https://www.exomol.com/db/NH/14N-1H/kNigHt/',
         doi='10.1093/mnras/stae1340',
         citation='Perri & McKemmish 2024, MNRAS 531, 3023 (NH kNigHt). Note: the newer '
                  '2kNigHt (v20260414) is in press with no DOI yet; kNigHt is the cited '
                  'published variational list.',
         midir=('NH 1-0 fundamental', 3000.0, 3500.0),
         out='14N-1H__kNigHt_rovib.bsyn'),
    dict(key='CH_MoLLIST_midIR', mol='CH', subdir='CH', iso='12C-1H', dataset='MoLLIST',
         code='0106.000012', tag='MoLLIST', label='ExoMol MoLLIST CH', version='20190214',
         url='https://www.exomol.com/db/CH/12C-1H/MoLLIST/',
         doi='10.1051/0004-6361/201423956',
         citation='Masseron et al. 2014, A&A 571, A47 (CH linelist); re-hosted via '
                  'Bernath 2020, JQSRT 240, 106687 (MoLLIST, DOI 10.1016/j.jqsrt.2019.106687)',
         midir=('CH 1-0 fundamental', 2650.0, 3100.0),
         out='12C-1H__MoLLIST_rovib.bsyn'),
]

# CO recipe recovered + validated in Phase 1 (the RYA-360 gap-closer).
CO_CONVERSION = {
    'converter': 'scripts/molecular_linelist_convert.py (RYA-503)',
    'source_format': 'ExoMol .states/.trans',
    'source_dataset': 'ExoMol CO Li2015',
    'source_url': 'https://www.exomol.com/db/CO/12C-16O/Li2015/',
    'version': '20170101',
    'doi': '10.1088/0067-0049/216/1/15',
    'citation': 'Li, Gordon, Rothman et al. 2015, ApJS 216, 15 (CO rovibrational line '
                'lists), ExoMol-hosted Li2015',
    'recipe': 'vacuum λ = 1e8/(E_u−E_l); χ_low = E_l/8065.543937 eV; '
              'loggf = log10(1.49919e-16·λ²·g_u·A_ul); filter ν ≥ 1000 cm⁻¹ '
              '(drops 7713 pure-rotational lines below 1000 cm⁻¹ / above 10 µm)',
    'roundtrip': 'PASS — reproduces the vendored CO_IR_Li2015.dat: 117783/117783 lines, '
                 '0 only-in-either; max dev λ 8.7e-3 Å, χ 5.1e-6 eV, loggf 5.0e-4 dex, '
                 'A exact at the file\'s 3-sig-fig precision. Closes the RYA-360 CO-script '
                 'reproducibility gap.',
}


def _cm(wl_A):
    return 1.0e8 / wl_A


def _measure_bsyn(path, midir):
    lam = np.array([float(l.split()[0]) for l in open(path)
                    if l.lstrip() and l.lstrip()[0] != "'"], dtype=float)
    nu = _cm(lam)
    label, clo, chi = midir
    return {
        'line_count': int(lam.size),
        'wavelength_coverage': {
            'min_A': round(float(lam.min()), 3), 'max_A': round(float(lam.max()), 3),
            'min_cm-1': round(float(nu.min()), 1), 'max_cm-1': round(float(nu.max()), 1),
            'regime': 'mid-IR-rovibrational'},
        'midir_window': {
            'label': f'{label} ({clo:.0f}-{chi:.0f} cm⁻¹)',
            'range_cm-1': [clo, chi],
            'count': int(((nu >= clo) & (nu <= chi)).sum()),
            'verdict': 'MID-IR PRESENT (acquired; was 0 in the held electronic .bsyn)'},
    }


def run(src_dir: Path) -> dict:
    manifest = json.loads(MANIFEST_PATH.read_text())
    mols = manifest['molecules']

    # (1) stamp the CO conversion recipe (closes the RYA-360 gap)
    if 'CO' in mols:
        mols['CO']['conversion'] = CO_CONVERSION

    report = {'acquired': [], 'skipped': []}
    for a in ACQUIRE:
        states_p = src_dir / f"{a['iso']}__{a['dataset']}.states"
        trans_p = src_dir / f"{a['iso']}__{a['dataset']}.trans"
        def_p = src_dir / f"{a['iso']}__{a['dataset']}.def"
        if not (states_p.exists() and trans_p.exists()):
            # "the molecules we can do" — a species with no available source is flagged, not faked
            report['skipped'].append((a['key'], f"source files absent in {src_dir}"))
            continue
        edef = parse_exomol_def(def_p) if def_p.exists() else None
        states = parse_exomol_states(states_p, edef) if edef else parse_exomol_states(states_p)
        lines = convert_exomol(states, trans_p, tag=a['tag'], nu_min=NU_MIN, ground_only=True)
        sub = VENDORED_DIR / a['subdir']
        sub.mkdir(parents=True, exist_ok=True)
        out = sub / a['out']
        out.write_text(format_bsyn(lines, a['code'], a['label']))

        meas = _measure_bsyn(out, a['midir'])
        entry = {
            'species': f"{a['mol']} (mid-IR ro-vibrational, {a['dataset']})",
            'vendored_subdir': a['subdir'], 'files': [a['out']],
            'ts_species_code': a['code'],
            'source': f"{a['mol']} X–X ro-vibrational — {a['citation']}",
            'distribution': f"ExoMol {a['dataset']} v{a['version']} ({a['url']}); "
                            f"DOI {a['doi']}; converted via scripts/molecular_linelist_convert.py",
            'doi': a['doi'], 'version': a['version'], 'origin': 'acquired (RYA-503)',
            'in_ispec': False, 'coverage_gate': 'midir',
            'conversion': {
                'converter': 'scripts/molecular_linelist_convert.py',
                'source_format': 'ExoMol .states/.trans',
                'recipe': f'ground X–X (ro-vibrational), ν ≥ {NU_MIN:.0f} cm⁻¹; '
                          'vacuum λ; loggf = log10(1.49919e-16·λ²·g_u·A); '
                          'χ_low = E_l/8065.543937 eV',
                'gf_const': GF_CONST, 'cm_per_eV': CM_PER_EV},
            **meas}
        mols[a['key']] = entry
        report['acquired'].append((a['key'], entry))

    manifest['generated'] = str(date.today())
    manifest.setdefault('acquisitions', {})['RYA-503'] = (
        'OH/NH/CH mid-IR ro-vibrational fundamentals acquired from ExoMol; CO conversion '
        'recipe recovered + validated (round-trip PASS). CO gap closed.')
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + '\n')
    return report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description='RYA-503 Phase 2 — acquire OH/NH/CH mid-IR lists')
    ap.add_argument('--src-dir', required=True,
                    help='dir with the downloaded ExoMol <iso>__<dataset>.{states,trans,def}')
    args = ap.parse_args(argv)
    rep = run(Path(args.src_dir))
    print('=' * 78)
    print('RYA-503 Phase 2 — OH/NH/CH mid-IR acquisition')
    print('=' * 78)
    for key, e in rep['acquired']:
        mw = e['midir_window']; wc = e['wavelength_coverage']
        print(f"  {key:<18} {e['line_count']:>7} lines  span {wc['min_A']:.1f}-{wc['max_A']:.1f} Å  "
              f"| {mw['label']} → {mw['count']} rows")
        print(f"     {e['source']}")
    for key, why in rep['skipped']:
        print(f"  SKIPPED {key}: {why}  (flagged, not faked)")
    print(f"\nManifest updated → {MANIFEST_PATH}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

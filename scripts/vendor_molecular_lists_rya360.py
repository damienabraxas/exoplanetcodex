#!/usr/bin/env python3
"""
scripts/vendor_molecular_lists_rya360.py
========================================
RYA-360 Steps 0+1 — verify the C/N/O molecular line lists at the iSpec path,
MEASURE them (HARPS-window headline counts + mid-IR row spans), and VENDOR them into
the git-tracked repo with a provenance manifest.

Reconciliation (the RYA-237 false-absence): the usable Turbospectrum `.bsyn` lists
live in the iSpec TOOL dir `ispec/input/linelists/turbospectrum/molecules/`, NOT the
raw-download dir `data/linelists/molecular/exomol/` (`.bz2`) the RYA-237 recon looked
in. This script reads from the iSpec dir (the RYA-236 audit's authoritative location).

Run once to (re)generate the vendored copies + manifest:
    python scripts/vendor_molecular_lists_rya360.py            # vendor + write manifest
    python scripts/vendor_molecular_lists_rya360.py --measure-only   # Step-0 report, no copy
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import date
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))

from pipeline.molecular_lists import (  # noqa: E402
    ISPEC_MOLECULES_DIR, VENDORED_DIR, MANIFEST_PATH,
    count_bsyn_lines, bsyn_wavelengths, source_label)

# cm⁻¹ ↔ Å (vacuum wavenumber convention; coarse — used only for span reporting +
# the mid-IR absence check, which is decisive by orders of magnitude regardless).
def _cm(wl_A):
    return 1.0e8 / wl_A


# ── the C/N/O molecule registry (what RYA-236 acquired; the ONLY in-scope species) ──
# Per molecule: isotopologue file prefixes, the underlying line-data source (as
# recorded by RYA-236 / cno_synthesis), the distribution/compilation origin, and the
# HARPS headline diagnostic window (the RYA-236 count to reproduce). OH/NH/CH also get
# the RYA-499 mid-IR fundamental windows measured.
_HARPS_RANGE = (3800.0, 6900.0)     # HARPS-VIS arm (RYA-237 keystone)

# RYA-499 mid-IR ro-vibrational fundamental windows (cm⁻¹) → measured per OH/NH/CH file.
_MIDIR_WINDOWS_CM = {
    'OH': ('OH 1-0 fundamental', 2600.0, 3600.0),
    'NH': ('NH 1-0 fundamental', 3000.0, 3500.0),
    'CH': ('CH 1-0 fundamental', 2650.0, 3100.0),
}

_DISTRIB = ("Masseron compilation, VALD-based, distributed with iSpec/Turbospectrum "
            "(Gerber et al. 2023, A&A 669, A43; uploaded 2022 to keeper.mpdl.mpg.de, "
            "split per wavelength region) — see molecules-dir README.md")

MOLECULES = {
    'CH': dict(
        prefixes=['12CH', '13CH'],
        source='CH A-X: Masseron et al. 2014 (MNRAS 442, 1489) / 2022 update',
        distribution=_DISTRIB,
        harps=('CH A-X G-band', 4290.0, 4315.0),   # RYA-236 baseline ≈ 583
        midir=True),
    'CN': dict(
        prefixes=['12C14N', '13C14N', '12C15N'],
        source='CN red A-X: Brooke et al. 2014 / Sneden et al. 2014 (isotopologues)',
        distribution=_DISTRIB,
        harps=('CN red A-X', 6000.0, 6200.0),       # RYA-236 baseline ≈ 3534
        midir=False),
    'C2': dict(
        prefixes=['12C12C', '12C13C', '13C13C'],
        source='C2 Swan d³Πg–a³Πu (Brooke et al. 2013 / Masseron compilation)',
        distribution=_DISTRIB,
        harps=('C2 Swan (0,0)', 5130.0, 5170.0),    # RYA-236 baseline ≈ 1019
        midir=False),
    'OH': dict(
        prefixes=['16OH'],
        source='OH A-X electronic (Masseron compilation, VALD-based)',
        distribution=_DISTRIB,
        harps=None,
        midir=True),
    'NH': dict(
        prefixes=['14NH'],
        source='NH A-X electronic (Masseron compilation, VALD-based)',
        distribution=_DISTRIB,
        harps=None,
        midir=True),
    'CO': dict(
        files=['CO_IR_Li2015.dat'],
        source='CO ExoMol line list, Li et al. 2015 (ApJS 216, 15)',
        distribution='RYA-236 conversion of the ExoMol Li2015 CO list to the '
                     'Turbospectrum babsma .dat format (converted to species code '
                     '0608.012016; conversion script external to this repo — RYA-236)',
        harps=None,
        midir=None),   # CO is the mid-IR exception (see measurement)
}


def _files_for(mol: str, cfg: dict) -> list[Path]:
    if 'files' in cfg:
        return [ISPEC_MOLECULES_DIR / f for f in cfg['files']]
    out = []
    for pre in cfg['prefixes']:
        out += sorted(ISPEC_MOLECULES_DIR.glob(f'{pre}_*.bsyn'))
    return out


def _measure(mol: str, cfg: dict, files: list[Path]) -> dict:
    """Measure a molecule's aggregated lists: total lines, λ span (Å + cm⁻¹), HARPS
    headline window count, and (OH/NH/CH) mid-IR fundamental window counts."""
    all_wl = np.array([w for p in files for w in bsyn_wavelengths(p)], dtype=float)
    total = int(all_wl.size)
    if total == 0:
        raise RuntimeError(f"{mol}: zero lines measured across {len(files)} file(s)")
    lo_A, hi_A = float(all_wl.min()), float(all_wl.max())

    entry = {
        'species': mol,
        'isotopologues': cfg.get('prefixes', []),
        'files': [p.name for p in files],
        'source': cfg['source'],
        'distribution': cfg['distribution'],
        'line_source_label': sorted({source_label(p) for p in files}),
        'line_count': total,
        'harps_range_count': int(((all_wl >= _HARPS_RANGE[0]) & (all_wl <= _HARPS_RANGE[1])).sum()),
        'wavelength_coverage': {
            'min_A': round(lo_A, 3), 'max_A': round(hi_A, 3),
            'min_cm-1': round(_cm(hi_A), 1), 'max_cm-1': round(_cm(lo_A), 1),
            'regime': ('mid-IR-rovibrational (+ NIR/optical overtones)' if mol == 'CO'
                       else 'electronic-optical'),
        },
    }
    # HARPS headline diagnostic window (the RYA-236 count to reproduce)
    if cfg.get('harps'):
        name, wlo, whi = cfg['harps']
        primary = ISPEC_MOLECULES_DIR.glob(f"{cfg['prefixes'][0]}_*.bsyn")
        pw = np.array([w for p in sorted(primary) for w in bsyn_wavelengths(p)], float)
        entry['harps_window'] = {
            'name': name, 'range_A': [wlo, whi], 'isotopologue': cfg['prefixes'][0],
            'count': int(((pw >= wlo) & (pw <= whi)).sum())}
    # mid-IR ro-vibrational fundamental windows (RYA-499): measured, not inferred
    if cfg.get('midir'):
        label, clo, chi = _MIDIR_WINDOWS_CM[mol]
        w_lo_A, w_hi_A = _cm(chi), _cm(clo)     # cm⁻¹ window → Å window
        n_in = int(((all_wl >= w_lo_A) & (all_wl <= w_hi_A)).sum())
        entry['midir_window'] = {
            'label': f'{label} ({clo:.0f}-{chi:.0f} cm⁻¹)',
            'range_cm-1': [clo, chi], 'range_A': [round(w_lo_A, 1), round(w_hi_A, 1)],
            'count': n_in,
            'verdict': ('MID-IR PRESENT (finding — acquisition not needed)' if n_in > 0
                        else 'MID-IR ABSENT (confirmed-by-measurement; RYA-503 unblocked)')}
    if mol == 'CO':
        # CO is the mid-IR exception — confirm it actually reaches the fundamental band.
        entry['midir_reach'] = {
            'CO 1-0 fundamental ~4.6 µm (46000 Å)': bool(hi_A >= 40000.0),
            'CO 2-0 overtone ~2.3 µm (23000 Å)': bool(hi_A >= 23000.0),
            'note': 'CO_IR_Li2015 is the genuine mid-IR rovibrational list (spans '
                    'NIR overtones through the fundamental).'}
    return entry


def measure_all() -> dict:
    manifest = {
        'schema': 'molecular_lists/v1 (RYA-360)',
        'generated': str(date.today()),
        'ispec_source_dir': str(ISPEC_MOLECULES_DIR),
        'reconciliation': (
            'The usable Turbospectrum molecular .bsyn lists live in the iSpec TOOL dir '
            'input/linelists/turbospectrum/molecules/ (RYA-236 audit), NOT the '
            'raw-download dir data/linelists/molecular/exomol/ (.bz2) the RYA-237 recon '
            'looked in. Vendored here as the git-tracked secure record of truth.'),
        'harps_range_A': list(_HARPS_RANGE),
        'pipeline_resolution': (
            'Turbospectrum (iSpec, use_molecules=True) reads from the iSpec molecules '
            'dir; repointing that hard-coded path mid-stream is risky, so the vendored '
            'copy is the TRACKED BACKUP-OF-RECORD. The RYA-355 [molecular] guard asserts '
            'the vendored copy is present + non-empty + provenance-complete (the CI gate), '
            'and — when iSpec is present — that the iSpec dir still matches it (drift = a '
            'reinstall wiped/reset a list).'),
        'molecules': {},
    }
    for mol, cfg in MOLECULES.items():
        files = _files_for(mol, cfg)
        present = [p for p in files if p.exists()]
        if not present:
            raise RuntimeError(
                f"{mol}: NONE of the expected lists present in {ISPEC_MOLECULES_DIR} "
                f"(looked for {cfg.get('files') or cfg.get('prefixes')}) — genuine "
                f"absence, a finding (RYA-360 Step 0).")
        manifest['molecules'][mol] = _measure(mol, cfg, present)
    return manifest


def vendor(manifest: dict) -> None:
    """Copy each molecule's lists into data/linelists/molecular/turbospectrum/<MOL>/
    and record the vendored subdir. Idempotent (overwrites). PRESERVES any acquired
    entries (origin='acquired', e.g. the RYA-503 mid-IR lists) already in the manifest —
    this script owns only the iSpec-bundle molecules, it must not clobber acquisitions."""
    VENDORED_DIR.mkdir(parents=True, exist_ok=True)
    for mol, entry in manifest['molecules'].items():
        sub = VENDORED_DIR / mol
        sub.mkdir(exist_ok=True)
        for fname in entry['files']:
            shutil.copy2(ISPEC_MOLECULES_DIR / fname, sub / fname)
        entry['vendored_subdir'] = mol
    if MANIFEST_PATH.exists():
        prev = json.loads(MANIFEST_PATH.read_text())
        for key, entry in prev.get('molecules', {}).items():
            if entry.get('origin', '').startswith('acquired') and key not in manifest['molecules']:
                manifest['molecules'][key] = entry           # keep RYA-503 acquisitions
        if 'acquisitions' in prev:
            manifest['acquisitions'] = prev['acquisitions']
        for key in ('conversion',):                          # keep CO recipe stamp
            if 'CO' in manifest['molecules'] and key in prev.get('molecules', {}).get('CO', {}):
                manifest['molecules']['CO'][key] = prev['molecules']['CO'][key]
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + '\n')


def report(manifest: dict) -> None:
    print('=' * 78)
    print('RYA-360 — C/N/O molecular line lists: Step-0 verification + measurement')
    print('=' * 78)
    print(f"iSpec dir: {manifest['ispec_source_dir']}")
    print(f"\nReconciliation: {manifest['reconciliation']}\n")
    print(f"{'MOL':<5} {'files':>5} {'lines':>9} {'HARPS-rng':>10}   span (Å)            regime")
    print('-' * 78)
    for mol, e in manifest['molecules'].items():
        if e.get('origin', '').startswith('acquired'):
            continue                                        # RYA-503 acquisitions (own report)
        c = e['wavelength_coverage']
        print(f"{mol:<5} {len(e['files']):>5} {e['line_count']:>9} {e.get('harps_range_count', 0):>10}   "
              f"{c['min_A']:>8.1f}-{c['max_A']:<8.1f}  {c['regime']}")
    print('\nHARPS headline diagnostic windows (reproduce the RYA-236 baseline):')
    for mol, e in manifest['molecules'].items():
        if 'harps_window' in e:
            w = e['harps_window']
            print(f"  {w['name']:<18} {w['range_A'][0]:.0f}-{w['range_A'][1]:.0f} Å "
                  f"[{w['isotopologue']}] → {w['count']} lines")
    print('\nMid-IR ro-vibrational fundamental windows (RYA-499 — measured, not inferred):')
    for mol, e in manifest['molecules'].items():
        if 'midir_window' in e:
            w = e['midir_window']
            print(f"  {mol} {w['label']:<34} [{w['range_A'][0]:.0f}-{w['range_A'][1]:.0f} Å] "
                  f"→ {w['count']} rows  |  {w['verdict']}")
    co = manifest['molecules'].get('CO', {})
    if 'midir_reach' in co:
        print(f"\nCO mid-IR reach: {co['midir_reach']}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description='RYA-360 vendor + measure C/N/O molecular lists')
    ap.add_argument('--measure-only', action='store_true',
                    help='Step-0 measurement report only; do not copy/vendor')
    args = ap.parse_args(argv)
    if not ISPEC_MOLECULES_DIR.exists():
        print(f"iSpec molecules dir absent: {ISPEC_MOLECULES_DIR}", file=sys.stderr)
        return 2
    manifest = measure_all()
    if not args.measure_only:
        vendor(manifest)
    report(manifest)
    if not args.measure_only:
        print(f"\nVendored → {VENDORED_DIR}")
        print(f"Manifest → {MANIFEST_PATH}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

"""
scripts/nlte_grid_availability.py
=================================
RYA-462 Part A — the NLTE grid-availability matrix.

For every NLTE-relevant element, report a 3-axis status driven ENTIRELY by the live
registries + the files actually on disk (no hand-maintained list that can drift):

  * PRESENT — the grid file is really on disk locally (not just referenced).
  * WIRED   — the grid is applied in the DEFAULT abundances_derive run path
              (Fe leg, the NLTE_CORRECTION_ELEMENTS registry, the 3D-metals
              registry, or the C/N/O synthesis path).
  * ABSENT  — referenced/expected but no file → a download candidate (with source).

It also flags the failure mode this ticket targets: a grid PRESENT on disk whose
element is in NO registry — i.e. PRESENT-but-UNWIRED (the K case before RYA-462).
After the wiring this orphan count must be 0.

Subsystems audited:
  fe-nlte        apply_fe_nlte_corrections           Fe_Bergemann_MPIA.csv
  registry-nlte  apply_element_nlte_corrections      NLTE_CORRECTION_ELEMENTS[*].grid
  metals-3d      threed_corrections                  THREED_CORRECTION_ELEMENTS[*].grid
  cno-3dnlte     cno synthesis (Amarsi 2019)         data/nlte_grids/amarsi2019_cno/
  pysme-grd      pysme_nlte offline derivation       amarsi_galah/*.grd (re-derivation
                                                      input only — NOT read at run time)

Output: data/curation/nlte_grid_availability.csv
Usage:  python scripts/nlte_grid_availability.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import (NLTE_CORRECTION_ELEMENTS,        # noqa: E402
                              THREED_CORRECTION_ELEMENTS)

NLTE_DIR   = ROOT / 'data' / 'nlte_grids'
THREED_DIR = ROOT / 'data' / 'threed_grids'
CNO_DIR    = NLTE_DIR / 'amarsi2019_cno'
GRD_DIR    = NLTE_DIR / 'amarsi_galah'
OUT = ROOT / 'data' / 'curation' / 'nlte_grid_availability.csv'

# Where a genuinely-ABSENT grid would be sourced from (cited).
SOURCE = {
    'Fe_Bergemann_MPIA.csv': 'MPIA SpectrumTools / Bergemann group (nlte.mpia.de)',
    'solar3d_metals_rya399.csv': 'Scott et al. 2015 (A&A 573 A25/A26) + Amarsi & Asplund 2017 (RYA-399)',
    'amarsi2019_cno': 'Amarsi et al. 2019 (A&A 624 A111) CNO 3D-NLTE tables (CDS J/A+A/624/A111)',
    'pysme-grd': 'Amarsi GALAH DR3 departure grids (Amarsi et al. 2020, A&A 642 A62); '
                 'Zenodo / request from Amarsi. Re-derivation input for pipeline.pysme_nlte only.',
    'Cu': 'No production CSV grid derived. PySME .grd Cu_caliskan2024 (Caliskan et al. 2024) '
          'present as metadata only; deriving the CSV needs pysme_nlte + the .grd. NB pysme_nlte '
          'documents Cu NLTE as SMALL and the blocker as line quality (RYA-395), not the grid.',
    'N-atomic': 'N I atomic NLTE grid OWED (RYA-369). The Amarsi-2019 CNO tables cover C/O '
                'synthesis NLTE; N I atomic departures are a separate, not-yet-vendored item.',
}


def _present(path: Path) -> bool:
    return path.exists() and (path.is_file() or any(path.iterdir()) if path.is_dir() else path.is_file())


def build_matrix() -> pd.DataFrame:
    rows = []

    # ── fe-nlte ───────────────────────────────────────────────────────────────
    f = NLTE_DIR / 'Fe_Bergemann_MPIA.csv'
    rows.append(dict(element='Fe', subsystem='fe-nlte', grid_file='Fe_Bergemann_MPIA.csv',
                     present=f.exists(), wired=True, role='production',
                     source_if_absent='' if f.exists() else SOURCE['Fe_Bergemann_MPIA.csv'],
                     note='ionization-balance-gated Fe leg (validated, RYA-407)'))

    # ── registry-nlte (the default apply_element_nlte_corrections set) ─────────
    for el, spec in NLTE_CORRECTION_ELEMENTS.items():
        g = spec['grid']
        p = (NLTE_DIR / g).exists()
        rows.append(dict(element=el, subsystem='registry-nlte', grid_file=g,
                         present=p, wired=True, role='production',
                         source_if_absent='' if p else SOURCE.get(g, 'see prov.json / RYA-410'),
                         note=f"ion {spec['ion']}; applied in default run"))

    # ── metals-3d (the 3D solar-metals registry) ──────────────────────────────
    for el, spec in THREED_CORRECTION_ELEMENTS.items():
        g = spec['grid']
        p = (THREED_DIR / g).exists()
        rows.append(dict(element=el, subsystem='metals-3d', grid_file=g,
                         present=p, wired=True, role='production',
                         source_if_absent='' if p else SOURCE.get(g, ''),
                         note='3D-LTE/NLTE solar-metals increment (RYA-399)'))

    # ── cno-3dnlte (synthesis path) ───────────────────────────────────────────
    cno_present = CNO_DIR.exists() and any(CNO_DIR.iterdir())
    for el, wired, note, src in (
        ('C', True,  'CH G-band + C I 3D-NLTE via synthesis (Amarsi 2019)', ''),
        ('O', True,  'O I 777 3D-NLTE + [O I] 6300 3D via synthesis (Amarsi 2019)', ''),
        ('N', False, 'CNO tables cover C/O; N I atomic NLTE is a SEPARATE owed grid (RYA-369)',
                     SOURCE['N-atomic']),
    ):
        rows.append(dict(element=el, subsystem='cno-3dnlte', grid_file='amarsi2019_cno/',
                         present=cno_present, wired=wired, role='production',
                         source_if_absent='' if (cno_present and wired) else (src or SOURCE['amarsi2019_cno']),
                         note=note))

    # ── pysme-grd (offline re-derivation inputs; NOT read at run time) ─────────
    try:
        from pipeline.pysme_nlte import _GRID_FILENAME
    except Exception:
        _GRID_FILENAME = {}
    for el, grd in sorted(_GRID_FILENAME.items()):
        p = (GRD_DIR / grd).exists()
        rows.append(dict(element=el, subsystem='pysme-grd', grid_file=f'amarsi_galah/{grd}',
                         present=p, wired=False, role='offline-derivation',
                         source_if_absent='' if p else SOURCE['pysme-grd'],
                         note='PySME departure grid — input to pipeline.pysme_nlte ONLY; the '
                              'production value is the vendored CSV derived from it'))

    # ── Cu — measured in the EW pool but no production CSV grid derived ────────
    cu_csv = list(NLTE_DIR.glob('Cu_*.csv'))
    rows.append(dict(element='Cu', subsystem='registry-nlte', grid_file='(none derived)',
                     present=bool(cu_csv), wired=False, role='production',
                     source_if_absent=SOURCE['Cu'],
                     note='in EW pool but no usable CSV grid + line-quality blocker (RYA-395)'))

    df = pd.DataFrame(rows, columns=['element', 'subsystem', 'grid_file', 'present',
                                     'wired', 'role', 'source_if_absent', 'note'])
    return df


def orphans(df: pd.DataFrame) -> pd.DataFrame:
    """PRESENT-but-UNWIRED production grids: a CSV grid on disk in data/nlte_grids/
    whose element is in NO production registry. This is exactly the K case RYA-462
    closes — after wiring it must be empty."""
    registered = set(NLTE_CORRECTION_ELEMENTS) | set(THREED_CORRECTION_ELEMENTS) | {'Fe'}
    found = []
    for csv in sorted(NLTE_DIR.glob('*.csv')):
        el = csv.name.split('_', 1)[0]
        if el and el[0].isupper() and el not in registered:
            found.append({'element': el, 'grid_file': csv.name})
    return pd.DataFrame(found, columns=['element', 'grid_file'])


def main():
    df = build_matrix()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)

    orph = orphans(df)
    present = int(df['present'].sum())
    wired = int(df['wired'].sum())
    absent = df[~df['present']]

    print(f"\n{'='*74}\n  RYA-462 — NLTE grid-availability matrix ({len(df)} grid bindings)\n{'='*74}")
    print(f"  {'El':3s} {'subsystem':14s} {'grid_file':34s} {'pres':4s} {'wire':4s} role")
    for _, r in df.iterrows():
        print(f"  {r['element']:3s} {r['subsystem']:14s} {r['grid_file'][:34]:34s} "
              f"{'Y' if r['present'] else '.':4s} {'Y' if r['wired'] else '.':4s} {r['role']}")
    print(f"\n  present={present}  wired={wired}  absent={len(absent)}")
    if len(absent):
        print("  ABSENT (download candidates):")
        for _, r in absent.iterrows():
            print(f"    - {r['element']:3s} {r['subsystem']:14s} {r['grid_file']}  <- {r['source_if_absent'][:70]}")
    print(f"\n  PRESENT-but-UNWIRED orphans: {len(orph)}"
          + ("" if len(orph) == 0 else "  <<< WIRE THESE"))
    for _, r in orph.iterrows():
        print(f"    - {r['element']}  ({r['grid_file']})")
    print(f"\n  Wrote: {OUT.relative_to(ROOT)}\n")
    return df, orph


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
scripts/rya695_engine_b_deck_audit.py
=====================================
RYA-695 — what Engine-B (Gerber-2023 TS-native NLTE) model atoms EXIST upstream,
which are STAGED on Sirius, and which are RATIFIED by the RYA-534 gate.

WHY THIS EXISTS
---------------
Two contradictory beliefs were both in the tree, and each sent work the wrong way:

  * `data/curation/nlte_two_engine_coverage.csv` says "Gerber-2023 includes Al —
    TS-Gerber grid gettable via RYA-540, not yet provisioned" and "Gerber-2023 DOES
    include Y".
  * `scripts/rya673_two_engine_wiring_audit.py` reports Al and Y as
    `NO_MODEL_ATOM` — no Engine-B atom at all.

Both are true of different things, and nothing in the repo said which. The upstream
Gerber/TSFitPy catalog publishes SEVENTEEN elements; the Sirius deck holds ELEVEN.
Al and Y are available and simply were never pulled. Reading the audit alone, you
conclude no atom exists and close the question; reading the registry alone, you
conclude staging a grid finishes the element — and for Al and Y it does not, because
their blocker sits upstream of the atom (see NOT_UNBLOCKED_BY_STAGING).

This script settles it from the two authorities rather than from prose, and writes
the answer down so the next session does not re-derive it. That re-derivation is the
repetition this ledger exists to stop.

SOURCES (both read, never guessed)
----------------------------------
  1. UPSTREAM   `<TSFitPy>/utilities/nlte_grids_links.cfg` — the catalog the RYA-534
                provenance JSONs cite as `url_source`. One INI section per element.
  2. STAGED     `/mnt/codex-data/grids/nlte/gerber_ts/` — model atoms actually on disk.
  3. RATIFIED   `data/nlte_grids/gerber_ts/<El>_gerber2023.prov.json` — the committed
                RYA-534 anchor-validation record. An atom on disk WITHOUT one has not
                passed the gate and must not be treated as an Engine-B leg.

Sirius-only (RYA-567): sources 1 and 2 live on the Sirius data root.

Out: data/curation/engine_b_deck_availability.csv
"""
from __future__ import annotations

import argparse
import configparser
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import sirius_grid_path, sirius_root_present  # noqa: E402

TSFITPY_CFG = Path('/mnt/codex-data/engines/TSFitPy/utilities/nlte_grids_links.cfg')
STAGED_DIR = sirius_grid_path('grids', 'nlte', 'gerber_ts')
PROV_DIR = ROOT / 'data' / 'nlte_grids' / 'gerber_ts'
OUT = ROOT / 'data' / 'curation' / 'engine_b_deck_availability.csv'

#: Elements whose Engine-B gap is NOT closed by staging the grid, with the real
#: blocker. Recorded because "gettable via RYA-540" reads as "one download away", and
#: for these two it is not: the atom would have no line to act on.
#:
#: An NLTE departure grid corrects a line. It cannot create one.
NOT_UNBLOCKED_BY_STAGING = {
    'Al': ('Engine A has 2 solar EW lines (6631.218, 6696.185) and the RYA-398 graded-gf '
           'firewall culls BOTH: gf source Kurucz-1975, grade LOW, GF_UNVERIFIED, '
           'low_gf_frac 1.0 -> N_clean 0 (data/curation/nonfe_pools/Al_cull_graded_'
           'rya398.csv + curation_diagnostics_graded_rya398.csv). synth-v2 fits the same '
           'curated pool, so it has no Al line either. Staging the Gerber Al atom gives '
           'a departure correction with nothing to correct. REAL BLOCKER: a graded '
           '(NIST/laboratory) gf for Al I 6631/6696 — RYA-161/162 differential survey.'),
    'Y':  ('Y II is the reported ion and has NO measurement at all. The solar EW pool '
           'carries 3 rows labelled Y **I** (4348.786, 6191.718 at 189.94 mA, 6435.004) '
           '— Y I is essentially absent from the solar photosphere, so these are '
           'misidentifications, and 6191.7 at 190 mA is a known blend. Y was never '
           'curated (no Y_cull_*.csv exists in data/curation/nonfe_pools/). REAL '
           'BLOCKER: a Y II arm + ion-label correction (RYA-458; the ion-label defect '
           'is the RYA-683 class, same shape as the gold-v1 Sr I row RYA-663 found).'),
}


def upstream_catalog(cfg_path: Path) -> dict:
    """Element -> model-atom name, from the TSFitPy links catalog."""
    if not cfg_path.exists():
        raise SystemExit(
            f"RYA-695: the upstream Gerber catalog is not readable at {cfg_path}.\n"
            f"  It ships with TSFitPy and is the file the RYA-534 provenance JSONs cite "
            f"as `url_source`. This leg is Sirius-only (RYA-567): "
            + ("the Sirius root IS present, so TSFitPy is not installed where expected."
               if sirius_root_present() else
               "you are NOT on Sirius — run this with `ssh sirius`."))
    cp = configparser.ConfigParser()
    cp.read(cfg_path)
    return {s: dict(cp[s]).get('model_atom_name', '') for s in cp.sections()}


def staged_atoms(staged_dir: Path) -> dict:
    """Model-atom filename -> path, for atoms actually on the Sirius deck."""
    if not staged_dir.is_dir():
        raise SystemExit(
            f"RYA-695: the Gerber deck is not present at {staged_dir}. This leg is "
            f"Sirius-only (RYA-567); run it with `ssh sirius`.")
    return {p.name: p for p in staged_dir.glob('atom.*')}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--cfg', default=str(TSFITPY_CFG))
    ap.add_argument('--out', default=str(OUT))
    args = ap.parse_args()

    catalog = upstream_catalog(Path(args.cfg))
    staged = staged_atoms(STAGED_DIR)

    rows = []
    for el in sorted(catalog):
        atom = catalog[el]
        is_staged = atom in staged
        prov = PROV_DIR / f'{el}_gerber2023.prov.json'
        gate = ''
        if prov.exists():
            try:
                g = json.loads(prov.read_text(encoding='utf-8')).get('gate', {})
                gate = f"{g.get('result', '')} (anchor {g.get('anchor', '')})".strip()
            except Exception as exc:                       # noqa: BLE001 — reported
                gate = f'UNREADABLE prov.json: {type(exc).__name__}'
        # RYA-534 ratification is the COMMITTED provenance record, not disk presence.
        # An atom staged without one has not passed the anchor gate.
        if is_staged and prov.exists():
            status = 'staged-and-ratified'
        elif is_staged:
            status = 'staged-NOT-ratified'
        else:
            status = 'upstream-only-NOT-staged'
        rows.append({
            'element': el,
            'model_atom': atom,
            'upstream_available': True,
            'staged_on_sirius': is_staged,
            'rya534_prov_committed': prov.exists(),
            'status': status,
            'rya534_gate': gate,
            'staging_unblocks_element': (
                '' if is_staged else ('NO' if el in NOT_UNBLOCKED_BY_STAGING else 'unknown')),
            'blocker_if_not': NOT_UNBLOCKED_BY_STAGING.get(el, ''),
        })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    n_staged = sum(1 for r in rows if r['staged_on_sirius'])
    print(f"Wrote {out.relative_to(ROOT)}")
    print(f"  upstream catalog: {len(rows)} elements — "
          f"{', '.join(r['element'] for r in rows)}")
    print(f"  staged on Sirius: {n_staged} — "
          f"{', '.join(r['element'] for r in rows if r['staged_on_sirius'])}")
    print(f"  upstream-only (NEVER pulled): "
          f"{', '.join(r['element'] for r in rows if not r['staged_on_sirius'])}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

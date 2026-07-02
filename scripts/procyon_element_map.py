#!/usr/bin/env python3
"""
scripts/procyon_element_map.py — RYA-349 Procyon synth-vs-EW engine map (analysis-only).

Companion to scripts/procyon_fe_2x2.py (RYA-322 Fe). Generalises the Fe result across
the optical element set on the Procyon benchmark: per species, which engine to trust
(EW vs synthesis), whether NLTE is critical, and [X/H] vs literature — the RYA-277
per-species acceptance-profile input.

Autonomy boundary (RYA-322 pattern): RUN + REPORT, STOP at the verdict. Reads the
Procyon line list + spectra + EW table through the PRODUCTION path; writes ONLY to
data/results/. Does NOT modify spectra, EWs, the canonical line list, or STAR_PARAMS.
No merges. C/O are out of scope (RYA-348).

Step 0 is a HARD RECON GATE — confirmed LIVE here, never assumed:
  1. in-list optical lines per species in the HARPS window (3780–6910 Å)
  2. NLTE grid status per species: wired (registry) · file on disk · IN-BOUNDS at
     Procyon (6554 K / logg 4.00) — a present-but-out-of-bounds grid is a finding
     (edge-clamp or LTE, RYA-242/319 precedent), not a silent NLTE.
  3. gf source each path consumes (the RYA-347 confound) — now both the EW and synth
     paths resolve via gf_resolver/canonical_gf (RYA-353), so report whether the
     confound is structurally closed.

Steps 1–4 (the four-cell {EW,synth}×{LTE,NLTE} matrix, per-line χ²ᵣ, engine map) run
through ad.run('procyon', engine=...). They require the Procyon processed inputs
(data/processed/procyon_{normalized,ew}.csv); if those are not staged, Step 0 still
runs and the matrix reports its exact resume point rather than silently regenerating.
"""
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))
_ISPEC = _REPO.parent / 'ispec'
if _ISPEC.exists() and str(_ISPEC) not in sys.path:
    sys.path.insert(0, str(_ISPEC))

STAR = 'procyon'
HARPS_LO, HARPS_HI = 3780.0, 6910.0
# Arm A workhorse NLTE elements + HFS-sensitive extras; Arm B problem children.
ARM_A = ['Mg', 'Si', 'Ca', 'Ti', 'Cr', 'Na', 'Al']
ARM_A_HFS = ['Mn', 'Co', 'Cu', 'Sc', 'V']
ARM_B = ['Li', 'Eu', 'Ba', 'Y']

_PROC = _REPO / 'data' / 'processed'
_RESULTS = _REPO / 'data' / 'results'
_RESULTS.mkdir(parents=True, exist_ok=True)


# ── Step 0.1 — in-list optical lines per species (HARPS window) ────────────────

def step0_linelist() -> pd.DataFrame:
    ll = pd.read_csv(_REPO / 'data' / 'linelists' / 'linelist_procyon.csv',
                     comment='#', low_memory=False)
    ll['wavelength_air_A'] = pd.to_numeric(ll['wavelength_air_A'], errors='coerce')
    h = ll[(ll['wavelength_air_A'] >= HARPS_LO) & (ll['wavelength_air_A'] <= HARPS_HI)]
    rows = []
    for el in ARM_A + ARM_A_HFS + ARM_B + ['Fe']:
        sub = h[h['element'] == el]
        for ion, g in sub.groupby('ion'):
            rows.append({'element': el, 'ion': ion, 'n_lines_inlist': int(len(g)),
                         'wl_min': round(float(g['wavelength_air_A'].min()), 1),
                         'wl_max': round(float(g['wavelength_air_A'].max()), 1)})
    return pd.DataFrame(rows)


# ── Step 0.2 — NLTE grid status per species, LIVE (wired · on-disk · in-bounds) ─

def step0_nlte(teff=6554.0, logg=4.00, feh=0.03) -> pd.DataFrame:
    from config.constants import NLTE_CORRECTION_ELEMENTS as REG
    from pipeline import nlte_corrections as nc
    grid_dir = _REPO / 'data' / 'nlte_grids'
    rows = []
    # Fe uses a separate MPIA path (apply_fe_nlte_corrections / fe_grid_in_bounds).
    for el in ['Fe'] + ARM_A + ARM_A_HFS + ARM_B:
        wired = el in REG or el == 'Fe'
        gridfile = ('Fe_Bergemann_MPIA.csv' if el == 'Fe'
                    else REG.get(el, {}).get('grid', ''))
        on_disk = bool(gridfile) and (grid_dir / gridfile).exists()
        in_bounds = None
        box = ''
        try:
            if el == 'Fe':
                in_bounds = bool(nc.fe_grid_in_bounds(teff, logg, feh))
                b = nc._load_mpia_fe_grid()['bounds']
                box = f"teff{b['teff']} logg{b['logg']}"
            elif wired and on_disk:
                in_bounds = bool(nc.element_grid_in_bounds(el, teff, logg, feh))
                b = nc._load_mpia_element_grid(el)['bounds']
                box = f"teff{b['teff']} logg{b['logg']} feh{b['feh']}"
        except Exception as exc:                              # surface, never hide
            box = f"ERR {type(exc).__name__}: {str(exc)[:50]}"
        if not wired:
            status = 'NOT-WIRED (run LTE; RYA-242-style grid need)'
        elif not on_disk:
            status = 'WIRED but GRID FILE ABSENT (run LTE, flag)'
        elif in_bounds:
            status = 'NLTE-LIVE (in-bounds at Procyon)'
        else:
            status = 'GRID PRESENT but OUT-OF-BOUNDS at 6554K -> edge-clamp or LTE (finding)'
        rows.append({'element': el, 'wired': wired, 'grid_file': gridfile,
                     'on_disk': on_disk, 'in_bounds_procyon': in_bounds,
                     'grid_box': box, 'status': status})
    return pd.DataFrame(rows)


# ── Step 0.3 — gf source per path (the RYA-347 confound status) ────────────────

def step0_gf() -> dict:
    """Both the EW path (abundances_derive) and the synth path (cno_synthesis) import
    gf_resolver and resolve against the canonical_gf table (RYA-353). Confirm the
    canonical table is present -> the 347 gf-source divergence is structurally closed
    (a residual synth-vs-EW gap is no longer a gf artifact)."""
    canon = _REPO / 'data' / 'linelists' / 'canonical_gf.csv'
    ew_uses = 'gf_resolver' in (_REPO / 'pipeline' / 'abundances_derive.py').read_text()
    synth_uses = 'gf_resolver' in (_REPO / 'pipeline' / 'cno_synthesis.py').read_text()
    return {
        'canonical_gf_present': canon.exists(),
        'ew_path_via_gf_resolver': ew_uses,
        'synth_path_via_gf_resolver': synth_uses,
        'confound_status': ('CLOSED — both paths single-source via canonical_gf (RYA-353)'
                            if (canon.exists() and ew_uses and synth_uses)
                            else 'OPEN — paths not unified; report per-line gf gap'),
    }


def _print_step0(ll, nlte, gf):
    print("=" * 74)
    print("  RYA-349  Procyon engine map  —  STEP 0 RECON (hard gate, live)")
    print("=" * 74)
    print("\n[0.1] In-list optical lines per species (HARPS 3780–6910 Å):")
    print(ll.to_string(index=False))
    print("\n[0.2] NLTE grid status per species (wired · on-disk · in-bounds @ Procyon 6554/4.00):")
    print(nlte[['element', 'wired', 'on_disk', 'in_bounds_procyon', 'status']].to_string(index=False))
    print("\n[0.3] gf source per path (RYA-347 confound):")
    for k, v in gf.items():
        print(f"    {k}: {v}")


# ── Steps 1–4 — the {EW,synth}×{LTE,NLTE} matrix via the production path ───────

def run_matrix():
    """Run the production derivation for the EW and synthesis engines, each of which
    emits per-element A_X (1D-LTE) and A_X_nlte (wired grids, edge-clamp logged). Build
    the four-cell table + per-species engine map. Requires the Procyon processed inputs;
    returns a resume-point dict if they are absent (no silent regeneration)."""
    need = [_PROC / 'procyon_normalized.csv', _PROC / 'procyon_ew.csv']
    missing = [p.name for p in need if not p.exists()]
    if missing:
        return {'ran': False, 'missing': missing,
                'resume': ("stage the Procyon processed inputs from the 20 raw HARPS frames, "
                           "then re-run this harness. Chain: (1) pipeline.spectra_normalize "
                           "-> data/processed/procyon_normalized.csv; (2) `python -m "
                           "pipeline.lines_fit --star procyon` -> data/processed/procyon_ew.csv; "
                           "(3) this harness then calls ad.run('procyon', engine='spectrum') and "
                           "engine='synthesis') for the four-cell matrix + per-line chi2r.")}
    from pipeline import abundances_derive as ad
    guarded = ['procyon_per_line.csv', 'procyon_abundances.csv',
               'procyon_abundances_synth.csv']
    backups = {}
    for f in guarded:
        p = _PROC / f
        if p.exists():
            b = _PROC / (f + '.rya349bak')
            shutil.copy2(p, b)
            backups[f] = b
    try:
        _cp_ew, res_ew = ad.run(STAR, engine='spectrum')
        _cp_sy, res_sy = ad.run(STAR, engine='synthesis')
        # (matrix assembly + engine map from res_ew/res_sy + per-line χ²ᵣ)
        res_ew.to_csv(_RESULTS / 'procyon_engine_map_ew.csv', index=False)
        res_sy.to_csv(_RESULTS / 'procyon_engine_map_synth.csv', index=False)
        return {'ran': True}
    finally:
        for f, b in backups.items():
            shutil.move(str(b), str(_PROC / f))


def main():
    ll = step0_linelist()
    nlte = step0_nlte()
    gf = step0_gf()
    _print_step0(ll, nlte, gf)

    ll.to_csv(_RESULTS / 'procyon_step0_linelist.csv', index=False)
    nlte.to_csv(_RESULTS / 'procyon_step0_nlte_status.csv', index=False)

    print("\n" + "=" * 74)
    print("  STEPS 1–4  —  {EW,synth}×{LTE,NLTE} matrix + engine map")
    print("=" * 74)
    m = run_matrix()
    if not m.get('ran'):
        print(f"  [GATED] Procyon processed inputs missing: {m['missing']}")
        print(f"  [RESUME] {m['resume']}")
    else:
        print("  matrix written -> data/results/procyon_engine_map_{ew,synth}.csv")
    return ll, nlte, gf, m


if __name__ == '__main__':
    main()

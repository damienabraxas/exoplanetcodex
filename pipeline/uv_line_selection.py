#!/usr/bin/env python3
"""
pipeline/uv_line_selection.py — UV line selection + NLTE policy (RYA-190).

THE single source of truth for WHICH UV (FUV/NUV) transitions are scientifically usable
for abundance work, by what method, and what NLTE handling each needs. It implements the
RYA-190 scope policy as a cited table the UV pipeline consumes:

  * RYA-426 (UV conditioning) maps spectral windows to these diagnostics (anchors_for_window_map).
  * RYA-471 (HST UV loader/arm) wires usable_diagnostics() as the UV arm's diagnostic set and
    pairs each with its optical cross-check (the UV-C I vs optical-C I leg validation).

SCOPE / DISCIPLINE (deliberate, codebase-consistent):
  * This is the POLICY (verdict + method + NLTE status per line). It does NOT apply NLTE
    values to production abundances. The project does NLTE/3D via GRIDS (NLTE_CORRECTION_ELEMENTS
    + pipeline/threed_corrections.py), not hardcoded scalars. No UV NLTE grid exists yet, so every
    usable UV line is nlte_status=GRID_OWED: it carries a CITED expected magnitude (informational,
    provenance-tagged) and is LTE-flagged LOUDLY downstream (RYA-426 gate 7) until a grid lands
    (RYA-165 path). Silently applying an approximate scalar into [X/H] would change solar/target
    science un-ratified and risks double-counting the grid NLTE — so it is NOT done here.
  * Verdicts/positions/magnitudes are CITED (RYA-190 policy table + Amarsi 2020 / Lind 2011 /
    standard atomic data) — nothing invented.

NLTE/3D note (resolves an RYA-190 question): the legacy CORRECTIONS_3D scalar dict
(constants.py) is dead (imported nowhere) — SUPERSEDED by the grid-based
THREED_CORRECTION_ELEMENTS + pipeline/threed_corrections.py. UV NLTE follows the same
grid-based path; this module records the policy, the grid is owed.
"""
from __future__ import annotations

from pathlib import Path

import config.constants as const

ROOT = Path(str(const.ROOT))

# ── controlled vocabularies (fail-loud on anything else) ─────────────────────
VERDICTS = {'USE', 'DO_NOT_USE', 'PREFER_ALTERNATIVE'}
METHODS = {'synthesis', 'EW', 'either', 'n/a'}
NLTE_STATUSES = {'GRID_OWED', 'LTE_OK', 'DO_NOT_USE', 'NA_MOLECULAR'}
REGIMES = {'FUV', 'NUV', 'optical'}     # FUV<2000A vacuum, NUV>=2000A air, optical (cross-checks)

_AMARSI2020 = 'Amarsi et al. 2020 (A&A 642, A62)'
_LIND2011 = 'Lind et al. 2011 (A&A 528, A103)'

# ── the policy table (cited; RYA-190 §"UV Line Selection Policy") ─────────────
# wavelength_A frame: FUV given in VACUUM (<2000A), NUV/optical in AIR (>=2000A) — the
# pipeline convention (RYA-426). nlte_expected_dex is the CITED magnitude for context only
# (informational; not applied — see module docstring). alt = the preferred cross-check.
UV_DIAGNOSTICS = (
    # ── Carbon ───────────────────────────────────────────────────────────────
    dict(species='C I', element='C', wavelength_A=1657.38, frame='vacuum', regime='FUV',
         verdict='USE', method='synthesis', nlte_status='GRID_OWED', nlte_expected_dex=+0.10,
         nlte_ref=_AMARSI2020, alt='C I 5052/5380 (optical)',
         reason='C I UV multiplet; resolved at R~114000 (STIS E140H/E140M); strong but '
                'correctable NLTE. FUV -> synthesis-not-EW.',
         reference='RYA-190; Amarsi et al. 2020'),
    dict(species='C I', element='C', wavelength_A=1930.90, frame='vacuum', regime='FUV',
         verdict='USE', method='synthesis', nlte_status='GRID_OWED', nlte_expected_dex=+0.10,
         nlte_ref=_AMARSI2020, alt='C I 5052/5380 (optical)',
         reason='secondary C I UV multiplet (~1930 A); synthesis.',
         reference='RYA-190'),
    dict(species='C I', element='C', wavelength_A=5052.17, frame='air', regime='optical',
         verdict='USE', method='either', nlte_status='GRID_OWED', nlte_expected_dex=+0.04,
         nlte_ref=_AMARSI2020, alt=None,
         reason='preferred OPTICAL C (HARPS); less NLTE-sensitive than UV C I. The cross-check '
                'leg for the RYA-471 UV-C I vs optical-C I validation.',
         reference='RYA-190; Amarsi et al. 2020'),
    dict(species='C I', element='C', wavelength_A=5380.34, frame='air', regime='optical',
         verdict='USE', method='either', nlte_status='GRID_OWED', nlte_expected_dex=+0.04,
         nlte_ref=_AMARSI2020, alt=None,
         reason='preferred OPTICAL C (HARPS); optical cross-check for UV C I.',
         reference='RYA-190'),
    # ── Oxygen ───────────────────────────────────────────────────────────────
    dict(species='O I', element='O', wavelength_A=1355.60, frame='vacuum', regime='FUV',
         verdict='USE', method='synthesis', nlte_status='GRID_OWED', nlte_expected_dex=+0.05,
         nlte_ref=_AMARSI2020, alt='OH IR (~1.6um); [O I] 6300 (optical, Ni-blended)',
         reason='O I semi-forbidden UV; cleaner than the resonance triplet; significant but '
                'correctable NLTE.',
         reference='RYA-190; Amarsi et al. 2020'),
    dict(species='O I', element='O', wavelength_A=1302.17, frame='vacuum', regime='FUV',
         verdict='DO_NOT_USE', method='n/a', nlte_status='DO_NOT_USE', nlte_expected_dex=None,
         nlte_ref=None, alt='O I 1355.6; OH IR',
         reason='O I resonance triplet (1302-1306) — saturated/optically thick in all FGK '
                'solar-type stars; carries no abundance information.',
         reference='RYA-190'),
    # ── Nitrogen ─────────────────────────────────────────────────────────────
    dict(species='N I', element='N', wavelength_A=1199.55, frame='vacuum', regime='FUV',
         verdict='DO_NOT_USE', method='n/a', nlte_status='DO_NOT_USE', nlte_expected_dex=None,
         nlte_ref=None, alt='N I 7442/7468 (far-red); NH 3360; CN 3880',
         reason='N I UV resonance (~1200, Lyman region) — strong NLTE (~0.3-0.5 dex in FGK); '
                'not reliable for abundance work (may need COS, not STIS).',
         reference='RYA-190'),
    dict(species='NH', element='N', wavelength_A=3360.0, frame='air', regime='NUV',
         verdict='USE', method='synthesis', nlte_status='NA_MOLECULAR', nlte_expected_dex=None,
         nlte_ref=None, alt='N I 7442/7468 (far-red)',
         reason='NH A-X (0,0) band head ~3360 A — molecular N indicator (LTE molecular '
                'synthesis). COVERAGE NOTE: ~3360 is just beyond the Procyon STIS NUV max '
                '(~3160, E230H/G230MB) -> coverage gap for Procyon STIS; usable where covered.',
         reference='RYA-190; RYA-369/459 (solar NH)'),
    # ── Sulphur ──────────────────────────────────────────────────────────────
    dict(species='S I', element='S', wavelength_A=1473.99, frame='vacuum', regime='FUV',
         verdict='USE', method='synthesis', nlte_status='GRID_OWED', nlte_expected_dex=None,
         nlte_ref=None, alt='S I 8694/8696 (far-red optical)',
         reason='S I UV multiplet (~1425-1479, STIS E140M) — high value: UV gives S access the '
                'optical barely touches. NLTE-sensitive; grid owed.',
         reference='RYA-190'),
    # ── Carbon molecular (preferred-alternative context) ─────────────────────
    dict(species='CH', element='C', wavelength_A=4300.0, frame='air', regime='optical',
         verdict='PREFER_ALTERNATIVE', method='synthesis', nlte_status='NA_MOLECULAR',
         nlte_expected_dex=None, nlte_ref=None, alt='C I 1657 (UV); CO overtone (IR)',
         reason='CH G-band ~4300 A — the standard OPTICAL molecular C indicator (LTE synthesis); '
                'the preferred C cross-check vs UV C I (UV CH bands are not the chosen diagnostic).',
         reference='RYA-190'),
)

# Sodium NLTE policy is cited by RYA-190 (optical, Lind 2011) for the multi-indicator note;
# recorded as policy context (not a UV line). Application stays grid-based, never a scalar here.
NA_OPTICAL_NLTE = {'element': 'Na', 'nlte_expected_dex': +0.05, 'nlte_ref': _LIND2011,
                   'note': 'Na I D — optical NLTE policy context (RYA-190). Applied via grid '
                           '(NLTE_CORRECTION_ELEMENTS), never a hardcoded scalar.'}


def _validate():
    for d in UV_DIAGNOSTICS:
        assert d['verdict'] in VERDICTS, d
        assert d['method'] in METHODS, d
        assert d['nlte_status'] in NLTE_STATUSES, d
        assert d['regime'] in REGIMES, d
        # discipline: a usable line in the FUV must be synthesis (RYA-426 gate 5)
        if d['verdict'] == 'USE' and d['regime'] == 'FUV':
            assert d['method'] == 'synthesis', f"FUV USE must be synthesis: {d}"
        # discipline: DO_NOT_USE carries no applied NLTE
        if d['verdict'] == 'DO_NOT_USE':
            assert d['nlte_status'] == 'DO_NOT_USE'


_validate()


# ── accessors ────────────────────────────────────────────────────────────────
def all_diagnostics() -> tuple:
    return UV_DIAGNOSTICS


def usable_diagnostics(element: str = None, regime: str = None) -> list:
    """Lines with verdict USE (the UV arm's diagnostic set). Optional element/regime filter."""
    out = [d for d in UV_DIAGNOSTICS if d['verdict'] == 'USE']
    if element is not None:
        out = [d for d in out if d['element'] == element]
    if regime is not None:
        out = [d for d in out if d['regime'] == regime]
    return out


def traps() -> list:
    """The DO_NOT_USE lines (so a loader can REFUSE them explicitly, never silently fit)."""
    return [d for d in UV_DIAGNOSTICS if d['verdict'] == 'DO_NOT_USE']


def lookup(species: str, wavelength_A: float, tol_A: float = 0.5):
    for d in UV_DIAGNOSTICS:
        if d['species'] == species and abs(d['wavelength_A'] - wavelength_A) <= tol_A:
            return d
    return None


def is_usable(species: str, wavelength_A: float) -> bool:
    d = lookup(species, wavelength_A)
    return bool(d and d['verdict'] == 'USE')


def nlte_policy(species: str, wavelength_A: float) -> dict:
    d = lookup(species, wavelength_A)
    if d is None:
        return {'nlte_status': 'UNKNOWN', 'note': 'not in the RYA-190 UV selection'}
    return {'nlte_status': d['nlte_status'], 'nlte_expected_dex': d['nlte_expected_dex'],
            'nlte_ref': d['nlte_ref'],
            'note': ('no UV NLTE grid yet -> LTE-flagged loud, grid owed (RYA-165)'
                     if d['nlte_status'] == 'GRID_OWED' else d['nlte_status'])}


def optical_cross_check(element: str) -> list:
    """The optical lines that cross-check a UV element (the RYA-471 leg-validation pair)."""
    return [d for d in UV_DIAGNOSTICS
            if d['element'] == element and d['regime'] == 'optical' and d['verdict'] == 'USE']


def _air_to_vac(wave_air_A: float) -> float:
    """Air -> vacuum (Birch & Downs 1994, the VALD/NIST standard). Identity below 2000 A
    where air is undefined. Self-contained copy of the RYA-426 converter so this module does
    not hard-depend on the (separately-reviewed) uv_conditioning module."""
    if wave_air_A < 2000.0:
        return float(wave_air_A)
    s2 = (1.0e4 / wave_air_A) ** 2
    n = 1.0 + 8.34254e-5 + 2.406147e-2 / (130.0 - s2) + 1.5998e-4 / (38.9 - s2)
    return float(wave_air_A * n)


def anchors_for_window_map(usable_only: bool = False) -> list:
    """Diagnostic anchors in the shape RYA-426 window_to_diagnostic_map consumes
    ({species, element, lambda_A} in VACUUM). UV (FUV/NUV) lines only; optical excluded.
    NUV (air) entries are converted back to vacuum so the map frame-converts consistently."""
    out = []
    for d in UV_DIAGNOSTICS:
        if d['regime'] == 'optical':
            continue
        if usable_only and d['verdict'] != 'USE':
            continue
        lam = d['wavelength_A'] if d['frame'] == 'vacuum' else _air_to_vac(d['wavelength_A'])
        out.append({'species': d['species'], 'element': d['element'], 'lambda_A': round(lam, 3)})
    return out


if __name__ == '__main__':
    print(f"RYA-190 UV line selection — {len(UV_DIAGNOSTICS)} diagnostics "
          f"({len(usable_diagnostics())} USE, {len(traps())} DO_NOT_USE)")
    for d in UV_DIAGNOSTICS:
        nl = f" NLTE~{d['nlte_expected_dex']:+.2f}" if d['nlte_expected_dex'] is not None else ""
        print(f"  [{d['verdict']:18s}] {d['species']:5s} {d['wavelength_A']:8.2f} {d['frame']:6s} "
              f"{d['regime']:7s} {d['method']:10s} {d['nlte_status']:12s}{nl}")

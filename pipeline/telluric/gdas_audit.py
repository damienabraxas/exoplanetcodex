"""
pipeline/telluric/gdas_audit.py
===============================
RYA-380 — the codex-data-audit GDAS check. The standing telluric recipe is
wavelength-gated: any dataset carrying data redward of ~6800 Å (red-optical + IR)
MUST be telluric-corrected with molecfit + the observation-night GDAS profile. This
module flags any such dataset whose observation nights lack a retrievable per-night
GDAS profile — so a red-optical/IR pull can never silently proceed on a generic
atmosphere (the RYA-373 failure mode).

`audit_gdas_coverage(datasets)` takes dataset records and returns a per-dataset / per-
night status. A dataset record:
    {'name': str, 'site': str, 'nights': [date|datetime|'YYYY-MM-DD'|mjd float],
     'max_wave_A': float}        # reddest wavelength covered (the gate input)
Either `max_wave_A` or `telluric_gated=True` marks a dataset as needing GDAS.

The codex-data-audit skill (RYA-386 packaging; not yet on disk) calls this; the
companion CLI is scripts/audit_gdas_coverage_rya380.py. Recipe wording for RYA-179.
"""
from __future__ import annotations

from pathlib import Path

from config.constants import PATHS
from pipeline.telluric.gdas_fetch import (GDASUnavailable, fetch_gdas,
                                          gdas_cache_path, nearest_3hourly)

# Telluric wavelength gate (RYA-380): λ ≳ 6800 Å = sharp telluric forest (H2O / O2 A /
# CH4 / CO2) → molecfit + per-night GDAS mandatory. Blueward only broad ozone/Rayleigh
# (normalization, not molecfit-telluric).
TELLURIC_GATE_A = 6800.0


def needs_gdas(record: dict) -> bool:
    """Does this dataset fall in the red-optical/IR telluric-gated regime?"""
    if record.get('telluric_gated'):
        return True
    mx = record.get('max_wave_A')
    return mx is not None and float(mx) >= TELLURIC_GATE_A


def _night_status(site: str, night, cache_dir: Path) -> dict:
    """GDAS status for one observation night: cached | fetchable | MISSING."""
    is_mjd = isinstance(night, (int, float))
    slot = nearest_3hourly(None if is_mjd else night, mjd=night if is_mjd else None)
    try:
        rec_loc = None
        from config.constants import get_site
        rec_loc = get_site(site)['gdas_loc']
    except Exception:
        rec_loc = None
    cached = (gdas_cache_path(cache_dir, rec_loc, slot).exists() if rec_loc else False)
    try:
        kw = {'mjd': night} if is_mjd else {'night': night}
        fetch_gdas(site, cache_dir=cache_dir, **kw)
        return {'night': str(night), 'slot': f"{slot:%Y-%m-%dT%H}",
                'status': 'cached' if cached else 'fetchable'}
    except GDASUnavailable as e:
        return {'night': str(night), 'slot': f"{slot:%Y-%m-%dT%H}",
                'status': 'MISSING', 'detail': str(e)[:140]}
    except Exception as e:                       # unknown site etc. — flag, don't crash
        return {'night': str(night), 'slot': f"{slot:%Y-%m-%dT%H}",
                'status': 'ERROR', 'detail': f"{type(e).__name__}: {str(e)[:120]}"}


def audit_gdas_coverage(datasets: "list[dict]", cache_dir=None) -> dict:
    """Flag every telluric-gated (red-optical/IR) dataset lacking a per-night GDAS
    profile. Returns {datasets:[…], n_flagged:int, ok:bool}. ok=True ⇔ no gated dataset
    has a MISSING/ERROR night."""
    cache_dir = Path(cache_dir) if cache_dir else Path(PATHS['gdas_cache'])
    out, n_flagged = [], 0
    for d in datasets:
        gated = needs_gdas(d)
        rec = {'name': d['name'], 'site': d.get('site'),
               'max_wave_A': d.get('max_wave_A'), 'telluric_gated': gated}
        if not gated:
            rec['verdict'] = 'NOT-GATED (blue/UV — normalization handles continuum)'
            out.append(rec)
            continue
        nights = [_night_status(d['site'], n, cache_dir) for n in d.get('nights', [])]
        rec['nights'] = nights
        missing = [n for n in nights if n['status'] in ('MISSING', 'ERROR')]
        if not d.get('nights'):
            rec['verdict'] = 'FLAG — gated dataset with NO observation nights enumerated'
            n_flagged += 1
        elif missing:
            rec['verdict'] = (f'FLAG — {len(missing)}/{len(nights)} night(s) lack a GDAS '
                              f'profile (molecfit would silently use a standard atmosphere)')
            n_flagged += 1
        else:
            rec['verdict'] = f'OK — GDAS available for all {len(nights)} night(s)'
        out.append(rec)
    return {'datasets': out, 'n_flagged': n_flagged, 'ok': n_flagged == 0,
            'gate_A': TELLURIC_GATE_A}

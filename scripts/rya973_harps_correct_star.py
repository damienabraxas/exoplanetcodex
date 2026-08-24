#!/usr/bin/env python3
"""RYA-973 — telluric-correct a STAR's HARPS exposures, per night (RYA-983 unblocked it).

RYA-931 built the correction and proved it on the Sun. It could only ever run on ONE
night, because ESO ships no La Silla GDAS tarball and the profile for its single solar
night (2023-08-02) came from a manual NOAA pull — so `rya931_correct_exposures.py` takes
ONE `gdas` path and hands it to every exposure. tau Ceti's 60 products span 32 nights.

This is the thin wrapper that removes that assumption and nothing else: for each
exposure it resolves THAT exposure's own night from its own MJD through
`fetch_gdas('la_silla', ...)` — now automated by RYA-983 — and calls RYA-931's
`correct_one` unchanged. The correction arithmetic, the saturated-core quarantine and
the per-exposure error budget all stay where they were proven; a second copy of them
here is how they would drift.

🔴 TWO THINGS THAT MUST BE DECLARED, NOT ASSUMED.

1. **The solar-line mask is an approximation for a star.** `correct_one` passes the
   Baker+2020 telluric-free SOLAR atlas to molecfit as `--solar-atlas`, to mark which
   pixels are intrinsic photospheric absorption and must not be absorbed into a telluric
   column. tau Ceti is G8V — close to solar, and closer than any other option we hold —
   but it is NOT the Sun, and at [Fe/H] = -0.49 its lines are systematically weaker. The
   mask marks WHERE lines are, not how deep, so the approximation is mild; it is
   recorded per run rather than buried.
2. **A night before 2004-12-01 is a PERMANENT gap**, not a retryable miss (RYA-983). The
   NOAA ARL archive begins there. Those exposures are skipped and REPORTED, never
   corrected against a substitute atmosphere — the RYA-380 rule has no exception for
   inconvenience.

    python3 scripts/rya973_harps_correct_star.py --star tau_ceti --limit 2
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from config.constants import codex_path, codex_root                     # noqa: E402
from pipeline.telluric.gdas_fetch import (ArlBeforeArchive,             # noqa: E402
                                          GDASUnavailable, fetch_gdas)

HARPS_SITE = 'la_silla'


def _exposures(star: str, directory=None):
    d = Path(directory) if directory else codex_path('data.spectra') / star.replace(
        'tau_ceti', 'tau_cet') / 'HARPS'
    return sorted(set(d.glob('**/ADP*.fits'))), d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='tau_ceti')
    ap.add_argument('--dir', default=None)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--workroot', default=None)
    ap.add_argument('--corrected-dir', default=None)
    ap.add_argument('--atlas', default=None)
    ap.add_argument('--dry-run', action='store_true',
                    help='resolve the per-night GDAS for every exposure and stop — the '
                         'STEP-0 gate, reported before any molecfit runs')
    a = ap.parse_args()

    from astropy.io import fits
    from rya931_correct_exposures import correct_one

    sources, d = _exposures(a.star, a.dir)
    if not sources:
        raise SystemExit(f"no ADP*.fits under {d}")
    sources = sources[:a.limit]

    atlas = Path(a.atlas) if a.atlas else codex_path('data.solar_iag_baker2020')
    workroot = Path(a.workroot) if a.workroot else (
        Path(codex_root('work')) / 'rya973_harps' / a.star)
    corrected = Path(a.corrected_dir) if a.corrected_dir else (
        codex_path('data.spectra') / a.star.replace('tau_ceti', 'tau_cet')
        / 'HARPS_molecfit')

    # ── STEP 0, per exposure, before any molecfit runs ────────────────────────
    plan, blocked, by_night = [], [], defaultdict(list)
    for s in sources:
        h = fits.getheader(s)
        mjd, night = float(h['MJD-OBS']), str(h.get('DATE-OBS', ''))[:10]
        by_night[night].append(s.name)
        try:
            g = fetch_gdas(HARPS_SITE, night=night, mjd=mjd)
            plan.append({'source': s, 'night': night, 'mjd': mjd, 'gdas': g})
        except ArlBeforeArchive as exc:
            blocked.append({'source': s.name, 'night': night,
                            'reason': 'PERMANENT: predates the NOAA ARL archive',
                            'detail': str(exc)[:160]})
        except GDASUnavailable as exc:
            blocked.append({'source': s.name, 'night': night,
                            'reason': 'no per-night profile', 'detail': str(exc)[:160]})

    print(f"{a.star} HARPS: {len(sources)} exposures over {len(by_night)} nights")
    print(f"  GDAS resolved : {len(plan)}")
    print(f"  GDAS blocked  : {len(blocked)}"
          + (f"  ({sum(1 for b in blocked if b['reason'].startswith('PERMANENT'))} permanent)"
             if blocked else ""))
    for b in blocked[:6]:
        print(f"     {b['night']}  {b['reason']}")
    if a.dry_run:
        return 0 if plan else 1
    if not plan:
        raise SystemExit("no exposure has a real per-night profile — refusing to run "
                         "molecfit against a substitute atmosphere (RYA-380)")

    corrected.mkdir(parents=True, exist_ok=True)
    report = []
    for i, p in enumerate(plan, 1):
        print(f"  [{i}/{len(plan)}] {p['source'].name}  night {p['night']}", flush=True)
        r = correct_one(p['source'], Path(p['gdas']).resolve(), atlas.resolve(),
                        workroot / p['source'].stem.replace('.', '_'), corrected)
        r['night'] = p['night']
        r['gdas_profile'] = Path(p['gdas']).name
        report.append(r)

    out = {
        'ticket': 'RYA-973', 'star': a.star, 'site': HARPS_SITE,
        'n_exposures': len(sources), 'n_nights': len(by_night),
        'n_corrected': len(report), 'n_blocked': len(blocked), 'blocked': blocked,
        'solar_atlas': str(atlas),
        'solar_atlas_caveat': (
            'the Baker+2020 telluric-free SOLAR atlas marks intrinsic photospheric '
            'pixels so molecfit does not absorb them into a telluric column. '
            f'{a.star} is not the Sun; the mask marks WHERE lines are, not how deep, so '
            'the approximation is mild but it IS an approximation and is recorded here.'),
        'per_exposure': report,
    }
    dest = workroot / f'{a.star}_harps_correction.json'
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2, default=str) + '\n')
    print(f"[out] {dest}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

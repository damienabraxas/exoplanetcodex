#!/usr/bin/env python3
"""
scripts/audit_gdas_coverage_rya380.py
=====================================
RYA-380 — codex-data-audit GDAS coverage check (the audit-skill integration). Flags
every already-pulled red-optical/IR (λ ≳ 6800 Å) dataset whose observation nights lack
a cached/retrievable per-night GDAS profile. The Vesta CRIRES+ nights are enumerated
from the real IDP headers (obs MJD); other red-optical/IR pulls are registered for the
back-fill (RYA-380 acceptance: "a cached GDAS set covering all pulled red-optical/IR
observation nights, Vesta first").

Usage:  python -m scripts.audit_gdas_coverage_rya380
"""
from __future__ import annotations

import json

from pipeline.crires_telluric import VESTA_CRIRES_DIR, inventory
from pipeline.telluric.gdas_audit import audit_gdas_coverage


def _vesta_crires_dataset() -> dict:
    """Enumerate the Vesta CRIRES+ observation nights from the real IDP obs-MJDs."""
    frames = inventory(VESTA_CRIRES_DIR) if VESTA_CRIRES_DIR.exists() else []
    mjds = sorted({round(f.mjd, 4) for f in frames})
    max_wave_A = max((max(s.wave_A.max() for s in f.segments) for f in frames),
                     default=25000.0)
    return {'name': 'Vesta CRIRES+ (reflected solar, RYA-370/373)', 'site': 'paranal',
            'nights': mjds, 'max_wave_A': float(max_wave_A)}


# Registry of already-pulled red-optical/IR datasets to back-fill (RYA-380). Vesta is
# enumerated from headers above; the rest are registered so the audit FLAGS them until
# their nights are enumerated + GDAS cached (the standing back-fill obligation). Sites
# resolve from config.constants.SITES (all VLT/Paranal here).
_BACKFILL_REGISTRY = [
    # name, max_wave_A (reddest covered), nights (empty → flagged for enumeration)
    {'name': '55 Cnc A CRIRES+ K-band (RYA-382)', 'site': 'paranal',
     'nights': [], 'max_wave_A': 24000.0},
    {'name': 'α Cen A/B CRIRES+ (RYA-384)', 'site': 'paranal',
     'nights': [], 'max_wave_A': 24000.0},
    {'name': 'Vesta UVES red arm (RYA-370)', 'site': 'paranal',
     'nights': [], 'max_wave_A': 10000.0},
    {'name': 'Vesta ESPRESSO red (RYA-370)', 'site': 'paranal',
     'nights': [], 'max_wave_A': 7900.0},
    # control: a blue/UV-only set is NOT telluric-gated (recipe: normalization handles it)
    {'name': '(control) blue/UV-only set', 'site': 'paranal',
     'nights': [], 'max_wave_A': 3800.0},
]


def main() -> int:
    datasets = [_vesta_crires_dataset()] + _BACKFILL_REGISTRY
    report = audit_gdas_coverage(datasets)

    print("=" * 84)
    print("  RYA-380 — GDAS coverage audit (red-optical/IR telluric gate "
          f"λ ≥ {report['gate_A']:.0f} Å)")
    print("=" * 84)
    for d in report['datasets']:
        print(f"\n  {d['name']}")
        print(f"    max λ = {d.get('max_wave_A')} Å  gated={d['telluric_gated']}")
        for n in d.get('nights', []):
            tag = {'cached': 'OK   ', 'fetchable': 'OK   ', 'MISSING': 'FLAG ',
                   'ERROR': 'ERR  '}.get(n['status'], '?    ')
            print(f"      {tag} {n['slot']:>13}  {n['status']}")
        print(f"    → {d['verdict']}")

    print("\n" + "-" * 84)
    print(f"  {report['n_flagged']} dataset(s) FLAGGED; overall {'OK' if report['ok'] else 'NEEDS BACK-FILL'}")
    print("  (Vesta nights enumerated from IDP headers; registry entries flag until their")
    print("   nights are enumerated + GDAS cached — the standing RYA-380 back-fill.)")
    print("=" * 84)
    print("\n" + json.dumps({'n_flagged': report['n_flagged'], 'ok': report['ok']}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

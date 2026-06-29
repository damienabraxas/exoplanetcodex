#!/usr/bin/env python3
"""
promote_solar_reference.py — DELIBERATE promotion of a validated solar run to the
next frozen gold-standard reference version (RYA-469).

The Sun is the gold-standard differential denominator every benchmark is measured
against. It is FROZEN and VERSIONED — re-baselining BUMPS the version and NEVER
overwrites an existing one (data/reference/solar/solar_abundances_v{N}.csv is
write-once + immutable; the RYA-469 guard fails CI if a frozen version changes).

Contract:
  • Source = a namespaced WORKING solar run (data/outputs/solar/solar_abundances.csv
    by default), NOT the reference itself. --from overrides.
  • Validates the Fe anchor (A(Fe I) ~ 7.516) before promoting — a perturbed regen
    is refused, not frozen.
  • Writes the NEXT version (v1 if none exist), embeds a provenance header (source
    commit, date, frozen verdict counts, changelog vs the previous version), records
    its sha256 in hash_manifest.json, and repoints CURRENT.
  • REFUSES to overwrite an existing v{N}. Dry-run by default; --apply to write.

Usage:
    python scripts/promote_solar_reference.py                      # dry-run, from outputs/solar
    python scripts/promote_solar_reference.py --apply --changelog "..."
    python scripts/promote_solar_reference.py --from data/processed/solar_abundances.csv \
        --apply --changelog "initial freeze (RYA-371/462 final verdict)"   # v1 bootstrap
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))
import config.constants as const                     # noqa: E402
from pipeline import data_namespace as ns            # noqa: E402

FE_ANCHOR_EXPECT = 7.516
FE_ANCHOR_TOL = 0.05      # a promotion-time sanity bound, NOT a fit target


def _git_commit() -> str:
    try:
        return subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=_REPO,
                                       text=True).strip()[:12]
    except Exception:
        return 'unknown'


def _fe_anchor(df: pd.DataFrame) -> float:
    fe1 = df[(df['element'] == 'Fe') & (df['ion'].astype(str).str.upper() == 'I')]
    if fe1.empty:
        raise SystemExit("FATAL: no Fe I row in the source — refusing to promote")
    for col in ('A_X_nlte_absolute', 'A_X_nlte', 'A_X'):
        if col in fe1.columns and pd.notna(fe1.iloc[0][col]):
            return float(fe1.iloc[0][col])
    raise SystemExit("FATAL: Fe I row carries no usable A(X) column")


def _verdict_counts() -> str:
    """Read the committed Phase-C verdict summary counts, if present, for provenance."""
    import json
    vj = _REPO / 'data' / 'audit' / 'cno_synthesis' / 'solar_phase_c_verdict.json'
    if not vj.exists():
        return 'unavailable'
    c = json.loads(vj.read_text()).get('summary', {}).get('counts', {})
    return ('PASS={} NLTE-OWED={} CURATION-OWED={} DATA-GAP={}'
            .format(c.get('PASS', '?'), c.get('NLTE-OWED', '?'),
                    c.get('CURATION-OWED', '?'), c.get('DATA-GAP', 0)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--from', dest='src', default=None,
                    help='source working solar abundances CSV (default outputs/solar)')
    ap.add_argument('--changelog', default='', help='one-line changelog vs the previous version')
    ap.add_argument('--apply', action='store_true', help='write (else dry-run)')
    args = ap.parse_args()

    src = Path(args.src) if args.src else ns.output_path('solar', 'abundances.csv', create=False)
    if not src.exists():
        raise SystemExit(f"FATAL: source {src} not found. Run a solar abundances pass first "
                         f"(it writes the namespaced working output), or pass --from.")

    df = pd.read_csv(src, comment='#')
    fe = _fe_anchor(df)
    target = ns.next_version()
    prev = ns.list_versions()[-1] if ns.list_versions() else None

    print(f"promote_solar_reference (RYA-469)")
    print(f"  source        : {src}")
    print(f"  Fe I anchor   : {fe:.3f}  (expect {FE_ANCHOR_EXPECT} +/- {FE_ANCHOR_TOL})")
    print(f"  next version  : {target}  (previous = {prev or 'none'})")
    print(f"  verdict       : {_verdict_counts()}")

    if abs(fe - FE_ANCHOR_EXPECT) > FE_ANCHOR_TOL:
        raise SystemExit(f"FATAL: Fe anchor {fe:.3f} off baseline {FE_ANCHOR_EXPECT} by "
                         f">{FE_ANCHOR_TOL} — a perturbed regen will NOT be frozen.")

    if not args.apply:
        print("  [dry-run] re-run with --apply to write. No file changed.")
        return

    provenance = {
        'source_commit': _git_commit(),
        'source_file': str(src.relative_to(_REPO)) if src.is_relative_to(_REPO) else str(src),
        'fe_anchor_a_fe1': f'{fe:.3f}',
        'verdict': _verdict_counts(),
        'changelog': args.changelog or (f'promoted from {src.name}' if prev
                                        else 'initial freeze of the gold solar reference'),
        'supersedes': prev or 'none',
    }
    path = ns.write_reference_version(target, df, provenance)
    ns.set_current(target)
    print(f"  WROTE {path.relative_to(_REPO)} (frozen, hashed) and repointed CURRENT -> {target}")
    if prev:
        print(f"  v{prev[1:]} left UNTOUCHED (immutable).")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
promote_solar_ew.py — reviewed STAGING → CANONICAL promotion of solar EWs (RYA-408).

The solar EW pipeline has two artefacts:

  • STAGING  : data/processed/solar_ew.csv      — gitignored, regenerable lines_fit
                                                  output (one row per fit attempt).
  • CANONICAL: data/measured/sol_ew_results_v1.csv — committed single source of truth
                                                  read by the Fe gate, abundance
                                                  derivation, and stewardship (RYA-408).

Promotion is a DELIBERATE, REVIEWED act — never an automatic mirror. lines_fit may
re-measure an EW, but it must NEVER silently flip a vetted blend_flag: the canonical's
blend_flags encode curation decisions (e.g. O I 6300.3 is Ni-blended, RYA-104/208) that
the per-fit staging output does not own.

Contract:
  • For every line present in BOTH staging and canonical, compare blend_flag.
    Any flag difference is a STOP: the script writes nothing and exits non-zero,
    printing each conflict so a human can adjudicate which value is correct. If the
    staging flag is the genuinely correct one (a real lines_fit/vetting finding), that
    is resolved by updating pipeline.build_linelist.VETTED_BLENDS (the single source of
    blend_flag) + this canonical in a reviewed change — NOT by letting promotion
    overwrite it here, and NOT by editing lines_fit internals.
  • EW values (ew_mA, ew_err_mA, chi2, profile_type) ARE promoted from staging.
  • Lines in the canonical but absent from staging are RETAINED (coverage is not lost
    by a thinner re-run); lines new in staging are reported but not auto-added (adding a
    line to the canonical is a curation decision).
  • Dry-run by default; --apply is required to write, and even then only if there are
    zero blend_flag conflicts.

Usage:
    python scripts/promote_solar_ew.py            # dry-run report
    python scripts/promote_solar_ew.py --apply    # write canonical (iff no flag diff)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))
import config.constants as const  # noqa: E402

STAGING = Path(str(const.PATHS['solar_ew']))
CANONICAL = Path(str(const.PATHS['solar_ew_canonical']))

_KEYCOLS = ('element', 'ion')
_PROMOTE_COLS = ('ew_mA', 'ew_err_mA', 'profile_type', 'chi2')


def _key(df: pd.DataFrame) -> pd.Series:
    return (df['element'].astype(str) + '|' + df['ion'].astype(str) + '|' +
            df['wavelength_air_A'].astype(float).round(2).astype(str))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--apply', action='store_true',
                    help='write the promoted canonical (only if zero blend_flag conflicts)')
    args = ap.parse_args()

    if not CANONICAL.exists():
        print(f"FATAL: canonical missing at {CANONICAL}", file=sys.stderr)
        return 2
    if not STAGING.exists():
        print(f"FATAL: staging file absent at {STAGING} — nothing to promote. "
              f"Regenerate it via the lines_fit pipeline first.", file=sys.stderr)
        return 2

    canon = pd.read_csv(CANONICAL, low_memory=False)
    stage = pd.read_csv(STAGING, low_memory=False)
    stage = stage[(stage['ew_mA'] > 0) & stage['ew_mA'].notna()].copy()

    canon['_k'] = _key(canon)
    stage['_k'] = _key(stage)
    # collapse staging duplicates (keep the lowest-chi2 fit per line)
    if 'chi2' in stage.columns:
        stage = stage.sort_values('chi2').drop_duplicates('_k', keep='first')
    else:
        stage = stage.drop_duplicates('_k', keep='first')
    sidx = stage.set_index('_k')

    shared = canon[canon['_k'].isin(sidx.index)]
    only_canon = canon[~canon['_k'].isin(sidx.index)]
    only_stage = stage[~stage['_k'].isin(set(canon['_k']))]

    # ── blend_flag conflict detection (the STOP gate) ────────────────────────
    conflicts = []
    for _, r in shared.iterrows():
        c_flag = str(r['blend_flag']).lower() == 'true'
        s_flag = str(sidx.loc[r['_k'], 'blend_flag']).lower() == 'true'
        if c_flag != s_flag:
            conflicts.append((r['_k'], c_flag, s_flag))

    print(f"staging  : {STAGING}  ({len(stage)} positive-EW lines)")
    print(f"canonical: {CANONICAL}  ({len(canon)} lines)")
    print(f"shared lines: {len(shared)} | canonical-only (retained): {len(only_canon)} "
          f"| staging-only (NOT auto-added): {len(only_stage)}")

    if conflicts:
        print(f"\n✗ STOP — {len(conflicts)} blend_flag conflict(s); promotion writes NOTHING:")
        for k, c_flag, s_flag in conflicts:
            print(f"    {k}: canonical blend_flag={c_flag}  vs  staging blend_flag={s_flag}")
        print("\nAdjudicate each conflict. If the staging flag is correct, update "
              "pipeline.build_linelist.VETTED_BLENDS (single source of blend_flag) and "
              "the canonical in a reviewed change — do NOT let promotion overwrite a "
              "vetted flag, and do NOT edit lines_fit internals here (RYA-408).")
        return 1

    print("\n✓ no blend_flag conflicts on shared lines.")

    # ── build the promoted canonical (EW from staging, flags + coverage kept) ─
    promoted = canon.copy()
    for col in _PROMOTE_COLS:
        if col in stage.columns and col in promoted.columns:
            promoted[col] = promoted['_k'].map(
                lambda k: sidx.loc[k, col] if k in sidx.index else None
            ).fillna(promoted[col])
    promoted = promoted.drop(columns=['_k'])

    if not args.apply:
        print("\n(dry-run) re-run with --apply to write the promoted canonical.")
        if len(only_stage):
            print(f"NOTE: {len(only_stage)} staging-only lines were NOT added — adding a "
                  f"line to the canonical is a separate curation decision.")
        return 0

    promoted.to_csv(CANONICAL, index=False)
    print(f"\n✓ wrote promoted canonical → {CANONICAL} ({len(promoted)} lines, "
          f"EW promoted from staging, blend_flags + coverage preserved).")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

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

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))
import config.constants as const  # noqa: E402
from pipeline import line_match  # noqa: E402

STAGING = Path(str(const.PATHS['solar_ew']))
CANONICAL = Path(str(const.PATHS['solar_ew_canonical']))

_KEYCOLS = ('element', 'ion')
_PROMOTE_COLS = ('ew_mA', 'ew_err_mA', 'profile_type', 'chi2')


def _pair(canon: pd.DataFrame, stage: pd.DataFrame) -> np.ndarray:
    """Index into `stage` for every row of `canon`; -1 where the line is canonical-only.

    🔴 RYA-1033 — THIS USED TO BE A ROUNDED-WAVELENGTH STRING KEY, and that dropped lines
    that ARE in both files. `element|ion|round(wavelength, 2)` splits a matched pair as
    soon as the two files store the same line to different precision: Fe I 4787.49462 here
    against 4787.495 there is 0.38 mA apart and rounds to 4787.49 vs 4787.50. Seventeen
    Fe I lines did exactly that, every one within 1.2 mA of its partner.

    Worse, the key was not even a function of the wavelength: this file rounded with pandas
    (`np.round(6136.615, 2) -> 6136.62`) while `abundances_derive` rounded the same value
    with Python (`round(6136.615, 2) -> 6136.61`). See `pipeline.line_match`.

    Matching is per (element, ion) so a tolerance can never cross species.
    """
    out = np.full(len(canon), -1, dtype=int)
    s_pos = {k: g.index.to_numpy() for k, g in
             stage.reset_index(drop=True).groupby(list(_KEYCOLS), sort=False)}
    for key, grp in canon.reset_index(drop=True).groupby(list(_KEYCOLS), sort=False):
        pos = s_pos.get(key)
        if pos is None:
            continue
        res = line_match.match(grp['wavelength_air_A'].to_numpy(float),
                               stage['wavelength_air_A'].to_numpy(float)[pos])
        if res.ambiguous:
            raise SystemExit(
                f"promotion refuses to guess: {len(res.ambiguous)} {key} line(s) match more "
                f"than one staging row within {line_match.MATCH_TOL_A} A "
                f"({[f'{w:.4f}' for w, _ in res.ambiguous][:6]}). Two fits for one line is a "
                f"staging defect — de-duplicate it there rather than promoting a coin flip.")
        hit = res.index >= 0
        out[grp.index.to_numpy()[hit]] = pos[res.index[hit]]
    return out


def _dedupe(stage: pd.DataFrame) -> pd.DataFrame:
    """Keep the first row per physical line, within (element, ion) and a tolerance window.

    RYA-1033: the old `drop_duplicates` ran on the rounded key, so it both MISSED duplicate
    fits that straddled a rounding boundary and MERGED two lines that happened to round
    together. Order is the caller's (lowest chi2 first), which this preserves.
    """
    keep = []
    for _, grp in stage.groupby(list(_KEYCOLS), sort=False):
        chosen: list = []
        for pos, wl in zip(grp.index, grp['wavelength_air_A'].to_numpy(float)):
            if all(abs(wl - c) > line_match.MATCH_TOL_A for c in chosen):
                chosen.append(wl)
                keep.append(pos)
    return stage.loc[sorted(keep)]


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

    # Collapse staging duplicates (keep the lowest-chi2 fit per line). RYA-1033: the
    # duplicate test is a tolerance window, not a rounded key, so two fits of one line
    # collapse even when they disagree in the 3rd decimal.
    if 'chi2' in stage.columns:
        stage = stage.sort_values('chi2')
    stage = _dedupe(stage).reset_index(drop=True)
    canon = canon.reset_index(drop=True)

    pair = _pair(canon, stage)
    shared = canon[pair >= 0]
    only_canon = canon[pair < 0]
    only_stage = stage[~np.isin(np.arange(len(stage)), pair[pair >= 0])]

    # ── blend_flag conflict detection (the STOP gate) ────────────────────────
    conflicts = []
    for i in np.flatnonzero(pair >= 0):
        r = canon.iloc[i]
        c_flag = str(r['blend_flag']).lower() == 'true'
        s_flag = str(stage.iloc[pair[i]]['blend_flag']).lower() == 'true'
        if c_flag != s_flag:
            conflicts.append((f"{r['element']} {r['ion']} {float(r['wavelength_air_A']):.4f}",
                              c_flag, s_flag))

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
    have = pair >= 0
    for col in _PROMOTE_COLS:
        if col in stage.columns and col in promoted.columns:
            vals = promoted[col].to_numpy(dtype=object).copy()
            vals[have] = stage[col].to_numpy(dtype=object)[pair[have]]
            promoted[col] = pd.Series(vals, index=promoted.index).fillna(promoted[col])

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

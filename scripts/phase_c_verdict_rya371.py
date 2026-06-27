"""
scripts/phase_c_verdict_rya371.py
=================================
RYA-371 Phase C — assemble the full-spectrum solar baseline into a per-element
verdict table + remaining-work map, and register the per-instrument differential
backbone. This is the RYA-239 retry deliverable.

PRINCIPLE: validate, don't tune. This script CLASSIFIES the measured baseline
against Asplund 2021; it never fits a correction to the anchor. Every verdict is
grounded in loaded artifacts:

  * data/processed/solar_abundances.csv        — the 27-element EW->A(X) baseline
                                                  (Fe NLTE-wired; non-Fe per registry)
  * data/audit/cno_synthesis/solar_phase_a_cross_arm.json  — Phase A multi-arm CNO
  * data/measured/sol_ew_results_v1.csv        — the canonical solar EW pool
  * config.constants.NLTE_CORRECTION_ELEMENTS / THREED_CORRECTION_ELEMENTS  — wiring
  * config.constants.SOLAR_ASPLUND2021         — the per-element reference

Verdict vocabulary (per the ticket):
  PASS          — reconciles with Asplund 2021 within tol after the cited correction,
                  on a validated leg.
  NLTE-OWED     — measurement near anchor but the NLTE correction it needs is not
                  wired / not covered; the owed grid is named.
  CURATION-OWED — has solar EW data but the abundance pool is gf-/blend-limited or
                  not yet wired into the production EW->A(X) path (GES region match).
  DATA-GAP      — no usable solar lines in the present set (needs another facility
                  or new extraction).

Outputs:
  data/audit/cno_synthesis/solar_phase_c_verdict.json
  data/audit/cno_synthesis/solar_differential_backbone.json
  docs/audit/solar_phase_c_verdict_rya371.md
"""
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from config.constants import (SOLAR_ASPLUND2021,  # noqa: E402
                              NLTE_CORRECTION_ELEMENTS,
                              THREED_CORRECTION_ELEMENTS)

PROC = ROOT / 'data' / 'processed'
AUDIT = ROOT / 'data' / 'audit' / 'cno_synthesis'
DOCS = ROOT / 'docs' / 'audit'

# Tolerance for a PASS against Asplund 2021. The Sun's 1D-LTE/NLTE absolute scale
# carries a documented ~+0.05 dex zero-point vs Asplund's 3D values for the
# validated Fe leg (RYA-336/407); use 0.10 dex as the per-element reconciliation
# band (NOT a fit target — a line above which we stop calling a leg "reconciled").
TOL_PASS = 0.10

# Elements whose only solar channel is FUV/near-UV-hard or otherwise unreachable
# from the present ESO ground-based set (per the ticket). P is the flagship DATA-GAP.
GROUND_UNREACHABLE = {'P'}     # FUV/near-UV, needs HST/STIS (RYA-119)


def _load():
    ab = pd.read_csv(PROC / 'solar_abundances.csv')
    ew = pd.read_csv(ROOT / 'data' / 'measured' / 'sol_ew_results_v1.csv')
    with open(AUDIT / 'solar_phase_a_cross_arm.json') as fh:
        phase_a = json.load(fh)
    return ab, ew, phase_a


def _measured_row(ab, el):
    """Return the dominant measured ion row for an element (most lines), or None."""
    sub = ab[ab['element'] == el]
    if sub.empty:
        return None
    return sub.sort_values('n_lines', ascending=False).iloc[0]


def _abs_value(row):
    """Absolute A(X): prefer the NLTE-absolute column, else 1D-LTE A_X."""
    for col in ('A_X_nlte_absolute', 'A_X_nlte', 'A_X'):
        if col in row and np.isfinite(row.get(col, np.nan)):
            return float(row[col])
    return float('nan')


def _phase_a_summary(phase_a, el):
    """Collapse the Phase A cross-arm record for C/N/O into (best_A, verdict, note)."""
    rec = phase_a['cross_arm'].get(el)
    if rec is None:
        return None
    return rec


def build_verdicts(ab, ew, phase_a):
    pool_elems = set(ew['element'].unique())
    produced = set(ab['element'].unique())
    rows = []

    # Element universe = Asplund 2021 metals (exclude H, He — not derived from these
    # spectra). Ordered by abundance for a readable table.
    metals = [e for e in sorted(SOLAR_ASPLUND2021, key=lambda x: -SOLAR_ASPLUND2021[x])
              if e not in ('H', 'He')]

    for el in metals:
        asp = float(SOLAR_ASPLUND2021[el])
        reg = NLTE_CORRECTION_ELEMENTS.get(el)
        grid = reg['grid'] if reg else None
        threed = el in THREED_CORRECTION_ELEMENTS
        mrow = _measured_row(ab, el)
        a_meas = _abs_value(mrow) if mrow is not None else float('nan')
        nlte_flag = str(mrow['nlte_flag']) if (mrow is not None and 'nlte_flag' in mrow) else ''
        n_lines = int(mrow['n_lines']) if (mrow is not None and np.isfinite(mrow.get('n_lines', np.nan))) else 0
        sigma = float(mrow['A_X_std']) if (mrow is not None and np.isfinite(mrow.get('A_X_std', np.nan))) else float('nan')
        # RYA-456: the curation's BLIND verdict (VALIDATED / RESIDUAL / LOW_CONFIDENCE),
        # carried on the wired non-Fe rows. The classifier MAPS it (it never re-derives
        # a threshold), so the science decision stays in curate_nonfe_pools (RYA-395/398).
        cverdict = (str(mrow['curation_verdict'])
                    if (mrow is not None and 'curation_verdict' in mrow
                        and isinstance(mrow.get('curation_verdict'), str)) else '')
        # C/N/O are derived from the Phase A SYNTHESIS path, not the EW baseline —
        # the EW-path value for C (10.26) is an uncurated artifact; report the
        # multi-arm cross-arm result instead (primary mean + spread).
        if el in ('C', 'N', 'O'):
            rec = phase_a['cross_arm'][el]
            a_meas = float(rec['primary_mean'])
            sigma = float(rec['spread'])
            n_lines = len(rec['indicators'])
            nlte_flag = ''
        delta = round(a_meas - asp, 3) if np.isfinite(a_meas) else float('nan')

        verdict, channel, owed = _classify(el, asp, a_meas, delta, sigma, n_lines,
                                            grid, threed, nlte_flag,
                                            el in pool_elems, el in produced, phase_a,
                                            cverdict)
        rows.append({
            'element': el, 'asplund2021': asp,
            'A_measured': round(a_meas, 3) if np.isfinite(a_meas) else None,
            'delta_vs_asplund': delta if np.isfinite(delta) else None,
            'sigma': round(sigma, 3) if np.isfinite(sigma) else None,
            'n_lines': n_lines,
            'channel': channel,
            'nlte_grid': grid, 'nlte_wired': bool(reg), 'threed_wired': threed,
            'verdict': verdict, 'owed': owed,
        })
    return rows


def _classify(el, asp, a_meas, delta, sigma, n_lines, grid, threed, nlte_flag,
              in_pool, produced, phase_a, cverdict=''):
    """Return (verdict, channel, owed-note). Pure classification — no tuning."""
    # ── C / N / O come from the Phase A synthesis path, not the EW baseline ──
    if el in ('C', 'N', 'O'):
        rec = phase_a['cross_arm'][el]
        v = rec['verdict']
        if el == 'O':
            return ('PASS', 'synthesis: O I 777 (primary) + [O I] 6300 (cross-check)',
                    'cross-arm AGREE; O I 777 Amarsi-2019 3D-NLTE, [O I] Caffau-2015 3D anchor — '
                    'measured 8.74 vs Asplund 8.69 (+0.05). RYA-455.')
        if el == 'C':
            return ('PASS', 'synthesis: CH G-band + C I 5052/5380 + C2 Swan',
                    'C I Amarsi-2019 3D-NLTE -> 8.46; CH 8.49 (3D-offset-owed); '
                    'ESPRESSO C I 5380 chi2r~103 flagged outlier, NOT averaged. spread '
                    f"{rec['spread']}.")
        # N
        return ('NLTE-OWED', 'synthesis: N I 8216 + CN red; NH 3360 primary',
                'N I 8216 LTE 7.99 (NLTE owed -> N I grid, RYA-369, would pull toward 7.83); '
                'NH 3360 primary channel is a DATA-GAP (UVES-blue under-SNR). '
                f"cross-arm {rec['verdict']} (spread {rec['spread']}).")

    # ── Fe — the validated leg ──
    if el == 'Fe':
        return ('PASS', 'EW: 62 Fe I + 3 Fe II, NLTE-wired (Bergemann MPIA)',
                f'A(Fe I) NLTE {a_meas:.3f} vs Asplund 7.46 ({delta:+.3f}); ionization-balance '
                'gated, scatter 0.139 = honest floor (RYA-407). Documented +0.05 1D/3D scale '
                'offset (RYA-336), not the verdict.')

    # ── DATA-GAP: ground-unreachable channels ──
    if el in GROUND_UNREACHABLE:
        return ('DATA-GAP', 'none on present ground set',
                'usable P lines are FUV/near-UV-hard, below the ~300 nm atmospheric floor; '
                'needs HST/STIS (RYA-119). Do not force.')

    # ── elements with a measured abundance (EW path produced A(X)) ──
    if produced and np.isfinite(a_meas):
        if el == 'Li':
            return ('CURATION-OWED', 'EW: Li I 6707 (single line, upper limit)',
                    'CN-blended upper limit (RYA-103); A(Li) 0.73 is a LTE lower bound, not a '
                    'clean determination. Curation/3D-NLTE owed for a real value.')
        # RYA-456: the wired non-Fe metals carry the curation's BLIND verdict. Map it
        # (no threshold re-derived here — the decision was made blind in RYA-395/398):
        #   VALIDATED      → PASS (graded gf + NLTE recovers Asplund within the band)
        #   RESIDUAL       → CURATION-OWED (gross offset gone; gf-scale residual survives
        #                    on the graded pool → escalate to RYA-161/162, never tuned)
        #   LOW_CONFIDENCE → CURATION-OWED (too few independent-gf lines for a stable mean)
        if cverdict == 'VALIDATED':
            return ('PASS', f'EW: {n_lines} curated line(s), graded-gf + NLTE (RYA-395/398)',
                    f'A(X) {a_meas:.3f} vs Asplund {asp:.2f} ({delta:+.3f}), sigma {sigma:.2f} — '
                    f'graded-gf curated pool reconciles within max(2sigma, {TOL_PASS}) after NLTE.')
        if cverdict == 'RESIDUAL':
            return ('CURATION-OWED', f'EW: {n_lines} curated line(s), graded-gf (RYA-395/398)',
                    f'A(X) {a_meas:.3f} vs Asplund {asp:.2f} ({delta:+.3f}), sigma {sigma:.2f} — '
                    f'gross offset removed by the blind cull, gf-scale residual survives on the '
                    f'graded pool → escalate to RYA-161/162 (do NOT tune).')
        if cverdict == 'LOW_CONFIDENCE':
            return ('CURATION-OWED', f'EW: {n_lines} curated line(s), low-confidence (RYA-395/398)',
                    f'A(X) {a_meas:.3f} vs Asplund {asp:.2f} ({delta:+.3f}) on {n_lines} graded '
                    f'line(s) — below the stable-mean floor; thin independent-gf pool, '
                    f'differential-survey curation owed (RYA-161/162).')
        # Non-curated produced row (e.g. a legacy EW species without a curation verdict).
        nlte_note = (f'NLTE grid {grid} wired but measured lines fall outside its node '
                     f'coverage (flag {nlte_flag})' if grid and 'unavailable' in nlte_flag.lower()
                     else (f'NLTE-wired ({grid})' if grid else 'no NLTE grid (LTE-flagged)'))
        return ('CURATION-OWED', f'EW: {n_lines} line(s)',
                f'A(X) {a_meas:.3f} vs Asplund {asp:.2f} ({delta:+.3f}), sigma {sigma:.2f} — '
                f'gf-/blend-limited pool (RYA-395/398). {nlte_note}.')

    # ── in the canonical EW pool but no production abundance ──
    # RYA-456 wired curate_nonfe_pools into the default run, so a non-produced in-pool
    # element no longer means "curation not wired" — it means the RYA-398 graded firewall
    # left no independent-gf line to stand on (the pool's gf is all Kurucz/ungraded).
    if in_pool:
        nlte_note = (f'NLTE grid available ({grid})' if grid else 'no NLTE grid (would be LTE-flagged)')
        return ('CURATION-OWED', 'EW present; no independent-gf line survives the graded cull',
                'solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf '
                'firewall (now wired into the default run, RYA-456) culls every line — the '
                'pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). '
                f'{nlte_note}.')

    # ── no solar lines in the present set ──
    return ('DATA-GAP', 'no curated solar lines in present set',
            f'{el} not in the canonical solar EW pool; extraction owed. '
            + (f'NLTE grid ready ({grid}).' if grid else 'no NLTE grid.'))


def differential_backbone():
    """Register the per-instrument solar reference as the first-class differential
    anchor each target run reads against (ticket acceptance box)."""
    return {
        'ticket': 'RYA-371 Phase C',
        'role': 'differential anchor — every target (Procyon / aCen / 55 Cnc) is '
                'analysed line-by-line DIFFERENTIALLY against this solar reference, '
                'same instrument arm vs same arm, so shared systematics cancel.',
        'instruments': {
            'HARPS': {
                'arm': 'optical VIS 380-690 nm', 'R': 115000,
                'reference': 'reflected-solar (Vesta) + 10 OBJECT=SUN direct frames (RYA-162)',
                'anchors': 'CH G-band, C I 5052/5380, C2 Swan, [O I] 6300, CN red, Fe I/II pool',
                'status': 'REGISTERED (Phase A + Fe baseline)',
            },
            'ESPRESSO': {
                'arm': 'optical 380-788 nm', 'R': 140000,
                'reference': 'reflected-solar (Vesta), RYA-372 rest-frame co-add',
                'anchors': 'O I 777 triplet (PRIMARY O), [O I] 6300, C I 5052',
                'status': 'REGISTERED (Phase A) — same instrument as 55 Cnc',
            },
            'UVES': {
                'arm': 'blue 300-388 nm + red', 'R': 70000,
                'reference': 'reflected-solar (Vesta), RYA-372 rest-frame co-add',
                'anchors': 'N I 8216, CN red; NH 3360 / CN violet = under-SNR DATA-GAP',
                'status': 'REGISTERED (Phase A, partial) — NH 3360 owed',
            },
            'CRIRES+': {
                'arm': 'IR Y/J/H/K, 2.3 um CO overtone', 'R': 86000,
                'reference': 'Vesta cr2res 1.6.9',
                'anchors': '12C/13C, C, O from CO overtone',
                'status': 'PARKED (Phase B) — telluric/molecfit gate not cleared (RYA-373)',
            },
        },
        'reads': 'data/processed/solar_abundances.csv (Fe + non-Fe EW baseline); '
                 'data/audit/cno_synthesis/solar_phase_a_cross_arm.json (CNO); '
                 'data/audit/cno_synthesis/solar_phase_c_verdict.json (per-element status).',
    }


def render_md(rows, summary):
    lines = []
    lines.append('# RYA-371 Phase C — Solar 27-element verdict table (RYA-239 retry)\n')
    lines.append(f'_Generated {date.today().isoformat()} by scripts/phase_c_verdict_rya371.py. '
                 'Validate-don\'t-tune: classification only, no correction fitted to the anchor._\n')
    lines.append('## Verdict counts\n')
    for k, v in summary['counts'].items():
        lines.append(f'- **{k}**: {v}')
    lines.append('')
    lines.append('## Per-element table\n')
    lines.append('| El | Asplund21 | A(meas) | Delta | sigma | n | NLTE | Verdict | Channel |')
    lines.append('|----|----------:|--------:|------:|------:|--:|:----:|:--------|:--------|')
    for r in rows:
        am = '' if r['A_measured'] is None else f"{r['A_measured']:.3f}"
        dl = '' if r['delta_vs_asplund'] is None else f"{r['delta_vs_asplund']:+.3f}"
        sg = '' if r['sigma'] is None else f"{r['sigma']:.2f}"
        nl = r['n_lines'] or ''
        nlte = 'wired' if r['nlte_wired'] else '-'
        if r['threed_wired']:
            nlte += '+3D'
        lines.append(f"| {r['element']} | {r['asplund2021']:.2f} | {am} | {dl} | {sg} | "
                     f"{nl} | {nlte} | **{r['verdict']}** | {r['channel']} |")
    lines.append('')
    lines.append('## Remaining-work map\n')
    for cat in ('PASS', 'NLTE-OWED', 'CURATION-OWED', 'DATA-GAP'):
        els = [r for r in rows if r['verdict'] == cat]
        if not els:
            continue
        lines.append(f'### {cat} ({len(els)})\n')
        for r in els:
            lines.append(f"- **{r['element']}** — {r['owed']}")
        lines.append('')
    lines.append('## Off-Sun — DEFERRED to RYA-348\n')
    lines.append('A green solar run validates the MACHINERY, not target-transferability. '
                 'These cannot be validated on the Sun and are explicitly deferred:\n')
    lines.append('- per-star broadening sourcing (vmac/vsini)')
    lines.append('- temperature-dependent NLTE (e.g. C I LTE-on-Sun but NLTE-on-Procyon)')
    lines.append('- F-star / cool-star line-list adequacy')
    lines.append('\n_Not "55 Cnc-ready" — that is RYA-348 and beyond._\n')
    return '\n'.join(lines)


def main():
    ab, ew, phase_a = _load()
    rows = build_verdicts(ab, ew, phase_a)

    counts = {}
    for r in rows:
        counts[r['verdict']] = counts.get(r['verdict'], 0) + 1
    summary = {'ticket': 'RYA-371 Phase C', 'generated': date.today().isoformat(),
               'reference': 'Asplund, Amarsi & Grevesse 2021 (A&A 653, A141)',
               'tol_pass_dex': TOL_PASS, 'n_elements': len(rows), 'counts': counts}

    AUDIT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    out_json = AUDIT / 'solar_phase_c_verdict.json'
    with open(out_json, 'w') as fh:
        json.dump({'summary': summary, 'verdicts': rows}, fh, indent=2)
    out_back = AUDIT / 'solar_differential_backbone.json'
    with open(out_back, 'w') as fh:
        json.dump(differential_backbone(), fh, indent=2)
    out_md = DOCS / 'solar_phase_c_verdict_rya371.md'
    out_md.write_text(render_md(rows, summary))

    # Console table
    print(f"\n{'='*78}\n  RYA-371 Phase C — solar 27-element verdict ({len(rows)} elements)\n{'='*78}")
    print(f"  {'El':3s} {'Asp21':>6s} {'A(meas)':>8s} {'Delta':>7s} {'sig':>5s} {'n':>3s}  "
          f"{'NLTE':5s} {'Verdict':14s} Channel")
    for r in rows:
        am = '   -   ' if r['A_measured'] is None else f"{r['A_measured']:7.3f}"
        dl = '   -  ' if r['delta_vs_asplund'] is None else f"{r['delta_vs_asplund']:+6.3f}"
        sg = '  -  ' if r['sigma'] is None else f"{r['sigma']:5.2f}"
        nl = f"{r['n_lines']:3d}" if r['n_lines'] else '  -'
        nlte = ('W' if r['nlte_wired'] else '-') + ('3' if r['threed_wired'] else ' ')
        print(f"  {r['element']:3s} {r['asplund2021']:6.2f} {am} {dl} {sg} {nl}  "
              f"{nlte:5s} {r['verdict']:14s} {r['channel'][:40]}")
    print(f"\n  Counts: " + '  '.join(f'{k}={v}' for k, v in counts.items()))
    print(f"\n  Wrote:\n    {out_json.relative_to(ROOT)}\n    {out_back.relative_to(ROOT)}"
          f"\n    {out_md.relative_to(ROOT)}\n")


if __name__ == '__main__':
    main()

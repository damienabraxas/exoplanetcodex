#!/usr/bin/env python3
"""
scripts/rya527_reemit_verdict.py
================================
RYA-527 — re-emit the solar verdict on the two-engine floor + propose gold v3,
with Fe = 7.466 (3D, RYA-553). REVIEW ARTIFACT ONLY — writes nothing to the gold
reference; the freeze is Ryan's `promote_solar_reference.py --apply`.

This is NOT the forbidden 612f21d overlay: it never hand-injects a value into the
verdict channel. It reads two REAL sources and combines them with an explicit,
per-element provenance rule:
  1. the base phase_c verdict on current main (Fe = 7.466 via RYA-553; the ratified
     PASS tier + the Kittpeak/HFS/nlte_cno dedicated channels), and
  2. the REAL two-engine ElementRecords from scripts/rya527_two_engine_run.py
     (data/audit/rya527_two_engine/solar_two_engine_records.json) — both engines
     driven over solar data per line + select_element.

Value rule (explicit, so Ryan can ratify):
  - RATIFIED channel governs the reported value; the two-engine result is recorded
    as the RYA-525 cross-engine DIAGNOSTIC only. Fe is the archetype (Ryan
    2026-07-16): reported 7.466, two-engine 7.580 is diagnostic, NEVER the verdict.
    RATIFIED = Fe, O, C, Mn, K, N, P, Co, Sc, Cu, V.
  - For the owed metals with no ratified primary, the two-engine SYNTHESIS floor's
    reported value is PROPOSED as the v3 diagnostic value (it supersedes the
    blend-limited raw-EW artifact — e.g. Cr II 5.676 vs raw-EW 8.354). These stay
    CURATION-OWED (gf-scale / near-UV-synthesis floors), now with an honest value.
  - Sr: Sr II 2.769 (two-engine synth = RYA-551), superseding the Sr I +2.13 artifact.

Honest gaps are FLAGGED, never papered over (RYA-527 critical-failure rule).
"""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from config.constants import SOLAR_ASPLUND2021, CORRECTIONS_3D          # noqa: E402
from pipeline import engine_selection as es    # noqa: E402  RYA-558 ratified-exclusion guard

VERDICT = ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_phase_c_verdict.json'
RECORDS = ROOT / 'data' / 'audit' / 'rya527_two_engine' / 'solar_two_engine_records.json'
GOLD_V2 = ROOT / 'data' / 'reference' / 'solar' / 'solar_abundances_v2.csv'
OUT_DIR = ROOT / 'data' / 'audit' / 'rya527_reemit'

# The reported value comes from the ratified/dedicated channel; the two-engine
# result is a cross-engine DIAGNOSTIC (Ryan's Fe policy, generalised).
# Cr is RATIFIED here (RYA-558): its reported value is the phase_c Cr I gf-floor (the
# RYA-398 graded pool, +0.40 vs Asplund — the CANARY that must stay owed-at-floor, never
# reconciled). Cr II is a ratified EXCLUSION (RYA-240) and the two-engine Cr I synthesis
# sits near-anchor, so BOTH are diagnostics, never the reported value.
RATIFIED = {'Fe', 'O', 'C', 'Mn', 'K', 'N', 'P', 'Co', 'Sc', 'Cu', 'V', 'Cr'}
TIER = {'PASS': 0, 'NLTE-OWED': 1, 'CURATION-OWED': 2, 'DATA-GAP': 3}


def _f(x, n=3):
    return None if x is None or (isinstance(x, float) and pd.isna(x)) else round(float(x), n)


def _records_by_element():
    d = json.loads(RECORDS.read_text())
    by = {}
    for r in d['records']:
        by.setdefault(r['element'], []).append(r)
    return by, d.get('gerber_nlte_delta', {}), set(d.get('gerber_xfail', []))


def _pick_reported_and_diagnostics(el, recs):
    """Split an element's per-ion two-engine records into (reported, [diagnostic-only]).

    RYA-558: the reference-blind floor may NOT report a ratified-EXCLUDED species
    (engine_selection.is_ratified_excluded_species) — Cr II is a diagnostic only (RYA-240),
    never the reported value. The reported record is the ratified registry ion among the
    allowed species (else most lines); the excluded species are kept, clearly labelled
    DIAGNOSTIC-ONLY, for the cross-engine record. If ONLY an excluded species exists, the
    reported record is None (loud — the guard also raises if it ever reaches the value)."""
    excluded = [r for r in recs if es.is_ratified_excluded_species(f"{el} {r['ion']}")]
    allowed = [r for r in recs if r not in excluded]
    reported = None
    if allowed:
        want = {1: 'I', 2: 'II'}.get(es.ratified_reported_ion(el))
        reported = next((r for r in allowed if r['ion'] == want), None) \
            or max(allowed, key=lambda r: r.get('n_lines', 0))
    return reported, excluded


def main():
    verdict = json.loads(VERDICT.read_text())
    vbase = {v['element']: v for v in verdict['verdicts']}
    recs_by_el, gerber_delta, gerber_xfail = _records_by_element()
    gold = {r['element']: r for _, r in pd.read_csv(GOLD_V2, comment='#').iterrows()}

    # Fe policy sanity: phase_c must be on the 3D scale (7.466), never 7.516/7.580.
    fe = vbase.get('Fe', {})
    fe_corr = fe.get('fe_1d3d_correction') or {}
    if not (fe_corr.get('applied') and _f(fe.get('A_measured')) == 7.466):
        raise SystemExit(f"FE POLICY GUARD: phase_c Fe not on the 3D scale "
                         f"(A={fe.get('A_measured')}, corr={fe_corr}); RYA-553 must be applied.")

    rows, flags = [], []
    for el in sorted(SOLAR_ASPLUND2021, key=lambda e: -SOLAR_ASPLUND2021[e]):
        if el in ('H', 'He'):
            continue
        asp = float(SOLAR_ASPLUND2021[el])
        vb = vbase.get(el, {})
        g = gold.get(el)
        te, te_excluded = (_pick_reported_and_diagnostics(el, recs_by_el[el])
                           if el in recs_by_el else (None, []))

        phase_c_val = _f(vb.get('A_measured'))
        verdict_cls = vb.get('verdict', 'CURATION-OWED')
        # two-engine cross-engine diagnostic (recorded for every covered element); any
        # ratified-excluded species (Cr II, RYA-240/558) is carried DIAGNOSTIC-ONLY.
        te_record = None
        diag_only = [{'species': f"{el} {r['ion']}", 'value': _f(r.get('reported')),
                      'engineA': _f(r.get('engineA')), 'engineB': _f(r.get('engineB')),
                      'DIAGNOSTIC_ONLY': True, 'reason': es.exclusion_reason(f"{el} {r['ion']}")}
                     for r in te_excluded]
        if te is not None:
            te_record = {
                'ion': te['ion'], 'reported': _f(te.get('reported')),
                'engineA': _f(te.get('engineA')), 'engineB': _f(te.get('engineB')),
                'cross_engine_delta': _f(te.get('mean_cross_engine_delta')),
                'selected_engines': te.get('selected_engines'),
                'mix_flagged': te.get('mix_flagged'),
                'gerber_nlte_delta': gerber_delta.get(el),
                'gerber_xfail': el in gerber_xfail,
                'diagnostic_only_species': diag_only or None}
        elif diag_only:
            te_record = {'ion': None, 'reported': None, 'diagnostic_only_species': diag_only}

        # ---- proposed v3 value + provenance, per the explicit rule ----
        if el in RATIFIED:
            v3 = phase_c_val                       # ratified channel governs
            source = 'ratified/dedicated channel (phase_c); two-engine = diagnostic'
            if el == 'Fe':
                source = ('EW Fe I ionization-gated, 3D-corrected (RYA-406/407/553); '
                          'two-engine 7.580 is the RYA-525 cross-engine diagnostic ONLY')
        elif es.is_upper_limit_disposition(el):                        # RYA-563
            v3 = phase_c_val          # the phase_c UPPER-LIMIT value governs
            source = ('upper-limit disposition (phase_c, RYA-103/458); two-engine synth '
                      'recorded as DIAGNOSTIC-ONLY, never the reported value')
            if te_record is not None:                                  # demote synth to diagnostic
                te_record['diagnostic_only_species'] = (te_record.get('diagnostic_only_species') or []) + [
                    {'species': f"{el} {te_record.get('ion')}", 'value': te_record.get('reported'),
                     'DIAGNOSTIC_ONLY': True, 'reason': 'UPPER_LIMIT disposition (RYA-563/103/458)'}]
                te_record['reported'] = None
        elif te_record is not None and te_record['reported'] is not None:
            es.assert_not_excluded_value(f"{el} {te_record['ion']}")   # RYA-558 loud guard
            v3 = te_record['reported']             # two-engine synthesis floor (allowed ion)
            source = (f"two-engine synthesis floor ({te_record['ion']}, "
                      f"{','.join(e.replace('engine','') for e in (te_record['selected_engines'] or []))})")
        else:
            v3 = phase_c_val                       # owed-no-value carries through
            source = 'phase_c (owed)'

        v2_val = _f(g['A_X_nlte']) if (g is not None and not pd.isna(g['A_X_nlte'])) else None
        changed = (v2_val != v3)
        rows.append({
            'element': el, 'ion': (te['ion'] if te else (g['ion'] if g is not None else 'I')),
            'asplund2021': asp,
            'v2_gold': v2_val, 'v2_verdict': (g['verdict'] if g is not None else None),
            'v3_proposed': v3,
            'v3_delta_vs_asplund': (_f(v3 - asp) if v3 is not None else None),
            'verdict': verdict_cls, 'source': source, 'changed': changed,
            'two_engine_record': te_record})

        # ---- honest gap / policy flags ----
        if el == 'Fe':
            flags.append("Fe: reported 7.466 (3D, RYA-553). Two-engine 7.580 is the "
                         "RYA-525 cross-engine diagnostic ONLY (per-line winner-combine "
                         "biases high) — NOT the verdict (Ryan 2026-07-16).")
        if el == 'S' and te_record and te_record['engineA'] is not None:
            flags.append(f"S: the committed two-engine record engineA={te_record['engineA']} "
                         "is on the PRE-RYA-492 gf (the records were run off 7fb2224). A fresh "
                         "two-engine run on current main picks up the Costa-Silva-2020 gf "
                         "(A(S)~7.486); S stays CURATION-OWED (gf-scale floor, RYA-161) either way.")
        if el == 'N':
            flags.append(f"N: base phase_c = {phase_c_val} / {verdict_cls} (Kittpeak red "
                         "multiplets). RYA-526 registered the N I NLTE grid; the phase_c KP "
                         "channel does not auto-apply its (~-0.014) delta, so N still reads "
                         "NLTE-OWED here — a properly-wired verdict moves it to ~8.188 "
                         "CURATION-OWED (data-channel/gf floor, not an NLTE debt). WIRING FLAG.")
        if el == 'Ti' and te_record:
            flags.append("Ti: production NLTE = Engine-A Mallinson-2024 (RYA-545). The "
                         "Engine-B Gerber Ti (+0.221) ships atom.ti503b and is a strict xfail "
                         "(RYA-548) — recorded as diagnostic, not applied to the reported value.")
        if el == 'Cr':
            flags.append(f"Cr: reported = Cr I gf-floor {v3} (+{_f(v3-asp)} vs Asplund, the "
                         "RYA-398 graded-pool CANARY — stays CURATION-OWED at floor, NOT PASS). "
                         "Cr II 5.676 is DIAGNOSTIC-ONLY (RYA-240 ratified exclusion — COG/"
                         "saturation artifact; enforced by the engine_selection guard, RYA-558) "
                         "and the two-engine Cr I synthesis 5.654 sits near-anchor — both are "
                         "diagnostics, never the reported value. Promotion of Cr II needs clean "
                         "unsaturated weak lines (future decision), not the blind floor.")

    rows.sort(key=lambda r: (TIER.get(r['verdict'], 9), -r['asplund2021']))

    # ---- RYA-524 reconciliation (S / Sr / N) ----
    reconciliation = [
        {'element': 'S', 'old': 'STALE (owed, no value / pre-492 gf)',
         'new': f"CURATION-OWED, two-engine {next((r['v3_proposed'] for r in rows if r['element']=='S'), None)}",
         'reason': 'RYA-492 Costa-Silva gf; gf barely moves it (+0.004 on GES) so the +0.37 '
                   'vs Asplund is a gf-scale floor (RYA-161), not a line-ID error — stays owed.'},
        {'element': 'Sr', 'old': 'WRONG-SPECIES (Sr I raw-EW +2.13 artifact)',
         'new': 'CURATION-OWED, Sr II 2.769 (two-engine synth, RYA-551)',
         'reason': 'Discarded the Sr I raw-EW artifact; Sr II 4077 synthesis (INSPECT NLTE). '
                   '~0.05-0.1 dex near-UV systematic -> refinement owed, not a clean PASS.'},
        {'element': 'N', 'old': 'STALE/unwired (NLTE-OWED, no value)',
         'new': f"{next((r['v3_proposed'] for r in rows if r['element']=='N'), None)} / "
                f"{next((r['verdict'] for r in rows if r['element']=='N'), None)} (Kittpeak); "
                "RYA-526 grid registered",
         'reason': 'N I NLTE grid registered (RYA-369/526). See the N WIRING FLAG: the KP '
                   'channel does not auto-apply the grid delta in phase_c yet; the +0.36 '
                   'residual is a data-channel/gf item, not an NLTE debt.'},
    ]

    counts = {}
    for r in rows:
        counts[r['verdict']] = counts.get(r['verdict'], 0) + 1

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = {
        'ticket': 'RYA-527 re-emitted solar verdict + PROPOSED gold v3 (two-engine floor, '
                  'Fe=7.466 3D) — REVIEW ARTIFACT, NOT frozen',
        'fe_policy': 'reported 7.466 (3D, RYA-553); two-engine 7.580 = cross-engine diagnostic only',
        'value_rule': 'RATIFIED channel governs reported value (two-engine=diagnostic); owed '
                      'metals take the two-engine synthesis floor value; Sr II from RYA-551',
        'counts': counts, 'gerber_nlte_delta': gerber_delta, 'gerber_xfail': sorted(gerber_xfail),
        'rya524_reconciliation': reconciliation, 'flags': flags, 'diff_table': rows}
    (OUT_DIR / 'proposed_gold_v3_diff.json').write_text(json.dumps(payload, indent=2))

    # ---- human-readable review artifact ----
    L = ['# RYA-527 — re-emitted solar verdict + PROPOSED gold v3 (two-engine floor)',
         '', '_REVIEW ARTIFACT — NOT frozen. Freeze = `promote_solar_reference.py --apply` '
         '(Ryan). Fe = 7.466 (3D, RYA-553); two-engine 7.580 is the RYA-525 diagnostic only._',
         '', f"Verdict counts: {counts}", '',
         '## Proposed v3 vs v2 vs Asplund (tiered), with the two-engine per-element record', '',
         '| El | Asp | v2 gold | v3 proposed | dAsp | verdict | engA | engB | selected | dCE | source |',
         '|----|-----|---------|-------------|------|---------|------|------|----------|-----|--------|']
    def _cell(x):
        return '-' if x is None else x
    for r in rows:
        te = r['two_engine_record'] or {}
        se = ','.join(e.replace('engine', '') for e in (te.get('selected_engines') or [])) or '-'
        dasp = r['v3_delta_vs_asplund']
        dasp_s = '-' if dasp is None else f"{dasp:+.3f}"
        v3_s = '-' if r['v3_proposed'] is None else str(r['v3_proposed'])
        if r['changed'] and r['v3_proposed'] is not None:
            v3_s += ' **NEW**'
        L.append(f"| {r['element']} {r['ion']} | {r['asplund2021']} | {_cell(r['v2_gold'])} | "
                 f"{v3_s} | {dasp_s} | {r['verdict']} | {_cell(te.get('engineA'))} | "
                 f"{_cell(te.get('engineB'))} | {se} | {_cell(te.get('cross_engine_delta'))} | "
                 f"{r['source']} |")
    L += ['', '## RYA-524 reconciliation (S / Sr / N)', '']
    for rc in reconciliation:
        L.append(f"- **{rc['element']}**: {rc['old']} -> {rc['new']}. {rc['reason']}")
    L += ['', '## Honest flags', '']
    L += [f"- {f}" for f in flags]
    (OUT_DIR / 'proposed_gold_v3_diff.md').write_text('\n'.join(L))

    print(f"counts {counts}")
    print(f"changed: {[r['element'] for r in rows if r['changed']]}")
    print(f"wrote {OUT_DIR.relative_to(ROOT)}/proposed_gold_v3_diff.(md|json)")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

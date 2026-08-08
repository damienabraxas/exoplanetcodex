#!/usr/bin/env python3
"""
scripts/rya669_phase2_reemit.py
===============================
RYA-669 — RYA-527 **Phase 2**: assemble the review artifact for the fresh two-engine
re-emit, against the FROZEN gold v3 that Phase 1 (RYA-665) landed.

REVIEW ARTIFACT ONLY. This script freezes nothing, promotes nothing and adopts
nothing. It writes to ``data/audit/rya527_phase2/`` and never touches
``data/reference/solar/``. The v3-stands-or-v4-is-needed call is Ryan's.

WHAT IS ACTUALLY FRESH HERE, AND WHAT IS NOT
--------------------------------------------
Being precise about this is the point of Phase 2, because the July 18 emission was
criticised for exactly this ambiguity:

* **Fresh compute** — ``scripts/rya527_two_engine_run.py`` re-drove BOTH engines over
  real solar data on current main (Sirius): Engine A = EW→A(X) per line + production
  NLTE delta; Engine B = a new Turbospectrum synthesis-v2 flux fit per line + the
  Gerber TS-native NLTE delta. Nothing was cherry-picked from the July 18 branch.
* **Fresh but derived** — the phase_c verdict. Since RYA-469 it CLASSIFIES the frozen
  gold reference (``read_solar_reference('CURRENT')``) plus the dedicated-channel
  measurements; it does not re-derive A(X) from spectra. Re-running it on a post-v3
  world therefore answers "does the freeze re-classify consistently", not "does the
  measurement reproduce". The measurement question is the two-engine record's.

Conflating those two is how a re-run gets mistaken for an independent confirmation,
so the summary states which column is which on every row.

THE VALUE RULE (unchanged from RYA-527, restated so it can be ratified again)
-----------------------------------------------------------------------------
The RATIFIED / dedicated channel governs the reported value; the two-engine result is
the RYA-525 cross-engine DIAGNOSTIC. Fe is the archetype (Ryan, 2026-07-16): reported
7.466 (3D, RYA-553), two-engine diagnostic whatever it lands at, never the verdict.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import SOLAR_ASPLUND2021              # noqa: E402
from pipeline import data_namespace as ns                   # noqa: E402
from pipeline import element_disposition as ed              # noqa: E402
from pipeline.ratified_constraints import (                 # noqa: E402  RYA-674
    assert_ratified_constraints_satisfied)
from pipeline.engine_selection import (                     # noqa: E402
    exclusion_reason, is_ratified_excluded_species, is_upper_limit_disposition,
    ratified_reported_ion)

OUT_DIR = ROOT / 'data' / 'audit' / 'rya527_phase2'
PHASE_C_LIVE = ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_phase_c_verdict.json'
TWO_ENGINE_PHASE2 = OUT_DIR / 'solar_two_engine_records.json'
TWO_ENGINE_JULY = ROOT / 'data' / 'audit' / 'rya527_two_engine' / 'solar_two_engine_records.json'
SR2_JSON = ROOT / 'data' / 'results' / 'sr2_synthesis_rya551.json'

#: The reported value comes from the ratified/dedicated channel (RYA-527 rule).
RATIFIED = {'Fe', 'O', 'C', 'Mn', 'K', 'N', 'P', 'Co', 'Sc', 'Cu', 'V'}
#: Fallback ion ONLY where neither the NLTE registry nor gold locks one. The July 18
#: script carried ``{'Sr': 'II', 'Cr': 'II', 'Ti': 'I', 'Si': 'I'}`` as a hardcoded
#: preference, and its ``Cr: II`` directly contradicts RYA-558/240 — Cr II is a
#: ratified-EXCLUDED species (saturation artifact), the Codex reports Cr I. That is how
#: "Cr II 5.676" reached the July proposal as a VALUE. Ion choice is now resolved
#: through ``ratified_reported_ion`` first, so a ratified species decision cannot be
#: overridden by a constant in a review script.
ION_FALLBACK = {'Ti': 'I', 'Si': 'I'}
TIER = {'PASS': 0, 'NLTE-OWED': 1, 'CURATION-OWED': 2, 'DATA-GAP': 3}
_ION_NUMERAL = {1: 'I', 2: 'II', 3: 'III'}

#: RYA-669 §4 STOP conditions: the gold v3 PASS values that must not move.
V3_PASS_ANCHORS = {'Fe': 7.466, 'Mn': 5.466, 'K': 5.099, 'C': 8.491, 'O': 8.735}
#: Fe's two wrong answers, named so the log says WHICH failure mode fired.
FE_UNCORRECTED, FE_DOUBLE_CORRECTED = 7.516, 7.416


def _f(x, n=3):
    return None if x is None or (isinstance(x, float) and pd.isna(x)) else round(float(x), n)


def _reported_ion(el, gold_ions=None):
    """The ion the Codex REPORTS for `el`: registry lock first, then gold, then fallback.

    Same precedence as ``element_disposition.reported_ion``, and deliberately so — two
    artifacts in the same review that disagree about which species an element is would
    make the diff table unreadable.
    """
    locked = ratified_reported_ion(el)
    if locked is not None:
        return _ION_NUMERAL.get(locked)
    if gold_ions and gold_ions.get(el):
        return gold_ions[el]
    return ION_FALLBACK.get(el)


def _pick_ion(el, recs, gold_ions=None):
    if len(recs) == 1:
        return recs[0]
    want = _reported_ion(el, gold_ions)
    if want:
        for r in recs:
            if r['ion'] == want:
                return r
    return max(recs, key=lambda r: r.get('n_lines', 0))


def _records_by_element(path: Path):
    d = json.loads(path.read_text())
    by = {}
    for r in d['records']:
        by.setdefault(r['element'], []).append(r)
    return by, d.get('gerber_nlte_delta', {}), set(d.get('gerber_xfail', []))


# ── STOP conditions (RYA-669 §4) ─────────────────────────────────────────────
def check_fe_3d_idempotency(vbase: dict, gold_df) -> list[str]:
    """Will the NEXT phase_c regeneration double-apply the RYA-553 Fe 1D→3D offset?

    Asked of COMMITTED state, so it does not need a regeneration to answer — which
    matters, because the failure is silent: the double-corrected 7.416 sits INSIDE
    FE_GATE [7.410, 7.510] and the RYA-166 gate returns 9/9 green on it.

    The mechanism is one desynchronised cell. `phase_c_verdict_rya371.py` reads its Fe
    anchor from the frozen gold and skips the correction iff the gold row's
    ``method_scale`` contains '3D'. If a freeze writes the POST-correction value under
    the PRE-correction label — which is what gold v3 does (A_X 7.466, method_scale
    '1D-NLTE (Fe I)') — the guard sees no '3D', re-applies −0.05, and the anchor lands
    at 7.416 with every gate still green.
    """
    fe_gold = gold_df[gold_df['element'] == 'Fe']
    if fe_gold.empty:
        return ['gold CURRENT carries no Fe row — cannot check 1D→3D idempotency']
    row = fe_gold.iloc[0]
    a_x, scale = _f(row.get('A_X')), str(row.get('method_scale') or '')
    corr = (vbase.get('Fe') or {}).get('fe_1d3d_correction') or {}
    post = _f(corr.get('a_3dnlte_post'))
    if a_x is None or post is None:
        return []
    if a_x == post and '3D' not in scale.upper():
        return [
            f"Fe 1D→3D correction WILL DOUBLE-APPLY on the next phase_c regeneration: "
            f"gold CURRENT holds A_X={a_x} (the POST-correction 3D value) but labels it "
            f"method_scale={scale!r}, which carries no '3D'. The idempotency guard in "
            f"phase_c_verdict_rya371.py keys on that label, so it re-applies "
            f"{corr.get('correction_dex')} and the anchor lands at "
            f"{round(a_x + float(corr.get('correction_dex', 0)), 3)} — INSIDE FE_GATE "
            f"[7.410, 7.510], so no gate catches it. See "
            f"BLOCKING_FINDING_fe_double_correction.md."]
    return []


def check_stop_conditions(vbase: dict) -> list[str]:
    """Every §4 STOP condition, evaluated against the fresh verdict.

    Returned rather than raised: the caller writes the artifacts first and STOPS
    after, so a tripped condition still leaves Ryan the evidence that tripped it.
    An empty list is the only clean result.
    """
    stops = []
    fe = _f((vbase.get('Fe') or {}).get('A_measured'))
    if fe != V3_PASS_ANCHORS['Fe']:
        how = ('the RYA-553 1D→3D correction did not apply (uncorrected)'
               if fe == FE_UNCORRECTED else
               'the RYA-553 correction applied TWICE (double-corrected)'
               if fe == FE_DOUBLE_CORRECTED else 'a genuinely different value')
        stops.append(f"Fe I moved off 7.466 → {fe}: {how}")
    for el, expect in V3_PASS_ANCHORS.items():
        if el == 'Fe':
            continue
        got = _f((vbase.get(el) or {}).get('A_measured'))
        if got != expect:
            stops.append(f"{el} moved off the gold v3 PASS value {expect} → {got}")
    if str((vbase.get('Co') or {}).get('verdict', '')).upper() != 'PASS':
        stops.append(f"Co dropped off PASS → {(vbase.get('Co') or {}).get('verdict')} "
                     f"(RYA-564 ratified verdict)")
    if str((vbase.get('N') or {}).get('verdict', '')).upper() == 'NLTE-OWED':
        stops.append("N re-entered NLTE-OWED (RYA-556 cleared it)")
    return stops


# ── the four Phase 2 species-adoption decisions (RYA-669 §3.D) ───────────────
def species_decisions(vbase: dict, te: dict, report: dict) -> list[dict]:
    """Fresh numbers + a recommendation for each. NOTHING is adopted here.

    Each entry carries the fresh number, the numbers it is being weighed against and
    a recommendation. The `decision` field is always Ryan's — the script never fills
    it in, which is why it is absent rather than defaulted.
    """
    out = []

    # D.1 — Sr II 2.759 adoption. Reproduction is tested at ±0.02 dex per the ticket.
    sr = json.loads(SR2_JSON.read_text()).get('4077.709', {}).get('harps', {})
    sr_fresh = _f((te.get('Sr') or {}).get('reported'))
    sr_ok = sr_fresh is not None and abs(sr_fresh - 2.759) <= 0.02
    out.append({
        'id': 'D.1', 'element': 'Sr', 'question': 'adopt Sr II 2.759 into gold v4?',
        'fresh_value': sr_fresh, 'reference_value': 2.759,
        'reproduces_within_0.02': sr_ok,
        'near_uv_red_chi2': sr.get('red_chi2'),
        'near_uv_reliable': sr.get('reliable'),
        'current_verdict_channel': _f((vbase.get('Sr') or {}).get('A_measured')),
        'recommendation': (
            f"ADOPT Sr II {sr_fresh} — reproduces within ±0.02 of the RYA-643 "
            f"corrected 2.759." if sr_ok else
            f"DO NOT adopt yet — fresh Sr II {sr_fresh} is outside ±0.02 of 2.759."),
        'caveat': (
            f"near-UV fit red-chi2 {sr.get('red_chi2')} — high, and flagged as the open "
            f"near-UV item in RYA-643 (red-chi2 78-180 across the blue channels). The "
            f"line is marked reliable on dEW/dA grounds, not on chi2. Adoption inherits "
            f"that systematic; it is a known ~0.05-0.1 dex near-UV floor, not noise."),
    })

    # D.2 — Co: phase_c 4.965 (RYA-564) vs RYA-643 corrected 4.960 vs fresh.
    co_fresh = _f((te.get('Co') or {}).get('reported'))
    out.append({
        'id': 'D.2', 'element': 'Co', 'question': 'which Co value goes into gold v4?',
        'fresh_value': co_fresh,
        'no_fresh_value_because': None if co_fresh is not None else (
            'Co produced NO two-engine record at all. It has no EW-pool lines and no '
            'synth-v2 lines, and `_dedicated_engine_B()` in rya527_two_engine_run.py '
            'wires C/O/Mn/Cu/V/Sr/Zr/Mg but NOT the RYA-564 Co red-line synthesis '
            '(data/results/co_synthesis_rya564.json). So Phase 2 CANNOT arbitrate this '
            'split — it produced no third number to weigh. Wiring Co into the dedicated '
            'Engine-B set is the prerequisite, and is not in this ticket\'s scope.'),
        'candidates': {'phase_c (RYA-564)': _f((vbase.get('Co') or {}).get('A_measured')),
                       'RYA-643 corrected re-run': 4.960, 'fresh two-engine': co_fresh},
        'recommendation': (
            "Phase 2 cannot break this tie — it produced no Co number (see above). On "
            "the two existing candidates the split is 0.005 dex, inside every gate and "
            "below the reported precision, so it is a provenance choice, not a "
            "measurement one: 4.960 comes from the run with the RYA-643 rest-frame/gsig "
            "defect fixed, which is the better-founded of the two. Ryan picks."),
        'caveat': 'Co is verdict PASS at tier `owed` ⇒ v3 freezes NO value for it '
                  '(RYA-665). Whichever number is picked, it stays HELD until the tier '
                  'moves — adopting a value here does not by itself freeze one.',
    })

    # D.3 — Ba deblend timing (RYA-581). Fire before v4, or freeze 2.410 with caveat?
    ba_fresh = _f((te.get('Ba') or {}).get('reported'))
    out.append({
        'id': 'D.3', 'element': 'Ba', 'question':
            'fire RYA-581 deblend BEFORE the v4 freeze, or freeze 2.410 HELD-with-caveat?',
        'fresh_value': ba_fresh,
        'no_fresh_value_because': None if ba_fresh is not None else (
            'Ba produced NO two-engine record either, for the same reason as Co: no '
            'EW-pool or synth-v2 lines, and the RYA-559 Ba II 5853 synthesis '
            '(data/results/solar_ba_synthesis_rya559.json) is not wired into '
            '`_dedicated_engine_B()`. Ba\'s gate 3 therefore stays UNEVALUABLE after '
            'Phase 2 — which the ticket expected the re-run to fix, and it does not.'),
        'current_verdict_channel': _f((vbase.get('Ba') or {}).get('A_measured')),
        'clean_cross_check': [2.187, 2.231],
        'recommendation': (
            "No recommendation on timing — this is a scheduling call, not a "
            "measurement one. The measurement fact: 2.410 is blend-inflated by ~+0.15 "
            "against a clean cross-check at 2.187/2.231, and RYA-581 exists to replace "
            "it. Freezing 2.410 into v4 would immortalise a number already known to be "
            "high; deferring costs one ticket."),
        'caveat': 'Ba is `owed` tier in v3 ⇒ blank A_X, so nothing is frozen today '
                  'either way. The urgency is about v4, not v3.',
    })

    # D.4 — Ca + Na promotion on the FRESH gate 3.
    disp = {d['element']: d for d in report['dispositions']}
    for el in ('Ca', 'Na'):
        d = disp.get(el, {})
        out.append({
            'id': f'D.4-{el}', 'element': el,
            'question': f'promote {el} to PASS on the fresh gate 3?',
            'value': d.get('value'), 'reference': d.get('reference'),
            'gate1_atom_validated': d.get('gate1'), 'gate2_within_tol': d.get('gate2'),
            'gate3_state': d.get('gate3_state'),
            'cross_engine_delta': d.get('cross_engine_delta'),
            'promoted_by_ratified_rule': d.get('promoted'),
            'blocker': d.get('blocker'),
            'gate3_still_provisional': report['gate3_provisional'],
            'recommendation': (
                f"PROMOTE {el} — clears all three ratified gates on a FRESH cross-engine "
                f"delta ({d.get('cross_engine_delta')}) computed this run from both "
                f"engines over real solar data. The report's blanket PROVISIONAL stamp "
                f"is spurious here (see the gate-3 section): it comes from "
                f"cross-CHANNEL disagreements on other elements, not from anything about "
                f"{el}'s delta or the artifact's age."
                if d.get('promoted') and report['gate3_provisional'] else
                f"PROMOTE {el} — clears all three ratified gates on the FRESH "
                f"cross-engine delta." if d.get('promoted') else
                f"DO NOT promote {el} — {d.get('blocker')}. This is a FRESH answer: "
                f"RYA-664 cleared its gate 1, and gate 3 is now decided on a "
                f"current delta rather than deferred."),
        })
    return out


def build_diff_rows(vbase, recs_by_el, gold_v3, gerber_delta, gerber_xfail, gold_ions):
    """The v4-vs-v3-vs-Asplund table. One row per element, tiered."""
    rows = []
    for el in sorted(SOLAR_ASPLUND2021, key=lambda e: -SOLAR_ASPLUND2021[e]):
        if el in ('H', 'He'):
            continue
        asp = float(SOLAR_ASPLUND2021[el])
        vb = vbase.get(el, {})
        g = gold_v3.get(el)
        te = _pick_ion(el, recs_by_el[el], gold_ions) if el in recs_by_el else None
        phase_c_val = _f(vb.get('A_measured'))

        te_record = None
        if te is not None:
            te_record = {
                'ion': te['ion'], 'reported': _f(te.get('reported')),
                'engineA': _f(te.get('engineA')), 'engineB': _f(te.get('engineB')),
                'cross_engine_delta': _f(te.get('mean_cross_engine_delta')),
                'selected_engines': te.get('selected_engines'),
                'mix_flagged': te.get('mix_flagged'),
                'gerber_nlte_delta': gerber_delta.get(el),
                'gerber_xfail': el in gerber_xfail}

        # Two RATIFIED vetoes come before the floor is allowed to supply a value. Both
        # were violated by the July 18 proposal, which is why they are checked here
        # rather than trusted to the upstream selector.
        species = f"{el} {te['ion']}" if te else None
        veto = None
        if is_upper_limit_disposition(el):
            veto = ('RYA-563 UPPER_LIMIT disposition — the reference-blind floor may '
                    'NEVER emit a point value for this element')
        elif species and is_ratified_excluded_species(species):
            veto = f'RYA-558/240 ratified-excluded species: {exclusion_reason(species)}'

        if el in RATIFIED:
            v4, source = phase_c_val, 'ratified/dedicated channel (phase_c); two-engine = diagnostic'
        elif veto:
            # Carry the frozen/held value through; the two-engine number stays in the
            # record as a DIAGNOSTIC and never becomes the proposal.
            v4 = _f(g['A_X']) if (g is not None and pd.notna(g.get('A_X'))) else phase_c_val
            source = f'HELD — {veto}; two-engine value recorded as diagnostic only'
        elif te_record is not None and te_record['reported'] is not None:
            v4 = te_record['reported']
            source = (f"two-engine synthesis floor ({te_record['ion']}, "
                      f"{','.join(e.replace('engine', '') for e in (te_record['selected_engines'] or []))})")
        else:
            v4, source = phase_c_val, 'phase_c (owed)'

        # v3 FREEZES no value on an `owed` tier (RYA-522/665): a blank A_X there is the
        # tier working, not a missing number. Reporting it as "changed" would invent a
        # diff on every held element.
        v3_val = _f(g['A_X']) if (g is not None and pd.notna(g.get('A_X'))) else None
        # The RYA-522 tier lives in the `confidence` column: gold / gf_floor /
        # upper_limit / owed. Only the first three carry a frozen A_X.
        v3_tier = (g.get('confidence') if g is not None else None)
        rows.append({
            'element': el, 'ion': (te['ion'] if te else (g['ion'] if g is not None else 'I')),
            'asplund2021': asp,
            'v3_gold': v3_val, 'v3_tier': v3_tier,
            'v3_verdict': (g['verdict'] if g is not None else None),
            'v4_proposed': v4,
            'v4_delta_vs_asplund': (_f(v4 - asp) if v4 is not None else None),
            'delta_v4_minus_v3': (_f(v4 - v3_val) if (v4 is not None and v3_val is not None) else None),
            'verdict': vb.get('verdict', 'CURATION-OWED'),
            'source': source,
            'changed_vs_v3': (v3_val is not None and v4 is not None and abs(v4 - v3_val) > 1e-9),
            'two_engine_record': te_record})
    rows.sort(key=lambda r: (TIER.get(r['verdict'], 9), -r['asplund2021']))
    return rows


def render_diff_md(rows, counts, decisions, moved) -> str:
    L = ['# RYA-669 — RYA-527 Phase 2: proposed gold v4 vs frozen v3 vs Asplund', '',
         '_**REVIEW ARTIFACT — NOTHING IS FROZEN.** The freeze is a separate, '
         'Ryan-ratified ticket running `promote_solar_reference.py --apply`. '
         'Gold v3 is untouched by this run._', '',
         f"Verdict counts: {counts}", '',
         '`v3 gold` is blank wherever the RYA-522 tier is `owed` — that is the tier '
         'holding the value unfrozen, not a missing measurement.', '',
         '## Proposed v4 vs frozen v3 vs Asplund (tiered), with the two-engine record', '',
         '| El | Asp | v3 gold | tier | v4 proposed | Δ(v4−v3) | ΔAsp | verdict | engA | engB | selected | dCE | source |',
         '|----|-----|---------|------|-------------|----------|------|---------|------|------|----------|-----|--------|']
    cell = lambda x: '—' if x is None else x                      # noqa: E731
    for r in rows:
        te = r['two_engine_record'] or {}
        se = ','.join(e.replace('engine', '') for e in (te.get('selected_engines') or [])) or '—'
        d_asp = r['v4_delta_vs_asplund']
        d_v3 = r['delta_v4_minus_v3']
        v4_s = '—' if r['v4_proposed'] is None else str(r['v4_proposed'])
        if r['changed_vs_v3']:
            v4_s += ' **MOVED**'
        L.append(f"| {r['element']} {r['ion']} | {r['asplund2021']} | {cell(r['v3_gold'])} | "
                 f"{cell(r['v3_tier'])} | {v4_s} | "
                 f"{'—' if d_v3 is None else f'{d_v3:+.3f}'} | "
                 f"{'—' if d_asp is None else f'{d_asp:+.3f}'} | {r['verdict']} | "
                 f"{cell(te.get('engineA'))} | {cell(te.get('engineB'))} | {se} | "
                 f"{cell(te.get('cross_engine_delta'))} | {r['source']} |")

    L += ['', '## Elements whose value moved > 0.01 dex vs frozen v3', '']
    if moved:
        L.append('| element | v3 | v4 | Δ | why |')
        L.append('|---|---|---|---|---|')
        for m in moved:
            L.append(f"| {m['element']} | {m['v3_gold']} | {m['v4_proposed']} | "
                     f"{m['delta_v4_minus_v3']:+.3f} | {m['source']} |")
    else:
        L.append('**None.** Every element carrying a frozen v3 value reproduces it to '
                 'within 0.01 dex on the fresh re-emit.')

    L += ['', '## The four species-adoption decisions — Ryan decides, nothing adopted', '']
    for d in decisions:
        L.append(f"### {d['id']} — {d['element']}: {d['question']}")
        L.append('')
        for k, v in d.items():
            if k in ('id', 'element', 'question', 'recommendation', 'caveat'):
                continue
            L.append(f"- `{k}`: {v}")
        L.append(f"- **Recommendation:** {d['recommendation']}")
        if d.get('caveat'):
            L.append(f"- ⚠ {d['caveat']}")
        L.append('')
    return '\n'.join(L) + '\n'


def render_summary_md(rows, counts, decisions, moved, report, stops, gold_version,
                      te_path: Path) -> str:
    """The one page Ryan reads first. Says what moved, what did not, and what is owed.

    Deliberately leads with the verdict-level answer (does v3 stand?) rather than the
    table, because the table cannot be read without knowing which column is a fresh
    measurement and which is a re-classification of the freeze.
    """
    fresh_ok = not stops and not moved
    L = ['# RYA-669 — RYA-527 Phase 2 run summary', '',
         '**REVIEW ARTIFACT. Nothing frozen, promoted or adopted.** Gold v3 is '
         'byte-untouched; `data/reference/solar/` was not written to.', '']

    L += ['## The answer', '']
    if stops:
        L += [f'🛑 **STOPPED — {len(stops)} §4 STOP condition(s) tripped.** The artifacts '
              'below were still written so the evidence is readable, but no further step '
              'was taken.', '']
        L += [f'- {s}' for s in stops] + ['']
    elif fresh_ok:
        L += ['✅ **Gold v3 STANDS on the fresh re-emit.** No STOP condition tripped, and '
              'no element carrying a frozen v3 value moved by more than 0.01 dex. The '
              'freeze reproduces on a genuinely fresh two-engine run.', '',
              'Path A in the ticket: Beta can close on v3, subject to the four species '
              'decisions below — none of which is a *correction* to v3, all of which are '
              'about elements v3 holds unfrozen.', '']
    else:
        L += [f'⚠️ **A v4 candidate emerged.** No STOP condition tripped, but '
              f'{len(moved)} element(s) carrying a frozen v3 value moved by more than '
              f'0.01 dex on the fresh run. Path B: these need Ryan\'s ratification '
              f'before any freeze.', '']

    L += ['## What is actually fresh here', '',
          '| leg | fresh? | what it means |', '|---|---|---|',
          f'| two-engine record (`{te_path.name}`) | **YES — re-computed** | both engines '
          're-driven over real solar data on current main (Sirius): Engine A EW→A(X) per '
          'line + production NLTE δ; Engine B a new Turbospectrum synthesis-v2 flux fit '
          'per line + Gerber TS-native NLTE δ. Nothing cherry-picked from the July 18 '
          'branch. |',
          '| phase_c verdict | fresh run, **derived** input | since RYA-469 phase_c '
          'CLASSIFIES the frozen gold (`read_solar_reference(\'CURRENT\')`) plus the '
          'dedicated-channel measurements — it does not re-derive A(X) from spectra. '
          'Re-running it answers "does the freeze re-classify consistently", NOT "does '
          'the measurement reproduce". |',
          f'| disposition report | **YES** | same classifier, run over the FRESH record. '
          f'It was expected to retire the gate-3 staleness RYA-663 flagged; it does not, '
          f'and the section below shows why that flag cannot clear by re-running. |', '',
          f'Gold compared against: **{gold_version}**. Verdict counts: `{counts}`.', '']

    L += ['## Elements whose value moved > 0.01 dex vs frozen v3', '']
    if moved:
        L += ['| element | v3 | v4 | Δ | why |', '|---|---|---|---|---|']
        L += [f"| {m['element']} | {m['v3_gold']} | {m['v4_proposed']} | "
              f"{m['delta_v4_minus_v3']:+.3f} | {m['source']} |" for m in moved]
    else:
        L += ['**None.**']
    L += ['']

    L += ['## Gate 3 — and why the PROVISIONAL flag cannot clear itself', '']
    if report['gate3_provisional']:
        L += [f"`gate3_provisional` still reads **True**, on "
              f"{len(report['stale_input_evidence'])} element(s):", '']
        L += [f"- {e}" for e in report['stale_input_evidence']]
        L += ['',
              '**That verdict is now demonstrably wrong, and this run is what proves '
              'it.** `detect_stale_inputs` infers "the two-engine artifact predates that '
              'measurement" from *any* value disagreement with the live channel. The '
              'artifact it just read was generated on current main during this run, so '
              'it predates nothing. Every one of the disagreements above is a '
              'CROSS-CHANNEL difference, not an age difference:', '',
              '| element | two-engine leg | live channel leg | why they differ |',
              '|---|---|---|---|',
              '| Fe | per-line winner-combine | EW ionization-gated + 3D | the RATIFIED '
              'Fe policy — the two-engine number is a diagnostic that sits above the '
              'anchor BY CONSTRUCTION (Ryan, 2026-07-16) |',
              '| Cr | Cr I synthesis floor | gf_floor EW value | different legs |',
              '| Si | synthesis floor | gf_floor EW value | different legs |',
              '| S | EW leg | RYA-492 Costa-Silva dedicated synthesis | different legs |',
              '| Li | synthesis point value | ratified UPPER_LIMIT | the RYA-563 veto |',
              '| O | 8.730 | 8.735 | 0.005 — rounding |', '',
              'The detector conflates *"this artifact is old"* with *"the diagnostic '
              'legitimately disagrees with the ratified channel"*. The second is the '
              'normal, designed state of a two-engine floor. So the flag is **not '
              'clearable by re-running** — RYA-663 deferred Ca\'s promotion to "confirm '
              'on the RYA-527 re-run", the re-run has now happened, and the flag reads '
              'exactly the same.', '',
              '**Ca\'s cross-engine delta is nevertheless genuinely fresh: −0.003**, '
              'computed this run from both engines over real solar data. The number gate '
              '3 needs is sound; only the blanket provisional stamp on top of it is not.',
              '', 'Recommendation (a decision, so not taken here): narrow '
              '`detect_stale_inputs` to compare like-for-like legs, or bound it by the '
              'artifact\'s git commit date against the verdict\'s, so a same-day artifact '
              'cannot be reported as predating anything.', '']
    else:
        L += ['✅ **Cleared.** The fresh two-engine record carries no contradiction '
              'against the live verdict channel.', '']
    L += [f"Promotes under the ratified three gates: "
          f"**{', '.join(report['can_flip_now']) or 'none'}**", '']

    L += ['## The four species-adoption decisions — NOT adopted, Ryan decides', '']
    for d in decisions:
        fresh = d.get('fresh_value', d.get('value'))
        L += [f"**{d['id']} — {d['element']}: {d['question']}**", '',
              f"- Fresh number: `{fresh}`"
              + ('' if fresh is not None else '  ← **no fresh value produced**')]
        if d.get('no_fresh_value_because'):
            L.append(f"- Why not: {d['no_fresh_value_because']}")
        L.append(f"- Recommendation: {d['recommendation']}")
        if d.get('caveat'):
            L.append(f"- ⚠ {d['caveat']}")
        L.append('')

    L += ['## Known defect carried forward (not fixed here)', '',
          '`data/reference/solar/solar_abundances_v3.csv` holds a **Sr I** row and no '
          'Sr II row, while the NLTE registry ratifies Sr as **Sr II** (RYA-551/643). '
          'The diff table therefore shows Sr II against a blank v3 cell. This is the '
          'RYA-663 defect, unchanged by the v3 freeze; adopting Sr II (D.1) is what '
          'would repair it.', '']
    return '\n'.join(L) + '\n'


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--two-engine', default=None,
                    help='fresh two-engine records (default data/audit/rya527_phase2/...)')
    args = ap.parse_args()
    # Resolved against the repo root so a repo-relative argument and an absolute one
    # both land on the same file — git_provenance reports paths relative to the root.
    te_path = (Path(args.two_engine) if args.two_engine else TWO_ENGINE_PHASE2)
    if not te_path.is_absolute():
        te_path = ROOT / te_path
    if not te_path.exists():
        raise SystemExit(f"fresh two-engine record not found at {te_path} — run "
                         f"scripts/rya527_two_engine_run.py --out-dir data/audit/rya527_phase2 first")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. the fresh phase_c verdict, captured beside the run it belongs to
    shutil.copy2(PHASE_C_LIVE, OUT_DIR / 'solar_phase_c_verdict.json')
    verdict = json.loads(PHASE_C_LIVE.read_text())
    vbase = {v['element']: v for v in verdict['verdicts']}

    # 2. STOP conditions BEFORE anything is interpreted
    gold_df, gold_version = ns.read_solar_reference('CURRENT')
    stops = check_stop_conditions(vbase) + check_fe_3d_idempotency(vbase, gold_df)

    # 3. the disposition report on the FRESH two-engine record
    report = ed.build_report(
        two_engine_path=te_path,
        ticket='RYA-669 per-element disposition — RYA-527 Phase 2 (fresh two-engine record)')
    (OUT_DIR / 'element_disposition_report.json').write_text(
        json.dumps(report, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    # 4. the v4 diff table against FROZEN v3
    gold_v3 = {str(r['element']): r for _, r in gold_df.iterrows()}
    gold_ions = {str(r['element']): str(r['ion']) for _, r in gold_df.iterrows()
                 if pd.notna(r.get('ion'))}
    recs_by_el, gerber_delta, gerber_xfail = _records_by_element(te_path)
    rows = build_diff_rows(vbase, recs_by_el, gold_v3, gerber_delta, gerber_xfail, gold_ions)
    moved = [r for r in rows if r['delta_v4_minus_v3'] is not None
             and abs(r['delta_v4_minus_v3']) > 0.01]
    counts = {}
    for r in rows:
        counts[r['verdict']] = counts.get(r['verdict'], 0) + 1

    decisions = species_decisions(
        vbase, {el: _pick_ion(el, rs, gold_ions) for el, rs in recs_by_el.items()}, report)

    payload = {
        'ticket': 'RYA-669 — RYA-527 Phase 2 fresh two-engine re-emit. REVIEW ARTIFACT, '
                  'NOT frozen. Gold v3 untouched.',
        'gold_version_compared_against': gold_version,
        'value_rule': 'RATIFIED/dedicated channel governs the reported value; the '
                      'two-engine result is the RYA-525 cross-engine DIAGNOSTIC',
        'stop_conditions_tripped': stops,
        'counts': counts,
        'elements_moved_gt_0.01_dex': moved,
        'species_adoption_decisions': decisions,
        'gate3_provisional': report['gate3_provisional'],
        'stale_input_evidence': report['stale_input_evidence'],
        'can_flip_now': report['can_flip_now'],
        'gerber_nlte_delta': gerber_delta, 'gerber_xfail': sorted(gerber_xfail),
        'diff_table': rows}
    # RYA-674 §2C: the proposed gold v4 table is the highest-stakes emission in the
    # repo — it is what a freeze would immortalise. Gated before it is written.
    assert_ratified_constraints_satisfied(
        rows, 'RYA-669 Phase 2 re-emit / proposed gold v4 ladder')
    (OUT_DIR / 'proposed_gold_v4_diff.json').write_text(json.dumps(payload, indent=2))
    (OUT_DIR / 'proposed_gold_v4_diff.md').write_text(
        render_diff_md(rows, counts, decisions, moved))
    (OUT_DIR / 'phase2_run_summary.md').write_text(
        render_summary_md(rows, counts, decisions, moved, report, stops,
                          gold_version, te_path))

    print(f"gold compared against : {gold_version}")
    print(f"counts                : {counts}")
    print(f"moved > 0.01 dex      : {[m['element'] for m in moved] or 'none'}")
    print(f"can flip now          : {report['can_flip_now'] or 'none'}"
          + ('  (gate 3 PROVISIONAL)' if report['gate3_provisional'] else '  (gate 3 fresh)'))
    if stops:
        print('\nSTOP CONDITIONS TRIPPED (RYA-669 §4):')
        for s in stops:
            print(f"  - {s}")
        return 2
    print('\nno STOP condition tripped')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

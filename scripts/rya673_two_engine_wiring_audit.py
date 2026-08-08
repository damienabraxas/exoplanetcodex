#!/usr/bin/env python3
"""
scripts/rya673_two_engine_wiring_audit.py
=========================================
RYA-673 — per-species Engine-A / Engine-B WIRING audit across the canonical set.

DISCOVERY ONLY. This script changes no wiring, no engine, no verdict and no state
surface. It reads the orchestrator's own coverage functions and reports what they
actually cover, so the per-element wiring tickets are filed against a known map
instead of a guess.

WHY AN AUDIT WAS NEEDED
-----------------------
RYA-669 found that `_dedicated_engine_B()` never wired the Co and Ba synthesis
harnesses, so Phase 2 could not arbitrate Co and left Ba's gate 3 UNEVALUABLE. That
raised the real question: for how many species has the RYA-525 "two-engine floor"
actually been a ONE-engine floor? Nobody knew, because nothing ever asked.

HOW IT ASKS — the orchestrator's own functions, not a reimplementation
----------------------------------------------------------------------
The audit imports `scripts/rya527_two_engine_run.py` and calls its three coverage
producers directly:

    _engine_A_perline()    Engine A — EW pool -> A(X) per line + production NLTE delta
    _engine_B_perline()    Engine B route 1 — the synthesis-v2 flux fit, per line
    _dedicated_engine_B()  Engine B route 2 — the dedicated synthesis harnesses

A second copy of the wiring logic would be worse than no audit: it could report a
coverage the pipeline does not have. So "wired" here means *that function returned a
value for this species on real solar data*, which is the only definition that cannot
drift from the pipeline.

TWO THINGS THE TICKET ASSUMED THAT ARE NOT TRUE OF THE CODE
------------------------------------------------------------
Recorded here because they change what the audit can even mean:

1. There is **no `_dedicated_engine_A()`**, and the dedicated engines do not live in
   `pipeline/engine_selection.py` — `_dedicated_engine_B()` is in
   `scripts/rya527_two_engine_run.py`. `engine_selection.py` is the *selector*
   (it picks between two values it is handed); it does no wiring at all.
2. The two engines are not symmetric, so "wired" does not mean the same thing on
   each side. Engine A is POOL-driven: a species is covered iff the curated EW pool
   has lines for it. Engine B has TWO independent routes (synth-v2, dedicated
   harness) and a species can be covered by either.

REASON TAXONOMY
---------------
The ticket's five classes, plus one the ticket's list cannot express. An Engine-A
gap caused by an empty EW pool is neither "no atom" nor "no grid" nor "no harness
call" — it is a measurement gap, and calling it UNKNOWN would hide a cause we know
exactly. `NO_EW_POOL` is therefore added rather than mislabelling those rows.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import NLTE_CORRECTION_ELEMENTS, TARGET_ELEMENTS  # noqa: E402
import pipeline.problem_children as pc                                  # noqa: E402
from pipeline.engine_selection import (                                 # noqa: E402
    is_upper_limit_disposition, nlte_atom_validation, ratified_reported_ion)

CSV_OUT = ROOT / 'data' / 'audit' / 'two_engine_wiring_audit.csv'
MD_OUT = ROOT / 'docs' / 'audit' / 'two_engine_wiring_report.md'
ORCHESTRATOR = ROOT / 'scripts' / 'rya527_two_engine_run.py'
PHASE_C = ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_phase_c_verdict.json'


def phase_c_channels() -> dict[str, tuple]:
    """element -> (A_measured, channel) from the live verdict.

    THE THIRD CHANNEL. Neither engine is the whole story: `phase_c` also reads
    dedicated channels the orchestrator never sees — chiefly the RYA-460 Kitt Peak
    atlas. K I is the sharpest case: it is a gold-tier PASS at 5.099, and it is
    `neither`-wired here. Without this column the audit would read as "K is
    unmeasured", which is false and would send someone to fix the wrong thing.

    The real finding for those species is subtler and worse: they are measured on
    exactly ONE channel that the two-engine floor cannot see, so the floor can never
    cross-check them. That is an architectural gap, not a measurement gap.
    """
    if not PHASE_C.exists():
        return {}
    doc = json.loads(PHASE_C.read_text(encoding='utf-8'))
    return {v['element']: (v.get('A_measured'), str(v.get('channel') or ''))
            for v in doc.get('verdicts', [])}

# ── reason classes ───────────────────────────────────────────────────────────
NO_MODEL_ATOM = 'NO_MODEL_ATOM'
NO_NLTE_GRID = 'NO_NLTE_GRID'
NO_HARNESS_INVOCATION = 'NO_HARNESS_INVOCATION'
DELIBERATELY_SKIPPED = 'DELIBERATELY_SKIPPED'
NO_EW_POOL = 'NO_EW_POOL'
UNKNOWN = 'UNKNOWN'

BOTH, A_ONLY, B_ONLY, NEITHER = 'both', 'A_only', 'B_only', 'neither'

_ION_NUMERAL = {1: 'I', 2: 'II', 3: 'III'}

#: Dedicated Engine-B synthesis harness RESULTS that exist on disk, per element, with
#: the ticket that produced them. A file present here but not referenced by the
#: orchestrator source is the NO_HARNESS_INVOCATION class — a measurement the Codex
#: already owns and does not use. Detection is by reading the orchestrator's source,
#: so wiring one up flips its row automatically with no edit here.
HARNESS_RESULTS = {
    'C':  ('data/audit/cno_synthesis/solar_phase_a_cross_arm.json', 'RYA-237/491'),
    'O':  ('data/audit/cno_synthesis/solar_phase_a_cross_arm.json', 'RYA-237/491'),
    'Mn': ('data/audit/mn_hfs_synthesis/solar_mn_hfs_synthesis_rya473.json', 'RYA-473'),
    'Cu': ('data/audit/cu_v_hfs_synthesis/solar_cu_v_hfs_synthesis_rya466.json', 'RYA-466'),
    'V':  ('data/audit/cu_v_hfs_synthesis/solar_cu_v_hfs_synthesis_rya466.json', 'RYA-466'),
    'Sr': ('data/results/sr2_synthesis_rya551.json', 'RYA-551/643'),
    'Zr': ('data/results/zr2_synthesis_rya560.json', 'RYA-560'),
    'Mg': ('data/results/mg_5528_synthesis_rya592.json', 'RYA-592'),
    'Co': ('data/results/co_synthesis_rya564.json', 'RYA-564'),
    'Ba': ('data/results/solar_ba_synthesis_rya559.json', 'RYA-559'),
    'S':  ('data/results/solar_s_costasilva_rya492.json', 'RYA-492'),
}


def _load_orchestrator():
    """Import the orchestrator as a module so its coverage functions can be called."""
    spec = importlib.util.spec_from_file_location('rya527_two_engine_run', ORCHESTRATOR)
    mod = importlib.util.module_from_spec(spec)
    sys.modules['rya527_two_engine_run'] = mod
    spec.loader.exec_module(mod)
    return mod


def canonical_species() -> list[tuple[str, str]]:
    """The canonical audit rows, DERIVED — never a hardcoded list.

    One row per element on its ratified/registry-reported ion, plus **Fe II**: Fe is
    the one element the Codex reports on two ions, because Fe II is the RYA-406
    ionization-balance arbiter and is wired independently of Fe I. That is what makes
    the row count 27 against 26 canonical elements — the ticket's "27 elements" counts
    Fe twice, and this is where the extra row comes from.
    """
    rows = []
    for el in TARGET_ELEMENTS:
        locked = ratified_reported_ion(el)
        if locked is not None:
            ion = _ION_NUMERAL.get(locked, 'I')
        else:
            disp = pc.disposition_for(el) or {}
            sp = str(disp.get('species') or '').split()
            ion = sp[1] if len(sp) > 1 else 'I'
        rows.append((el, ion))
        if el == 'Fe':
            rows.append(('Fe', 'II'))
    return rows


def synthesis_required(mod) -> set[str]:
    """The RYA-520 set, taken from the orchestrator's own rule, not restated.

    These elements have their Engine-A raw-EW leg DELIBERATELY suppressed: a
    saturated / HFS / blended species must not report a raw-EW abundance (that is how
    C=10.26 happened). Their Engine B is therefore not optional — it is the only leg
    they have.
    """
    return {'C', 'N', 'O'} | {
        el for el in TARGET_ELEMENTS
        if (d := pc.disposition_for(el)) and
        d.get('required_treatment') in ('synthesis', 'HFS_sum')}


def audit(mod) -> list[dict]:
    src = ORCHESTRATOR.read_text(encoding='utf-8')
    params = mod._solar_params()
    a_cov = mod._engine_A_perline(params)      # {(el, ion): {...lines}}
    b_cov = mod._engine_B_perline()            # {(el, ion): {...lines}}
    ded_b = mod._dedicated_engine_B()          # {(el, ion): (value, source)}
    synth_req = synthesis_required(mod)
    pcc = phase_c_channels()

    rows = []
    for el, ion in canonical_species():
        key = (el, ion)
        disp = pc.disposition_for(el) or {}
        treatment = str(disp.get('required_treatment') or 'none')
        is_synth_req = el in synth_req

        # ── Engine A ────────────────────────────────────────────────────────
        a_lines = len(a_cov.get(key, {}))
        if is_synth_req:
            # Not a gap. The suppression is the RYA-520 fix working as designed.
            a_wired, a_reason = False, DELIBERATELY_SKIPPED
            a_note = (f'raw-EW leg suppressed by RYA-520 (required_treatment='
                      f'{treatment}) — Engine B is this species\' only valid leg')
        elif a_lines:
            a_wired, a_reason, a_note = True, '', f'{a_lines} EW-pool line(s)'
        else:
            a_wired, a_reason = False, NO_EW_POOL
            a_note = 'no line survives the curated EW pool for this species'

        # Engine-A NLTE grid presence is a QUALITY axis, not a wiring one: no grid
        # means the leg runs in LTE, it does not mean the leg is missing.
        a_nlte = el in NLTE_CORRECTION_ELEMENTS

        # ── Engine B ────────────────────────────────────────────────────────
        b_synth = len(b_cov.get(key, {}))
        b_ded = key in ded_b
        harness_path, harness_ticket = HARNESS_RESULTS.get(el, (None, None))
        # "Referenced" = the orchestrator's source actually names the result file.
        harness_referenced = bool(harness_path) and Path(harness_path).name in src
        harness_exists = bool(harness_path) and (ROOT / harness_path).exists()

        atom_ok, atom_citation = nlte_atom_validation(el)

        if b_ded or b_synth:
            b_wired, b_reason = True, ''
            route = 'dedicated harness' if b_ded else 'synth-v2 per-line'
            b_note = f'Engine B via {route}'
            if b_ded and b_synth:
                b_note = 'Engine B via dedicated harness (synth-v2 also covers it)'
            # A ratified dedicated measurement that exists and is NOT called is still
            # a wiring gap even when synth-v2 happens to cover the species.
            if harness_exists and not harness_referenced:
                b_note += (f'; ⚠ the ratified {harness_ticket} harness result '
                           f'({harness_path}) is NOT referenced by the orchestrator')
        elif harness_exists and not harness_referenced:
            b_wired, b_reason = False, NO_HARNESS_INVOCATION
            b_note = (f'the {harness_ticket} synthesis result already exists at '
                      f'{harness_path} — the orchestrator never reads it')
        elif not atom_ok and el not in HARNESS_RESULTS:
            b_wired, b_reason = False, NO_MODEL_ATOM
            b_note = f'no validated Engine-B NLTE atom: {atom_citation}'
        elif b_ded and not b_synth:
            b_wired, b_reason = False, UNKNOWN
            b_note = 'wired but produced nothing — see the gated-harness note'
        else:
            b_wired, b_reason = False, UNKNOWN
            b_note = ('no synth-v2 lines and no dedicated harness result on record — '
                      'cause not determinable from the code; needs Ryan review')

        # A harness that is wired AND gated shut (Zr's reliability floor, Mg's
        # concordance band) is a wired leg producing a null, not an unwired one —
        # "wire it" is the wrong follow-on for those. But this may only downgrade a
        # species that has NO other Engine-B route: Mg is covered by synth-v2
        # independently of its gated 5528 harness, and an earlier version of this
        # check reported Mg as `neither` on the strength of the gated harness alone,
        # hiding a working Engine-B leg.
        if (el in ('Zr', 'Mg') and harness_referenced and not b_ded and not b_synth):
            b_wired, b_reason = False, DELIBERATELY_SKIPPED
            b_note = (f'harness IS wired ({harness_ticket}) and is gated shut: it '
                      f'returns nothing until its reliability/concordance condition '
                      f'is met. Wiring is NOT the fix — the measurement is.')

        if is_upper_limit_disposition(el):
            b_note += '; ratified UPPER_LIMIT (RYA-563) — may never carry a point value'

        status = (BOTH if (a_wired and b_wired) else
                  A_ONLY if a_wired else
                  B_ONLY if b_wired else NEITHER)

        # The RYA-525 loud-fail was written to make exactly this unrepresentable.
        blocks_beta = is_synth_req and not b_wired

        # Is this species measured at all, and if so on a channel neither engine sees?
        pc_val, pc_channel = pcc.get(el, (None, ''))
        measured = pc_val is not None
        off_channel = measured and status == NEITHER
        if off_channel:
            b_note += (f'; ⚠ MEASURED OFF-ORCHESTRATOR: phase_c carries A={pc_val} via '
                       f'"{pc_channel[:60]}" — a channel the two-engine floor never '
                       f'sees, so this value has NO cross-engine confirmation and '
                       f'never can until it is wired')

        rows.append({
            'element': el, 'ion': ion,
            'engine_a_wired': a_wired, 'engine_a_reason': a_reason,
            'engine_b_wired': b_wired, 'engine_b_reason': b_reason,
            'wiring_status': status,
            'required_treatment': treatment,
            'engine_a_nlte_grid': a_nlte,
            'engine_b_atom_validated': atom_ok,
            'engine_a_lines': a_lines, 'engine_b_synthv2_lines': b_synth,
            'dedicated_harness_exists': harness_exists,
            'dedicated_harness_wired': harness_referenced,
            'harness_ticket': harness_ticket or '',
            'synthesis_required': is_synth_req,
            'blocks_beta': blocks_beta,
            'phase_c_measured': measured,
            'phase_c_value': pc_val,
            'measured_off_orchestrator': off_channel,
            'phase_c_channel': pc_channel,
            'notes': f'A: {a_note}. B: {b_note}',
        })
    return rows


def render_md(rows: list[dict]) -> str:
    counts = {s: sum(1 for r in rows if r['wiring_status'] == s)
              for s in (BOTH, A_ONLY, B_ONLY, NEITHER)}
    blockers = [r for r in rows if r['blocks_beta']]
    unwired_b = [r for r in rows if not r['engine_b_wired']]
    silent = [r for r in rows if r['wiring_status'] == NEITHER]

    L = ['# Two-engine wiring audit — all canonical species (RYA-673)', '',
         '**GENERATED — do not hand-edit.** Regenerate with '
         '`python scripts/rya673_two_engine_wiring_audit.py`.', '',
         '**Discovery only.** No wiring was changed, no engine touched, no verdict '
         'regenerated. Every "wired" below means *the orchestrator\'s own coverage '
         'function returned a value for this species on real solar data* — not that a '
         'code path appears to exist.', '',
         f"`both` **{counts[BOTH]}** · `A_only` **{counts[A_ONLY]}** · "
         f"`B_only` **{counts[B_ONLY]}** · `neither` **{counts[NEITHER]}** "
         f"— {len(rows)} species", '']

    L += ['## The headline', '',
          f'**{len(unwired_b)} of {len(rows)} species have no Engine B.** Of those, '
          f'**{len(blockers)} are synthesis-required** — elements whose raw-EW leg is '
          'deliberately suppressed by RYA-520, so Engine B is not their second opinion, '
          'it is their *only* leg. Those rows are reporting on one engine or on nothing.',
          '']
    if blockers:
        L += ['| species | treatment | why Engine B is missing |', '|---|---|---|']
        L += [f"| {r['element']} {r['ion']} | `{r['required_treatment']}` | "
              f"`{r['engine_b_reason']}` |" for r in blockers]
        L.append('')

    L += ['## ⚠ The RYA-525 loud-fail does not cover this', '',
          'RYA-525 added a guard whose stated job is to refuse a synthesis-required '
          'element with no Engine-B value:', '',
          '```python',
          'species = sorted(set(a_pl) | set(b_pl) | set(ded_b), ...)',
          'for (el, ion) in species:',
          '    ...',
          '    if synth_required and not has_B:',
          '        loud.append(...)          # -> SystemExit',
          '```', '',
          'The guard iterates the **union of the three coverage sources**. An element '
          'absent from all three never enters the loop, so it is never tested — it is '
          'skipped in silence. The guard catches a *partially* covered species and is '
          'blind to a *completely* uncovered one, which is the more serious case.', '']
    if silent:
        L += [f"Species that produce **no two-engine record at all**, and therefore "
              f"never reach the guard: **{', '.join(r['element'] + ' ' + r['ion'] for r in silent)}**.",
              '']

    off = [r for r in rows if r['measured_off_orchestrator']]
    unmeasured = [r for r in rows if r['wiring_status'] == NEITHER
                  and not r['phase_c_measured']]
    L += ['## `neither` splits in two, and the halves need opposite responses', '',
          'A species wired to no engine is not necessarily unmeasured. `phase_c` reads '
          'dedicated channels the orchestrator never sees — chiefly the RYA-460 Kitt '
          'Peak atlas. Reading `neither` as "unmeasured" would send someone to fix the '
          'wrong thing entirely.', '']
    if off:
        L += [f'### Measured, but invisible to the floor ({len(off)})', '',
              'These carry a real value on exactly ONE channel, and the two-engine floor '
              'cannot see it. **They have no cross-engine confirmation and cannot '
              'acquire one until they are wired** — for Beta\'s "best of abilities on '
              'all engines" bar, this is the important class, and it is invisible in '
              'every existing report.', '',
              '| species | value | channel |', '|---|---|---|']
        L += [f"| {r['element']} {r['ion']} | {r['phase_c_value']} | "
              f"{r['phase_c_channel'][:70]} |" for r in off]
        L.append('')
    if unmeasured:
        L += [f'### Genuinely unmeasured ({len(unmeasured)})', '',
              f"**{', '.join(r['element'] + ' ' + r['ion'] for r in unmeasured)}** — no "
              f"engine and no value anywhere. These need a measurement, not wiring.", '']

    L += ['## Per-species', '',
          '| species | status | A | B | A reason | B reason | treatment | blocks Beta |',
          '|---|---|---|---|---|---|---|---|']
    tick = lambda b: '✓' if b else '✗'                                  # noqa: E731
    for r in rows:
        L.append(f"| {r['element']} {r['ion']} | `{r['wiring_status']}` | "
                 f"{tick(r['engine_a_wired'])} | {tick(r['engine_b_wired'])} | "
                 f"{r['engine_a_reason'] or '—'} | {r['engine_b_reason'] or '—'} | "
                 f"`{r['required_treatment']}` | {'**YES**' if r['blocks_beta'] else '—'} |")

    L += ['', '## What each reason class means, and the follow-on it implies', '', ]
    for cls, meaning in (
        (NO_HARNESS_INVOCATION,
         'the measurement **already exists** in the repo and the orchestrator simply '
         'never reads it. Cheapest possible class — a wiring ticket, no new science.'),
        (NO_MODEL_ATOM,
         'no validated Engine-B NLTE atom for the species. Needs atom sourcing + an '
         'RYA-534/548-style anchor validation BEFORE wiring is meaningful.'),
        (NO_NLTE_GRID,
         'no 1D delta CSV for the Engine-A leg. Affects quality (the leg runs LTE), '
         'not wiring — no row is marked unwired for this alone.'),
        (NO_EW_POOL,
         'the curated EW pool has no surviving line for the species. A measurement '
         'gap, not a plumbing gap — line-pool / gf work, not a wiring ticket.'),
        (DELIBERATELY_SKIPPED,
         'a documented, ratified decision — the RYA-520 raw-EW suppression, or a '
         'harness that is wired and gated shut. **Not** a gap; do not file a ticket.'),
        (UNKNOWN,
         'the audit could not determine the cause from the code. Deferred to Ryan by '
         'design rather than guessed.'),
    ):
        members = [f"{r['element']} {r['ion']}" for r in rows
                   if cls in (r['engine_a_reason'], r['engine_b_reason'])]
        L.append(f"- **`{cls}`** ({', '.join(members) if members else 'none'}) — {meaning}")

    L += ['', '## Recommended follow-on per non-`both` species', '',
          'Recommendations only. **No tickets filed** — per §4, Ryan directs which get '
          'filed.', '',
          '| species | recommended action |', '|---|---|']
    for r in rows:
        if r['wiring_status'] == BOTH:
            continue
        reasons = {r['engine_a_reason'], r['engine_b_reason']} - {''}
        if r['blocks_beta'] and r['engine_b_reason'] == NO_HARNESS_INVOCATION:
            act = ('**FILE A WIRING TICKET** — result exists, orchestrator does not '
                   'read it. Blocks Beta: this species currently reports on no engine.')
        elif r['blocks_beta']:
            act = (f"**BLOCKS BETA** — synthesis-required with no Engine B "
                   f"(`{r['engine_b_reason']}`). Needs the underlying gap closed first.")
        elif reasons == {DELIBERATELY_SKIPPED}:
            act = 'None — ratified decision, working as designed.'
        elif r['engine_b_reason'] == NO_HARNESS_INVOCATION:
            act = ('File a wiring ticket (low cost, result already exists). Not '
                   'Beta-blocking — the species has another leg.')
        elif r['measured_off_orchestrator']:
            act = (f"Measured off-orchestrator (A={r['phase_c_value']}) with NO "
                   f"cross-engine confirmation. Wiring it into the floor is what makes "
                   f"that value confirmable — a Beta-quality question, not a "
                   f"measurement one.")
        elif NO_EW_POOL in reasons:
            act = 'Line-pool / gf work, not wiring. No wiring ticket.'
        elif r['engine_b_reason'] == NO_MODEL_ATOM:
            act = 'Atom sourcing + anchor validation first; wiring only after.'
        else:
            act = 'Needs Ryan review (UNKNOWN cause).'
        L.append(f"| {r['element']} {r['ion']} | {act} |")

    L += ['', '## Per-species narrative', '']
    for r in rows:
        L += [f"### {r['element']} {r['ion']} — `{r['wiring_status']}`", '',
              f"{r['notes']}", '']
    return '\n'.join(L) + '\n'


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--check', action='store_true',
                    help='verify the committed audit matches a fresh run; write nothing')
    args = ap.parse_args()

    mod = _load_orchestrator()
    rows = audit(mod)

    if args.check:
        if not CSV_OUT.exists():
            print(f'MISSING: {CSV_OUT.relative_to(ROOT)}', file=sys.stderr)
            return 1
        committed = pd.read_csv(CSV_OUT).to_dict('records')
        fresh = pd.DataFrame(rows).to_dict('records')
        if json.dumps(committed, sort_keys=True, default=str) != json.dumps(
                fresh, sort_keys=True, default=str):
            print(f'STALE: {CSV_OUT.relative_to(ROOT)} does not match a fresh run',
                  file=sys.stderr)
            return 1
        print('Two-engine wiring audit is up to date.')
        return 0

    CSV_OUT.parent.mkdir(parents=True, exist_ok=True)
    MD_OUT.parent.mkdir(parents=True, exist_ok=True)
    with CSV_OUT.open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    MD_OUT.write_text(render_md(rows), encoding='utf-8')

    counts = {s: sum(1 for r in rows if r['wiring_status'] == s)
              for s in (BOTH, A_ONLY, B_ONLY, NEITHER)}
    print(f'Wrote {CSV_OUT.relative_to(ROOT)} and {MD_OUT.relative_to(ROOT)}')
    print(f'  {len(rows)} species: ' + ' · '.join(f'{k}={v}' for k, v in counts.items()))
    blockers = [f"{r['element']} {r['ion']}" for r in rows if r['blocks_beta']]
    print(f'  synthesis-required with NO Engine B (Beta-blocking): '
          f'{", ".join(blockers) if blockers else "none"}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

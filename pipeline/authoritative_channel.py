#!/usr/bin/env python3
"""
pipeline/authoritative_channel.py  (RYA-521)
============================================
ONE authoritative solar-abundance channel: the **phase_c verdict**
(synthesis + HFS + atlas + NLTE + curation + the RYA-463 disposition source).

`run()`'s raw EW output (`solar_abundances.csv`) is **DIAGNOSTIC-ONLY** — it is a
per-line EW cross-check, never the reported/consumed abundance for any element.
Everything downstream (gold reference / RYA-469 freeze, website, C/O ratio, the
CODEX state register) must read the abundance from the verdict via
`load_verdict_abundances()`, NOT from the raw EW file.

Why this exists: the two channels can disagree. On the flagship they did by 1.77
dex — raw EW C = 10.260 (saturated C I 5380) vs verdict C = 8.491 (CH synthesis),
RYA-520. Until there is one documented source, "the abundance of X" is ambiguous
and that class recurs. `channel_divergence()` is the guard that surfaces it.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# The single source of truth. Consumers assert against this.
AUTHORITATIVE_CHANNEL = 'phase_c_verdict'
DIVERGENCE_TOL_DEX = 0.10   # per-element raw-EW vs verdict tolerance (RYA-371 TOL_PASS)

DIAGNOSTIC_ONLY_HEADER = (
    "# DIAGNOSTIC-ONLY (RYA-521): raw per-line EW abundances. NOT the authoritative\n"
    "# solar abundance for any element — the authoritative value is the phase_c\n"
    "# verdict (pipeline/authoritative_channel.load_verdict_abundances). Synthesis-\n"
    "# required elements (C/N/O, HFS heavies) are suppressed here (RYA-520).\n"
)


def _verdict_path(star: str = 'solar', path=None) -> Path:
    """Canonical per-star verdict json (RYA-469 namespaced), audit copy as fallback."""
    if path is not None:
        return Path(path)
    cand = [
        ROOT / 'data' / 'outputs' / star / f'{star}_verdict.json',
        ROOT / 'data' / 'audit' / 'cno_synthesis' / f'{star}_phase_c_verdict.json',
        ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_phase_c_verdict.json',
    ]
    for c in cand:
        if c.exists():
            return c
    raise FileNotFoundError(
        f"No phase_c verdict json for star={star!r} (looked: {[str(c) for c in cand]}). "
        f"Run scripts/phase_c_verdict_rya371.py first — the verdict is the authoritative "
        f"channel (RYA-521).")


def load_verdict_abundances(star: str = 'solar', path=None) -> dict:
    """{element: {'A': float|None, 'verdict': str, 'channel': str, 'owed': str}} from
    the authoritative phase_c verdict. `A` is the verdict's A_measured (None where the
    element is owed and carries no point value). THE authoritative reader."""
    vj = json.loads(_verdict_path(star, path).read_text())
    out = {}
    for r in vj.get('verdicts', []):
        a = r.get('A_measured')
        out[str(r['element'])] = {
            'A': (float(a) if a is not None else None),
            'verdict': r.get('verdict'),
            'channel': r.get('channel'),
            'owed': r.get('owed'),
        }
    return out


def load_raw_ew_abundances(path) -> dict:
    """{element: A} from the raw DIAGNOSTIC-ONLY run output, comparable to the verdict:
    the dominant/neutral ion (the verdict reports one abundance per element, from the
    workhorse ion) and the NLTE-applied value where present (else 1D-LTE A_X) — so the
    guard compares like-with-like and does not false-positive on LTE-vs-NLTE or an Fe I
    vs Fe II ion mismatch. Rows suppressed by RYA-520 (authoritative=False) are skipped."""
    df = pd.read_csv(path, comment='#')
    out = {}
    for el, g in df.groupby(df['element'].astype(str)):
        if 'authoritative' in g.columns:
            g = g[g['authoritative'] != False]
        if g.empty:
            continue
        if 'ion' in g.columns:
            neutral = g[g['ion'].astype(str).str.strip().isin(('I', '1', '1.0'))]
            g = neutral if not neutral.empty else g
        r = g.iloc[0]
        a = r.get('A_X_nlte')
        if a is None or pd.isna(a):
            a = r.get('A_X')
        if pd.notna(a):
            out[el] = float(a)
    return out


def channel_divergence(star: str = 'solar', raw_path=None, verdict_path=None,
                       tol: float = DIVERGENCE_TOL_DEX) -> list[dict]:
    """Per-element |raw-EW − verdict| beyond `tol`. Divergence is EXPECTED (raw EW is
    diagnostic), but a large delta on a PASS element is a signal (as C=10.260 was).
    Returns a list of dicts (sorted worst-first); each carries `pass_element` = the
    verdict is PASS, which escalates it. Compares only where BOTH channels have a value."""
    verdict = load_verdict_abundances(star, verdict_path)
    if raw_path is None:
        raw_path = ROOT / 'data' / 'outputs' / star / f'{star}_abundances.csv'
    raw = load_raw_ew_abundances(raw_path)
    flags = []
    for el, ra in raw.items():
        v = verdict.get(el)
        if v is None or v['A'] is None:
            continue
        delta = abs(ra - v['A'])
        if delta > tol:
            flags.append({
                'element': el, 'raw_ew': round(ra, 3), 'verdict': round(v['A'], 3),
                'delta': round(delta, 3), 'verdict_status': v['verdict'],
                'pass_element': v['verdict'] == 'PASS',
            })
    flags.sort(key=lambda d: (not d['pass_element'], -d['delta']))
    return flags


def assert_authoritative_is_verdict(source: str) -> None:
    """Consumer guard: a downstream that reports/consumes an abundance must source it
    from the verdict channel. Call with the channel you're reading; raises if it is the
    raw EW file. Makes 'read the verdict, not the raw EW' enforceable in code."""
    s = str(source).lower()
    if 'abundances.csv' in s and 'reference' not in s:
        raise RuntimeError(
            f"RYA-521: {source!r} is the DIAGNOSTIC-ONLY raw EW channel. The authoritative "
            f"solar abundance is the phase_c verdict ({AUTHORITATIVE_CHANNEL}) — read it via "
            f"authoritative_channel.load_verdict_abundances(), never the raw EW file.")

#!/usr/bin/env python3
"""
scripts/measure_cu_v_hfs_synthesis_rya466.py
============================================
RYA-466 — MEASURE solar Cu via HFS-resolved synthesis (+ RYA-402 b-factor NLTE);
V as an LTE stretch (the NLTE-void).

THE BLOCKER (from RYA-354 Batch 2, corrected): Cu and V are NOT gf-grade-limited.
They have ZERO measured lines reaching the EW curation pool — the single-profile EW
fit cannot measure a hyperfine-split line, so Cu/V never enter the pool gf-grading
operates on. Same root cause as Eu 6645 / Li 6707: the WRONG MEASUREMENT TOOL for an
HFS element. gf adjudication structurally cannot help when there is no measured line.

THE FIX (reuse, not rebuild): synthesise. The GES atomic line list the production
synthesis path already loads (`GESv6_atom_hfs_iso`, RYA-353) carries the Cu/V hyperfine
components NATIVELY — each feature is 4..10 resolved components, NOT a collapsed single
line — so a flux-space synthesis fit (`cno_synthesis._fit_element`, the RYA-338 synthesis
v2 machinery) handles the HFS the EW path could not. This routes solar Cu/V through that
existing path, varying only A(X) per line and reading χ²ᵣ against the observed HARPS solar
flux. No new synthesis engine; no new line data.

gf ADJUDICATION (Cu, the RYA-354 conflict — resolved here):
  RYA-354 flagged a genuine VALD-vs-GES disagreement on Cu I (median |Δ|~0.275 dex, all
  5 lines). The two scales are recorded IN-REPO in canonical_gf.csv (gf_synth_ges vs
  gf_linelist_vald). The GES column is the Kock & Richter gf (reference_code 'KR') — the
  NIST-graded lab authority RYA-354 named (Kock & Richter 1968, Z. Astrophys.) — and it is
  the gf the synthesis ACTUALLY USES (the synthesis reads the GES list). The adjudication
  therefore ADOPTS GES = Kock & Richter and SUPERSEDES the VALD3 value that canonical_gf
  records for 4 of the 5 lines (the 5th, 5782, already carries KR). This is cited, never
  invented: the value is read from the in-repo GES list, not transcribed from a paper
  (the WebFetch-hallucination rule). We write adjudication_status ONLY (value/grade/ref
  UNTOUCHED — the synthesis already uses KR), the RYA-354 Batch-2 precedent.

NLTE:
  * Cu — RYA-402 built the b-factor (departure-coefficient) Cu NLTE model. The solar Cu
    correction it derived is SMALL and POSITIVE: single-component +0.001, HFS-resolved
    +0.003 dex (consistent with Shi et al. 2014, small positive optical Cu). We TRY the
    live PySME derivation (pipeline.pysme_nlte.nlte_delta('Cu')); if the multi-GB .grd
    departure grid is not on disk in this worktree we APPLY the RYA-402-derived value as a
    vendored correction, LOUDLY FLAGGED — never silently LTE, never re-fitted to Asplund.
  * V — the NLTE-VOID (no model atom anywhere; RYA-463 registry). LTE-only, flagged,
    lower confidence. Named gf sources for a future graded pull: V I Lawler+2014 ApJS
    215,20; V II Wood+2014 ApJS 214,18.

VALIDATE-DON'T-TUNE: gf cited/adjudicated (KR, in-repo), NLTE vendored (RYA-402), no value
pulled toward Asplund. A synthesised Cu that still sits high is a FINDING, not a failure.

Out:
  data/audit/cu_v_hfs_synthesis/solar_cu_v_hfs_synthesis_rya466.json  (phase_c folds this)
  data/audit/cu_v_hfs_synthesis/solar_cu_v_hfs_synthesis_rya466.csv   (per-line table)
  data/linelists/canonical_gf.csv  (adjudication_status on the 5 Cu I lines; value frozen)

    python -m scripts.measure_cu_v_hfs_synthesis_rya466            # full run
    python -m scripts.measure_cu_v_hfs_synthesis_rya466 --quick    # Cu only, no V
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import warnings
from datetime import date
from pathlib import Path

import numpy as np

warnings.filterwarnings('ignore')

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from config.constants import ISPEC_DIR, HARPS_R, SOLAR_ASPLUND2021, get_star_params  # noqa: E402

sys.path.insert(0, str(ISPEC_DIR))
import ispec  # noqa: E402

from pipeline.abundances_derive import (  # noqa: E402
    _load_atmosphere, _load_synth_resources, _load_observed_spectrum,
    _ISPEC_SOLAR_ABUND_FILE,
)
from pipeline.cno_synthesis import _atom_codes, _solar_A, _fit_element  # noqa: E402

CANON = _REPO / 'data' / 'linelists' / 'canonical_gf.csv'
OUT_DIR = _REPO / 'data' / 'audit' / 'cu_v_hfs_synthesis'
TMP_DIR = '/tmp/ispec_cu466'

# ── Diagnostic lines (air Å) and per-line fit half-windows (Å) ────────────────
# Cu I — the 5 solar diagnostic lines RYA-354 flagged (all HFS-split in GES).
CU_LINES = [
    (5105.537, 0.45),
    (5218.198, 0.45),
    (5220.066, 0.45),
    (5700.237, 0.45),
    (5782.122, 0.45),
]
# V I — clean low-EP solar lines (HFS-split in GES). V II is dominated by weak
# blue/UV lines in the Sun (not reachable cleanly in the HARPS-VIS optical) → V I only.
V_LINES = [
    (5727.048, 0.30),
    (6039.722, 0.30),
    (6081.441, 0.30),
    (6090.214, 0.30),
    (6251.827, 0.30),
    (6531.415, 0.30),
]

# RYA-402 derived solar Cu b-factor NLTE correction (HFS-resolved 10-component value);
# small + positive, consistent with Shi et al. 2014. Applied as the vendored value when
# the live .grd grid is not on disk (loudly flagged). Validate-don't-tune.
CU_NLTE_DELTA_DOC = 0.003
CU_NLTE_REF = ('RYA-402 b-factor (departure-coefficient) Cu model, HFS-resolved solar '
               'delta +0.003 dex; cross-check Shi et al. 2014 (small positive optical Cu).')

_GF_COLS = None  # filled by _canon_index()


def _solar_params() -> dict:
    rec = get_star_params('solar')
    return {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
            'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}, rec


# ── HFS verification (the feature is resolved, not a collapsed single line) ────

def verify_hfs(ll, element_ion: str, lines) -> dict:
    """For each feature, count the GES hyperfine components within ±0.25 Å. Returns
    {wl: n_components}. The whole premise is that these are >1 (HFS-resolved); a feature
    that collapses to 1 component would mean the EW path could have measured it."""
    notes = np.array([str(x) for x in ll['element']])
    w = np.asarray(ll['wave_A'], float)
    out = {}
    for wl, _ in lines:
        m = np.where((notes == element_ion) & (np.abs(w - wl) < 0.25))[0]
        out[round(wl, 3)] = int(len(m))
    return out


# ── gf adjudication (Cu): adopt GES = Kock & Richter, supersede VALD3 ──────────

def _canon_index():
    global _GF_COLS
    rows = list(csv.reader(open(CANON)))
    hdr = rows[0]
    _GF_COLS = {n: i for i, n in enumerate(hdr)}
    return hdr, rows


def adjudicate_cu_gf(apply: bool = True) -> list:
    """Resolve the Cu I VALD-vs-GES gf conflict to the NIST-graded Kock & Richter value
    (already the synthesis gf, GES column), recording adjudication_status on the 5 lines.
    Value/grade/reference are FROZEN (the synthesis already uses KR; we never invent or
    transcribe a gf). Returns the per-line audit (GES vs VALD Δ).

    BYTE-SAFE textual line-edit (the RYA-354 Batch-2 precedent): we read the file as raw
    lines, replace ONLY the adjudication_status field on the 5 target rows, and write the
    bytes back — every other row is preserved byte-for-byte. (A csv.writer round-trip would
    rewrite all 145k line endings; the canonical_gf mixed-type columns also corrupt under
    pandas to_csv — hence the textual edit.)"""
    hdr, rows = _canon_index()      # parsed view, for locating + auditing
    c = _GF_COLS
    raw = CANON.read_text().splitlines(keepends=True)
    audit, edits = [], []
    for wl, _ in CU_LINES:
        match = None
        for i, r in enumerate(rows[1:], start=1):
            if r[c['species']] == 'Cu I' and abs(float(r[c['wavelength_air_A']]) - wl) < 0.02:
                match = i
                break
        if match is None:
            raise SystemExit(f"FATAL: Cu I {wl} not in canonical_gf — refusing to proceed")
        r = rows[match]
        # refuse the textual edit if the row has an embedded comma (field count drift)
        if len(raw[match].rstrip('\n').split(',')) != len(hdr) or len(r) != len(hdr):
            raise SystemExit(f"Cu I {wl}: embedded comma (ncols mismatch) — refusing textual edit")
        ges = r[c['gf_synth_ges']]
        vald = r[c['gf_linelist_vald']] or r[c['gf_regions_vald']]
        try:
            delta = float(ges) - float(vald)
        except ValueError:
            delta = float('nan')
        already_kr = 'KR' in r[c['loggf_reference']]
        audit.append(dict(line=round(wl, 3), gf_GES_KR=ges, gf_VALD3=vald,
                          delta_GES_minus_VALD=round(delta, 3),
                          canonical_adopted=r[c['log_gf']],
                          canonical_ref=r[c['loggf_reference']],
                          synthesis_uses='GES = Kock & Richter (reference_code KR)',
                          disposition=('already KR' if already_kr else
                                       'VALD3 SUPERSEDED by GES=KR (Kock & Richter)')))
        edits.append(match)

    if apply:
        for i in edits:
            fields = raw[i].rstrip('\n').split(',')
            assert len(fields) == len(hdr), "field-count drift — refusing edit"
            fields[c['adjudication_status']] = 'adjudicated_kr1968_rya466'
            assert not any(',' in f for f in fields), "new field contains a comma"
            eol = '\n' if raw[i].endswith('\n') else ''
            raw[i] = ','.join(fields) + eol
        CANON.write_text(''.join(raw))
    return audit


# ── flux-space synthesis fit (single free element) ────────────────────────────

def _measure(ow, of, atm, params, element, lines, ll, iso, sab, broad, verbose=True):
    codes = _atom_codes((element,), _CHEM, sab)
    A0 = _solar_A((element,), _CHEM, sab)
    state = dict(A0)
    per_line = []
    for wl, hw in lines:
        t0 = time.time()
        r = _fit_element(ow, of, atm, params, element, state, codes,
                         ((wl - hw, wl + hw),), False, broad,
                         A0[element] - 1.2, A0[element] + 1.2,
                         ll, iso, sab, TMP_DIR)
        r['line'] = round(wl, 3)
        r['wall_s'] = round(time.time() - t0, 1)
        per_line.append(r)
        if verbose:
            print(f"    {element} {wl:9.3f}  A={r['A_X']}  χ²ᵣ={r['red_chi2']}  "
                  f"σfit={r['sigma_fit']}  npix={r['n_pix']}  [{r['status']}]  ({r['wall_s']}s)")
    return per_line


def _summarise(per_line):
    vals = [r['A_X'] for r in per_line if isinstance(r['A_X'], float) and np.isfinite(r['A_X'])]
    if not vals:
        return dict(n=0, mean=None, median=None, scatter=None)
    return dict(n=len(vals), mean=round(float(np.mean(vals)), 3),
                median=round(float(np.median(vals)), 3),
                scatter=round(float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0, 3))


def cu_nlte_delta() -> dict:
    """Try the live RYA-402 PySME b-factor Cu delta; fall back to the RYA-402-derived
    vendored value if the .grd grid is offline. Never silent, never re-fitted."""
    try:
        from pipeline.pysme_nlte import nlte_delta
        res = nlte_delta('Cu')
        return dict(delta=round(float(res['delta_median']), 3), live=True,
                    source='RYA-402 PySME b-factor Cu (live, HFS-resolved); '
                           + CU_NLTE_REF)
    except Exception as exc:                                    # grid offline / no pysme
        return dict(delta=CU_NLTE_DELTA_DOC, live=False,
                    source=f'{CU_NLTE_REF} (live grid unavailable in worktree: '
                           f'{type(exc).__name__}; vendored value applied, FLAGGED).')


# ── main ──────────────────────────────────────────────────────────────────────

_CHEM = None


def main(argv=None):
    global _CHEM
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true', help='Cu only (skip V)')
    ap.add_argument('--no-adjudicate', action='store_true',
                    help='do not write adjudication_status to canonical_gf')
    args = ap.parse_args(argv)

    Path(TMP_DIR).mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    params, rec = _solar_params()
    print(f"\n{'='*74}\n  RYA-466 — solar Cu (+V) via HFS-resolved synthesis "
          f"(Teff={params['teff_K']:.0f} logg={params['logg']:.2f} "
          f"[Fe/H]={params['feh']:+.2f} ξ={params['vturb_kms']:.2f})\n{'='*74}")

    ll, iso, _CHEM = _load_synth_resources()
    sab = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    ow, of = _load_observed_spectrum('solar')
    atm = _load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    broad = (float(HARPS_R), float(rec.get('vmac', 1.6) or 1.6), float(rec.get('vsini', 1.9) or 1.9))

    # ── 1. gf adjudication (Cu) ───────────────────────────────────────────────
    print("\n  ── Cu I gf adjudication (VALD-vs-GES → Kock & Richter) ──")
    gf_audit = adjudicate_cu_gf(apply=not args.no_adjudicate)
    absd = [abs(a['delta_GES_minus_VALD']) for a in gf_audit if np.isfinite(a['delta_GES_minus_VALD'])]
    for a in gf_audit:
        print(f"    {a['line']:9.3f}  GES(KR)={a['gf_GES_KR']:>8}  VALD3={a['gf_VALD3']:>8}  "
              f"Δ={a['delta_GES_minus_VALD']:+.3f}  → {a['disposition']}")
    print(f"    median |Δ(GES−VALD)| = {np.median(absd):.3f} dex over {len(absd)} lines; "
          f"adjudicated to GES = Kock & Richter (KR) — the gf the synthesis uses.")

    # ── 2. HFS verification ───────────────────────────────────────────────────
    cu_hfs = verify_hfs(ll, 'Cu 1', CU_LINES)
    print("\n  ── HFS resolution check (GES components per feature) ──")
    print(f"    Cu I: {cu_hfs}  (all >1 → hyperfine-resolved, NOT collapsed)")
    assert all(n > 1 for n in cu_hfs.values()), "a Cu feature collapsed to <=1 component"

    # ── 3. Cu HFS synthesis fit ───────────────────────────────────────────────
    print("\n  ── Cu I HFS-resolved synthesis fit (flux space) ──")
    cu_pl = _measure(ow, of, atm, params, 'Cu', CU_LINES, ll, iso, sab, broad)
    cu_sum = _summarise(cu_pl)
    nd = cu_nlte_delta()
    a_lte = cu_sum['median']
    a_nlte = round(a_lte + nd['delta'], 3) if a_lte is not None else None
    asp_cu = SOLAR_ASPLUND2021['Cu']
    print(f"    A(Cu)_LTE  median {a_lte} (mean {cu_sum['mean']}, σ {cu_sum['scatter']}, n={cu_sum['n']})")
    print(f"    NLTE Δ {nd['delta']:+.3f} ({'live' if nd['live'] else 'vendored RYA-402'}) "
          f"→ A(Cu)_NLTE {a_nlte}  ({a_nlte - asp_cu:+.3f} vs Asplund {asp_cu:.2f})")

    # ── 4. V LTE fit (NLTE-void) ──────────────────────────────────────────────
    v_pl, v_sum = [], dict(n=0, mean=None, median=None, scatter=None)
    if not args.quick:
        v_hfs = verify_hfs(ll, 'V 1', V_LINES)
        print("\n  ── V I LTE synthesis fit (NLTE-VOID — no model atom; flagged) ──")
        print(f"    V I HFS components per feature: {v_hfs}")
        v_pl = _measure(ow, of, atm, params, 'V', V_LINES, ll, iso, sab, broad)
        v_sum = _summarise(v_pl)
        asp_v = SOLAR_ASPLUND2021['V']
        print(f"    A(V)_LTE  median {v_sum['median']} (mean {v_sum['mean']}, σ {v_sum['scatter']}, "
              f"n={v_sum['n']})  ({(v_sum['median'] or 0) - asp_v:+.3f} vs Asplund {asp_v:.2f}) "
              f"— LTE-only, V NLTE-void")

    # ── 5. write products ─────────────────────────────────────────────────────
    payload = {
        'ticket': 'RYA-466', 'star': 'solar', 'generated': date.today().isoformat(),
        'engine': 'pipeline.cno_synthesis._fit_element (synthesis v2, RYA-338) on '
                  'GESv6_atom_hfs_iso (HFS-resolved); Turbospectrum 1D-LTE flux fit',
        'params': params,
        'Cu': {
            'method': 'HFS-resolved synthesis (flux-space χ²ᵣ fit)',
            'gf_provenance': 'GES = Kock & Richter (reference_code KR; the NIST-graded '
                             'authority RYA-354 named) — VALD3 superseded, '
                             f'median |Δ(GES−VALD)|={round(float(np.median(absd)),3)} dex',
            'gf_adjudication': gf_audit,
            'hfs_components': cu_hfs,
            'A_lte_median': a_lte, 'A_lte_mean': cu_sum['mean'],
            'scatter': cu_sum['scatter'], 'n_lines': cu_sum['n'],
            'nlte_delta': nd['delta'], 'nlte_live': nd['live'], 'nlte_source': nd['source'],
            'A_nlte': a_nlte, 'asplund2021': asp_cu,
            'delta_vs_asplund': round(a_nlte - asp_cu, 3) if a_nlte is not None else None,
            'per_line': [{k: r[k] for k in ('line', 'A_X', 'red_chi2', 'sigma_fit',
                                            'n_pix', 'status')} for r in cu_pl],
        },
        'V': {
            'method': 'LTE synthesis (NLTE-VOID — no model atom; RYA-463 registry)',
            'nlte_void': True,
            'gf_sources_for_future_pull': 'V I Lawler+2014 ApJS 215,20; V II Wood+2014 ApJS 214,18',
            'A_lte_median': v_sum['median'], 'A_lte_mean': v_sum['mean'],
            'scatter': v_sum['scatter'], 'n_lines': v_sum['n'],
            'asplund2021': SOLAR_ASPLUND2021['V'],
            'delta_vs_asplund': round((v_sum['median'] or 0) - SOLAR_ASPLUND2021['V'], 3)
                                if v_sum['median'] is not None else None,
            'per_line': [{k: r[k] for k in ('line', 'A_X', 'red_chi2', 'sigma_fit',
                                            'n_pix', 'status')} for r in v_pl],
        },
    }
    out_json = OUT_DIR / 'solar_cu_v_hfs_synthesis_rya466.json'
    out_json.write_text(json.dumps(payload, indent=2, default=float))

    out_csv = OUT_DIR / 'solar_cu_v_hfs_synthesis_rya466.csv'
    with open(out_csv, 'w', newline='') as fh:
        wr = csv.writer(fh)
        wr.writerow(['element', 'line_A', 'A_X_lte', 'red_chi2', 'sigma_fit', 'n_pix', 'status'])
        for el, pl in (('Cu', cu_pl), ('V', v_pl)):
            for r in pl:
                wr.writerow([el, r['line'], r['A_X'], r['red_chi2'], r['sigma_fit'],
                             r['n_pix'], r['status']])

    print(f"\n  wrote {out_json.relative_to(_REPO)}\n        {out_csv.relative_to(_REPO)}")
    print(f"\n  VERDICT MOVE: Cu/V no-value → MEASURED. "
          f"A(Cu)_NLTE {a_nlte} (+{a_nlte - asp_cu:.3f} vs Asplund) n={cu_sum['n']}; "
          f"A(V)_LTE {v_sum['median']} n={v_sum['n']} (NLTE-void).")
    return payload


if __name__ == '__main__':
    main()

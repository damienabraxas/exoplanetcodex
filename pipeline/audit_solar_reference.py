#!/usr/bin/env python3
"""
pipeline/audit_solar_reference.py
=================================
RYA-459 (under RYA-162) — audit + coverage verification for the solar reference
library, and the provenance gate (UV must be cited-composite, never measured).

This is the SMOKE TEST / verification entry point:

    python -m pipeline.audit_solar_reference --verify

It does NOT re-download anything. It:
  * reports what is staged under data/solar_reference/ (and what is missing),
  * prints the COVERAGE MATRIX — for each key diagnostic region, which atlas covers
    it, at what sampling, MEASURED or CITED,
  * asserts the provenance flags from config.constants.SOLAR_REFERENCE_SPECTRA
    (every UV/composite source is provenance=cited-composite; no UV is 'measured').

The "actual point" of RYA-459 is the coverage matrix + the cited-vs-measured
discipline, not "files downloaded".
"""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from config.constants import SOLAR_REFERENCE_SPECTRA  # noqa: E402

# Key diagnostic regions the library must be judged against (air Angstrom).
# Mirrors the RYA-459 spec list + the RYA-369 N strategy.
KEY_DIAGNOSTICS = [
    ('NH 3360',          3360.0, 'N'),
    ('CN violet 3883',   3883.0, 'N'),
    ('N I 7442/7468',    7455.0, 'N'),
    ('N I 8216/8223',    8219.5, 'N'),
    ('N I 8680-8718',    8699.0, 'N'),
    ('[O I] 6300',       6300.3, 'O'),
    ('O I 777 triplet',  7773.0, 'O'),
    ('P I 10581/10596', 10589.0, 'P'),
    ('K I 7665/7699',    7682.0, 'K'),
    ('Co I 3845',        3845.5, 'Co'),
    ('Sc II 4246',       4246.8, 'Sc'),
]


def _load_provenance(spec):
    """Load the per-source provenance JSON if present; tolerate absence."""
    d = ROOT / spec['path']
    if not d.exists():
        return None, []
    provs = sorted(d.glob('*provenance*.json'))
    out = []
    for p in provs:
        try:
            out.append(json.loads(p.read_text()))
        except Exception:
            pass
    return d, out


def audit():
    print(f"\n{'='*78}\n  RYA-459 — solar reference library audit\n{'='*78}")
    staged, missing = [], []
    kp_diags = None
    uv_range_nm = None
    for key, spec in SOLAR_REFERENCE_SPECTRA.items():
        d, provs = _load_provenance(spec)
        present = d is not None and any(d.glob('*.csv'))
        flag = 'STAGED' if present else ('DEFERRED' if spec.get('status') == 'deferred' else 'MISSING')
        (staged if present else missing).append(key)
        print(f"\n  [{flag:8s}] {key}  ({spec['provenance']})")
        print(f"     path     : {spec['path']}")
        print(f"     coverage : {spec['wavelength_coverage_nm']} nm   res: {spec['resolution']}")
        print(f"     units    : {spec['flux_units']}")
        print(f"     cite     : {spec.get('citation', '?')}")
        if present:
            n_csv = len(list(d.glob('*.csv')))
            print(f"     files    : {n_csv} csv segment(s)")
        for pr in provs:
            if 'extracted_diagnostics' in pr:
                kp_diags = pr['extracted_diagnostics']
            if 'wl_vac_A_range' in pr:
                uv_range_nm = [pr['wl_vac_A_range'][0] / 10.0, pr['wl_vac_A_range'][1] / 10.0]
    return staged, missing, kp_diags, uv_range_nm


def coverage_matrix(kp_diags, uv_range_nm):
    print(f"\n{'='*78}\n  COVERAGE MATRIX — diagnostic region vs reference atlas\n{'='*78}")
    print(f"  {'Diagnostic':18s} {'line_A':>8s}  {'Kitt Peak (MEASURED)':28s}  {'CALSPEC UV (CITED)':20s}")
    print(f"  {'-'*18} {'-'*8}  {'-'*28}  {'-'*20}")
    kp_by_line = {}
    if kp_diags:
        for e in kp_diags:
            kp_by_line[round(float(e['line_A']), 1)] = e
    rows = []
    for name, line_A, el in KEY_DIAGNOSTICS:
        # Kitt Peak: match the extracted diagnostic by nearest line, else range test
        kp = None
        for e in (kp_diags or []):
            lo, hi = e['window_A']
            if lo - 1 <= line_A <= hi + 1:
                kp = e
                break
        if kp and kp.get('covered'):
            dl = kp.get('median_dlambda_A')
            kp_txt = f"YES  n={kp['n_points']} dl={dl}A"
            kp_cov = True
        elif 2960.0 <= line_A <= 13000.0:
            kp_txt = "in-range (not extracted)"
            kp_cov = True
        else:
            kp_txt = "NO (out of 296-1300nm)"
            kp_cov = False
        # CALSPEC composite range (cited)
        uv_cov = bool(uv_range_nm and uv_range_nm[0] <= line_A / 10.0 <= uv_range_nm[1])
        uv_txt = "YES (cited-composite)" if uv_cov else "NO"
        print(f"  {name:18s} {line_A:8.1f}  {kp_txt:28s}  {uv_txt:20s}")
        rows.append({'diagnostic': name, 'line_A': line_A, 'element': el,
                     'kitt_peak_measured': kp_cov, 'calspec_cited': uv_cov})
    # the headline unblock
    n_cov = sum(r['kitt_peak_measured'] for r in rows)
    print(f"\n  -> {n_cov}/{len(rows)} diagnostics have a MEASURED (Kitt Peak) reference.")
    n_unblocked = [r['diagnostic'] for r in rows
                   if r['element'] == 'N' and r['kitt_peak_measured']]
    print(f"  -> solar-N unblock: {len(n_unblocked)} N channels measured "
          f"({', '.join(n_unblocked)}).")
    return rows


def assert_provenance():
    """Gate: every UV/composite source is cited-composite; no UV tagged measured."""
    print(f"\n{'='*78}\n  PROVENANCE GATE\n{'='*78}")
    ok = True
    for key, spec in SOLAR_REFERENCE_SPECTRA.items():
        prov = spec['provenance']
        is_uv = 'uv' in key.lower() or 'composite' in prov
        # the cardinal rule: a UV / composite source must NEVER be tagged 'measured'
        if is_uv and prov == 'measured':
            print(f"  FAIL  {key}: UV/composite source tagged 'measured' — forbidden.")
            ok = False
        else:
            print(f"  ok    {key}: provenance={prov}")
        if prov not in ('measured', 'cited-composite', 'model'):
            print(f"  FAIL  {key}: unknown provenance {prov!r}")
            ok = False
    # explicit: the UV composite is cited-composite
    uv = SOLAR_REFERENCE_SPECTRA.get('uv_composite')
    if uv and uv['provenance'] != 'cited-composite':
        print(f"  FAIL  uv_composite must be 'cited-composite', got {uv['provenance']!r}")
        ok = False
    print(f"\n  PROVENANCE GATE: {'PASS' if ok else 'FAIL'}")
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true', help='run the full audit + gate')
    args = ap.parse_args(argv)
    staged, missing, kp_diags, uv_range_nm = audit()
    rows = coverage_matrix(kp_diags, uv_range_nm)
    gate_ok = assert_provenance()
    print(f"\n{'='*78}")
    print(f"  STAGED: {staged}")
    print(f"  MISSING/DEFERRED: {missing}")
    print(f"{'='*78}\n")
    if args.verify and not gate_ok:
        sys.exit(1)
    return rows, gate_ok


if __name__ == '__main__':
    main()

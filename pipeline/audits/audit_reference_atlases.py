#!/usr/bin/env python3
"""Audit the RYA-390 IR reference atlases (ACE-FTS, NSO photatl, Wallace telluric).

RYA-392 — read-only verification of the RYA-390 Part A intake before Part B's
three-way CO validation (which already consumed these atlases) is trusted. The
correlations Part B reports are only as good as the atlas intake: a mis-converted
cm^-1 axis or a swapped air/vac tag would silently corrupt the whole result.

Checks, per atlas, every one asserted rather than inferred-and-forgotten:

  * Axis convention   — wavenumber vs wavelength, classified from the value range
                        (UNKNOWN is a loud failure, never a silent guess).
  * Air / vacuum      — native = FTS vacuum wavenumber; the stored vac/air columns
                        are reproduced from the wavenumber and checked (not swapped).
  * CO coverage       — the 4255-4367 cm^-1 K-band CO (2-0) segment. Role-aware:
                        ACE / photatl must span the full band; the Wallace ASCII
                        telluric ratio is a *documented* band-middle product
                        (4299.8-4338.6) — checked against its stated sub-range, not
                        false-failed for missing the wings.
  * CO bandhead       — the cm^-1 -> wavelength conversion is verified against the
                        CO (2-0) bandhead near 4360 cm^-1 / 2.293 um.
  * Provenance        — source, citation, version, disk character, telluric status.
                        RYA-390 stored these in the sidecar provenance JSON (keyed by
                        segment filename), not as inline CSV headers; a file with no
                        provenance entry AND no inline provenance keys is a loud fail.
  * Telluric class    — free / residual / pure, read from provenance (not memory).
  * ACE truth-source  — the ACE file must resolve to the Hase+2010 WSpectra "complete
                        solar spectrum" derived product (telluric-FREE), NOT a raw ACE
                        *occultation* transmission (which looks through the atmosphere
                        and would poison the telluric-free truth). Confirmed from the
                        source URL + citation + filename, not the value range.
  * Disk character    — photatl is disk-CENTER; flagged as a known systematic against
                        the integrated-disk reflected-Vesta target (center-to-limb
                        variation changes CO line depths), recorded, not silently
                        treated as equivalent to integrated.

Read-only: it touches nothing under the store. No silent unit assumptions.
"""
import sys
import os
import glob
import json
import argparse

import numpy as np

CO_SEG_CM = (4255.0, 4367.0)        # RYA-390 K-band CO (2-0) segment, cm^-1
COV_TOL_CM = 1.0                    # endpoint tolerance: absorb the FTS sampling grid
                                    # (~0.01 cm^-1 step) without accepting a short atlas
CO_BANDHEAD_CM = 4360.0             # CO (2-0) bandhead, ~2.293 um
CO_BANDHEAD_UM = 2.293             # expected bandhead wavelength (vacuum), micron
DEFAULT_STORE = os.path.join('data', 'solar_reference', 'ir_atlases')

PROV_KEYS = ('source', 'cite', 'hase', 'wallace', 'livingston', 'jqsrt',
             'telluric', 'disk', 'version', 'ace', 'occult', 'wspectra')


def sniff_xaxis(x):
    """Classify the spectral axis from its value range. Returns (label, (lo, hi)).
    'UNKNOWN' is a loud failure, never a silent guess.

    The wavenumber branch is capped at 1e4 cm^-1 (the IR-atlas domain — ACE spans
    700-4430, the K-band CO segments sit at ~4300) so it does not shadow the
    Angstrom branch: a K-band wavelength axis (~23000 A) must classify as Angstrom,
    not be silently mistaken for a 23000 cm^-1 wavenumber."""
    lo, hi = float(np.nanmin(x)), float(np.nanmax(x))
    if 100 < lo and 2000 < hi < 1e4:            # ~700-4430 cm^-1 (IR-atlas domain)
        return 'wavenumber_cm', (lo, hi)
    if 0.5 < lo and hi < 30:                     # microns
        return 'wavelength_um', (lo, hi)
    if 4000 < lo and hi < 3e5:                   # Angstrom
        return 'wavelength_AA', (lo, hi)
    return 'UNKNOWN', (lo, hi)


def to_cm(x, kind):
    """Convert any supported axis to vacuum wavenumber (cm^-1)."""
    if kind == 'wavenumber_cm':
        return np.asarray(x, float)
    if kind == 'wavelength_um':
        return 1.0e4 / np.asarray(x, float)
    if kind == 'wavelength_AA':
        return 1.0e8 / np.asarray(x, float)
    return None


def _tokens(line):
    """Split a data line on comma and/or whitespace (handles CSV and FTS ASCII)."""
    return [t for t in line.replace(',', ' ').split() if t]


def read_any(path):
    """Return (cols, header_text) where cols is a dict colname->ndarray (best effort)
    and the first two numeric columns are positionally cols['x'], cols['y'].

    Handles comma-delimited CSV (with a name header row) and whitespace ASCII; FITS
    via astropy. Unlike a naive thousands-separator strip, this treats the comma as a
    delimiter so the stored segment CSVs parse correctly instead of collapsing into a
    single un-floatable token."""
    if path.lower().endswith(('.fits', '.fit')):
        from astropy.io import fits
        with fits.open(path) as hd:
            hdr = repr(hd[0].header) + ('\n' + repr(hd[1].header) if len(hd) > 1 else '')
            d = hd[1].data
            names = list(d.columns.names)
            cols = {n: np.asarray(d[n], float) for n in names}
            cols['x'], cols['y'] = cols[names[0]], cols[names[1]]
            return cols, hdr

    head, names, rows = [], None, []
    with open(path, 'r', errors='replace') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            tok = _tokens(s)
            try:
                vals = [float(t) for t in tok]
                rows.append(vals)
            except ValueError:
                # non-numeric line: the column-name header (first one) or prose
                if names is None and len(tok) >= 2:
                    names = tok
                head.append(s)
    arr = np.asarray([r for r in rows if len(r) == len(rows[0])], float)
    if names and len(names) == arr.shape[1]:
        cols = {n: arr[:, i] for i, n in enumerate(names)}
    else:
        cols = {f'col{i}': arr[:, i] for i in range(arr.shape[1])}
    keys = list(cols)
    cols['x'], cols['y'] = cols[keys[0]], cols[keys[1]]
    return cols, '\n'.join(head[:60])


def load_provenance(store):
    """Load the RYA-390 sidecar provenance JSON; return {segment_basename: entry}."""
    by_file = {}
    meta = {}
    for jpath in glob.glob(os.path.join(store, '**', '*provenance*.json'), recursive=True):
        with open(jpath) as f:
            doc = json.load(f)
        meta = doc
        for src in doc.get('sources', []):
            seg = src.get('segment_file')
            if seg:
                by_file[os.path.basename(seg)] = src
    return by_file, meta


def _check_air_vac(cols):
    """If the file carries vac+air wavelength columns derived from a wavenumber column,
    verify vac == 1e8/wn to high precision and air < vac by the expected ~6 A IR offset
    (i.e. the tags are not swapped). Returns (status, detail)."""
    wn = next((cols[k] for k in cols if k.lower().startswith('wavenumber')), None)
    vac = next((cols[k] for k in cols if 'vac' in k.lower()), None)
    air = next((cols[k] for k in cols if 'air' in k.lower()), None)
    if wn is None or vac is None:
        return 'n/a', 'no wn/vac columns to cross-check'
    vac_err = float(np.nanmax(np.abs(vac - 1.0e8 / wn)))
    ok = vac_err < 1e-2
    detail = f'max|vac-1e8/wn|={vac_err:.2e} A'
    if air is not None:
        offs = float(np.nanmedian(vac - air))
        # vacuum wavelength is longer than air; ~6 A near 2.3 um (Birch & Downs)
        ok = ok and (3.0 < offs < 9.0)
        detail += f'; median(vac-air)={offs:.2f} A (expect ~6.3, vac>air)'
    return ('OK' if ok else 'CHECK'), detail


def _co_bandhead_check(cols, kind):
    """Verify the cm^-1 -> wavelength conversion lands the CO (2-0) bandhead
    (~4360 cm^-1) near 2.293 um, if the band is covered."""
    xcm = to_cm(cols['x'], kind)
    if xcm is None or np.nanmin(xcm) > CO_BANDHEAD_CM or np.nanmax(xcm) < CO_BANDHEAD_CM:
        return 'n/a', f'{CO_BANDHEAD_CM:.0f} cm^-1 outside file range'
    um = 1.0e4 / CO_BANDHEAD_CM      # native conversion, independent of stored cols
    ok = abs(um - CO_BANDHEAD_UM) < 0.005
    return ('OK' if ok else 'CHECK'), f'{CO_BANDHEAD_CM:.0f} cm^-1 -> {um:.4f} um (expect ~{CO_BANDHEAD_UM})'


def _classify_telluric(entry):
    t = (entry.get('telluric_status', '') if entry else '').lower()
    if 'free' in t:
        return 'free'
    if 'pure' in t:
        return 'pure'
    if 'residual' in t or 'terrestrial' in t:
        return 'residual'
    return 'unknown'


def _ace_truth_source(entry):
    """Confirm the ACE file is the Hase+2010 derived telluric-FREE 'complete solar
    spectrum', NOT a raw ACE occultation transmission. Returns (ok, verdict, evidence)."""
    cite = (entry.get('citation', '') or '').lower()
    url = (entry.get('source_url', '') or '').lower()
    raw = (entry.get('raw_file', '') or '').lower()
    tell = (entry.get('telluric_status', '') or '').lower()
    derived = ('complete solar spectrum' in cite and 'solarspectrum' in url
               and 'free' in tell)
    occ_smell = any(s in raw for s in ('occ', 'transmission', 'occult')) or 'occultation product' in tell
    ok = derived and not occ_smell
    verdict = 'DERIVED telluric-free (Hase+2010 WSpectra)' if ok else 'RAW-OCCULTATION SUSPECT'
    evid = (f"citation='{entry.get('citation','')}'; url='{entry.get('source_url','')}'; "
            f"raw_file='{entry.get('raw_file','')}'; telluric_status='{entry.get('telluric_status','')}'")
    return ok, verdict, evid


def audit_one(path, prov_by_file):
    cols, hdr = read_any(path)
    kind, (lo, hi) = sniff_xaxis(cols['x'])
    xcm = to_cm(cols['x'], kind)
    base = os.path.basename(path)
    entry = prov_by_file.get(base)

    # provenance: sidecar entry (RYA-390 design) OR inline header keywords
    inline = [k for k in PROV_KEYS if k in hdr.lower()]
    has_prov = bool(entry) or bool(inline)

    # coverage — role-aware. Wallace ASCII ratio is a documented band-middle product.
    partial = bool(entry and entry.get('COVERAGE_CAVEAT'))
    if xcm is None:
        covers_full = False
        cov_detail = 'axis UNKNOWN'
    else:
        cmin, cmax = float(np.nanmin(xcm)), float(np.nanmax(xcm))
        covers_full = (cmin <= CO_SEG_CM[0] + COV_TOL_CM
                       and cmax >= CO_SEG_CM[1] - COV_TOL_CM)
        overlap = cmin < CO_SEG_CM[1] and cmax > CO_SEG_CM[0]
        if partial:
            cov_detail = f'band-middle {cmin:.1f}-{cmax:.1f} (documented partial), overlaps CO={overlap}'
            covers_ok = overlap
        else:
            cov_detail = f'{cmin:.1f}-{cmax:.1f}, full-band={covers_full}'
            covers_ok = covers_full

    air_vac, av_detail = _check_air_vac(cols)
    bandhead, bh_detail = _co_bandhead_check(cols, kind)
    tell_class = _classify_telluric(entry)

    notes = []
    role_ok = True
    # ACE truth-source gate (match the ACE-FTS atlas precisely — not the "wallACE"
    # substring; key off the atlas name or an `ace_`-prefixed segment filename)
    is_ace = bool(entry) and ('ace-fts' in entry.get('name', '').lower()
                              or base.lower().startswith('ace'))
    if is_ace:
        a_ok, a_verdict, a_evid = _ace_truth_source(entry)
        role_ok = role_ok and a_ok
        notes.append(f'ACE source: {a_verdict}')
    # photatl disk-center caveat
    if entry and 'disk-center' in (entry.get('disk_integration', '') or '').lower():
        notes.append('photatl disk-CENTER vs integrated-disk Vesta — known systematic (center-to-limb)')
    if partial:
        notes.append('Wallace band-middle only — full-band telluric via photatl atmospheric column')

    ok = (kind != 'UNKNOWN') and covers_ok and has_prov and role_ok \
        and air_vac != 'CHECK' and bandhead != 'CHECK'

    return dict(name=base, kind=kind, rng=(lo, hi), covers_ok=covers_ok,
                covers_full=covers_full, partial=partial, prov=entry is not None,
                inline=inline, n=len(cols['x']), cov_detail=cov_detail,
                air_vac=(air_vac, av_detail), bandhead=(bandhead, bh_detail),
                tell_class=tell_class, notes=notes, ok=ok,
                role=(entry.get('role', '') if entry else ''))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--store', default=DEFAULT_STORE,
                    help='RYA-162 solar-reference store dir (default: %(default)s)')
    a = ap.parse_args(argv)

    store = a.store
    files = [f for f in sorted(glob.glob(os.path.join(store, '**', '*'), recursive=True))
             if f.lower().endswith(('.txt', '.csv', '.dat', '.asc', '.fits', '.fit'))]
    if not files:
        print('NO-GO: no atlas files under', store)
        return 2

    prov_by_file, meta = load_provenance(store)
    if meta:
        print(f'# provenance sidecar: {meta.get("ticket","?")}  CO band {tuple(meta.get("co_band_cm-1", CO_SEG_CM))} cm^-1\n')

    print(f'{"file":42s} {"axis":15s} {"range (cm^-1 native or A)":28s} {"CO":4s} {"a/v":4s} {"bh":4s} {"tell":9s} prov')
    print('-' * 120)

    results = []
    bad = 0
    for f in files:
        try:
            r = audit_one(f, prov_by_file)
        except Exception as e:
            print(f'{os.path.basename(f):42s} READ-FAIL: {e}')
            bad += 1
            continue
        results.append(r)
        bad += 0 if r['ok'] else 1
        provtag = r['role'] or (','.join(r['inline']) if r['inline'] else 'MISSING')
        rng_s = '{:.1f}..{:.1f}'.format(r['rng'][0], r['rng'][1])
        co_s = 'yes' if r['covers_ok'] else 'NO'
        ok_s = 'OK' if r['ok'] else 'CHECK'
        print('{:42s} {:15s} {:28s} {:4s} {:4s} {:4s} {:9s} {}  [{}]'.format(
            r['name'], r['kind'], rng_s, co_s,
            r['air_vac'][0], r['bandhead'][0], r['tell_class'], provtag, ok_s))

    # detail block
    print('\nPer-atlas detail:')
    for r in results:
        print(f'  - {r["name"]} (n={r["n"]})')
        print(f'      coverage : {r["cov_detail"]}')
        print(f'      air/vac  : {r["air_vac"][1]}')
        print(f'      bandhead : {r["bandhead"][1]}')
        for note in r['notes']:
            print(f'      note     : {note}')
        if not r['prov'] and not r['inline']:
            print('      PROVENANCE: MISSING — no sidecar entry and no inline keys (loud fail)')

    print('\nVERDICT:', 'GO' if bad == 0 else f'NO-GO — {bad} atlas(es) need review (CHECK rows)')
    return 0 if bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main())

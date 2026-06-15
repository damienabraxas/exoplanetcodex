#!/usr/bin/env python3
# =============================================================================
# THE EXOPLANET CODEX
# exoplanetcodex.org  |  github.com/damienabraxas/exoplanetcodex
# =============================================================================
#
# File:         vald_inbox_status.py
# Module:       scripts (VALD intake — inbox watcher)
# Description:  READ-ONLY status of the VALD delivery folders. For every archive
#               Ryan has dropped in ~/Documents/Exoplanet Codex/VALD/<star>/ it
#               peeks the header (wavelength band) and the footer model-atmosphere
#               metallicity (Castelli grid node) STRAIGHT FROM THE .gz — it never
#               decompresses to disk, never writes, never touches a committed
#               line list. It then runs the RYA-321 metallicity gate per archive
#               and tells you, per star, whether a build is READY or HELD.
#
#               Purpose: this is the "new-archive alert". Wired as a Claude Code
#               SessionStart hook, it surfaces any freshly-dropped delivery (and
#               whether its composition is correct) the moment a session starts,
#               so the corrected re-extractions (RYA-323) get picked up and built
#               "like we did before" without anyone having to remember to look.
#
#               It is a DETECTOR, not a builder: the actual extraction/merge is
#               still scripts/intake_build_alpha_cen.py (α Cen) / the 55 Cnc merge
#               driver (RYA-324). Run those once this says BUILD READY.
#
# Author:       Ryan Schmitt
# Contributors: Claude (Anthropic) via Claude Code
# Linear issue: RYA-321 — VALD intake metallicity gate
# =============================================================================

import gzip
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from data.linelists.vald_parse import (  # noqa: E402
    _CASTELLI_RE, _nearest_castelli_node)
from config.constants import STAR_PARAMS  # noqa: E402

MET_TOL = 0.02  # matches vald_parse.verify_metallicity default

VALD_ROOT = Path.home() / 'Documents' / 'Exoplanet Codex' / 'VALD'

# VALD delivery folder -> Codex STAR_PARAMS star_id.
FOLDER_STAR = {
    'Sun':              'solar',
    'Procyon':          'procyon',
    'Alpha Centauri A': 'alpha_cen_a',
    'Alpha Centauri B': 'alpha_cen_b',
    '55 Cancri':        '55cnc_a',
}

# Bands a complete per-system extraction must span (matches intake_build_alpha_cen).
BAND_RANGE = {
    'uv1':     (1150.0, 2000.0),
    'uv2':     (2000.0, 3780.0),
    'optical': (3780.0, 6910.0),
    'nir':     (6910.0, 17000.0),
}
RANGE_TOL = 5.0

# Stars whose extraction uses the 4-band scheme above (so "all 4 bands -> BUILD
# READY" applies). Anchors (solar/Procyon) use a different split and are done.
FOUR_BAND_STARS = {'alpha_cen_a', 'alpha_cen_b', '55cnc_a'}


def _open_text(path):
    """Open a delivery as text whether it is gzipped or plain."""
    if path.suffix == '.gz':
        return gzip.open(path, 'rt', errors='replace')
    return open(path, errors='replace')


def peek_band(path):
    """Map a delivery to a band by its header wavelength range (read-only)."""
    try:
        with _open_text(path) as f:
            line1 = f.readline().strip()
            meta = f.readline().strip() if line1.startswith('WARNING') else line1
        lo, hi = (float(x.strip()) for x in meta.split(',')[:2])
    except Exception:
        return None
    return next((b for b, (blo, bhi) in BAND_RANGE.items()
                 if abs(lo - blo) <= RANGE_TOL and abs(hi - bhi) <= RANGE_TOL), None)


def peek_metallicity(path):
    """Castelli [M/H] node VALD applied, read straight from the (gz or plain)
    footer — same token as vald_parse.parse_vald_model_metallicity, but able to
    stream a .gz so nothing is decompressed to disk."""
    try:
        with _open_text(path) as f:
            for line in f:
                m = _CASTELLI_RE.search(line)
                if m:
                    sign = -1.0 if m.group(1).lower() == 'm' else 1.0
                    return sign * int(m.group(2)) / 10.0
    except Exception:
        return None
    return None


def archives(folder_path):
    """Delivery files in a folder: every .gz, plus any extension-less file not
    shadowed by a same-stem .gz (some Procyon drops are uncompressed)."""
    gz = sorted(folder_path.glob('*.gz'))
    gz_stems = {p.stem for p in gz}
    plain = sorted(p for p in folder_path.iterdir()
                   if p.is_file() and p.suffix == '' and p.name not in gz_stems
                   and not p.name.startswith('.'))
    return gz + plain


def verdict_for(path, star_id):
    """(band, mh, (verdict, msg)) for one archive. Applies the same gate as
    vald_parse.verify_metallicity (nearest Castelli node to catalog feh_ref,
    tol MET_TOL), evaluated on the node peeked from the gz/plain delivery."""
    band = peek_band(path)
    mh = peek_metallicity(path)
    rec = STAR_PARAMS.get(star_id)
    if not rec or 'feh_ref' not in rec:
        return band, mh, ('REJECT', f'no STAR_PARAMS[{star_id}].feh_ref')
    if mh is None:
        return band, mh, ('REJECT', 'no model-atmosphere metallicity in delivery')
    catalog = float(rec['feh_ref'])
    expected = _nearest_castelli_node(catalog)
    if abs(mh - expected) > MET_TOL:
        return band, mh, ('REJECT',
                          f'M/H={mh:+.2f} != catalog {catalog:+.2f} (node {expected:+.2f})')
    return band, mh, ('ACCEPT', f'M/H={mh:+.2f} ok vs catalog {catalog:+.2f}')


# "Seen" ledger — archives this watcher has already reported. Kept beside the
# delivery folders (a stable path that survives git-worktree churn, unlike a
# committed file's reset mtime) so --quiet fires exactly once per NEW archive.
LEDGER = VALD_ROOT / '.inbox_seen.json'


def _ident(path):
    """Stable identity for a delivery: parent folder / name / byte size. A
    re-extraction lands under a new request id (new name) or differing size, so
    it reads as a new archive."""
    return f'{path.parent.name}/{path.name}:{path.stat().st_size}'


def _load_seen():
    import json
    try:
        return set(json.loads(LEDGER.read_text()))
    except Exception:
        return set()


def _save_seen(seen):
    import json
    try:
        LEDGER.write_text(json.dumps(sorted(seen)))
    except Exception:
        pass  # ledger is a convenience; never let it break the session hook


def scan():
    """Per-folder scan -> list of (folder, star_id, rows) where each row is
    (path, band, mh, (verdict, msg))."""
    out = []
    for folder, star_id in FOLDER_STAR.items():
        fp = VALD_ROOT / folder
        if not fp.is_dir():
            continue
        files = archives(fp)
        if files:
            out.append((folder, star_id, [(p,) + verdict_for(p, star_id) for p in files]))
    return out


def _star_summary(star_id, rows):
    """One-line verdict for a star's archive set."""
    accepts = [r for r in rows if r[3][0] == 'ACCEPT']
    bands_ok = {r[1] for r in accepts}
    n_rej = len(rows) - len(accepts)
    if n_rej:
        return ('HELD', f'{n_rej} archive(s) fail the metallicity gate — '
                're-extract at correct M/H (RYA-323).')
    if star_id not in FOUR_BAND_STARS:
        return ('OK', 'metallicity OK (anchor star, already correct).')
    if bands_ok == set(BAND_RANGE):
        return ('READY', 'all 4 bands metallicity-correct — run the per-star '
                'intake/build script, then commit.')
    missing = sorted(set(BAND_RANGE) - bands_ok)
    return ('PARTIAL', f'correct bands {sorted(b for b in bands_ok if b)}, '
            f'awaiting {missing}.')


def _print_star(folder, star_id, rows):
    print(f'\n  {folder}  ({star_id})')
    for path, band, mh, (v, msg) in rows:
        mark = '✓' if v == 'ACCEPT' else '✗'
        mhs = f'{mh:+.2f}' if mh is not None else 'none'
        print(f'    {mark} {path.name:24} {str(band):8} M/H={mhs:6} {v}  {msg}')
    state, note = _star_summary(star_id, rows)
    print(f'    -> {state}: {note}')


def main(argv):
    quiet = '--quiet' in argv
    scanned = scan()

    if quiet:
        # Fire only for archives never reported before; then remember them, so a
        # routine session stays silent and a fresh drop alerts exactly once.
        seen = _load_seen()
        fresh = [(folder, sid, rows) for folder, sid, rows in scanned
                 if any(_ident(r[0]) not in seen for r in rows)]
        if not fresh:
            return 0
        print('VALD inbox — NEW archive(s) detected (RYA-321 alert)')
        for folder, sid, rows in fresh:
            _print_star(folder, sid, rows)
        ready = [sid for _, sid, rows in fresh if _star_summary(sid, rows)[0] == 'READY']
        print(f"\n  Action: {'extract now (BUILD READY) — ' + ', '.join(ready) if ready else 'corrected re-extraction still pending (see above)'}.")
        for folder, sid, rows in scanned:
            seen.update(_ident(r[0]) for r in rows)
        _save_seen(seen)
        return 0

    print('VALD inbox — metallicity-gated status (RYA-321)')
    if not scanned:
        print('  (no deliveries in any VALD folder)')
    for folder, sid, rows in scanned:
        _print_star(folder, sid, rows)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

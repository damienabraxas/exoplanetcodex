#!/usr/bin/env python3
# =============================================================================
# RYA-269 Part C — 55 Cnc HFS-ON intake check
# Three new deliveries: UV-a (019516), UV-b (019517), NIR (019518).
# Primary test: split-group signature. Secondary: Cu I direct (UV-b), NIR cross-check.
# =============================================================================

import gzip
import shutil
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from data.linelists.vald_parse import read_vald_header, parse_vald_long, TRUNCATION_WARNING

LINELISTS   = REPO_ROOT / 'data' / 'linelists'
VALD_55CNC  = Path.home() / 'Documents' / 'Exoplanet Codex' / 'VALD' / '55 Cancri'

NEW_GZ = {
    'uv_a': VALD_55CNC / 'RyanSchmitt.019516.gz',
    'uv_b': VALD_55CNC / 'RyanSchmitt.019517.gz',
    'nir':  VALD_55CNC / 'RyanSchmitt.019518.gz',
    'uv_b_dup': VALD_55CNC / 'RyanSchmitt.019520.gz',  # check if duplicate of 017
}
OLD_HFS_OFF_NIR = LINELISTS / 'vald_55cnc_nir_raw.txt'   # 5000–30000, HFS-OFF
OLD_HFS_OFF_UV  = LINELISTS / 'vald_55cnc_uv_raw.txt'    # 019509 1150–3780, HFS-OFF

HFS_SPECIES = {
    ('Mn', 'I'), ('Mn', 'II'),
    ('Co', 'I'), ('Co', 'II'),
    ('Cu', 'I'),
    ('V',  'I'), ('V',  'II'),
    ('Sc', 'I'), ('Sc', 'II'),
    ('Ba', 'II'),
    ('Eu', 'II'),
}


def split_group_signature(records, window=0.2, min_cluster=3):
    hfs = [(r['element'], r['ion'], r['wavelength'])
           for r in records if (r['element'], r['ion']) in HFS_SPECIES]
    by_species = defaultdict(list)
    for elem, ion, wl in hfs:
        by_species[(elem, ion)].append(wl)

    in_groups = 0
    cluster_detail = {}
    for (elem, ion), wls in sorted(by_species.items()):
        wls_sorted = sorted(wls)
        species_groups = 0
        i = 0
        while i < len(wls_sorted):
            j = i
            while j < len(wls_sorted) and wls_sorted[j] - wls_sorted[i] <= window:
                j += 1
            if j - i >= min_cluster:
                species_groups += (j - i)
            i = j
        if species_groups:
            cluster_detail[f'{elem} {ion}'] = species_groups
        in_groups += species_groups

    total = len(hfs)
    frac = in_groups / total if total else 0.0
    return total, in_groups, frac, cluster_detail


def gunzip_to_tempfile(gz_path):
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.txt', mode='wb')
    with gzip.open(gz_path, 'rb') as f_in:
        shutil.copyfileobj(f_in, tmp)
    tmp.close()
    return Path(tmp.name)


def check_one(tag, gz_path, requested_range):
    print(f"\n{'='*60}")
    print(f"  {tag} — {requested_range} Å  ({gz_path.name})")
    print(f"{'='*60}")

    if not gz_path.exists():
        print(f"  ERROR: file not found — {gz_path}")
        return None, 'ERROR'

    tmp = gunzip_to_tempfile(gz_path)
    try:
        # 1. Truncation check
        hdr = read_vald_header(tmp)
        if hdr['truncated']:
            print(f"  TRUNCATION WARNING on line 1 — REJECT")
            return hdr, 'REJECT (truncated)'
        print(f"  Line 1: metadata header (not truncation warning) ✓")
        print(f"  Header: {hdr['wl_start']:.2f}–{hdr['wl_end']:.2f} Å, "
              f"{hdr['n_selected']} selected / {hdr['n_processed']} processed, "
              f"vmicro={hdr['vmicro']}")

        # 2. Parse
        records, report = parse_vald_long(tmp)
        n = report['n_parsed']
        print(f"  Parsed: {n}/{hdr['n_selected']} records "
              f"({'complete' if n == hdr['n_selected'] else 'INCOMPLETE'})")
        if report['n_failures']:
            print(f"  Parse failures: {report['n_failures']}")
            for ex in report['examples'][:3]:
                print(f"    {ex}")

        # 3. Coverage
        wls = [r['wavelength'] for r in records]
        print(f"  Delivered: {min(wls):.3f}–{max(wls):.3f} Å")

        # 4. Split-group signature (PRIMARY HFS test)
        total, in_grp, frac, detail = split_group_signature(records)
        pct = frac * 100
        if total == 0:
            hfs_verdict = 'INCONCLUSIVE (no HFS-capable lines in file)'
        elif frac >= 0.40:
            hfs_verdict = f'ON ✓'
        elif frac <= 0.05:
            hfs_verdict = f'OFF ✗  → REJECT'
        else:
            hfs_verdict = f'AMBIGUOUS — escalate'
        print(f"  Split-group signature: {in_grp}/{total} HFS-capable lines in groups "
              f"({pct:.0f}%) → {hfs_verdict}")
        top = sorted(detail.items(), key=lambda x: -x[1])[:6]
        if top:
            print(f"    Top species: " + ', '.join(f"{s}:{n}" for s, n in top))

        # REJECT if HFS is OFF
        if total > 0 and frac <= 0.05:
            return records, 'REJECT (HFS-OFF)'
        if total > 0 and frac < 0.40:
            return records, 'AMBIGUOUS'

        return records, 'ACCEPT'
    finally:
        tmp.unlink(missing_ok=True)


def cu_i_direct_check(records, tag):
    print(f"\n  Cu I direct check ({tag}):")
    for center, window in [(3247.5, 2.0), (3274.0, 2.0)]:
        nearby = [r for r in records
                  if r['element'] == 'Cu' and r['ion'] == 'I'
                  and abs(r['wavelength'] - center) <= window]
        wls = [f"{r['wavelength']:.3f}" for r in nearby]
        status = '✓ HFS cluster' if len(nearby) >= 3 else ('⚠ weak' if nearby else '✗ absent')
        print(f"    Cu I ~{center} Å: {len(nearby)} records {wls[:8]} → {status}")


def nir_cross_check(new_records, old_path, new_range=(6910, 17000)):
    """Compare new HFS-ON NIR vs old HFS-OFF NIR in the same wavelength slice."""
    print(f"\n  NIR cross-check (new HFS-ON vs old HFS-OFF, {new_range[0]}–{new_range[1]} Å):")
    if not old_path.exists():
        print(f"    Old HFS-OFF NIR not found ({old_path}) — skipping cross-check")
        return

    old_recs, _ = parse_vald_long(old_path)
    old_slice = [r for r in old_recs
                 if new_range[0] <= r['wavelength'] <= new_range[1]]
    new_slice = [r for r in new_records
                 if new_range[0] <= r['wavelength'] <= new_range[1]]

    print(f"    Old HFS-OFF ({old_path.name}), {new_range[0]}–{new_range[1]} Å slice: "
          f"{len(old_slice)} records")
    print(f"    New HFS-ON  (019518),           {new_range[0]}–{new_range[1]} Å:        "
          f"{len(new_slice)} records")
    delta = len(new_slice) - len(old_slice)
    print(f"    Delta: {delta:+d}")

    # Count HFS-capable records in each
    def hfs_count(recs):
        return sum(1 for r in recs if (r['element'], r['ion']) in HFS_SPECIES)
    old_hfs = hfs_count(old_slice)
    new_hfs = hfs_count(new_slice)
    print(f"    HFS-capable records: old={old_hfs}, new={new_hfs}, delta={new_hfs-old_hfs:+d}")

    if delta > 0 and (new_hfs - old_hfs) > 0:
        print(f"    Cross-check: new has MORE records, excess in HFS-capable species ✓")
    elif delta <= 0:
        print(f"    Cross-check: new has FEWER or equal records than old HFS-OFF — WARN")
    else:
        print(f"    Cross-check: more records overall but not concentrated in HFS-capable — WARN")


def main():
    print("RYA-269 Part C — 55 Cnc intake check (HFS-ON re-extractions)")
    print("Primary test: split-group signature\n")

    # Check the three new deliveries
    checks = [
        ('uv_a',    NEW_GZ['uv_a'],     '1150–2000'),
        ('uv_b',    NEW_GZ['uv_b'],     '2000–3780'),
        ('nir',     NEW_GZ['nir'],      '6910–17000'),
    ]

    results = {}
    all_records = {}
    for tag, gz, rng in checks:
        recs, verdict = check_one(tag, gz, rng)
        results[tag] = verdict
        if recs:
            all_records[tag] = recs

    # Extra: Cu I direct check on UV-b
    if 'uv_b' in all_records:
        cu_i_direct_check(all_records['uv_b'], 'uv_b')

    # Extra: NIR cross-check vs old HFS-OFF
    if 'nir' in all_records:
        nir_cross_check(all_records['nir'], OLD_HFS_OFF_NIR)

    # Check 019520 vs 019517 duplicate
    print(f"\n{'='*60}")
    print(f"  019520 duplicate check vs UV-b (019517)")
    print(f"{'='*60}")
    gz520 = NEW_GZ['uv_b_dup']
    if gz520.exists():
        tmp = gunzip_to_tempfile(gz520)
        try:
            hdr520 = read_vald_header(tmp)
            print(f"  019520 header: {hdr520['wl_start']:.2f}–{hdr520['wl_end']:.2f} Å, "
                  f"{hdr520['n_selected']} selected — "
                  f"{'identical to 019517 → duplicate' if hdr520['n_selected'] == 90537 else 'DIFFERENT'}")
        finally:
            tmp.unlink(missing_ok=True)
    else:
        print(f"  019520 not found")

    # Summary table
    print(f"\n{'='*60}")
    print(f"  INTAKE SUMMARY")
    print(f"{'='*60}")
    print(f"  {'File':<8} {'Verdict'}")
    print(f"  {'-'*8} {'-'*30}")
    for tag, verdict in results.items():
        print(f"  {tag:<8} {verdict}")

    all_accept = all(v == 'ACCEPT' for v in results.values())
    print(f"\n  {'ALL ACCEPT — proceed to merge' if all_accept else 'ONE OR MORE REJECT — do not merge'}")
    return 0 if all_accept else 1


if __name__ == '__main__':
    sys.exit(main())

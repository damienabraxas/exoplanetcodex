#!/usr/bin/env python3
"""RYA-963 — register the telluric-corrected α Cen A CRIRES+ holding.

A corrected SIBLING row, following the `solar_harps_molecfit_corrected` precedent
(RYA-931): the uncorrected `alpha_cen_a_crires_plus` EXTRACTC IDPs stay the archival
base and are NOT modified, and the corrected products get their own holding_id.

`telluric_applied` is **derived, not asserted** — `pipeline.telluric_intake` reads it off
the products' own MTRANS extension, and refuses `applied` if that transmission is all
unity. So a frozen molecfit fit (every column at its prior, transmission ≡ 1.0) cannot
reach the registry as a correction even if this script were run on it.

    python3 scripts/rya963_register_corrected.py --products <dir>          # report
    python3 scripts/rya963_register_corrected.py --products <dir> --write  # write

⚠️ The registry is edited with `csv`, never round-tripped through pandas: pandas rewrites
quoting and line endings across every untouched row and buries the real change (RYA-786).
"""
from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import io
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline.telluric_intake import from_many                      # noqa: E402

HOLDINGS = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"

#: One entry per corrected CRIRES+ stellar set. Keyed by the --set the driver used, so
#: the registration cannot drift from the run that produced the products.
SETS = {
    'alpha_cen_a_crires': {
        'holding_id': 'alpha_cen_a_crires_plus_molecfit',
        'base_holding': 'alpha_cen_a_crires_plus',
        'system_id': 'alpha_cen_a', 'ticket': 'RYA-963',
        'audit_dir': 'rya963_crires_telluric',
    },
    'tau_ceti_crires': {
        'holding_id': 'tau_cet_crires_plus_molecfit',
        'base_holding': 'tau_cet_crires_plus',
        'system_id': 'tau_ceti', 'ticket': 'RYA-973',
        'audit_dir': 'rya973_crires_telluric',
    },
}


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rows_from_products(products: list) -> list:
    """Read the manifest rows out of the PRODUCTS THEMSELVES.

    The first cut required the driver's run JSON. That is a side channel, and it went
    missing the moment a frame was re-run on its own: tau Ceti's Y1029 aborted in the set
    run, was retried alone with fewer fit windows, and wrote a perfectly good product
    that the set's JSON knew nothing about. Registering from the run record would have
    silently dropped it. The product is the artifact of record and carries everything the
    manifest needs in its own header, so it is the thing that is read."""
    from astropy.io import fits
    rows = []
    for path in products:
        h = fits.getheader(path)
        def g(k, d=''):
            v = h.get(k, d)
            return '' if v is None else v
        rows.append({
            'product': Path(path).name, 'sha256': sha256(path),
            'base_frame': g('BASEFILE'), 'wlen_id': g('WLEN'), 'band': g('BAND'),
            'date_obs': str(g('DATE-OBS'))[:10],
            'gdas_profile': g('GDAS'), 'gdas_md5': g('GDASMD5'),
            'star_id': g('STARID'),
            # STARIDQ is written only when the id rule is contested; its ABSENCE means
            # the product predates the card, not that the id is settled.
            'star_id_contested': (bool(h['STARIDQ']) if 'STARIDQ' in h else 'unrecorded'),
            'gate_before': f"{float(g('GATEBEF', float('nan'))):.5f}"
                           if isinstance(g('GATEBEF'), (int, float)) else '',
            'gate_after': f"{float(g('GATEAFT', float('nan'))):.5f}"
                          if isinstance(g('GATEAFT'), (int, float)) else '',
            'gate_passed': bool(g('GATEPASS', False)),
        })
    return rows


def build_notes(run: dict, products: list, cfg: dict) -> str:
    """The registry note, built from what is actually in hand. Reads the run record when
    one is given and the PRODUCTS otherwise, so a set assembled from a re-run frame still
    describes itself correctly."""
    from astropy.io import fits
    rows = rows_from_products(products)
    settings = ', '.join(sorted({r['wlen_id'] for r in rows if r['wlen_id']}))
    nights = sorted({r['date_obs'] for r in rows if r['date_obs']})
    gdas = sorted({r['gdas_profile'] for r in rows if r['gdas_profile']})
    gates = '; '.join(f"{r['wlen_id']} {r['gate_before']}->{r['gate_after']}"
                      f"{'' if r['gate_passed'] else ' FAIL'}" for r in rows)
    pwv = []
    for path in products:
        h = fits.getheader(path)
        if 'H2OCOLMM' in h and 'WLEN' in h:
            pwv.append(f"{h['WLEN']} {float(h['H2OCOLMM']):.3f}")
    n_pass = sum(1 for r in rows if r['gate_passed'])
    return (
        f"{cfg['ticket']} telluric-corrected sibling of {cfg['base_holding']}; that "
        f"holding is UNMODIFIED. {len(rows)} corrected CRIRES+ products over "
        f"{len(nights)} night(s) ({', '.join(nights)}), settings {settings}, each "
        f"corrected with ESO molecfit (molecfit_model on derived fit windows, then "
        f"molecfit_calctrans over every chip) against the REAL per-night Paranal GDAS "
        f"profile for ITS OWN night -- {len(gdas)} distinct profile(s): "
        f"{', '.join(gdas)}. No standard-atmosphere fallback is in play (the RYA-373 "
        f"failure mode). ONE atmosphere is fitted per exposure and evaluated over all "
        f"chips: one exposure saw one sky. THE FIT IS PROVEN TO HAVE MOVED, not merely "
        f"to have exited 0 -- molecfit writes every product even when mpfit takes its "
        f"status-4 exit with a zero Jacobian, leaving each column at its prior and "
        f"transmission 1.0 everywhere, so assert_fit_moved requires a FITTED molecular "
        f"column with a non-zero uncertainty and best_chi2 below initial_chi2. "
        f"D1 residual gate at telluric-dominated, stellar-clean pixels, before->after: "
        f"{gates} ({n_pass}/{len(rows)} PASS). Fitted PWV (mm): {'; '.join(pwv)}. "
        f"telluric_applied is DERIVED by pipeline.telluric_intake from each product's "
        f"own non-unity MTRANS extension, not asserted here; the registration script "
        f"refuses to write anything intake does not independently call corrected. "
        f"Evidence: data/audit/{cfg['audit_dir']}/; reproduce with "
        f"python3 -m pipeline.crires_stellar_telluric --set <set>."
    ).replace("\n", " ")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", dest="set_name", default="alpha_cen_a_crires",
                    choices=sorted(SETS))
    ap.add_argument("--products", required=True, help="directory of corrected products")
    ap.add_argument("--run-json", default=None,
                    help="optional run_set record; the manifest is derived from the "
                         "PRODUCT HEADERS regardless, so a frame re-run on its own is "
                         "not lost")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()
    cfg = SETS[a.set_name]
    run = json.loads(Path(a.run_json).read_text()) if a.run_json else {
        'frames': [], 'gdas_gate': {'nights': [], 'slots': {}}}
    products = sorted(glob.glob(str(Path(a.products) / "*_telluric.fits")))
    if not products:
        raise SystemExit(f"no corrected products under {a.products}")

    # DERIVE the state from the products; never assert it.
    state, evidence = from_many(products)
    print(f"products: {len(products)}")
    for p in products:
        print(f"  {Path(p).name}  sha256={sha256(p)}")
    print(f"derived telluric_applied = {state}")
    for e in evidence[:2]:
        print(f"  evidence: {e}")
    if state != "applied":
        raise SystemExit(
            f"REFUSING to register: intake derives telluric_applied={state!r} from the "
            f"products themselves. A holding is not corrected because this script says "
            f"so — if the transmission is all unity the fit never moved.")

    manifest_dir = ROOT / "data" / "audit" / cfg['audit_dir']
    manifest = manifest_dir / "corrected_manifest.csv"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    rows = rows_from_products(products)
    fields = ["product", "sha256", "base_frame", "wlen_id", "band", "date_obs",
              "gdas_profile", "gdas_md5", "star_id", "star_id_contested",
              "gate_before", "gate_after", "gate_passed"]
    with open(manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"[manifest] {manifest}")

    rows = list(csv.DictReader(open(HOLDINGS, newline="")))
    fields = list(rows[0].keys())
    row = {"holding_id": cfg['holding_id'], "system_id": cfg['system_id'],
           "instrument_id": "crires_plus",
           "manifest_path": str(manifest.relative_to(ROOT)),
           "evidence_state": "verified", "source_issue_ids": cfg["ticket"],
           "notes": build_notes(run, products, cfg), "telluric_applied": state}
    existing = [i for i, r in enumerate(rows) if r["holding_id"] == cfg["holding_id"]]
    if existing:
        rows[existing[0]] = row
        print(f"UPDATE row {cfg['holding_id']}")
    else:
        rows.append(row)
        print(f"APPEND row {cfg['holding_id']}")

    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=fields, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    if a.write:
        HOLDINGS.write_text(buf.getvalue())
        print(f"[written] {HOLDINGS}")
    else:
        print("(report only; pass --write to edit the registry)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

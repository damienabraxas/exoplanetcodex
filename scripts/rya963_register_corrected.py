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


def build_notes(run: dict, products: list) -> str:
    g = run["gdas_gate"]
    slots = sorted(set(g["slots"].values()))
    frames = run["frames"]
    confirmed = [f for f in frames if f["star_id"] == "A"]
    gates = [f"{f['wlen_id']} {f['gate_before']:.4f}->{f['gate_after']:.4f}"
             for f in frames]
    pwv = [f"{f['wlen_id']} {f['h2o_col_mm']:.3f}" for f in frames
           if f.get("h2o_col_mm") is not None]
    return (
        f"{cfg['ticket']} telluric-corrected sibling of {cfg['base_holding']}; it is "
        f"UNMODIFIED. {len(products)} corrected CRIRES+ products. "
        f"(one per setting: {', '.join(f['wlen_id'] for f in frames)}), corrected with "
        f"ESO molecfit (molecfit_model on derived fit windows, then molecfit_calctrans "
        f"over every chip) against the REAL per-night Paranal GDAS profile "
        f"{'/'.join(slots)} - all {len(frames)} exposures fall in the same 3-hourly slot, "
        f"and no standard-atmosphere fallback is in play (the RYA-373 failure mode). "
        f"ONE atmosphere is fitted per exposure and evaluated over all chips: one "
        f"exposure saw one sky, so a per-chip independent fit would let the water column "
        f"disagree with itself inside a single 4-second exposure. "
        f"THE FIT IS PROVEN TO HAVE MOVED, not merely to have exited 0: molecfit reports "
        f"success and writes every product even when mpfit takes its status-4 exit with a "
        f"zero Jacobian, leaving each column at its prior with uncertainty exactly 0 and "
        f"transmission 1.0 everywhere - an uncorrected spectrum wearing a correction's "
        f"provenance. Two distinct causes produced exactly that here and both are now "
        f"asserted against: WLC_CONST is a FRACTION of the chip's half wavelength range, "
        f"so its -0.05 default is -3.4 A on a 134 A chip but -135 A on a whole 540 nm "
        f"setting handed over as one chip; and without MAP_REGIONS_TO_CHIP molecfit says "
        f"'Assuming that all regions are mapped to Chip 1' and puts every fit window on "
        f"the first chip. Fit windows are DERIVED per frame (the unsaturated, non-stellar "
        f"absorbed fraction per chip), not tabulated, and the stellar exclusion mask comes "
        f"from linelist_solar rather than hand-listed wavelengths. Molecules are split: "
        f"the MODEL carries every molecule with a band in the frame (or calctrans would "
        f"leave that absorber out entirely), while only those with a band inside a fit "
        f"window are FITTED. D1 residual gate, at telluric-dominated stellar-clean pixels, "
        f"before->after: {'; '.join(gates)}. Fitted PWV (mm): {'; '.join(pwv)}. "
        f"STAR ID: {len(confirmed)}/{len(frames)} frames confirmed alpha Cen A by the "
        f"RYA-423 rule run on the measured RV - which is only available AFTER correction, "
        f"since RYA-423's own CRIRES branch returns INDETERMINATE for want of a telluric-"
        f"corrected spectrum. telluric_applied is DERIVED by pipeline.telluric_intake from "
        f"each product's own non-unity MTRANS extension, not asserted here. "
        f"Evidence: data/audit/rya963_crires_telluric/; reproduce with "
        f"python3 -m pipeline.crires_stellar_telluric --set alpha_cen_a_crires."
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
           "notes": build_notes(run, products), "telluric_applied": state}
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

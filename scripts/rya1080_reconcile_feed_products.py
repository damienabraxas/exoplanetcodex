#!/usr/bin/env python3
"""
RYA-1080 — make the repo hold what the feed already publishes.
==============================================================
The RYA-1015 products were produced on the Mac, ingested into `data/products/solar/Fe.json`,
and never copied into the repo. `provenance.copied_to` — the field designed to record that
copy — is `None` on 67 of 75 live products. So every check that reads
`data/results/band_products/` sees a stale tree and answers a question about the feed with
an answer about the disk. That cost two wrong diagnoses in one session ("HARPS has no graded
product"; "Fe II has never been re-run"), the second of which nearly triggered a from-scratch
re-run of work already done.

🔴 AND `copied_to` BEING SET IS NOT EVIDENCE OF BEING COMMITTED. The eight rows that DO carry
a `copied_to` point into `/private/tmp/g3d/` and `/private/tmp/sirius_orphans/` — scratch
directories, outside the repo, untracked, gone at the next reboot. Counting them as
reconciled is how the exposure looked like 67 instead of 75. The guard therefore asks
"is it committed IN THIS REPO?", never "is the field non-null?".

This script does NOT re-measure and does NOT touch a single abundance, grade or sigma. It
copies the artifact the feed measured, verified by the checksum the feed recorded, and
writes back a REPO-RELATIVE `copied_to`.

    mismatch => BLOCKED, never committed. A file whose sha256 has changed is not the file
    the feed measured, whatever its numbers say, and committing it under the feed's hash
    would launder exactly the drift this ticket exists to close.

Usage:
    python3 scripts/rya1080_reconcile_feed_products.py            # dry run, reports only
    python3 scripts/rya1080_reconcile_feed_products.py --apply    # copy + rewrite copied_to
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FEED_DIR = ROOT / "data" / "products" / "solar"
DEST = ROOT / "data" / "results" / "band_products"
OUT = ROOT / "data" / "results" / "rya1080"

#: Published fields worth diffing when a checksum fails, so a mismatch is characterised
#: rather than merely reported. feed key -> products.csv column.
PUBLISHED = {"A": "A", "sigma_stat": "stat_dex", "sigma_syst": "syst_dex",
             "n_lines": "n_lines", "n_excluded": "n_excluded",
             "dominant_term": "dominant"}

#: `route` is normalised to upper case ON INGEST — all 57 checksum-VERIFYING products differ
#: from their file on it too. Comparing it would report a feed-side transform as file drift.
FEED_SIDE_TRANSFORMS = ("route",)


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _source_for(prov: dict) -> Path | None:
    """The artifact the feed measured. `path` is where it was produced; `copied_to` is
    where someone put it afterwards. Either can be the surviving copy, so try both."""
    for key in ("path", "copied_to"):
        v = prov.get(key)
        if v and Path(v).exists():
            return Path(v)
    return None


def _diff_published(product: dict, src: Path) -> list[str]:
    """What actually differs, for a product whose checksum failed."""
    try:
        df = pd.read_csv(src)
    except Exception as e:
        return [f"unreadable: {type(e).__name__}"]
    m = df[df.get("treatment").astype(str) == str(product.get("treatment"))] \
        if "treatment" in df.columns else df
    if not len(m):
        return ["no row matches this product's treatment"]
    r, out = m.iloc[0], []
    for fk, ck in PUBLISHED.items():
        if fk in FEED_SIDE_TRANSFORMS or fk not in product or ck not in df.columns:
            continue
        a, b = product[fk], r[ck]
        try:
            same = abs(float(a) - float(b)) < 1e-9
        except (TypeError, ValueError):
            same = str(a) == str(b)
        if not same:
            out.append(f"{fk}: feed={a!r} file={b!r}")
    return out


def reconcile(apply: bool) -> pd.DataFrame:
    rows = []
    for feed in sorted(FEED_DIR.glob("*.json")):
        doc = json.loads(feed.read_text())
        changed = False
        for p in doc.get("products", []):
            prov = p["provenance"]
            ident = {"feed": feed.name, "element": p.get("element"), "ion": p.get("ion"),
                     "band": p.get("band"), "instrument": p.get("instrument"),
                     "tier": p.get("tier"), "treatment": p.get("treatment")}
            src = _source_for(prov)
            expect = prov.get("sha256")

            if src is None:
                rows.append({**ident, "status": "UNREACHABLE", "detail":
                             f"neither path nor copied_to exists here: {prov.get('path')}"})
                continue
            if not expect:
                rows.append({**ident, "status": "NO_RECORDED_SHA", "detail": str(src)})
                continue

            actual = sha256(src)
            if actual != expect:
                rows.append({**ident, "status": "BLOCKED_SHA_MISMATCH",
                             "source": str(src), "feed_sha256": expect,
                             "actual_sha256": actual,
                             "detail": "; ".join(_diff_published(p, src))
                                       or "no published field differs — bytes only"})
                continue

            dest = DEST / src.name
            rel = str(dest.relative_to(ROOT))
            if apply:
                DEST.mkdir(parents=True, exist_ok=True)
                if not dest.exists() or sha256(dest) != expect:
                    shutil.copy2(src, dest)
                # Re-verify the COPY, not the source. A copy is a second chance to differ.
                if sha256(dest) != expect:
                    raise SystemExit(f"CRITICAL: copy of {src} does not match the feed sha256")
                if prov.get("copied_to") != rel:
                    prov["copied_to"] = rel
                    changed = True
            rows.append({**ident, "status": "RECONCILED" if apply else "READY",
                         "source": str(src), "copied_to": rel, "feed_sha256": expect})
        if apply and changed:
            feed.write_text(json.dumps(doc, indent=2) + "\n")
            print(f"[feed] rewrote copied_to in {feed.name}")
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true",
                    help="copy the artifacts in and rewrite copied_to (default: dry run)")
    args = ap.parse_args()

    df = reconcile(args.apply)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "rya1080_reconciliation.csv", index=False)

    print(f"\n=== RYA-1080 reconciliation ({'APPLY' if args.apply else 'DRY RUN'}) ===")
    for s, n in df.status.value_counts().items():
        print(f"  {s:<22} {n}")
    blocked = df[df.status == "BLOCKED_SHA_MISMATCH"]
    if len(blocked):
        blocked.to_csv(OUT / "rya1080_blocked.csv", index=False)
        print(f"\n=== 🔴 BLOCKED — not committed, by design ({len(blocked)}) ===")
        for _, r in blocked.iterrows():
            print(f"  {r.ion:<3}{str(r.band):<12}{str(r.tier):<11}{str(r.treatment):<10} "
                  f"{r.feed_sha256[:10]} -> {r.actual_sha256[:10]}  {r.detail}")
        print("\n  These files changed AFTER the feed ingested them. Committing one under "
              "the\n  feed's recorded sha256 would record a file the feed never measured.")
    print(f"\n[out] {OUT}/rya1080_reconciliation.csv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
RYA-1084 — re-ingest the ten RYA-1080-blocked products, once their writer is cleared.
=====================================================================================
RYA-1080 blocked ten feed products whose source files no longer matched the sha256 the
feed recorded. Re-ingesting BEFORE identifying the writer would launder an unverified
process into the feed as truth, so RYA-1080 held them. This script is the other half, and
it refuses to run on anything it has not been given a cleared reason for.

WHAT THE INVESTIGATION FOUND (see docs/science/rya1084_blocked_writer.md):

* The writer is `scripts/derive_band_products.py`, run legitimately. The ten rewrites come
  in five base + ENGINE-A pairs 20-40 s apart, which is the shape of a derivation run, not
  a stray touch.
* It is DETERMINISTIC. The 2026-08-27 13:34 run reproduced the file RYA-1051 had already
  committed BYTE FOR BYTE. A nondeterministic writer does not do that.
* Six of the ten differ by a new `deck` COLUMN, added by RYA-1044/1045 in the same window.
  Schema evolution, values untouched.
* 🔴 The ticket's leading hypothesis — a Python `round()` vs `np.round()` half-tie — is
  REFUTED. They agree at every point tested. The one moved value (`sigma_stat` 0.0217 ->
  0.0218) is a ONE-ULP difference in the upstream RMS amplified by a decimal boundary:
  0.02175 as a float is 0.0217499999999999985, just below the tie. `stat_basis` prints five
  decimals and says "0.02175" for both, which is why it looked like a tie.
* The write site is now `pipeline.error_budget.round_dex`, one rule, half-up on the
  shortest round-tripping decimal — stable across exactly that ULP pair.

So the current files are what a cleared, deterministic generator produced, and the FEED is
the stale side. Re-ingest updates the feed's recorded sha256 (and the one published bar
that genuinely moved) to match them.

Usage:
    python3 scripts/rya1084_reingest_blocked.py            # dry run
    python3 scripts/rya1084_reingest_blocked.py --apply
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FEED_DIR = ROOT / "data" / "products" / "solar"
BLOCKED = ROOT / "data" / "results" / "rya1080" / "rya1080_blocked.csv"
OUT = ROOT / "data" / "results" / "rya1084"

#: feed key -> products.csv column, for the fields the feed republishes.
PUBLISHED = {"A": "A", "sigma_stat": "stat_dex", "sigma_syst": "syst_dex",
             "n_lines": "n_lines", "n_excluded": "n_excluded",
             "dominant_term": "dominant"}

CLEARED = ("RYA-1084: writer identified as scripts/derive_band_products.py, shown "
           "deterministic (the 2026-08-27 13:34 run reproduced RYA-1051's committed bytes "
           "exactly); diffs are the RYA-1044/1045 `deck` column plus one 1-ULP rounding "
           "boundary now pinned by pipeline.error_budget.round_dex")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    blocked = pd.read_csv(BLOCKED)
    by_feed: dict[str, pd.DataFrame] = {f: g for f, g in blocked.groupby("feed")}
    rows = []

    for feed_name, g in by_feed.items():
        feed = FEED_DIR / feed_name
        doc = json.loads(feed.read_text())
        for _, b in g.iterrows():
            src = Path(b.source)
            if not src.exists():
                rows.append({"status": "SOURCE_GONE", "source": str(src)})
                continue
            actual = hashlib.sha256(src.read_bytes()).hexdigest()
            if actual != b.actual_sha256:
                # The file moved again since RYA-1080 measured it. Re-ingesting now would
                # record a third state nobody has looked at.
                rows.append({"status": "MOVED_AGAIN", "source": str(src),
                             "expected": b.actual_sha256, "found": actual})
                continue

            match = [p for p in doc["products"]
                     if p["provenance"].get("sha256") == b.feed_sha256
                     and str(p.get("treatment")) == str(b.treatment)
                     and str(p.get("band")) == str(b.band)
                     and str(p.get("ion")) == str(b.ion)]
            if len(match) != 1:
                rows.append({"status": "AMBIGUOUS_FEED_ROW", "source": str(src),
                             "n_matched": len(match)})
                continue
            p = match[0]

            df = pd.read_csv(src)
            m = df[df.treatment.astype(str) == str(p["treatment"])]
            if len(m) != 1:
                rows.append({"status": "NO_UNIQUE_TREATMENT_ROW", "source": str(src)})
                continue
            r = m.iloc[0]

            moved = {}
            for fk, ck in PUBLISHED.items():
                if fk not in p or ck not in df.columns:
                    continue
                a, c = p[fk], r[ck]
                try:
                    same = abs(float(a) - float(c)) < 1e-12
                    newv = float(c)
                except (TypeError, ValueError):
                    same, newv = str(a) == str(c), c
                if not same:
                    moved[fk] = (a, newv)

            if args.apply:
                for fk, (_, newv) in moved.items():
                    p[fk] = newv
                p["provenance"]["sha256"] = actual
                p["provenance"]["reingested_by"] = "RYA-1084"
                p["provenance"]["reingest_reason"] = CLEARED
            rows.append({"status": "REINGESTED" if args.apply else "READY",
                         "ion": b.ion, "band": b.band, "treatment": b.treatment,
                         "source": str(src), "old_sha256": b.feed_sha256,
                         "new_sha256": actual,
                         "fields_moved": "; ".join(f"{k}: {v[0]} -> {v[1]}"
                                                   for k, v in moved.items()) or "none"})
        if args.apply:
            feed.write_text(json.dumps(doc, indent=2) + "\n")
            print(f"[feed] re-ingested into {feed_name}")

    out = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / "rya1084_reingest.csv", index=False)
    print(f"\n=== RYA-1084 re-ingest ({'APPLY' if args.apply else 'DRY RUN'}) ===")
    for s, n in out.status.value_counts().items():
        print(f"  {s:<24} {n}")
    moved_rows = out[out.get("fields_moved", "none") != "none"] if len(out) else out
    print(f"\n  published fields that actually moved: {len(moved_rows)} of {len(out)}")
    for _, r in moved_rows.iterrows():
        print(f"    {r.ion} {r.band} {r.treatment}: {r.fields_moved}")
    print(f"\n[out] {OUT}/rya1084_reingest.csv")


if __name__ == "__main__":
    main()

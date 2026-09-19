#!/usr/bin/env python3
"""RYA-1223 Step 1 — VERIFY. Diff the LIVE site against the current feed. Writes no page.

    python3 scripts/rya1223_site_reconcile.py --site <exoplanetcodex-site checkout> \
        [--output data/audit/rya1223_site_reconcile]

Ryan's 2026-09-17 rescope: this is a verify-first INCREMENTAL reconcile, not a rebuild.
"Catching a wrong live value is the priority." So every product currently rendered on the
live site is bucketed against the feed BEFORE anything is edited:

    CORRECT  live A and sigma_reported match the feed        -> leave untouched
    STALE    live value or uncertainty differs               -> needs correction
    NEW      in the feed, not yet on the site                -> needs adding
    ORPHAN   on the site, not in the feed                    -> FLAG, never delete

🔴 SCOPE IS FOUR BANDS (Ryan's 2026-09-19 Option B): VIS / red-optical / NIR / H. near-UV is
REMOVED from this ticket and publishes after RYA-1226 wires its measured components; it is
reported here under its own bucket as HELD so the split is visible rather than silent, and
so nobody later reads its absence as an oversight.

⚠️ THE SITE'S OWN SNAPSHOT IS THE LIVE TRUTH, NOT A LOCAL CLONE'S BRANCH. The RYA-1221 [J]
finding ("site HEAD 61cfb27, PRs never deployed") was a stale-clone read and was WRONG --
the site is live at #36. So the comparison reads `assets/data/fe-publication/Fe.json` from a
checkout of the site's REMOTE main, and records the SHA it read, so the claim can be
re-checked rather than believed.

⚠️ IDENTITY IS THE RYA-1127 KEY, NOT (ion, band, holding, treatment). The latter is not
unique over this feed -- GRADED and DEEPGRADED products share it -- so a dict on it
silently collapses pairs and then reports every loser as a difference.
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FEED = ROOT / "data/products/solar/Fe.json"
SITE_FEED = "assets/data/fe-publication/Fe.json"

#: Ryan's Option B scope. near-UV is tracked separately, never folded in.
BANDS = ("VIS", "red-optical", "NIR", "H")
HELD_BAND = "near-UV"

#: The fields whose disagreement makes a live row STALE. `A` is the value; the other two are
#: the bar. A bar that is too SMALL is the defect RYA-1213/587 closed, so the uncertainty is
#: compared as strictly as the abundance.
COMPARED = ("A", "sigma_reported", "sigma_xi", "xi_state")

KEY_FIELDS = ("element", "ion", "band", "instrument", "holding", "tier", "selector",
              "route", "treatment")


def key(p: dict) -> tuple:
    return tuple(str(p.get(k) or "") for k in KEY_FIELDS)


def label(k: tuple) -> str:
    el, ion, band, inst, hold, tier, sel, route, treat = k
    return f"{el} {ion} | {band} | {tier} | {inst} | {hold} | {treat}"


def _git(cwd: Path, *args: str) -> str:
    return subprocess.run(["git", *args], cwd=cwd, capture_output=True,
                          text=True).stdout.strip()


def reconcile(feed: dict, site: dict) -> dict:
    live = {key(p): p for p in site["products"]}
    cur = {key(p): p for p in feed["products"]}
    if len(live) != len(site["products"]):
        raise SystemExit("the site snapshot's identity key is not unique -- refusing to "
                         "bucket, because a collapsed pair reads as a difference")
    if len(cur) != len(feed["products"]):
        raise SystemExit("the feed's identity key is not unique -- refusing to bucket")

    out = {b: {"CORRECT": [], "STALE": [], "NEW": [], "ORPHAN": []}
           for b in (*BANDS, HELD_BAND)}
    other: dict = {"CORRECT": [], "STALE": [], "NEW": [], "ORPHAN": []}

    def bucket(band: str) -> dict:
        return out.get(band, other)

    for k, now in cur.items():
        band = now["band"]
        was = live.get(k)
        if was is None:
            bucket(band)["NEW"].append({
                "product": label(k),
                **{f"{f}": now.get(f) for f in COMPARED},
            })
            continue
        diffs = {f: {"live": was.get(f), "feed": now.get(f)}
                 for f in COMPARED if was.get(f) != now.get(f)}
        if diffs:
            bucket(band)["STALE"].append({"product": label(k), "differs": diffs})
        else:
            bucket(band)["CORRECT"].append({"product": label(k)})

    for k, was in live.items():
        if k not in cur:
            bucket(was["band"])["ORPHAN"].append({
                "product": label(k),
                **{f"{f}": was.get(f) for f in COMPARED},
            })

    counts = {b: {s: len(v) for s, v in d.items()} for b, d in out.items()}
    counts["_other_bands"] = {s: len(v) for s, v in other.items()}
    in_scope = {s: sum(counts[b][s] for b in BANDS) for s in
                ("CORRECT", "STALE", "NEW", "ORPHAN")}
    return {"counts_by_band": counts, "in_scope_totals": in_scope,
            "by_band": out, "other_bands": other}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--site", required=True, type=Path,
                    help="a checkout of exoplanetcodex-site at its REMOTE main")
    ap.add_argument("--output", type=Path,
                    default=ROOT / "data/audit/rya1223_site_reconcile")
    a = ap.parse_args()

    feed = json.loads(FEED.read_text())
    site_path = a.site / SITE_FEED
    if not site_path.is_file():
        raise SystemExit(f"no site feed snapshot at {site_path}")
    site = json.loads(site_path.read_text())

    rep = reconcile(feed, site)
    rep["ticket"] = "RYA-1223"
    rep["read_only"] = True
    rep["scope"] = {"bands": list(BANDS), "held": HELD_BAND,
                    "why": ("Ryan 2026-09-19 Option B: publish the complete bands now, "
                            "hold near-UV until RYA-1226 wires its measured components.")}
    rep["feed"] = {"path": str(FEED.relative_to(ROOT)), "version": feed.get("version"),
                   "products": len(feed["products"]),
                   "codex_commit": _git(ROOT, "rev-parse", "HEAD")}
    rep["site"] = {"checkout": str(a.site), "snapshot": SITE_FEED,
                   "feed_version_rendered": site.get("version"),
                   "products_rendered": len(site["products"]),
                   "site_commit": _git(a.site, "rev-parse", "HEAD"),
                   "site_describe": _git(a.site, "log", "--oneline", "-1")}
    rep["compared_fields"] = list(COMPARED)

    print(f"feed  : v{rep['feed']['version']}  {rep['feed']['products']} products  "
          f"(codex {rep['feed']['codex_commit'][:8]})")
    print(f"site  : v{rep['site']['feed_version_rendered']}  "
          f"{rep['site']['products_rendered']} rendered  "
          f"(site {rep['site']['site_commit'][:8]})")
    print()
    print(f"{'band':14}{'CORRECT':>9}{'STALE':>7}{'NEW':>6}{'ORPHAN':>8}")
    for b in (*BANDS, HELD_BAND):
        c = rep["counts_by_band"][b]
        tag = "   <- HELD (not this ticket)" if b == HELD_BAND else ""
        print(f"{b:14}{c['CORRECT']:>9}{c['STALE']:>7}{c['NEW']:>6}{c['ORPHAN']:>8}{tag}")
    t = rep["in_scope_totals"]
    print(f"{'IN SCOPE (4)':14}{t['CORRECT']:>9}{t['STALE']:>7}{t['NEW']:>6}{t['ORPHAN']:>8}")

    a.output.mkdir(parents=True, exist_ok=True)
    (a.output / "reconciliation.json").write_text(json.dumps(rep, indent=2) + "\n")
    print(f"\nwrote {(a.output / 'reconciliation.json').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

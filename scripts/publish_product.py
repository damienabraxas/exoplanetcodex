#!/usr/bin/env python3
"""THE source of truth for published element products — RYA-1034.

WHY THIS EXISTS. A real, current, four-arm 3D-NLTE result (A(Fe I) VIS = 7.511, n=50)
was produced on Sirius and then could not be found: it lived in `~/out_g3d_9e651a/`,
an ad-hoc directory, on the only machine that can run that leg. Meanwhile every Mac
copy of `rya817_3dnlte_products.csv` still held the superseded ungraded 7.604, and
`data/results/band_products/` is GITIGNORED, so nothing about any of it was under
version control. The value existed in a Linear ticket and nowhere else durable.

That is not a one-off. It is the third time this week a result outlived its artifact.

THE RULE THIS FILE ENFORCES: a product is published HERE or it does not exist.

  data/products/solar/<El>.json      <- git-TRACKED, one file per element

Tracked, not ignored: the whole failure above is what a gitignored product directory
buys you. These files are small (values and provenance, never per-line data), so they
cost nothing to version and they diff readably.

FOUR PROPERTIES, each of which was violated at least once this week:

  1. ONE KEY, ONE CURRENT VALUE. A product is identified by
     (element, ion, band, instrument, holding, tier, treatment). Publishing the same key
     twice UPDATES it. No glob decides which of five copies wins.
  2. NOTHING IS EVER DELETED. A changed value moves to `superseded` with the reason and
     the timestamp (RYA-711: quarantined, never culled). You can always answer "what did
     this say yesterday, and why did it change".
  3. A CHANGE NEEDS A STATED REASON. Overwriting a differing value without `--reason`
     is refused. "We only change it when there's a damn good reason" is a check, not an
     aspiration.
  4. PROVENANCE TRAVELS WITH THE NUMBER. host, absolute source path, sha256, and the
     artifact's own mtime. Copying a foreign artifact without this is laundering
     (RYA-772); with it, the row states exactly where it came from and can be re-verified.

Idempotent: re-publishing byte-identical content is a no-op and does NOT bump the
version, so a re-run in a loop cannot inflate the history.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import socket
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
STORE = ROOT / "data" / "products"
SCHEMA = "codex.element_product/1"

#: RYA-906 axes -> the display name. Derived, never stored by hand: the whole point of
#: RYA-906 is that the NAME follows the physics, so a renamed engine cannot strand a label.
_DISPLAY = {
    "1D-LTE":            "Synth · 1D-LTE",
    "ENGINE-A":          "EW · 1D-NLTE · Bergemann",
    "ENGINE-B":          "Synth · 1D-LTE",
    "ENGINE-B-NLTE":     "Synth · 1D-NLTE · Gerber",
    "ENGINE-A-3DNLTE":   "EW · 3D-NLTE · Amarsi",
}


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def display_name(treatment: str, *, gf: str | None, route: str | None) -> str:
    """`treatment` names the physics; the tier suffix names the gf pedigree (RYA-906)."""
    base = _DISPLAY.get(str(treatment), str(treatment))
    # The band route writes 1D-LTE for BOTH the EW and the synth leg; only the route
    # separates them, and reading a profile-fit row as "Synth" was a real mislabel.
    if str(treatment) == "1D-LTE" and str(route or "").upper().startswith("PROFILEFIT"):
        base = "EW · 1D-LTE"
    if str(gf or "").lower() in {"lab", "cited lab", "gf scale (cited lab)"}:
        base += " · lab-gf"
    return base


def normalise(df: pd.DataFrame, *, holding: str, tier: str, route: str | None,
              selector: str | None = None) -> list[dict]:
    """Both product schemas -> one row shape.

    🔴 THE TWO ROUTES DISAGREE ON COLUMN NAMES, and that is part of why 3D never joined
    the same picture: the band route writes `A`/`stat_dex`/`syst_dex`, the RYA-817 3D
    route writes `value`/`sigma`. Normalising here -- rather than teaching every reader
    both dialects -- is the single-sourcing rule (RYA-353) applied to a schema.
    """
    out = []
    for _, r in df.iterrows():
        has_A = "A" in df.columns and pd.notna(r.get("A"))
        A = r.get("A") if has_A else r.get("value")
        if pd.isna(A):
            continue                      # a cell with no value is a GAP, not a product
        stat = r.get("stat_dex") if "stat_dex" in df.columns else r.get("sigma")
        syst = r.get("syst_dex") if "syst_dex" in df.columns else None
        treatment = str(r.get("treatment"))
        out.append({
            "element": str(r.get("element")), "ion": str(r.get("ion")),
            "band": str(r.get("band")), "instrument": str(r.get("instrument")),
            "holding": holding, "tier": tier, "selector": selector,
            "treatment": treatment,
            "display": display_name(treatment, gf=r.get("gf") or r.get("dominant"),
                                    route=route),
            "A": round(float(A), 4),
            "sigma_stat": None if stat is None or pd.isna(stat) else round(float(stat), 4),
            "sigma_syst": None if syst is None or pd.isna(syst) else round(float(syst), 4),
            "n_lines": None if pd.isna(r.get("n_lines")) else int(r.get("n_lines")),
            "n_excluded": None if pd.isna(r.get("n_excluded")) else int(r.get("n_excluded")),
            "dominant_term": (None if pd.isna(r.get("dominant")) else str(r.get("dominant")))
                             if "dominant" in df.columns else None,
            "route": route,
        })
    return out


#: 🔴 ROUTE AND SELECTOR ARE PART OF PRODUCT IDENTITY, not metadata. Leaving them out
#: collided the EW and synthesis legs onto one key: IAG VIS GRADED ENGINE-A is 7.481 on
#: n=4 by profile fit and 7.484 on n=37 by synthesis -- two different measurements of two
#: different pools, and the store would have silently kept whichever was written last.
#: RYA-1026 lists `1D-LTE` and `1D-LTE Synth` as SEPARATE rows for exactly this reason,
#: and RYA-984 says two runs differing only in selector are two products. Caught by the
#: overwrite guard on the very first backfill, which is what the guard is for.
KEY_FIELDS = ("element", "ion", "band", "instrument", "holding",
              "tier", "selector", "route", "treatment")


def key_of(row: dict) -> str:
    return "|".join(str(row.get(k) or "") for k in KEY_FIELDS)


def load(element: str, star: str) -> dict:
    p = STORE / star / f"{element}.json"
    if p.exists():
        return json.loads(p.read_text())
    return {"schema": SCHEMA, "element": element, "star": star, "version": "1.0",
            "updated_at": None,
            #: CURRENT, publishable. What the tracker and the website read.
            "products": [],
            #: ARCHIVE — prior values of a key that later changed. Answers "what did this
            #: say yesterday, and why did it move".
            "archive": [],
            #: QUARANTINE — products WITHDRAWN from display but kept in full (RYA-711:
            #: quarantined, never culled). A withdrawn product still has to be defensible
            #: in the appendix, so deleting it destroys the evidence for its own rejection.
            "quarantine": []}


def bump(v: str) -> str:
    major, _, minor = v.partition(".")
    return f"{major}.{int(minor or 0) + 1}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="src", nargs="+",
                    help="one or more *_products.csv to publish. Several files in ONE "
                         "invocation is ONE version bump: a backfill of 23 artifacts is "
                         "a single change to the store, not 23 of them.")
    ap.add_argument("--holding",
                    help="the holding the spectrum came from. REQUIRED and never inferred "
                         "from the filename: RYA-1026 makes two holdings of one instrument "
                         "two different PRODUCTS, and a name containing 'corrected' is not "
                         "evidence that it is.")
    ap.add_argument("--tier",
                    choices=["GRADED", "DEEPGRADED", "CONSISTENT", "UNGRADED"],
                    help="the line-selection tier — Ryan's three sub-products plus ungraded")
    ap.add_argument("--route", default=None, help="PROFILEFIT | SYNTH | EW-3D")
    ap.add_argument("--selector", default=None,
                    help="the line-selection variant as the stem records it (GRADED, "
                         "DEEPGRADED, FROMEW, ...). Part of the KEY: RYA-984 makes two "
                         "runs differing only in selector two different products.")
    ap.add_argument("--star", default="solar")
    ap.add_argument("--host", default=None, help="where the artifact was PRODUCED "
                                                 "(default: this machine's hostname)")
    ap.add_argument("--origin-path", default=None,
                    help="the artifact's path ON THE HOST THAT PRODUCED IT. Required in "
                         "spirit whenever --host is not this machine: recording the path "
                         "of a local COPY would state that a Sirius-only result was "
                         "produced here, which is the laundering RYA-772 names. The "
                         "sha256 is of the bytes actually read, so the copy stays "
                         "verifiable against the origin.")
    ap.add_argument("--reason", default=None,
                    help="why an existing differing value is being replaced. REQUIRED to "
                         "overwrite; a value that changes without a stated reason is the "
                         "thing this store exists to prevent.")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--quarantine-where", default=None, metavar="FIELD=VALUE",
                    help="withdraw every CURRENT product whose FIELD equals VALUE "
                         "(e.g. holding=solar_kpno). Requires --reason. The rows move to "
                         "`quarantine` in full — never deleted — so the appendix can still "
                         "defend why they were withdrawn (RYA-711/844).")
    ap.add_argument("--element", default=None, help="element to act on for --quarantine-where")
    ap.add_argument("--quarantine-older-than", default=None, metavar="ISO8601",
                    help="withdraw every CURRENT product whose ARTIFACT was produced "
                         "before this instant. Keyed on `provenance.artifact_mtime` -- "
                         "when the measurement was made -- NOT on ingested_at, which only "
                         "says when it was copied into the store and would mark a "
                         "freshly-backfilled stale result as current.")
    a = ap.parse_args()

    if a.quarantine_where or a.quarantine_older_than:
        if not a.reason:
            print("REFUSING: --quarantine-where needs --reason. A product withdrawn "
                  "without a stated reason cannot be defended in the appendix, and is "
                  "indistinguishable from one that was lost.", file=sys.stderr)
            return 5
        if not a.element:
            print("REFUSING: --quarantine-where needs --element.", file=sys.stderr)
            return 5
        field = value = None
        if a.quarantine_where:
            field, _, value = a.quarantine_where.partition("=")
        doc = load(a.element, a.star)
        doc.setdefault("archive", []); doc.setdefault("quarantine", [])
        keep, pulled = [], []
        for row in doc["products"]:
            hit = (field is not None and str(row.get(field)) == value)
            if a.quarantine_older_than:
                mt = (row.get("provenance") or {}).get("artifact_mtime") or ""
                hit = hit or (mt < a.quarantine_older_than)
            if hit:
                r = dict(row)
                r["quarantined_at"] = _now(); r["quarantine_reason"] = a.reason
                doc["quarantine"].append(r); pulled.append(key_of(row))
            else:
                keep.append(row)
        if not pulled:
            print("no current product matches "
                  + (f"{field}={value} " if field else "")
                  + (f"older-than {a.quarantine_older_than}" if a.quarantine_older_than else ""))
            return 0
        doc["products"] = keep
        doc["version"] = bump(doc["version"]); doc["updated_at"] = _now()
        out = STORE / a.star / f"{a.element}.json"
        if a.dry_run:
            print(f"[dry-run] would quarantine {len(pulled)} row(s)")
            return 0
        out.write_text(json.dumps(doc, indent=2) + "\n")
        print(f"{out.relative_to(ROOT)}  ->  v{doc['version']}   "
              f"QUARANTINED {len(pulled)}, {len(doc['products'])} remain current")
        for k in pulled: print(f"    ! {k}")
        print(f"    reason: {a.reason}")
        return 0

    if not (a.src and a.holding and a.tier):
        print("publishing needs --from, --holding and --tier "
              "(or use --quarantine-where)", file=sys.stderr)
        return 2
    srcs = [Path(x).resolve() for x in a.src]
    for src in srcs:
        if not src.exists():
            print(f"no such product file: {src}", file=sys.stderr)
            return 2
    pending = []
    for src in srcs:
        df = pd.read_csv(src)
        rows = normalise(df, holding=a.holding, tier=a.tier, route=a.route,
                         selector=a.selector)
        if not rows:
            print(f"REFUSING: {src.name} carries no row with a value. An empty publish "
                  f"would read as a measurement of nothing (RYA-833) — if the cell is a "
                  f"gap, it belongs in the gap list, not here.", file=sys.stderr)
            return 3
        prov = {"host": a.host or socket.gethostname().split(".")[0],
                "path": a.origin_path or str(src),
                "copied_to": str(src) if a.origin_path else None,
                "sha256": _sha256(src),
                "artifact_mtime": _dt.datetime.fromtimestamp(
                    src.stat().st_mtime, _dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "ingested_at": _now()}
        for r in rows:
            r["provenance"] = prov
        pending.extend(rows)

    element = pending[0]["element"]
    doc = load(element, a.star)
    doc.setdefault("archive", []); doc.setdefault("quarantine", [])

    by_key = {key_of(r): i for i, r in enumerate(doc["products"])}
    added, updated, unchanged = [], [], []
    for row in pending:
        k = key_of(row)
        if k not in by_key:
            doc["products"].append(row); added.append(k); continue
        cur = doc["products"][by_key[k]]
        same = all(cur.get(f) == row.get(f) for f in ("A", "sigma_stat", "sigma_syst",
                                                      "n_lines", "n_excluded"))
        if same:
            unchanged.append(k); continue
        if not a.reason:
            print(f"REFUSING to change {k}\n"
                  f"    current {cur['A']} (n={cur['n_lines']}) -> new {row['A']} "
                  f"(n={row['n_lines']})\n"
                  f"    Pass --reason. A published value that moves without a stated "
                  f"reason is exactly what this store exists to prevent.", file=sys.stderr)
            return 4
        old = dict(cur)
        old["superseded_at"] = _now(); old["superseded_reason"] = a.reason
        doc.setdefault("archive", []).append(old)        # RYA-711: never deleted
        doc["products"][by_key[k]] = row
        updated.append(k)

    if not added and not updated:
        print(f"no change — {len(unchanged)} row(s) already current at v{doc['version']}")
        return 0

    doc["version"] = bump(doc["version"])
    doc["updated_at"] = _now()
    doc["products"].sort(key=key_of)
    out = STORE / a.star / f"{element}.json"
    if a.dry_run:
        print(f"[dry-run] would write {out} at v{doc['version']}: "
              f"+{len(added)} ~{len(updated)} ={len(unchanged)}")
        return 0
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"{out.relative_to(ROOT)}  ->  v{doc['version']}   "
          f"added {len(added)}, updated {len(updated)}, unchanged {len(unchanged)}, "
          f"total {len(doc['products'])}")
    for k in added:   print(f"    + {k}")
    for k in updated: print(f"    ~ {k}   ({a.reason})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

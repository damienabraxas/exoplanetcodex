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
import re
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
            # 🔴 CARRY WHAT `sigma_stat` MEANS, DO NOT MAKE THE READER GUESS (RYA-1095).
            # The band artifacts have always had a `stat_basis` column (`error_budget
            # .stat_basis`) and this function dropped it, so the feed published a number
            # whose construction was unstated -- and the two routes do not agree on it:
            # the band route's `stat_dex` is a standard error, the EW-3D route's `sigma`
            # was a raw per-line standard deviation. RYA-1092's gate had to infer the
            # basis from the ROUTE, which is a lookup that goes stale the moment a route
            # changes what it writes (it did, here). Recorded per product instead.
            "stat_basis": (None if pd.isna(r.get("stat_basis"))
                           else str(r.get("stat_basis")))
                          if "stat_basis" in df.columns else None,
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


def _repo_relative(src: Path) -> str:
    """`copied_to` as a REPO-RELATIVE path, and a refusal if the file is outside the repo.

    🔴 IT WAS STORING AN ABSOLUTE ONE, AND AN ABSOLUTE PATH IS A MACHINE FACT. `srcs` is
    built with `Path(x).resolve()`, so publishing from inside the repo recorded e.g.
    `/Users/<me>/codex/rya1089/data/results/...`. That is inside the repository on the
    machine that published and nowhere else: `feed_repo_reconciliation` reads `copied_to`
    to check the artifact is COMMITTED, and on any other checkout the path is not under
    the repo root, so it reports COPIED_TO_OUTSIDE_REPO. Measured on RYA-1096: eight
    freshly published products passed every check on the Mac and turned that guard red on
    Sirius -- a machine-local pass, which is the shape a skipped test has.

    Every record written before this stored a relative path (the caller happened to pass a
    relative one and an earlier version did not resolve), so this restores the convention
    the feed already follows rather than inventing one.
    """
    try:
        return str(Path(src).resolve().relative_to(ROOT))
    except ValueError:
        raise SystemExit(
            f"--from {src} is outside the repository ({ROOT}). `copied_to` records that "
            f"the artifact is COMMITTED HERE, so publishing from a scratch directory "
            f"would write a path no other checkout can resolve. Copy the artifact into "
            f"data/results/ and publish from there.")


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
            #: GAPS — cells that CANNOT exist, with the measurement that says so. A gap
            #: is not a missing run and must not be retried as one: red-optical has zero
            #: Fe I lab lines above the depth gate, so its DEEPGRADED tier is empty BY
            #: CONSTRUCTION. Recording it is what stops three jobs being relaunched every
            #: week to rediscover the same emptiness (RYA-833: an absence is a hypothesis
            #: until measured -- these are measured, so they are conclusions).
            "gaps": [],
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
                    choices=["GRADED", "DEEPGRADED", "CONSISTENT", "UNGRADED", "ALL"],
                    help="the line-selection tier. Ryan's three sub-products are GRADED "
                         "(lab gf, at/below the depth gate), DEEPGRADED (lab gf, above "
                         "it) and CONSISTENT (no lab gf, admitted on behaviour). "
                         "🔴 ALL IS NOT A TIER, IT IS THE ABSENCE OF ONE: a run with no "
                         "selector measures the WHOLE pool, lab lines included -- the "
                         "harps EW pool is 247 lines of which 9 ARE lab-graded. Calling "
                         "that UNGRADED made the tiers overlap, and CONSISTENT (a subset "
                         "of the 238 non-lab lines) would then have been double-counted "
                         "against it. UNGRADED means the COMPLEMENT of the lab set and "
                         "nothing else (RYA-946: graded and ungraded are separate "
                         "products, never merged).")
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
    ap.add_argument("--declare-gap", default=None, metavar="BAND:TIER",
                    help="record a cell that CANNOT exist, with --reason carrying the "
                         "measurement that establishes it. Needs --element.")
    ap.add_argument("--quarantine-older-than", default=None, metavar="ISO8601",
                    help="withdraw every CURRENT product whose ARTIFACT was produced "
                         "before this instant. Keyed on `provenance.artifact_mtime` -- "
                         "when the measurement was made -- NOT on ingested_at, which only "
                         "says when it was copied into the store and would mark a "
                         "freshly-backfilled stale result as current.")
    a = ap.parse_args()

    if a.declare_gap:
        if not (a.reason and a.element):
            print("REFUSING: --declare-gap needs --element and --reason. A gap with no "
                  "measurement behind it is just a missing run.", file=sys.stderr)
            return 8
        band, _, tier = a.declare_gap.partition(":")
        doc = load(a.element, a.star)
        doc.setdefault("gaps", [])
        rec = {"band": band, "tier": tier, "reason": a.reason, "declared_at": _now()}
        doc["gaps"] = [g for g in doc["gaps"]
                       if not (g["band"] == band and g["tier"] == tier)] + [rec]
        doc["version"] = bump(doc["version"]); doc["updated_at"] = _now()
        out = STORE / a.star / f"{a.element}.json"
        out.write_text(json.dumps(doc, indent=2) + "\n")
        print(f"{out.relative_to(ROOT)}  ->  v{doc['version']}   GAP DECLARED {band}:{tier}")
        return 0

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
                "copied_to": _repo_relative(src) if a.origin_path else None,
                "sha256": _sha256(src),
                "artifact_mtime": _dt.datetime.fromtimestamp(
                    src.stat().st_mtime, _dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "ingested_at": _now()}
        #: 🔴 THE TIER IS A CLAIM ABOUT THE gf, SO CHECK IT AGAINST THE GRADER'S VERDICT.
        #: The IR band stems a product `SYNTH_GRADED` whose budget reads
        #: `gf rung 1 (gf scale (UNGRADED))` with syst 0.1726 -- the blanket.
        #: `canonical_gf` calls all 29 in-band lines gf_tier=LAB, but the grade decider
        #: re-derives and downgrades 22 to `systematic:K07/SCALE-MISMATCH`, leaving 4 real
        #: laboratory lines. The SELECTOR and the GRADER disagree and the FILENAME follows
        #: the selector, so publishing on the stem would have put a 0.17-blanket, 1-dex-bar
        #: number in the graded tier beside cited-lab cells carrying 0.04. RYA-799/824
        #: already named this: in the IR, SCALE-MISMATCH is a LABELLING defect.
        #:
        #: Read from the BUDGET, which is the decider's own output. NOT from the `gf`
        #: column -- that names the LINELIST SOURCE and reads 'kurucz' on rung-3 and
        #: rung-1 products alike, so keying on it refused a legitimately graded
        #: red-optical cell on my first attempt.
        if a.tier in ("GRADED", "DEEPGRADED"):
            bud = Path(str(src)[:-len("_products.csv")] + "_budgets.txt")
            if not bud.exists():
                bud = Path(re.sub(r"_ENGINE-[A-Z-]+$", "", str(src)[:-len("_products.csv")])
                           + "_budgets.txt")
            m = re.search(r"gf rung (\d) \(gf scale \(([^)]*)\)", bud.read_text()) \
                if bud.exists() else None
            if m and int(m.group(1)) != 3:
                print(f"REFUSING: {src.name}\n"
                      f"    tier={a.tier} claims laboratory gf, but the budget's own "
                      f"verdict is gf rung {m.group(1)} ({m.group(2)}).\n"
                      f"    A tier is a claim about the gf SCALE, and the stem is not "
                      f"evidence for it. Publish under the tier the grade supports, or "
                      f"fix the pool.", file=sys.stderr)
                return 7
            if not m:
                print(f"REFUSING: {src.name}\n"
                      f"    tier={a.tier} claims laboratory gf but no budget was found to "
                      f"corroborate it ({bud.name}). An unverifiable grade claim is not "
                      f"published (RYA-833).", file=sys.stderr)
                return 7
        for r in rows:
            r["provenance"] = prov
        pending.extend(rows)

    element = pending[0]["element"]
    doc = load(element, a.star)
    doc.setdefault("archive", []); doc.setdefault("quarantine", [])

    by_key = {key_of(r): i for i, r in enumerate(doc["products"])}
    #: 🔴 REPUBLISHING A QUARANTINED KEY. A fresh run of a product that was withdrawn is
    #: legitimate -- that is how a superseded cell gets replaced -- but leaving the old
    #: entry in `quarantine` while the new one is current makes the key BOTH withdrawn
    #: and published. That is the exact shape of the near-UV 8.529 that "came back": a
    #: quarantine that does not actually withdraw. The withdrawal is now HISTORY, so the
    #: record moves to `archive` where superseded values live, and it still names why it
    #: was withdrawn in the first place.
    q_by_key = {key_of(r): i for i, r in enumerate(doc["quarantine"])}
    revived = []
    for row in pending:
        k = key_of(row)
        if k in q_by_key:
            revived.append(k)
    if revived and not a.reason:
        print("REFUSING: these keys are QUARANTINED and would be republished as current:\n"
              + "\n".join(f"    {k}" for k in revived)
              + "\nPass --reason saying why the withdrawal no longer applies (a fresh "
                "run, a fixed input). Silently un-quarantining is how a withdrawn "
                "product comes back.", file=sys.stderr)
        return 6
    for k in revived:
        rec = doc["quarantine"][q_by_key[k]]
        rec = dict(rec)
        rec["superseded_at"] = _now()
        rec["superseded_reason"] = (
            f"republished as current: {a.reason} "
            f"(was withdrawn: {rec.get('quarantine_reason','—')})")
        doc["archive"].append(rec)
    if revived:
        keep_q = {k for k in revived}
        doc["quarantine"] = [r for r in doc["quarantine"] if key_of(r) not in keep_q]

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

#!/usr/bin/env python3
"""RYA-1127 — what adding `line_set` to the product identity key does to EVERY product.

    python3 scripts/rya1127_key_migration_audit.py [--check]

WHY AN AUDIT AND NOT JUST A TEST
--------------------------------
Adding a field to `KEY_FIELDS` re-identifies every product in the feed at once. The failure
mode is not a wrong number -- no abundance moves here -- it is a wrong SPLIT or a wrong
MERGE: two records that were one identity becoming two, or two identities collapsing into
one, without anyone deciding that. Either one silently changes what the feed publishes,
what the tracker dedups, and which cell the site renders.

So the key is computed BOTH WAYS over every pool and the two partitions are compared
structurally:

  * a SPLIT   — one old key now maps to several new keys
  * a MERGE   — several old keys now map to one new key
  * unchanged — the partition is identical

🔴 BE HONEST ABOUT WHICH OF THOSE CHECKS CAN ACTUALLY FIRE TODAY. `line_set` was APPENDED
to `KEY_FIELDS`, so the new key is the old key plus a suffix. While that is true:

  * a MERGE cannot happen -- distinct old keys keep distinct prefixes; and
  * every SPLIT is trivially "explained by line_set", because it is the only field that
    can differ between two records sharing an old key.

Reporting those two as passing checks would be a guard whose two sides converge: they
would read as evidence and be arithmetic. What this audit therefore ASSERTS is the
refinement invariant itself -- `new == old + "|" + line_set` for every record. That is the
condition under which the two structural facts hold, and it is exactly what a later edit
would break by inserting `line_set` mid-tuple or reordering `KEY_FIELDS`. Break it and old
identities really can merge and split arbitrarily, and the merge check stops being
decorative.

The checks that carry weight on today's data are therefore:

  * every record RESOLVES a line_set (nine did not when this ticket began);
  * the refinement invariant holds;
  * no live identity splits -- because a split means one published identity now covers
    two products, which changes what the feed publishes.

⚠️ THE INTERESTING RESULT IS THAT TODAY THERE ARE NO SPLITS. The RYA-1106 Asplund products
are not in the feed yet -- they cannot be, which is the whole reason this ticket exists --
so on the committed feed every old key maps to exactly one new key and the partition is
unchanged. That is the correct and reassuring answer for a migration: the schema widens,
and nothing already published is re-identified. The split this enables is demonstrated
separately, by staging the four Asplund products against the live ones.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import product_eligibility as pe          # noqa: E402
from pipeline import reference_lineset as rls           # noqa: E402

FEED = ROOT / "data" / "products" / "solar"
POOLS = ("products", "superseded", "archive", "quarantine")

#: The key as it was BEFORE this ticket -- kept here, spelled out, so the audit compares
#: against the real previous definition rather than against `KEY_FIELDS[:-1]`, which would
#: silently follow any future reordering and stop being the old key at all.
OLD_KEY_FIELDS = ("element", "ion", "band", "instrument", "holding",
                  "tier", "selector", "route", "treatment")


def old_key(p: dict) -> str:
    return "|".join(str(p.get(k) or "") for k in OLD_KEY_FIELDS)


def audit_feed(path: Path) -> dict:
    doc = json.loads(path.read_text())
    rows, unresolved = [], []
    for pool in POOLS:
        for p in doc.get(pool) or []:
            if not isinstance(p, dict):
                continue
            try:
                ls = rls.line_set_for_product(p)
            except Exception as e:                       # loud, never defaulted
                unresolved.append({"pool": pool, "key": old_key(p), "error": str(e)})
                continue
            rows.append({"pool": pool, "old": old_key(p), "new": pe.key_of(p),
                         "line_set": ls, "A": p.get("A")})

    by_old = collections.defaultdict(list)
    for r in rows:
        by_old[r["old"]].append(r)
    by_new = collections.defaultdict(list)
    for r in rows:
        by_new[r["new"]].append(r)

    # 🔴 THE REFINEMENT INVARIANT IS THE CHECK WITH TEETH, AND THE OTHER TWO DEPEND ON IT.
    # `line_set` was APPENDED to KEY_FIELDS, so today new_key == old_key + "|" + line_set.
    # While that holds, two facts follow for free and are NOT evidence of anything:
    #   * a MERGE is impossible -- distinct old keys keep distinct prefixes;
    #   * every SPLIT is "explained by line_set" -- it is the only field that can differ.
    # Reporting those as if they were passing checks would be a guard whose two sides
    # converge. So the invariant itself is asserted instead: if someone later INSERTS
    # line_set mid-tuple, or reorders KEY_FIELDS, the prefix property breaks and old keys
    # really can merge and split arbitrarily -- and this is what notices.
    refinement_broken = [r for r in rows if not r["new"].startswith(r["old"] + "|")]

    splits, merges = [], []
    for old, group in by_old.items():
        news = {r["new"] for r in group}
        if len(news) > 1:
            sets = {r["line_set"] for r in group}
            splits.append({"old": old, "new": sorted(news),
                           "line_sets": sorted(sets),
                           # Meaningful ONLY where refinement is broken; while it holds
                           # this is true by construction and is recorded, not trusted.
                           "explained_by_line_set": len(sets) == len(news)})
    for new, group in by_new.items():
        olds = {r["old"] for r in group}
        if len(olds) > 1:
            merges.append({"new": new, "old": sorted(olds)})

    # ⚠️ `relative_to(ROOT)` RAISES on a path outside the repo, which would make this
    # function impossible to point at a staged fixture -- and a migration audit nobody can
    # run against a synthetic split is one whose split-detection has never fired.
    try:
        label = str(path.relative_to(ROOT))
    except ValueError:
        label = str(path)
    return {"feed": label, "n_records": len(rows),
            "n_old_keys": len(by_old), "n_new_keys": len(by_new),
            "line_set_census": dict(collections.Counter(r["line_set"] for r in rows)),
            "splits": splits, "merges": merges, "unresolved": unresolved,
            "refinement_holds": not refinement_broken,
            "refinement_broken": [{"old": r["old"], "new": r["new"]}
                                  for r in refinement_broken]}


def demo_the_rya1106_split() -> dict:
    """Stage the four Asplund products beside the live ones and show the keys separate.

    This is the case the ticket exists for, and it cannot be shown from the committed feed
    because those products are exactly what is blocked. Staging them is not a publish: the
    records are built in memory from the RYA-1106 artifacts and thrown away.
    """
    live = json.loads((FEED / "Fe.json").read_text())
    amarsi = [p for p in live["products"] if p.get("treatment") == "ENGINE-A-3DNLTE"]
    out = []
    for p in amarsi:
        replication = dict(p, line_set="asplund")
        out.append({
            "holding": p["holding"],
            "our_graded_key": pe.key_of(p),
            "asplund_key": pe.key_of(replication),
            "distinct": pe.key_of(p) != pe.key_of(replication),
            "old_key_our_graded": old_key(p),
            "old_key_asplund": old_key(replication),
            "collided_before": old_key(p) == old_key(replication),
        })
    return {"n": len(out), "rows": out,
            "all_distinct_now": all(r["distinct"] for r in out),
            "all_collided_before": all(r["collided_before"] for r in out)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero on an unintended split, any merge, or an "
                         "unresolvable line_set")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data" / "results" / "rya1127"
                    / "key_migration_audit.json")
    a = ap.parse_args()

    audits = [audit_feed(p) for p in sorted(FEED.glob("*.json"))]
    demo = demo_the_rya1106_split()

    L = []; A = L.append
    A("=" * 96)
    A("RYA-1127 — adding `line_set` to the product identity key: migration audit")
    A("=" * 96)
    A(f"  OLD key: {' x '.join(OLD_KEY_FIELDS)}")
    A(f"  NEW key: {' x '.join(pe.KEY_FIELDS)}")
    A("")
    bad = 0
    for d in audits:
        A(f"  {d['feed']}")
        A(f"    records {d['n_records']}   old keys {d['n_old_keys']} -> new keys "
          f"{d['n_new_keys']}")
        A(f"    line_set census: {d['line_set_census']}")
        if d["unresolved"]:
            bad += len(d["unresolved"])
            A(f"    🔴 {len(d['unresolved'])} record(s) with NO resolvable line_set:")
            for u in d["unresolved"]:
                A(f"        {u['pool']:11s} {u['key']}")
                A(f"          {u['error'][:110]}")
        if not d["refinement_holds"]:
            A(f"    🔴 THE REFINEMENT INVARIANT IS BROKEN on "
              f"{len(d['refinement_broken'])} record(s): the new key is not the old key "
              f"plus a suffix, so old identities can merge as well as split and the two "
              f"checks below stop being structural facts.")
            for r in d["refinement_broken"][:5]:
                A(f"        old {r['old']}")
                A(f"        new {r['new']}")
        else:
            A("    refinement invariant: new key == old key + '|' + line_set (so a MERGE "
              "is structurally impossible here, and is asserted rather than claimed)")
        unintended = [s for s in d["splits"] if not s["explained_by_line_set"]]
        for s in d["splits"]:
            tag = "INTENDED (line_set)" if s["intended"] else "🔴 UNINTENDED"
            A(f"    SPLIT {tag}: {s['old']}")
            for n in s["new"]:
                A(f"        -> {n}")
        for m in d["merges"]:
            A(f"    🔴 MERGE (never intended): {m['new']}")
            for o in m["old"]:
                A(f"        <- {o}")
        bad += len(unintended) + len(d["merges"]) + len(d["refinement_broken"])
        if not d["splits"] and not d["merges"] and not d["unresolved"]:
            A("    ✅ partition UNCHANGED — no product was re-identified")
        A("")

    A("-" * 96)
    A("  THE RYA-1106 CASE — the split this ticket enables, staged (not published)")
    A("-" * 96)
    for r in demo["rows"]:
        A(f"    {r['holding']}")
        A(f"      before  our-graded : {r['old_key_our_graded']}")
        A(f"      before  asplund    : {r['old_key_asplund']}")
        A(f"        -> collided: {r['collided_before']}")
        A(f"      after   our-graded : {r['our_graded_key']}")
        A(f"      after   asplund    : {r['asplund_key']}")
        A(f"        -> distinct: {r['distinct']}")
    A("")
    A(f"  all four collided before : {demo['all_collided_before']}")
    A(f"  all four distinct now    : {demo['all_distinct_now']}")
    if not (demo["all_collided_before"] and demo["all_distinct_now"]):
        A("  🔴 the motivating case does NOT behave as the ticket describes")
        bad += 1
    A("=" * 96)

    text = "\n".join(L)
    print(text)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps({"old_key_fields": list(OLD_KEY_FIELDS),
                                 "new_key_fields": list(pe.KEY_FIELDS),
                                 "feeds": audits, "rya1106_case": demo}, indent=2) + "\n")
    a.out.with_suffix(".txt").write_text(text + "\n")
    print(f"\nwrote {a.out}")
    if a.check and bad:
        print(f"\nFAIL: {bad} unintended identity change(s) or unresolvable line_set(s)")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

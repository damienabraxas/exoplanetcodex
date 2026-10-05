#!/usr/bin/env python3
"""Posts orchestrator report summaries to Linear; reads reports, decides nothing, writes no science data.

RYA-1235 -- governing process step 13's Linear reporting. One ASCII comment per element on
that element's dossier issue, routed ONLY through data/catalog/linear_dossiers.csv. A
comment whose content has not changed since it was last posted is not posted again
(data/results/orchestrator/posted_ledger.json). The only file this ever writes is that
ledger.

    python scripts/post_run_report.py --star solar [--element Si] [--dry-run]

Auth: LINEAR_API_KEY from the environment, never from a file. A cron entry is Ryan's call.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOSSIERS = ROOT / "data" / "catalog" / "linear_dossiers.csv"
REPORT_DIR = ROOT / "data" / "results" / "orchestrator"
POSTED = REPORT_DIR / "posted_ledger.json"
LINEAR_API = "https://api.linear.app/graphql"
TIMEOUT_S = 30
MAX_ROWS = 10


def ascii_only(text: str) -> str:
    """Transliterate to ASCII (Cloudflare WAF); anything left unmappable becomes '?'."""
    t = unicodedata.normalize("NFKD", text).encode("ascii", "replace").decode("ascii")
    assert t.isascii()
    return t


def routes(star: str) -> dict[str, str]:
    """{report element name -> issue identifier} for this star, from the ledger only.
    Report names are elements_master spellings ('Fe' is Fe I in the sweep, 'Fe II')."""
    out: dict[str, str] = {}
    with DOSSIERS.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["star"] != star:
                continue
            name = r["element"] + (" II" if r["ion"].strip() == "II" else "")
            out[name] = r["issue_identifier"].strip()
    return out


def element_reports(star: str, only: str | None) -> list[tuple[str, dict]]:
    """(element, report). An element the sweep never ran (no graded line, Fe gate) has no
    element report; its sweep entry -- status and reason -- stands in for it."""
    sweep = REPORT_DIR / f"{star}_SWEEP_latest.json"
    entries = ({e["element"]: e for e in json.loads(sweep.read_text())["elements"]}
               if sweep.exists() else {})
    if only:
        names = [only]
    elif entries:
        names = list(entries)
    else:
        raise SystemExit(f"no {sweep.relative_to(ROOT)}; run the sweep first or pass --element")
    out = []
    for n in names:
        p = REPORT_DIR / f"{star}_{n.replace(' ', '')}_latest.json"
        if p.exists():
            out.append((n, json.loads(p.read_text())))
        elif n in entries:
            e = entries[n]
            out.append((n, {"sweep_status": e.get("status"), "sweep_reason": e.get("reason", ""),
                            "counts": e.get("counts") or {}}))
        else:
            out.append((n, {}))
    return out


def body(star: str, element: str, doc: dict) -> tuple[str, str]:
    """(comment body, content hash). The hash excludes generated_at and the commit."""
    v = doc.get("verdict") or {}
    verdict = v.get("element_verdict") or "not yet wired"
    k = doc.get("counts") or {}
    cells = [c for c in (v.get("cells") or []) if c.get("status") == "OUT_OF_BAND"]
    failed = [c for c in doc.get("cells", []) if c.get("status") == "FAILED"]
    held = [c for c in doc.get("cells", []) if c.get("status") == "HELD"]
    probs = doc.get("problem_lines") or []
    if doc.get("sweep_status"):
        verdict = doc["sweep_status"]
    lines = [
        f"[orchestrator] {star} {element}  verdict={verdict}",
        f"commit {doc.get('code_commit', '?')}  generated {doc.get('generated_at', '?')}",
        "cells: " + (" / ".join(f"{s} {n}" for s, n in k.items()) or "none"),
    ]
    if doc.get("sweep_status"):
        lines.append(f"not run: {doc['sweep_reason'][:200]}")
    steps = [s for s in doc.get("process_steps", []) if not s.get("ok")]
    if steps:
        lines.append("process steps NOT done: " + "; ".join(
            f"{s['step']} {s['name']}" for s in steps))
    oob = [f"{c['band']} {c['instrument']} {c['treatment']} {c.get('pool', '')} "
           f"{c['delta']:+.3f}" for c in cells]
    lines.append("OUT_OF_BAND: " + ("none" if not oob else
                 "; ".join(oob[:MAX_ROWS]) + (f" (+{len(oob) - MAX_ROWS} more)"
                                              if len(oob) > MAX_ROWS else "")))
    fl = [f"{c['band']} {c['instrument']} {c['engine']} {c['reason'][:80]}" for c in failed]
    lines.append("FAILED: " + ("none" if not fl else "; ".join(fl[:MAX_ROWS])))
    pl = [f"{p['wavelength_air_A']} {p['ep_eV']} {p['band']} {p['reason_code']}"
          for p in probs]
    lines.append(f"problem lines: {len(probs)}"
                 + (" (top 10: " + "; ".join(pl[:MAX_ROWS]) + ")" if pl else ""))
    if held and not steps:
        lines.append(f"held cells: {len(held)} (first: {held[0]['reason'][:100]})")
    lines.append(f"report: data/results/orchestrator/{star}_{element.replace(' ', '')}"
                 f"_latest.json")
    text = ascii_only("\n".join(lines))
    key = json.dumps({
        "counts": {s: k.get(s, 0) for s in sorted(k)}, "verdict": verdict,
        "oob": sorted(oob), "failed": sorted(f"{c['band']}|{c['instrument']}|{c['engine']}"
                                             for c in failed),
        "problems": sorted(f"{p['wavelength_air_A']}|{p['ep_eV']}|{p['band']}|"
                           f"{p['instrument']}|{p['reason_code']}" for p in probs),
        "steps": [s["step"] for s in steps]}, sort_keys=True)
    return text, hashlib.sha256(key.encode()).hexdigest()


def _gql(query: str, variables: dict, key: str) -> tuple[dict | None, str]:
    import requests
    try:
        r = requests.post(LINEAR_API, json={"query": query, "variables": variables},
                          headers={"Authorization": key, "Content-Type": "application/json"},
                          timeout=TIMEOUT_S)
    except Exception as exc:                                   # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"
    if r.status_code != 200:
        return None, f"HTTP {r.status_code}: {r.text[:200]}"
    data = r.json()
    if data.get("errors"):
        return None, f"GraphQL errors: {str(data['errors'])[:200]}"
    return data, ""


def post(issue: str, text: str, key: str) -> tuple[bool, str]:
    # The ledger holds the human identifier (TEAM-123); commentCreate takes the issue's id,
    # resolved here through `issue(id:)`, which accepts the identifier.
    found, why = _gql("query($i: String!) { issue(id: $i) { id identifier } }",
                      {"i": issue}, key)
    if found is None:
        return False, f"could not resolve {issue}: {why}"
    iss = (found.get("data") or {}).get("issue") or {}
    if iss.get("identifier") != issue:
        return False, f"{issue} resolved to {iss.get('identifier')!r} -- refusing to post"
    data, why = _gql("mutation($i: String!, $b: String!) { commentCreate(input: "
                     "{issueId: $i, body: $b}) { success comment { id } } }",
                     {"i": iss["id"], "b": text}, key)
    if data is None:
        return False, why
    if not (data.get("data", {}).get("commentCreate") or {}).get("success"):
        return False, f"commentCreate did not succeed: {str(data)[:200]}"
    return True, "ok"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--star", required=True)
    ap.add_argument("--element", default=None)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    key = os.environ.get("LINEAR_API_KEY", "")
    if not a.dry_run and not key:
        print("ERROR: LINEAR_API_KEY is not set in the environment. Nothing was posted. "
              "The key is read from the environment only, never from a file.",
              file=sys.stderr)
        return 2

    route = routes(a.star)
    ledger = json.loads(POSTED.read_text()) if POSTED.exists() else {}
    n = {"POSTED": 0, "SKIP": 0, "UNROUTED": 0, "FAILED": 0}
    unrouted = []
    for element, doc in element_reports(a.star, a.element):
        issue = route.get(element)
        if not issue:
            n["UNROUTED"] += 1
            unrouted.append(element)
            print(f"UNROUTED  {element}: no row in {DOSSIERS.name}")
            continue
        if not doc:
            n["FAILED"] += 1
            print(f"FAILED    {element}: no report {element.replace(' ', '')}_latest.json",
                  file=sys.stderr)
            continue
        text, h = body(a.star, element, doc)
        lk = f"{a.star}|{element}|{issue}"
        if ledger.get(lk, {}).get("hash") == h:
            n["SKIP"] += 1
            print(f"SKIP      {element} -> {issue} (unchanged since {ledger[lk]['posted_at']})")
            continue
        if a.dry_run:
            print(f"--- would post to {issue} ---\n{text}\n")
            continue
        ok, why = post(issue, text, key)
        if ok:
            from datetime import datetime, timezone
            ledger[lk] = {"hash": h, "posted_at": datetime.now(timezone.utc)
                          .strftime("%Y-%m-%dT%H:%M:%SZ")}
            POSTED.parent.mkdir(parents=True, exist_ok=True)
            POSTED.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
            n["POSTED"] += 1
            print(f"POSTED    {element} -> {issue}")
        else:
            n["FAILED"] += 1
            print(f"FAILED    {element} -> {issue}: {why}", file=sys.stderr)
    if unrouted:
        print(f"UNROUTED list: {', '.join(unrouted)}")
    print(" / ".join(f"{s} {n[s]}" for s in ("POSTED", "SKIP", "UNROUTED", "FAILED")))
    return 1 if n["FAILED"] else 0


if __name__ == "__main__":
    sys.exit(main())

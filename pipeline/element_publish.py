"""Governing-process steps 11-16 as ONE orchestrator stage: budget, assemble, publish, check.

    python run_pipeline.py --star solar --element Si --publish [--jobs 8]

RYA-1233. On Si every one of these steps was a hand-run per-ticket script with a hand-written
unit list and hand-filtered reports, and each hand step broke something once (ENGINE-B
duplicates, a stale pre-cull stage, missing `copied_to`, stale display labels). This stage is
those same, already-proven tools composed, with every filter the hand runs needed made a rule:

  1  units      the derive invocation of every synthesis cell the last run BUILT (DONE /
                UNPUBLISHED / SKIP), taken from the orchestrator's own resolution
  2  legs       scripts/rya1230_cno_budget_legs.py -- the RYA-587 paired legs, N in parallel
  3  assemble   scripts/rya1230_cno_budget.py -- the budget per product; holds stay holds
  4  filter     drop ENGINE-B (on the ts-lte deck it is the 1D-LTE fit relabelled), products
                with no line in the aggregate, and held budgets -- every drop RECORDED
  5  publish    scripts/rya1230_publish_cno.py --apply --origin-dir (copied_to recorded),
                through publish_product's own gates
  6  check      element_verdict.judge_feed (literature) + coverage_report (definition of done)

Everything is written under data/results/orchestrator/budget/<star>_<El>/: legs/ (ignored),
stage/ + report.json + publish_summary.json (committed: the feed's provenance points here).
"""
from __future__ import annotations

import json
import shlex
import subprocess
import sys
from pathlib import Path

from pipeline import run_matrix as rm

ROOT = Path(__file__).resolve().parents[1]
BUILT = (rm.DONE, rm.UNPUBLISHED, rm.SKIP)


def budget_dir(star: str, element: str) -> Path:
    return rm.REPORT_DIR / "budget" / f"{star}_{element.replace(' ', '')}"


def units_for(star: str, element: str, *, report: dict, interpreter: str) -> list[tuple[str, str]]:
    """(cell label, derive command line) for every synthesis cell the run built."""
    from pipeline.run_descriptor import resolve
    symbol, _ = rm.split_symbol(element)
    built = {(c["band"], c["holding"], c["ion"], c.get("pool") or "")
             for c in report.get("cells", []) if c.get("status") in BUILT}
    out = []
    for d in rm.expand(star, symbol, engines=["ts-lte"], methods=["synthesis"],
                       pools=["set", "codex", "deep"]):
        if (rm._band_of(d, star).name_declared, d.holding, d.ion, rm.pool_tier(d)) not in built:
            continue
        r = resolve(d, interpreter=interpreter or None, ispec_dir=None)
        step = next(s for s in r.steps if s["name"] == "derive_products")
        out.append((f"{d.ion} {d.band} {d.holding} {d.pool}",
                    "python3 scripts/derive_band_products.py " + shlex.join(step["args"])))
    return out


def publish(star: str, element: str, *, report: dict, jobs: int = 8,
            interpreter: str | None = None, reason: str | None = None) -> dict:
    py = interpreter or sys.executable
    work = budget_dir(star, element)
    legs, stage = work / "legs", work / "stage"
    work.mkdir(parents=True, exist_ok=True)
    units = units_for(star, element, report=report, interpreter=py)
    summary = {"star": star, "element": element, "units": [u[0] for u in units],
               "dropped": [], "published": [], "held": []}
    if not units:
        summary["note"] = "no built synthesis cell in the run report -- nothing to budget"
        (work / "publish_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
        return summary
    ufile = work / "units.txt"
    ufile.write_text("\n".join(u[1] for u in units) + "\n")
    # A changed unit list invalidates leg directories keyed by unit index (the legs runner
    # skips a leg whose DONE exists, so unchanged units are reused across invocations).
    prev = work / "units.prev.txt"
    if prev.exists() and prev.read_text() != ufile.read_text():
        import shutil
        shutil.rmtree(legs, ignore_errors=True)
    prev.write_text(ufile.read_text())
    subprocess.run([py, "scripts/rya1230_cno_budget_legs.py", "--units", str(ufile),
                    "--out", str(legs), "--jobs", str(jobs)], cwd=ROOT, check=True)
    import shutil
    shutil.rmtree(stage, ignore_errors=True)
    rep = work / "report.json"
    subprocess.run([py, "scripts/rya1230_cno_budget.py", "--legs", str(legs), "--stage", str(stage),
                    "--report", str(rep)], cwd=ROOT, check=True)
    r = json.loads(rep.read_text())
    keep = []
    for o in r["products"]:
        tag = f"{o.get('band')} {o.get('holding')} {o.get('selector')} {o.get('key_treatment')}"
        if o.get("key_treatment") == "ENGINE-B":
            summary["dropped"].append({"product": tag, "why": "ENGINE-B on ts-lte = 1D-LTE relabelled"})
        elif o.get("skip"):
            summary["held"].append({"product": tag, "why": o["skip"]})
        elif o.get("verdict"):
            summary["held"].append({"product": tag, "why": o["verdict"]})
        else:
            keep.append(o)
    r["products"] = keep
    pub = work / "publish.json"
    pub.write_text(json.dumps(r, indent=1, default=str) + "\n")
    why = reason or (f"orchestrator --publish ({rm._code_commit()[:12]}): RYA-587 budget from paired "
                     f"legs on each product's own pool; published through publish_product's gates")
    # publish_product directly, one staged product at a time, through its own gates --
    # with --origin-path so `copied_to` is recorded (feed_repo_reconciliation, RYA-1080).
    done = set()
    for o in keep:
        src = stage / o["prod_stem"]
        if src in done or not src.exists():
            continue
        done.add(src)
        cmd = [py, "scripts/publish_product.py", "--from", str(src), "--holding", o["holding"],
               "--tier", "ALL", "--route", "SYNTH", "--element", o["element"], "--star", star,
               "--reason", why, "--origin-path", str(src)]
        if o.get("selector"):
            cmd += ["--selector", o["selector"]]
        res = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        tag = f"{o.get('band')} {o.get('holding')} {o.get('selector')} {o.get('key_treatment')} A={o.get('A')}"
        if res.returncode == 0:
            summary["published"].append(tag)
        else:
            summary["held"].append({"product": tag, "why": "publish_product refused: "
                                    + " | ".join((res.stdout + res.stderr).strip().splitlines()[-2:])})
    from pipeline import coverage_report, element_verdict
    summary["verdict"] = {k: v for k, v in element_verdict.judge_feed(star, element).items() if k != "cells"}
    cov = coverage_report.coverage(star, element)
    summary["coverage"] = {"counts": cov["counts"], "done": cov["done"]}
    (work / "coverage.json").write_text(json.dumps(cov, indent=1, default=str) + "\n")
    (work / "publish_summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
    print(coverage_report.render(cov))
    return summary

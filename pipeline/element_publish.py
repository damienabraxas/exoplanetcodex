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
  4  filter     keep only what each cell's deck EMITS (run_matrix.DECK_EMITS: ENGINE-B on
                ts-lte is the 1D-LTE fit relabelled; a Gerber deck's 1D-LTE leg duplicates), products
                with no line in the aggregate, and held budgets -- every drop RECORDED
  5  publish    scripts/rya1230_publish_cno.py --apply --origin-dir (copied_to recorded),
                through publish_product's own gates
  6  check      element_verdict.judge_feed (literature) + coverage_report (definition of done)

Everything is written under data/results/orchestrator/budget/<star>_<El>/: legs/ (ignored),
stage/ + report.json + publish_summary.json (committed: the feed's provenance points here).
"""
from __future__ import annotations

import json
import re
import shlex
import subprocess
import sys
from pathlib import Path

from pipeline import run_matrix as rm

ROOT = Path(__file__).resolve().parents[1]
BUILT = (rm.DONE, rm.UNPUBLISHED, rm.SKIP)


def budget_dir(star: str, element: str) -> Path:
    return rm.REPORT_DIR / "budget" / f"{star}_{element.replace(' ', '')}"


def units_for(star: str, element: str, *, report: dict,
              interpreter: str) -> list[tuple[str, str, str]]:
    """(cell label, derive command line, deck) for every synthesis cell the run built -- on
    EVERY deck it built (ts-lte and, where the element has them, the Gerber decks)."""
    from pipeline.run_descriptor import resolve
    symbol, _ = rm.split_symbol(element)
    built = {(c["band"], c["holding"], c["ion"], c.get("engine") or "ts-lte", c.get("pool") or "")
             for c in report.get("cells", []) if c.get("status") in BUILT}
    # RYA-1251: the cells the report names at all, whatever their status. A SET cell the run
    # HELD is not the GRADED cell beside it -- aliasing it budgeted 4 Na NIR cells with no
    # product (36 legs, every one rc=1).
    reported = {(c["band"], c["holding"], c["ion"], c.get("engine") or "ts-lte", c.get("pool") or "")
                for c in report.get("cells", [])}
    out = []
    for d in rm.expand(star, symbol, methods=["synthesis"], pools=["set", "codex", "deep"]):
        identity = (rm._band_of(d, star).name_declared, d.holding, d.ion,
                    d.engine_deck, rm.pool_tier(d))
        # Older run reports represented a published line-set cell as the generic
        # GRADED selector. Accept that spelling when matching the descriptor's
        # SET-* selector so artifacts already built by the orchestrator reach the
        # budget/publication stage. Canonical and deep pools remain exact matches.
        if identity not in built and not (
                identity not in reported
                and d.pool and d.pool.startswith(rm.SET_POOL_PREFIX)
                and (*identity[:-1], "GRADED") in built):
            continue
        r = resolve(d, interpreter=interpreter or None, ispec_dir=None)
        step = next(s for s in r.steps if s["name"] == "derive_products")
        out.append((f"{d.ion} {d.band} {d.holding} {d.engine_deck} {d.pool}",
                    "python3 scripts/derive_band_products.py " + shlex.join(step["args"]),
                    d.engine_deck, rm.ledger_key(star, d, identity[0])))
    return out


def _observed_conditioning(o: dict) -> str | None:
    """The single `observed_conditioning` value of the product's own per-line artifact."""
    import pandas as pd
    stem = o.get("stem") or ""
    # A correction leg (ENGINE-A ...) is built from the SAME observed artifact as the unit's
    # 1D-LTE leg, whose per-line rows carry the column.
    lte = re.sub(r"_[A-Za-z0-9.-]+_lines\.csv$", "_1D-LTE_lines.csv", stem)
    for f in [ROOT / "data/results/band_products" / s for s in dict.fromkeys([stem, lte])] + \
            [g for s in dict.fromkeys([stem, lte])
             for g in sorted((ROOT / "data/results/band_products").glob(f"deck_*/{s}"))]:
        if stem and f.exists():
            d = pd.read_csv(f)
            if "observed_conditioning" in d.columns:
                vals = sorted(set(d["observed_conditioning"].dropna().astype(str)))
                if len(vals) == 1:
                    return vals[0]
                if len(vals) > 1:
                    return None              # ambiguous: refuse rather than pick
                # empty on a correction leg: fall through to the unit's 1D-LTE artifact
    return None


def _identity(row: dict) -> tuple:
    return tuple(str(row.get(k) or "") for k in
                 ("ion", "band", "instrument", "holding", "selector", "treatment"))


def _built_from(rec: dict | None) -> str | None:
    """A unit's legs are keyed on what its cell was built from: the inputs hash, plus the
    forced-rebuild stamp where `--force` rebuilt it at an unchanged hash (RYA-1315)."""
    rec = rec or {}
    h = rec.get("inputs_hash")
    return f"{h} forced {rec['forced_at']}" if rec.get("forced_at") else h


def _scope_cells(report: dict) -> dict[tuple, set]:
    """(ion, band, holding, selector) -> treatments, for every cell a narrowed run BUILT."""
    out: dict[tuple, set] = {}
    for c in report.get("cells", []):
        if c.get("status") in BUILT:
            out.setdefault((c["ion"], c["band"], c["holding"], c.get("pool") or ""), set()).update(
                rm.DECK_EMITS.get(c.get("engine") or "ts-lte", ()))
    return out


def _withdraw_superseded(star: str, element: str, *, stage: Path, keep: list,
                         held: list, why: str, started: str,
                         scope: dict[tuple, set] | None = None) -> list[str]:
    """Move every live SYNTH product this run did NOT republish to `archive`.

    RYA-1233 Si: a re-measured cell whose new budget HOLDS (or whose cell no longer exists --
    a narrowed window) left its previous value live, pointing at a stage file this stage had
    just deleted, and coverage still counted it published. This stage owns the element's
    synthesis products: the current run is the only source of truth for them, so a product
    it did not publish is superseded -- archived (never deleted) with the run's reason.
    """
    import pandas as pd
    sys.path.insert(0, str(ROOT / "scripts"))
    import publish_product as pp
    current = set()
    for o in keep:
        src = stage / o["prod_stem"]
        if src.exists():
            for r in pd.read_csv(src).to_dict("records"):
                current.add(_identity({**r, "holding": o["holding"],
                                       "selector": o.get("selector")}))
    why_held = {h["product"]: h["why"] for h in held}
    path = pp.STORE / star / f"{element}.json"
    if not path.exists():
        return []
    doc = json.loads(path.read_text())
    # A rerun the eligibility gate sent to quarantine (RYA-1097) leaves the earlier record
    # live; the current run's verdict is the one that stands, so that record goes too.
    quarantined_now = {pp.key_of(q): q for q in doc.get("quarantine") or []
                       if (q.get("quarantined_at") or "") >= started}
    gone = [p for p in doc["products"] if p.get("route") == "SYNTH"
            and p.get("star", star) == star
            and (_identity(p) not in current or pp.key_of(p) in quarantined_now)]
    if scope is not None:
        # RYA-1315: a narrowed run withdraws only what it re-measured and could not stand
        # behind -- a product in its scope that it HELD or quarantined. One it did not
        # re-emit at all is left live; out-of-scope products are never considered.
        def _held(p):
            return (f"{p.get('band')} {p.get('holding')} {p.get('selector')} {p.get('treatment')}"
                    in why_held or pp.key_of(p) in quarantined_now)
        gone = [p for p in gone if str(p.get("treatment")) in scope.get(
                    (str(p.get("ion")), str(p.get("band")), str(p.get("holding")),
                     str(p.get("selector") or "")), ()) and _held(p)]
    if not gone:
        return []
    now = pp._now()
    for p in gone:
        tag = f"{p.get('band')} {p.get('holding')} {p.get('selector')} {p.get('treatment')}"
        doc.setdefault("archive", []).append(dict(
            p, superseded_at=now, superseded_reason=(
                f"{why}. Superseded: the current run "
                + ("QUARANTINED its rerun ("
                   + (quarantined_now[pp.key_of(p)].get("quarantine_reason") or "")[-200:] + ")"
                   if pp.key_of(p) in quarantined_now else
                   f"HOLDS this product ({why_held[tag]})" if tag in why_held
                   else "no longer builds this cell") + ", so the earlier value is withdrawn")))
    doc["products"] = [p for p in doc["products"] if p not in gone]
    doc["version"] = pp.bump(doc["version"])
    doc["updated_at"] = now
    pp.write_feed(path, doc)
    return [f"{pp.key_of(p)} A={p.get('A')}" for p in gone]


def _preserve_history(star: str, element: str, *, work: Path, prev_stage: Path) -> list[str]:
    """Keep the artifact of every archived / quarantined record that names one.

    The stage is rebuilt on every run, so a superseded record's `copied_to` pointed at a file
    the rebuild deleted -- the history could no longer be checked against its own bytes
    (feed_repo_reconciliation, RYA-1080/1096). The file is copied to `superseded/` from the
    previous stage (or, failing that, the last commit holding it), verified by sha256, and the
    record repointed. Measured fields are never touched.
    """
    import hashlib
    import shutil
    sys.path.insert(0, str(ROOT / "scripts"))
    import publish_product as pp
    path = pp.STORE / star / f"{element}.json"
    if not path.exists():
        return []
    doc = json.loads(path.read_text())
    keep_dir = work / "superseded"
    done = []
    for sec in ("archive", "quarantine"):
        for rec in doc.get(sec) or []:
            prov = rec.get("provenance") or {}
            ct, sha = prov.get("copied_to"), prov.get("sha256")
            if not ct or (ROOT / ct).exists():
                continue
            name = Path(ct).name
            dst = keep_dir / f"{(sha or 'nosha')[:12]}_{name}"
            if not dst.exists():
                data = None
                cand = prev_stage / name
                if cand.exists():
                    data = cand.read_bytes()
                if data is None or (sha and hashlib.sha256(data).hexdigest() != sha):
                    data = None
                    revs = subprocess.run(["git", "log", "--all", "--format=%H", "--", ct], cwd=ROOT,
                                          capture_output=True, text=True).stdout.split()
                    for rev in revs:
                        got = subprocess.run(["git", "show", f"{rev}:{ct}"], cwd=ROOT,
                                             capture_output=True)
                        if got.returncode == 0 and (not sha or
                                                    hashlib.sha256(got.stdout).hexdigest() == sha):
                            data = got.stdout
                            break
                if data is None:
                    continue                 # left as is: reconciliation keeps naming it
                keep_dir.mkdir(parents=True, exist_ok=True)
                dst.write_bytes(data)
            rec["provenance"] = {**prov, "copied_to": str(dst.relative_to(ROOT)),
                                 "copied_to_was": ct}
            done.append(f"{sec}: {name}")
    if done:
        doc["version"] = pp.bump(doc["version"])
        doc["updated_at"] = pp._now()
        pp.write_feed(path, doc)
    return done


def unit_signature(command: str, inputs_hash: str | None) -> str:
    """What a unit's legs were built from: its derive command line and its cell's inputs hash."""
    return f"{command}\n# {inputs_hash}"


def reuse_legs(legs: Path, prev_file: Path, sigs: list[str]) -> list[int]:
    """Keep each unit's legs while its signature is unchanged, renumbered to the unit's new
    index. Returns the new indices whose legs were kept.

    A previous unit absent from this run is FILED under its signature (`legs.store/`), never
    deleted, and comes back when its signature does: RYA-1250 ran publish on a one-band run
    report and lost 132 of 172 Mg units' legs (~30 h). Only a unit whose inputs changed is
    stale, and it stays filed under its OLD signature, where nothing current matches it.

    No per-unit record (a first run, or the pre-RYA-1250 whole-element `units.prev.txt`)
    proves nothing current, so those legs go. A signature listed twice is ambiguous and is
    never reused."""
    import hashlib
    import shutil

    def _file(s: str) -> str:
        return hashlib.sha256(s.encode()).hexdigest()[:32]
    prev = json.loads(prev_file.read_text()) if prev_file.exists() else None
    store = legs.parent / "legs.store"
    kept: list[int] = []
    if prev is None:
        shutil.rmtree(legs, ignore_errors=True)
    else:
        old = {s: j for j, s in enumerate(prev) if prev.count(s) == 1}
        hold = legs.parent / "legs.reuse"
        shutil.rmtree(hold, ignore_errors=True)
        hold.mkdir(parents=True)
        for i, s in enumerate(sigs):
            if sigs.count(s) != 1:
                continue
            j = old.pop(s, None)
            src = legs / f"unit{j:02d}" if j is not None else store / _file(s)
            if src.is_dir():
                src.rename(hold / f"unit{i:02d}")
                kept.append(i)
        store.mkdir(parents=True, exist_ok=True)
        for s, j in old.items():                       # not in this run: file, never delete
            src, dst = legs / f"unit{j:02d}", store / _file(s)
            if src.is_dir():
                shutil.rmtree(dst, ignore_errors=True)
                src.rename(dst)
        shutil.rmtree(legs, ignore_errors=True)
        hold.rename(legs)
    prev_file.parent.mkdir(parents=True, exist_ok=True)
    prev_file.write_text(json.dumps(sigs, indent=0) + "\n")
    return kept


def publish(star: str, element: str, *, report: dict, jobs: int = 8,
            interpreter: str | None = None, reason: str | None = None) -> dict:
    py = interpreter or sys.executable
    # RYA-1250: a one-band report once made publish treat 132 of 172 Mg units as gone.
    # RYA-1315: a NARROWED report publishes its own cells only -- legs are filed by signature
    # (never deleted), the stage is written into rather than rebuilt (every live product's
    # `copied_to` points there), and only products INSIDE the scope can be withdrawn.
    scoped = bool(report.get("scope"))
    if report.get("dry_run"):
        raise ValueError("publish refuses a dry-run report: it built nothing to price")
    if report.get("star") not in (None, star) or report.get("element") not in (None, element):
        raise ValueError(f"publish: report is for {report.get('star')} {report.get('element')}, "
                         f"not {star} {element}")
    import datetime as _dt
    started = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
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
    # Leg directories are keyed by unit index and the legs runner skips a leg whose DONE
    # exists, so legs are reused across invocations -- ONLY while the unit's command line and
    # its own cell's inputs hash are unchanged. Keyed on the command line alone, a canonical_gf
    # or synthesis-list change re-built the cells and then priced them from stale legs.
    # RYA-1250: keyed on the whole element, culling one K line re-priced all 172 Mg units.
    ledger = json.loads(rm.LEDGER.read_text()) if rm.LEDGER.exists() else {}
    reuse_legs(legs, work / "units.prev.json",
               [unit_signature(u[1], _built_from(ledger.get(u[3]))) for u in units])
    # RYA-1251: a leg whose every line failed SYNTHESIS measured nothing -- an engine failure
    # the signature above cannot see (a library module is not in `inputs_hash`). Re-run it.
    for done in legs.glob("unit*/*/DONE"):
        arts = [{"path": str(p.resolve())} for p in done.parent.glob("*_lines.csv")]
        if rm._every_line_failed_synthesis(arts):
            done.unlink()
    # Legs write DONE only on success; a leg killed mid-run (RYA-1233 Si: 16 legs, logs cut off
    # with no traceback, under memory pressure) is simply re-run -- once -- and whatever is
    # still missing is RECORDED, never silently assembled around.
    for _attempt in range(2):
        subprocess.run([py, "scripts/rya1230_cno_budget_legs.py", "--units", str(ufile),
                        "--out", str(legs), "--jobs", str(jobs)], cwd=ROOT, check=True)
        missing = sorted(str(d.relative_to(legs)) for d in legs.glob("unit*/*")
                         if d.is_dir() and not (d / "DONE").exists())
        if not missing:
            break
    summary["legs_incomplete"] = missing
    import shutil
    prev_stage = work / "stage.prev"
    shutil.rmtree(prev_stage, ignore_errors=True)
    if scoped:
        stage.mkdir(parents=True, exist_ok=True)
        shutil.copytree(stage, prev_stage)  # the out-of-scope files stay where the feed points
    elif stage.exists():
        stage.rename(prev_stage)        # kept until superseded records are preserved
    rep = work / "report.json"
    subprocess.run([py, "scripts/rya1230_cno_budget.py", "--legs", str(legs), "--stage", str(stage),
                    "--report", str(rep)], cwd=ROOT, check=True)
    r = json.loads(rep.read_text())
    keep = []
    for o in r["products"]:
        tag = f"{o.get('band')} {o.get('holding')} {o.get('selector')} {o.get('key_treatment')}"
        deck = units[o["unit"]][2] if isinstance(o.get("unit"), int) and o["unit"] < len(units) else "ts-lte"
        if o.get("key_treatment") not in rm.DECK_EMITS.get(deck, ()):
            # A cell publishes only what its deck EMITS (run_matrix.DECK_EMITS): ENGINE-B on
            # ts-lte is the 1D-LTE fit relabelled, and a Gerber deck's 1D-LTE / ENGINE-A legs
            # duplicate the ts-lte cell's.
            summary["dropped"].append({"product": tag, "why": f"not emitted by deck {deck} "
                                       f"(it emits {', '.join(rm.DECK_EMITS.get(deck, ()))})"})
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
        tier = o.get("selector") if o.get("selector") in ("GRADED", "DEEPGRADED") else "ALL"
        cmd = [py, "scripts/publish_product.py", "--from", str(src), "--holding", o["holding"],
               "--tier", tier, "--route", "SYNTH", "--element", o["element"], "--star", star,
               "--reason", why, "--origin-path", str(src)]
        if o.get("selector"):
            cmd += ["--selector", o["selector"]]
        cond = _observed_conditioning(o)
        if cond:
            # RYA-1176: Al publishing requires the observed-conditioning identity, read from
            # the measured artifact itself (derive writes it into every per-line row).
            cmd += ["--observed-conditioning", cond]
        res = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        tag = f"{o.get('band')} {o.get('holding')} {o.get('selector')} {o.get('key_treatment')} A={o.get('A')}"
        if res.returncode == 0:
            summary["published"].append(tag)
        else:
            summary["held"].append({"product": tag, "why": "publish_product refused: "
                                    + " | ".join((res.stdout + res.stderr).strip().splitlines()[-2:])})
    summary["withdrawn"] = _withdraw_superseded(
        star, element, stage=stage, keep=keep, held=summary["held"], why=why, started=started,
        scope=_scope_cells(report) if scoped else None)
    if scoped:
        summary["scope"] = report["scope"]
    summary["preserved"] = _preserve_history(star, element, work=work, prev_stage=prev_stage)
    shutil.rmtree(prev_stage, ignore_errors=True)
    from pipeline import coverage_report, element_verdict
    summary["verdict"] = {k: v for k, v in element_verdict.judge_feed(star, element).items() if k != "cells"}
    cov = coverage_report.coverage(star, element)
    summary["coverage"] = {"counts": cov["counts"], "done": cov["done"]}
    (work / "coverage.json").write_text(json.dumps(cov, indent=1, default=str) + "\n")
    (work / "publish_summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
    print(coverage_report.render(cov))
    return summary

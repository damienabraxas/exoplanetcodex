"""The definition of DONE for one element: every cell, every treatment, an outcome (RYA-1233).

Ryan, 2026-10-08: "did we hit every band? every instrument? every engine? with uncertainty?"
-- a question the orchestrator must answer itself. This table is that answer, built from the
orchestrator's own expansion (`run_matrix.expand`), its process holds, the last run report and
the published feed. Every (band x holding x ion x grade-pool x treatment) cell gets exactly one
outcome:

  PUBLISHED   a product in data/products/<star>/<El>.json, with its A, sigma and n
  HELD        the cell, or that treatment, cannot be produced -- with the stated reason
  NOT_RUN     applicable and not yet measured -- work is owed
  STALE       published, but its inputs (lines, gf, spectrum, culls, NLTE table) changed
              since the build -- work is owed
  CODE_DRIFT  published; only the code changed since the build, no input did -- reported,
              not owed (Ryan 2026-11-03: run only the cells a confirmed defect touches)

RYA-1349 (Ryan 2026-11-07: "Did we do Engine A? Did we do CRIRES?"): a missing ENGINE-A
product is HELD only when its grid really holds none of the pool's lines, a published
ENGINE-A product that corrected no line is owed, and products of the same lines on
different holdings that disagree are FLAGGED for a decision (`disagreements`).

A cell's treatments are its deck's emissions (run_matrix.DECK_EMITS); ENGINE-B on the
ts-lte deck is the 1D-LTE fit relabelled and is not a separate cell.
An element is DONE when no cell is NOT_RUN or STALE. ELEMENT_PROTOCOL (docs/) "three honest outcomes"
map here: RESOLVED -> PUBLISHED, BOUNDED/UNRESOLVED -> HELD with its reason.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from pipeline import band_policy
from pipeline import run_matrix as rm

#: Band order is the band policy's own (blue to red), never typed here.
BAND_ORDER = [p.name for p in band_policy.POLICIES]

#: A cell's treatments are what its deck EMITS (run_matrix.DECK_EMITS): ts-lte -> 1D-LTE +
#: ENGINE-A; each Gerber deck -> its own token. ENGINE-B on ts-lte is not a cell.
#: Orchestrator pool -> the feed's `selector` token.
POOL_SELECTOR = {"codex": "GRADED", "deep": "DEEPGRADED", "reference": "REFERENCE"}
#: RYA-1213 / Ryan 2026-10-21: Reference Grade = the published (Asplund) set PLUS every
#: lab-gf line -- the `reference` pool is the second half of Reference Grade, not a pool of
#: its own and not Codex + Deep.
GRADE_OF_POOL = {"codex": "Codex Grade", "deep": "Deep Grade", "reference": "Reference Grade"}


#: Ryan 2026-11-02/03 (docs/ORCHESTRATOR_RUNBOOK.md "only Deep is dropped"): no NEW Deep
#: product; an existing one is corrected, never created. An unbuilt Deep cell owes nothing.
DEEP_POLICY = "Deep dropped (Ryan 2026-11-03): no new Deep product; existing ones are only corrected"
#: The canonical_gf fields a fit or its synthesis reads. A source tag or gf sigma edit
#: (RYA-1311 cited Ca II 5001.48 to Asplund, log gf unchanged) moves none of them.
GF_PHYSICS = ("physical_id", "species", "wavelength_air_A", "excitation_potential_eV",
              "log_gf", "hfs_n_components")
#: Products of the same ion/band/grade/treatment on different holdings: one is FLAGGED when
#: it sits further than this many of the OTHER holdings' (median) sigma from their median.
#: Judged on the agreeing holdings' precision, so an outlier's own inflated sigma cannot
#: excuse it (Ti II NIR KP-Kurucz: 5.30 +/- 0.38 against 4.81 / 4.84, RYA-1315).
DISAGREE_K = 2.0


def _git_text(commit: str | None, rel: str, _cache: dict = {}) -> str | None:
    """`rel` as committed at `commit` (None: the working tree); None if unreadable."""
    import subprocess
    key = (commit, rel)
    if key not in _cache:
        try:
            if commit is None:
                _cache[key] = (rm.ROOT / rel).read_text(encoding="utf-8")
            else:
                res = subprocess.run(["git", "show", f"{commit}:{rel}"], cwd=rm.ROOT,
                                     capture_output=True, text=True)
                _cache[key] = res.stdout if res.returncode == 0 else None
        except OSError:
            _cache[key] = None
    return _cache[key]


def _git_table(commit: str | None, rel: str, _cache: dict = {}):
    """`_git_text` parsed as strings; None if unreadable."""
    import io
    import pandas as pd
    key = (commit, rel)
    if key not in _cache:
        text = _git_text(commit, rel)
        try:
            _cache[key] = (None if text is None else
                           pd.read_csv(io.StringIO(text), dtype=str, keep_default_na=False))
        except Exception:
            _cache[key] = None
    return _cache[key]


def _built_table(kind: str, rel: str, commit: str | None, digest: str | None,
                 lo: float, hi: float):
    """The table a build read, recovered from its recorded commit -- only when re-hashing that
    commit's file with the build's own scoped fingerprint gives the RECORDED digest. A build
    from a dirty tree records `<sha>-dirty` (RYA-1315 Ti: 7b28c94f-dirty); its data files
    are the commit's exactly when this proof holds. None when it does not."""
    import tempfile
    if not commit or not digest:
        return None
    sha = commit.removesuffix("-dirty")
    text = _git_text(sha, rel)
    if text is None:
        return None
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / Path(rel).name
        f.write_text(text, encoding="utf-8")
        got = (rm._rows_fingerprint(f, "wavelength_air_A", lo, hi) if kind == "ledger_scoped"
               else rm.registry_fingerprint(round(lo), round(hi), f))
    return _git_table(sha, rel) if got == digest else None


def _input_moved(kind: str, name: str, commit: str | None, digest: str | None, symbol: str,
                 window: tuple[float, float], waves: list[float] = ()) -> bool | None:
    """Whether a scoped ledger input changed in a way the cell can see, between the build and
    now. True/False; None when it cannot be told (then the cell is STALE).

    canonical_gf is scoped to the window on purpose -- another species' line there is blend
    opacity (RYA-1233) -- but only its GF_PHYSICS fields reach the synthesis. A cull in
    problem_children changes the pool of its OWN species only. With the cell's own line
    wavelengths, a gf row counts only inside the span a fit reads around one of them (the
    RYA-1230 continuum envelope): Ca II 6456.88's AGSS21 gf sits 51 A from any Ti line."""
    import re
    import pandas as pd
    rel = name.rsplit("[", 1)[0]
    lo, hi = window
    old, new = _built_table(kind, rel, commit, digest, lo, hi), _git_table(None, rel)
    if old is None or new is None:
        return None
    if kind == "ledger_scoped":
        def rows(t):
            w = pd.to_numeric(t["wavelength_air_A"], errors="coerce")
            return set(map(tuple, t[(w >= lo) & (w <= hi)][list(GF_PHYSICS)].values.tolist()))
        moved = rows(old) ^ rows(new)
        if waves:
            from pipeline.local_continuum import ENV_HALF_WIDTH_A
            span = ENV_HALF_WIDTH_A + 0.5          # derive_band_products._observed's read
            moved = {r for r in moved if any(abs(float(r[2]) - w) <= span for w in waves)}
        return bool(moved)

    def culls(t):
        keep = set()
        for r in t.values.tolist():
            m = re.match(r"\s*([0-9]+\.?[0-9]*)", str(r[1]))
            if (m is None or round(lo) <= float(m.group(1)) <= round(hi)) \
                    and str(r[0]).split(" ")[0] == symbol:
                keep.add(tuple(r))
        return keep
    return culls(old) != culls(new)


def _staleness(rec: dict, now: list[dict], symbol: str, window: tuple[float, float],
               waves: list[float] = ()) -> tuple[str, str]:
    """('', '') current | ('STALE', what moved) | ('CODE_DRIFT', which code)."""
    old = {(r["kind"], r["name"]): r.get("digest") for r in rec.get("inputs") or []}
    if not old:
        return "STALE", "the build record lists no inputs"
    cur = {(r["kind"], r["name"]): r.get("digest") for r in now}
    data, code = [], []
    for key in sorted(set(old) | set(cur)):
        if old.get(key) == cur.get(key):
            continue
        kind, name = key
        if kind == "stage_script":
            code.append(name)
        elif kind in ("ledger_scoped", "registry_scoped") and key in old and key in cur:
            moved = _input_moved(kind, name, rec.get("code_commit"), old[key], symbol,
                                 window, waves)
            if moved is not False:
                data.append(name)
        else:
            data.append(name)
    if data:
        return "STALE", "inputs changed since the build: " + ", ".join(data)
    if code:
        return "CODE_DRIFT", "code changed since the build (no input did): " + ", ".join(code)
    return "", ""


def lines_artifact(product: dict) -> Path | None:
    """The per-line artifact of a published product, found from the file its own
    `provenance.copied_to` names (RYA-1229: read what the feed points at, never a folder
    glob). `<stem>_products.csv` -> `<stem>_lines.csv`; a ts-lte deck's un-suffixed products
    file carries its treatment in the lines name (`<stem>_1D-LTE_lines.csv`)."""
    ct = (product.get("provenance") or {}).get("copied_to") or ""
    if not ct.endswith("_products.csv"):
        return None
    stem = rm.ROOT / ct[: -len("_products.csv")]
    for cand in (Path(f"{stem}_lines.csv"), Path(f"{stem}_{product.get('treatment')}_lines.csv")):
        if cand.exists():
            return cand
    return None


def correction_applied(product: dict) -> tuple[int, int] | None:
    """(lines the NLTE/3D correction moved, lines in the aggregate) read from the product's own
    per-line artifact; None when it has none. An ENGINE-A product that moved no line is the
    1D-LTE fit relabelled -- the Ca miss RYA-1349 is here to catch."""
    f = lines_artifact(product)
    if f is None:
        return None
    import pandas as pd
    d = pd.read_csv(f, low_memory=False)
    if "in_aggregate" in d.columns:
        d = d[d["in_aggregate"].astype(str).str.lower().isin(("true", "1"))]
    if "nlte_delta_dex" not in d.columns:
        return 0, len(d)
    delta = pd.to_numeric(d["nlte_delta_dex"], errors="coerce")
    return int((delta.notna() & (delta != 0)).sum()), len(d)


def _pool_waves(record: dict | None) -> list[float]:
    """The air wavelengths of a product's own pool (its budget's indicator ids)."""
    if not record:
        return []
    ids = set(record.get("uncertainty_indicator_ids") or [])
    gf = _git_table(None, "data/linelists/canonical_gf.csv")
    if not ids or gf is None:
        return []
    return [float(w) for w in gf[gf["physical_id"].isin(ids)]["wavelength_air_A"]]


def correction_reach(symbol: str, ion: str, treatment: str, waves: list[float]) -> tuple[int, str]:
    """How many of `waves` the treatment's correction source holds, and that source's name.
    ENGINE-A: the element's registered PySME grid. ENGINE-A-3DNLTE: the cited solar tables
    (data/reference/solar_nlte_tables.yaml) for this element and ion."""
    import pandas as pd
    tol = 0.01
    if treatment == "ENGINE-A":
        from pipeline.element_prepare import engine_a_table
        from pipeline.nlte_corrections import nlte_ions
        table = engine_a_table(symbol)
        if not table.exists():
            return 0, f"no registered ENGINE-A grid for {symbol}"
        stage = rm._ION_ORDER.index(ion) + 1
        if stage not in nlte_ions(symbol):
            return 0, f"{table.name} is not applied to {symbol} {ion} (nlte_ions)"
        g = pd.read_csv(table, comment="#", usecols=["ion", "wave_A"])
        g = g[g["ion"] == stage]
        held = g["wave_A"].astype(float).unique()
        src = table.name
    else:
        import yaml
        reg = rm.ROOT / "data/reference/solar_nlte_tables.yaml"
        tables = [t for t in (yaml.safe_load(reg.read_text()).get("tables") or [])
                  if t["element"] == symbol and ion in t["ions"]]
        if not tables:
            return 0, f"no cited 3D NLTE table for {symbol} {ion} (solar_nlte_tables.yaml)"
        held = []
        for t in tables:
            d = pd.read_csv(rm.ROOT / t["source"])
            d = d[d[t["species_column"]].astype(str).str.strip() == f"{symbol} {ion}"]
            held += d[t["wavelength_column"]].astype(float).tolist()
        src = ", ".join(t["source"] for t in tables)
    return sum(any(abs(w - h) <= tol for h in held) for w in waves), src


def disagreements(products: list[dict], k: float = DISAGREE_K) -> list[dict]:
    """Products the other holdings contradict: same ion, band, selector and treatment on at
    least three holdings, and this one further than k x the others' median sigma from their
    median A. A FLAG for a decision, never a block."""
    import statistics as st
    from collections import defaultdict
    groups = defaultdict(list)
    for p in products:
        if p.get("A") is not None and p.get("sigma_reported"):
            groups[(str(p.get("ion")), p.get("band"), p.get("selector") or "",
                    p.get("treatment"))].append(p)
    out = []
    for (ion, band, sel, t), ps in sorted(groups.items()):
        if len({p.get("holding") for p in ps}) < 3:
            continue
        for p in ps:
            others = [q for q in ps if q.get("holding") != p.get("holding")]
            med = st.median(float(q["A"]) for q in others)
            lim = k * st.median(float(q["sigma_reported"]) for q in others)
            if abs(float(p["A"]) - med) > lim:
                out.append({"ion": ion, "band": band, "selector": sel, "treatment": t,
                            "holding": p.get("holding"), "A": p["A"],
                            "sigma": p["sigma_reported"], "others_median": round(med, 3),
                            "limit_dex": round(lim, 3),
                            "others": {q.get("holding"): q["A"] for q in others}})
    return out


def _selector(pool: str | None) -> str:
    if pool and pool.startswith(rm.SET_POOL_PREFIX):
        return "SET-" + pool[len(rm.SET_POOL_PREFIX):]
    return POOL_SELECTOR.get(pool or "", "")


def _grade(pool: str | None) -> str:
    return "Reference Grade" if pool and pool.startswith(rm.SET_POOL_PREFIX) else GRADE_OF_POOL.get(pool or "", "")


def coverage(star: str, element: str, *, report_dir: Path | None = None) -> dict:
    symbol, _ = rm.split_symbol(element)
    feed = rm.load_feed(star, symbol) or {}
    products = feed.get("products", [])
    index = {}
    for p in products:
        index[(p.get("band"), p.get("holding"), str(p.get("ion")), p.get("selector") or "",
               p.get("treatment"))] = p
    last = {}
    rep = (report_dir or rm.REPORT_DIR) / f"{star}_{element.replace(' ', '')}_latest.json"
    run_doc = json.loads(rep.read_text()) if rep.exists() else {}
    if not run_doc.get("dry_run"):
        # RYA-1349: a dry run's WOULD_RUN / HELD statuses are a plan, not the state of the
        # cells (Mg's committed `_latest` is the 10-27 dry run; Ti's was until RYA-1349).
        for c in run_doc.get("cells", []):
            last[(c["band"], c["holding"], c["ion"], c.get("engine") or "ts-lte", c.get("pool") or "")] = c
    # The budget stage's own verdicts (element_publish): a product built but held by its
    # RYA-587 budget, or skipped as empty, is HELD with that reason -- not NOT_RUN.
    budget_held = {}
    bud = rm.REPORT_DIR / "budget" / f"{star}_{element.replace(' ', '')}" / "report.json"
    if bud.exists():
        for o in json.loads(bud.read_text()).get("products", []):
            why = o.get("skip") or o.get("verdict")
            if why:
                budget_held[(o.get("band"), o.get("holding"), o.get("selector") or "",
                             o.get("key_treatment"))] = f"budget: {why}"
    # A product the publisher's eligibility gate sent to quarantine (RYA-1097) is a verdict on
    # the cell, not missing work: HELD with the gate's reason (RYA-1233 Si CRIRES+ ELGUETA).
    quarantined = {}
    for q in feed.get("quarantine", []):
        quarantined[(q.get("band"), q.get("holding"), str(q.get("ion")), q.get("selector") or "",
                     q.get("treatment"))] = "quarantined: " + (q.get("quarantine_reason") or "")[-200:]
    # A PUBLISHED cell is only done while its inputs are the ones it was built on: the ledger
    # records the inputs hash each build ran against; recompute it now (run()'s own hash).
    ledger = json.loads(rm.LEDGER.read_text()) if rm.LEDGER.exists() else {}
    manifests = rm._manifest_paths(star)

    def stale(d, product: dict) -> tuple[str, str]:
        rec = ledger.get(rm.ledger_key(star, d, rm._band_of(d, star).name_declared))
        if rec is None:
            # Published without an orchestrator build record (a hand run, or before the
            # ledger): nothing proves it was built on today's inputs.
            return "STALE", "no orchestrator build record for this cell"
        try:
            from pipeline.run_descriptor import resolve
            now = rm.input_fingerprints(d, resolve(d, interpreter=None, ispec_dir=None),
                                        manifest_path=manifests.get(d.holding))
        except Exception:
            return "", ""
        if hashlib.sha256(json.dumps(now, sort_keys=True).encode()).hexdigest() == rec.get("inputs_hash"):
            return "", ""
        # RYA-1349: say WHAT moved. A window-wide hash moves on any species' metadata edit
        # and on any code edit; neither is work owed.
        window = (float(d.lo_A) - rm.GF_SCOPE_PAD_A, float(d.hi_A) + rm.GF_SCOPE_PAD_A)
        return _staleness(rec, now, symbol, window, _pool_waves(product))

    quarantine_rec = {(q.get("band"), q.get("holding"), str(q.get("ion")), q.get("selector") or "",
                       q.get("treatment")): q for q in feed.get("quarantine", [])}

    rows = []
    for ion in rm.graded_ions(symbol) or ["I"]:
        # The THREE grades (Ryan): Reference, Codex, Deep. Reference Grade has TWO pools
        # (RYA-1213): the published line set(s) AND the `reference` pool -- every lab-gf line,
        # no depth gate. Both are expected cells.
        for d in rm.expand(star, symbol, ions=[ion], methods=["synthesis"],
                           pools=["set", "reference", "codex", "deep"]):
            sel = _selector(d.pool)
            hold = rm.cell_process_hold(d)
            # The run report stores the pool as the feed token (pool_tier: GRADED / DEEPGRADED /
            # SET-<NAME>), not the orchestrator's pool name.
            run = last.get((d.band, d.holding, ion, d.engine_deck, rm.pool_tier(d) or ""))
            for t in rm.DECK_EMITS[d.engine_deck]:
                p = index.get((d.band, d.holding, ion, sel, t))
                applied = None
                base = index.get((d.band, d.holding, ion, sel, "1D-LTE"))
                base_q = quarantine_rec.get((d.band, d.holding, ion, sel, "1D-LTE"))
                if p is not None:
                    out, why = stale(d, p)
                    out = out or "PUBLISHED"
                    if t in ("ENGINE-A", "ENGINE-A-3DNLTE"):
                        applied = correction_applied(p)
                        if applied is not None and applied[0] == 0:
                            out, why = "NOT_RUN", (f"{t} product corrects 0 of {applied[1]} "
                                                   "lines: the 1D-LTE fit relabelled")
                elif hold:
                    out, why = "HELD", hold
                elif d.pool == "deep":
                    # after the process holds, so a real reason (no deck, no lines) still shows
                    out, why = "HELD", DEEP_POLICY
                elif run is not None and run.get("status") in (rm.FAILED, rm.NOT_READY):
                    # a failure or a cell owed on another host is WORK OWED, never a hold
                    out, why = "NOT_RUN", f"{run['status']}: {run.get('reason', '')}"
                elif run is not None and run.get("status") in (rm.HELD, rm.BLOCKED):
                    out, why = "HELD", f"{run['status']}: {run.get('reason', '')}"
                elif (d.band, d.holding, sel, t) in budget_held:
                    out, why = "HELD", budget_held[(d.band, d.holding, sel, t)]
                elif (d.band, d.holding, ion, sel, t) in quarantined:
                    out, why = "HELD", quarantined[(d.band, d.holding, ion, sel, t)]
                elif t in ("ENGINE-A", "ENGINE-A-3DNLTE") and base is None and base_q is not None:
                    # the pool itself is quarantined (Ca I near-UV Reference, RYA-1311): a
                    # correction of it would be too, so nothing is owed
                    out, why = "HELD", ("the cell's 1D-LTE is quarantined: "
                                        + (base_q.get("quarantine_reason") or "")[-160:])
                elif t in ("ENGINE-A", "ENGINE-A-3DNLTE") and base is not None:
                    # RYA-1349: HELD only when the correction source really holds none of the
                    # pool's lines. A grid that covers them and no product is work owed --
                    # the hold used to be granted on the 1D-LTE product's existence alone.
                    waves = _pool_waves(base)
                    n, src = correction_reach(symbol, ion, t, waves)
                    if n == 0:
                        out, why = "HELD", f"no {t} correction for this pool's lines: {src}"
                    else:
                        out, why = "NOT_RUN", (f"{src} covers {n} of {len(waves)} pool lines "
                                               f"and no {t} product exists")
                elif t in ("ENGINE-A", "ENGINE-A-3DNLTE") and \
                        (d.band, d.holding, sel, "1D-LTE") in budget_held:
                    out, why = "HELD", budget_held[(d.band, d.holding, sel, "1D-LTE")]
                else:
                    out, why = "NOT_RUN", ""
                rows.append({"band": d.band, "holding": d.holding, "ion": ion,
                             "grade": _grade(d.pool), "pool": d.pool, "treatment": t,
                             "outcome": out, "reason": why,
                             "A": p.get("A") if p else None,
                             "sigma": p.get("sigma_reported") if p else None,
                             "n_lines": p.get("n_lines") if p else None,
                             "correction_applied": list(applied) if applied else None})
    counts = {k: sum(r["outcome"] == k for r in rows)
              for k in ("PUBLISHED", "CODE_DRIFT", "HELD", "NOT_RUN", "STALE")}
    return {"star": star, "element": element, "feed_version": feed.get("version"),
            "counts": counts, "done": counts["NOT_RUN"] == 0 and counts["STALE"] == 0,
            "checklist": checklist(rows),
            "disagreements": disagreements(products),
            "asplund": asplund_comparison(products, symbol),
            "asplund_mu1": _asplund_mu1(rm.REPORT_DIR / "budget" / f"{star}_{symbol}"),
            "rows": rows}


#: Asplund, Amarsi & Grevesse 2021, A&A 653, A141, Table 2 (value, sigma), as the site's
#: scripts/asplund2021.py transcribes it.
ASPLUND2021 = {
    'Li': (0.96, 0.06), 'C': (8.46, 0.04), 'N': (7.83, 0.07), 'O': (8.69, 0.04),
    'Na': (6.22, 0.03), 'Mg': (7.55, 0.03), 'Al': (6.43, 0.03), 'Si': (7.51, 0.03),
    'P': (5.41, 0.03), 'S': (7.12, 0.03), 'K': (5.07, 0.03), 'Ca': (6.30, 0.03),
    'Sc': (3.14, 0.04), 'Ti': (4.97, 0.05), 'V': (3.90, 0.08), 'Cr': (5.62, 0.04),
    'Mn': (5.42, 0.06), 'Fe': (7.46, 0.04), 'Co': (4.94, 0.05), 'Ni': (6.20, 0.04),
    'Cu': (4.18, 0.05), 'Zn': (4.56, 0.05), 'Sr': (2.83, 0.06), 'Y': (2.21, 0.05),
    'Zr': (2.59, 0.04), 'Ba': (2.27, 0.05), 'Eu': (0.52, 0.04),
}
#: RYA-1311 (Ca, 2026-11-01), settled for every element: Asplund fits the disk-centre
#: intensity (mu=1); we fit the integrated flux, the geometry of every point-source target.
#: In flux, stronger lines read high by up to ~0.2 dex; at mu=1 the same code landed on A21.
#: An offset inside that is expected and owes no work -- never re-measure it at mu=1
#: (Ryan 2026-11-09). Beyond it, the flux geometry does not explain it: look at the fit.
GEOMETRY_DEX = 0.20


def _asplund_mu1(budget_dir: Path) -> dict | None:
    """The element's Asplund line set measured as Asplund does (scripts/asplund_mu1.py)."""
    f = budget_dir / "asplund_mu1.json"
    return json.loads(f.read_text()) if f.exists() else None


def asplund_comparison(products: list[dict], symbol: str) -> dict | None:
    """Our Asplund line-set products (SET-<El>_AGSS21) per treatment against A21 Table 2,
    with the settled reading of the offset, so no run reopens the flux-vs-disk question."""
    ref = ASPLUND2021.get(symbol)
    sel = f"SET-{symbol.upper()}_AGSS21"
    ps = [p for p in products if p.get("selector") == sel and p.get("A") is not None]
    if ref is None or not ps:
        return None
    by_t: dict[str, list[float]] = {}
    for p in ps:
        by_t.setdefault(p["treatment"], []).append(float(p["A"]))
    rows = []
    for t, vals in by_t.items():
        med = sorted(vals)[len(vals) // 2] if len(vals) % 2 else \
            sum(sorted(vals)[len(vals) // 2 - 1:len(vals) // 2 + 1]) / 2
        off = med - ref[0]
        rows.append({"treatment": t, "n_products": len(vals), "median_A": round(med, 3),
                     "offset_dex": round(off, 3),
                     "reading": ("expected: flux vs Asplund's disk-centre fit (RYA-1311); "
                                 "nothing owed, do not re-measure at mu=1")
                     if abs(off) <= GEOMETRY_DEX else
                     f"beyond the {GEOMETRY_DEX} dex flux geometry explains: look at the fit"})
    return {"asplund": ref[0], "asplund_sigma": ref[1], "selector": sel, "rows": rows}


def checklist(rows: list[dict]) -> list[dict]:
    """Ryan's questions, one line each: per treatment and per instrument family, how many
    cells are done (PUBLISHED / CODE_DRIFT), held with a reason, and owed (NOT_RUN / STALE)."""
    def tally(name, sel):
        rs = [r for r in rows if sel(r)]
        done = sum(r["outcome"] in ("PUBLISHED", "CODE_DRIFT") for r in rs)
        owed = sum(r["outcome"] in ("NOT_RUN", "STALE") for r in rs)
        held = sum(r["outcome"] == "HELD" for r in rs)
        # a published NLTE/3D product whose per-line artifact is not kept cannot show it
        # corrected anything (Al, Si: their stage predates the lines files)
        unverified = sum(r["outcome"] in ("PUBLISHED", "CODE_DRIFT") and not r.get("correction_applied")
                         and r["treatment"] in ("ENGINE-A", "ENGINE-A-3DNLTE") for r in rs)
        return {"item": name, "done": done, "held": held, "owed": owed,
                "unverified": unverified, "ok": owed == 0 and done > 0}
    items = [tally(f"treatment {t}", lambda r, t=t: r["treatment"] == t)
             for t in dict.fromkeys(r["treatment"] for r in rows)]
    items += [tally(f"instrument {fam}", lambda r, fam=fam: r["holding"].startswith(prefix))
              for fam, prefix in (("CRIRES+", "solar_crires_plus"), ("HARPS", "solar_harps"),
                                  ("Kitt Peak", "solar_kpno"), ("IAG", "solar_iag"))]
    return [i for i in items if i["done"] + i["held"] + i["owed"]]


def render(cov: dict) -> str:
    lines = [f"COVERAGE  {cov['star']} {cov['element']}  feed v{cov['feed_version']}  "
             + "  ".join(f"{k} {v}" for k, v in cov["counts"].items())
             + ("  -> DONE" if cov["done"] else "  -> NOT DONE (NOT_RUN / STALE cells owe work)")]
    lines.append("  CHECKLIST  [x] done, nothing owed   [-] every cell held, with its reason   [ ] work owed")
    for i in cov.get("checklist", []):
        mark = "x" if i["ok"] else ("-" if not i["owed"] else " ")   # '-': all held, each with a reason
        lines.append(f"    [{mark}] {i['item']:40} done {i['done']:4}  "
                     f"held {i['held']:4}  owed {i['owed']:4}"
                     + (f"  ({i['unverified']} with no per-line evidence of the correction)"
                        if i.get("unverified") else ""))
    flags = cov.get("disagreements", [])
    lines.append(f"  HOLDINGS THAT DISAGREE ({len(flags)}; flagged for a decision, not blocked)")
    for f in flags:
        lines.append(f"    {f['ion']:2} {f['band']:11} {f['selector'][:14]:14} {f['treatment'][:20]:20} "
                     f"{f['holding'][:30]:30} A={f['A']:.3f} vs others {f['others_median']:.3f} "
                     f"(limit {f['limit_dex']:.3f})")
    asp = cov.get("asplund")
    if asp:
        lines.append(f"  ASPLUND 2021  {asp['asplund']:.2f} +/- {asp['asplund_sigma']:.2f} "
                     f"against our {asp['selector']} products (median over holdings and bands)")
        for a in asp["rows"]:
            lines.append(f"    {a['treatment'][:32]:32} {a['median_A']:.3f} ({a['offset_dex']:+.3f}, "
                         f"n={a['n_products']})  {a['reading']}")
        mu = cov.get("asplund_mu1")
        if mu and mu.get("mean_A_mu1") is not None:
            lines.append(f"    [x] measured as Asplund does (mu=1, Liege): {mu['mean_A_mu1']:.3f} "
                         f"({mu['mean_A_mu1'] - asp['asplund']:+.3f}) on {mu['n_lines']} lines; "
                         f"per line vs A21 1D MARCS {mu['mean_minus_A21_marcs_lte']:+.3f} "
                         f"(n={mu['n_paired']})" if mu.get("mean_minus_A21_marcs_lte") is not None
                         else f"    [x] measured as Asplund does (mu=1, Liege): {mu['mean_A_mu1']:.3f}")
        else:
            lines.append("    [ ] measured as Asplund does (mu=1, Liege): owed -- "
                         "run_pipeline.py --element <El> --asplund-mu1")
    for r in sorted(cov["rows"], key=lambda r: (BAND_ORDER.index(r["band"]) if r["band"] in BAND_ORDER else 99,
                                                r["holding"], r["ion"], r["grade"], r["treatment"])):
        val = (f"{r['A']:.3f} +/- {r['sigma']:.3f} (n={r['n_lines']})" if r["A"] is not None
               and r["sigma"] is not None else "")
        if r.get("correction_applied"):
            val += " corrects {}/{} lines".format(*r["correction_applied"])
        lines.append(f"  {r['outcome']:10} {r['band']:11} {r['holding'][:30]:30} {r['ion']:2} "
                     f"{r['grade'][:22]:22} {r['treatment']:8} {val or r['reason'][:70]}")
    return "\n".join(lines)

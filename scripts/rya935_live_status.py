#!/usr/bin/env python3
"""Derive the live tracker's status from products on disk (RYA-935).

The tracker's JSON was HAND-TYPED. Two things follow from that, and both had
already happened by the time this was written:

  * it violated RYA-686 -- a result artifact landed with no generating harness,
    which is the RYA-559 hole that convention exists to close, and CI said so;
  * it went stale silently. The committed copy listed RYA-931 as BACKLOG on the
    day RYA-931 merged, and showed HARPS as "correction owed" after the
    corrected holding existed. A dashboard that is typed is a dashboard that
    lies eventually.

So this derives every cell it can from artifacts that already exist, and is
element-agnostic on purpose: Fe appears the moment its products land, with no
edit here.

What it deliberately does NOT emit: the ticket pipeline. That is Linear state,
not repository state, and a copy of it in a committed file is guaranteed to
drift -- which is exactly the failure above. The dashboard shows what the repo
can prove.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import product_eligibility as pe  # noqa: E402  RYA-1097
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

#: `<El><Ion>_<lo>_<hi>_<instrument>_<holding>_<HANDLER>_products.csv`. The holding
#: is in the stem because RYA-933/934 put it there -- before that, two holdings of
#: one instrument wrote the same filename and the second overwrote the first.
STEM = re.compile(r"^(?P<el>[A-Z][a-z]?)(?P<ion>I+|IV|VI*)_(?P<lo>\d+)_(?P<hi>\d+)_"
                  r"(?P<rest>.+?)_(?P<handler>PROFILEFIT|SYNTH)"
                  # 🔴 RYA-1031: ONE OR MORE selector segments. This line ended at the
                  # handler, so a stem carrying ANY selector -- _GRADED, _DEEPGRADED,
                  # _FROMEW, _LOCALRENORM, _ENGINE-A -- did not match, `parse_stem`
                  # returned None, and the product was silently absent from the page.
                  # That is why HARPS and IAG have never appeared and why every graded
                  # arm was invisible: not a missing run, an unparseable filename.
                  r"(?P<selector>(?:_[A-Za-z0-9<>-]+)*)_products\.csv$")


def species_display(species: str) -> str:
    """`FeI` -> `Fe`, `ScII` -> `Sc II`. The ion is shown only when there IS one.

    🔴 RYA-1031 (Ryan): a neutral IS the base element and carries no indicator. Writing
    every neutral as `FeI`/`CI`/`VI` put a roman numeral on 25 of 26 roster entries where
    it said nothing, and on one it said something WRONG -- `VI` is vanadium-neutral and
    reads as roman six. Only an actual ion (II and up) earns the suffix.

    Underlying keys are NOT renamed: `reference` and the product rows stay keyed by
    species, because that is how they join. This is the rendered name only.
    """
    m = re.fullmatch(r"([A-Z][a-z]?)(I+|IV|VI*)", species or "")
    if not m:
        return species
    element, ion = m.groups()
    return element if ion == "I" else f"{element} {ion}"


def _source_path(path: Path) -> str:
    """Repo-relative when the product is inside the repo, absolute when it is not."""
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def parse_stem(name: str, instruments: set[str], holdings: set[str]) -> dict | None:
    m = STEM.match(name)
    if not m:
        return None
    rest = m.group("rest")
    instrument = next((i for i in sorted(instruments, key=len, reverse=True)
                       if rest == i or rest.startswith(i + "_")), None)
    if instrument is None:
        return None
    tail = rest[len(instrument):].lstrip("_")
    # No holding in the stem = a product from BEFORE RYA-933/934, when the stem
    # keyed on instrument alone. Say so; do not guess. Every product committed
    # before that change is attributable to an instrument and NOT to a holding,
    # and an instrument can serve several -- including a corrected and an
    # uncorrected one. Inferring which would be exactly the collapse the stem
    # change was made to prevent.
    if tail and tail in holdings:
        holding, source = tail, "filename"
    elif tail:
        holding, source = tail, "filename (unregistered holding)"
    else:
        holding, source = None, "absent -- product predates RYA-933/934"
    return {"element": m.group("el"), "ion": m.group("ion"),
            "lo_A": float(m.group("lo")), "hi_A": float(m.group("hi")),
            "instrument": instrument, "holding": holding,
            "holding_source": source, "handler": m.group("handler"),
            # Carried, not discarded: two runs differing only in selector are two
            # different products (RYA-984) and must not collapse into one row.
            "selector": (m.group("selector") or "").lstrip("_") or "default"}


#: Sub-paths that mark a DIAGNOSTIC variant rather than a headline product --
#: RYA-877's before/after control pair, RYA-847's gated sweep. They are real
#: products and must not be deleted from the page; they are also not what the
#: matrix is asking about, so they are labelled and hidden behind a toggle.
VARIANT_MARKERS = ("control", "gated")

#: Extra results roots handed in via `--products-root`, so `run_context` can anchor a
#: product that lives outside this checkout instead of raising on it (RYA-1031).
_EXTRA_RESULT_ROOTS: list[Path] = []


#: The dates the telluric-corrected holdings first existed in the repository.
#: A product committed before its instrument had a corrected holding CANNOT have
#: used one -- this is provenance, not inference, and it is the first question
#: anyone asks of a number on this page.
CORRECTED_HOLDING_BORN = {
    "harps": "2026-08-20",              # RYA-931, commit 4d8abf8
    "kpno_solar_atlas": "2026-08-21",   # RYA-940, commit c0465b1
}


def telluric_state_of(row: dict, committed: str | None) -> dict:
    """Was this product made before or after its instrument had a corrected holding?

    Where the row names a holding, the holding answers it outright. Where it does
    not -- every product predating RYA-933/934 -- the date still answers it: the
    corrected holdings did not exist, so nothing could have used them.
    """
    if row.get("holding"):
        return {"telluric_basis": "named holding", "telluric_epoch": None}
    born = CORRECTED_HOLDING_BORN.get(row["instrument"])
    if born and committed and committed < born:
        return {"telluric_basis": "PRE-correction (provenance)",
                "telluric_epoch": f"committed {committed}; {row['instrument']} had no "
                                  f"corrected holding until {born}"}
    return {"telluric_basis": "unknown", "telluric_epoch": None}


def run_context(path: Path) -> dict:
    """WHICH RUN produced this row. Part of a product's identity, not decoration.

    Six identities appear more than once across ticket output directories, and
    they are NOT duplicates: rya877/control/before gives Fe II 7.568 where
    rya877/control/after gives 7.542. That pair is the whole point of a control.
    Deduplicating on (species, holding, band, engine) would silently keep one and
    discard the other -- picking a winner between two numbers that were produced
    to be compared.
    """
    import subprocess
    # RYA-1031: a product from another `--products-root` is not under this repo, so
    # anchor on whichever results root actually contains it. This said
    # relative_to(ROOT/"data"/"results") unconditionally and raised, which is the
    # second reason the multi-root flag could never be used across checkouts.
    for _base in (ROOT / "data" / "results", *_EXTRA_RESULT_ROOTS):
        try:
            rel = path.relative_to(_base).parent
            break
        except ValueError:
            continue
    else:
        rel = Path(path.parent.name)
    parts = rel.parts
    try:
        committed = subprocess.run(
            ["git", "log", "--format=%ad", "--date=short", "-1", "--", str(path)],
            cwd=ROOT, capture_output=True, text=True, timeout=20).stdout.strip() or None
    except Exception:                                           # noqa: BLE001
        committed = None
    return {"run_context": str(rel) if parts else "(root)",
            "ticket_dir": parts[0] if parts else None,
            "committed": committed,
            "is_variant": any(m in parts for m in VARIANT_MARKERS)}


def display_name(row) -> str:
    """The physics-axis name, DERIVED from the stored axes (RYA-906).

    Not `treatment`. "ENGINE-A" / "ENGINE-B" are letters, not physics: they say
    nothing about route, scale or model, and the same letter means different
    things on different rows. RYA-906 stored the five axes precisely so the name
    could be derived, and the tracker was still showing the legacy labels.

    Route comes from the HANDLER where the row carries one -- never from the
    label. RYA-906 measured this over 2153 committed rows: `1D-LTE` and
    `1D-LTE-LABGF` each pair with BOTH ProfileFitHandler and SynthesisHandler, so
    on those labels the legacy string is not merely lossy, it is FALSE.
    """
    from pipeline.treatment_axes import Axes
    stored = {k: row.get(k) for k in ("route", "scale", "model", "atmos", "gf")}
    if all(v is not None and not pd.isna(v) for v in stored.values()):
        return Axes(route=str(stored["route"]), scale=str(stored["scale"]),
                    model=str(stored["model"]), atmos=str(stored["atmos"]),
                    gf=str(stored["gf"]),
                    route_basis=str(row.get("route_basis") or "stored")).display
    # No stored axes = a row that predates RYA-906. Derive from the legacy label
    # plus whatever route evidence the row carries, and fail loudly on an unknown
    # label rather than defaulting -- a silent default is how RYA-869 published
    # four wrong systematics.
    from pipeline.treatment_axes import display_for, UnknownTreatment
    try:
        return display_for(str(row.get("treatment")), handler=row.get("handler"))
    except (UnknownTreatment, ValueError):
        return f"{row.get('treatment')} (unresolved axes)"


#: `gf rung: gf rung N (term): reason` -- the reason text states the graded count in one
#: of two shapes, and both are parsed here rather than recomputed. Recomputing would mean
#: re-grading every line against the line list inside the tracker, which is a SECOND
#: implementation of membership (the RYA-845 two-homes shape); the budget file is the
#: artifact the product was actually CHARGED on, so it is the honest source.
_RUNG_MIXED = re.compile(r"MIXED POOL:\s*(\d+)\s+of\s+(\d+)\s")
_RUNG_ALL = re.compile(r"every one of the\s+(\d+)\s")
_RUNG_HEAD = re.compile(r"gf rung:\s*gf rung\s*(\d+)\s*\(([^)]*)\)")


def graded_counts(products_csv: Path) -> dict:
    """(n_graded, n_pool, rung) per treatment, read from the sibling *_budgets.txt.

    Returns {} when there is no budget file -- an older artifact predates the gf rung and
    must read as UNKNOWN, never as zero. Zero graded lines is a real, different statement
    from "this product was written before we recorded the rung".
    """
    b = products_csv.with_name(products_csv.name.replace("_products.csv", "_budgets.txt"))
    if not b.exists():
        return {}
    out, treatment = {}, None
    for line in b.read_text(errors="replace").splitlines():
        t = line.strip()
        # Budget blocks open with the cell header, e.g. "Fe . VIS . n=6"; the gf rung line
        # belongs to the block above it. Track the most recent non-indented header.
        if t and not line.startswith(" ") and "gf rung" not in t:
            treatment = t
        if "gf rung:" not in t:
            continue
        head = _RUNG_HEAD.search(t)
        rung = int(head.group(1)) if head else None
        m = _RUNG_MIXED.search(t)
        if m:
            n_graded, n_pool = int(m.group(1)), int(m.group(2))
        else:
            m2 = _RUNG_ALL.search(t)
            if not m2:
                continue
            n_graded = n_pool = int(m2.group(1))
        # KEYED BY POOL SIZE, not by treatment name. The budget block header is
        # "Fe · VIS · n=148" -- it carries the cell and the pool size but NOT the
        # treatment string the products table uses, and the two vocabularies do not
        # match (RYA-906: `display` is derived, `treatment` is the legacy label). The
        # header's n IS the pool the budget was charged on, so it joins to the row's
        # n_lines exactly.
        out[int(n_pool)] = {"n_graded": n_graded, "n_pool": n_pool, "gf_rung": rung,
                            "budget_cell": treatment}
    return out


def collect_products(roots: list[Path], instruments: set[str],
                     holdings: set[str]) -> list[dict]:
    rows: list[dict] = []
    for root in roots:
        if not root.exists():
            continue
        for path in sorted(root.rglob("*_products.csv")):
            meta = parse_stem(path.name, instruments, holdings)
            if meta is None:
                continue
            try:
                frame = pd.read_csv(path)
            except Exception:                                   # noqa: BLE001
                continue
            ctx = run_context(path)
            graded = graded_counts(path)
            for _, r in frame.iterrows():
                rows.append({
                    **meta,
                    "band": r.get("band"), "treatment": r.get("treatment"),
                    "display": display_name(r),
                    "A": None if pd.isna(r.get("A")) else float(r["A"]),
                    "sigma_stat": None if pd.isna(r.get("stat_dex")) else float(r["stat_dex"]),
                    "sigma_syst": None if pd.isna(r.get("syst_dex")) else float(r["syst_dex"]),
                    "n_lines": None if pd.isna(r.get("n_lines")) else int(r["n_lines"]),
                    # RYA-946 — how many of THIS engine's lines carry a primary-lab gf.
                    # None means the artifact predates the gf rung, which is not 0.
                    **((graded.get(int(r["n_lines"]))
                        if not pd.isna(r.get("n_lines")) else None)
                       or {"n_graded": None, "n_pool": None, "gf_rung": None}),
                    # RYA-1031: `--products-root` takes MULTIPLE roots, but this said
                    # relative_to(ROOT) and raised on any product outside the repo --
                    # so the multi-root flag could never actually be used with one.
                    # Relative where it can be, absolute where it cannot; a path that
                    # cannot be expressed is still a path, and dropping the row would
                    # hide a product rather than locate it.
                    # 🔴 RYA-1031: WHEN THIS RUN ACTUALLY WROTE. `committed` is the git
                    # date and is EMPTY for every band product, because band_products/ is
                    # gitignored (RYA-469) -- so the page could not distinguish a number
                    # measured an hour ago from one measured two days ago. The file mtime
                    # is the only honest answer available, and a stale number that LOOKS
                    # current is the failure this whole page exists to prevent.
                    "measured_at": datetime.fromtimestamp(
                        path.stat().st_mtime, timezone.utc).isoformat(timespec="minutes"),
                    "source": _source_path(path),
                    **ctx,
                })
    return rows


def collect_feed_products(store: Path, instruments: set[str],
                         holdings: set[str]) -> tuple[list[dict], list[dict]]:
    """The LIVE products, read from the element feed -- the same source the public site reads.

    🔴 THE TRACKER USED TO CLIMB IN THE WINDOW. RYA-1092 gates `products[]` in
    `data/products/<star>/<El>.json`, and RYA-1091's website reads that feed -- but this
    page globbed `data/results/` for `*_products.csv` and built its own list. So a product
    the gate had just quarantined (syst-incomplete, pre-continuum-fix, an uninterpretable
    error bar) still rendered here as live. The gate protected the public front door while
    the internal page came in another way, and an INTERNAL page that disagrees with the
    public one is worse than either being wrong alone: it is the surface we check the
    public one against.

    Reading the feed also makes the two move together. A product is live here the moment it
    is published and gated, and gone the moment it is withdrawn -- no second scan to fall
    out of step, and no path where a file on disk is enough to appear.

    ⚠️ THE RUN CONTEXT IS ENRICHED FROM THE COMMITTED ARTIFACT, NOT RE-DERIVED. The feed
    carries the published value and its provenance; the per-run detail this page shows
    (band edges, handler, gf rung, when it was measured) lives in the artifact. It is read
    through `provenance.copied_to`, which RYA-1096 made repo-relative for exactly this kind
    of cross-machine read. A product whose artifact is absent still appears -- with the
    enrichment fields None and `enrichment` saying why -- because dropping it would hide a
    live product behind a missing sidecar (RYA-833).

    Returns (live, withheld). `withheld` carries the feed's own quarantine/archive records
    so a diagnostics view can show them EXPLICITLY, labelled, never mixed into the live set.
    """
    live: list[dict] = []
    withheld: list[dict] = []
    for feed in sorted(Path(store).glob("*/[A-Z]*.json")):
        doc = json.loads(feed.read_text())
        star = str(doc.get("star") or feed.parent.name)
        for pool in ("products", "quarantine", "superseded", "archive"):
            for rec in doc.get(pool) or []:
                row = _feed_row(rec, star=star, feed=feed, pool=pool,
                                instruments=instruments, holdings=holdings)
                (live if pool == "products" else withheld).append(row)
    return live, withheld


def _feed_row(rec: dict, *, star: str, feed: Path, pool: str,
              instruments: set[str], holdings: set[str]) -> dict:
    prov = rec.get("provenance") or {}
    copied = prov.get("copied_to")
    artifact = (ROOT / copied) if copied and not Path(copied).is_absolute() else None
    row = {
        "star": star, "feed": str(feed.relative_to(ROOT)), "pool": pool,
        "element": rec.get("element"), "ion": rec.get("ion"),
        "band": rec.get("band"), "instrument": rec.get("instrument"),
        "holding": rec.get("holding"),
        # 🔴 `holding_source` is 'feed', not 'filename'. The old scanner PARSED the
        # holding out of the stem and had to say so, because a stem from before
        # RYA-933/934 carries none. A published record states it as a field, so the
        # provenance of the provenance changes and the page must not claim otherwise.
        "holding_source": "feed record",
        "tier": rec.get("tier"), "selector": rec.get("selector"),
        "route": rec.get("route"), "treatment": rec.get("treatment"),
        "display": rec.get("display"),
        "A": rec.get("A"), "sigma_stat": rec.get("sigma_stat"),
        "sigma_syst": rec.get("sigma_syst"),
        # RYA-1095: what `sigma_stat` MEANS, carried per product. A bar whose
        # construction is unstated is not comparable with the one beside it.
        "stat_basis": rec.get("stat_basis"),
        "n_lines": rec.get("n_lines"), "n_excluded": rec.get("n_excluded"),
        "dominant_term": rec.get("dominant_term"),
        "identity": pe.key_of(rec),
        "measured_at": prov.get("artifact_mtime"),
        "ingested_at": prov.get("ingested_at"),
        "source": copied or prov.get("path"),
    }
    if pool != "products":
        row["withdrawn_reason"] = (rec.get("quarantine_reason")
                                   or rec.get("superseded_reason"))
        row["withdrawn_codes"] = rec.get("quarantine_codes")
        row["withdrawn_at"] = rec.get("quarantined_at") or rec.get("superseded_at")
    if artifact is not None and artifact.exists():
        g = graded_counts(artifact)
        n = rec.get("n_lines")
        row.update((g.get(int(n)) if n is not None else None)
                   or {"n_graded": None, "n_pool": None, "gf_rung": None})
        row.update(run_context(artifact))
        row["enrichment"] = "from the committed artifact"
    else:
        row.update({"n_graded": None, "n_pool": None, "gf_rung": None})
        row["enrichment"] = (
            f"artifact not readable from this checkout ({copied!r}); the published value "
            f"stands, the per-run detail is unavailable here" if copied else
            "the record carries no copied_to, so no committed artifact to enrich from")
    return row


#: The order fields are tried when disambiguating a shared display label. Holding first
#: because it is what actually differs between most collisions and what a reader needs;
#: `treatment` last because a treatment that renders the same string as another is a
#: LABELLING defect one level up (`publish_product._DISPLAY` maps both `1D-LTE` and
#: `ENGINE-B` to "Synth · 1D-LTE"), and naming it here is a workaround, not the fix.
_DISAMBIGUATION_ORDER = ("holding", "band", "tier", "selector", "route", "instrument",
                         "treatment")


def _label_for(row: dict, group: list, candidates: list) -> str:
    """The display label plus the FEWEST fields that make this row unique in its group.

    🔴 NOT "every field that differs somewhere in the group". That was the first version
    and it produced `Synth · 1D-LTE · holding=... · instrument=... · band=... · route=... ·
    treatment=...` -- five clauses, most of which distinguish some OTHER pair. A label a
    reader will not read is not a disambiguation.

    Greedy in `_DISAMBIGUATION_ORDER`, stopping the moment the row is unique, so a
    collision that only differs by holding says only the holding.
    """
    used: list = []
    for f in candidates:
        used.append(f)
        key = tuple(str(row.get(x)) for x in used)
        if sum(1 for o in group
               if tuple(str(o.get(x)) for x in used) == key) == 1:
            break
    return " · ".join([str(row.get("display"))]
                      + [f"{f}={row.get(f)}" for f in used])


def identity_report(rows: list[dict]) -> dict:
    """True duplicates (a bug) vs distinct products that share a display label (a labelling
    problem). The two are not the same finding and must not be reported as one.

    A TRUE DUPLICATE is two live records with the same identity key -- the nine fields
    `publish_product.KEY_FIELDS` names, `route` and `selector` included. Those two are NOT
    decoration: IAG VIS GRADED ENGINE-A is 7.481 on n=4 by profile fit and 7.484 on n=37 by
    synthesis, so a key without `route` collides two different measurements of two
    different pools. RYA-1097's spec lists the key without `route`; it is kept here because
    the store's own guard caught that collision on its first backfill.

    A SHARED LABEL is several distinct products rendering the same string. The remedy is to
    surface the field that distinguishes them, which `disambiguated` gives per row, so a
    reader never sees two identical-looking lines that are different measurements.
    """
    by_identity: dict = {}
    by_label: dict = {}
    for r in rows:
        by_identity.setdefault(r["identity"], []).append(r)
        by_label.setdefault(str(r.get("display")), []).append(r)

    true_dups = {k: [_brief(x) for x in v] for k, v in by_identity.items() if len(v) > 1}
    shared: dict = {}
    for label, group in by_label.items():
        if len(group) < 2:
            continue
        # What actually differs across the group -- reported, not guessed at.
        fields = [f for f in _DISAMBIGUATION_ORDER
                  if len({str(r.get(f)) for r in group}) > 1]
        shared[label] = {"n": len(group), "distinguished_by": fields,
                         "rows": [_brief(x) for x in group]}
        for r in group:
            r["disambiguated"] = _label_for(r, group, fields)
    for r in rows:
        r.setdefault("disambiguated", r.get("display"))
    return {
        "key_fields": list(pe.KEY_FIELDS),
        "n_live": len(rows),
        "n_unique_identities": len(by_identity),
        "true_duplicates": true_dups,
        "shared_display_labels": shared,
        "note": ("A true duplicate is a BUG -- one identity, one live record. A shared "
                 "label is a LABELLING problem and is resolved by `disambiguated`, which "
                 "appends the fields that actually differ."),
    }


def _brief(r: dict) -> dict:
    return {k: r.get(k) for k in ("holding", "instrument", "band", "tier", "selector",
                                  "route", "treatment", "A", "n_lines", "source")}


def collect_instruments() -> list[dict]:
    """Coverage and telluric state, read from the registries -- never retyped."""
    import measure_band_ew as M
    from pipeline.telluric_policy import applied_state
    out = []
    for instrument, specs in M._INSTRUMENT_HOLDINGS.items():
        for spec in specs:
            try:
                telluric = applied_state(spec.holding_id)
            except KeyError:
                telluric = "unregistered"
            out.append({
                "instrument": instrument, "holding": spec.holding_id,
                "telluric_applied": telluric,
                "pre_normalised": spec.pre_normalised,
                "coverage_A": list(spec.span_A) if spec.span_A else None,
                "coverage_note": ("declared" if spec.span_A else
                                  "discovered by the reader; not declarable as one "
                                  "interval"),
            })
    return out


def collect_reference(root: Path) -> dict:
    """Literature anchor per species, from the FROZEN gold reference.

    Read, not typed. The previous hand-written page carried Al 6.43 as a literal;
    it is right, but a literal cannot follow the reference when it is re-frozen,
    and this project's whole discipline is that a value cites its source. The
    pointer file `data/reference/solar/CURRENT` names the live version.
    """
    current = (root / "data" / "reference" / "solar" / "CURRENT")
    version = current.read_text().strip() if current.exists() else "v5"
    table = root / "data" / "reference" / "solar" / f"solar_abundances_{version}.csv"
    if not table.exists():
        return {}
    frame = pd.read_csv(table, comment="#")
    out = {}
    for _, r in frame.iterrows():
        if pd.isna(r.get("asplund2021")):
            continue
        out[f"{r['element']}{r['ion']}"] = {
            "asplund2021": float(r["asplund2021"]),
            "codex_A_X": None if pd.isna(r.get("A_X")) else float(r["A_X"]),
            "verdict": str(r.get("verdict")),
            "source": f"data/reference/solar/solar_abundances_{version}.csv",
            "sigma_external": None, "band": None, "best_external": None,
            "scale": None, "deviate_beyond": None,
        }

    # The gold table carries the literature VALUE but no uncertainty. litscan does,
    # and it is the ratified comparator: best-external +/- sigma_external, with the
    # source named. Take the band from there wherever an element has a litscan.
    #
    # NOTE the band is AGREEMENT WITH THE LITERATURE, not a pass/fail gate. litscan's
    # own basis text warns against conflating it with the FE_GATE policy window
    # ([7.41, 7.51], RYA-166) -- they answer different questions.
    try:
        from pipeline import litscan
        for element in litscan.available_elements():
            rng = litscan.literature_range(element)
            if rng is None or rng.sigma_external is None:
                continue
            for key in [k for k in out if k.startswith(element)
                        and k[len(element):].strip("IV") == ""]:
                out[key].update({
                    "asplund2021": rng.central,
                    "sigma_external": rng.sigma_external,
                    "band": [rng.min, rng.max],
                    "deviate_beyond": rng.deviate_beyond,
                    "best_external": rng.best_external,
                    "scale": rng.scale,
                    "band_meaning": ("agreement with the literature (best external +/- "
                                     "sigma_ext) -- NOT a pass/fail gate"),
                    "source": f"pipeline/litscan.py :: {element}.yaml",
                })
    except Exception:                                           # noqa: BLE001
        pass

    # A SECOND comparator, kept separate rather than averaged: Lodders, Bergemann &
    # Palme 2025 Table 6. Its PRESENT column is the photospheric-era value; the
    # proto-solar column runs ~0.09 dex higher and quoting it against a photospheric
    # measurement would manufacture a discrepancy.
    #
    # Comparators are a LIST, and each may be scoped to particular bands. A source
    # whose determination only covers the infrared must not be drawn across the
    # optical as if it applied there.
    lodders = root / "data" / "reference" / "solar" / "lodders2025_table6.csv"
    if lodders.exists():
        frame = pd.read_csv(lodders, comment="#")
        by_element = {r["element"]: r for _, r in frame.iterrows()}
        for key, entry in out.items():
            element = "".join(c for c in key if not c.isupper() or c == key[0])
            element = key.rstrip("IV") or key
            row = by_element.get(element)
            if row is None:
                continue
            entry.setdefault("comparators", []).append({
                "name": "Lodders+ 2025", "value": float(row["A_present"]),
                "sigma": float(row["sigma_present"]), "colour": "gold",
                "bands": None,          # applies everywhere
                "note": "Table 6, present-day Sun (proto-solar is "
                        f"{float(row['A_protosolar']):.2f}, ~0.09 dex higher)",
                "source": "data/reference/solar/lodders2025_table6.csv",
            })
    return out


def why_no_product(element: str, ion: str, holding: str, instrument: str,
                   lo: float, hi: float) -> str:
    """Why this cell is empty -- from the SAME resolver that plans the runs.

    An empty cell with no explanation is the RYA-833 shape: "we do not hold this"
    becomes indistinguishable from "nobody looked". Every blank on this grid has
    to say which it is.
    """
    from pipeline.run_descriptor import RunDescriptor, resolve
    d = RunDescriptor(element, ion, instrument, holding, lo, hi)
    # Interpreter and engine dir are supplied so that only REAL blockers -- coverage,
    # wiring, the telluric gate -- surface here. Whether a run host happens to have
    # the right numpy is not a fact about the science and does not belong on the grid.
    r = resolve(d, interpreter="(host)", ispec_dir="(host)")
    if r.blocked_reason:
        return r.blocked_reason
    return "no run yet"


def collect_graded(root: Path) -> list[dict]:
    """The GRADED lab-gf cells -- RYA-850's primary reported value.

    These were missing from the tracker entirely, which is why the page appeared to
    show error bars that had grown. They had not: the page was showing only the
    MIXED pools, which carry the blanket 0.17 dex ungraded-gf placeholder. A fully
    graded cell carries its own CITED pool sigma instead, and the systematic drops
    to 0.061-0.113.

    They live in `rya850_summary.json`, not in a `*_products.csv`, so the product
    glob never saw them.

    RYA-851's reporting contract: GRADED is primary, UNGRADED is secondary, and the
    ungraded value is never the headline. Both are shown -- more lines buy a wider
    gf floor, which is a trade, not a defect.
    """
    summary = root / "data" / "results" / "rya850" / "rya850_summary.json"
    if not summary.exists():
        return []
    doc = json.loads(summary.read_text())
    out = []
    for cell in doc.get("graded_cells", []):
        out.append({
            "element": "Fe", "ion": cell.get("ion", "I"), "band": cell["band"],
            "engine": cell.get("engine"), "A": cell["A"], "n_lines": cell["n_lines"],
            "sigma_stat": cell["stat_dex"], "sigma_syst": cell["syst_dex"],
            "total_dex": cell.get("total_dex"),
            "ungraded_total_dex": cell.get("ungraded_total_dex"),
            "graded_beats_ungraded": cell.get("graded_beats_ungraded"),
            "gf_term": doc.get("gf_term_published"),
            "role": "PRIMARY (graded lab-gf)",
            "source": "data/results/rya850/rya850_summary.json",
        })
    return out


def collect_telluric(audit_root: Path) -> list[dict]:
    """Before/after residuals from the correction tickets' own evidence."""
    out = []
    for manifest in sorted(audit_root.glob("rya940_kp1984_telluric/*/fit_manifest.json")):
        d = json.loads(manifest.read_text())
        c = d.get("correction")
        if not c:
            out.append({"metric": "pct_below_0.5", "window": f"{int(d['band_A'][0])}-{int(d['band_A'][1])} A",
                        "product": "solar_kpno_molecfit_corrected",
                        "before_pct_below_0.5": None, "after_pct_below_0.5": None,
                        "note": "no admissible fit; band NOT corrected"})
            continue
        out.append({
            "metric": "pct_below_0.5",
            "window": f"{int(d['band_A'][0])}-{int(d['band_A'][1])} A",
            "product": "solar_kpno_molecfit_corrected",
            "before_pct_below_0.5": c["before"]["pct_below_0.5"],
            "after_pct_below_0.5": c["after"]["pct_below_0.5"],
            "externally_validated": d.get("externally_validated"),
        })
    for row in _stellar_crires(audit_root):
        out.append(row)
    for row in _harps_state(audit_root):
        out.append(row)
    verification = audit_root / "rya931_molecfit_runtime" / "verification.json"
    if verification.exists():
        v = json.loads(verification.read_text())["o2b_gate"]
        out.append({"metric": "pct_below_0.5", "window": "6867-6884 A", "product": "solar_harps_molecfit_corrected",
                    "before_pct_below_0.5": v["o2b_before"]["pct_below_0.5"],
                    "after_pct_below_0.5": v["o2b_after"]["pct_below_0.5"],
                    "externally_validated": True})
    return out


#: The STELLAR telluric legs score with a DIFFERENT metric from the solar ones, and the
#: two must never share a column. RYA-940/931 report `pct_below_0.5` in a registered
#: telluric band; RYA-963/973 report the D1 residual — the median |1 - continuum| at
#: pixels molecfit calls telluric-DOMINATED and the star's own line list calls CLEAN.
#: They answer different questions and are not comparable as numbers, so each row names
#: the metric it carries (RYA-873: report a value under its DERIVED name).
#: Which corrected holding each ticket's evidence directory belongs to. Keyed here
#: because the manifest records products, and the HOLDING is what the registry and the
#: reader think in.
_AUDIT_DIR_HOLDING = {
    'rya963_crires_telluric': 'alpha_cen_a_crires_plus_molecfit',
    'rya973_crires_telluric': 'tau_cet_crires_plus_molecfit',
}


def _holding_for_audit_dir(name: str) -> str:
    return _AUDIT_DIR_HOLDING.get(name, name)


def _stellar_crires(audit_root: Path) -> list[dict]:
    """Per-frame D1 residuals from the CRIRES+ stellar corrections (RYA-963, RYA-973)."""
    import csv
    rows = []
    for manifest in sorted(audit_root.glob("rya*_crires_telluric/corrected_manifest.csv")):
        ticket = manifest.parent.name.split("_")[0].upper().replace("RYA", "RYA-")
        for r in csv.DictReader(open(manifest, newline="")):
            try:
                before, after = float(r["gate_before"]), float(r["gate_after"])
            except (KeyError, TypeError, ValueError):
                continue
            rows.append({
                "window": f"{r.get('wlen_id', '?')} ({r.get('band', '?')} band)",
                # Name the INSTRUMENT and the HOLDING. The audit directory is where the
                # evidence lives, not what was corrected, and a reader scanning for
                # "CRIRES" found nothing because every row said rya963_crires_telluric.
                "instrument": "crires_plus",
                "holding": _holding_for_audit_dir(manifest.parent.name),
                "product": manifest.parent.name, "ticket": ticket,
                "star_id": r.get("star_id", ""), "date_obs": r.get("date_obs", ""),
                # NOT-contested and NOBODY-RECORDED-IT are different states and the
                # dashboard must not collapse them: alpha Cen's A/B verdict rests on a
                # branch assignment RYA-963 left contested, and a manifest written before
                # that column existed says nothing about it either way.
                "star_id_contested": (r["star_id_contested"].strip().lower() == "true"
                                      if "star_id_contested" in r else "unrecorded"),
                "metric": "d1_residual",
                "before_d1_residual": before, "after_d1_residual": after,
                "gate_passed": str(r.get("gate_passed", "")).lower() == "true",
                "gdas_profile": r.get("gdas_profile", ""),
            })
    return rows


def _harps_state(audit_root: Path) -> list[dict]:
    """HARPS telluric STATE determinations (RYA-973). A determination is not a
    correction, and the tracker must not let one read as the other: these rows carry
    `state`, never a before/after, because nothing has been corrected."""
    rows = []
    for path in sorted(audit_root.glob("rya*_harps_telluric/*_harps_telluric_state.json")):
        d = json.loads(path.read_text())
        rows.append({
            "window": f"{int(d['o2b_window_A'][0])}-{int(d['o2b_window_A'][1])} A",
            "product": f"{d['star']}_harps", "ticket": d.get("ticket", ""),
            "metric": "state_only__NOT_corrected",
            "header_state": d.get("header_state"),
            "flux_state": d.get("flux_state"),
            "o2b_median_frac_below": d.get("o2b_median_frac_below"),
            "control_median_frac_below": d.get("control_median_frac_below"),
            "excess_ratio": d.get("excess_ratio"),
            "n_products": d.get("n_products"),
            "note": d.get("disposition", ""),
        })
    return rows


def telluric_summary(rows: list[dict], instruments: list[dict]) -> dict:
    """The narrative state of telluric work, DERIVED from the rows the collectors found.

    The page needs prose, not just a table — a reader should be able to see at a glance
    what is corrected, what is merely determined, and what is blocked and why. But prose
    with numbers typed into it is the exact failure RYA-935 exists to prevent, so every
    figure here is computed from the same evidence the tables render, and the page only
    formats it."""
    from collections import defaultdict
    corrected = defaultdict(lambda: {'products': 0, 'pass': 0, 'fail': 0,
                                     'bands': set(), 'settings': set(), 'tickets': set()})
    for r in rows:
        if r.get('metric') != 'd1_residual':
            continue
        k = (r.get('instrument', '?'), r.get('holding', '?'))
        e = corrected[k]
        e['products'] += 1
        e['pass' if r.get('gate_passed') else 'fail'] += 1
        w = str(r.get('window', ''))
        if '(' in w:
            e['bands'].add(w.split('(')[1].split()[0])
        e['settings'].add(w.split(' (')[0])
        if r.get('ticket'):
            e['tickets'].add(r['ticket'])

    solar = [r for r in rows if r.get('metric', 'pct_below_0.5') == 'pct_below_0.5']
    solar_corrected = [r for r in solar if r.get('after_pct_below_0.5') is not None]
    determined = [{'holding': r.get('product'), 'window': r.get('window'),
                   'header_state': r.get('header_state'), 'flux_state': r.get('flux_state'),
                   'excess_ratio': r.get('excess_ratio'), 'ticket': r.get('ticket')}
                  for r in rows if str(r.get('metric', '')).startswith('state_only')]
    uncorrected = [i for i in instruments if i.get('telluric_applied') == 'not-applied']

    return {
        'corrected_stellar': [
            {'instrument': k[0], 'holding': k[1], 'products': v['products'],
             'gates_pass': v['pass'], 'gates_fail': v['fail'],
             'bands': sorted(v['bands']), 'settings': sorted(v['settings']),
             'tickets': sorted(v['tickets'])}
            for k, v in sorted(corrected.items())],
        'corrected_solar_bands': len(solar_corrected),
        'solar_bands_not_corrected': len(solar) - len(solar_corrected),
        'determined_not_corrected': determined,
        'holdings_still_not_applied': [i.get('holding') for i in uncorrected],
        # Facts about the world, each carrying the ticket that established it. Declared
        # rather than derived because no artifact in this repo can measure them.
        'blockers': [
            {'what': 'HARPS / any La Silla instrument',
             'why': ('ESO ships a GDAS tarball for PARANAL ONLY, so a per-night La Silla '
                     'profile had to be pulled by hand. RYA-983 automated the NOAA ARL '
                     'archive (byte-range, 10.7 MB per night rather than 599 MB), so this '
                     'is unblocked but not yet run.'),
             'ticket': 'RYA-983'},
            {'what': 'any night before 2004-12-01',
             'why': ('the NOAA ARL GDAS1 archive begins there. Under the RYA-380 '
                     'no-standard-atmosphere rule those frames can never be corrected by '
                     'this route — a PERMANENT gap, not a backlog item. Three of tau '
                     "Ceti's 32 HARPS nights sit below it."),
             'ticket': 'RYA-983'},
            {'what': 'scoring a corrected stellar spectrum',
             'why': ('measure_band_ew._INSTRUMENT_HOLDINGS still contains no non-solar '
                     'holding and HoldingSpec carries no star, so there is nothing to '
                     'point the harness at. RYA-985 made the synthesis star-generic; the '
                     'holding half is open.'),
             'ticket': 'RYA-985'},
        ],
        'caveats': [
            {'what': 'the K-band telluric anchor',
             'note': ('Y/J/H close within 1.93 km/s on every frame; K2148 and K2192 rail '
                      'on tau Ceti and K2148 failed on alpha Cen. The CORRECTION is '
                      'unaffected (the fit is topocentric and never uses the zero-point) '
                      'but an absolute RV from a K frame is not trustworthy.'),
             'ticket': 'RYA-973'},
            {'what': "alpha Cen's A/B star id",
             'note': ('rests on the acen_orbit branch assignment RYA-963 left CONTESTED: '
                      'the CRIRES pair says it is inverted, RYA-384\'s NIRPS anchor says '
                      'it is not, and both cannot be right.'),
             'ticket': 'RYA-963'},
        ],
    }



def collect_model_matrix() -> dict:
    """RYA-1015 element x model-type availability + engines + molecules.

    Loud-fails visibly (an `error` key the page renders) rather than dropping the
    section, so a blind Sirius scan can never look like "no grids".
    """
    try:
        from pipeline.model_availability_matrix import (
            build_matrix, build_engine_matrix, build_molecule_matrix)
        m = build_matrix()
        m["engines"] = build_engine_matrix(m)
        m["molecules"] = build_molecule_matrix()
        from pipeline.model_availability_matrix import write_findings_csv
        from pathlib import Path as _P
        write_findings_csv(m, m["engines"], m["molecules"],
                           _P(__file__).resolve().parents[1] / "data" / "results"
                           / "rya935" / "model_availability_findings.csv")
        m["findings_csv"] = "model_availability_findings.csv"
        return m
    except Exception as exc:                      # noqa: BLE001 - surfaced on the page
        return {"error": f"{type(exc).__name__}: {exc}", "cells": [],
                "engines": [], "molecules": [], "problem_count": None}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--products-root", type=Path, action="append", default=None)
    ap.add_argument("--audit-root", type=Path, default=ROOT / "data" / "audit")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "results" / "rya935"
                    / "live_status.json")
    ap.add_argument("--refresh-seconds", type=int, default=5)
    ap.add_argument("--products-store", type=Path,
                    default=ROOT / "data" / "products",
                    help="the element feed the PUBLIC site reads. Its gated products[] is "
                         "this page's source (RYA-1097).")
    ap.add_argument("--from-results", action="store_true",
                    help="⚠️ fall back to scanning data/results/ for *_products.csv -- the "
                         "PRE-RYA-1097 behaviour, which shows band products that were "
                         "never published and so were never gated. For inspecting an "
                         "unpublished run, never for the live page.")
    args = ap.parse_args()

    import measure_band_ew as M
    instruments = set(M._INSTRUMENT_HOLDINGS)
    holdings = {s.holding_id for specs in M._INSTRUMENT_HOLDINGS.values() for s in specs}
    roots = args.products_root or [ROOT / "data" / "results"]
    _EXTRA_RESULT_ROOTS.extend(Path(r).resolve() for r in roots)

    # 🔴 RYA-1097 — THE FEED IS THE SOURCE. `--from-results` keeps the old directory scan
    # for a run that has not been published yet, and it is NOT the default: the whole
    # defect was that this page could show what the gate had just withdrawn.
    if args.from_results:
        products = collect_products(roots, instruments, holdings)
        feed_withheld: list[dict] = []
    else:
        products, feed_withheld = collect_feed_products(
            args.products_store, instruments, holdings)
    status = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "scripts/rya935_live_status.py",
        "refresh_seconds": args.refresh_seconds,
        "derivation_note": ("Every value here is read from a product or a registry. "
                            "Nothing is typed. The ticket pipeline is deliberately "
                            "absent: it is Linear state, and a committed copy drifts."),
        # EVERY species in the frozen gold reference, not merely those with products.
        # This page is the progress framework for the whole solar calibration and
        # then for every star after it, so early on the absences ARE the content --
        # listing only what is finished would hide the work that remains.
        "elements": None,   # filled below, once the reference is read
        "elements_with_products": sorted({p["element"] + p["ion"] for p in products}),
        "bands": ["near-UV", "VIS", "red-optical", "NIR"],
        "instruments": (_inst := collect_instruments()),
        "products": products,
        "telluric": (_tel := collect_telluric(args.audit_root)),
        "telluric_summary": None,       # filled below; needs both lists
        "reference": collect_reference(ROOT),
        # 🔴 RYA-1031: SUPERSEDED, and it was rendering as "PRIMARY". This is the older
        # RYA-850 graded lineage (Fe I VIS n=9, A=7.445) read from rya850_summary.json.
        # It sat at the TOP of the forest labelled PRIMARY, above the RYA-933 graded
        # products measured on the same band with 67 lab-gf lines. Two graded lineages on
        # one plot with the SMALLER one called primary is worse than showing neither --
        # the reader cannot tell which pool the headline refers to. Kept in the artifact
        # under `graded_superseded`; off the plot.
        "graded": [],
        "graded_superseded": collect_graded(ROOT),
        "model_matrix": collect_model_matrix(),
        "reporting_contract": {
            "primary": "graded lab-gf pool, on its own CITED pool sigma (RYA-850)",
            "secondary": "ungraded all-lines pool, on the 0.17 dex gf placeholder",
            "headline_rule": "the ungraded value is NEVER the headline (RYA-851)",
            "bars": "statistical SOLID, systematic WIREFRAME -- never summed; "
                    "error_budget.py deliberately provides no combined()",
        },
    }
    # 🔴 THE ROSTER STAYS KEYED BY SPECIES. The page joins products to the roster with
    # `p.element + p.ion === sp`, so rewriting "FeI" to "Fe" here silently broke every
    # join and the page rendered nothing at all -- a display change is not allowed to
    # move a join key (RYA-906). The label travels ALONGSIDE the key instead.
    status["elements"] = sorted(status["reference"])
    status["element_labels"] = {sp: species_display(sp) for sp in status["elements"]}
    status["telluric_summary"] = telluric_summary(_tel, _inst)
    status["system"] = "solar"   # the framework is per-star; this build is the Sun

    # WHY a cell is empty. Computed once per (holding, band) and NOT per element,
    # because every blocker the resolver reports -- wiring, coverage, the telluric
    # gate -- is a property of the holding and the band. Recomputing it per element
    # would be 26x the work for identical answers, and would invite someone to read
    # element-specific meaning into a reason that has none.
    from pipeline import band_policy
    reachability = {}
    for inst in status["instruments"]:
        for policy in band_policy.POLICIES:
            try:
                reason = why_no_product("Fe", "I", inst["holding"], inst["instrument"],
                                        policy.lo_A, policy.hi_A)
            except Exception as exc:                            # noqa: BLE001
                reason = f"could not resolve: {exc}"
            reachability[f"{inst['holding']}|{policy.name}"] = reason
    status["reachability"] = reachability
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(status, indent=2) + "\n")
    have = ", ".join(status["elements_with_products"]) or "(none)"
    print(f"{len(products)} product rows across {have}; "
          f"{len(status['elements'])} species tracked -> {args.out}")
    unattributed = sum(1 for p in products if p["holding"] is None)
    status["unattributed_products"] = unattributed
    status["variant_products"] = sum(1 for p in products if p["is_variant"])
    status["run_contexts"] = sorted({p["run_context"] for p in products})
    for row in products:
        row.update(telluric_state_of(row, row.get("committed")))
    status["pre_correction_products"] = sum(
        1 for p in products if p["telluric_basis"].startswith("PRE"))
    # 🔴 RYA-1026 display gate. Applied HERE, after every row is built, so the split is
    # visible in the artifact rather than done silently in the page's javascript.
    # 🔴 RYA-1026's RATIFIED gate, not a second one. This called a display_class()
    # hand-rolled here, which duplicated `pipeline.telluric_display_policy` and would
    # have drifted from it the first time the policy changed -- two implementations of
    # one rule is how a ratified decision quietly stops being enforced.
    #
    # `display_state` is used rather than `assert_displayable`: the tracker must SHOW a
    # blocked product with its reason, and raising would kill the whole page over one
    # bad row. An omission reads as "no data", which is absence-as-conclusion (RYA-833).
    from pipeline.telluric_display_policy import display_state, anomaly
    _science, _withheld = [], []
    for row in products:
        holding = row.get("holding")
        state = display_state(holding) if holding else "UNREGISTERED"
        row["display_state"] = state
        # CLEAN_WITH_ANOMALY is not CLEAN, and the difference is the anomaly text.
        # Carrying the state without it would render a caveated holding as unqualified.
        row["telluric_anomaly"] = anomaly(holding) if holding else None
        # The gf-graded tier is what the page reports (RYA-1026). DEEPGRADED is rung 3
        # too and is still withheld: it is a DIFFERENT line selection (the 109 saturated
        # lines above the EW depth gate), not a second opinion on the same 67, and
        # mixing selections in one view reads a selection difference as a measurement
        # difference (RYA-842/984).
        sel = str(row.get("selector") or "")
        if not sel.startswith("GRADED"):
            row["not_displayed_because"] = (
                f"selector {sel!r} is not the graded tier; RYA-1026 reports the "
                f"gf-graded product. Kept on disk and listed here, not deleted")
            _withheld.append(row)
        elif state in ("CLEAN", "CLEAN_WITH_ANOMALY"):
            _science.append(row)
        else:
            row["not_displayed_because"] = (
                f"holding {holding!r} is {state}; RYA-1026 displays telluric-corrected "
                f"input only")
            _withheld.append(row)
    # 🔴 RYA-1097 — THE FEED'S OWN WITHDRAWALS ARE LISTED, EXPLICITLY LABELLED, AND NEVER
    # MIXED INTO THE LIVE VIEW. `_withheld` above is this page's DISPLAY policy (RYA-1026
    # selector/telluric); `feed_withheld` is the eligibility gate's verdict, carried with
    # the code and reason it recorded. Two different questions, kept apart: "we chose not
    # to show this" and "this is not publishable".
    status["products"] = _science
    status["products_withheld"] = _withheld + feed_withheld
    status["withheld_summary"] = {
        "n_withheld_by_display_policy": len(_withheld),
        "n_withdrawn_in_the_feed": len(feed_withheld),
        "source": ("the gated products[] of data/products/<star>/<El>.json -- the same "
                   "feed the public site reads (RYA-1097)" if not args.from_results else
                   "⚠️ data/results/ scan (--from-results): these were never gated"),
        "rule": "RYA-1026: displayed science is telluric-corrected input only; the "
                "KP2005-vs-KP1984 control pair is the whitelisted exception",
        "note": "Withheld rows are NOT deleted -- they are on disk and listed here with "
                "the reason, so a gap in the page reads as 'owed a corrected re-run' "
                "rather than as 'nothing measured' (RYA-833).",
    }
    # 🔴 A TRUE DUPLICATE IS A BUG; A SHARED LABEL IS NOT. Two live records with the SAME
    # identity mean the page cannot say which one IS the product for that cell. Two records
    # that merely render the same string are DISTINCT products the display collapses --
    # e.g. "Synth · 1D-LTE" on molecfit and on kurucz2005 -- and the fix there is to say
    # what differs, not to hide one.
    status["identity"] = identity_report(_science)
    args.out.write_text(json.dumps(status, indent=2) + "\n")
    for i in status["instruments"]:
        have = {p["band"] for p in products if p["holding"] == i["holding"]}
        print(f"  {i['holding']:<34} telluric={i['telluric_applied']:<12} "
              f"bands: {sorted(have) or '—'}")
    if unattributed:
        print(f"\n  {unattributed} product rows predate RYA-933/934 and name no "
              f"holding in their filename.\n  They are attributable to an INSTRUMENT "
              f"only. Re-run to attribute them.")


if __name__ == "__main__":
    main()

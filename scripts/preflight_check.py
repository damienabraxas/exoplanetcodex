#!/usr/bin/env python3
"""Readiness reconciliation, run BEFORE an element is measured — RYA-905.

    python3 scripts/preflight_check.py --star sun --element Fe
    python3 scripts/preflight_check.py --star sun --element Al --json out.json

WHAT NIGHT THIS COMES FROM
--------------------------
The Fe/HARPS night (RYA-896/897/898) cost a working day to a class of defect that never
touched an abundance value: **things silently MISSING that looked present**. The direct
solar HARPS arm was configured in `scripts/measure_band_ew.py` — it is a key of
`PRE_NORMALISED` to this day — but the `load_window` dispatch never got the branch. The
holding was verified, the config named it, and no code path could read it. Nothing
raised. It surfaced only when a human eyeballed the rendered page and counted
instruments by hand.

The failure was SILENCE, not error. So the defence is a readiness check that runs before
each element and turns silent absence into a visible warning.

ADVISORY, NEVER A GATE
----------------------
This exits 0 whatever it finds. A survey legitimately lacks data — no UVES for a given
star is normal, not a showstopper — and a gate that fires on the normal state of a survey
is a gate that gets disabled. It WARNS, and it suggests tickets. Ryan reads the report,
confirms the expected absences, and files the rest.

THE TWO SEVERITIES ARE THE ENTIRE POINT
---------------------------------------
    INFO   expected-absence.  A holding / line / grid we simply do not have for this
                              element or star. The normal survey state. No alarm.
    WARN   silent-gap.        Something we DO have — verified, reachable-should-be —
                              that the pipeline cannot actually see.

"We're missing HARPS data" (we weren't) and "we HAVE verified HARPS and the code cannot
reach it" (the bug) were INDISTINGUISHABLE from the rendered page. A check that cannot
tell those two apart is useless, so every finding below carries the discriminator that
separated them.

A third severity, ERROR, is reserved for the check itself being broken — a registry that
will not load, a dispatch this module can no longer read. An ERROR means the report's
absences are UNTRUSTWORTHY and must not be read as evidence of anything. That distinction
exists because the cheapest way to make a reconciliation report look clean is to break its
own ability to see (RYA-833: *an absence is a hypothesis, never a conclusion*).

CONTROLS
--------
Two of the six checks assert a NEGATIVE ("this is not reachable", "this line is not
covered"), and a negative needs a positive control or it is unfalsifiable. So:

  * the loader-dispatch reader must find `kpno_solar_atlas` — the one instrument every
    committed band product in this repo was measured on. If the reader cannot see the
    arm we KNOW is wired, its report that some other arm is unwired is worthless, and the
    check returns ERROR instead of a page of false WARNs.
  * the same reader must NOT report a sentinel instrument name that cannot exist. A
    reader that says yes to everything is as broken as one that says no to everything.

READ-ONLY, AND SINGLE-SOURCED
-----------------------------
This module WRITES NOTHING except its own report (and only when `--json` asks). It does
not re-verify anything: the intake framework already produced the ground truth, and every
registry below is read from the ONE place that owns it — never a copy, never a value
restated here:

    holdings + telluric_applied   data/catalog/holdings_manifest_registry.csv   RYA-806
    instrument axis               data/catalog/instrument_catalog.csv           RYA-786
    telluric verdict              pipeline.telluric_policy.gate_holding         RYA-806
    band harness dispatch         scripts/measure_band_ew.py (read, not run)    RYA-897
    band regimes                  pipeline.band_policy                          RYA-713
    NLTE grids                    data/nlte_grids/ + nlte_grid_availability.csv
    line-list coverage            data/audit/vald_inventory/coverage_matrix.csv RYA-376
    the anchor's declared chain   data/audit/cno_synthesis/solar_phase_c_verdict.json
    the frozen anchor value       data/reference/solar/elements/<El>_<ion>.json RYA-814

WHY THE DISPATCH IS READ AND NOT IMPORTED
-----------------------------------------
`scripts/measure_band_ew.py` resolves the Kitt Peak atlas AT IMPORT and raises SystemExit
when it is absent, so importing it to ask a question about reachability would make the
preflight impossible on any machine without the spectra staged — i.e. it would fail
hardest exactly where a readiness check is most needed. The dispatch is therefore read
statically, from the AST of the harness's own source. That is a read of the dispatch, not
a duplicate of it: nothing here restates which instruments are wired.
"""
from __future__ import annotations

import argparse
import ast
import json
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path

# Standalone-script bootstrap (RYA-313): repo root on sys.path BEFORE importing
# config/pipeline, so this runs from any cwd. Derived from __file__, never cwd.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pandas as pd  # noqa: E402

from pipeline import band_policy  # noqa: E402
from pipeline import telluric_policy  # noqa: E402
from pipeline.species import parse_ion  # noqa: E402  RYA-345 canonical ion normalizer

# ── The registries. One line each, and each is the owner of its fact ──────────
HOLDINGS_CSV = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"
INSTRUMENTS_CSV = ROOT / "data" / "catalog" / "instrument_catalog.csv"
SYSTEMS_CSV = ROOT / "data" / "catalog" / "system_catalog.csv"
HARNESS_PY = ROOT / "scripts" / "measure_band_ew.py"
GRID_DIR = ROOT / "data" / "nlte_grids"
GRID_AVAILABILITY_CSV = ROOT / "data" / "curation" / "nlte_grid_availability.csv"
VALD_COVERAGE_CSV = ROOT / "data" / "audit" / "vald_inventory" / "coverage_matrix.csv"
LITSCAN_DIR = ROOT / "data" / "reference" / "litscan"
ANCHOR_VERDICT_JSON = ROOT / "data" / "audit" / "cno_synthesis" / "solar_phase_c_verdict.json"
ANCHOR_ELEMENT_DIR = ROOT / "data" / "reference" / "solar" / "elements"
CANONICAL_GF_CSV = ROOT / "data" / "linelists" / "canonical_gf.csv"
#: RYA-709's per-line accounting ledger — every usable line, which arms reach it, and
#: whether it made the measured pool. A line ON this ledger is not invisible: it is
#: counted, and `summary.csv` carries the per-element totals. That is what separates
#: "not measured" from "nobody knows this line exists".
LINE_ACCOUNTING_CSV = ROOT / "data" / "audit" / "line_accounting" / "per_line.csv"

#: Star aliases. The registries key on `system_id` / `star_params_key`; a human types
#: "sun". Resolved here rather than in six places.
STAR_ALIASES = {"sun": "solar", "sol": "solar", "the sun": "solar"}

#: Evidence states that mean "we hold this, and someone checked". A holding at a weaker
#: state is not a silent gap — it is an intake step that has not finished, which is a
#: different (INFO) thing.
VERIFIED_STATES = ("verified",)

#: The gf grades this project treats as a graded laboratory measurement. D/E/F are the
#: cull tiers (RYA-398) and an ungraded blank is Kurucz-theoretical (RYA-161).
GRADED_NIST = ("AAA", "AA", "A+", "A", "B+", "B", "C+", "C")

#: Wavelength match tolerance, in Angstrom, when joining a line across two catalogues.
#: NOT a fudge factor: catalogues disagree about a line's air wavelength at the ~0.1 A
#: level (Al I 6696.185 in canonical_gf vs 6696.023 in the Amarsi-2020 grid is the SAME
#: line), while the measured separability floor is 0.30 A (RYA-761). A tolerance below
#: the catalogue disagreement MANUFACTURES an absence — it reports a covered line as
#: uncovered — and one above the separability floor merges two real lines. 0.20 A sits
#: between the two, and every match reports the delta it used so the join is auditable.
WAVE_TOL_A = 0.20

#: A name no instrument can have. The negative control for the dispatch reader.
DISPATCH_SENTINEL = "__preflight_no_such_instrument__"

#: The RYA-1030 registry column and its vocabulary. IMPORTED, never re-spelled: a second
#: copy of a state name drifts silently and this check would then pass on a typo.
from pipeline.normalization_intake import NORMALISED as NI_NORMALISED  # noqa: E402
from pipeline.normalization_intake import UNKNOWN as NI_UNKNOWN        # noqa: E402

COLUMN_NORMALISATION = "normalization_state"

OK, INFO, WARN, ERROR = "OK", "INFO", "WARN", "ERROR"
_RANK = {OK: 0, INFO: 1, WARN: 2, ERROR: 3}


# ── Findings ──────────────────────────────────────────────────────────────────

@dataclass
class Finding:
    """One statement about readiness.

    `discriminator` is mandatory prose, not decoration: it records WHAT separated this
    from the other severity. A WARN whose author could not say why it is not an
    expected-absence has not made the distinction this whole module exists to make.
    """
    check: str
    severity: str
    subject: str
    message: str
    discriminator: str
    suggested_ticket: str | None = None

    def line(self) -> str:
        head = f"  [{self.severity:<5}] {self.subject}: {self.message}"
        if self.severity in (WARN, ERROR):
            head += f"\n           why not expected-absence: {self.discriminator}"
        if self.suggested_ticket:
            head += f"\n           SUGGEST TICKET: {self.suggested_ticket}"
        return head


@dataclass
class CheckResult:
    number: int
    name: str
    findings: list[Finding] = field(default_factory=list)

    @property
    def severity(self) -> str:
        return max((f.severity for f in self.findings), key=lambda s: _RANK[s], default=OK)


# ── Reading the band-harness dispatch (check 1's instrument) ──────────────────

@dataclass
class Dispatch:
    """What `measure_band_ew.load_window` can actually serve, read from its source."""
    instruments: tuple[str, ...]
    #: instrument -> the holdings that arm can serve, IN PREFERENCE ORDER.
    #: 🔴 WAS `dict[str, str]` — one holding per instrument — until RYA-904 proved that
    #: shape wrong in the harness itself. `crires_plus` has three solar holdings and the
    #: single-valued map named the ONE the telluric gate refuses, so both corrected
    #: holdings read as unreachable. A reader that mirrors a broken shape reproduces the
    #: broken answer, so this is a tuple.
    served_holdings: dict[str, tuple[str, ...]]
    configured: tuple[str, ...]         # instruments named by harness config tables
    controls_ok: bool
    control_note: str


def read_dispatch(harness_py: Path | None = None) -> Dispatch:
    """Which instruments `load_window` has a branch for — statically, with controls.

    Reads three things out of the harness's own AST:

      * the instrument literals `load_window` compares against — the dispatch itself;
      * `_INSTRUMENT_HOLDINGS`, the harness's own declaration (RYA-806/904) of WHICH
        HOLDINGS each arm serves, because "the instrument is wired" and "this holding is
        reachable" are not the same claim. ⚠️ READ BOTH SHAPES: RYA-904 replaced the
        single-valued `_LOADER_HOLDING` with a preference-ordered table of `HoldingSpec`
        entries, and the legacy name is still parsed so this reader keeps working against
        a harness that predates it (and against the synthetic fixtures the controls use);
      * the instrument keys of the harness's config tables (`PRE_NORMALISED` and friends).
        The RYA-897 signature is precisely an instrument that is CONFIGURED but NOT
        DISPATCHED, and naming that signature is more useful than "not wired".
    """
    # Resolved from the module attribute at CALL time, not bound as a default at def
    # time: a default would freeze the path this reader looks at, which is exactly what
    # makes a reader untestable — the control that proves it reacts to a harness WITH the
    # loader could never point it at one.
    harness_py = Path(harness_py) if harness_py is not None else HARNESS_PY
    tree = ast.parse(harness_py.read_text(), filename=str(harness_py))

    dispatch: list[str] = []
    served: dict[str, tuple[str, ...]] = {}
    configured: set[str] = set()

    for node in ast.walk(tree):
        # the dispatch: `if instrument == "kpno_solar_atlas": ...`
        if isinstance(node, ast.FunctionDef) and node.name == "load_window":
            for sub in ast.walk(node):
                if not isinstance(sub, ast.Compare):
                    continue
                names = {n.id for n in ast.walk(sub.left) if isinstance(n, ast.Name)}
                if "instrument" not in names:
                    continue
                for comp in sub.comparators:
                    for const in ast.walk(comp):
                        if isinstance(const, ast.Constant) and isinstance(const.value, str):
                            dispatch.append(const.value)
        # the module-level tables.
        # ⚠️ `ast.AnnAssign` TOO — RYA-904. `_INSTRUMENT_HOLDINGS: dict[...] = {...}` is
        # an ANNOTATED assignment, which is a different node type, and an `ast.Assign`-only
        # walk sees nothing at all. That is not a soft failure: the reader returned an
        # EMPTY dispatch, which would have read as "no instrument is wired". The only
        # reason it was caught is that RYA-905's positive control refuses to report
        # absences when it cannot see kpno_solar_atlas — the absence would otherwise have
        # been a page of confident false WARNs (RYA-833).
        _is_assign = isinstance(node, (ast.Assign, ast.AnnAssign))
        if _is_assign and isinstance(getattr(node, "value", None), ast.Dict):
            _tgts = node.targets if isinstance(node, ast.Assign) else [node.target]
            targets = [t.id for t in _tgts if isinstance(t, ast.Name)]
            keys = [k.value for k in node.value.keys
                    if isinstance(k, ast.Constant) and isinstance(k.value, str)]
            if "_INSTRUMENT_HOLDINGS" in targets:
                # RYA-904 shape: instrument -> (HoldingSpec("id", ...), ...). The holding
                # id is each spec's FIRST argument, positional or `holding_id=`. Read
                # positionally rather than by evaluating the call, because this reader is
                # deliberately static -- it must work against a harness it cannot import
                # (measure_band_ew resolves the Kitt Peak atlas at module import).
                for k, v in zip(keys, node.value.values):
                    ids: list[str] = []
                    for elt in getattr(v, "elts", []):
                        if not isinstance(elt, ast.Call):
                            continue
                        first = (elt.args[0] if elt.args else next(
                            (kw.value for kw in elt.keywords if kw.arg == "holding_id"),
                            None))
                        if isinstance(first, ast.Constant) and isinstance(first.value, str):
                            ids.append(first.value)
                    if ids:
                        served[k] = tuple(ids)
                        dispatch.append(k)
            elif "_LOADER_HOLDING" in targets:
                # LEGACY single-valued shape (pre-RYA-904). Widened to a 1-tuple here so
                # everything downstream speaks one language.
                vals = [v.value if isinstance(v, ast.Constant) else None
                        for v in node.value.values]
                served.update({k: (v,) for k, v in zip(keys, vals) if isinstance(v, str)})
            elif targets:
                configured.update(keys)

    instruments = tuple(sorted(set(dispatch)))

    # ── the controls ─────────────────────────────────────────────────────────
    # POSITIVE: every committed band product in this repo was measured on the Kitt Peak
    # atlas, so `kpno_solar_atlas` IS wired. A reader that cannot see it is broken, and
    # its report that some other arm is unwired would be a page of false WARNs.
    # NEGATIVE: a reader that reports an impossible name is equally broken.
    pos = "kpno_solar_atlas" in instruments
    neg = DISPATCH_SENTINEL not in instruments
    if pos and neg:
        note = (f"positive control: kpno_solar_atlas found among {len(instruments)} "
                f"dispatched arms; negative control: sentinel absent")
    else:
        note = (f"CONTROL FAILED (positive={pos}, negative={neg}) — the dispatch reader "
                f"parsed {instruments!r} out of {harness_py.name}. Its absences prove "
                f"nothing until this is fixed.")
    return Dispatch(instruments, served, tuple(sorted(configured)), pos and neg, note)


# ── Loading the state this reconciles ─────────────────────────────────────────

@dataclass
class State:
    star: str
    element: str
    ion: str
    species: str
    holdings: pd.DataFrame          # rows for THIS star
    all_holdings: pd.DataFrame
    instruments: pd.DataFrame
    dispatch: Dispatch
    roots: tuple[Path, ...]
    ew_runs: list[dict]             # measured EW artifacts found, parsed
    products: pd.DataFrame          # rendered product rows for this species
    product_sources: tuple[str, ...]
    best_lines: pd.DataFrame        # wavelength_air_A, basis
    best_line_basis: str


def _resolve_star(raw: str, holdings: pd.DataFrame) -> str:
    key = STAR_ALIASES.get(raw.strip().lower(), raw.strip().lower())
    known = set(holdings.system_id.astype(str))
    if SYSTEMS_CSV.exists():
        sys_df = pd.read_csv(SYSTEMS_CSV, comment="#")
        known |= {str(x) for x in sys_df.get("star_params_key", pd.Series(dtype=str)).dropna()}
    if key not in known:
        raise SystemExit(
            f"star {raw!r} resolves to {key!r}, which is in neither {HOLDINGS_CSV.name} "
            f"nor {SYSTEMS_CSV.name}. Known: {', '.join(sorted(known))}. Refusing to "
            f"report 'no holdings' for a star id that does not exist — that would be a "
            f"MANUFACTURED absence (RYA-833), not a finding.")
    return key


def _parse_ew_name(path: Path) -> dict | None:
    """`FeI_3800_6910_kpno_solar_atlas_PROFILEFIT_ew.csv` -> its run coordinates.

    The harness names its own outputs; parsing that name is how a reconciliation learns
    WHICH BAND WAS ACTUALLY RUN, which is the difference between "this line was dropped"
    and "nothing has been measured here yet".
    """
    stem = path.name
    for suffix in ("_ew.csv",):
        if not stem.endswith(suffix):
            return None
        stem = stem[: -len(suffix)]
    parts = stem.split("_")
    if len(parts) < 4:
        return None
    species_tok, lo_tok, hi_tok = parts[0], parts[1], parts[2]
    try:
        lo, hi = float(lo_tok), float(hi_tok)
    except ValueError:
        return None
    handler = parts[-1]
    instrument = "_".join(parts[3:-1])
    return {"path": path, "species_token": species_tok, "lo_A": lo, "hi_A": hi,
            "instrument": instrument, "handler": handler}


def _species_token(element: str, ion: str) -> str:
    return f"{element}{ion}"


def _collect_ew_runs(roots, element, ion) -> list[dict]:
    tok = _species_token(element, ion)
    out = []
    for root in roots:
        d = Path(root) / "data" / "measured" / "band_ew"
        if not d.is_dir():
            continue
        for p in sorted(d.glob(f"{tok}_*_ew.csv")):
            meta = _parse_ew_name(p)
            if meta is None:
                continue
            try:
                meta["frame"] = pd.read_csv(p)
            except Exception as exc:          # a corrupt artifact is not an absence
                meta["frame"] = None
                meta["error"] = str(exc)
            out.append(meta)
    return out


def _collect_products(roots, element, ion) -> tuple[pd.DataFrame, tuple[str, ...]]:
    """Every rendered product row for this species, from every root inspected.

    Two shapes are read, because the project has two: the per-band `*_products.csv` the
    harness writes, and the committed cross-band `*_product_matrix.csv` that the site and
    the tickets quote. A reconciliation that read only one of them would call a band
    missing because it looked in the wrong file.
    """
    frames, sources = [], []
    for root in roots:
        root = Path(root)
        for pat in ("data/results/band_products/**/*_products.csv",
                    "data/results/**/*product_matrix*.csv"):
            for p in sorted(root.glob(pat)):
                try:
                    df = pd.read_csv(p)
                except Exception:
                    continue
                if "element" in df.columns:
                    df = df[df.element.astype(str) == element]
                    if "ion" in df.columns:
                        df = df[df.ion.astype(str) == ion]
                else:
                    # the per-band file carries the species in its NAME, not a column
                    if not p.name.startswith(_species_token(element, ion)):
                        continue
                    df = df.copy()
                    df["element"], df["ion"] = element, ion
                if not len(df):
                    continue
                df = df.copy()
                df["_artifact"] = str(p)
                if "instrument" not in df.columns:
                    df["instrument"] = _instrument_from_name(p.name)
                frames.append(df)
                sources.append(str(p))
    if not frames:
        return pd.DataFrame(columns=["element", "ion", "instrument", "_artifact"]), ()
    return pd.concat(frames, ignore_index=True), tuple(sources)


def _instrument_from_name(name: str) -> str | None:
    """Recover the instrument from a harness-named product file, or None.

    None, deliberately, rather than a guess: a product whose instrument cannot be
    established is not a product measured on some default arm (that is exactly how a
    number gets labelled with an instrument it was not measured on — the reason
    `load_window` is loud on an unknown instrument).
    """
    stem = name.rsplit(".", 1)[0]
    for tail in ("_products", "_per_line"):
        if stem.endswith(tail):
            stem = stem[: -len(tail)]
    parts = stem.split("_")
    for i, tok in enumerate(parts):
        if tok.isdigit() and i + 1 < len(parts):
            rest = parts[i + 1:]
            # drop a trailing handler token (PROFILEFIT / SYNTH / LABGF ...)
            if rest and rest[-1].isupper():
                rest = rest[:-1]
            return "_".join(rest) or None
    rest = parts[1:]
    if rest and rest[-1].isupper():
        rest = rest[:-1]
    return "_".join(rest) or None


def _load_best_lines(element, ion, star, ew_runs, products) -> tuple[pd.DataFrame, str]:
    """This element's BEST lines — the subject of checks 2, 3 and 5.

    THE UNION, NOT A CHOICE (RYA-823). Two independent things make a line "best" here and
    neither subsumes the other:

      a) a line the project has already ACCEPTED into an aggregate (`in_aggregate`). It
         survived feature verification, the blend test and the fit-quality test on real
         data — proof by measurement.
      b) a line carrying a GRADED laboratory gf (NIST A..C). Proof by atomic data. It may
         never have been measured; that is what makes it a candidate rather than a
         result, and exactly what check 2 is for.

    A `best_lines:` block in the litscan projection (data/reference/litscan/<El>.yaml)
    OVERRIDES both when present, because the dossier is the source of truth for what the
    literature considers this element's best lines. No such block exists today for any
    element; the hook is declared so the dossier can drive this without a code change.
    """
    lit = LITSCAN_DIR / f"{element}.yaml"
    if lit.exists():
        declared = _litscan_best_lines(lit)
        if declared:
            return (pd.DataFrame({"wavelength_air_A": declared, "origin": "litscan"}),
                    f"litscan dossier projection ({lit.name} best_lines)")

    rows, origins = [], []
    for run in ew_runs:
        df = run.get("frame")
        if df is None or "wavelength_air_A" not in df.columns:
            continue
        if "in_aggregate" in df.columns:
            keep = df[df.in_aggregate.astype(str).str.lower() == "true"]
        else:
            keep = df
        for w in keep.wavelength_air_A.astype(float):
            rows.append(w)
            origins.append("in-aggregate (measured)")
    for col in ("wave", "wavelength_air_A"):
        if col in products.columns:
            for w in pd.to_numeric(products[col], errors="coerce").dropna():
                rows.append(float(w))
                origins.append("in-aggregate (product per-line)")
            break

    n_measured = len(rows)
    graded = _graded_lines(element, ion)
    rows.extend(graded)
    origins.extend(["graded gf (NIST A-C)"] * len(graded))

    if not rows:
        star_list = ROOT / "data" / "linelists" / f"linelist_{star}.csv"
        if star_list.exists():
            df = pd.read_csv(star_list)
            sel = df[(df.element.astype(str) == element) & (df.ion.astype(str) == ion)]
            waves = sorted(set(round(float(w), 4) for w in sel.wavelength_air_A))
            return (pd.DataFrame({"wavelength_air_A": waves,
                                  "origin": ["star linelist (UNGRADED)"] * len(waves)}),
                    f"star linelist {star_list.name} — every catalogued {element} {ion} "
                    f"line, UNGRADED: neither measured nor lab-graded. The weakest basis; "
                    f"findings below inherit that weakness")
        return pd.DataFrame(columns=["wavelength_air_A", "origin"]), "NO SOURCE FOUND"

    df = pd.DataFrame({"wavelength_air_A": rows, "origin": origins})
    df = df.sort_values("wavelength_air_A").drop_duplicates("wavelength_air_A")
    basis = (f"UNION of {n_measured} in-aggregate measured line(s) and {len(graded)} "
             f"graded-gf line(s) -> {len(df)} distinct")
    return df.reset_index(drop=True), basis


def _litscan_best_lines(path: Path) -> list[float]:
    """`best_lines:` from the litscan projection, without requiring a YAML dependency.

    Deliberately narrow: a flat list of numbers under a top-level `best_lines:` key. If
    the dossier ever grows a richer shape this returns nothing and the union basis is
    used instead — which is a visible, stated fallback, not a silent one.
    """
    out, in_block = [], False
    for raw in path.read_text().splitlines():
        if raw.startswith("best_lines:"):
            in_block = True
            continue
        if in_block:
            s = raw.strip()
            if not s.startswith("-"):
                if s and not raw.startswith(" "):
                    break
                continue
            try:
                out.append(float(s.lstrip("- ").split("#")[0].strip()))
            except ValueError:
                continue
    return out


def _graded_lines(element: str, ion: str) -> list[float]:
    if not CANONICAL_GF_CSV.exists():
        return []
    df = pd.read_csv(CANONICAL_GF_CSV,
                     usecols=["species", "wavelength_air_A", "nist_grade"])
    sel = df[(df.species.astype(str) == f"{element} {ion}")
             & (df.nist_grade.astype(str).str.strip().isin(GRADED_NIST))]
    return sorted(round(float(w), 4) for w in sel.wavelength_air_A)


def _nearest(target: float, pool) -> tuple[float | None, float]:
    best, best_d = None, float("inf")
    for w in pool:
        d = abs(float(w) - target)
        if d < best_d:
            best, best_d = float(w), d
    return best, best_d


def load_state(star_raw: str, element: str, ion: str, roots) -> State:
    holdings = pd.read_csv(HOLDINGS_CSV)
    instruments = pd.read_csv(INSTRUMENTS_CSV, comment="#")
    star = _resolve_star(star_raw, holdings)
    mine = holdings[holdings.system_id.astype(str) == star].copy()
    dispatch = read_dispatch()
    ew_runs = _collect_ew_runs(roots, element, ion)
    products, product_sources = _collect_products(roots, element, ion)
    best, basis = _load_best_lines(element, ion, star, ew_runs, products)
    return State(star=star, element=element, ion=ion, species=f"{element} {ion}",
                 holdings=mine, all_holdings=holdings, instruments=instruments,
                 dispatch=dispatch, roots=tuple(Path(r) for r in roots),
                 ew_runs=ew_runs, products=products, product_sources=product_sources,
                 best_lines=best, best_line_basis=basis)


# ── CHECK 1 — instrument reachability (the RYA-897 class) ─────────────────────

def check_instrument_reachability(st: State) -> CheckResult:
    """Can the band harness actually LOAD every holding we have verified for this star?

    This is the RYA-897 guard, run per element. The discriminator between the two
    severities is possession: an instrument we do not hold is an INFO, and an instrument
    we hold, verified, that no `load_window` branch can read is the WARN — that is the
    exact shape of the night this check exists to prevent.
    """
    res = CheckResult(1, "instrument reachability (band-harness loader)")
    if not st.dispatch.controls_ok:
        res.findings.append(Finding(
            "1", ERROR, "dispatch reader", st.dispatch.control_note,
            "the reader that would decide reachability cannot be trusted, so NOTHING "
            "below it is evidence — neither its OKs nor its absences",
            "FIX: preflight_check.read_dispatch can no longer read measure_band_ew's "
            "load_window dispatch; repair the reader before trusting any preflight"))
        return res

    served = st.dispatch.served_holdings
    # RYA-904 — a gap someone has already WRITTEN DOWN is still a gap, and still WARNs;
    # what changes is that the reader can say who owns it instead of leaving the reader
    # of the report to go and find out. `DECLARED_GAPS` is the existing declaration, not
    # a second one made here.
    try:
        from pipeline.loader_coverage import DECLARED_GAPS as _declared
    except Exception:
        _declared = {}
    for _, h in st.holdings.iterrows():
        hid, inst, state = str(h.holding_id), str(h.instrument_id), str(h.evidence_state)
        _decl = (f" DECLARED in pipeline.loader_coverage.DECLARED_GAPS: {_declared[hid]}"
                 if hid in _declared else "")
        if inst not in st.dispatch.instruments:
            configured = inst in st.dispatch.configured
            sev = WARN if state in VERIFIED_STATES else INFO
            sig = (" It IS named by the harness's own config tables while the dispatch "
                   "has no branch for it — the exact RYA-897 signature." if configured
                   else "")
            sig += _decl
            res.findings.append(Finding(
                "1", sev, hid,
                f"held on {inst} (evidence_state={state}) but "
                f"`measure_band_ew.load_window` has no branch for {inst}: the band "
                f"harness cannot read this holding.{sig}",
                f"we HOLD this and intake marked it {state} — this is not missing data, "
                f"it is data no code path can reach",
                (f"BUILD: band-harness loader for {inst} (holding {hid}) — verified at "
                 f"intake, unreachable by load_window (RYA-897 class)")
                if sev == WARN else None))
            continue

        targets = served.get(inst)
        if targets and hid not in targets:
            res.findings.append(Finding(
                "1", WARN, hid,
                f"{inst} IS dispatched, but the harness declares that arm serves "
                f"{', '.join(targets)} — this holding is not among them, so nothing "
                f"reads it.{_decl}",
                f"the instrument is wired, so a per-instrument check would have called "
                f"this reachable; it is the per-HOLDING join (RYA-806/904) that shows it "
                f"is not",
                f"AUDIT: holding {hid} ({inst}) is verified but the {inst} band-harness "
                f"branch serves {', '.join(targets)} — confirm expected, or wire it"))
        else:
            res.findings.append(Finding(
                "1", OK, hid,
                f"{inst} dispatched by load_window"
                + (" and declared to serve this holding" if targets and hid in targets
                   else ""),
                "held, verified, reachable"))

    held = set(st.holdings.instrument_id.astype(str))
    catalogued = set(st.instruments.instrument_id.astype(str))
    absent = sorted(catalogued - held)
    if absent:
        res.findings.append(Finding(
            "1", INFO, "not held",
            f"{len(absent)} catalogued instrument(s) have no {st.star} holding: "
            f"{', '.join(absent)}",
            "a survey legitimately lacks arms; nothing claims we have these"))
    return res


# ── CHECK 2 — line coverage ──────────────────────────────────────────────────

def check_line_coverage(st: State) -> CheckResult:
    """Are this element's best lines present, and slated to be measured?

    The discriminator here is whether a RUN COVERED THE WAVELENGTH. A best line nobody
    has run a band over is an expected absence. A best line inside a band that WAS run,
    which does not appear in that run's output at all, was dropped without a word — and a
    silent drop is the thing RYA-429's rejection ledger exists to forbid.
    """
    res = CheckResult(2, "line coverage (best lines present + slated)")
    if not len(st.best_lines):
        res.findings.append(Finding(
            "2", INFO, st.species,
            f"no best-line set could be built ({st.best_line_basis})",
            "nothing to reconcile — scoped to the artifacts inspected, not a claim "
            "about the element"))
        return res

    res.findings.append(Finding(
        "2", OK, "best-line basis", st.best_line_basis, "declares what 'best' means here"))

    if not st.ew_runs:
        res.findings.append(Finding(
            "2", INFO, st.species,
            f"no measured-EW artifact for {st.species} under the inspected root(s) "
            f"{', '.join(str(r) for r in st.roots)} — no band has been run here, so no "
            f"line can have been dropped by one",
            "the absence is scoped to these roots, not to the project: pass "
            "--artifact-root to widen it"))
        return res

    measured: dict[float, list[dict]] = {}
    for run in st.ew_runs:
        df = run.get("frame")
        if df is None or "wavelength_air_A" not in df.columns:
            res.findings.append(Finding(
                "2", ERROR, Path(run["path"]).name,
                f"artifact unreadable or has no wavelength column ({run.get('error','')})",
                "an unreadable artifact is not an empty one",
                "FIX: unreadable band-EW artifact breaks preflight line accounting"))
            continue
        for _, r in df.iterrows():
            w = float(r.wavelength_air_A)
            measured.setdefault(round(w, 4), []).append(
                {"run": run, "in_aggregate": str(r.get("in_aggregate", "")).lower() == "true",
                 "reason": str(r.get("excluded_reason", "") or "").strip()})

    ledger = _line_accounting(st)
    catalogued = _star_linelist_waves(st)
    silent_drops, unreasoned, accounted, expected = [], [], 0, 0
    on_ledger = 0
    for w in st.best_lines.wavelength_air_A.astype(float):
        hit, delta = _nearest(w, measured.keys())
        if hit is not None and delta <= WAVE_TOL_A:
            for m in measured[round(hit, 4)]:
                if m["in_aggregate"]:
                    accounted += 1
                elif m["reason"]:
                    accounted += 1
                else:
                    unreasoned.append((w, Path(m["run"]["path"]).name))
            continue
        # NOT IN THE RUN'S OUTPUT — but is it on the accounting ledger? RYA-709 records
        # every usable line and whether it made the measured pool, and `summary.csv`
        # totals them per element. A line ON that ledger has been counted and can be
        # found; that is an unmeasured line, which is the survey's normal state. A line
        # on NEITHER the run nor the ledger is the one nothing knows about.
        near_ledger, dl = _nearest(w, ledger)
        if near_ledger is not None and dl <= WAVE_TOL_A:
            on_ledger += 1
            continue
        covering = [r for r in st.ew_runs if r["lo_A"] <= w <= r["hi_A"]]
        if covering:
            near_ll, dll = _nearest(w, catalogued)
            where = ("catalogued in the star line list" if near_ll is not None
                     and dll <= WAVE_TOL_A else "not in the star line list either")
            silent_drops.append((w, ", ".join(Path(r["path"]).name for r in covering), where))
        else:
            expected += 1

    if accounted:
        res.findings.append(Finding(
            "2", OK, st.species,
            f"{accounted} best line(s) measured and accounted — in the aggregate, or "
            f"excluded WITH a stated reason",
            "an exclusion that states its reason is not a silent gap"))
    if on_ledger:
        res.findings.append(Finding(
            "2", INFO, st.species,
            f"{on_ledger} best line(s) are not in any run's output but ARE on the RYA-709 "
            f"line-accounting ledger ({LINE_ACCOUNTING_CSV.relative_to(ROOT)}) as usable "
            f"and unmeasured — counted, not lost",
            "the project's own ledger holds these lines with their reaching instruments; "
            "an unmeasured line that is counted is the survey's normal state"))
    if expected:
        spans = ", ".join("{:.0f}-{:.0f}".format(r["lo_A"], r["hi_A"]) for r in st.ew_runs)
        res.findings.append(Finding(
            "2", INFO, st.species,
            f"{expected} best line(s) lie outside every band that has been run "
            f"({spans} A)",
            "no run covered the wavelength, so no run could have dropped it"))
    if silent_drops:
        in_list = sum(1 for _, _, where in silent_drops if where.startswith("catalogued"))
        shown = "; ".join(f"{w:.3f} A ({where}; run {n})" for w, n, where in silent_drops[:6])
        res.findings.append(Finding(
            "2", WARN, st.species,
            f"{len(silent_drops)} best line(s) sit INSIDE a band that was run, appear "
            f"nowhere in that run's output, AND are absent from the line-accounting "
            f"ledger — {in_list} of them ARE catalogued in the star line list. "
            f"{shown}" + (" ..." if len(silent_drops) > 6 else ""),
            "the band was measured and the ledger that counts unmeasured lines does not "
            "hold these either — a filter upstream of both removed them, and the lines "
            "it removed are recorded nowhere (RYA-429 class)",
            f"AUDIT: {len(silent_drops)} {st.species} best line(s) fall inside a measured "
            f"band but appear in neither its EW artifact nor the RYA-709 line-accounting "
            f"ledger — name the filter that drops them and ledger its rejections"))
    if unreasoned:
        shown = "; ".join(f"{w:.3f} A ({n})" for w, n in unreasoned[:8])
        res.findings.append(Finding(
            "2", WARN, st.species,
            f"{len(unreasoned)} best line(s) excluded from the aggregate with an EMPTY "
            f"reason: {shown}" + (" ..." if len(unreasoned) > 8 else ""),
            "the line WAS measured and then excluded — an exclusion with no stated cause "
            "is indistinguishable from tuning (RYA-844)",
            f"FIX: {st.species} lines excluded with a blank excluded_reason — every "
            f"exclusion must state its cause"))
    return res


# ── CHECK 3 — NLTE grid reach ────────────────────────────────────────────────

def check_grid_reach(st: State) -> CheckResult:
    """Does the NLTE grid cover this element's BEST lines — not merely SOME lines?

    The RYA-773 Al gap is the shape: the Amarsi-2020 departure grid served 6696/6698 and
    not the headline 7835/8772, so the best lines fell back to LTE while a grid sat right
    there looking present. No grid at all is an INFO; a grid that misses the best lines
    is the WARN.
    """
    res = CheckResult(3, "NLTE grid reach (over the BEST lines)")
    grids = sorted(GRID_DIR.glob(f"{st.element}_*.csv")) + \
        sorted(GRID_DIR.glob(f"*/{st.element}_*.csv"))
    if not grids:
        res.findings.append(Finding(
            "3", INFO, st.element,
            f"no NLTE departure grid for {st.element} under {GRID_DIR.relative_to(ROOT)} "
            f"— the element runs LTE, which is a stated treatment, not a gap",
            "there is no grid to fail to reach"))
        return res

    if GRID_AVAILABILITY_CSV.exists():
        avail = pd.read_csv(GRID_AVAILABILITY_CSV)
        rows = avail[avail.element.astype(str) == st.element]
        for _, r in rows.iterrows():
            present = str(r.get("present", "")).strip().lower() == "true"
            wired = str(r.get("wired", "")).strip().lower() == "true"
            if present and not wired:
                res.findings.append(Finding(
                    "3", WARN, f"{st.element} / {r.get('grid_file','')}",
                    f"the grid is registered present but NOT wired "
                    f"({GRID_AVAILABILITY_CSV.name}: role={r.get('role','')})",
                    "we HAVE the grid; nothing consumes it, so every line silently runs "
                    "LTE with an NLTE grid on disk",
                    f"BUILD: wire the {st.element} NLTE grid {r.get('grid_file','')} — "
                    f"registered present, wired=False"))

    # THE ION KEY IS NORMALISED, NOT COMPARED AS TEXT. The grids disagree about how they
    # spell an ionisation stage — the Amarsi-2020 Al grid writes `1`, the register writes
    # `I` — and a string comparison between the two returns an empty selection, i.e. a
    # MANUFACTURED absence dressed up as "the grid does not cover this ion". Routed
    # through the project's canonical normaliser (RYA-345) so the two spellings are one.
    want_ion = parse_ion(st.ion)
    covered_waves: list[float] = []
    any_wave_col = False
    for g in grids:
        try:
            df = pd.read_csv(g)
        except Exception:
            continue
        if "wave_A" not in df.columns:
            continue
        any_wave_col = True
        sel = df
        if "ion" in df.columns:
            keep = []
            for v in df.ion:
                try:
                    keep.append(parse_ion(v) == want_ion)
                except ValueError:
                    keep.append(False)
            sel = df[pd.Series(keep, index=df.index)]
        covered_waves.extend(float(w) for w in sel.wave_A.unique())
    covered_waves = sorted(set(covered_waves))

    if not covered_waves:
        detail = ("carries a wave_A column but no row for this ion" if any_wave_col
                  else "carries no wave_A column at all")
        res.findings.append(Finding(
            "3", WARN, st.species,
            f"{len(grids)} grid file(s) exist for {st.element} and the set "
            f"{detail}: {', '.join(g.name for g in grids)}",
            "the grid file is present, so this is not a missing grid — it is a grid that "
            "serves no line of the ion being measured",
            f"AUDIT: {st.element} NLTE grid(s) present but serve no {st.species} "
            f"wavelength — confirm the ion split"))
        return res

    # TWO SUBJECTS, TWO SEVERITIES. A line ALREADY IN AN AGGREGATE that the grid misses is
    # live: the published product reads as NLTE while that line is not, which is the
    # RYA-773 Al shape exactly. A best line that has never been measured and that the grid
    # also misses is a warning about the future, not a defect in a number that exists —
    # INFO, so the WARN count stays a count of things that are wrong right now.
    is_agg = st.best_lines.origin.astype(str).str.startswith("in-aggregate")
    for subject, label, sev in (
            (st.best_lines[is_agg], "line(s) already in an aggregate", WARN),
            (st.best_lines[~is_agg], "candidate best line(s) not yet in an aggregate", INFO)):
        if not len(subject):
            continue
        hit, miss, max_delta = [], [], 0.0
        for w in subject.wavelength_air_A.astype(float):
            near, d = _nearest(w, covered_waves)
            if near is not None and d <= WAVE_TOL_A:
                hit.append((w, near, d))
                max_delta = max(max_delta, d)
            else:
                miss.append((w, near, d))
        if hit:
            res.findings.append(Finding(
                "3", OK, st.species,
                f"{len(hit)}/{len(subject)} {label} served by the grid "
                f"({', '.join(g.name for g in grids)}); largest wavelength join delta "
                f"{max_delta:.3f} A against a {WAVE_TOL_A} A tolerance",
                "reports the join delta so the match is auditable, not asserted"))
        if miss:
            shown = _by_band(miss)
            stated = _stated_engine_a_coverage(st) if sev == WARN else None
            if stated:
                # ALREADY STATED IS NOT SILENT. The product does not quietly serve these
                # lines as LTE under an NLTE label — it drops them from the Engine-A
                # aggregate and prints the count. Re-raising that as a WARN would report
                # a working guard as a defect.
                res.findings.append(Finding(
                    "3", INFO, st.species,
                    f"{len(miss)} of {len(subject)} {label} are not grid-served — "
                    f"{shown}. Not silent: {stated}",
                    "the reduced coverage is declared by the product itself, so the gap "
                    "is visible to anyone reading it — which is the whole difference "
                    "from the RYA-773 Al shape"))
                continue
            res.findings.append(Finding(
                "3", sev, st.species,
                f"{len(miss)} of {len(subject)} {label} are NOT covered by the "
                f"{st.element} grid and fall back to LTE — {shown}",
                "a grid for this element EXISTS and is consumed for other lines, and no "
                "product row declares the reduced coverage — so the product reads as "
                "NLTE while these lines are not (the RYA-773 Al shape)"
                if sev == WARN else
                "no product depends on these yet; recorded so the gap is known before "
                "they are measured, not after",
                (f"BUILD: extend the {st.element} NLTE grid to {len(miss)} line(s) that "
                 f"are IN a published aggregate and not grid-served (level-reach check "
                 f"first, RYA-773)") if sev == WARN else None))
    return res


# ── CHECK 4 — anchor consistency (the RYA-898 lesson) ────────────────────────

def check_anchor_consistency(st: State) -> CheckResult:
    """Is the solar anchor on the SAME instrument/method chain as the measurements?

    RYA-898: the [X/H] cancellation is what a differential abundance is FOR. It only
    cancels if the solar reference and the program-star measurement travelled the same
    chain. A solar anchor sitting on a literature EW pool while the star products sit on
    a different arm does not cancel anything — and nothing in the frozen record makes
    that visible, because the anchor row carries a value and a scale, not a chain.

    The discriminator: an element with no anchor yet is an INFO (nothing to be
    inconsistent with). An anchor whose DECLARED chain is disjoint from the chain the
    products were measured on is the WARN.
    """
    res = CheckResult(4, "anchor consistency (anchor chain vs product chain)")
    if st.star != "solar":
        res.findings.append(Finding(
            "4", INFO, st.star,
            "the solar anchor check applies to the solar reference itself; for a program "
            "star the relevant question is whether ITS chain matches the solar anchor's, "
            "which this check reports when run with --star sun",
            "not applicable rather than absent"))

    if not ANCHOR_VERDICT_JSON.exists():
        res.findings.append(Finding(
            "4", ERROR, "anchor channel",
            f"{ANCHOR_VERDICT_JSON.relative_to(ROOT)} is missing — the anchor's declared "
            f"chain cannot be read",
            "a missing channel is not a consistent one",
            "FIX: the phase-C verdict artifact preflight reads for anchor provenance "
            "is absent"))
        return res

    verdicts = json.loads(ANCHOR_VERDICT_JSON.read_text()).get("verdicts", [])
    row = next((v for v in verdicts if str(v.get("element")) == st.element), None)
    if row is None:
        res.findings.append(Finding(
            "4", INFO, st.element,
            f"{st.element} has no row in {ANCHOR_VERDICT_JSON.name} — no anchor has been "
            f"declared for it yet",
            "no anchor exists, so no chain can disagree with the products"))
        return res

    declared = " ".join(str(row.get(k, "")) for k in ("provenance", "channel"))
    anchor_instruments = _instruments_named_in(declared, st.instruments)
    product_instruments = _product_instruments(st)

    res.findings.append(Finding(
        "4", OK, f"{st.element} anchor",
        f"declared chain: {declared.strip()!r} (verdict {row.get('verdict')}, "
        f"n_lines={row.get('n_lines')})",
        "quotes the declaration rather than paraphrasing it"))

    if not anchor_instruments:
        res.findings.append(Finding(
            "4", WARN, f"{st.element} anchor",
            f"the anchor's declared provenance names NO instrument in "
            f"{INSTRUMENTS_CSV.name}, so whether it shares a chain with the products "
            f"cannot be established from the record",
            "we HAVE an anchor and we HAVE products; what is missing is the one field "
            "that would let the [X/H] cancellation be checked at all",
            f"AUDIT: the {st.element} solar anchor declares no instrument chain — record "
            f"it so anchor/product chain consistency is decidable (RYA-898)"))
        return res

    if not product_instruments:
        res.findings.append(Finding(
            "4", INFO, f"{st.element} anchor",
            f"anchor chain {sorted(anchor_instruments)}; no rendered product for "
            f"{st.species} under the inspected roots to compare it against",
            "scoped to the artifacts inspected — no product, so no disagreement"))
        return res

    shared = anchor_instruments & product_instruments
    if shared:
        res.findings.append(Finding(
            "4", OK, f"{st.element} anchor",
            f"anchor chain {sorted(anchor_instruments)} shares "
            f"{sorted(shared)} with the product chain {sorted(product_instruments)}",
            "the differential cancellation has a common arm"))
    else:
        res.findings.append(Finding(
            "4", WARN, f"{st.element} anchor",
            f"anchor chain {sorted(anchor_instruments)} is DISJOINT from the product "
            f"chain {sorted(product_instruments)} — a differential built on these two "
            f"cancels nothing",
            "both sides exist and both are declared; they simply are not the same "
            "measurement chain, which is invisible in a table of values",
            f"RE-DECIDE: the {st.element} solar anchor and its band products sit on "
            f"different instrument chains — the [X/H] cancellation does not hold "
            f"(RYA-898)"))
    return res


def _instruments_named_in(text: str, instruments: pd.DataFrame) -> set[str]:
    """Which catalogued instruments a provenance string names.

    Matched against the CATALOG (id and display name), never a list written here — so a
    new instrument becomes matchable by being registered, not by editing this module.
    """
    low = f" {text.lower()} "
    found = set()
    for _, r in instruments.iterrows():
        iid = str(r.instrument_id).strip().lower()
        name = str(r.get("instrument_name", "")).strip().lower()
        for tok in {iid, iid.replace("_", " "), iid.replace("_", "-"), name}:
            if tok and tok in low:
                found.add(str(r.instrument_id))
                break
    return found


def _product_instruments(st: State) -> set[str]:
    if "instrument" not in st.products.columns:
        return set()
    return {str(x) for x in st.products.instrument.dropna().unique() if str(x) != "None"}


# ── CHECK 5 — rendered-output reconciliation ─────────────────────────────────

def check_rendered_output(st: State) -> CheckResult:
    """Does the rendered product account for every verified holding and every band?

    This is the Fe-page class — the one only an eyeball caught. The registries say what we
    have; the product says what was used; nobody had ever subtracted one from the other.

    The discriminator is a STATED REASON. A verified holding absent from the product
    because no loader can read it (check 1), because the telluric gate withholds it
    (check 6), or because its wavelength range does not reach this element's lines, is
    accounted for — INFO, with the pointer. A verified holding absent for NO stated
    reason is the WARN.
    """
    res = CheckResult(5, "rendered-output reconciliation (product vs registries)")
    if not len(st.products):
        res.findings.append(Finding(
            "5", INFO, st.species,
            f"no rendered product for {st.species} under the inspected root(s) "
            f"{', '.join(str(r) for r in st.roots)} — nothing to reconcile yet",
            "scoped to these roots: it says no artifact was found here, not that no "
            "product exists"))
        return res

    product_instruments = _product_instruments(st)
    res.findings.append(Finding(
        "5", OK, st.species,
        f"{len(st.products)} product row(s) from {len(st.product_sources)} artifact(s) "
        f"on instrument(s) {sorted(product_instruments) or ['<undeclared>']}",
        "counts what the product actually contains"))

    verified = st.holdings[st.holdings.evidence_state.astype(str).isin(VERIFIED_STATES)]
    lines = st.best_lines.wavelength_air_A.astype(float).tolist()

    for _, h in verified.iterrows():
        hid, inst = str(h.holding_id), str(h.instrument_id)
        if inst in product_instruments:
            continue
        reason = _stated_reason(st, hid, inst, lines)
        if reason:
            res.findings.append(Finding(
                "5", INFO, hid,
                f"verified on {inst} and absent from the {st.species} product — "
                f"accounted for: {reason}",
                "absent WITH a stated reason is not a silent gap"))
        else:
            res.findings.append(Finding(
                "5", WARN, hid,
                f"verified on {inst}, reachable by the band harness, telluric-clear, and "
                f"within wavelength reach of {st.species} — yet it appears nowhere in the "
                f"rendered product",
                "every other explanation for its absence has been checked and none "
                "applies; the product simply does not use a holding we have",
                f"AUDIT: {st.species} product omits verified holding {hid} ({inst}) with "
                f"no stated reason — measure it or state why not"))

    if "band" in st.products.columns:
        present = {str(b) for b in st.products.band.dropna().unique()}
        expected = set()
        for w in lines:
            try:
                expected.add(band_policy.resolve(w).name)
            except Exception:
                continue
        for b in sorted(expected - present):
            ran = [r for r in st.ew_runs
                   if any(r["lo_A"] <= w <= r["hi_A"] for w in lines
                          if _band_name(w) == b)]
            if ran:
                res.findings.append(Finding(
                    "5", WARN, f"band {b}",
                    f"a run covered this band ({', '.join(Path(r['path']).name for r in ran)}) "
                    f"but no {st.species} product row reports it",
                    "the measurement exists; the rendered product drops it",
                    f"AUDIT: {st.species} {b} band was measured but carries no product row"))
            else:
                res.findings.append(Finding(
                    "5", INFO, f"band {b}",
                    f"{st.species} has best lines in {b} and no product there; no run "
                    f"covers it either",
                    "unmeasured, not dropped"))
    return res


def _line_accounting(st: State) -> list:
    """Wavelengths this element has on the RYA-709 accounting ledger, across every root.

    Read across the SAME roots as the products, because the ledger is regenerated by the
    harness and a checkout that has never run one carries the committed copy only.
    """
    waves = []
    for root in st.roots:
        p = Path(root) / LINE_ACCOUNTING_CSV.relative_to(ROOT)
        if not p.exists():
            continue
        try:
            df = pd.read_csv(p, usecols=["element", "ion", "wave_air_A"])
        except Exception:
            continue
        sel = df[(df.element.astype(str) == st.element)
                 & (df.ion.astype(str) == st.ion)]
        waves.extend(float(w) for w in sel.wave_air_A)
    return sorted(set(waves))


def _star_linelist_waves(st: State) -> list:
    """This species' wavelengths in the star's own line list.

    The third tier of "where does this line exist". A best line missing from the run AND
    from the ledger reads very differently depending on whether the star's line list has
    it at all: catalogued-and-dropped names a filter, absent-everywhere names a line-list
    coverage gap. Saying which is the difference between a finding and a guess.
    """
    p = ROOT / "data" / "linelists" / f"linelist_{st.star}.csv"
    if not p.exists():
        return []
    try:
        df = pd.read_csv(p, usecols=["element", "ion", "wavelength_air_A"])
    except Exception:
        return []
    sel = df[(df.element.astype(str) == st.element) & (df.ion.astype(str) == st.ion)]
    return sorted(set(float(w) for w in sel.wavelength_air_A))


def _stated_engine_a_coverage(st: State) -> str | None:
    """What the ENGINE-A product already says about lines the grid does not serve.

    `derive_band_products` REFUSES to average an Engine-A product over lines MPIA does not
    serve — it reports reduced coverage with the count instead (RYA-783). So the grid gap
    is, for those lines, already stated in the product. Returning that statement is what
    keeps check 3 from re-reporting a documented behaviour as a silent one; returning None
    means nothing states it, and then it IS silent.
    """
    df = st.products
    if not len(df):
        return None
    label = None
    for col in ("treatment", "product"):
        if col in df.columns:
            label = df[col].astype(str)
            break
    if label is None:
        return None
    rows = df[label.str.upper().str.contains("ENGINE-A", na=False)]
    if not len(rows):
        return None
    n_lines = pd.to_numeric(rows.get("n_lines"), errors="coerce").sum()
    n_excl = pd.to_numeric(rows.get("n_excluded"), errors="coerce").sum()
    return (f"the ENGINE-A product reports its own reduced coverage across "
            f"{len(rows)} band row(s): n_lines={n_lines:.0f}, n_excluded={n_excl:.0f}")


def _by_band(miss) -> str:
    """Summarise uncovered lines BY BAND, with one example each.

    A flat list sorted by wavelength shows six near-UV lines and hides that the gap is a
    whole regime: the Fe graded pool is 550 near-UV lines against a grid that is 251 VIS
    lines, and "3000.451, 3000.948, 3003.030 ..." reads like a handful of stragglers.
    Counting by band is what makes the shape of the hole visible.
    """
    by: dict[str, list] = {}
    for w, n, d in miss:
        by.setdefault(_band_name(w) or "<outside every declared band>", []).append((w, n, d))
    parts = []
    for band in sorted(by, key=lambda b: min(x[0] for x in by[b])):
        rows = by[band]
        w, n, d = rows[0]
        ex = (f"e.g. {w:.3f} A, nearest grid line {n:.3f} ({d:.3f} A away)"
              if n is not None else f"e.g. {w:.3f} A, grid empty")
        parts.append(f"{band}: {len(rows)} ({ex})")
    return "; ".join(parts)


def _band_name(w: float) -> str | None:
    try:
        return band_policy.resolve(w).name
    except Exception:
        return None


def _stated_reason(st: State, holding_id: str, instrument: str, lines) -> str | None:
    """Why this verified holding is legitimately absent from the product — or None.

    Every branch names the check that owns the reason, so a reader can go and read the
    finding itself rather than taking this one's word for it.
    """
    if instrument not in st.dispatch.instruments:
        return (f"the band harness has no loader for {instrument} (check 1 WARNs on it; "
                f"not double-counted here)")
    served = st.dispatch.served_holdings.get(instrument)
    if served and holding_id not in served:
        return (f"the {instrument} harness branch serves {', '.join(served)} "
                f"(check 1)")
    try:
        ok, why = telluric_policy.gate_holding(holding_id, instrument)
        if not ok:
            return f"the telluric gate withholds it — {why} (check 6)"
    except Exception as exc:
        return f"the telluric gate refuses it: {type(exc).__name__} (check 6)"
    row = st.instruments[st.instruments.instrument_id.astype(str) == instrument]
    if len(row) and lines:
        lo_nm = pd.to_numeric(row.iloc[0].get("wavelength_min_nm"), errors="coerce")
        hi_nm = pd.to_numeric(row.iloc[0].get("wavelength_max_nm"), errors="coerce")
        if pd.notna(lo_nm) and pd.notna(hi_nm):
            if not any(lo_nm * 10.0 <= w <= hi_nm * 10.0 for w in lines):
                return (f"{instrument} covers {lo_nm:.0f}-{hi_nm:.0f} nm, which reaches "
                        f"none of this element's best lines")
    return None


# ── CHECK 6 — telluric state (the correct version of the gate) ───────────────

def check_telluric_state(st: State) -> CheckResult:
    """Is telluric correction RECOGNISED as done, not merely required?

    Both directions of this were live defects. Uncorrected IR that gets measured anyway is
    the RYA-805 hole; CORRECTED IR that the gate fails to recognise and drops anyway is
    the CRIRES+ class — we paid for the correction and then threw the data away.

    The verdict is `telluric_policy.gate_holding`, called, not reimplemented. Uncorrected
    and withheld is the expected absence. Corrected but missing from the product, and
    verified-but-undetermined, are the silent gaps.
    """
    res = CheckResult(6, "telluric state (per-holding, RYA-806)")
    product_instruments = _product_instruments(st)
    lines = st.best_lines.wavelength_air_A.astype(float).tolist()

    for _, h in st.holdings.iterrows():
        hid, inst = str(h.holding_id), str(h.instrument_id)
        state = str(h.get("telluric_applied", "")).strip()
        try:
            ok, why = telluric_policy.gate_holding(hid, inst)
        except Exception as exc:
            res.findings.append(Finding(
                "6", WARN, hid,
                f"telluric_applied={state or '<blank>'} and the gate REFUSES the holding: "
                f"{type(exc).__name__}. {exc}",
                "we HOLD this product and it is verified; what blocks it is an "
                "undetermined fact about it, not a missing observation — and RYA-806 "
                "forbids defaulting it either way",
                f"BUILD: determine telluric_applied for holding {hid} ({inst}) from its "
                f"own headers/flux and record it on {HOLDINGS_CSV.name}"))
            continue

        if not ok:
            res.findings.append(Finding(
                "6", INFO, hid,
                f"telluric_applied={state} — correctly withheld: {why}",
                "we hold it uncorrected and the gate says so; withholding it is the "
                "policy working, not a gap"))
            continue

        # THE CRIRES+ CLASS, AND ONLY IT. Scoped to instruments the catalog says NEED a
        # correction: those are the arms whose data gets withheld for tellurics, so they
        # are the only ones where "we paid for the correction and the pipeline dropped it
        # anyway" is the telluric axis's own failure. A corrected arm that needs no
        # correction stage (IAG, basis 'corrected') can also be missing from a product,
        # but that is an ACCOUNTING absence and check 5 owns it — reporting it here as
        # well would double-count one gap as two.
        if (state == "applied" and telluric_policy.requires_correction(inst)
                and inst not in product_instruments and len(st.products)):
            in_reach = _reaches(st, inst, lines)
            # ONE GAP, ONE WARN. If check 1 already found that no harness path reaches
            # this holding, THAT is why it contributes nothing, and it is already
            # counted. Re-raising it here would turn one defect into two tickets.
            unreachable = _stated_reason(st, hid, inst, lines)
            if in_reach and unreachable and "check 1" in unreachable:
                res.findings.append(Finding(
                    "6", INFO, hid,
                    f"telluric_applied=applied and the gate passes it, but it contributes "
                    f"nothing to the {st.species} product for a reason the telluric axis "
                    f"does not own: {unreachable}",
                    "the telluric state is fine; the absence is a reachability defect "
                    "already counted by check 1"))
                continue
            res.findings.append(Finding(
                "6", WARN if in_reach else INFO, hid,
                f"telluric_applied=applied (the gate passes it: {why.split('.')[0]}) but "
                f"{inst} — an arm the catalog registers telluric_required=yes — "
                f"contributes nothing to the {st.species} product"
                + ("" if in_reach else ", and its wavelength range reaches none of this "
                                       "element's best lines"),
                "the correction was applied and the gate recognises it, so the data is "
                "measurable — it is simply not being measured (the CRIRES+ class)"
                if in_reach else "out of wavelength reach",
                (f"AUDIT: corrected holding {hid} ({inst}) passes the telluric gate but "
                 f"contributes no {st.species} product row") if in_reach else None))
            continue

        res.findings.append(Finding(
            "6", OK, hid, f"telluric_applied={state}: {why.split('.')[0]}",
            "gate passes and nothing is being withheld"))

    ir = [w for w in lines if w >= 9199.0]
    if ir:
        res.findings.append(Finding(
            "6", OK, st.species,
            f"{len(ir)} best line(s) at or beyond 9199 A — the regime where the telluric "
            f"axis is load-bearing",
            "states the subject the IR half of this check applies to"))
    return res


def _reaches(st: State, instrument: str, lines) -> bool:
    row = st.instruments[st.instruments.instrument_id.astype(str) == instrument]
    if not len(row) or not lines:
        return True
    lo = pd.to_numeric(row.iloc[0].get("wavelength_min_nm"), errors="coerce")
    hi = pd.to_numeric(row.iloc[0].get("wavelength_max_nm"), errors="coerce")
    if pd.isna(lo) or pd.isna(hi):
        return True
    return any(lo * 10.0 <= w <= hi * 10.0 for w in lines)


def check_normalisation_state(st: State) -> CheckResult:
    """Does the recorded `normalization_state` agree with the DECLARED pre-normalisation?

    The seventh conditioning check, and the newest axis (RYA-1030). It mirrors check 6:
    read the per-holding state the registry carries, call the module that owns the rule,
    and never re-derive either here.

    THE DEFECT IT CATCHES. A declared flag and a MIS-ROUTED file agree with each other
    perfectly and are both wrong, so no amount of cross-reading DECLARATIONS finds it --
    only the flux can. KP2005 ran that way for months: `pre_normalised=False` (RYA-929)
    while the reader opened `irradthu.dat`, absolute irradiance, and the harness fitted
    its own continuum on top, biasing A(Fe I) low by 0.022 dex with a wavelength trend.

    SEVERITIES, on this module's own expected-vs-silent split:
      INFO  the scan could not speak -- below the blue edge, inside a telluric band, or
            the product is not reachable from this machine. An expected absence.
      WARN  the registry records a state that CONTRADICTS the holding's `pre_normalised`
            flag. That is the silent gap: nothing downstream refuses on it today, and a
            continuum stage will happily run.
      WARN  the holding carries NO recorded state at all AND a continuum stage could run
            on it -- the undeclared case the ticket exists for.
    """
    res = CheckResult(7, "normalisation state (per-holding, RYA-1030)")
    declared = _declared_pre_normalised()

    if COLUMN_NORMALISATION not in st.holdings.columns:
        res.findings.append(Finding(
            "7", WARN, HOLDINGS_CSV.name,
            f"the registry has no `{COLUMN_NORMALISATION}` column, so no product's "
            f"normalisation state is recorded anywhere",
            "the column IS this check (RYA-1030); without it every holding is undeclared "
            "and a continuum stage can run on any of them unchallenged",
            "BUILD: run scripts/rya1030_backfill_normalisation.py --write"))
        return res

    for _, h in st.holdings.iterrows():
        hid = str(h.holding_id)
        state = str(h.get(COLUMN_NORMALISATION, "") or "").strip()
        flag = declared.get(hid)

        if flag is None:
            continue          # not served by the measurement harness; check 1 owns that

        if not state:
            res.findings.append(Finding(
                "7", WARN, hid,
                f"declared pre_normalised={flag} but the registry records NO scanned "
                f"normalisation state",
                "we hold this product and the harness has a reader for it, so the flux "
                "CAN be scanned -- an undeclared state is the case RYA-1030 exists for, "
                "and it is silent because a continuum stage will simply run",
                f"BUILD: scan {hid} with pipeline.normalization_intake and record "
                f"{COLUMN_NORMALISATION}"))
            continue

        if state == NI_UNKNOWN:
            res.findings.append(Finding(
                "7", INFO, hid,
                f"declared pre_normalised={flag}; the scan returned `unknown` -- the flux "
                f"could not speak (below the blue edge, inside a telluric band, or the "
                f"product is not reachable here)",
                "an expected absence: `unknown` is a real answer, not a failed scan, and "
                "the module reports it rather than defaulting either way (RYA-833)"))
            continue

        scanned_normalised = state == NI_NORMALISED
        if scanned_normalised != flag:
            res.findings.append(Finding(
                "7", WARN, hid,
                f"THE FLUX AND THE FLAG DISAGREE: declared pre_normalised={flag}, but the "
                f"recorded scan says {state}",
                "a declared flag and a mis-routed file agree with each other perfectly "
                "and are both wrong, so only the flux could have caught this -- and "
                "nothing downstream refuses on it, so a continuum stage runs anyway "
                "(the KP2005 class, RYA-929 -> RYA-933/1026)",
                f"FIX: check WHICH FILE {hid}'s reader opens before touching the flag"))
            continue

        res.findings.append(Finding(
            "7", OK, hid,
            f"{COLUMN_NORMALISATION}={state} agrees with pre_normalised={flag}",
            "the flux and the declaration were checked against each other and match"))
    return res


def _declared_pre_normalised() -> dict[str, bool]:
    """holding_id -> `pre_normalised`, parsed STATICALLY from the harness source.

    Parsed rather than imported because `measure_band_ew` resolves the Kitt Peak atlas at
    import time and `SystemExit`s when it is not staged -- importing it would make this
    check pass or fail on whether a data drive happens to be mounted.
    """
    import ast as _ast
    tree = _ast.parse((ROOT / "scripts" / "measure_band_ew.py").read_text())
    out: dict[str, bool] = {}
    for node in _ast.walk(tree):
        if isinstance(node, _ast.Call) and getattr(node.func, "id", None) == "HoldingSpec":
            for kw in node.keywords:
                if kw.arg == "pre_normalised":
                    out[node.args[0].value] = kw.value.value
    return out


# ── Report ───────────────────────────────────────────────────────────────────

CHECKS = (
    check_instrument_reachability,
    check_line_coverage,
    check_grid_reach,
    check_anchor_consistency,
    check_rendered_output,
    check_telluric_state,
    check_normalisation_state,
)


def run(star: str, element: str, ion: str, roots) -> tuple[State, list[CheckResult]]:
    st = load_state(star, element, ion, roots)
    return st, [fn(st) for fn in CHECKS]


def render(st: State, results: list[CheckResult]) -> str:
    out = []
    out.append("=" * 78)
    out.append(f"  PRE-FLIGHT READINESS — star={st.star}  element={st.species}   (RYA-905)")
    out.append("=" * 78)
    out.append(f"  advisory: this NEVER blocks. INFO = expected absence, WARN = silent gap.")
    out.append(f"  roots inspected : {', '.join(str(r) for r in st.roots)}")
    out.append(f"  holdings        : {len(st.holdings)} for {st.star} "
               f"({(st.holdings.evidence_state.astype(str).isin(VERIFIED_STATES)).sum()} verified)")
    out.append(f"  dispatch control: {st.dispatch.control_note}")
    out.append(f"  best lines      : {len(st.best_lines)} — {st.best_line_basis}")
    for res in results:
        out.append("")
        out.append(f"CHECK {res.number}. {res.name}   [{res.severity}]")
        for f in res.findings:
            out.append(f.line())

    warns = [f for r in results for f in r.findings if f.severity == WARN]
    errors = [f for r in results for f in r.findings if f.severity == ERROR]
    infos = [f for r in results for f in r.findings if f.severity == INFO]
    out.append("")
    out.append("-" * 78)
    out.append(f"SUMMARY  WARN={len(warns)}  INFO={len(infos)}  ERROR={len(errors)}   "
               f"(advisory — exit status is 0 either way)")
    if errors:
        out.append("  ERROR present: this report's absences are NOT evidence until fixed.")
    if warns:
        out.append("")
        out.append("SUGGESTED TICKETS (one line each; Ryan confirms or discards):")
        for i, f in enumerate(warns, 1):
            if f.suggested_ticket:
                out.append(f"  {i}. [check {f.check}] {f.suggested_ticket}")
                out.append(f"       why: {f.discriminator}")
    else:
        out.append("  No silent gaps found. Expected absences above are the survey's "
                   "normal state.")
    out.append("-" * 78)
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Readiness reconciliation before an element run (RYA-905). "
                    "Advisory: warns, never blocks.")
    ap.add_argument("--star", required=True, help="star id or alias, e.g. sun / solar / procyon")
    ap.add_argument("--element", required=True, help="element symbol, e.g. Fe")
    ap.add_argument("--ion", default="I", help="ionisation stage (default I)")
    ap.add_argument("--artifact-root", action="append", default=[],
                    help="extra repo-shaped root to search for measured/ and results/ "
                         "artifacts (repeatable). The default is this checkout alone, so "
                         "every absence is scoped to what was actually inspected.")
    ap.add_argument("--json", dest="json_out", default=None,
                    help="also write the findings as JSON to this path")
    a = ap.parse_args(argv)

    roots = [ROOT] + [Path(r) for r in a.artifact_root]
    st, results = run(a.star, a.element, a.ion, roots)
    print(render(st, results))

    if a.json_out:
        payload = {
            "ticket": "RYA-905",
            "star": st.star, "element": st.element, "ion": st.ion,
            "roots": [str(r) for r in st.roots],
            "dispatch_control": st.dispatch.control_note,
            "best_line_basis": st.best_line_basis,
            "checks": [{"number": r.number, "name": r.name, "severity": r.severity,
                        "findings": [asdict(f) for f in r.findings]} for r in results],
        }
        Path(a.json_out).write_text(json.dumps(payload, indent=2))
        print(f"\nwrote {a.json_out}")

    # ADVISORY BY CONSTRUCTION. A survey lacking an arm is normal, and a check that
    # fails the run on the normal state is a check somebody turns off.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

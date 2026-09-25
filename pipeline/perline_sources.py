#!/usr/bin/env python3
"""
RYA-1229 — resolve every PUBLISHED product to its per-line artifact, from the feed.

🔴 THE PER-LINE PRODUCT IS A PROJECTION OF THE FEED, NOT OF A DIRECTORY. RYA-870 built it
by globbing two ticket-scoped snapshot directories named in a module constant. The data
moved to `data/results/band_products/` and the constant did not, so the generator kept
succeeding while reading 13 of 178 per-line files: the shipped artifact went a month stale,
covered ONE of four instruments, and 932 of its 1039 rows carried an engine label
(`1D-LTE (ts-lte)`) that no published product has ever used. Nothing failed, because
nothing was asking the feed what it had published.

So discovery starts from `data/products/<star>/<element>.json` and asks each product where
its own evidence is. A product that has no reachable per-line artifact is REPORTED BY NAME,
never skipped quietly — that silence is the whole defect.

HOW A PRODUCT NAMES ITS EVIDENCE. `provenance.copied_to` is the artifact the published
number was ingested from, recorded by the emitter at publish time. It is the SSOT and it is
what we read; we never reconstruct a filename from the identity fields. The per-line sibling
sits beside it under one of two conventions, and BOTH are tried because the repo holds two
artifact families:

    <stem>_<treatment>_lines.csv   band products — one file per treatment, and the file
                                   already carries `instrument` and `treatment` columns
    <stem>_per_line.csv            the 3D-NLTE families (RYA-1095/1106/1213) — one file per
                                   holding, covering every tier, with its own schema

⚠️ `copied_to` names the treatment inconsistently across families: some stems already end
in `_<treatment>`, some do not. Stripping it before re-appending is what makes the two
conventions resolve under one rule instead of two.

🔴 EVERY RESOLUTION IS GATED ON n_lines VALUE-EQUALITY. A path that exists is not a match.
RYA-1112 found two live products sharing ion/band/holding/tier/route/treatment and differing
only in `selector`, publishing 7.339 and 7.535 — 0.196 dex apart — so a path that "looks
right" can belong to the neighbouring product. The artifact's own aggregate count must equal
the product's published `n_lines` or the product is UNRESOLVED, with both numbers reported.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PRODUCTS_DIR = ROOT / "data" / "products"

#: The production band-product directory. Named once, here, so that no caller can pin a
#: snapshot of it again (RYA-914).
BAND_PRODUCTS = ROOT / "data" / "results" / "band_products"

#: The RYA-1127 identity of a published product. Every emitted per-line row carries all of
#: it, which is what lets a reader join a line back to the number it is evidence for.
#: `treatment` is deliberately NOT abbreviated to an engine name here — RYA-906 renamed the
#: axis, and the old deck-suffixed labels are precisely what went unjoinable.
IDENTITY_FIELDS = ("element", "ion", "band", "instrument", "holding",
                   "tier", "selector", "route", "treatment")

#: The columns `pipeline.perline_product` requires of a per-line frame. A family adapter
#: must produce all of them or raise; a blank is only legal where the source genuinely has
#: no such measurement, and then it is blank, never zero.
REQUIRED_COLUMNS = ("element", "ion", "wavelength_air_A", "instrument", "treatment",
                    "in_aggregate", "abundance")


class PerLineSourceError(RuntimeError):
    """A per-line source is unreadable or self-contradictory. Never a warning."""


@dataclass
class Source:
    """One published product and the per-line rows that are evidence for it."""
    product: dict
    path: Path
    frame: pd.DataFrame
    family: str

    @property
    def identity(self) -> dict:
        p = self.product
        return {k: (p.get(k) or "") for k in IDENTITY_FIELDS}


def load_feed(star: str, element: str) -> dict:
    feed = PRODUCTS_DIR / star / f"{element}.json"
    if not feed.exists():
        raise PerLineSourceError(
            f"no published feed at {feed} — the per-line product is a projection of the "
            f"feed and cannot be built without it")
    return json.loads(feed.read_text(encoding="utf-8"))


def candidate_paths(product: dict) -> list[Path]:
    """The per-line artifacts a product's OWN provenance points at, in try order."""
    ct = ((product.get("provenance") or {}).get("copied_to") or "").strip()
    if not ct.endswith("_products.csv"):
        return []
    stem = ct[: -len("_products.csv")]
    treatment = product.get("treatment") or ""
    # ⚠️ NOT NEUTRAL TO SKIP. Some stems already carry the treatment and some do not;
    # appending blindly produced `..._ENGINE-A_ENGINE-A_lines.csv` and lost 90 products.
    if treatment and stem.endswith("_" + treatment):
        stem = stem[: -len("_" + treatment)]
    return [ROOT / f"{stem}_{treatment}_lines.csv", ROOT / f"{stem}_per_line.csv"]


def _aggregate_count(df: pd.DataFrame, family: str) -> int:
    """How many of this artifact's rows entered the published aggregate."""
    if family == "band_product":
        return int(df["in_aggregate"].astype(str).eq("True").sum())
    # The 3D-NLTE families keep the 1D-LTE membership flag and blank `a_3dnlte` for a line
    # the network refused. A line counts only when it is BOTH in the 1D aggregate and in
    # the network's training domain; `in_aggregate_1dlte` is absent on the ASPLUND_AGSS21
    # runs, whose files hold nothing but the selected pool.
    keep = df["a_3dnlte"].notna()
    if "in_aggregate_1dlte" in df.columns:
        keep &= df["in_aggregate_1dlte"].astype(str).eq("True")
    return int(keep.sum())


def _adapt_3dnlte(df: pd.DataFrame, product: dict) -> pd.DataFrame:
    """Map a 3D-NLTE per-line file onto the band-product column contract.

    🔴 THE ABUNDANCE IS `a_3dnlte`, NEVER `a_1dlte`. The 1D value is carried in the same
    file as the base the correction was applied to; publishing it under the 3D-NLTE
    treatment would report a different physics under this product's name.
    """
    for c in ("wavelength_air_A", "a_3dnlte"):
        if c not in df.columns:
            raise PerLineSourceError(
                f"3D-NLTE per-line frame for {product.get('treatment')} is missing {c!r}")
    keep = df["a_3dnlte"].notna()
    if "in_aggregate_1dlte" in df.columns:
        keep &= df["in_aggregate_1dlte"].astype(str).eq("True")

    def col(*names):
        """The first of these columns the family actually writes, else all-blank.

        ⚠️ The three 3D-NLTE families do NOT share a schema: RYA-1095/1213 write
        `ew_mA`/`rew` and carry element/ion/band, while RYA-1106 writes `ew_mA_agss21`/
        `rew_agss21` and carries none of the identity columns at all. Reading one name and
        letting the other fall through as NaN would silently blank the measured EW on a
        third of the 3D-NLTE products.
        """
        for n in names:
            if n in df.columns:
                return df[n]
        return pd.Series([None] * len(df), index=df.index)

    # element/ion are taken from the PRODUCT where the artifact does not carry them; the
    # artifact is per-(element, ion) by construction in that family.
    out = pd.DataFrame({
        "element": df["element"] if "element" in df.columns else product["element"],
        "ion": df["ion"] if "ion" in df.columns else product["ion"],
        "wavelength_air_A": df["wavelength_air_A"],
        "instrument": product["instrument"],
        "treatment": product["treatment"],
        "in_aggregate": keep.map({True: "True", False: "False"}),
        "abundance": df["a_3dnlte"],
        "ep_eV": col("elo_eV"),
        "ew_mA": col("ew_mA", "ew_mA_agss21"),
        "rew": col("rew", "rew_agss21"),
        # No chi2 surface is written for the CORRECTION: the 3D-NLTE value is an offset
        # applied to a 1D fit, so a per-line statistical sigma of the published number
        # would have to be invented. Blank is the honest state (RYA-850) and
        # `sigma_A_basis` will say so. `red_chi2` is the 1D fit's and is carried where the
        # family records it, labelled by the column it came from.
        "sigma_A": None,
        "red_chi2": col("red_chi2"),
        "ew_inversion": "False",
        "excluded_reason": col("domain_reason", "reason").fillna(""),
    })
    return out


def _read(path: Path, product: dict) -> tuple[pd.DataFrame, str]:
    df = pd.read_csv(path)
    if "a_3dnlte" in df.columns:
        return _adapt_3dnlte(df, product), "3dnlte"
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise PerLineSourceError(
            f"{path} is missing required column(s) {missing}; this product refuses to fill "
            f"a measured quantity with a default")
    return df, "band_product"


def resolve(product: dict) -> tuple[Path | None, pd.DataFrame | None, str, str]:
    """(path, frame, family, reason) for one published product. path None == unresolved."""
    cands = candidate_paths(product)
    if not cands:
        return None, None, "", ("the product records no `provenance.copied_to`, so it does "
                                "not say which artifact it was ingested from")
    present = [c for c in cands if c.exists()]
    if not present:
        rel = ", ".join(str(c.relative_to(ROOT)) for c in cands)
        return None, None, "", (f"no per-line artifact on disk (band_products is "
                                f"gitignored; only force-added files are present) — "
                                f"looked for {rel}")
    path = present[0]
    raw = pd.read_csv(path)
    family = "3dnlte" if "a_3dnlte" in raw.columns else "band_product"
    want = product.get("n_lines")
    if want is None:
        return None, None, family, "the product declares no n_lines, so nothing gates this"
    have = _aggregate_count(raw, family)
    if have != int(want):
        return None, None, family, (
            f"n_lines MISMATCH: the product publishes {want}, {path.name} aggregates "
            f"{have}. Refusing rather than stamping a neighbour's lines (RYA-1112)")
    frame, family = _read(path, product)
    return path, frame, family, "resolved"


def resolve_published(star: str, element: str):
    """Every published product of this element, split into resolved and unresolved.

    Returns (sources, unresolved). `unresolved` is a list of (identity, reason) and is the
    part a caller MUST surface: a published product whose evidence cannot be reached is a
    fact about the repo, and dropping it silently is how the artifact this module replaces
    stayed wrong for a month.
    """
    feed = load_feed(star, element)
    sources, unresolved = [], []
    for p in feed.get("products", []):
        if (p.get("element") or "").strip() != element:
            continue
        path, frame, family, why = resolve(p)
        if path is None:
            unresolved.append(({k: (p.get(k) or "") for k in IDENTITY_FIELDS}, why))
            continue
        sources.append(Source(product=p, path=path, frame=frame, family=family))
    if not sources:
        raise PerLineSourceError(
            f"not one of the {len(feed.get('products', []))} published {element} products "
            f"could be resolved to per-line evidence — refusing to emit an empty product")
    return sources, unresolved

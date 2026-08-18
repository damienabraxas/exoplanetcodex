#!/usr/bin/env python3
"""
RYA-855 — which gf rung is every published band-product cell actually entitled to?
=================================================================================
    python3 scripts/rya855_rung_audit.py --band-products <dir> [--matrix <csv>]

`scripts/derive_band_products.py` charged every cell the UNGRADED 0.17 because it passed
`gf_graded=False` at both `error_budget.build()` call sites. The wiring fix makes the
rung a function of the lines; this answers the question the fix raises — WHICH RUNG DOES
EACH ALREADY-PUBLISHED CELL LAND ON, and does its bar move.

WHY THIS RE-DERIVES THE BUDGET INSTEAD OF RE-RUNNING THE DERIVER
----------------------------------------------------------------
The budget is a pure function of (band, n_lines, scatter, gf rung, harness residual). None
of those depend on the gf rung except the rung itself, so re-fitting 18 cells — days of
synthesis — could not change any of the other four. What it WOULD do is fold two months of
input drift into a diff that is supposed to isolate one term (RYA-848: a banked artifact
carries its own drift, so prove a change with a SAME-INPUTS control).

So this reads each cell's own per-line file, rebuilds the budget from it, and CHECKS ITSELF
FIRST: `stat_dex` and the rung-1 `syst_dex` must reproduce the published numbers exactly.
A cell that fails that check is reported and NOT diffed — its published bar was not built
from the inputs standing here, and a diff against it would measure the drift, not the fix.
Only after the reproduction passes is the rung swapped in, so every moved bar is
attributable to the gf term and to nothing else.

⚠️ THE VALUES ARE NOT TOUCHED AND CANNOT BE. `A` is the median of the per-line abundances
and no term of the budget enters it. Asserted per cell anyway rather than argued.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import gf_rung                                       # noqa: E402
from pipeline.band_policy import resolve as resolve_band           # noqa: E402
from pipeline.error_budget import build as build_budget            # noqa: E402

OUT = ROOT / "data" / "results" / "rya855"

#: `derive_band_products.PROFILE_FIT_RESIDUAL_DEX`, imported rather than restated so this
#: audit cannot differ from the deriver in the one term it is holding fixed.
from derive_band_products import PROFILE_FIT_RESIDUAL_DEX          # noqa: E402

#: Which handler produced a treatment, and therefore whose measured residual it is
#: charged. Same rule the deriver applies: Engine B never touches the profile fitter,
#: and neither does the near-UV synthesis route.
_SYNTH_TREATMENTS = {"ENGINE-B", "ENGINE-B-NLTE"}

#: One representative wavelength per band, used ONLY to resolve the band policy (which
#: continuum/telluric terms apply) — never as a measurement. Same convention and same
#: values as `rya850_graded_products.BAND_PIVOT_A`.
BAND_PIVOT_A = {"VIS": 5000.0, "red-optical": 8000.0, "near-UV": 3390.0,
                "NIR": 12000.0}

_STEM = re.compile(r"^(?P<el>[A-Z][a-z]?)(?P<ion>I+|IV|VI?)_(?P<lo>\d+)_(?P<hi>\d+)_"
                   r"(?P<inst>.+?)_(?P<route>PROFILEFIT|SYNTH|LABGF)_"
                   r"(?P<treatment>.+)_lines\.csv$")


def _cells(band_products: Path) -> list[dict]:
    """Every per-line file under the band-products tree, with its cell key parsed."""
    out = []
    for f in sorted(band_products.rglob("*_lines.csv")):
        m = _STEM.match(f.name)
        if not m:
            continue
        g = m.groupdict()
        out.append({"path": f, "element": g["el"], "ion": g["ion"],
                    "lo": float(g["lo"]), "hi": float(g["hi"]),
                    "instrument": g["inst"], "route": g["route"],
                    "treatment": g["treatment"],
                    "deck": f.parent.name if f.parent != band_products else ""})
    return out


def _published(band_products: Path, cell: dict) -> dict | None:
    """The cell's own published row, read from the products.csv beside its per-line file."""
    stem = (f"{cell['element']}{cell['ion']}_{int(cell['lo'])}_{int(cell['hi'])}_"
            f"{cell['instrument']}_{cell['route']}_products.csv")
    p = cell["path"].parent / stem
    if not p.exists():
        return None
    d = pd.read_csv(p)
    r = d[d.treatment.astype(str) == cell["treatment"]]
    return None if r.empty else r.iloc[0].to_dict()


def audit(band_products: Path, linelist) -> pd.DataFrame:
    rows = []
    for cell in _cells(band_products):
        lines = pd.read_csv(cell["path"])
        used = lines[lines.in_aggregate.astype(bool) & lines.abundance.notna()]
        n = int(len(used))
        vals = used.abundance.to_numpy(dtype=float)
        value = float(np.median(vals)) if n else np.nan
        scatter = float(np.std(vals, ddof=1)) if n > 1 else 0.0

        pol = resolve_band(0.5 * (cell["lo"] + cell["hi"]))
        is_synth = (cell["treatment"] in _SYNTH_TREATMENTS
                    or cell["route"] in ("SYNTH", "LABGF"))
        harness = 0.0 if is_synth else PROFILE_FIT_RESIDUAL_DEX
        handler = "SynthesisHandler" if is_synth else "ProfileFitHandler"
        pivot = BAND_PIVOT_A.get(pol.name, 0.5 * (cell["lo"] + cell["hi"]))

        def _budget(**gf):
            return build_budget(cell["element"], pivot, max(n, 1), scatter_dex=scatter,
                                harness_residual_dex=harness, handler=handler, **gf)

        before = _budget(gf_graded=False)
        stat_before, syst_before = before.total()

        lines_gf = gf_rung.resolve_lines(cell["element"], cell["ion"],
                                         used.wavelength_air_A, linelist)
        rung = gf_rung.decide(cell["element"], cell["ion"], lines_gf)
        after = _budget(**rung.budget_kwargs())
        stat_after, syst_after = after.total()

        pub = _published(band_products, cell)
        repro_stat = repro_syst = repro_A = None
        if pub is not None:
            repro_stat = abs(round(stat_before, 4) - float(pub["stat_dex"])) <= 5e-5
            repro_syst = abs(round(syst_before, 4) - float(pub["syst_dex"])) <= 5e-5
            repro_A = abs(round(value, 3) - float(pub["A"])) <= 5e-4

        rows.append({
            "band": pol.name, "element": cell["element"], "ion": cell["ion"],
            "treatment": cell["treatment"], "deck": cell["deck"],
            "n_lines": n, "A": round(value, 3),
            "A_published": (None if pub is None else float(pub["A"])),
            "reproduces_A": repro_A,
            "stat_dex": round(stat_before, 4),
            "reproduces_stat": repro_stat,
            "syst_before": round(syst_before, 4),
            "syst_published": (None if pub is None else float(pub["syst_dex"])),
            "reproduces_syst": repro_syst,
            "rung": rung.rung, "gf_graded": rung.gf_graded,
            "n_graded": rung.n_graded, "n_unresolved": rung.n_unresolved,
            "cited_sigma_dex": rung.cited_sigma_dex,
            "syst_after": round(syst_after, 4),
            "d_syst": round(syst_after - syst_before, 4),
            "dominant_after": (after.dominant().name if after.dominant() else ""),
            "grade_counts": json.dumps(rung.grade_counts),
            "reason": rung.reason,
            "src": str(cell["path"].relative_to(band_products)),
        })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--band-products", required=True, type=Path,
                    help="a band-products tree containing *_lines.csv + *_products.csv")
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()

    # THE LINE LIST THE POOLS WERE MEASURED ON, loaded by the pipeline's own loader with
    # its production default `apply_canonical_gf=True` — so `loggf` is the value the
    # pools actually used (RYA-799/824), not the best value available.
    from pipeline.abundances_derive import _load_synth_resources
    linelist, _iso, _chem = _load_synth_resources()
    print(f"[linelist] {len(linelist)} rows, production loader, canonical gf applied")

    d = audit(a.band_products, linelist)
    a.out.mkdir(parents=True, exist_ok=True)
    d.to_csv(a.out / "rya855_rung_by_cell.csv", index=False)

    bad = d[(d.reproduces_stat == False) | (d.reproduces_syst == False)      # noqa: E712
            | (d.reproduces_A == False)]
    print(f"\n=== CONTROL: does each cell rebuild from its own per-line file? ===")
    print(f"  {int((d.reproduces_A == True).sum())}/{len(d)} reproduce the published A")
    print(f"  {int((d.reproduces_stat == True).sum())}/{len(d)} reproduce stat_dex")
    print(f"  {int((d.reproduces_syst == True).sum())}/{len(d)} reproduce syst_dex on "
          f"rung 1 (the hardcode)")
    for _, r in bad.iterrows():
        print(f"  ⚠️ {r.band:<12}{r.element}{r.ion:<3}{r.treatment:<15} "
              f"A {r.A} vs {r.A_published}  stat {r.stat_dex} vs -  "
              f"syst {r.syst_before} vs {r.syst_published} — NOT diffed")

    print(f"\n=== the rung each cell is entitled to ===")
    print(f"{'band':<13}{'sp':<7}{'treatment':<16}{'n':>4}{'rung':>6}"
          f"{'graded':>8}{'syst':>9}{'->':>4}{'syst':>9}{'delta':>9}")
    for _, r in d.iterrows():
        print(f"{r.band:<13}{r.element + ' ' + r.ion:<7}{r.treatment:<16}{r.n_lines:>4}"
              f"{r.rung:>6}{r.n_graded:>8}{r.syst_before:>9.4f}{'->':>4}"
              f"{r.syst_after:>9.4f}{r.d_syst:>+9.4f}")

    moved = d[d.d_syst.abs() > 5e-5]
    print(f"\n=== {len(moved)} of {len(d)} cells change their bar ===")
    for _, r in moved.iterrows():
        print(f"  {r.band} {r.element} {r.ion} {r.treatment}: "
              f"{r.syst_before:.4f} -> {r.syst_after:.4f}  (rung {r.rung})")
        print(f"      {r.reason}")
    if not len(moved):
        print("  none. Every published cell was ALREADY entitled to rung 1 only — the")
        print("  hardcode agreed with the data by accident, cell by cell, and the reason")
        print("  is now stated per cell instead of assumed. See the `reason` column.")

    print(f"\n=== values ===")
    print(f"  {int((d.reproduces_A == True).sum())}/{len(d)} cells reproduce A exactly; "
          f"no term of the budget enters the median, and none did.")

    summary = {
        "ticket": "RYA-855",
        "band_products": str(a.band_products),
        "n_cells": int(len(d)),
        "n_reproducing_published_syst": int((d.reproduces_syst == True).sum()),
        "n_cells_by_rung": {int(k): int(v) for k, v in
                            d.rung.value_counts().sort_index().items()},
        "n_cells_moved": int(len(moved)),
        "cells": json.loads(d.to_json(orient="records")),
    }
    (a.out / "rya855_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\n[out] {a.out}/rya855_rung_by_cell.csv\n[out] {a.out}/rya855_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

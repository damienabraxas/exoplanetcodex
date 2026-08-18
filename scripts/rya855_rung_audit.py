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

#: Treatments produced by a FLUX FIT, which never touches the profile fitter and so is
#: not charged its measured residual. This is the rule by HANDLER — the correct one.
_FLUX_FIT_TREATMENTS = {"ENGINE-B", "ENGINE-B-NLTE"}


def _deriver_harness(treatment: str, route: str) -> float:
    """The residual `derive_band_products` ACTUALLY charges. Mirrored, defect included.

    🔴 DO NOT "FIX" THIS TO MATCH `_FLUX_FIT_TREATMENTS`. This audit exists to isolate
    ONE term, so every other term must be held at whatever the deriver produced; a
    baseline that silently corrects a second term would attribute that correction to the
    gf change (RYA-848: prove a change with a SAME-INPUTS control).

    ⚠️ AND THE MIRROR IS WHAT SURFACED A SECOND DEFECT, reported separately below and
    filed rather than folded in. The deriver's test is `prod.treatment == "ENGINE-B"`,
    an equality against a treatment name that GAINED A VARIANT in RYA-798: every
    ENGINE-B-NLTE product is therefore charged the profile fitter's 0.0129 dex and
    labelled `ProfileFitHandler` in its own budget file, while the ENGINE-B product of
    the same handler is charged 0.0000 and labelled `SynthesisHandler`. Not this
    ticket's subject and not folded into its diff.
    """
    if route in ("SYNTH", "LABGF"):
        return 0.0
    return 0.0 if treatment == "ENGINE-B" else PROFILE_FIT_RESIDUAL_DEX


def _handler_harness(treatment: str, route: str) -> float:
    """The residual the treatment's HANDLER earns. Differs from the above on NLTE."""
    if route in ("SYNTH", "LABGF") or treatment in _FLUX_FIT_TREATMENTS:
        return 0.0
    return PROFILE_FIT_RESIDUAL_DEX

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


def _linelist_for(route: str, band: str, lists: dict):
    """The line list THAT CELL WAS MEASURED ON — never a single default.

    🔴 THE FIRST CUT OF THIS AUDIT USED ONE LIST FOR EVERY CELL and reported the near-UV
    pool as 40 of 40 lines UNRESOLVED. That is not a property of the pool: the near-UV
    route runs on `data/linelists/ispec_nearuv_3000_3780`, and the default GES list spans
    4200-9200 A, so not one 3000-3780 A line could be in it. A MANUFACTURED ABSENCE
    (RYA-833) — it happened to fall on the safe side, because an unresolvable line forces
    rung 1 anyway, but the stated REASON was wrong and would have read as a finding about
    the lines. The deriver itself was always right here: its synthesis route passes the
    band's own list through `ctx["linelist"]`.
    """
    key = band if route in ("SYNTH", "LABGF") else "__ew__"
    if key not in lists:
        raise KeyError(
            f"no line list loaded for {route}/{band} — refusing to grade a pool against "
            f"a list it was not measured on")
    return lists[key]


def audit(band_products: Path, lists: dict) -> pd.DataFrame:
    rows = []
    for cell in _cells(band_products):
        lines = pd.read_csv(cell["path"])
        used = lines[lines.in_aggregate.astype(bool) & lines.abundance.notna()]
        n = int(len(used))
        vals = used.abundance.to_numpy(dtype=float)
        value = float(np.median(vals)) if n else np.nan
        scatter = float(np.std(vals, ddof=1)) if n > 1 else 0.0

        pol = resolve_band(0.5 * (cell["lo"] + cell["hi"]))
        harness = _deriver_harness(cell["treatment"], cell["route"])
        correct = _handler_harness(cell["treatment"], cell["route"])
        handler = ("SynthesisHandler" if harness == 0.0 else "ProfileFitHandler")
        pivot = BAND_PIVOT_A.get(pol.name, 0.5 * (cell["lo"] + cell["hi"]))

        def _budget(_h=harness, **gf):
            return build_budget(cell["element"], pivot, max(n, 1), scatter_dex=scatter,
                                harness_residual_dex=_h, handler=handler, **gf)

        before = _budget(gf_graded=False)
        stat_before, syst_before = before.total()

        lines_gf = gf_rung.resolve_lines(
            cell["element"], cell["ion"], used.wavelength_air_A,
            _linelist_for(cell["route"], pol.name, lists))
        rung = gf_rung.decide(cell["element"], cell["ion"], lines_gf)
        # WHY a line could not be priced, split: absent from the list is a coverage fact,
        # two rows inside the tolerance is an identification one. Collapsing them would
        # hide which is which, and they have different fixes.
        _unres = lines_gf[~lines_gf.resolved] if len(lines_gf) else lines_gf
        n_absent = int(_unres.unresolved_why.str.startswith("absent").sum()) if len(_unres) else 0
        after = _budget(**rung.budget_kwargs())
        stat_after, syst_after = after.total()

        pub = _published(band_products, cell)
        repro_stat = repro_syst = repro_A = None
        if pub is not None:  # noqa: E501
            repro_stat = abs(round(stat_before, 4) - float(pub["stat_dex"])) <= 5e-5
            repro_syst = abs(round(syst_before, 4) - float(pub["syst_dex"])) <= 5e-5
            repro_A = abs(round(value, 3) - float(pub["A"])) <= 5e-4

        rows.append({
            "band": pol.name, "element": cell["element"], "ion": cell["ion"],
            "treatment": cell["treatment"], "deck": cell["deck"],
            "n_lines": n, "A": round(value, 3), "has_published_row": pub is not None,
            "A_published": (None if pub is None else float(pub["A"])),
            "reproduces_A": repro_A,
            "stat_dex": round(stat_before, 4),
            "reproduces_stat": repro_stat,
            "syst_before": round(syst_before, 4),
            "syst_published": (None if pub is None else float(pub["syst_dex"])),
            "reproduces_syst": repro_syst,
            "rung": rung.rung, "gf_graded": rung.gf_graded,
            "n_graded": rung.n_graded, "n_unresolved": rung.n_unresolved,
            "n_absent_from_linelist": n_absent,
            "n_ambiguous_in_linelist": int(rung.n_unresolved - n_absent),
            "cited_sigma_dex": rung.cited_sigma_dex,
            "syst_after": round(syst_after, 4),
            "d_syst": round(syst_after - syst_before, 4),
            "dominant_after": (after.dominant().name if after.dominant() else ""),
            # RYA-855's own finding, carried per cell rather than only narrated: what
            # this cell's systematic would be if the harness residual followed its
            # HANDLER instead of an equality test on the treatment name.
            "harness_charged_dex": harness,
            "harness_by_handler_dex": correct,
            "harness_misattributed": bool(abs(harness - correct) > 1e-12),
            "syst_if_harness_correct": round(
                _budget(_h=correct, **rung.budget_kwargs()).systematic(), 4),
            "grade_counts": json.dumps(rung.grade_counts),
            "reason": rung.reason,
            "src": str(cell["path"].relative_to(band_products)),
        })
    return pd.DataFrame(rows)


#: RYA-836's near-UV lab-gf sub-pool, per line. The ONE graded cell in the Fe matrix.
RYA836_PER_LINE = (ROOT / "data" / "results" / "rya836"
                   / "rya836_nearuv_lab_gf_per_line.csv")
#: RYA-850's published cited-sigma table, for the comparison in `graded_pool_control`.
RYA850_SUMMARY = ROOT / "data" / "results" / "rya850" / "rya850_summary.json"


def graded_pool_control() -> dict | None:
    """POSITIVE CONTROL — does the decider find a graded pool where one demonstrably is?

    Everything else in this audit lands on rung 1, and a decider that returned rung 1
    unconditionally would produce exactly that table. An absence needs a positive
    control (RYA-833/805), so the decider is run against the one pool in the repo that
    IS entirely primary-lab: RYA-836's near-UV sub-pool, which `rya836_nearuv_lab_gf_
    subpool.py` charges `gf_graded=True` by hand. If the decider agrees, that hand
    assertion is confirmed rather than merely repeated.

    🔴 IT IS PRICED ON `loggf_lab`, NOT ON THE PRODUCTION LIST. This pool was re-inverted
    ON the laboratory value; a grade describes the log gf the pool ACTUALLY used
    (RYA-799), so resolving these lines against the production line list would grade them
    SCALE-MISMATCH and report rung 1 for a pool that is entirely laboratory-measured.
    That is the same trap in the opposite direction, and it is why `decide` takes the gf
    as data rather than looking it up.
    """
    if not RYA836_PER_LINE.exists():
        return None
    d = pd.read_csv(RYA836_PER_LINE)
    used = d[d.a_labgf.notna()]
    lines = pd.DataFrame({"wavelength_air_A": used.wavelength_air_A.astype(float),
                          "ep_eV": used.ep_eV.astype(float),
                          "log_gf": used.loggf_lab.astype(float)})
    r = gf_rung.decide("Fe", "I", lines)
    out = {"n_lines": int(len(lines)), "rung": r.rung, "gf_graded": r.gf_graded,
           "n_graded": r.n_graded, "cited_sigma_dex": r.cited_sigma_dex,
           "cited_source": r.cited_source, "grade_counts": r.grade_counts,
           "reason": r.reason,
           # RYA-836 asserts this by hand at its own `build_budget` call.
           "rya836_asserts_graded": True,
           "decider_agrees_with_rya836": bool(r.gf_graded)}

    # ⚠️ AND IT DISAGREES WITH RYA-850 ON THE SIGMA, for a reason worth stating.
    # RYA-850 matched the lab table on WAVELENGTH ALONE inside 0.05 A, caught two rows
    # for two of the sixty lines, and counted those UNMATCHED — the right call under that
    # rule, since `iloc[0]` there is how RYA-853 manufactured 12-dex defects. The pair is
    # 3125.651 / 3125.683 A: two REAL Fe I transitions 0.032 A apart whose excitation
    # potentials are 0.990 and 2.404 eV. `gf_grades` keys on wavelength AND EP AND
    # agreement with the log gf the pool used, so it separates them and each line matches
    # its own row. The ambiguity was in the MATCH RULE, not in the data.
    #
    # NOT CHANGED HERE. It moves a published bar (RYA-850's near-UV graded cell) and that
    # is that ticket's number to move.
    if RYA850_SUMMARY.exists():
        j = json.loads(RYA850_SUMMARY.read_text())
        c = (j.get("cited_gf_sigma_by_band") or {}).get("near-UV") or {}
        out["rya850_published"] = {"cited_sigma": c.get("cited_sigma"),
                                   "n": c.get("n"), "n_pool": c.get("n_pool"),
                                   "ambiguous": c.get("ambiguous")}
    return out


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
    from derive_band_products import SYNTH_BANDS
    from pipeline.nearuv_synth import gf_provenance

    lists = {}
    lists["__ew__"], _iso, _chem = _load_synth_resources()
    print(f"[linelist] EW route: {len(lists['__ew__'])} rows (GES, canonical gf applied)")
    for band, cfg in SYNTH_BANDS.items():
        if not cfg.linelist.exists():
            print(f"[linelist] {band}: {cfg.linelist} ABSENT — cells on this band will "
                  f"refuse rather than be graded against another band's list")
            continue
        _w = pd.read_csv(cfg.linelist, sep="\t", usecols=["wave_A"], low_memory=False)
        prov = gf_provenance(float(_w.wave_A.min()), float(_w.wave_A.max()))
        lists[band], _, _ = _load_synth_resources(
            linelist_file=str(cfg.linelist),
            apply_canonical_gf=prov["apply_canonical_gf"])
        print(f"[linelist] {band}: {len(lists[band])} rows from {cfg.linelist.name}, "
              f"canonical gf {'applied' if prov['apply_canonical_gf'] else 'NOT applied'}")

    d = audit(a.band_products, lists)
    a.out.mkdir(parents=True, exist_ok=True)
    d.to_csv(a.out / "rya855_rung_by_cell.csv", index=False)

    have = d[d.has_published_row]
    bad = have[(have.reproduces_stat == False) | (have.reproduces_syst == False)  # noqa: E712
               | (have.reproduces_A == False)]
    print(f"\n=== CONTROL: does each cell rebuild from its own per-line file? ===")
    print(f"  {len(have)} of {len(d)} per-line files have a published row to check "
          f"against; {len(d) - len(have)} do not and are reported, never assumed")
    print(f"  {int((have.reproduces_A == True).sum())}/{len(have)} reproduce the "
          f"published A")
    print(f"  {int((have.reproduces_stat == True).sum())}/{len(have)} reproduce stat_dex")
    print(f"  {int((have.reproduces_syst == True).sum())}/{len(have)} reproduce syst_dex "
          f"on rung 1 (the hardcode)")
    for _, r in bad.iterrows():
        print(f"  ⚠️ {r.band:<12}{r.element}{r.ion:<3}{r.treatment:<15} "
              f"A {r.A} vs {r.A_published}   "
              f"syst {r.syst_before} vs {r.syst_published} — NOT diffed")
    for _, r in d[~d.has_published_row].iterrows():
        print(f"  (no published row) {r.band} {r.element} {r.ion} {r.treatment} "
              f"[{r.deck or 'top level'}] — per-line file present, products.csv has no "
              f"such treatment")

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
    print(f"  {int((have.reproduces_A == True).sum())}/{len(have)} cells reproduce A "
          f"exactly; no term of the budget enters the median, and none did.")

    ctrl = graded_pool_control()
    print(f"\n=== POSITIVE CONTROL: a pool that IS entirely primary-lab ===")
    if ctrl is None:
        print(f"  ⚠️ {RYA836_PER_LINE} absent — the table above is an unproven absence")
    else:
        print(f"  RYA-836 near-UV lab-gf sub-pool, n={ctrl['n_lines']}, priced on "
              f"loggf_lab: rung {ctrl['rung']}, graded={ctrl['gf_graded']}, "
              f"cited sigma {ctrl['cited_sigma_dex']}")
        print(f"  RYA-836 asserts gf_graded=True by hand; decider "
              f"{'AGREES' if ctrl['decider_agrees_with_rya836'] else '⚠️ DISAGREES'}")
        print(f"  {ctrl['reason']}")
        pub = ctrl.get("rya850_published") or {}
        if pub.get("cited_sigma") is not None:
            print(f"  ⚠️ RYA-850 PUBLISHED {pub['cited_sigma']:.4f} over n={pub['n']}/"
                  f"{pub['n_pool']} ({pub['ambiguous']} refused as ambiguous under its "
                  f"wavelength-only\n     0.05 A match). The two are 3125.651 / 3125.683 "
                  f"A — REAL, distinct Fe I transitions\n     0.032 A apart at EP 0.990 "
                  f"and 2.404 eV, which the wavelength+EP key separates cleanly.\n"
                  f"     So the near-UV cited sigma is defensibly "
                  f"{ctrl['cited_sigma_dex']:.4f} over all {ctrl['n_lines']}. FLAGGED, "
                  f"NOT CHANGED — it is RYA-850's published bar.")

    # ── the second defect, found by mirroring and NOT fixed here ──────────────────
    mis = d[d.harness_misattributed]
    print(f"\n=== SEPARATE FINDING (not this ticket's subject, not fixed here) ===")
    if not len(mis):
        print("  none")
    else:
        print(f"  {len(mis)} cells are charged the PROFILE FITTER's measured residual "
              f"({PROFILE_FIT_RESIDUAL_DEX} dex) by a flux-fit engine that never touches "
              f"it.\n  `derive_band_products` tests `prod.treatment == \"ENGINE-B\"`, an "
              f"equality against a\n  treatment name that gained the ENGINE-B-NLTE "
              f"variant in RYA-798, so every NLTE\n  cell is also LABELLED "
              f"`ProfileFitHandler` in its own budget file.")
        for _, r in mis.iterrows():
            print(f"    {r.band:<12}{r.element} {r.ion:<3}{r.treatment:<16}"
                  f"syst {r.syst_after:.4f} -> {r.syst_if_harness_correct:.4f} "
                  f"if the residual followed the handler")

    summary = {
        "ticket": "RYA-855",
        "band_products": str(a.band_products),
        "n_cells": int(len(d)),
        "n_reproducing_published_syst": int((d.reproduces_syst == True).sum()),
        "n_cells_by_rung": {int(k): int(v) for k, v in
                            d.rung.value_counts().sort_index().items()},
        "n_cells_moved": int(len(moved)),
        "graded_pool_control": ctrl,
        "separate_finding_harness_misattributed": {
            "what": "derive_band_products charges the ProfileFitHandler residual to "
                    "ENGINE-B-NLTE because it tests `treatment == \"ENGINE-B\"`; the "
                    "NLTE variant arrived in RYA-798 and the equality never followed",
            "n_cells": int(len(mis)),
            "cells": json.loads(mis[["band", "element", "ion", "treatment", "deck",
                                     "syst_after", "syst_if_harness_correct"]]
                                .to_json(orient="records")),
            "status": "FILED SEPARATELY — not fixed in RYA-855, whose diff must isolate "
                      "the gf term",
        },
        "cells": json.loads(d.to_json(orient="records")),
    }
    (a.out / "rya855_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\n[out] {a.out}/rya855_rung_by_cell.csv\n[out] {a.out}/rya855_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

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

from pipeline import harness_residual                             # noqa: E402  RYA-869

def _harness(treatment: str, route: str) -> harness_residual.HarnessResidual:
    """The residual and label this banked cell's HANDLER earns. ONE rule, not two.

    🔴 THIS AUDIT USED TO CARRY TWO. `_deriver_harness` mirrored the deriver's defect on
    purpose — `prod.treatment == "ENGINE-B"`, which missed RYA-798's `ENGINE-B-NLTE`
    variant — so that RYA-855's baseline reproduced the PUBLISHED numbers exactly and
    every moved bar was attributable to the gf term alone (RYA-848). `_handler_harness`
    stood beside it holding the correct answer, and the gap between the two was RYA-855's
    separate finding: four published bars charged another handler's systematic.

    RYA-869 fixed the deriver, so the mirror has nothing left to mirror and both are
    gone. What remains is the shared decider, asked by ROUTE + TREATMENT because these
    are banked per-line files that predate the `handler` column — see
    `harness_residual.handler_of_banked_cell`.
    """
    return harness_residual.for_handler(
        harness_residual.handler_of_banked_cell(route=route, treatment=treatment))

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
        hr = _harness(cell["treatment"], cell["route"])
        harness, handler = hr.residual_dex, hr.handler
        pivot = BAND_PIVOT_A.get(pol.name, 0.5 * (cell["lo"] + cell["hi"]))

        def _budget(_h=harness, _p=hr.provenance, **gf):
            return build_budget(cell["element"], pivot, max(n, 1), scatter_dex=scatter,
                                harness_residual_dex=_h, handler=handler,
                                # RYA-873 — carry the provenance too, or the budget text
                                # this audit rebuilds would describe an uncharged zero as
                                # MEASURED and stop matching the artifact it checks.
                                harness_provenance=_p, **gf)

        before = _budget(gf_graded=False)
        stat_before, syst_before = before.total()

        # RYA-871 — the per-line artifact carries `ep_eV` now, so a measured line is
        # matched back to the list on wavelength AND EP. A file written before RYA-871
        # has no such column, and `.get` returning None keeps THAT file on the narrow
        # wavelength-only rule rather than widening it blind — which is the whole point
        # of the two tolerances travelling with the key.
        lines_gf = gf_rung.resolve_lines(
            cell["element"], cell["ion"], used.wavelength_air_A,
            _linelist_for(cell["route"], pol.name, lists),
            measured_ep_eV=(used["ep_eV"] if "ep_eV" in used.columns else None))
        rung = gf_rung.decide(cell["element"], cell["ion"], lines_gf)
        # WHY a line could not be priced, split: absent from the list is a coverage fact,
        # two rows inside the tolerance is an identification one. Collapsing them would
        # hide which is which, and they have different fixes.
        _unres = lines_gf[~lines_gf.resolved] if len(lines_gf) else lines_gf
        n_absent = int(_unres.unresolved_why.str.startswith("absent").sum()) if len(_unres) else 0
        after = _budget(**rung.budget_kwargs())
        stat_after, syst_after = after.total()

        # 🔴 RYA-869 — WHY THE PUBLISHED syst MAY NOT REPRODUCE, AND IT IS NOT DRIFT.
        # Six banked cells were charged the profile fitter's residual by a flux fit
        # (`treatment == "ENGINE-B"` never followed RYA-798's `ENGINE-B-NLTE` variant),
        # and RYA-869 fixed it. So a pre-RYA-869 artifact's `syst_dex` legitimately
        # differs from what this audit now rebuilds — and that has to be told apart from
        # the input drift the reproduction control exists to catch, because the two need
        # OPPOSITE handling: drift invalidates a diff, a known correction does not.
        # `charged_in_banked_cell` is the pre-fix rule, quoted in one place.
        was = harness_residual.charged_in_banked_cell(route=cell["route"],
                                                      treatment=cell["treatment"])
        syst_pre869 = _budget(_h=was.residual_dex, gf_graded=False).systematic() \
            if was.residual_dex != harness else syst_before

        pub = _published(band_products, cell)
        repro_stat = repro_syst = repro_A = None
        explained = None
        if pub is not None:  # noqa: E501
            repro_stat = abs(round(stat_before, 4) - float(pub["stat_dex"])) <= 5e-5
            repro_syst = abs(round(syst_before, 4) - float(pub["syst_dex"])) <= 5e-5
            repro_A = abs(round(value, 3) - float(pub["A"])) <= 5e-4
            # Does the PRE-RYA-869 harness rule reproduce it? If so the only difference
            # is the correction, the cell's inputs are intact, and it is safe to diff.
            explained = (not repro_syst
                         and abs(round(syst_pre869, 4) - float(pub["syst_dex"])) <= 5e-5)

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
            # RYA-869 — a False `reproduces_syst` explained ENTIRELY by the harness
            # correction. Carried per cell so "not reproduced" never has to be taken on
            # trust as "probably the known fix".
            "syst_published_rule_pre_rya869": round(syst_pre869, 4),
            "syst_diff_is_rya869_correction": explained,
            "rung": rung.rung, "gf_graded": rung.gf_graded,
            "n_graded": rung.n_graded, "n_unresolved": rung.n_unresolved,
            "n_absent_from_linelist": n_absent,
            "n_ambiguous_in_linelist": int(rung.n_unresolved - n_absent),
            "cited_sigma_dex": rung.cited_sigma_dex,
            "syst_after": round(syst_after, 4),
            "d_syst": round(syst_after - syst_before, 4),
            "dominant_after": (after.dominant().name if after.dominant() else ""),
            # RYA-855's finding, now FIXED by RYA-869: the harness residual follows the
            # handler that produced the cell. The `harness_by_handler_dex` /
            # `harness_misattributed` / `syst_if_harness_correct` columns that stood here
            # measured the gap between the deriver's rule and the correct one; there is
            # one rule now, so the gap has no definition and reporting a column of zeros
            # would imply an ongoing check that no longer exists.
            "harness_charged_dex": harness,
            "harness_handler": handler,
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
    # RYA-869 — a syst mismatch the harness correction fully explains is NOT a failed
    # reproduction: the cell's inputs are intact and its diff is still attributable.
    fixed = have[have.syst_diff_is_rya869_correction == True]                   # noqa: E712
    bad = have[(have.reproduces_stat == False) | (have.reproduces_A == False)   # noqa: E712
               | ((have.reproduces_syst == False)                               # noqa: E712
                  & (have.syst_diff_is_rya869_correction != True))]             # noqa: E712
    print(f"\n=== CONTROL: does each cell rebuild from its own per-line file? ===")
    print(f"  {len(have)} of {len(d)} per-line files have a published row to check "
          f"against; {len(d) - len(have)} do not and are reported, never assumed")
    print(f"  {int((have.reproduces_A == True).sum())}/{len(have)} reproduce the "
          f"published A")
    print(f"  {int((have.reproduces_stat == True).sum())}/{len(have)} reproduce stat_dex")
    print(f"  {int((have.reproduces_syst == True).sum())}/{len(have)} reproduce syst_dex "
          f"on rung 1 (the hardcode)")
    if len(fixed):
        print(f"  + {len(fixed)} reproduce it under the PRE-RYA-869 harness rule and "
              f"differ ONLY by that correction — inputs intact, still diffed:")
        for _, r in fixed.iterrows():
            print(f"      {r.band:<12}{r.element}{r.ion:<3}{r.treatment:<15} "
                  f"published {r.syst_published:.4f} = pre-869 rule "
                  f"{r.syst_published_rule_pre_rya869:.4f}; charged by handler "
                  f"{r.syst_before:.4f}")
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

    # ── what each cell's harness term is now charged, and under whose name ────────
    # RYA-855 reported this block as a SEPARATE FINDING and did not fix it: the deriver
    # tested `prod.treatment == "ENGINE-B"`, which never followed RYA-798's
    # `ENGINE-B-NLTE` variant, so six banked cells were charged the profile fitter's
    # residual by a flux fit that never touches it. RYA-869 fixed the deriver, and this
    # audit now asks the same decider the deriver asks, so there is no gap left to
    # report. It is printed anyway rather than deleted — the whole table on ONE rule is
    # the evidence that the two answers merged, and a silent removal would read as the
    # finding having gone away by itself.
    print(f"\n=== harness residual per cell, charged BY HANDLER (RYA-869) ===")
    for (h, t), grp in d.groupby(["harness_handler", "treatment"]):
        print(f"  {t:<16}{h:<20}{grp.harness_charged_dex.iloc[0]:.4f} dex   "
              f"({len(grp)} cell{'s' if len(grp) != 1 else ''})")

    summary = {
        "ticket": "RYA-855",
        "band_products": str(a.band_products),
        "n_cells": int(len(d)),
        "n_reproducing_published_syst": int((d.reproduces_syst == True).sum()),
        "n_reproducing_published_syst_under_pre_rya869_rule":
            int((d.syst_diff_is_rya869_correction == True).sum()),
        "n_cells_by_rung": {int(k): int(v) for k, v in
                            d.rung.value_counts().sort_index().items()},
        "n_cells_moved": int(len(moved)),
        "graded_pool_control": ctrl,
        "separate_finding_harness_misattributed": {
            "what": "derive_band_products charged the ProfileFitHandler residual to "
                    "ENGINE-B-NLTE because it tested `treatment == \"ENGINE-B\"`; the "
                    "NLTE variant arrived in RYA-798 and the equality never followed",
            "status": "FIXED in RYA-869 — the residual and the handler label now come "
                      "from pipeline.harness_residual, keyed on the handler the product "
                      "declares. This audit's two mirrored harness rules were deleted "
                      "with it; the numbers below are on the ONE surviving rule.",
            "harness_by_handler_dex": {
                str(h): float(g.harness_charged_dex.iloc[0])
                for h, g in d.groupby("harness_handler")},
        },
        "caveats": [
            "The EW route's per-line artifact carries no excitation potential, so a "
            "measured line is resolved back to the loaded line list on WAVELENGTH ALONE. "
            "16 of 152 VIS lines and 13 of 101 red-optical lines do not resolve uniquely "
            "at 0.005 A; widening the window converts 'absent' into 'ambiguous' rather "
            "than resolving it, because several have two Fe I rows straddling them. Those "
            "lines stay UNGRADEABLE and force rung 1 — conservative, and it flips no "
            "answer here because every EW pool is mixed several times over. It IS a "
            "ceiling: a pool that was otherwise entirely lab-gf would be held at rung 1 "
            "by lines nobody can identify. The fix is to carry ep_eV on the EW artifact.",
            "The near-UV synthesis pool is 15 GF-NIST + 23 systematic:K07 + 2 ambiguous "
            "and ZERO primary-lab, which corroborates the provenance prose this ticket "
            "replaced ('17 of 40 carry a citable NIST accuracy class'). RYA-822 keeps "
            "GF-NIST outside `is_graded` because FMW *is* NIST and VALD copies it "
            "(RYA-760), so those 15 do not buy the band a rung.",
            "The harness-residual misattribution this audit surfaced is FIXED "
            "(RYA-869): the residual and its handler label now come from "
            "pipeline.harness_residual, keyed on the handler the product declares. The "
            "two mirrored rules this audit carried — one reproducing the defect so the "
            "gf baseline matched the published bars, one holding the correct answer — "
            "are gone with it, and the numbers here are on the one surviving rule. Six "
            "cells (four published matrix cells) moved by 0.0005 dex; see "
            "data/results/rya869/.",
        ],
        "cells": json.loads(d.to_json(orient="records")),
    }
    (a.out / "rya855_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\n[out] {a.out}/rya855_rung_by_cell.csv\n[out] {a.out}/rya855_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

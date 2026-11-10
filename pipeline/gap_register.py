"""Atomic Data Gap Register (RYA-1378): every line of an element the grading keeps out of the
reference pool, as a list for future laboratory or theoretical study.

Ryan 2026-11-09: "all lines that are held or whatever due to grading should be added to a
list for future study". One generator for every element (`run_pipeline.py --element El
--gap-register`, and at the end of every element run); the Cr template was RYA-728's
scripts/rya728_cr_gap_candidates.py. Nothing here changes a gf, a grade, a pool or a product.

Population: every canonical_gf row of the element (I, II; inside a band) that is NOT a reference (LAB-tier) line,
plus every curated lab line with no canonical_gf row at all. The pools use the ORCH's own
definitions (run_matrix._canonical_species): reference = LAB tier; graded = any published gf
uncertainty (pipeline.gf_grades.is_gf_graded); a registry-excluded line is in neither.

held_reason (why the grading keeps the line out of the reference pool; one per row):
  REGISTRY_EXCLUDED  the line registry excludes it for this star
  UNGRADED_GF        no published gf uncertainty: in no graded pool
  GRADED_NOT_LAB     graded (stored sigma / NIST A-C) but no laboratory gf: codex pool only
  NOT_IN_CANONICAL   a laboratory line with no canonical_gf row

Flags (RYA-1378 vocabulary; may coexist):
  UNRATED_GF         no primary laboratory gf attached
  MISSING_GF         no canonical_gf row at all
  PROVENANCE_GAP     gf from a compilation (MFW, CSE, PhFl) or an opaque VALD3 tag
  GF_DISCREPANCY     curated: a lab gf and the list gf disagree (data/reference/atomic_gap_known.csv)
  SOLAR_DISCREPANCY  curated, or a band diagnostic fit > SOLAR_DEX from the litscan solar
                     value / non-physical. Fit evidence is only valid from runs on or after
                     01f51846 (RYA-1350 continuum gate); regenerate after a rerun.

Priority (higher = more useful to a laboratory or theoretical effort):
  1  the species has NO graded line anywhere in the band: every cell there is HELD
  2  a GF_DISCREPANCY or SOLAR_DISCREPANCY line in a band with a graded route
  3  in the synthesis list   4  not in it
Priority 1 does not key on in_synth: IR bands are fitted from their own band lists.
Observability (depth, blending) is NOT established here.

The Codex fit evidence (orchestrator problem lines) is attached as evidence, not as a flag:
a fit failure is the empirical axis, not a statement about the gf.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "data/linelists/canonical_gf.csv"
KNOWN = ROOT / "data/reference/atomic_gap_known.csv"
PROBLEMS = ROOT / "data/results/orchestrator/{star}_problem_lines.csv"
DIAG = ROOT / "data/audit/orchestrator_prepare/{star}_{el}_band_diagnostics.csv"
OUT_DIR = ROOT / "data/audit/atomic_gap_register"

COMPILATIONS = {"MFW", "CSE", "1968PhFl.."}
MATCH_A = 0.02
SOLAR_DEX = 0.5
#: the ions the orchestrator measures; higher stages have no solar photospheric line
IONS = ("I", "II")
#: identity + the register's own columns. gf, its source, tier and uncertainty are NOT copied:
#: canonical_gf stays the one authority, joined by line_id (RYA-1378: no rival values).
COLS = ["line_id", "physical_id", "species", "wavelength_air_A", "excitation_potential_eV",
        "band", "in_synth", "held_reason", "no_graded_route_in_band", "flags", "priority",
        "evidence", "codex_fit_evidence", "review_state"]


def band_of(w: float) -> str:
    """The narrowest band_policy band holding w (H and K sit inside NIR)."""
    from pipeline import band_policy
    hits = [p for p in band_policy.POLICIES if p.lo_A <= w < p.hi_A]
    return min(hits, key=lambda p: p.hi_A - p.lo_A).name if hits else "outside"


def _solar_value(symbol: str) -> float | None:
    from pipeline.litscan import load_litscan
    rng = (load_litscan(symbol) or {}).get("range") or {}
    v = rng.get("central") if isinstance(rng, dict) else None
    return float(v) if isinstance(v, (int, float)) else None


def _excluded_by_registry(symbol: str, rows: pd.DataFrame) -> pd.Series:
    from pipeline import problem_children as pc
    t = pc.line_dispositions()
    return pd.Series([pc.aggregate_action(pc.disposition_for_line(
        symbol, str(sp).split()[1], float(w), table=t)) == "exclude"
        for sp, w in zip(rows.species, rows.wavelength_air_A)], index=rows.index)


def build(symbol: str, star: str = "solar") -> pd.DataFrame:
    from pipeline.gf_grades import is_gf_graded
    canon = pd.read_csv(CANON, low_memory=False)
    every = canon[canon.species.astype(str).isin([f"{symbol} {i}" for i in IONS])].copy()
    every["band"] = every.wavelength_air_A.astype(float).map(band_of)
    every = every[every.band != "outside"]
    lab = every.gf_tier.astype(str).str.contains("LAB", na=False)
    culled = _excluded_by_registry(symbol, every) if len(every) else pd.Series(dtype=bool)
    graded = is_gf_graded(every) if len(every) else pd.Series(dtype=bool)
    # a (species, band) with no graded line: every cell there is HELD (no measured route)
    route = set(zip(every.species[graded & ~culled], every.band[graded & ~culled]))

    c = every[~lab | culled].copy()
    c["held_reason"] = "GRADED_NOT_LAB"
    c.loc[~graded.reindex(c.index, fill_value=False), "held_reason"] = "UNGRADED_GF"
    c.loc[culled.reindex(c.index, fill_value=False), "held_reason"] = "REGISTRY_EXCLUDED"
    flags = {i: [] for i in c.index}
    evidence = {i: [] for i in c.index}
    for i in c.index:
        if not lab[i]:
            flags[i].append("UNRATED_GF")
    ref = c.loggf_reference.astype(str)
    for i in c.index[ref.isin(COMPILATIONS) | (ref == "VALD3")]:
        flags[i].append("PROVENANCE_GAP")
        evidence[i].append(f"gf via {ref[i]}: source paper not attached")

    def near(table, sp, w):
        m = (table.species == sp) & ((table.wavelength_air_A.astype(float) - w).abs() <= MATCH_A)
        return list(table.index[m])

    absent = []
    if KNOWN.exists():
        k = pd.read_csv(KNOWN)
        for r in k[k.element == symbol].itertuples():
            fl = [f for f in str(r.flags).split(";") if f]
            hits = near(c, r.species, float(r.wavelength_air_A))
            for i in hits:
                flags[i] += [f for f in fl if f not in flags[i]]
                evidence[i].append(f"{r.evidence} [{r.source}]")
            # absence is judged against EVERY row, LAB included: a curated line that has a
            # LAB row is in the reference pool, not missing
            if not near(every, r.species, float(r.wavelength_air_A)):
                w = float(r.wavelength_air_A)
                absent.append({"species": r.species, "wavelength_air_A": w, "band": band_of(w),
                               "in_synth": False, "held_reason": "NOT_IN_CANONICAL",
                               "flags": ";".join(["MISSING_GF"] + fl),
                               "evidence": f"{r.evidence} [{r.source}] | no row in canonical_gf"})

    diag = Path(str(DIAG).format(star=star, el=symbol))
    solar = _solar_value(symbol)
    if diag.exists():
        d = pd.read_csv(diag)
        for r in d.itertuples():
            a = pd.to_numeric(getattr(r, "A_1d_lte", None), errors="coerce")
            off = solar is not None and pd.notna(a) and abs(float(a) - solar) > SOLAR_DEX
            if bool(r.in_aggregate) and not off:
                continue
            for i in near(c, f"{symbol} {r.ion}", float(r.wavelength_air_A)):
                if "SOLAR_DISCREPANCY" not in flags[i]:
                    flags[i].append("SOLAR_DISCREPANCY")
                evidence[i].append(f"band diagnostic {r.holding}: A={a}"
                                   + ("" if bool(r.in_aggregate) else " (excluded)"))

    fit = {i: [] for i in c.index}
    prob = Path(str(PROBLEMS).format(star=star))
    if prob.exists():
        p = pd.read_csv(prob)
        p = p[p.element == symbol]
        for r in p.itertuples():
            for i in near(c, f"{symbol} {r.ion}", float(r.wavelength_air_A)):
                tag = f"{r.band}/{r.instrument}/{r.engine}:{r.reason_code}"
                if tag not in fit[i]:
                    fit[i].append(tag)

    c["flags"] = [";".join(flags[i]) for i in c.index]
    c["evidence"] = [" | ".join(evidence[i]) for i in c.index]
    c["codex_fit_evidence"] = [" | ".join(fit[i]) for i in c.index]
    if absent:
        c = pd.concat([c, pd.DataFrame(absent)], ignore_index=True)
    c["no_graded_route_in_band"] = [(s, b) not in route for s, b in zip(c.species, c.band)]
    disc = c["flags"].str.contains("DISCREPANCY")
    c["priority"] = 3
    c.loc[disc, "priority"] = 2
    c.loc[~c.in_synth.fillna(False).astype(bool) & ~disc, "priority"] = 4
    c.loc[c.no_graded_route_in_band, "priority"] = 1
    c["review_state"] = "candidate"
    for col in COLS:
        if col not in c:
            c[col] = None
    return c[COLS].sort_values(["priority", "species", "wavelength_air_A"]).reset_index(drop=True)


def summary(t: pd.DataFrame, symbol: str) -> dict:
    return {"element": symbol, "n_rows": int(len(t)),
            "held_reason": {k: int(v) for k, v in t.held_reason.value_counts().items()},
            "flags": {f: int(t["flags"].str.contains(f).sum()) for f in
                      ("UNRATED_GF", "MISSING_GF", "PROVENANCE_GAP", "GF_DISCREPANCY",
                       "SOLAR_DISCREPANCY")},
            "priority": {str(k): int(v) for k, v in t.priority.value_counts().sort_index().items()},
            "no_graded_route": sorted({f"{s} {b}" for s, b in zip(
                t.species[t.no_graded_route_in_band], t.band[t.no_graded_route_in_band])}),
            "with_codex_fit_evidence": int((t.codex_fit_evidence.fillna("") != "").sum())}


def write(symbol: str, star: str = "solar", out_dir: Path | None = None) -> dict:
    out = out_dir or OUT_DIR
    out.mkdir(parents=True, exist_ok=True)
    t = build(symbol, star)
    t.to_csv(out / f"{star}_{symbol}.csv", index=False)
    s = summary(t, symbol)
    s["candidate_for"] = "theoretical grade (RYA-1351) or laboratory measurement"
    s["gf_source"] = "data/linelists/canonical_gf.csv, joined by line_id"
    (out / f"{star}_{symbol}.json").write_text(json.dumps(s, indent=1) + "\n")
    return s


def render(s: dict) -> str:
    return (f"[gap-register] {s['element']}: {s['n_rows']} lines held out of the reference pool "
            f"{s['held_reason']}; flags {s['flags']}; priority {s['priority']}; "
            f"no graded route: {', '.join(s['no_graded_route']) or 'none'}")

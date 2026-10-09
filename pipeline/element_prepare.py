"""Governing-process steps 1-8 as ONE orchestrator stage, BEFORE any measurement (RYA-1233).

    python run_pipeline.py --star solar --element Al --prepare           # report only
    python run_pipeline.py --star solar --element Al --prepare --apply   # write the approved plan

Every check that RYA-1233 did by hand on Si -- and that, done late, cost a measurement rerun
each time -- run up front, into one report Ryan reviews once:

  A  published sets   every solar-literature line: in canonical_gf? the authors' gf vs the
                      canonical gf, the authors' sigma -> an ADOPTION plan (Ryan's rule: the
                      published set, with the gf the authors used)
  B  gf sources       the gf source of every graded line classified in gf_error_model.csv?
  C  cull candidates  reduced EW on a reference spectrum (RYA-458 ceiling -4.90) x membership
                      in a solar-literature set -> PROPOSED culls (saturated + unused)
  D  NLTE coverage    graded lines the element's ENGINE-A table does not serve
  E  holdings         each corrected holding serves the graded lines in its window, at rest
                      (line core within REST_TOL_KMS of the line)
  F  literature       litscan present; every declared bibliography key resolves

Written to data/results/orchestrator/prepare/<star>_<El>.json (+ .md). `--apply` writes ONLY
the adoption plan (authors' gf + sigma into canonical_gf, the old value kept in the audit) and
the proposed culls (problem_children, exclude/active, observed_in the star) and rebuilds the
element's line sets -- the one human checkpoint is reading the report before passing --apply.
"""
from __future__ import annotations

import csv
import io
import json
import math
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "results" / "orchestrator" / "prepare"
CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"
REGISTRY_PC = ROOT / "data" / "registry" / "problem_children.csv"
REW_CEILING = -4.90          # RYA-458 linear ceiling (pipeline.problem_children.REW_LINEAR_CEILING)
SAT_MIN_DEPTH = 0.30         # below this a line is on the linear part of the curve of growth
REST_TOL_KMS = 3.0
ADOPT_TIER = "PUBLISHED-SET"  # graded through its stated sigma (pipeline.gf_grades.is_gf_graded)
STAR_NAME = {"solar": "Sun"}


def _spec_rows(element: str, *, whole_selection: bool = False) -> list[tuple[dict, dict]]:
    """The declared sets' rows. `whole_selection` drops the build-time wavelength cut: a
    set's line SELECTION is literature membership even where its gf is not adopted
    (Deshmukh's optical lines are used, only their gf differs from Asplund's)."""
    from pipeline import line_sets as ls
    out = []
    for spec in ls.load_sources():
        if spec["element"] != element:
            continue
        s = {k: v for k, v in spec.items() if not (whole_selection and k == "min_wavelength_A")}
        for r in ls._source_rows(s):
            out.append((spec, r))
    return out


def _measure_rew(element: str, ion: str, w: float) -> dict:
    """Depth + reduced EW of the line itself on a reference spectrum (KP molecfit to 13000 A,
    CRIRES+ H beyond): the cull protocol's saturation test."""
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from measure_band_ew import load_window_ex
    hold = (("kpno_solar_atlas", "solar_kpno_molecfit_corrected") if w < 13000
            else ("crires_plus", "solar_crires_plus_h_rya1094"))
    try:
        x = load_window_ex(hold[0], w, 3.0, holding=hold[1])
    except Exception as exc:                                   # noqa: BLE001
        return {"holding": hold[1], "error": f"{type(exc).__name__}: {str(exc)[:80]}"}
    wl, f = np.asarray(x[0], float), np.asarray(x[1], float)
    g = np.isfinite(f)
    wl, f = wl[g], f[g]
    i = int(np.argmin(np.abs(wl - w)))
    depth = float(1 - f[max(0, i - 3):i + 4].min())
    hw = 0.35 * w / 6000 + 0.15
    s = np.abs(wl - w) <= hw
    ew = float(np.trapezoid(np.clip(1 - f[s], 0, None), wl[s])) * 1000
    rew = math.log10(ew / 1000 / w) if ew > 0 else float("nan")
    core = float(wl[max(0, i - 3) + int(np.argmin(f[max(0, i - 3):i + 4]))])
    # Saturated = past the RYA-458 REW ceiling AND deep enough to be on the flat part of the
    # curve of growth. A shallow line cannot be saturated: its window-integrated EW is
    # neighbours' absorption (Si 8435, depth 0.03, read rew -4.76 in a crowded window).
    return {"holding": hold[1], "depth": round(depth, 3), "ew_mA": round(ew, 1),
            "rew": round(rew, 3), "saturated": bool(rew > REW_CEILING and depth >= SAT_MIN_DEPTH),
            "core_dv_kms": round((core - w) / w * 299792.458, 2)}


def prepare(star: str, element: str) -> dict:
    from pipeline import gf_error_model, line_sets as ls, problem_children as pc
    from pipeline import run_matrix as rm
    from pipeline.gf_grades import is_gf_graded
    canon = pd.read_csv(CANONICAL, low_memory=False)
    sp_all = canon[canon.species.astype(str).str.split().str[0] == element]
    graded = sp_all[is_gf_graded(sp_all)]
    table = pc.line_dispositions()
    culled = {(s, round(w, 2)) for (s, w), _c in ls.culled_by_star().get(star, {}).items()}
    rep: dict = {"star": star, "element": element, "generated": str(date.today()),
                 "sets": [], "adopt": [], "missing_lines": [], "gf_sources": [],
                 "cull_candidates": [], "nlte_missing": [], "holdings": [], "literature": {}}

    # A -- published sets, the authors' gf, the adoption plan
    lit_waves = set()
    for spec, r in _spec_rows(element):
        c = spec["columns"]
        sp = ls._species(spec, r)
        w = float(r[c["wavelength"]])
        mt = spec.get("match") or {}
        if "canonical_id_col" in mt and str(r.get(mt["canonical_id_col"], "")).strip():
            # the set's OWN join (RYA-1169 / RYA-1058), not a wavelength window
            hit = canon[canon.line_id.astype(str) == str(r[mt["canonical_id_col"]]).strip()]
        else:
            tol = float(mt.get("canonical_tol_A", 0.012))
            hit = canon[(canon.species == sp) & ((canon.wavelength_air_A - w).abs() <= tol)]
        a_gf = float(r[c["loggf"]]) if str(r.get(c["loggf"], "")).strip() not in ("", "nan") else None
        a_sig = (float(r[c["loggf_sigma"]]) if "loggf_sigma" in c
                 and str(r.get(c["loggf_sigma"], "")).strip() not in ("", "nan") else None)
        row = {"set": spec["name"], "species": sp, "wavelength_A": w,
               "solar_literature": bool(spec.get("solar_literature")),
               "authors_loggf": a_gf, "authors_sigma_dex": a_sig}
        if spec.get("solar_literature"):
            lit_waves.add((sp, round(w, 1)))
        if len(hit) != 1:
            row["canonical"] = "MISSING" if len(hit) == 0 else f"AMBIGUOUS ({len(hit)})"
            rep["missing_lines"].append(row)
        else:
            h = hit.iloc[0]
            row.update(canonical_wavelength_A=float(h.wavelength_air_A),
                       canonical_loggf=float(h.log_gf), canonical_tier=str(h.gf_tier),
                       canonical_sigma=(None if pd.isna(h.gf_sigma_dex) else float(h.gf_sigma_dex)),
                       d_loggf=(None if a_gf is None else round(a_gf - float(h.log_gf), 3)))
            if (spec.get("solar_literature") and a_gf is not None and a_sig is not None
                    and (abs(a_gf - float(h.log_gf)) > 0.0005 or row["canonical_sigma"] != a_sig)):
                rep["adopt"].append({**row, "line_id": str(h.line_id)})
        rep["sets"].append(row)

    for spec, r in _spec_rows(element, whole_selection=True):
        if spec.get("solar_literature"):
            lit_waves.add((ls._species(spec, r), round(float(r[spec["columns"]["wavelength"]]), 1)))

    # B -- gf sources of graded lines, classified?
    model = gf_error_model.load()
    tags = graded.lab_source_tag.fillna("").astype(str)
    for t, n in tags.value_counts().items():
        rep["gf_sources"].append({"lab_source_tag": t or "(none: NIST grade / stored sigma)",
                                  "n_graded_lines": int(n), "classified": t in model})
    for sname in sorted({a["set"] for a in rep["adopt"]}):
        rep["gf_sources"].append({"lab_source_tag": f"{sname} (to be adopted)",
                                  "n_graded_lines": sum(a["set"] == sname for a in rep["adopt"]),
                                  "classified": sname in model})

    # C -- cull candidates: graded, not yet culled, saturated on the reference spectrum, and
    #      not a line the solar literature measured the Sun on
    for _, g in graded.iterrows():
        sp, w = str(g.species), float(g.wavelength_air_A)
        if (sp, round(w, 2)) in culled:
            continue
        m = _measure_rew(element, sp.split()[1], w)
        if m.get("saturated") and (sp, round(w, 1)) not in lit_waves:
            rep["cull_candidates"].append({"species": sp, "wavelength_A": w,
                                           "problem_class": "SATURATION_COG", **m})

    # D -- NLTE coverage of the graded, un-culled lines
    nlte = ROOT / "data" / "nlte_grids" / f"{element}_Amarsi2020_PySME.csv"
    served = set(pd.read_csv(nlte).wave_A.round(1)) if nlte.exists() else set()
    for _, g in graded.iterrows():
        sp, w = str(g.species), float(g.wavelength_air_A)
        if sp.endswith(" I") and (sp, round(w, 2)) not in culled and round(w, 1) not in served:
            rep["nlte_missing"].append({"species": sp, "wavelength_A": w})
    rep["nlte_table"] = str(nlte.relative_to(ROOT)) if nlte.exists() else None

    # E -- every corrected holding serves its window's graded lines, at rest
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from measure_band_ew import load_window_ex
    seen = set()
    for d in rm.expand(star, element, engines=["ts-lte"], methods=["synthesis"], pools=["codex"]):
        key = (d.holding, d.instrument, d.lo_A, d.hi_A)
        if key in seen:
            continue
        seen.add(key)
        inwin = graded[(graded.wavelength_air_A >= d.lo_A) & (graded.wavelength_air_A <= d.hi_A)]
        res = {"holding": d.holding, "band": d.band, "lo_A": d.lo_A, "hi_A": d.hi_A,
               "graded_lines": len(inwin), "served": 0, "tested": 0, "off_rest": [],
               "not_served": []}
        for _, g in inwin.head(5).iterrows():
            w = float(g.wavelength_air_A)
            res["tested"] += 1
            try:
                x = load_window_ex(d.instrument, w, 0.6, holding=d.holding)
            except Exception as exc:                           # noqa: BLE001
                res["not_served"].append({"wavelength_A": w,
                                          "why": f"{type(exc).__name__}: {str(exc)[:90]}"})
                continue
            wl, f = np.asarray(x[0], float), np.asarray(x[1], float)
            ok = np.isfinite(f)
            if ok.sum() < 10:
                res["not_served"].append({"wavelength_A": w, "why": "fewer than 10 finite pixels"})
                continue
            res["served"] += 1
            # The line's OWN core: within 3x the rest tolerance of the line, never the deepest
            # pixel of the whole window (that is a neighbour's -- 6696.185 sits beside the much
            # stronger 6696.023). A line too weak to locate is not tested for frame.
            near = ok & (np.abs(wl - w) <= w * REST_TOL_KMS / 299792.458 + 0.02)
            if near.sum() < 3 or (1 - np.nanmin(f[near])) < 0.03:
                continue
            core = float(wl[near][int(np.nanargmin(f[near]))])
            dv = (core - w) / w * 299792.458
            if abs(dv) > REST_TOL_KMS:
                res["off_rest"].append({"wavelength_A": w, "dv_kms": round(dv, 1)})
        # A holding is broken when it serves NOTHING it should, or serves lines off rest. A
        # single line in a coverage gap (KP 1984 has none at 11253 A) is a coverage note.
        res["ok"] = res["graded_lines"] == 0 or (res["served"] > 0 and not res["off_rest"])
        rep["holdings"].append(res)

    # F -- literature
    lit = ROOT / "data" / "reference" / "litscan" / f"{element}.yaml"
    bib = {r["key"] for r in csv.DictReader((ROOT / "data/refs/bibliography.csv").open())}
    keys = {spec["bib"] for spec in ls.load_sources() if spec["element"] == element}
    rep["literature"] = {"litscan": lit.exists(), "set_bib_keys": sorted(keys),
                         "missing_bib_keys": sorted(keys - bib)}

    rep["ready"] = (rep["literature"]["litscan"] and not rep["literature"]["missing_bib_keys"]
                    and not rep["missing_lines"] and not rep["adopt"] and not rep["cull_candidates"]
                    and all(h["ok"] for h in rep["holdings"]))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{star}_{element}.json").write_text(json.dumps(rep, indent=1, default=str) + "\n")
    (OUT / f"{star}_{element}.md").write_text(render(rep) + "\n")
    return rep


def render(rep: dict) -> str:
    L = [f"PREPARE {rep['star']} {rep['element']}  -> "
         + ("READY to measure" if rep["ready"] else "NOT READY -- review, then --apply")]
    lit = rep["literature"]
    L.append(f"  F literature: litscan {'yes' if lit['litscan'] else 'MISSING'}; "
             f"set bibliography keys {lit['set_bib_keys']}"
             + (f"; MISSING {lit['missing_bib_keys']}" if lit["missing_bib_keys"] else ""))
    L.append(f"  A published sets: {len(rep['sets'])} lines; adopt authors' gf on "
             f"{len(rep['adopt'])}; missing from canonical_gf: {len(rep['missing_lines'])}")
    for a in rep["adopt"]:
        L.append(f"      adopt {a['set']} {a['species']} {a['wavelength_A']}: canonical "
                 f"{a['canonical_loggf']} -> authors {a['authors_loggf']} +/- {a['authors_sigma_dex']}"
                 f"  (d {a['d_loggf']:+.3f})")
    for m in rep["missing_lines"]:
        L.append(f"      MISSING {m['set']} {m['species']} {m['wavelength_A']}: {m['canonical']} "
                 f"-- add it to canonical_gf (VALD/NIST extraction) before it can be measured")
    for s in rep["gf_sources"]:
        L.append(f"  B gf source {s['lab_source_tag']}: {s['n_graded_lines']} graded lines, "
                 f"{'classified' if s['classified'] else 'NOT in gf_error_model (UNREVIEWED)'}")
    L.append(f"  C cull candidates (saturated, not a solar-literature line): {len(rep['cull_candidates'])}")
    for c in rep["cull_candidates"]:
        L.append(f"      {c['species']} {c['wavelength_A']}: rew {c.get('rew')} depth {c.get('depth')}")
    L.append(f"  D NLTE: {len(rep['nlte_missing'])} graded line(s) not in {rep['nlte_table']}")
    bad = [h for h in rep["holdings"] if not h["ok"]]
    L.append(f"  E holdings: {len(rep['holdings'])} windows; {len(bad)} not serving / off rest")
    for h in bad:
        L.append(f"      {h['holding']} {h['band']}: served {h['served']}/{h['tested']}"
                 + (f"; off rest {h['off_rest']}" if h["off_rest"] else "")
                 + (f"; NOT served {h['not_served']}" if h["not_served"] else ""))
    return "\n".join(L)


def apply(star: str, element: str, rep: dict) -> dict:
    """Write the reviewed plan: adoptions into canonical_gf, culls into problem_children."""
    from pipeline import line_sets as ls
    done = {"adopted": 0, "culled": 0}
    if rep["adopt"]:
        cg = pd.read_csv(CANONICAL, low_memory=False)
        audit = []
        for a in rep["adopt"]:
            i = cg.index[cg.line_id.astype(str) == a["line_id"]]
            if len(i) != 1:
                continue
            i = i[0]
            audit.append({"line_id": a["line_id"], "species": a["species"],
                          "wavelength_A": a["canonical_wavelength_A"], "old_loggf": cg.at[i, "log_gf"],
                          "old_tier": cg.at[i, "gf_tier"], "new_loggf": a["authors_loggf"],
                          "new_sigma_dex": a["authors_sigma_dex"], "set": a["set"]})
            cg.at[i, "log_gf"] = a["authors_loggf"]
            cg.at[i, "gf_sigma_dex"] = a["authors_sigma_dex"]
            cg.at[i, "gf_tier"] = ADOPT_TIER
            cg.at[i, "lab_source_tag"] = a["set"]
            cg.at[i, "adjudication_status"] = f"adopted_rya1233_{a['set'].lower()}"
            done["adopted"] += 1
        cg.to_csv(CANONICAL, index=False)
        ad = ROOT / "data" / "audit" / "orchestrator_prepare" / f"{star}_{element}_adoption.csv"
        ad.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(audit).to_csv(ad, index=False)
    if rep["cull_candidates"]:
        b = REGISTRY_PC.read_bytes()
        s = io.StringIO()
        w = csv.writer(s, lineterminator="\r\n")
        for c in rep["cull_candidates"]:
            w.writerow([c["species"], f"{c['wavelength_A']:.3f}", c["problem_class"], "exclude",
                        STAR_NAME.get(star, star), "Teff↓ logg↑", "high", "1233,458", "active",
                        "curated",
                        f"Saturated: reduced EW log(W/lambda) = {c.get('rew')} > the RYA-458 ceiling "
                        f"{REW_CEILING} ({c.get('holding')}, depth {c.get('depth')}); not a line any "
                        f"declared solar-literature set measured the Sun on. Culled for this star by "
                        f"orchestrator --prepare --apply (reviewed); kept in every line list."])
            done["culled"] += 1
        REGISTRY_PC.write_bytes(b + s.getvalue().encode())
    ls.build(element, quiet=True)
    return done

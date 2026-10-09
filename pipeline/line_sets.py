"""The published reference line sets, built from their DECLARATIONS (RYA-1233).

    python -m pipeline.line_sets [--element Si]

Every set in data/reference/line_sets/SOURCES.yaml becomes, per species:
  <set>_<species>.csv              every published line, kept for good (with gf_graded and
                                   culled_in = star:class)
  <set>_<species>_graded.csv       the lines whose gf carries a published uncertainty
                                   (pipeline.gf_grades.is_gf_graded) -- for any star
  <set>_<species>_graded_<star>.csv  that star's problem_children culls removed (written only
                                   where the star culls something)
and one row in REGISTRY.csv, which the orchestrator reads (run_matrix.line_sets_for,
run_descriptor.star_graded_csv). Replaces scripts/rya1233_build_si_line_sets.py, whose sources
were typed into the script -- so Al's published set (RYA-1173) could not be registered
without writing another script.

`match_tol_A` is the window a set line is matched into the synthesis list with: half a unit of
the source's last printed decimal, floored at RYA-1171's 0.005 A dual-key tolerance.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "reference" / "line_sets"
SOURCES = OUT / "SOURCES.yaml"
REGISTRY = OUT / "REGISTRY.csv"
CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"
FLOOR_A = 0.005
EV_PER_CM = 1 / 8065.544
#: problem_children `observed_in` (the literature's name) -> the orchestrator's star id.
STAR_ID = {"Sun": "solar"}


def load_sources() -> list[dict]:
    return yaml.safe_load(SOURCES.read_text())["sets"]


def _tol(text: str) -> float:
    dec = len(text.split(".")[1]) if "." in text else 0
    return max(FLOOR_A, 0.5 * 10 ** (-dec))


def _source_rows(spec: dict) -> list[dict]:
    with (ROOT / spec["source"]).open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    for k, v in (spec.get("where") or {}).items():
        rows = [r for r in rows if str(r.get(k, "")).strip() == str(v)]
    lo = spec.get("min_wavelength_A")
    if lo is not None:
        rows = [r for r in rows if float(r[spec["columns"]["wavelength"]]) >= float(lo)]
    return rows


def _species(spec: dict, r: dict) -> str:
    s = spec["species"]
    if isinstance(s, str):
        return s
    v = str(r[s["column"]]).strip()
    return (s.get("map") or {}).get(v, v)


def rows_for(spec: dict, canon: pd.DataFrame) -> list[tuple]:
    """(set, species, wavelength to MATCH on, ep, published log gf, band, published wavelength).

    The match wavelength is the CANONICAL one wherever the transcription already joined the
    line to its canonical id, or the unique canonical row within `canonical_tol_A`: a source
    printed to 0.01 A sits up to 0.012 A from the list (Amarsi's 6741.64 vs 6741.628), which
    a printed-precision window misses. The published wavelength is kept beside it."""
    c, m = spec["columns"], spec.get("match") or {}
    by_id = dict(zip(canon.line_id.astype(str), canon.wavelength_air_A.astype(str)))
    out = []
    for r in _source_rows(spec):
        sp = _species(spec, r)
        pub = str(r[c["wavelength"]]).strip()
        if "canonical_id_col" in m:
            match = by_id.get(str(r.get(m["canonical_id_col"], "")).strip(), pub)
        elif "canonical_tol_A" in m:
            w = float(pub)
            hit = canon[(canon.species == sp) & ((canon.wavelength_air_A - w).abs() <= float(m["canonical_tol_A"]))]
            if len(hit) > 1:
                raise ValueError(f"{spec['name']} {sp} {pub}: {len(hit)} canonical rows within "
                                 f"{m['canonical_tol_A']} A -- refusing a guessed match")
            # NONE: the published line is not in canonical_gf yet (a step-6 gap the prepare
            # report names). Kept at its published wavelength; it cannot be graded until added.
            match = str(hit.wavelength_air_A.iloc[0]) if len(hit) else pub
        else:
            match = pub
        ep = (str(r[c["ep"]]).strip() if "ep" in c
              else f"{float(r[c['elow_cm']]) * EV_PER_CM:.3f}")
        band = c.get("band_value") or str(r[c["band"]]).strip()
        out.append((spec["name"], sp, match, ep, str(r[c["loggf"]]).strip(), band, pub))
    return out


def culled_by_star() -> dict:
    """star -> {(species, wavelength): problem_class}: the registry's `exclude` + `active` rows
    (RYA-807's discriminator). A cull is per star; the set itself keeps every line."""
    out: dict = {}
    with (ROOT / "data/registry/problem_children.csv").open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["required_treatment"].strip() != "exclude" or r["status"].strip() != "active":
                continue
            mm = re.match(r"\s*([0-9]+\.[0-9]+)", r["lambda_or_scope"])
            if not mm:
                continue
            for obs in r["observed_in"].split(";"):
                star = STAR_ID.get(obs.strip())
                if star:
                    out.setdefault(star, {})[(r["species"].strip(), float(mm.group(1)))] = \
                        r["problem_class"].strip()
    return out


def _cull_of(culled: dict, species: str, wave: float) -> str:
    return next((cl for (sp, w), cl in culled.items() if sp == species and abs(w - wave) < 0.01), "")


def _write(path: Path, head: list, rows: list) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(head)
        w.writerows(rows)


def build(element: str | None = None, *, quiet: bool = False) -> list[dict]:
    from pipeline.gf_grades import is_gf_graded
    canon = pd.read_csv(CANONICAL, low_memory=False)
    graded_mask = is_gf_graded(canon)
    culled = culled_by_star()
    specs = load_sources()
    groups: dict = {}
    meta = {s["name"]: s for s in specs}
    for spec in specs:
        if element and spec["element"] != element:
            continue
        for row in rows_for(spec, canon):
            groups.setdefault((row[0], row[1]), []).append(row)

    def graded(species: str, wave: float) -> bool:
        m = (canon.species == species) & ((canon.wavelength_air_A - wave).abs() < 0.01)
        return bool(m.sum() == 1 and graded_mask[m].iloc[0])

    reg = []
    for (name, species), rows in sorted(groups.items()):
        spec = meta[name]
        tol = min(_tol(r[6]) for r in rows)
        stem = f"{name.lower()}_{species.replace(' ', '')}"
        fn, gfn = OUT / f"{stem}.csv", OUT / f"{stem}_graded.csv"
        head = ["line_set", "species", "wavelength_air_A", "wavelength_published_A", "ep_eV",
                "loggf_published", "band_published", "match_tol_A", "source", "doi",
                "gf_graded", "culled_in"]
        table = []
        for r in sorted(rows, key=lambda r: float(r[2])):
            g = graded(species, float(r[2]))
            cin = ";".join(f"{st}:{cl}" for st in sorted(culled)
                           if (cl := _cull_of(culled[st], species, float(r[2]))))
            table.append([r[0], r[1], r[2], r[6], r[3], r[4], r[5], tol, spec["citation"],
                          spec["doi"], g, cin])
        _write(fn, head, table)
        gr = [t for t in table if t[10]]
        _write(gfn, head, gr)
        per_star = {}
        for st in sorted(culled):
            keep = [t for t in gr if not _cull_of(culled[st], species, float(t[2]))]
            sf = OUT / f"{stem}_graded_{st}.csv"
            if len(keep) != len(gr):
                _write(sf, head, keep)
                per_star[st] = len(keep)
            elif sf.exists():
                sf.unlink()
        el, ion = species.split()
        ws = [float(r[2]) for r in rows]
        reg.append({"element": el, "ion": ion, "set_name": name, "csv": str(fn.relative_to(ROOT)),
                    "n_lines": len(rows), "graded_csv": str(gfn.relative_to(ROOT)),
                    "n_graded": len(gr), "n_graded_solar": per_star.get("solar", len(gr)),
                    "lo_A": min(ws), "hi_A": max(ws), "bibliography_key": spec["bib"],
                    "doi": spec["doi"], "source": spec["citation"], "gf_basis": spec["gf_basis"],
                    "solar_literature": bool(spec.get("solar_literature")), "ticket": "RYA-1233"})
        if not quiet:
            print(f"  {name:<17} {species:<6} {len(rows):3d} lines ({len(gr)} graded, "
                  f"{per_star.get('solar', len(gr))} for solar)  {min(ws):8.1f}-{max(ws):8.1f} A")
    # The registry carries every element: rebuilding one element keeps the others' rows.
    if element and REGISTRY.exists():
        with REGISTRY.open(newline="") as fh:
            reg = [r for r in csv.DictReader(fh) if r["element"] != element] + reg
        reg.sort(key=lambda r: (r["element"], r["set_name"], r["ion"]))
    with REGISTRY.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(reg[0]))
        w.writeheader()
        w.writerows(reg)
    return reg


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--element", default=None)
    a = ap.parse_args(argv)
    reg = build(a.element)
    print(f"registry: {REGISTRY.relative_to(ROOT)} ({len(reg)} sets)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

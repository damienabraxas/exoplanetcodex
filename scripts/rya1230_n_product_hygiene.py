"""RYA-1230 Part A -- what does each N.json product actually measure, and which are duplicates.

A1. SPECIES COMPOSITION. For every product, count the lines inside its OWN fit windows by
    species: atomic rows from the band's synthesis line list, molecular rows from the
    vendored Turbospectrum .bsyn lists (RYA-360; the files the synthesis reads). The
    windows are read through `rya1214_cnosynth_to_products._windows_for`, the same call
    that sized the published `n_lines`, so the count is reproduced, not re-derived.

A2. DUPLICATES. Group products on the feed's own identity key (`plot_grid.key_fields`), with
    the selector NORMALISED by `pipeline.band_products.TREATMENTS` -- the normalisation
    RYA-1214's xi-cost script already applies, because an earlier publish split the artifact
    stem on `_SYNTH_` and kept `SET-AGSS21_ENGINE-A` while `treatment` already said
    ENGINE-A. A group of >1 is a duplicate when the rows are value-identical AND point at
    the same artifact bytes.

Read-only by default. `--apply` retires the redundant member of each duplicate group into
`archive` -- the publisher's own place for SUPERSEDED rows (RYA-711: never deleted) -- and
NOT into `quarantine`, which means WITHDRAWN: a duplicate is not a rejected measurement,
and its correctly keyed twin, carrying the identical value, stays live. The write goes
through `publish_product.write_feed`, so the RYA-587 publication gate and the plot-grid
rebuild run exactly as on any publish. No value, gf or line list changes.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

MOLDIR = REPO / "data/linelists/molecular/turbospectrum"
BAND_PRODUCTS = REPO / "data/results/band_products"
VALUE_FIELDS = ("A", "sigma_stat", "sigma_syst", "n_lines", "n_excluded")


def _bsyn_species_in(windows) -> Counter:
    """Every molecular line inside `windows`, keyed by isotopologue (from the filename)."""
    out = Counter()
    for f in sorted(MOLDIR.rglob("*.bsyn")):
        m = re.match(r"(.*?)_(\d+)-(\d+)\.bsyn$", f.name)
        if m:
            lo_nm, hi_nm = float(m.group(2)), float(m.group(3))
            if not any(lo_nm <= lo / 10 <= hi_nm or lo_nm <= hi / 10 <= hi_nm
                       or (lo / 10 <= lo_nm and hi / 10 >= hi_nm) for lo, hi in windows):
                continue
            species = m.group(1)
        else:
            species = f.stem
        with f.open() as fh:
            for ln in fh:
                p = ln.split()
                if not p:
                    continue
                try:
                    w = float(p[0])
                except ValueError:
                    continue
                if any(lo <= w <= hi for lo, hi in windows):
                    out[species] += 1
    return out


def _atomic_species_in(linelist_tsv: Path, windows) -> Counter:
    out = Counter()
    with linelist_tsv.open() as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            w = float(r["wave_A"]) if "wave_A" in r else float(r["wave_nm"]) * 10.0
            if any(lo <= w <= hi for lo, hi in windows):
                out[r["element"]] += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true",
                    help="retire each duplicate group's mis-keyed member into archive")
    ap.add_argument("--out", type=Path,
                    default=REPO / "data/output/rya1230/n_product_hygiene.json")
    args = ap.parse_args()

    from pipeline.band_products import TREATMENTS
    from pipeline.abundances_derive import ISPEC_DIR
    import rya1214_cnosynth_to_products as shim

    doc = json.loads((REPO / "data/products/solar/N.json").read_text())
    key_fields = doc["plot_grid"]["key_fields"]
    ges = (ISPEC_DIR / "input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm"
           / "atomic_lines.tsv")
    ir = REPO / "data/linelists/ispec_ir_9200_13000/atomic_lines.tsv"

    def norm_selector(sel: str) -> str:
        for t in TREATMENTS:
            if sel.endswith("_" + t):
                return sel[: -(len(t) + 1)]
        return sel

    products = []
    for i, p in enumerate(doc["products"]):
        sel = p.get("selector") or ""
        rec = dict(index=i, selector=sel, treatment=p["treatment"], band=p["band"],
                   holding=p["holding"], ion=p["ion"], A=p["A"], n_lines=p["n_lines"],
                   artifact=p["provenance"].get("copied_to"),
                   sha256=p["provenance"].get("sha256"))
        if sel.startswith("MOL-"):
            diag = sel[4:]
            region = {("CN_red", "solar_harps_molecfit_corrected"): "vis",
                      ("CN_AX_IR", "solar_iag"): "nir_cn_iag",
                      ("CN_AX_IR", "solar_kpno_molecfit_corrected"): "nir_cn_kp"}[
                (diag, p["holding"])]
            windows = shim._windows_for(region, diag)
            mol = _bsyn_species_in(windows)
            atom = _atomic_species_in(ges if p["band"] == "VIS" else ir, windows)
            molecule = shim.DIAG[diag][2]
            recount = shim._count_molecular_lines(molecule, windows)
            rec.update(kind="MOLECULAR", diagnostic=diag, region=region,
                       windows_A=[list(w) for w in windows],
                       window_span_A=round(sum(hi - lo for lo, hi in windows), 3),
                       molecular_lines_by_species=dict(mol.most_common()),
                       atomic_lines_by_species_top=dict(atom.most_common(8)),
                       n_I_lines_in_windows=int(atom.get("N 1", 0)),
                       n_lines_is=f"{molecule} lines in the fit windows",
                       n_lines_recount=recount,
                       n_lines_reproduced=(recount == p["n_lines"]),
                       measures=f"A(N) from the {molecule} (CN A-X) band at fixed A(C); "
                                f"NOT an N I measurement")
        else:
            lines_csv = BAND_PRODUCTS / Path(rec["artifact"]).name.replace(
                "_products.csv",
                "_lines.csv" if p["treatment"] != "1D-LTE" else "_1D-LTE_lines.csv")
            if not lines_csv.exists():
                lines_csv = BAND_PRODUCTS / Path(rec["artifact"]).name.replace(
                    "_products.csv", "_1D-LTE_lines.csv")
            rows = list(csv.DictReader(lines_csv.open()))
            rec.update(kind="ATOMIC", lines_csv=str(lines_csv.relative_to(REPO)),
                       species=sorted({f"{r['element']} {r['ion']}" for r in rows}),
                       lines=[float(r["wavelength_air_A"]) for r in rows],
                       in_aggregate=[float(r["wavelength_air_A"]) for r in rows
                                     if r["in_aggregate"] == "True"],
                       measures="A(N) from N I atomic lines (AGSS21 five-line set)")
        rec["key"] = {k: (norm_selector(sel) if k == "selector" else p.get(k))
                      for k in key_fields}
        products.append(rec)

    groups = defaultdict(list)
    for rec in products:
        groups[json.dumps(rec["key"], sort_keys=True)].append(rec)
    dups = []
    for k, g in groups.items():
        if len(g) < 2:
            continue
        vals = {json.dumps({f: doc["products"][r["index"]].get(f) for f in VALUE_FIELDS})
                for r in g}
        dups.append(dict(key=json.loads(k), indices=[r["index"] for r in g],
                         selectors=[r["selector"] for r in g],
                         value_identical=len(vals) == 1,
                         same_artifact_bytes=len({r["sha256"] for r in g}) == 1,
                         keep=min((r for r in g), key=lambda r: r["selector"] != r["key"]["selector"])["index"]))
    report = dict(ticket="RYA-1230", part="A hygiene", feed_version=doc["version"],
                  n_products=len(products), key_fields=key_fields, products=products,
                  duplicate_groups=dups, n_duplicate_rows=sum(len(d["indices"]) - 1 for d in dups))
    if args.apply and dups:
        import publish_product as pp
        drop = set()
        for d in dups:
            if not (d["value_identical"] and d["same_artifact_bytes"]):
                raise SystemExit(f"refusing: group {d['indices']} is not a pure duplicate")
            drop |= set(d["indices"]) - {d["keep"]}
        now = pp._now()
        for i in sorted(drop):
            old = dict(doc["products"][i])
            twin = next(d["keep"] for d in dups if i in d["indices"])
            old["superseded_at"] = now
            old["superseded_reason"] = (
                "RYA-1230 hygiene: DUPLICATE, not a withdrawn measurement. Same artifact "
                "sha256 and same A/sigma/n as the live row keyed selector="
                f"{doc['products'][twin]['selector']} treatment={doc['products'][twin]['treatment']}. "
                f"This row's selector {old['selector']} is the mis-parse RYA-1214's xi-cost "
                "script documents (an earlier publish split the artifact stem on _SYNTH_ and "
                "kept the treatment suffix). The correctly keyed twin stays live.")
            doc.setdefault("archive", []).append(old)
        doc["products"] = [p for i, p in enumerate(doc["products"]) if i not in drop]
        doc["version"] = pp.bump(doc["version"])
        doc["updated_at"] = now
        pp.write_feed(REPO / "data/products/solar/N.json", doc)
        report["applied"] = dict(archived_indices=sorted(drop), new_version=doc["version"],
                                 n_live=len(doc["products"]))
        print(f"APPLIED: {len(drop)} duplicate row(s) -> archive; N.json v{doc['version']}, "
              f"{len(doc['products'])} live")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=1, default=str) + "\n")
    for r in products:
        tag = (f"{r['diagnostic']}: n_lines={r['n_lines']} recount={r['n_lines_recount']} "
               f"N_I_in_windows={r['n_I_lines_in_windows']} mol={r['molecular_lines_by_species']}"
               if r["kind"] == "MOLECULAR" else f"{r['species']} {r['lines']}")
        print(f"[{r['index']:2d}] {r['selector']:20s} {r['treatment']:9s} {r['holding']:34s} "
              f"A={r['A']}  {tag}")
    for d in dups:
        print("DUP", d["indices"], d["selectors"], "value_identical", d["value_identical"],
              "same_bytes", d["same_artifact_bytes"], "keep", d["keep"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

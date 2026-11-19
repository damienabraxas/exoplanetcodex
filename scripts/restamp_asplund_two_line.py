#!/usr/bin/env python3
"""Restamp every live product the budget priced from the scatter of TWO lines -- RYA-470.

Ryan 2026-11-19: "We do all uncertainty the Asplund way." Asplund+2021 Sect. 2.1: "In cases
when only one or two lines can be used, we estimated the statistical error from the goodness of
the profile fitting". The budget assembler priced two resolved features by their scatter
(SD/sqrt 2; 0.000 when the two agree) until 4eb34bde. This replaces that measurement term with
the fit-priced one (`rya1230_cno_budget.fit_stat`) and rebuilds the document through the
contract's own `assemble`; `publish_product.write_feed` restamps sigma_asplund and gates it.

NOTHING IS RE-MEASURED. Each product's own per-line file supplies the two lines' chi2-curvature
sigmas, and a file is accepted only if its accepted features REPRODUCE the stored scatter
(raw_sigma) and the product's A -- proof it is the pool the product was measured on. A product
with no such file, or a line with no recorded sigma, is reported and left untouched.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]

import publish_product as pp  # noqa: E402
import rya1230_cno_budget as bud  # noqa: E402
from pipeline.uncertainty_contract import assemble, fit_measurement  # noqa: E402

REASON = ("RYA-470 (Ryan 2026-11-19, 'We do all uncertainty the Asplund way'): a two-line "
          "product's statistical term is priced from the goodness of the two profile fits "
          "(Asplund+2021 Sect. 2.1), not from the scatter of two numbers; no line re-measured")


def is_two_line_scatter(p: dict) -> bool:
    m = next((c for c in (p.get("uncertainty") or {}).get("components", [])
              if c["name"] == "measurement"), None)
    ev = (m or {}).get("evidence") or {}
    return ev.get("method") == "line_scatter" and ev.get("n_lines") == 2


def _rel(f: Path) -> str:
    try:
        return str(f.relative_to(ROOT))
    except ValueError:
        return str(f)


def _candidates(p: dict, tracked: list[str], extra: list[str] = ()) -> list[Path]:
    ct = (p.get("provenance") or {}).get("copied_to") or ""
    t = p.get("treatment")
    out = [ct.replace("_products.csv", f"_{t}_lines.csv"), ct.replace("_products.csv", "_lines.csv")]
    stem = Path(ct).name[: -len("_products.csv")] if ct.endswith("_products.csv") else ""
    stem = re.sub(r"_(ENGINE-A-3DNLTE|ENGINE-A|ENGINE-B[^_]*|1D-LTE)$", "", stem)
    if stem:
        out += [f for f in list(tracked) + list(extra) if Path(f).name == f"{stem}_{t}_lines.csv"]
    return [ROOT / f for f in dict.fromkeys(out) if f and (ROOT / f).exists()]


def _sigmas(lines: Path, acc: pd.DataFrame) -> dict:
    """Per-line chi2-curvature sigma, falling back to the 1D-LTE fit beside an engine file --
    the same rule the assembler uses (a departure leg adds a delta to that SAME fit)."""
    sa = (pd.to_numeric(acc["sigma_A"], errors="coerce") if "sigma_A" in acc
          else pd.Series(np.nan, index=acc.index))
    lte = lines.parent / re.sub(r"_(ENGINE-A-3DNLTE|ENGINE-A|ENGINE-B[^_]*)_lines\.csv$",
                                "_1D-LTE_lines.csv", lines.name)
    l1 = pd.read_csv(lte) if lte.exists() and lte != lines else None
    out = {}
    for w, s in zip(acc["wavelength_air_A"].astype(float), sa):
        if not np.isfinite(s) and l1 is not None and "sigma_A" in l1:
            m = l1[(l1.wavelength_air_A - w).abs() < 1e-3]
            s = float(pd.to_numeric(m["sigma_A"], errors="coerce").iloc[0]) if len(m) else np.nan
        out[round(w, 3)] = float(s)
    return out


def new_measurement(p: dict, tracked: list[str], extra: list[str] = ()) -> tuple[dict | None, str]:
    m = next(c for c in p["uncertainty"]["components"] if c["name"] == "measurement")
    raw = float(m["evidence"]["raw_sigma"])
    for f in _candidates(p, tracked, extra):
        df = pd.read_csv(f)
        if "in_aggregate" not in df:
            continue
        acc = df[df["in_aggregate"] == True]  # noqa: E712
        feats = bud._features(acc)
        if len(feats) != 2:
            continue
        fv = [float(np.mean([v for _, v in g])) for g in feats]
        if not (math.isclose(float(np.std(fv, ddof=1)), raw, rel_tol=1e-6, abs_tol=1e-9)
                and abs(float(np.median(fv)) - float(p["A"])) <= 6e-4):
            continue                                   # not this product's pool
        by_w = _sigmas(f, acc)
        fs = [float(np.mean([by_w[round(w, 3)] for w, _ in g])) for g in feats]
        if not all(np.isfinite(fs)):
            return None, f"{_rel(f)}: a line records no curvature sigma"
        npx = int(pd.to_numeric(acc.get("n_pix"), errors="coerce").fillna(0).sum()) \
            if "n_pix" in acc else 0
        comp = fit_measurement(
            bud.fit_stat(fs),
            source=("Asplund+2021 Sect. 2.1: one or two lines -> statistical error from the "
                    "goodness of the profile fit; each line's chi2-curvature sigma, rescaled to "
                    "red_chi2 = 1; two lines -> sqrt(s1^2 + s2^2) / 2, the sigma of their mean "
                    "(= the median the product reports). Restamped from the product's own "
                    f"per-line file {_rel(f)} (RYA-470)"),
            likelihood="synthesis chi2 over the fit window, curvature at the best abundance",
            correlation_treatment=("red_chi2 rescaling absorbs pixel correlation and model "
                                   "inadequacy; the worse-constrained side is taken"),
            n_pixels=max(npx, 1))
        comp["evidence"].update(n_lines=2, line_sigma_dex=fs, replaced_scatter_dex=m["sigma_dex"],
                                lines_file=str(_rel(f)))
        return comp, str(_rel(f))
    return None, "no per-line file reproduces this product's two-line pool"


def rebuild(p: dict, comp: dict) -> dict:
    u = p["uncertainty"]
    comps = [comp if c["name"] == "measurement" else c for c in u["components"]]
    cov = [list(r) for r in u["covariance"]]
    i = u["covariance_order"].index("measurement")
    if any(cov[i][j] or cov[j][i] for j in range(len(cov)) if j != i):
        raise ValueError("measurement is correlated with another term; refusing to rescale it")
    cov[i][i] = comp["sigma_dex"] ** 2
    return assemble(u["scope"], comps, covariance=cov, covariance_source=u["covariance_source"],
                    assumptions=u["assumptions"], required_extra=u.get("required_extra", ()))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--element", action="append", required=True)
    ap.add_argument("--star", default="solar")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--extra-lines", default=None,
                    help="file listing further per-line CSVs (absolute paths) to try; each is "
                         "still accepted only if it reproduces the product's pool")
    a = ap.parse_args()
    tracked = [f for f in subprocess.run(["git", "ls-files", "*_lines.csv"], cwd=ROOT,
                                         capture_output=True, text=True).stdout.splitlines()]
    extra = (Path(a.extra_lines).read_text().split() if a.extra_lines else [])
    rc = 0
    for el in a.element:
        doc = pp.load(el, a.star)
        done, left = [], []
        for p in doc["products"]:
            if not is_two_line_scatter(p):
                continue
            comp, why = new_measurement(p, tracked, extra)
            if comp is None:
                left.append((pp.key_of(p), why))
                continue
            old = (p.get("sigma_stat"), p.get("sigma_asplund"))
            u = rebuild(p, comp)
            p["uncertainty"] = u
            p["sigma_stat"] = round(comp["sigma_dex"], 4)
            p["sigma_reported"] = u["sigma_reported"]
            p["restamped_at"] = pp._now()
            p["restamp_reason"] = REASON
            done.append((pp.key_of(p), old, comp["sigma_dex"]))
        print(f"== {el}: {len(done)} restamped, {len(left)} left")
        for k, (s0, a0), s1 in done:
            print(f"   {k}\n      stat {s0} -> {s1:.4f}  (asplund was {a0})")
        for k, why in left:
            print(f"   LEFT {k}: {why}")
        if left:
            rc = 1
        if done and not a.dry_run:
            doc["version"] = pp.bump(doc["version"])
            doc["updated_at"] = pp._now()
            out = pp.STORE / a.star / f"{el}.json"
            pp.write_feed(out, doc)
            new = json.loads(out.read_text())
            print(f"   wrote {out.relative_to(ROOT)} v{new['version']}")
    return rc


if __name__ == "__main__":
    sys.exit(main())

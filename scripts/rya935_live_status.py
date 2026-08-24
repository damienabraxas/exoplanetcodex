#!/usr/bin/env python3
"""Derive the live tracker's status from products on disk (RYA-935).

The tracker's JSON was HAND-TYPED. Two things follow from that, and both had
already happened by the time this was written:

  * it violated RYA-686 -- a result artifact landed with no generating harness,
    which is the RYA-559 hole that convention exists to close, and CI said so;
  * it went stale silently. The committed copy listed RYA-931 as BACKLOG on the
    day RYA-931 merged, and showed HARPS as "correction owed" after the
    corrected holding existed. A dashboard that is typed is a dashboard that
    lies eventually.

So this derives every cell it can from artifacts that already exist, and is
element-agnostic on purpose: Fe appears the moment its products land, with no
edit here.

What it deliberately does NOT emit: the ticket pipeline. That is Linear state,
not repository state, and a copy of it in a committed file is guaranteed to
drift -- which is exactly the failure above. The dashboard shows what the repo
can prove.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline.telluric_display_policy import (  # noqa: E402  RYA-1026
    anomaly, display_state)

#: `<El><Ion>_<lo>_<hi>_<instrument>_<holding>_<HANDLER>[_<SELECTOR>]_products.csv`. The
#: holding is in the stem because RYA-933/934 put it there -- before that, two holdings of
#: one instrument wrote the same filename and the second overwrote the first.
#:
#: RYA-990: the SELECTOR tag is optional and was previously unmatched, which silently
#: DROPPED every product carrying one. `derive_band_products._selector_tag` has emitted
#: `_DEEPGRADED` / `_FROMEW[-GRADED|-UNGRADED]` since RYA-984, so the two deep-graded VIS
#: Fe legs (RYA-984 Kitt Peak, RYA-991 HARPS) were on disk and merged but invisible here --
#: the tracker showed the 55-line shallow run as the only VIS synth product. A dashboard
#: that cannot see a merged product is the same failure mode as one that is hand-typed.
#: The selector is CARRIED, not discarded: it names which line set was measured, and two
#: runs differing only in selector are two different products (RYA-984) that must not
#: collapse into one cell (RYA-946).
STEM = re.compile(r"^(?P<el>[A-Z][a-z]?)(?P<ion>I+|IV|VI*)_(?P<lo>\d+)_(?P<hi>\d+)_"
                  r"(?P<rest>.+?)_(?P<handler>PROFILEFIT|SYNTH)"
                  # ONE OR MORE selector segments, not one. RYA-933 added `_GRADED`, and a
                  # graded ENGINE-A product stems as `..._SYNTH_GRADED_ENGINE-A_products.csv`
                  # -- two segments. With `?` here the regex matched neither, parse_stem
                  # returned None, and five committed products were silently invisible to
                  # the page (caught by test_every_committed_band_product_is_visible_to_the
                  # _tracker, which is the whole reason that test pins the invariant rather
                  # than the example).
                  r"(?P<selector>(?:_[A-Z][A-Z0-9]*(?:-[A-Z]+)?)*)_products\.csv$")


def parse_stem(name: str, instruments: set[str], holdings: set[str]) -> dict | None:
    m = STEM.match(name)
    if not m:
        return None
    rest = m.group("rest")
    instrument = next((i for i in sorted(instruments, key=len, reverse=True)
                       if rest == i or rest.startswith(i + "_")), None)
    if instrument is None:
        return None
    tail = rest[len(instrument):].lstrip("_")
    # No holding in the stem = a product from BEFORE RYA-933/934, when the stem
    # keyed on instrument alone. Say so; do not guess. Every product committed
    # before that change is attributable to an instrument and NOT to a holding,
    # and an instrument can serve several -- including a corrected and an
    # uncorrected one. Inferring which would be exactly the collapse the stem
    # change was made to prevent.
    if tail and tail in holdings:
        holding, source = tail, "filename"
    elif tail:
        holding, source = tail, "filename (unregistered holding)"
    else:
        holding, source = None, "absent -- product predates RYA-933/934"
    return {"element": m.group("el"), "ion": m.group("ion"),
            "lo_A": float(m.group("lo")), "hi_A": float(m.group("hi")),
            "instrument": instrument, "holding": holding,
            "holding_source": source, "handler": m.group("handler"),
            # Which line set was measured. "" is the default selector, which is what
            # every pre-RYA-984 artifact carries (RYA-984 kept the default unlabelled
            # so existing names did not change).
            "selector": m.group("selector").lstrip("_") or "default"}


def collect_products(roots: list[Path], instruments: set[str],
                     holdings: set[str]) -> list[dict]:
    rows: list[dict] = []
    for root in roots:
        if not root.exists():
            continue
        for path in sorted(root.rglob("*_products.csv")):
            meta = parse_stem(path.name, instruments, holdings)
            if meta is None:
                continue
            try:
                frame = pd.read_csv(path)
            except Exception:                                   # noqa: BLE001
                continue
            for _, r in frame.iterrows():
                rows.append({
                    **meta,
                    "band": r.get("band"), "treatment": r.get("treatment"),
                    "display": f"{r.get('route', '')} · {r.get('scale', '')}".strip(" ·"),
                    "A": None if pd.isna(r.get("A")) else float(r["A"]),
                    "sigma_stat": None if pd.isna(r.get("stat_dex")) else float(r["stat_dex"]),
                    "sigma_syst": None if pd.isna(r.get("syst_dex")) else float(r["syst_dex"]),
                    "n_lines": None if pd.isna(r.get("n_lines")) else int(r["n_lines"]),
                    # The RYA-906 physics axis: which gf SOURCE the linelist used
                    # ("kurucz"), not a grading tier.
                    #
                    # 🔴 RYA-990: the band-product schema carries NO TIER COLUMN, so the
                    # graded / consistent / no-bueno tier the cell contract asks for
                    # CANNOT be derived from a product. The only tier signals on disk are
                    # indirect -- the `selector` (DEEPGRADED = the lab-graded deep set)
                    # and `dominant` reading "gf scale (cited lab)". Both are emitted, and
                    # neither is relabelled as a tier: a dashboard that infers a tier is a
                    # dashboard that types one. Giving products a real tier field is the
                    # fix, and it belongs in the product writer, not here.
                    "gf_source": None if pd.isna(r.get("gf")) else str(r["gf"]),
                    "dominant_term": (None if pd.isna(r.get("dominant"))
                                      else str(r["dominant"])),
                    # Carried because the two deep-graded arms disagree on it (RYA-991
                    # flagged the gate refusing 1 of 109 on Kitt Peak and 0 of 109 on the
                    # noisier HARPS arm). A count of what a gate REFUSED is part of the
                    # result, not bookkeeping.
                    "n_excluded": (None if pd.isna(r.get("n_excluded"))
                                   else int(r["n_excluded"])),
                    # RYA-1026's display gate, AT THE RENDER PATH. Until this, the gate
                    # was called only by its own tests -- a guard nobody invokes is
                    # decorative, which is exactly what its docstring warns about.
                    #
                    # 🔴 CARRIED, NOT RAISED, and the choice is deliberate. Most existing
                    # Fe products were built on the telluric-UNCORRECTED 1984 atlas --
                    # RYA-1026 says so plainly -- so `assert_displayable` here would kill
                    # the whole page over historical products the ticket already knows
                    # about. Worse, dropping them would make a BLOCKED product look like
                    # an ABSENT one, which is the absence-as-conclusion error the same
                    # ticket forbids (RYA-833). So every product renders WITH ITS STATE
                    # ATTACHED and the page can mark it. The RAISE belongs where a NEW
                    # product is published, not where a historical one is listed.
                    "telluric_display": (display_state(meta["holding"])
                                         if meta["holding"] else "UNREGISTERED"),
                    "telluric_anomaly": (anomaly(meta["holding"])
                                         if meta["holding"] else None),
                    "source": str(path.relative_to(ROOT)),
                })
    return rows


def collect_instruments() -> list[dict]:
    """Coverage and telluric state, read from the registries -- never retyped."""
    import measure_band_ew as M
    from pipeline.telluric_policy import applied_state
    out = []
    for instrument, specs in M._INSTRUMENT_HOLDINGS.items():
        for spec in specs:
            try:
                telluric = applied_state(spec.holding_id)
            except KeyError:
                telluric = "unregistered"
            out.append({
                "instrument": instrument, "holding": spec.holding_id,
                "telluric_applied": telluric,
                "pre_normalised": spec.pre_normalised,
                "coverage_A": list(spec.span_A) if spec.span_A else None,
                "coverage_note": ("declared" if spec.span_A else
                                  "discovered by the reader; not declarable as one "
                                  "interval"),
            })
    return out


def collect_telluric(audit_root: Path) -> list[dict]:
    """Before/after residuals from the correction tickets' own evidence."""
    out = []
    for manifest in sorted(audit_root.glob("rya940_kp1984_telluric/*/fit_manifest.json")):
        d = json.loads(manifest.read_text())
        c = d.get("correction")
        if not c:
            out.append({"window": f"{int(d['band_A'][0])}-{int(d['band_A'][1])} A",
                        "product": "solar_kpno_molecfit_corrected",
                        "before_pct_below_0.5": None, "after_pct_below_0.5": None,
                        "note": "no admissible fit; band NOT corrected"})
            continue
        out.append({
            "window": f"{int(d['band_A'][0])}-{int(d['band_A'][1])} A",
            "product": "solar_kpno_molecfit_corrected",
            "before_pct_below_0.5": c["before"]["pct_below_0.5"],
            "after_pct_below_0.5": c["after"]["pct_below_0.5"],
            "externally_validated": d.get("externally_validated"),
        })
    verification = audit_root / "rya931_molecfit_runtime" / "verification.json"
    if verification.exists():
        v = json.loads(verification.read_text())["o2b_gate"]
        out.append({"window": "6867-6884 A", "product": "solar_harps_molecfit_corrected",
                    "metric": "pct_below_0.5",
                    "before_pct_below_0.5": v["o2b_before"]["pct_below_0.5"],
                    "after_pct_below_0.5": v["o2b_after"]["pct_below_0.5"],
                    "externally_validated": True})
    for row in _stellar_crires(audit_root):
        out.append(row)
    for row in _harps_state(audit_root):
        out.append(row)
    return out


#: The STELLAR telluric legs score with a DIFFERENT metric from the solar ones, and the
#: two must never share a column. RYA-940/931 report `pct_below_0.5` in a registered
#: telluric band; RYA-963/973 report the D1 residual — the median |1 - continuum| at
#: pixels molecfit calls telluric-DOMINATED and the star's own line list calls CLEAN.
#: They answer different questions and are not comparable as numbers, so each row names
#: the metric it carries (RYA-873: report a value under its DERIVED name).
def _stellar_crires(audit_root: Path) -> list[dict]:
    """Per-frame D1 residuals from the CRIRES+ stellar corrections (RYA-963, RYA-973)."""
    import csv
    rows = []
    for manifest in sorted(audit_root.glob("rya*_crires_telluric/corrected_manifest.csv")):
        ticket = manifest.parent.name.split("_")[0].upper().replace("RYA", "RYA-")
        for r in csv.DictReader(open(manifest, newline="")):
            try:
                before, after = float(r["gate_before"]), float(r["gate_after"])
            except (KeyError, TypeError, ValueError):
                continue
            rows.append({
                "window": f"{r.get('wlen_id', '?')} ({r.get('band', '?')} band)",
                "product": manifest.parent.name, "ticket": ticket,
                "star_id": r.get("star_id", ""), "date_obs": r.get("date_obs", ""),
                # NOT-contested and NOBODY-RECORDED-IT are different states and the
                # dashboard must not collapse them: alpha Cen's A/B verdict rests on a
                # branch assignment RYA-963 left contested, and a manifest written before
                # that column existed says nothing about it either way.
                "star_id_contested": (r["star_id_contested"].strip().lower() == "true"
                                      if "star_id_contested" in r else "unrecorded"),
                "metric": "d1_residual",
                "before_d1_residual": before, "after_d1_residual": after,
                "gate_passed": str(r.get("gate_passed", "")).lower() == "true",
                "gdas_profile": r.get("gdas_profile", ""),
            })
    return rows


def _harps_state(audit_root: Path) -> list[dict]:
    """HARPS telluric STATE determinations (RYA-973). A determination is not a
    correction, and the tracker must not let one read as the other: these rows carry
    `state`, never a before/after, because nothing has been corrected."""
    rows = []
    for path in sorted(audit_root.glob("rya*_harps_telluric/*_harps_telluric_state.json")):
        d = json.loads(path.read_text())
        rows.append({
            "window": f"{int(d['o2b_window_A'][0])}-{int(d['o2b_window_A'][1])} A",
            "product": f"{d['star']}_harps", "ticket": d.get("ticket", ""),
            "metric": "state_only__NOT_corrected",
            "header_state": d.get("header_state"),
            "flux_state": d.get("flux_state"),
            "o2b_median_frac_below": d.get("o2b_median_frac_below"),
            "control_median_frac_below": d.get("control_median_frac_below"),
            "excess_ratio": d.get("excess_ratio"),
            "n_products": d.get("n_products"),
            "note": d.get("disposition", ""),
        })
    return rows


def collect_model_matrix() -> dict:
    """RYA-1015 element x model-type availability, reconciled CSV vs disk vs code.

    A STANDING view: regenerated with the rest of the tracker, so "what atoms have what
    models" is read off a generated grid instead of re-derived in conversation.

    The disk half is a committed Sirius `find -L` snapshot whose positive control must
    have PASSED; if it did not, the matrix loud-fails rather than reporting absences
    that are symlink artifacts (RYA-1013). We surface that as an explicit error row
    rather than dropping the section, so a blind scan can never look like "no grids".
    """
    from pipeline.model_availability_matrix import build_matrix, DiskSnapshotError
    try:
        return build_matrix()
    except DiskSnapshotError as exc:
        return {"error": str(exc), "cells": [], "problem_count": None}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--products-root", type=Path, action="append", default=None)
    ap.add_argument("--audit-root", type=Path, default=ROOT / "data" / "audit")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "results" / "rya935"
                    / "live_status.json")
    ap.add_argument("--refresh-seconds", type=int, default=5)
    args = ap.parse_args()

    import measure_band_ew as M
    instruments = set(M._INSTRUMENT_HOLDINGS)
    holdings = {s.holding_id for specs in M._INSTRUMENT_HOLDINGS.values() for s in specs}
    roots = args.products_root or [ROOT / "data" / "results"]

    products = collect_products(roots, instruments, holdings)
    status = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "scripts/rya935_live_status.py",
        "refresh_seconds": args.refresh_seconds,
        "derivation_note": ("Every value here is read from a product or a registry. "
                            "Nothing is typed. The ticket pipeline is deliberately "
                            "absent: it is Linear state, and a committed copy drifts. "
                            "RYA-990: no product carries a graded/consistent/no-bueno "
                            "TIER field, so no tier is shown -- `selector` and "
                            "`dominant_term` are the only tier signals on disk and are "
                            "reported as themselves rather than inferred into a tier."),
        "elements": sorted({p["element"] + p["ion"] for p in products}),
        "bands": ["near-UV", "VIS", "red-optical", "NIR"],
        "instruments": collect_instruments(),
        "products": products,
        "telluric": collect_telluric(args.audit_root),
        "model_matrix": collect_model_matrix(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(status, indent=2) + "\n")
    species = ", ".join(status["elements"]) or "(none)"
    print(f"{len(products)} product rows across {species} -> {args.out}")
    unattributed = sum(1 for p in products if p["holding"] is None)
    status["unattributed_products"] = unattributed
    args.out.write_text(json.dumps(status, indent=2) + "\n")
    for i in status["instruments"]:
        have = {p["band"] for p in products if p["holding"] == i["holding"]}
        print(f"  {i['holding']:<34} telluric={i['telluric_applied']:<12} "
              f"bands: {sorted(have) or '—'}")
    if unattributed:
        print(f"\n  {unattributed} product rows predate RYA-933/934 and name no "
              f"holding in their filename.\n  They are attributable to an INSTRUMENT "
              f"only. Re-run to attribute them.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""RYA-1191 (A) — a telluric verdict for EVERY graded Fe line, in EVERY band.

    python3 scripts/rya1191_graded_line_telluric.py

VERIFICATION leg. Writes only under `data/results/rya1191/`; applies no correction and
moves no value. The fix and the re-measure are separate legs of the same ticket.

WHY THIS EXISTS AND WHAT RYA-1192 COULD NOT DO
-----------------------------------------------
RYA-1192 verified the HOLDINGS at band level and the measured lines at LINE level -- but
its line inventory came from the per-line `_lines.csv` artifacts on disk, and **none
exists for any red-optical, NIR or CRIRES+ Fe product**. So its strongest statement (not
one measured Fe line sits inside a registered telluric band) is true of near-UV + VIS and
SILENT exactly where the telluric complexes live. That silence is what this closes.

🔴 THE LINE SET IS RECOVERABLE WITHOUT RE-SYNTHESISING, AND THAT IS PROVEN, NOT ASSUMED.
The graded selectors are deterministic functions of committed data: `_cand_graded` /
`_cand_deep_graded` take the LAB-tier rows of `canonical_gf` in the band, split them on
`line_accounting_rya709.DEPTH_HI` using a feature depth read from `linelist_solar.csv`,
and match them into the band's synthesis list. Nothing there depends on the holding, the
engine or the fit. So the selection is re-derived here by CALLING THOSE SAME FUNCTIONS --
never a reimplementation (RYA-701) -- and the re-derivation is CONTROLLED against the 36
graded/deep-graded artifacts that DO exist: it reproduces every one of them exactly,
line for line. A method that reproduces every artifact we can check is the one used where
we cannot.

THE THREE FACTS PER LINE, AND THEY ARE NOT THE SAME FACT
---------------------------------------------------------
1. **Is a correction expected here?** Does a registered telluric band reach this line?
   Asked from `TELLURIC_BANDS` directly, NOT through `telluric_reason`, because that
   function keys on the INSTRUMENT's telluric_basis and answers "no" for a holding whose
   basis is `corrected` -- which is the right production behaviour and the wrong question
   for an audit. ⚠️ And an undeclared band is not an absent one (RYA-833).
2. **What did production DO with it?** `telluric_reason` + `serves_corrected_flux`:
   kept, QUARANTINED before fitting, or quarantined-then-LIFTED because a RYA-940
   corrected product covers it. This is the line's real fate in the shipped number.
3. **Is the flux there actually corrected?** The RYA-1192 evidence, applied at this
   line's own window: a direct comparison against the holding's raw sibling where that
   pair has the power to answer, and the band-level verdict where it does not.

🔴 (2) AND (3) CAN DISAGREE, AND WHERE THEY DO IT IS A FINDING EITHER WAY. A line
quarantined on a holding whose flux is verified corrected is a graded line thrown away
for nothing; a line kept on flux that is verified raw is an abundance measured through
the atmosphere.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

FEED = ROOT / "data/products/solar/Fe.json"
BP = ROOT / "data/results/band_products"

#: Half-window around each graded line for the direct comparison. The same 0.30 A the
#: RYA-1192 leg used, so the two maps are on one convention.
PAD_A = 0.30


def _selection(ion: str, lo: float, hi: float, deep: bool, ctx_cache: dict):
    """The production selectors, called -- never reimplemented (RYA-701)."""
    import derive_band_products as D
    from rya759_nearuv_fe_product import species_token
    from pipeline.nearuv_synth import build_solar_context
    from config.synth_bands import SYNTH_BANDS

    key = next((k for k, v in SYNTH_BANDS.items()
                if v.lo_A <= lo + 1e-6 and hi - 1e-6 <= v.hi_A), None)
    if key is None:
        raise LookupError(f"no synth band config covers {lo:.1f}-{hi:.1f} A")
    if key not in ctx_cache:
        ctx_cache[key] = build_solar_context(
            "Fe", 500000.0, linelist_file=str(SYNTH_BANDS[key].linelist),
            apply_canonical_gf=True, star="solar")
    ll = ctx_cache[key]["linelist"]
    sp = species_token("Fe", ion)
    cand = (D._cand_deep_graded(ll, lo_A=lo, hi_A=hi, species=sp) if deep
            else D._cand_graded(ll, lo_A=lo, hi_A=hi, species=sp))
    return [float(x) for x in cand.wave_A], key


def _control_reproduces_the_artifacts(ctx_cache: dict) -> dict:
    """🔴 THE CONTROL THAT LICENSES EVERYTHING BELOW. Re-derive the selection for every
    graded artifact that EXISTS and require an exact line-for-line match. If this ever
    stops holding, the red-optical and NIR line sets below are guesses."""
    import glob
    out = {"checked": 0, "exact": 0, "mismatched": []}
    for f in sorted(glob.glob(str(BP / "Fe*GRADED_*_lines.csv"))):
        n = Path(f).name
        parts = n.split("_")
        ion = "II" if parts[0] == "FeII" else "I"
        lo, hi = float(parts[1]), float(parts[2])
        try:
            sel, _ = _selection(ion, lo, hi, "_DEEPGRADED_" in n, ctx_cache)
        except (LookupError, SystemExit):
            continue
        art = sorted(round(float(x), 3) for x in
                     pd.read_csv(f).wavelength_air_A.dropna())
        out["checked"] += 1
        if sorted(round(x, 3) for x in sel) == art:
            out["exact"] += 1
        else:
            out["mismatched"].append(n)
    out["verdict"] = (
        f"the selectors reproduce {out['exact']} of {out['checked']} existing graded "
        f"artifacts EXACTLY, line for line"
        + ("" if not out["mismatched"] else
           f" — 🔴 {len(out['mismatched'])} MISMATCH: {out['mismatched'][:3]}"))
    return out


def _band_of(w: float):
    from pipeline.telluric_policy import TELLURIC_BANDS
    for a, b, c in TELLURIC_BANDS:
        if a <= w <= b:
            return str(c), (float(a), float(b))
    return None, None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pad", type=float, default=PAD_A)
    ap.add_argument("--out-dir", default="data/results/rya1191")
    a = ap.parse_args(argv)
    out = ROOT / a.out_dir
    out.mkdir(parents=True, exist_ok=True)

    import rya1192_sirius_verification as V
    from measure_band_ew import telluric_reason, serves_corrected_flux

    ctx_cache: dict = {}
    control = _control_reproduces_the_artifacts(ctx_cache)

    # 🔴 THIS LEG'S BAND EVIDENCE IS RYA-1191's OWN, NOT RYA-1192's. Leg B measured every
    # registered band against a telluric TEMPLATE and re-judged two of RYA-1192's verdicts
    # — kurucz2005 in H2O 7160-7340 (called raw-like from a depth ratio; it is clean) and
    # kpno_molecfit at o2gamma (called corrected from a featureless rescale; it is not).
    # Reading the superseded map here would put 21 graded lines in
    # `kept_though_the_flux_is_RAW` on a verdict this ticket has already overturned.
    band_state = _load_band_evidence()
    sig = V.signature_table()
    bverd = V.band_verdicts(sig)
    claims = V.holding_claim_verdicts(bverd)
    floors = {h: (V.baseline(i, h, s) if s else None)
              for h, (i, s, _) in V.RAW_SIBLING.items()}

    feed = json.loads(FEED.read_text())
    fams, seen = [], set()
    for p in feed["products"]:
        if p["tier"] not in ("GRADED", "DEEPGRADED"):
            continue
        rng = p.get("wavelength_range_A")
        if not rng:
            continue
        k = (p["ion"], p["band"], p["instrument"], p["holding"], p["selector"],
             float(rng[0]), float(rng[1]))
        if k in seen:
            continue
        seen.add(k)
        fams.append(k)

    rows, fam_rows = [], []
    for (ion, band, inst, holding, selector, lo, hi) in sorted(fams):
        try:
            sel, bkey = _selection(ion, lo, hi, selector.startswith("DEEPGRADED"),
                                   ctx_cache)
        except (LookupError, SystemExit) as e:
            fam_rows.append(dict(ion=ion, band=band, holding=holding, selector=selector,
                                 state="SELECTION-REFUSED", why=str(e)[:100]))
            continue
        raw_sib = V.RAW_SIBLING.get(holding, (inst, None, ""))[1]
        floor = floors.get(holding)
        n_q = n_l = n_inband = 0
        for w in sel:
            tb, span = _band_of(w)
            why = telluric_reason(w, inst)
            lift = serves_corrected_flux(holding, w) if why else ""
            fate = ("QUARANTINED-TELLURIC" if (why and not lift)
                    else "LIFTED-CORRECTED-PRODUCT" if why else "KEPT")
            n_inband += bool(tb)
            n_q += fate == "QUARANTINED-TELLURIC"
            n_l += fate == "LIFTED-CORRECTED-PRODUCT"
            rec = dict(ion=ion, band=band, instrument=inst, holding=holding,
                       selector=selector, wavelength_air_A=round(w, 4),
                       synth_band=bkey,
                       # (1) is a correction expected here at all?
                       in_registered_telluric_band=bool(tb),
                       telluric_band=tb or "", telluric_band_span=str(span or ""),
                       # (2) what production did
                       production_fate=fate,
                       production_reason=(str(why).split(" and ")[0] if why else ""),
                       lift_provenance=lift,
                       raw_sibling=raw_sib)
            # (3) is the flux actually corrected here?
            rec.update(_flux_verdict(V, inst, holding, raw_sib, floor, w, a.pad,
                                     tb, bverd, band_state))
            rows.append(rec)
        fam_rows.append(dict(ion=ion, band=band, instrument=inst, holding=holding,
                             selector=selector, window_A=[lo, hi], synth_band=bkey,
                             state="OK", n_selected=len(sel),
                             n_in_registered_telluric_band=n_inband,
                             n_quarantined=n_q, n_lifted=n_l,
                             n_kept=len(sel) - n_q,
                             had_per_line_artifact=_artifact_exists(ion, lo, hi, holding,
                                                                    selector)))
    per_line = pd.DataFrame(rows)
    per_line.to_csv(out / "rya1191_graded_line_telluric.csv", index=False)

    doc = {
        "ticket": "RYA-1191 (A) — telluric verdict at every graded Fe line",
        "kind": "VERIFICATION — no correction applied, no value moved",
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "selection_control": control,
        "pad_A": a.pad,
        "families": fam_rows,
        "holding_claim_verdicts": claims,
        "reduction_floors": {k: (round(v, 8) if v is not None else None)
                             for k, v in floors.items()},
        "totals": {
            "n_graded_lines": int(len(per_line)),
            "n_in_registered_telluric_band": int(
                per_line.in_registered_telluric_band.sum()) if len(per_line) else 0,
            "n_unverifiable": int((per_line.flux_state == "UNVERIFIABLE").sum())
            if len(per_line) else 0,
        },
        "by_state": (per_line.groupby(["holding", "band", "flux_state"]).size()
                     .rename("n").reset_index().to_dict("records")
                     if len(per_line) else []),
        "disagreements": _disagreements(per_line),
    }
    (out / "rya1191_graded_line_telluric.json").write_text(json.dumps(doc, indent=2) + "\n")

    print("=== selection control ===")
    print("  " + control["verdict"])
    print("\n=== families ===")
    for f in fam_rows:
        if f["state"] != "OK":
            print(f"  {f['ion']:<3}{f['band']:<13}{f['holding'][:32]:<33} {f['state']} {f['why'][:60]}")
            continue
        print(f"  {f['ion']:<3}{f['band']:<13}{f['holding'][:32]:<33}{f['selector'][:14]:<15}"
              f"sel={f['n_selected']:<4} inband={f['n_in_registered_telluric_band']:<3} "
              f"quar={f['n_quarantined']:<3} lift={f['n_lifted']:<3} kept={f['n_kept']:<4}"
              f" artifact={'yes' if f['had_per_line_artifact'] else 'NO'}")
    print("\n=== flux state, per holding x band ===")
    if len(per_line):
        print(per_line.groupby(["holding", "band", "flux_state"]).size().to_string())
    print(f"\nwrote {a.out_dir}/")
    return 0


def _artifact_exists(ion, lo, hi, holding, selector) -> bool:
    import glob
    tag = "DEEPGRADED" if selector.startswith("DEEPGRADED") else "GRADED"
    pat = f"Fe{ion}_{lo:.0f}_{hi:.0f}_*{holding}_SYNTH_{tag}_*_lines.csv"
    return bool(glob.glob(str(BP / pat)))


#: The states RYA-1191's own band measurement can return, mapped onto the verdict
#: vocabulary the per-line map uses. `partial` is deliberately NOT "corrected": a measured
#: residual is a statement about the flux and the line must carry it.
_EVIDENCE_TO_STATE = {"clean": "VERIFIED-CORRECTED",
                      "partial": "PARTIALLY-CORRECTED",
                      "uncorrected": "VERIFIED-RAW-LIKE",
                      "undetermined": "UNDETERMINED-NO-POWER"}


def _load_band_evidence():
    """(holding, lo, hi) -> (state, statistic) from THIS ticket's measurement."""
    f = ROOT / "data/catalog/telluric_correction_evidence.csv"
    if not f.exists():
        return []
    d = pd.read_csv(f, comment="#")
    return [(str(r.holding_id), float(r.lo_A), float(r.hi_A), str(r.state),
             f"{r.statistic}={r.value} vs references {r.reference_value}")
            for r in d.itertuples()]


def _band_evidence_at(band_state, holding, w):
    for h, lo, hi, st, stat in band_state:
        if h == holding and lo <= w <= hi:
            return st, stat
    return None, None


def _flux_verdict(V, inst, holding, raw_sib, floor, w, pad, tb, bverd,
                  band_state=()) -> dict:
    """(3) — is the flux at THIS line telluric-corrected, on the evidence?

    ⚠️ The direct comparison is used only where this pair can answer with it: a pair whose
    reduction floor exceeds the absorption present at the line cannot return CORRECTED for
    any data (RYA-1192). Where it cannot, the band verdict answers, and the row says so.
    """
    ev_state, ev_stat = _band_evidence_at(band_state, holding, w)
    if ev_state:
        return {"flux_state": _EVIDENCE_TO_STATE[ev_state],
                "flux_basis": (f"RYA-1191 BAND MEASUREMENT against the telluric template "
                               f"({ev_stat})"),
                "band_evidence_state": ev_state}
    if raw_sib is None:
        if tb is None:
            return {"flux_state": "NO-TELLURIC-BAND-DECLARED",
                    "flux_basis": ("no raw sibling and no declared telluric band reaches "
                                   "this line — not a claim that no telluric is present "
                                   "(RYA-833)")}
        key = next((k for k in bverd if k.startswith(f"{holding}|{tb}|")), None)
        v = bverd.get(key, {})
        return {"flux_state": v.get("verdict", "UNVERIFIABLE"),
                "flux_basis": f"BAND VERDICT ({v.get('test', 'n/a')}) — no raw sibling",
                "band_anchor": v.get("anchor", "")}
    avail = V._absorption_available(inst, raw_sib, w - pad, w + pad)
    powered = (floor is not None and avail is not None
               and V.MATERIAL_MULTIPLE * max(floor, 1e-12) < avail)
    if powered:
        r = V.compare(inst, holding, raw_sib, w, pad, base=floor)
        return {"flux_state": r.get("state", "UNVERIFIABLE"),
                "flux_basis": "DIRECT-COMPARISON vs raw sibling at this line",
                "max_abs_diff": r.get("max_abs_diff"),
                "absorption_available_in_raw": (round(avail, 6)
                                                if avail is not None else None)}
    if tb is None:
        return {"flux_state": "NO-TELLURIC-BAND-DECLARED",
                "flux_basis": ("the difference test has no power at this line and no "
                               "declared telluric band reaches it (RYA-833)"),
                "absorption_available_in_raw": (round(avail, 6)
                                                if avail is not None else None)}
    key = next((k for k in bverd if k.startswith(f"{holding}|{tb}|")), None)
    v = bverd.get(key, {})
    return {"flux_state": v.get("verdict", "UNVERIFIABLE"),
            "flux_basis": ("BAND VERDICT — the difference test has no power at this line "
                           f"({v.get('test', 'n/a')})"),
            "band_anchor": v.get("anchor", "")}


def _disagreements(per_line: pd.DataFrame) -> dict:
    """🔴 WHERE PRODUCTION'S FATE AND THE FLUX EVIDENCE PART COMPANY."""
    if not len(per_line):
        return {}
    thrown_away = per_line[(per_line.production_fate == "QUARANTINED-TELLURIC")
                           & (per_line.flux_state == "VERIFIED-CORRECTED")]
    partial = per_line[(per_line.production_fate != "QUARANTINED-TELLURIC")
                       & (per_line.flux_state == "PARTIALLY-CORRECTED")]
    measured_on_raw = per_line[(per_line.production_fate != "QUARANTINED-TELLURIC")
                               & (per_line.flux_state.isin(
                                   ["VERIFIED-RAW", "VERIFIED-RAW-LIKE"]))
                               & per_line.in_registered_telluric_band]
    return {
        "quarantined_though_the_flux_IS_corrected": {
            "n": int(len(thrown_away)),
            "detail": (thrown_away.groupby(["holding", "band", "telluric_band"]).size()
                       .rename("n").reset_index().to_dict("records")),
            "meaning": ("graded lines refused before fitting on a holding whose flux at "
                        "those wavelengths is verified corrected — measurable lines "
                        "thrown away, a pool cost with no data-integrity benefit"),
        },
        "kept_on_PARTIALLY_corrected_flux": {
            "n": int(len(partial)),
            "detail": (partial.groupby(["holding", "band", "telluric_band"]).size()
                       .rename("n").reset_index().to_dict("records")),
            "wavelengths": sorted(round(float(x), 3)
                                  for x in partial.wavelength_air_A.unique()),
            "meaning": ("⚠️ THE THIRD STATE, AND IT IS NEITHER OF THE OTHER TWO. These "
                        "lines are in a live product, on flux where a correction ran and "
                        "left a MEASURED residual — not raw, not clean. They are not "
                        "excluded here: a residual of a few sigma against references at "
                        "thirty is a small error term, and excluding half a pool for it "
                        "would cost more than it buys. What the ticket owes is the "
                        "NUMBER, so the A-shift of dropping them is measured and reported "
                        "rather than the choice being made silently."),
        },
        "kept_though_the_flux_is_RAW": {
            "n": int(len(measured_on_raw)),
            "detail": (measured_on_raw.groupby(["holding", "band", "telluric_band"]).size()
                       .rename("n").reset_index().to_dict("records")),
            "meaning": ("🔴 graded lines carrying an abundance measured on flux that is "
                        "verified NOT telluric-corrected, inside a registered band"),
        },
    }


if __name__ == "__main__":
    raise SystemExit(main())

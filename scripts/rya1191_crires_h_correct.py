#!/usr/bin/env python3
"""RYA-1191 (C, fix) — actually telluric-correct the CRIRES+ H arm, on OUR Vesta IDPs.

    python3 scripts/rya1191_crires_h_correct.py --plan     # what would run, no molecfit
    python3 scripts/rya1191_crires_h_correct.py --run      # run molecfit + verify

🔴 WHY THIS EXISTS RATHER THAN ANOTHER MEASUREMENT
---------------------------------------------------
RYA-1192 flagged three 100 A windows of the CRIRES+ H arm. RYA-1191 (C) re-asked it per
graded line and could not settle 13 of the 25: the per-line window is the right UNIT for
an abundance fitted over 1.89 A, and it has too few telluric pixels for a shape test to
clear its own null. Measuring harder was not going to fix that — the answer is to CORRECT
the arm and then ask again, which is what the permanent IR rule wants anyway.

Ryan's ruling (2026-09-07) is that the SCIENCE instruments are the assignment: HARPS and
CRIRES+ are what run on Alpha Cen, and a line rescued on a solar-only atlas is calibrating
a path the science targets do not have. CRIRES+ H is the IR half of that.

WHAT MAKES THIS RUNNABLE, ALL FOUR CHECKED BEFORE A LINE OF IT WAS WRITTEN
--------------------------------------------------------------------------
* **esorex + molecfit are installed** at /srv/codex/eso/molecfit/bin (model, calctrans,
  correct). `pipeline.telluric.esorex_runtime.resolve_esorex` finds them.
* **A real per-night GDAS profile exists for the Vesta night.** `fetch_gdas('paranal',
  mjd=...)` returns the 2022-11-22T00 profile. ⚠️ This matters more than it looks: a
  SILENT standard-atmosphere fallback is the RYA-373 critical bug, and `_resolve_gdas`
  raises rather than falling back.
* **The molecfit step is already generic.** `_molecfit_segment` is documented as reused
  "for the science order, the RV-anchor order, and every alpha Cen chip" — it takes a
  frame, ONE segment and a molecule list. Only `_molecfit_driver` was K-specific, and
  only because it hardcodes `segment_at(CO_2_0_BANDHEAD_NM)`.
* **Runtime is minutes.** RYA-963's H1559 and H1582 alpha Cen runs took 2m42s and 5m47s
  end to end; this is six Vesta frames, not an overnight job.

WHAT THIS DOES *NOT* DO, SAID PLAINLY
--------------------------------------
⚠️ It corrects the arm in the TOPOCENTRIC frame and verifies it there. That is the right
frame for the telluric question (tellurics are at rest in it) and the same frame the
RYA-1191 template test uses. It does NOT RV-condition to the solar rest frame and does NOT
emit a new Fe product: the reflected-solar RV anchors in `crires_telluric` are K-band
lines, and there is no live CRIRES+ H Fe product to re-measure (RYA-1094 derived one but
it is in no pool of data/products/solar/Fe.json). Building the H product is the next step
and is not smuggled in here.
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

#: The graded-Fe span of the H arm (RYA-1094): the lab-graded Ruffoni pool runs
#: 15051.7-17277.5 A. Segments outside it are not corrected — there is no graded line
#: there to protect, and each segment is a molecfit run.
H_LO_A, H_HI_A = 15007.0, 17494.0
MOLECULES = ("H2O", "CH4", "CO2")


def h_frames(crires_dir=None):
    from pipeline.crires_telluric import inventory, VESTA_CRIRES_DIR
    frames = inventory(crires_dir or VESTA_CRIRES_DIR)
    return [f for f in frames if f.band == "H"]


def graded_lines() -> np.ndarray:
    import pandas as pd
    cg = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    m = ((cg.species.astype(str) == "Fe I")
         & cg.wavelength_air_A.between(H_LO_A, H_HI_A)
         & (cg.gf_tier.astype(str) == "LAB"))
    return np.sort(cg.loc[m, "wavelength_air_A"].to_numpy(float))


def segments_with_graded(frame, lines) -> list:
    """Only the chip segments that actually carry a graded line.

    Each segment is one molecfit fit, and a CRIRES+ frame has dozens. Correcting the ones
    with nothing in them would multiply the runtime for no science.
    """
    out = []
    for s in frame.segments:
        w = np.asarray(s.wave_A, float)
        w = w[np.isfinite(w)]
        if w.size < 50:
            continue
        lo, hi = float(w.min()), float(w.max())
        if hi < H_LO_A or lo > H_HI_A:
            continue
        n = int(((lines >= lo) & (lines <= hi)).sum())
        if n:
            out.append((s, n, lo, hi))
    return out


def covering_set(plan_pairs, lines):
    """A minimal set of segments covering every graded line — greedy set cover.

    ⚠️ WHY NOT CORRECT EVERYTHING. Six H frames across four settings carry ~50 segments
    holding a graded line, and each segment is one molecfit fit at 3-6 minutes. Most of
    those segments are re-observations of the SAME lines, so correcting all of them buys
    duplicate answers and costs hours. The question here is per LINE, so one good segment
    per line answers it; the covering set is chosen by how many graded lines a segment
    holds and then by how many pixels it has, never by which answer it gives.
    """
    need = {round(float(x), 4) for x in lines}
    ranked = sorted(plan_pairs,
                    key=lambda t: (-t[1], -int(np.isfinite(np.asarray(t[0].wave_A)).sum())))
    chosen = []
    for seg, n, lo, hi, frame in ranked:
        got = {x for x in need if lo <= x <= hi}
        if not got:
            continue
        chosen.append((seg, len(got), lo, hi, frame))
        need -= got
        if not need:
            break
    return chosen, sorted(need)


def verify(work: Path, out: Path) -> int:
    """Re-ask RYA-1191 (C)'s question on the CORRECTED flux, and settle the 13.

    🔴 THE POINT OF THE WHOLE LEG. (C) could not judge 13 of the 25 graded H lines: the
    per-line window is the right unit for an abundance fitted over 1.89 A and it holds too
    few telluric pixels for a shape test to clear its own null. That is a limit of asking,
    not of the data — so the arm was corrected and the question re-asked here.

    The statistic is the one molecfit itself fits: MTRANS, the transmission it removed.
    For each graded line, the RAW flux must show absorption where MTRANS is deep (the
    positive control, which is what (C) could not get) and the CORRECTED flux must not.
    """
    import glob
    files = sorted(glob.glob(str(work / "*_corrected.npz")))
    if not files:
        raise SystemExit(f"no corrected segments under {work} — run --run first")
    lines = graded_lines()
    rows = []
    for f in files:
        d = np.load(f)
        w, raw, cor, mt = (np.asarray(d[k], float)
                           for k in ("wave_A", "flux_raw", "flux_corr", "mtrans"))
        o = np.argsort(w); w, raw, cor, mt = w[o], raw[o], cor[o], mt[o]
        tag = Path(f).name.replace("_corrected.npz", "")
        for lam in lines[(lines >= w.min()) & (lines <= w.max())]:
            m = np.abs(w - lam) <= FIT_HALF_WIDTH_A
            if m.sum() < 20:
                rows.append({"wavelength_air_A": round(float(lam), 4), "segment": tag,
                             "state": "TOO-THIN", "n": int(m.sum())})
                continue
            deep, clean = mt[m] < 0.90, mt[m] > 0.995
            rec = {"wavelength_air_A": round(float(lam), 4), "segment": tag,
                   "n_pixels": int(m.sum()),
                   "mtrans_min": round(float(np.nanmin(mt[m])), 4),
                   "n_deep_pixels": int(deep.sum())}
            if deep.sum() < 5 or clean.sum() < 5:
                rec.update(state="NO-TELLURIC-IN-WINDOW",
                           why=(f"only {int(deep.sum())} pixel(s) below MTRANS 0.90 in "
                                f"this line's +/-{FIT_HALF_WIDTH_A} A window — molecfit "
                                f"found nothing here to remove"))
                rows.append(rec); continue
            def contrast(flux):
                fl = flux[m]
                c = np.nanmedian(fl[fl >= np.nanpercentile(fl, 80)])
                a = 1.0 - fl / c
                return float(np.nanmean(a[deep]) - np.nanmean(a[clean]))
            rec["contrast_raw"] = round(contrast(raw), 4)
            rec["contrast_corrected"] = round(contrast(cor), 4)
            # ⚠️ The control is the RAW frame at THIS line, not the arm. A line whose raw
            # flux shows no telluric contrast cannot be cleared by a quiet corrected one.
            if rec["contrast_raw"] < 0.02:
                rec.update(state="UNDETERMINED",
                           why=(f"the RAW flux shows only {rec['contrast_raw']:+.4f} "
                                f"contrast at the telluric pixels — no positive control"))
            elif rec["contrast_corrected"] > 0.25 * rec["contrast_raw"]:
                rec.update(state="RESIDUAL",
                           why=(f"corrected retains {rec['contrast_corrected']:+.4f} of "
                                f"the raw {rec['contrast_raw']:+.4f}"))
            else:
                rec.update(state="CORRECTED",
                           why=(f"telluric contrast {rec['contrast_raw']:+.4f} -> "
                                f"{rec['contrast_corrected']:+.4f}"))
            rows.append(rec)
    d = pd.DataFrame(rows)
    d.to_csv(out / "rya1191_crires_h_verified.csv", index=False)
    best = (d[d.state.isin(["CORRECTED", "RESIDUAL", "NO-TELLURIC-IN-WINDOW"])]
            .sort_values("state").drop_duplicates("wavelength_air_A", keep="first"))
    doc = {"ticket": "RYA-1191 (C, verify) — the H arm after correction",
           "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "n_graded_lines": int(len(lines)),
           "n_with_a_verdict": int(len(best)),
           "by_state": d.state.value_counts().to_dict(),
           "per_line_best": best.to_dict("records"),
           "note": ("⚠️ A line measured on several segments keeps its best verdict; the "
                    "per-segment rows are all in the CSV. NO-TELLURIC-IN-WINDOW is a "
                    "statement about molecfit's own fitted transmission, so it needs no "
                    "positive control — there is nothing there to remove.")}
    (out / "rya1191_crires_h_verified.json").write_text(json.dumps(doc, indent=2) + "\n")
    print(d.state.value_counts().to_string())
    print(f"\n  {len(best)} of {len(lines)} graded H lines have a verdict")
    print(best[["wavelength_air_A", "segment", "state", "contrast_raw",
                "contrast_corrected"]].to_string(index=False)
          if "contrast_raw" in best else best.to_string(index=False))
    print(f"\nwrote {out.name}/")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--plan", action="store_true", help="report what would run; no molecfit")
    g.add_argument("--run", action="store_true", help="run molecfit on every planned segment")
    g.add_argument("--verify", action="store_true",
                   help="re-ask the per-line telluric question on the CORRECTED flux")
    ap.add_argument("--all-segments", action="store_true",
                    help="correct EVERY segment carrying a graded line, not a covering set")
    ap.add_argument("--work", default="/tmp/rya1191_crires_h")
    ap.add_argument("--out-dir", default="data/results/rya1191")
    a = ap.parse_args(argv)
    out = ROOT / a.out_dir; out.mkdir(parents=True, exist_ok=True)
    work = Path(a.work); work.mkdir(parents=True, exist_ok=True)

    if a.verify:
        import pandas as pd  # noqa: F401  (used by verify)
        return verify(work, out)
    lines = graded_lines()
    frames = h_frames()
    plan, total = [], 0
    for f in frames:
        segs = segments_with_graded(f, lines)
        total += len(segs)
        plan.append({"file": f.path.name, "wlen_id": f.wlen_id, "mjd": float(f.mjd),
                     "specsys": f.specsys, "n_segments_with_graded": len(segs),
                     "segments": [{"order": int(s.order), "detector": int(s.detector),
                                   "lo_A": round(lo, 2), "hi_A": round(hi, 2),
                                   "n_graded": n} for s, n, lo, hi in segs]})
    print(f"=== CRIRES+ H frames: {len(frames)}  segments carrying a graded line: {total} ===")
    for p in plan:
        print(f"  {p['file']}  {p['wlen_id']}  mjd={p['mjd']:.5f}  "
              f"{p['n_segments_with_graded']} segment(s)")
        for s in p["segments"]:
            print(f"      order {s['order']:>3} det {s['detector']}  "
                  f"{s['lo_A']:.1f}-{s['hi_A']:.1f} A   {s['n_graded']} graded")
    served = set()
    for p in plan:
        for s in p["segments"]:
            served |= {float(x) for x in lines if s["lo_A"] <= x <= s["hi_A"]}
    print(f"\n  graded lines in the H arm: {len(lines)}; on a correctable segment: {len(served)}")
    missing = sorted(set(float(x) for x in lines) - served)
    if missing:
        print(f"  ⚠️ NOT on any segment ({len(missing)}): {[round(x,3) for x in missing]}")

    doc = {"ticket": "RYA-1191 (C, fix) — CRIRES+ H telluric correction",
           "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "molecules": list(MOLECULES), "h_span_A": [H_LO_A, H_HI_A],
           "n_graded_lines": int(len(lines)),
           "n_graded_on_a_correctable_segment": len(served),
           "graded_not_on_any_segment": [round(x, 4) for x in missing],
           "plan": plan}
    if a.plan:
        (out / "rya1191_crires_h_correction_plan.json").write_text(
            json.dumps(doc, indent=2) + "\n")
        print(f"\nwrote {a.out_dir}/rya1191_crires_h_correction_plan.json (PLAN ONLY)")
        return 0

    from pipeline.crires_telluric import _molecfit_segment
    pairs = [(sg, n, lo, hi, f) for f in frames
             for sg, n, lo, hi in segments_with_graded(f, lines)]
    if a.all_segments:
        todo, uncovered = pairs, []
    else:
        todo, uncovered = covering_set(pairs, lines)
        print(f"\n  covering set: {len(todo)} of {len(pairs)} segments cover all "
              f"{len(lines) - len(uncovered)} graded lines")
        if uncovered:
            print(f"  ⚠️ uncovered: {uncovered}")
    results = []
    for sg, n, lo, hi, f in todo:
        s = sg
        if True:
            tag = f"{f.wlen_id}_o{int(s.order)}d{int(s.detector)}"
            print(f"\n  [molecfit] {f.path.name} {tag} {lo:.1f}-{hi:.1f} A ({n} graded)")
            rec = {"file": f.path.name, "wlen_id": f.wlen_id, "tag": tag,
                   "order": int(s.order), "detector": int(s.detector),
                   "lo_A": round(lo, 2), "hi_A": round(hi, 2), "n_graded": n}
            try:
                r = _molecfit_segment(f, s, work / tag, MOLECULES)
            except Exception as e:
                rec.update(state="FAILED", error=f"{type(e).__name__}: {str(e)[:200]}")
                print(f"      FAILED {rec['error']}")
                results.append(rec); continue
            mt = np.asarray(r["mtrans"], float)
            rec.update(state="CORRECTED",
                       mtrans_min=round(float(np.nanmin(mt)), 5),
                       mtrans_frac_below_097=round(float(np.nanmean(mt < 0.97)), 4))
            rec["chi2"] = (round(float(r["chi2"]), 4) if r.get("chi2") is not None else None)
            rec["gdas"] = str(r.get("gdas"))
            rec["npix"] = int(mt.size)
            np.savez_compressed(work / f"{tag}_corrected.npz",
                                wave_A=r["lam_A"], flux_corr=r["corr"],
                                flux_raw=r["flux_raw"], mtrans=mt)
            print(f"      OK mtrans min={rec['mtrans_min']} "
                  f"frac<0.97={rec['mtrans_frac_below_097']} chi2={rec['chi2']}")
            results.append(rec)
    doc["results"] = results
    doc["n_corrected"] = sum(1 for r in results if r["state"] == "CORRECTED")
    doc["n_failed"] = sum(1 for r in results if r["state"] == "FAILED")
    (out / "rya1191_crires_h_correction.json").write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\n  corrected {doc['n_corrected']} / failed {doc['n_failed']}")
    print(f"wrote {a.out_dir}/rya1191_crires_h_correction.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

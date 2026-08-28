#!/usr/bin/env python3
"""
scripts/rya1098_validate_amarsi_models.py
=========================================
RYA-1098 — prove the Amarsi 3D-NLTE model load is correct, or say it is not.

    python3 scripts/rya1098_validate_amarsi_models.py

🔴 THE QUESTION IS NOT "DID IT LOAD". The vendored MLP pickles were written by
scikit-learn 1.0.2 and are loaded under 1.9.0. Cross-version unpickling of sklearn
estimators is unsupported and does not reliably raise -- a clean deserialize followed by
subtly wrong predictions is the documented failure mode, and these models output a PHYSICS
CORRECTION applied to Fe abundances. So a silent mis-load is an abundance error with no
stack trace.

WHAT THE TICKET ASKED FOR, AND WHAT IS ACTUALLY POSSIBLE
--------------------------------------------------------
The preferred route was to capture (input -> output) pairs in the ORIGINAL 1.0.2
environment and check the current path reproduces them. That route is CLOSED and saying so
is part of the answer: sklearn 1.0.2 predates Python 3.12 and is not installable in any
environment we hold (Mac venv 1.9.0, Mac system 1.6.1, Sirius venv312 1.9.0). Inventing a
substitute "1.0.2 environment" would be validating the load against itself.

So the validation uses ground truth that never passed through our environment at all:

  CHECK 1 — THE AUTHORS' OWN PUBLISHED PAIR. `vendor/1L-3NErrors/README.md` §2.4 prints a
  worked example computed under 1.0.2 by the model's authors.

  CHECK 2 — AN INDEPENDENT REIMPLEMENTATION. The learned parameters are plain numpy arrays
  and pickle round-trips those exactly regardless of sklearn version, so the PARAMETERS
  cannot have been corrupted -- only their USE could be. Reimplementing the documented
  forward pass in numpy and comparing to `predict()` tests exactly that.

  CHECK 3 — THE CONTAMINATION QUESTION, ON THE FEED'S OWN INPUTS. Check 2 is repeated on
  the feature vectors the published runs actually spanned: every in-domain line of every
  Amarsi product, crossed with the whole A(Fe;3N) axis range the per-line iteration can
  visit. A general agreement is reassuring; this one answers "are the numbers already in
  the feed affected".

Emits `data/results/rya1098/rya1098_amarsi_model_integrity.json`.
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.amarsi_model_integrity import (            # noqa: E402
    MODEL_FILES, VALIDATED_RUNTIMES, VENDOR_DIR, forward, load_models,
    runtime_version, validate, written_version)

#: The A(Fe;3N) grid axis. The per-line solve iterates on this, so the feed's corrections
#: were evaluated somewhere inside it and the contamination check must span all of it.
AXIS_LO, AXIS_HI, AXIS_STEP = 4.5, 7.5, 0.05


def feed_vectors(per_line_csvs: list[Path]) -> pd.DataFrame:
    frames = []
    for f in per_line_csvs:
        d = pd.read_csv(f)
        g = d[(d.element == "Fe") & (d.band == "VIS") & d.in_domain.map(bool)]
        frames.append(g.dropna(subset=["aberr"])[["ion", "elo_eV", "eup_eV", "loggf"]])
    if not frames:
        return pd.DataFrame(columns=["ion", "elo_eV", "eup_eV", "loggf"])
    return pd.concat(frames, ignore_index=True)


def contamination_check(lines: pd.DataFrame, *, teff=5772.0, logg=4.438,
                        vmic=1.0) -> dict:
    models = load_models(require_validated=False)
    axes = np.round(np.arange(AXIS_LO, AXIS_HI + 0.5 * AXIS_STEP, AXIS_STEP), 3)
    X = np.vstack([
        np.column_stack([np.full(len(lines), teff), np.full(len(lines), logg),
                         np.full(len(lines), a), np.full(len(lines), vmic),
                         lines.elo_eV.to_numpy(float), lines.eup_eV.to_numpy(float),
                         lines.loggf.to_numpy(float)])
        for a in axes])
    worst, per_net = 0.0, {}
    for key in ("lt02", "gt02"):
        sc, m = models[key]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sk = np.asarray(m.predict(sc.transform(X)), dtype=float).ravel()
        d = float(np.max(np.abs(sk - forward(sc, m, X))))
        per_net[key] = d
        worst = max(worst, d)
    return {"n_lines": int(len(lines)), "n_axis_values": int(len(axes)),
            "n_feature_vectors": int(len(X)), "per_network_max_abs_diff": per_net,
            "max_abs_diff": worst, "contaminated": bool(worst > 0.0)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-line", action="append", type=Path, default=None,
                    help="an Amarsi rya817_3dnlte_per_line.csv from a PUBLISHED run; "
                         "repeatable. Defaults to the committed RYA-817 artifact.")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data" / "results" / "rya1098"
                    / "rya1098_amarsi_model_integrity.json")
    a = ap.parse_args()

    rep = validate()
    print(f"\n=== RYA-1098 — Amarsi 3D-NLTE model load ===")
    print(f"  written by sklearn : {rep.written_sklearn}   (read from the PICKLE BYTES; "
          f"the loaded object reports the runtime version and is not asked)")
    print(f"  loaded under       : {rep.runtime_sklearn}"
          f"{'   <-- SKEW' if rep.skewed else ''}")
    print(f"  runtime validated  : {rep.validated_runtime}  "
          f"(VALIDATED_RUNTIMES={sorted(VALIDATED_RUNTIMES)})")
    print(f"\n  CHECK 1  authors' published pair : predicted {rep.reference_predicted:+.6f} "
          f"vs published {rep.reference_published:+.3f}  -> "
          f"{'PASS' if rep.reference_ok else 'FAIL'}")
    print(f"  CHECK 2  numpy vs predict()      : max|diff| {rep.forward_max_abs_diff:.3e} "
          f" -> {'PASS' if rep.forward_ok else 'FAIL'}")

    csvs = a.per_line or [ROOT / "data" / "results" / "rya817"
                          / "rya817_3dnlte_per_line.csv"]
    lines = feed_vectors([p for p in csvs if p.exists()])
    cc = contamination_check(lines) if len(lines) else {"n_lines": 0}
    if cc.get("n_lines"):
        print(f"  CHECK 3  the feed's own inputs   : {cc['n_feature_vectors']} vectors "
              f"({cc['n_lines']} lines x {cc['n_axis_values']} axis values), "
              f"max|diff| {cc['max_abs_diff']:.3e}  -> "
              f"{'CONTAMINATED' if cc['contaminated'] else 'NOT contaminated'}")
    else:
        print("  CHECK 3  the feed's own inputs   : SKIPPED -- no per-line artifact given")

    out = {"ticket": "RYA-1098",
           "written_sklearn": rep.written_sklearn, "runtime_sklearn": rep.runtime_sklearn,
           "skewed": rep.skewed, "runtime_in_validated_set": rep.validated_runtime,
           "model_files": list(MODEL_FILES),
           "written_version_per_file": {f: written_version(VENDOR_DIR / f)
                                        for f in MODEL_FILES},
           "check1_authors_published_pair": {
               "predicted": rep.reference_predicted,
               "published": rep.reference_published, "pass": rep.reference_ok},
           "check2_independent_forward_pass": {
               "max_abs_diff": rep.forward_max_abs_diff, "pass": rep.forward_ok},
           "check3_feed_inputs": cc,
           "verdict": ("the 1.0.2 -> 1.9.0 skew is BENIGN for these estimators; the "
                       "corrections already in the feed are NOT contaminated"
                       if rep.ok and not cc.get("contaminated", True)
                       else "NOT PROVEN BENIGN -- see the failing check")}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\n  VERDICT: {out['verdict']}")
    try:
        shown = a.out.relative_to(ROOT)
    except ValueError:
        shown = a.out          # an out-of-repo scratch path is legitimate for a dry look
    print(f"  -> {shown}")
    return 0 if rep.ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

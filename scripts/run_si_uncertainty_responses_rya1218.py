#!/usr/bin/env python3
"""Run the Si stellar-parameter response campaign on Sirius.

This harness is intentionally explicit about the execution host: the bundled MOOG
binary is Linux x86-64 and must be run on Sirius.  It records the aggregate response
diagnostic, while the RYA-587 contract still requires exact per-line responses before
admitting a stellar component.
"""
from __future__ import annotations

import contextlib
import io
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import abundances_derive as ad


BASE = {"teff_K": 5772.0, "logg": 4.438, "vturb_kms": 1.0, "feh": 0.0}
CAMPAIGN = (("teff_K", 100.0), ("vturb_kms", 0.10))


def run(output: Path) -> None:
    runs = []
    for parameter, step in CAMPAIGN:
        for side in ("minus", "plus"):
            params = dict(BASE)
            params[parameter] += -step if side == "minus" else step
            started = time.time()
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                _, frame = ad.run(star_id="solar", engine="spectrum", skip_convergence=True,
                                  stellar_params_override=params,
                                  uncertainty_per_line_out=str(output.parent / f"solar_per_line_{parameter}_{side}.csv"))
            si = frame[frame.element.astype(str).str.lower().eq("si")]
            runs.append({"parameter": parameter, "side": side, "step": step,
                         "params": params, "elapsed_s": time.time() - started,
                         "status": "OK", "rows": si.to_dict(orient="records")})
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"ticket": "RYA-1218", "contract": "RYA-587",
                                  "execution_host": "sirius", "nominal": BASE,
                                  "runs": runs}, indent=2) + "\n")


if __name__ == "__main__":
    run(Path("data/results/rya1218/uncertainty_runs_20260916/stellar_responses.json"))

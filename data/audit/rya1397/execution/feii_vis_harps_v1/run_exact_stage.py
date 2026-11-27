import os
import runpy
import sys
from pathlib import Path
from pipeline.line_patch import LinePatch, LinePatchError
from pipeline.run_descriptor import RunDescriptor, resolve
from pipeline.run_matrix import verify_numpy_ceiling

out = Path("rya1397-benchmark-output/exact-harps-reference-v1").resolve()
out.mkdir(parents=True, exist_ok=False)
prior = Path("data/results/orchestrator/budget/solar_Fe/stage").resolve()
stem = "FeII_4200_6908_harps_solar_harps_molecfit_corrected_SYNTH_REFERENCE"
expected = {4233.162, 4385.377, 4416.819, 4508.280, 4522.628, 4549.466, 4583.829}
interpreter = os.environ["CODEX_SYNTH_PYTHON"]
ok, reason = verify_numpy_ceiling(interpreter)
if not ok:
    raise SystemExit(reason)
descriptor = RunDescriptor("Fe", "II", "harps", "solar_harps_molecfit_corrected",
                           4200.0, 6908.0, "ts-lte", "synthesis", "reference")
plan = resolve(descriptor, interpreter=interpreter, ispec_dir=os.environ["ISPEC_DIR"])
if plan.blocked_reason:
    raise SystemExit(plan.blocked_reason)
step, = plan.steps
bind = LinePatch.bind
get_prior = LinePatch.prior

def strict_bind(self, *args, **kwargs):
    bind(self, *args, **kwargs)
    if self.stem != stem or self.not_in_pool or self.refit_wl != expected:
        raise LinePatchError(f"Benchmark binding differs: {self.stem}, {sorted(self.refit_wl)}, {self.not_in_pool}")

def strict_prior(self, treatment, wave, ep):
    if not (prior / f"{stem}_{treatment}_lines.csv").is_file():
        raise LinePatchError(f"No pinned prior fit table for {treatment}")
    row = get_prior(self, treatment, wave, ep)
    if row is None and round(float(wave), 3) not in expected:
        raise LinePatchError(f"No reusable fit for untouched line {wave:.4f} in {treatment}")
    return row

LinePatch.bind = strict_bind
LinePatch.prior = strict_prior
script = step["script"]
args = step["args"]
args[args.index("--out") + 1] = str(out)
sys.path.insert(0, str(Path(script).parent.resolve()))
sys.argv = [script, *args, "--refit-lines",
            "data/audit/rya1397/patch_inputs/solar_FeII_VIS_harps_reference.csv",
            "--reuse-fits-from", str(prior)]
print("Resolved and gated stage:", sys.argv, flush=True)
runpy.run_path(script, run_name="__main__")

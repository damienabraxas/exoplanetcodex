"""RYA-1178 A — per-run xi provenance, and a pairing guard that reads it.

🔴 WHY THIS EXISTS. `dA/dxi` is a CENTRAL DIFFERENCE over two synthesis runs, one at
`xi_nominal - step` and one at `xi_nominal + step`. Nothing in either run's output said
which xi produced it. The pairing was true only because each leg was written into its own
directory and the caller remembered which was which -- WORKTREE ISOLATION AS A PROOF,
which is not a proof at all. The near-UV worktree collision is the demonstration: two runs
wrote the same stem, the second silently overwrote the first, and the only reason those
numbers are known-good today is that somebody re-ran them. No downstream check would have
caught a mis-pairing, because no artifact carried the fact that could contradict it.

So every perturbation run stamps its own xi and its own leg into its output directory, and
the derivative REFUSES to compute unless the two stamps it is handed form the expected
`nominal -/+ step` pair. A guard that cannot fail is not a guard (RYA-853), so the refusal
is exercised by a test that hands it a deliberately mis-paired run.

⚠️ THIS IS THE DERIVATIVE-RUN LEVEL, not the product level. `xi_value` on a published
product (RYA-1178 Part 1) records the xi the PUBLISHED abundance was measured at -- the
pinned 1.0. This records the xi a PERTURBATION LEG was run at, which is deliberately not
1.0. The two must never be conflated: a product whose xi_value is 1.0 is correct; a
derivative leg whose stamped xi is 1.0 means the perturbation did not bite.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

STAMP_NAME = "xi_run.json"

#: The two legs of a central difference. Stored as a word, not a sign, so a stamp read
#: by eye cannot be misread and a flipped sign cannot masquerade as a valid pair.
LEG_MINUS, LEG_PLUS = "minus", "plus"
LEGS = (LEG_MINUS, LEG_PLUS)

#: xi is carried in km/s at 4 dp in the stamp; two legs must match their nominal to this.
XI_TOL_KMS = 1e-9


class XiPairingError(RuntimeError):
    """A derivative was asked for from runs that do not provably form a xi pair."""


def stamp_path(outdir: Path | str) -> Path:
    return Path(outdir) / STAMP_NAME


def write_stamp(outdir: Path | str, *, xi_kms: float, leg: str,
                xi_nominal: float, step_kms: float, unit: str = "",
                extra: dict | None = None) -> Path:
    """Record, IN THE RUN'S OWN OUTPUT DIRECTORY, the xi that produced it.

    Written by the leg itself, beside the products it emitted, so the fact travels with
    the artifact rather than living in the caller's memory.
    """
    if leg not in LEGS:
        raise XiPairingError(f"leg must be one of {LEGS}, got {leg!r}")
    expect = xi_nominal - step_kms if leg == LEG_MINUS else xi_nominal + step_kms
    if not math.isclose(xi_kms, expect, abs_tol=XI_TOL_KMS):
        raise XiPairingError(
            f"refusing to stamp leg {leg!r} with xi={xi_kms!r}: for nominal "
            f"{xi_nominal} and step {step_kms} that leg must be {expect}. A stamp that "
            f"disagrees with its own leg is worse than no stamp -- it would make a "
            f"mis-paired run look verified.")
    p = stamp_path(outdir)
    p.parent.mkdir(parents=True, exist_ok=True)
    doc = {"xi_kms": float(xi_kms), "leg": leg,
           "xi_nominal_kms": float(xi_nominal), "step_kms": float(step_kms),
           "unit": unit, "ticket": "RYA-1178"}
    if extra:
        doc["extra"] = extra
    p.write_text(json.dumps(doc, indent=2) + "\n")
    return p


def read_stamp(outdir: Path | str) -> dict:
    p = stamp_path(outdir)
    if not p.exists():
        raise XiPairingError(
            f"no {STAMP_NAME} in {outdir} -- this run does not record which xi produced "
            f"it, so its leg cannot be established. Refusing to guess from the directory "
            f"name: the directory is exactly the isolation this guard exists to stop "
            f"relying on.")
    try:
        d = json.loads(p.read_text())
    except json.JSONDecodeError as exc:
        raise XiPairingError(f"{p} is not readable JSON: {exc}") from exc
    for k in ("xi_kms", "leg", "xi_nominal_kms", "step_kms"):
        if k not in d:
            raise XiPairingError(f"{p} is missing required key {k!r}")
    if d["leg"] not in LEGS:
        raise XiPairingError(f"{p} carries leg {d['leg']!r}, not one of {LEGS}")
    return d


def assert_pair(minus_dir: Path | str, plus_dir: Path | str, *,
                xi_nominal: float | None = None,
                step_kms: float | None = None) -> tuple[dict, dict]:
    """Refuse unless these two runs provably are the -/+ legs of ONE central difference.

    Every failure mode below has to be its own refusal, because each one produces a
    NUMBER rather than an error when it is not checked -- a wrong derivative is
    indistinguishable from a right one downstream.
    """
    a, b = read_stamp(minus_dir), read_stamp(plus_dir)

    if a["leg"] == b["leg"]:
        raise XiPairingError(
            f"both runs are stamped leg {a['leg']!r}. Two copies of the same leg "
            f"difference to ~0 and would publish 'xi does not matter' with total "
            f"confidence -- the RYA-853 vacuous-guard shape.")
    if a["leg"] != LEG_MINUS or b["leg"] != LEG_PLUS:
        raise XiPairingError(
            f"legs are the wrong way round: minus_dir is stamped {a['leg']!r} and "
            f"plus_dir {b['leg']!r}. Swapped legs flip the SIGN of dA/dxi silently.")

    for name, k in (("nominal xi", "xi_nominal_kms"), ("step", "step_kms")):
        if not math.isclose(a[k], b[k], abs_tol=XI_TOL_KMS):
            raise XiPairingError(
                f"the two legs disagree on {name}: {a[k]} vs {b[k]}. They are not two "
                f"halves of one central difference.")

    nominal = a["xi_nominal_kms"] if xi_nominal is None else float(xi_nominal)
    step = a["step_kms"] if step_kms is None else float(step_kms)
    if xi_nominal is not None and not math.isclose(a["xi_nominal_kms"], nominal, abs_tol=XI_TOL_KMS):
        raise XiPairingError(
            f"stamps say nominal xi = {a['xi_nominal_kms']}, caller expected {nominal}")
    if step_kms is not None and not math.isclose(a["step_kms"], step, abs_tol=XI_TOL_KMS):
        raise XiPairingError(
            f"stamps say step = {a['step_kms']}, caller expected {step}")

    if step <= 0:
        raise XiPairingError(
            f"step is {step}: a zero or negative perturbation cannot produce a "
            f"derivative, it produces a division by ~0.")

    want_lo, want_hi = nominal - step, nominal + step
    if not math.isclose(a["xi_kms"], want_lo, abs_tol=XI_TOL_KMS):
        raise XiPairingError(
            f"minus leg was run at xi={a['xi_kms']}, but nominal {nominal} - step {step} "
            f"is {want_lo}. Refusing: the span the derivative divides by would not be "
            f"the span the runs actually spanned.")
    if not math.isclose(b["xi_kms"], want_hi, abs_tol=XI_TOL_KMS):
        raise XiPairingError(
            f"plus leg was run at xi={b['xi_kms']}, but nominal {nominal} + step {step} "
            f"is {want_hi}. Refusing for the same reason.")
    if math.isclose(a["xi_kms"], b["xi_kms"], abs_tol=XI_TOL_KMS):
        raise XiPairingError(
            f"both legs were run at xi={a['xi_kms']} -- there is no perturbation to "
            f"differentiate.")
    return a, b


def span_kms(minus_stamp: dict, plus_stamp: dict) -> float:
    """The central-difference span the two runs ACTUALLY spanned, from the stamps.

    Taken from the stamps rather than from the caller's `step` so the number the
    derivative divides by is the number the runs were made at.
    """
    return float(plus_stamp["xi_kms"]) - float(minus_stamp["xi_kms"])

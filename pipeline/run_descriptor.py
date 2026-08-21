"""One object that describes a RUN, and a resolver that says whether it can happen.

RYA-767 Seam 2. Owns: the identity of a run, and the deterministic resolution of
that identity against what data and tooling ACTUALLY exist. Does NOT own: any
science. It runs nothing, measures nothing, and decides no physics -- it answers
"is this run possible, in what order, and under what preconditions", and every
answer it gives is derived from a registry, never from a literal typed here.

WHY THIS EXISTS
---------------
Nothing in the repo represented "a run" as an object. The decisions were spread
across `engine_selection`, `band_policy`, the holdings registry, the drivers'
argparse, and -- mostly -- the operator's head. RYA-933/934 made that concrete:
driving 24 Fe runs by hand needed six pieces of knowledge that live in NO code
path, and getting any of them wrong produced a plausible wrong answer rather
than an error:

  1. which interpreter (numpy must be < 2.3, RYA-682, or iSpec writes a
     zero-row artifact and exits 0);
  2. ISPEC_DIR, or the synthesis resolves its line list beside the repo and dies
     on a file that was never there;
  3. which holdings a band is actually covered by -- two telluric-corrected
     Kitt Peak products were registered and reachable by no loader (RYA-904);
  4. that the EW step must precede the products step, and that the products step
     reads the EW table by a filename convention;
  5. that the filename convention keyed on INSTRUMENT, so two holdings of one
     instrument silently overwrote each other;
  6. that near-UV forbids profile-fit entirely (`band_policy`), so the same
     request means a different method there.

Every one of those is now a PRECONDITION on the descriptor rather than a thing
the operator has to remember. That is the whole point: the judgment moves out of
the operator and into data, so the executor can be something that makes no
judgments at all -- a shell loop, a CI job, or a local model on Sirius.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Literal

from pipeline import band_policy

#: RYA-682. Above this, `ispec/abundances.py:132` fails and synthesis writes a
#: ZERO-ROW artifact while exiting 0 -- a control run reports success on no data.
#: The ceiling is a property of the ENGINE, so it is declared with the engine and
#: checked before dispatch, not discovered when a product comes back empty.
NUMPY_CEILING = (2, 3)

Method = Literal["profile-fit", "synthesis", "interval-integration"]


class RunNotPossible(RuntimeError):
    """A run that cannot happen, carrying WHY.

    Distinct from a failure: this is answered BEFORE anything executes, and the
    reason names which precondition was not met. "We do not hold this" and "we
    hold it in a state we may not measure" are different answers and must not
    collapse into one (RYA-796 / RYA-833).
    """


@dataclass(frozen=True)
class RunDescriptor:
    """The identity of one run. Carries state; decides nothing.

    Ryan's Salesforce analogy taken literally (RYA-767): one trigger that holds
    state, many stateless handlers. Everything downstream is a pure function of
    this object, which is what makes a run serialisable, resumable, diffable --
    and dispatchable by something other than a person.
    """
    element: str
    ion: str
    instrument: str
    holding: str
    lo_A: float
    hi_A: float
    engine_deck: str = "ts-lte"

    @property
    def band(self) -> str:
        return band_policy.resolve(0.5 * (self.lo_A + self.hi_A)).name

    @property
    def species(self) -> str:
        return f"{self.element}{self.ion}"

    @property
    def key(self) -> str:
        """Stable identity, and the ONLY thing that should name an artifact.

        The holding is in here deliberately. RYA-933/934: the stem keyed on
        instrument, so `solar_harps` and `solar_harps_molecfit_corrected` wrote
        the same file and the second overwrote the first -- two products
        differing precisely by whether tellurics were removed, which is the one
        pair this project must never collapse.
        """
        return (f"{self.species}_{int(self.lo_A)}_{int(self.hi_A)}"
                f"_{self.instrument}_{self.holding}")

    def as_dict(self) -> dict:
        d = asdict(self)
        d.update(band=self.band, key=self.key)
        return d


@dataclass
class Precondition:
    """One thing that must be true, and how to find out. Never assumed."""
    name: str
    satisfied: bool
    detail: str


@dataclass
class ResolvedRun:
    """A descriptor plus the answer: can it run, in what steps, under what checks."""
    descriptor: RunDescriptor
    method: Method
    preconditions: list[Precondition] = field(default_factory=list)
    steps: list[dict] = field(default_factory=list)
    blocked_reason: str | None = None

    @property
    def runnable(self) -> bool:
        return self.blocked_reason is None and all(p.satisfied for p in self.preconditions)

    def as_dict(self) -> dict:
        return {"descriptor": self.descriptor.as_dict(), "method": self.method,
                "runnable": self.runnable, "blocked_reason": self.blocked_reason,
                "preconditions": [asdict(p) for p in self.preconditions],
                "steps": self.steps}


def method_for(descriptor: RunDescriptor) -> Method:
    """Which measurement method this band PERMITS -- from policy, never from taste.

    near-UV permits synthesis only; asking for a profile fit there is not a
    preference the caller gets to express.
    """
    policy = band_policy.resolve(0.5 * (descriptor.lo_A + descriptor.hi_A))
    if "profile-fit" in policy.permitted_methods:
        return "profile-fit"
    if "synthesis" in policy.permitted_methods:
        return "synthesis"
    raise RunNotPossible(
        f"band {policy.name} permits none of the methods this layer can dispatch "
        f"(permitted: {policy.permitted_methods})")


def _holding_spec(descriptor: RunDescriptor):
    """The HoldingSpec this run names, from the harness's own table.

    Read from `measure_band_ew` rather than restated here: a second list of
    holdings would be a second declaration of one fact, which is how RYA-845's
    double-count survived and what `loader_coverage` already refuses to repeat.
    """
    import sys
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    if str(root / "scripts") not in sys.path:
        sys.path.insert(0, str(root / "scripts"))
    from measure_band_ew import _INSTRUMENT_HOLDINGS
    for spec in _INSTRUMENT_HOLDINGS.get(descriptor.instrument, ()):
        if spec.holding_id == descriptor.holding:
            return spec
    raise RunNotPossible(
        f"{descriptor.holding!r} is not wired for instrument {descriptor.instrument!r}. "
        f"It may be registered and still unreachable -- that is the RYA-904 shape, and "
        f"an unreachable holding reads to every caller exactly like having no data. "
        f"Wired: {[s.holding_id for s in _INSTRUMENT_HOLDINGS.get(descriptor.instrument, ())]}")


def resolve(descriptor: RunDescriptor, *, interpreter: str | None = None,
            ispec_dir: str | None = None) -> ResolvedRun:
    """Answer whether this run can happen, and emit its ordered steps.

    Every precondition here is one an operator had to hold in their head during
    RYA-933/934. None of them is checked by the drivers themselves, which is why
    each one failed silently or late at least once.
    """
    checks: list[Precondition] = []
    blocked: str | None = None

    try:
        method = method_for(descriptor)
    except RunNotPossible as exc:
        return ResolvedRun(descriptor, "synthesis", blocked_reason=str(exc))

    # ── 1. is the holding reachable at all, and does it cover this band ──────
    try:
        spec = _holding_spec(descriptor)
        covers = spec.covers(0.5 * (descriptor.lo_A + descriptor.hi_A),
                             0.5 * (descriptor.hi_A - descriptor.lo_A))
        checks.append(Precondition(
            "holding_wired", True, f"{descriptor.holding} -> reader {spec.reader!r}"))
        checks.append(Precondition(
            "holding_covers_band", covers,
            f"span {spec.span_A} vs requested {descriptor.lo_A}-{descriptor.hi_A} A"))
        if not covers:
            blocked = (f"{descriptor.holding} declares span {spec.span_A} and does not "
                       f"cover {descriptor.lo_A}-{descriptor.hi_A} A. NOT a failure -- "
                       f"we do not hold this window in this product.")
    except RunNotPossible as exc:
        checks.append(Precondition("holding_wired", False, str(exc)))
        blocked = str(exc)
        spec = None

    # ── 2. the telluric gate is consulted, never inferred from the instrument ─
    try:
        from pipeline.telluric_policy import gate_holding
        may_run, reason = gate_holding(descriptor.holding, descriptor.instrument)
        checks.append(Precondition("telluric_gate", bool(may_run), reason[:200]))
        if not may_run and blocked is None:
            blocked = f"telluric gate refuses {descriptor.holding}: {reason[:200]}"
    except Exception as exc:                                   # noqa: BLE001
        checks.append(Precondition("telluric_gate", False, f"gate could not answer: {exc}"))

    # ── 3. the continuum contract is a property of the HOLDING (RYA-904/713) ─
    if spec is not None:
        checks.append(Precondition(
            "continuum_contract", True,
            f"pre_normalised={spec.pre_normalised} -- "
            + ("the product ships its own continuum and the harness consumes it"
               if spec.pre_normalised else
               "the product ships NO continuum, so the harness must place one")))

    # ── 4. the engine's numpy ceiling, checked BEFORE dispatch (RYA-682) ─────
    if interpreter:
        checks.append(Precondition(
            "numpy_below_ceiling", True,
            f"caller pinned {interpreter}; the executor MUST verify "
            f"numpy < {'.'.join(map(str, NUMPY_CEILING))} before running -- above it "
            f"iSpec writes a zero-row artifact and exits 0"))
    else:
        checks.append(Precondition(
            "numpy_below_ceiling", False,
            f"no interpreter pinned. Synthesis needs numpy < "
            f"{'.'.join(map(str, NUMPY_CEILING))} (RYA-682); the default interpreter is "
            f"NOT guaranteed to satisfy it and the failure mode is a silent empty product."))
        if blocked is None:
            blocked = "no interpreter pinned and the engine has a hard numpy ceiling"

    # ── 5. the engine's environment ─────────────────────────────────────────
    checks.append(Precondition(
        "ispec_dir_set", bool(ispec_dir),
        f"ISPEC_DIR={ispec_dir}" if ispec_dir else
        "ISPEC_DIR unset -- synthesis resolves its line list RELATIVE TO THE REPO and "
        "dies on a path that was never staged there"))
    if not ispec_dir and blocked is None:
        blocked = "ISPEC_DIR unset; the synthesis step cannot find its line list"

    # ── 6. the ordered steps, each with the postcondition that proves it ─────
    env = {"ISPEC_DIR": ispec_dir} if ispec_dir else {}
    common = ["--element", descriptor.element, "--ion", descriptor.ion,
              "--lo", f"{descriptor.lo_A:g}", "--hi", f"{descriptor.hi_A:g}",
              "--instrument", descriptor.instrument, "--holding", descriptor.holding]
    steps: list[dict] = []
    if method == "profile-fit":
        steps.append({
            "name": "measure_ew", "script": "scripts/measure_band_profilefit.py",
            "args": common, "env": env, "interpreter": interpreter,
            "produces": f"data/measured/band_ew/{descriptor.key}_PROFILEFIT_ew.csv",
            "postcondition": "the EW table exists and has at least one row",
            "why_ordered_first": ("derive_band_products READS this table by filename; "
                                  "it does not measure. Reversing the order fails with "
                                  "'no measured EWs', which reads like missing data "
                                  "rather than a missing step."),
        })
    steps.append({
        "name": "derive_products", "script": "scripts/derive_band_products.py",
        "args": common + ["--engine-b-deck", descriptor.engine_deck],
        "env": env, "interpreter": interpreter,
        "produces": f"{descriptor.key}_"
                    f"{'PROFILEFIT' if method == 'profile-fit' else 'SYNTH'}_products.csv",
        "postcondition": "the products table exists and every treatment carries a value "
                         "and an ErrorBudget; a zero-row product is a FAILURE, not a null",
    })
    return ResolvedRun(descriptor, method, checks, steps, blocked)

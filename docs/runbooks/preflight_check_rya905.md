# Pre-flight readiness check — RYA-905

`scripts/preflight_check.py` is the runnable form of the RYA-716 Al pre-flight comment:
the "is this element ready to measure" checklist, made executable and reusable for every
element and every star.

```bash
python3 scripts/preflight_check.py --star sun --element Fe
python3 scripts/preflight_check.py --star sun --element Al --ion I --json report.json
```

It **never blocks**. Exit status is 0 whatever it finds.

## Why it exists

The Fe/HARPS night (RYA-896/897/898) cost a working day to a defect that never touched an
abundance value: the direct-solar HARPS arm was named in the harness's own config table
and the `load_window` dispatch had no branch for it. The holding was verified. Nothing
raised. It surfaced only when a human eyeballed the rendered page and counted instruments
by hand.

The failure was **silence**, not error. This turns that silence into a line of output.

## The two severities

| | meaning | example |
|---|---|---|
| `INFO` | **expected absence** — a holding/line/grid we simply do not have | no UVES for the Sun |
| `WARN` | **silent gap** — something we DO have that the pipeline cannot see | verified HARPS, no loader |
| `ERROR` | the check itself is broken | the dispatch reader cannot parse the harness |

The INFO/WARN distinction is the whole point: *"we're missing HARPS data"* (we weren't)
and *"we HAVE verified HARPS and the code can't reach it"* (the bug) looked identical from
the rendered page. Every WARN carries a **discriminator** — the sentence saying why it is
not an expected absence — and a one-line **suggested-ticket stub**.

An `ERROR` means the report's absences are not evidence of anything. Read nothing into a
report that carries one.

## The six checks

1. **Instrument reachability** — every verified holding for this star, against the band
   harness's `load_window` dispatch. Held+verified but not loadable is the RYA-897 class.
   The dispatch is *read* from `scripts/measure_band_ew.py`'s AST, never imported: the
   harness resolves the Kitt Peak atlas at import and exits when it is absent, so
   importing it would make the readiness check impossible exactly where it is most needed.
2. **Line coverage** — the element's best lines against what was actually measured. A line
   inside a band that was run, absent from that run's output, was silently dropped.
3. **NLTE grid reach** — over the **best** lines, not merely some lines. The RYA-773 Al
   gap is the shape: the grid served 6696/6698 and not the headline 7835/8772.
4. **Anchor consistency** — the solar anchor's declared instrument chain against the chain
   the products were measured on. Disjoint chains do not cancel in `[X/H]` (RYA-898).
5. **Rendered-output reconciliation** — every verified holding and every band, counted
   against what the product actually contains. The Fe-page class.
6. **Telluric state** — `telluric_policy.gate_holding`, called, not reimplemented.
   Uncorrected and withheld is expected; corrected and dropped anyway is the CRIRES+ class.

**One gap yields one WARN.** Checks 5 and 6 defer to the check that owns an absence and
report it as INFO with a pointer, so the suggested-ticket list is not three tickets for
one fix.

## The controls

Two checks assert negatives, and a negative needs a control at both ends:

* **positive** — the dispatch reader must find `kpno_solar_atlas`, the arm every committed
  band product was measured on. A reader that cannot see the arm we know is wired would
  emit a page of false WARNs; the controls make that state `ERROR` instead.
* **negative** — the same reader must not report a sentinel name that cannot exist.

`tests/test_preflight_check_rya905.py` runs the reader against a synthetic harness *with*
a `harps` branch and asserts the WARN disappears. Without that, "it warns about HARPS" is
compatible with "it warns about everything".

## Scope of an absence

Every absence is scoped to **the artifacts inspected**, never to the project. The default
root is this checkout alone. Widen it explicitly:

```bash
python3 scripts/preflight_check.py --star sun --element Fe \
    --artifact-root /mnt/codex-data/codex/rya845
```

The report prints the roots it looked in, so "no measured-EW artifact for Fe I" always
reads as a statement about those roots.

## Calling it from the element sweep (RYA-872)

The sweep runs this **before** measuring each element and treats it as advisory:

```python
from scripts import preflight_check as pf

state, results = pf.run(star, element, ion, [pf.ROOT])
print(pf.render(state, results))
# proceed with the measurement regardless — a survey legitimately lacks data
```

Or via the CLI with `--json`, whose payload carries every finding with its severity,
discriminator and ticket stub.

Do not turn this into a gate. A survey lacking an arm is the normal state, and a gate that
fires on the normal state is a gate somebody switches off.

## What it does not do

* It **writes nothing** except its own report, and only when `--json` asks.
* It **re-verifies nothing**. The intake framework (RYA-301/265, RYA-806, RYA-376) already
  produced the ground truth; this reconciles those registries against what the harness and
  the grids can actually reach.
* It reads every registry from the one place that owns it — never a copy, and no value is
  restated in the module.

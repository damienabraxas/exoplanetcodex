# The Element Protocol — what "we have done this element" means

**RYA-709.** Ryan, 2026-08-09: *"we create 27 tickets, for each element, and we go through
each one, and verify, rerun, give it the Al treatment. If it fails, why, and document."*

This is that protocol. It exists because "the Al treatment" is ambiguous, and the
ambiguity is dangerous: aluminium was worked hard for a session and reached **4 of its 26
reachable lines**. Calling that done, and stamping it on twenty-six more elements, would
industrialise a 15% pass and call it a sweep.

---

## The definition of done

An element is **DONE** when every one of its reachable usable lines has an outcome. Not a
value — an *outcome*. A line that cannot be measured is done when the reason is recorded
and checkable.

```
usable          predicted central depth 0.05-0.60
reachable       inside an instrument we hold, per pipeline.coverage
DONE            every reachable usable line has an outcome
```

`scripts/line_accounting_rya709.py` is the scoreboard. An element is done when its
`unmeasured_reachable` is zero — every remaining line having an explicit, recorded
disposition rather than silence.

## The eight steps

Each is a gate. A step that cannot be completed is recorded as a stop, with its reason,
and the element carries that reason into the SPP appendix.

1. **ACCOUNT.** Run the line accounting. Know the four numbers before touching anything:
   usable, reachable, measured, unmeasured. Never start from the pool.
2. **COVER.** For every reachable line, which instruments see it. Never ask a loaded
   array — ask `pipeline.coverage`. A wavelength no instrument reaches is an acquisition
   target and leaves the element's scope.
3. **MEASURE.** Local continuum, window half-width from the **line separation** (not the
   FWHM — pair members must never share flux). Measure on **every** covering instrument,
   not the first one.
4. **GRADE.** A NIST ASD pull for every line that will carry a value. Unresolved
   components at one wavelength are **summed**, exactly as HFS is. Ungraded gf HOLDS the
   promotion — it does not block the measurement.
5. **REACH.** If the abundance path cannot see the line, author the line region from the
   graded pull: `loggf` / `Ei` / `Ek` / `J` sourced, `nlte` set honestly, fit columns
   **zeroed**, every inherited constant named as inherited.
6. **LADDER.** EW/1D-LTE, then NLTE, then synthesis — **executed, not narrated**. Each
   rung records its own outcome. Two states get named because they hide: *stopped early*
   and *escalated without cause*.
7. **CROSS-CHECK.** Where two instruments cover a line, measure both. Agreement is
   corroboration a single arm cannot give; disagreement concentrated on blended lines is
   the system working, not noise.
8. **REPORT.** Per (instrument × band), never one collapsed number. Every unresolved line
   defends its blank in SPP Appendix A with plots and measured evidence.

## The three honest outcomes

| outcome | meaning |
|---|---|
| **RESOLVED** | a value, per instrument and band, with its gf grades and its caveats |
| **BOUNDED** | an upper limit or a held value, with the bound's basis |
| **UNRESOLVED** | every applicable rung executed and failed, each failure recorded |

A fourth state is a **defect, not an outcome**: unresolved where a rung was never
attempted. That is a finding against the pipeline, not against the star.

## Aluminium, scored honestly against this

The worked example, and it is a **partial**:

| step | state |
|---|---|
| 1 ACCOUNT | done — 55 usable, 26 reachable, 26 unmeasured |
| 2 COVER | done — HARPS / Kitt Peak / IAG registered and verified |
| 3 MEASURE | **4 of 26** lines, on Kitt Peak; IAG depths only |
| 4 GRADE | done for those 4 (B/B+) and the optical pair (C+) |
| 5 REACH | done — four line regions authored |
| 6 LADDER | rung 1 done; **rung 2 blocked, no NLTE grid**; rung 3 not run |
| 7 CROSS-CHECK | done for 2 optical lines; **not done for the 4 IR** |
| 8 REPORT | Appendix A prototype published |

**A(Al) = 6.415 ± 0.037, 1D-LTE, Kitt Peak, NIR band.** Not promoted, nothing in the pool.
**Al is IN PROGRESS, not done** — 22 reachable lines remain untouched.

## What this protocol cost, measured

Aluminium's four lines took roughly one working session and produced **thirteen defects in
one adapted harness**, a false NO-DATA claim that reached the state register, a duplicated
instrument catalog, and six broken consumers from a CSV header. That is the honest rate.
Budget for it rather than being surprised by it.

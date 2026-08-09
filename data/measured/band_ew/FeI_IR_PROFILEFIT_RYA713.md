# Fe I in the IR, re-measured by profile fit — RYA-713

Ryan: *"rerun the IR band with the profile fitter"*

Run under the method that **passed** the optical control at −0.013 dex, with the corrected
continuum policy (Kitt Peak ships residual flux; unity is the continuum).

| | interval integration | **profile fit** |
|---|---|---|
| measured | 445 | **447** |
| in aggregate | 271 | **103** |
| `corr(strength, log EW)` | 0.740 | **0.761** |
| EW in aggregate | — | 10.1 – 108.5 mÅ |
| REW | — | −5.94 … −4.90, all unsaturated |

## The yield dropped and that is the point

271 → 103. The interval method admitted 168 lines that the profile fit rejects, and it did
so because integration always returns *a number* — it cannot fail. A fit can, and 185 lines
came back `FIT-PINNED`: the optimiser could not find a profile of plausible width, which
means the feature is not a single clean line. Those 185 previously entered the aggregate
carrying an integral over whatever was in the window.

| quarantine cause | n |
|---|---|
| FIT-PINNED | 185 |
| BLEND-DOMINATED | 105 |
| saturated (REW > −4.9) | 32 |
| GF-GHOST-ABSENT | 18 |
| GF-GHOST | 4 |

## The physical sanity check

`corr(intrinsic strength, log EW)` = **+0.761**, against +0.399 for the raw interval
method. For lines of one element in one ionisation stage this correlation must be strongly
positive — it is the most basic statement that what we measured behaves like the lines we
think they are. The IR now satisfies it.

It is *not* proof of accuracy. It is proof of internal consistency, which the interval
method never had.

## A real frontier cross-check exists and has not been run

**All 103 in-aggregate lines are also covered by the IAG FTS atlas** (4050–10650 Å). That is
the frontier analogue of the optical control: no reference *abundance* exists out here, but
a second independent instrument does. Two instruments with different normalisation
histories agreeing on the same lines would be genuine evidence; disagreeing would localise
the problem exactly as the Al 6696-vs-6698 test did.

IAG lives on Sirius only, so this was not run here. **It is the single highest-value
outstanding check on the IR measurement.**

## Still not an abundance

The gf blocker is untouched: none of these lines carries a NIST grade, and the root-cause
split for the band remains dominated by atomic data. A correct EW through an unknown gf is
an unknown abundance. What has changed across this sequence is that the *measurement* is no
longer the weak link — the atomic data is, unambiguously.

## Data-hygiene defect found and fixed

A `--max-lines 2` smoke test **overwrote the documented 445-line interval result**, and the
clobbered 2-row file was committed in 85b025d. Recovered from b46788f. Both drivers now
append `_SUBSET{n}` to the output stem whenever `--max-lines` is set, so a partial run can
never be mistaken for — or written over — a complete one.

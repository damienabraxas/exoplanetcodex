# The difficulty ordering is inverted — the IR is the hard band, not the UV

Ryan: *"next mission, we did Vis a control, now onto IR and UV for Fe."*

Before running either, one measurement changed the plan. **Line-list completeness is worst
exactly where the spectrum is cleanest.**

## Measured: features the Sun shows vs lines our catalogue carries

Six windows per band on the Kitt Peak atlas. Features counted as absorption minima with
prominence ≥ 2 % of the local continuum.

| band | observed features/Å | catalogued/Å | coverage |
|---|---|---|---|
| near-UV 3000–3800 | 5.10 | 3.95 | **77 %** |
| VIS 3900–6900 | 2.00 | 1.10 | 55 % |
| red-optical 6910–10000 | 1.20 | 0.45 | **37 %** |
| NIR 10000–12500 | 0.45 | 0.10 | **22 %** |

*Caveat that must travel with this:* in crowded spectrum many catalogued lines merge into
one resolvable feature, so coverage **saturates** and the near-UV figure is a lower bound.
In sparse spectrum the ratio is meaningful — a feature we can see and do not carry is a
genuinely missing line.

## Why this inverts the expected ordering

The intuition is that the near-UV is the hard band: 5.1 features/Å, median flux 0.283 at
3100 Å (72 % of the spectrum absorbed), no true continuum anywhere. All of that is true —
**observationally**. But the near-UV has been the workhorse of laboratory and classical
solar spectroscopy for a century, so it is the *best-catalogued* region we have.

The IR is the opposite: sparse, clean, continuum within 0.3 % of unity — and **we carry
between a fifth and a third of the lines that are actually there.**

## This explains the IR result exactly

The Fe I IR root-cause split found **98 of 174 failures were "a real solar line missing
from our list dominates the window."** That is not a coincidence or a harness quirk: at
37 % coverage, roughly two of every three features have no catalogue entry, so a window
centred on one of our lines is very likely to be dominated by one of theirs.

The two findings are the same fact measured two ways.

## What it means for each band

**IR — the measurement is fine, the catalogue is not.** The profile fitter passes its
control and produces internally consistent EWs (`corr(strength, log EW)` = +0.761). The
blocker is that we cannot know whether a clean-looking line is clean when we hold a third
of the line list. Extending the IR line list is worth more than any further measurement
work there.

**near-UV — the catalogue is adequate; the physics is the obstacle.** Synthesis is the only
valid method (median gap 0.146 Å < a strong line's wings), and synthesis is *exactly* the
method that needs a complete list — which the near-UV has. The obstacles there are the
pseudo-continuum (median flux 0.283–0.805, the true continuum never observed) and the sheer
blending, both of which synthesis models rather than avoids.

**So the near-UV is more tractable than the IR right now**, which is the opposite of how I
would have ordered them an hour ago.

## Status of the Fe legs

| band | method | state |
|---|---|---|
| VIS 3924–6905 | profile fit | **control PASSED**, −0.013 dex |
| red-optical 6910–9199 | profile fit | measured, 103 in aggregate |
| red-optical 9199–10000 | profile fit | measured, 12 in aggregate of 51 (40 skipped — atlas gaps) |
| **near-UV 2960–3800** | **synthesis** | **not started** — needs a synthesis harness and its own control |
| NIR 10000+ | synthesis + telluric | not started |

**The near-UV cannot reuse the validated harness.** Profile fitting is refused there by band
policy, so the UV leg needs a synthesis path — and per the control/frontier rule that path
must first reproduce the known optical answer before it is trusted in the UV. The control
discipline does not transfer between *methods*, only the requirement does.

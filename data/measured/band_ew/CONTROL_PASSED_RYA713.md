# The optical control PASSES — profile fit reproduces the known answer

Ryan: *"wire in the RYA-429 profile fitter and rerun the control"*

| method | n | median KP/HARPS | MAD | within 10 % | dex offset |
|---|---|---|---|---|---|
| interval integration | 146 | 0.773 | 0.490 | 15 % | **−0.1121** |
| profile fit, re-normalised continuum | 74 | 0.816 | 0.116 | 20 % | −0.0884 |
| **profile fit, atlas continuum** | **47** | **0.971** | **0.060** | **55 %** | **−0.0129** |

**−0.112 dex → −0.013 dex.** Scatter improved 8× (MAD 0.490 → 0.060). Weak lines
(HARPS 10–30 mÅ) give a ratio of **1.004**; 30–80 mÅ gives 0.971.

## Three separate faults, fixed in order

**1. Interval integration is not an equivalent width.** A separation-derived window clips
the wings of crowded strong lines and over-reaches on isolated weak ones. Fixed by fitting
a Voigt/Gaussian and integrating the *model* — `pipeline/lines_fit`, the RYA-429 fitter,
**imported not copied**, so there is one EW definition in the project rather than two that
drift apart.

**2. The line width is not the instrument's.** `_fit_profile` seeded σ from the HARPS value
0.025 Å with a floor at 0.005 Å. Kitt Peak's instrumental σ *is* 0.005 Å — the floor. But an
observed line is broadened by star **and** instrument in quadrature, and at R ≈ 500 000 the
star dominates completely:

| | HARPS σ_inst | Kitt Peak σ_inst | stellar σ | observed on KP |
|---|---|---|---|---|
| 5500 Å | 0.0203 | 0.0047 | 0.0312 | **0.0315** |

Seeding from the instrument alone started the fit 4× below the truth; the optimiser walked
to the lower bound and stayed. Fe I 4995 returned **20.6 mÅ against a pool value of 138**.
Now `σ_init` = quadrature sum, `σ_min` = the instrumental σ as a hard physical floor, and a
fit landing *on* that floor is **flagged FIT-PINNED and excluded** rather than trusted — it
is reporting a width it never measured.

`_fit_profile` gained these as arguments with the HARPS values as defaults, so the
production HARPS path is byte-identical in behaviour.

**3. The continuum policy was not applied here.** I established it for the interval method
and then re-normalised anyway in the profile-fit path. Kitt Peak ships **residual flux** —
unity *is* the continuum by construction. Re-fitting a local continuum through
percentile-filtered edge strips lands below unity in crowded spectrum, because those strips
contain lines; dividing by it makes every line shallower. That was a coherent deficit, not
scatter, and it was worth **−0.088 dex** on its own.

## Yield, and why it is low

**159 of 2191 measured lines enter the aggregate.** That is 7 %, and it is not a bug:

| quarantine reason | n |
|---|---|
| BLEND-DOMINATED | 923 |
| FIT-PINNED | 756 |
| saturated (REW > −4.9) | 200 |
| GF-GHOST-ABSENT | 113 |
| GF-GHOST | 40 |

The candidate list is *every* catalogued Fe I line in the band at usable predicted depth,
including thousands that are hopelessly blended in the crowded blue. The HARPS pool by
comparison holds 444 hand-curated Fe I lines. A 7 % yield from an uncurated list is the
expected shape — and every rejection carries a named cause and an appendix panel.

**FIT-PINNED at 756 is worth its own look.** A pinned fit is diagnostic, not merely failed:
it means the optimiser could not find a profile of plausible width there, which usually
means the feature is not a single clean line. Whether these are recoverable with a
two-component fit is open.

## What this licenses, and what it does not

**Licenses:** the harness reproduces a known answer to −0.013 dex, so it has earned the
right to report a number where no reference exists. Per the control/frontier rule, that
residual is now the *measured* harness systematic and belongs in the frontier error budget
rather than being assumed zero.

**Does not license:** any abundance yet. The gf blocker is untouched — 0 of the surviving
IR lines carry a NIST grade, and 123 of 174 IR failures were atomic-data faults. A correct
EW through an unknown gf is still an unknown abundance. What changed is that the
*measurement* is no longer the weak link.

## Owed

1. Re-run the IR band with the profile fitter under the corrected continuum policy.
2. Carry −0.013 dex into the frontier budget as the harness systematic.
3. Investigate the 756 FIT-PINNED lines — two-component fits, or genuinely unmeasurable.
4. Still owed: NIST gf grades before any abundance.

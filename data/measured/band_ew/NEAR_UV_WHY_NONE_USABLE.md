# Near-UV Fe I: 901 candidates, 0 usable — the evidence — RYA-713

Ryan: *"So why were none of the lines usable? This should be documented as well in our
appendix."*

Figure: `data/plots/band_appendix/FeI_3000_3800_kpno_solar_atlas_PROFILEFIT_NEAR_UV_WHY_NONE.png`
— one panel per failure mode, including the four **widest-gap** lines, i.e. the ones most
likely to succeed.

## The ledger

| reason | n |
|---|---|
| BAND-POLICY (nearest catalogued neighbour < 0.30 Å) | **845** |
| BLEND-DOMINATED | 31 |
| saturated (REW > −4.90) | 13 |
| FIT-PINNED | 8 |
| GF-GHOST | 3 |
| GF-GHOST-ABSENT | 1 |
| **usable** | **0** |

## Why: the crowding, measured

Nearest **catalogued** neighbour for all 901 candidates:

| gap | n | |
|---|---|---|
| < 0.10 Å | **502** | 55.7 % |
| 0.10–0.20 Å | 244 | 27.1 % |
| 0.20–0.30 Å | 99 | 11.0 % |
| 0.30–0.50 Å | 50 | 5.5 % |
| > 0.50 Å | 6 | 0.7 % |

**Median gap 0.086 Å.** A solar Fe line has σ ≈ 0.03 Å, so its wings run to ~0.15 Å —
**wider than the distance to its neighbour for 83 % of these lines.** There is no interval
that contains one profile and excludes the next, which is precisely what the band policy
asserted from the band-level median of 0.146 Å.

## The 56 that passed the isolation test still failed — and this is the important part

| | |
|---|---|
| passed gap ≥ 0.30 Å | 56 |
| measured EW | **107 – 922 mÅ** |
| REW | −4.51 to −3.59 |
| **above the −4.90 saturation ceiling** | **56 of 56** |

A single solar Fe line is 10–150 mÅ. **922 mÅ is not a line.**

The panels show why. In every one, the flux is a forest: Fe I 3393.856 swings between 0.1
and 0.9 five times across 1.7 Å; Fe I 3585.388 never rises above **0.36** anywhere in its
window. Every panel carries **side-bands ABSORBED** — 0.885, 0.764, 0.884, **0.351**,
0.623, 0.913, 0.833. There is no continuum in the window to normalise against, and the
"isolated" line sits inside a blend complex that the profile fit swallows whole.

## The root cause of the escape's failure

The gap test measures distance to the nearest **catalogued** neighbour, and near-UV
catalogue completeness is **~77 %**. About a quarter of the features actually present are
absent from our list, so a 0.30 Å catalogued gap routinely contains an uncatalogued line.

**Absence of a neighbour in the catalogue is not absence in the spectrum.** This is the
identical lesson the IR root-cause split produced from the other direction, where 98 of 174
failures were windows dominated by a real solar line we do not carry.

## Verdict

**The band-wide ban was correct.** `permits_profile_fit_for_line()` was a reasonable
hypothesis, it was tested on the full candidate set, and it is falsified for the near-UV.
The mechanism is retained for bands whose crowding is genuinely marginal and is recorded in
`band_policy.py` as tested-and-negative here, so it is not re-tried as a workaround.

**Near-UV Fe is blocked on exactly one thing: a level-identified line list covering
3000–4200 Å.** Not the handler, not the control, and not Engine A — Engine A serves 203 Fe I
+ 828 Fe II lines down to 2960 Å and is ready the moment something can measure them.

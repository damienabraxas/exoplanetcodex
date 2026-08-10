# The near-UV per-line escape was tested and falsified — RYA-713

Ryan: *"Engine A for UV it is."*

Engine A does reach the near-UV (203 Fe I + 828 Fe II lines served in 2960–4200 Å,
verified live). But Engine A supplies the *correction*; something still has to measure
A(LTE), and the band policy bans profile fitting in the near-UV on its **median** 0.146 Å
line gap.

Isolation is a per-line property, so `permits_profile_fit_for_line()` was added to let an
individual line escape a band-level ban when its own neighbour is far enough away
(`PROFILE_FIT_MIN_GAP_A = 0.30`). Reasonable hypothesis. **It was tested on all 901
near-UV Fe I lines at usable depth and it yields nothing.**

| | n |
|---|---|
| candidates | 901 |
| refused by the gap test | 845 |
| **passed the isolation escape** | **56** |
| blend-dominated / ghost / pinned | 43 |
| **usable** | **0** |

The 56 that passed came back with **EW 107–922 mÅ**, REW **−4.51 to −3.59** — every one
above the −4.90 saturation ceiling. A single solar Fe line is 10–150 mÅ. **922 mÅ is a
blended complex, not a line.**

## Why the gap test passes lines that are not isolated

It measures the distance to the nearest **catalogued** neighbour, and near-UV catalogue
completeness is **~77 %**. Roughly a quarter of the features actually present are not in
our list, so a 0.30 Å "gap" routinely contains an uncatalogued line.

**Absence of a neighbour in the catalogue is not absence in the spectrum** — the same
lesson the IR root-cause split produced from the other direction, where 98 of 174 failures
were windows dominated by a real solar line we do not carry.

## Conclusion

**The band-wide ban was correct.** Synthesis really is the only valid route in the near-UV,
and it remains blocked on a line list below 4200 Å. The escape mechanism is retained —
it is the right shape for a band whose crowding is genuinely marginal — but it is recorded
in `band_policy.py` as tested-and-negative for the near-UV, so it is not re-tried as a
workaround.

**Near-UV Fe is therefore blocked on ONE thing: a level-identified line list covering
3000–4200 Å.** Not the handler, not the control, not Engine A.

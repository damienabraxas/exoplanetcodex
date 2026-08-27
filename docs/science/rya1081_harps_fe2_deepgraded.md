# RYA-1081 — why HARPS Fe II VIS DEEPGRADED reads 7.966

**Bottom line: it is not a few bad lines. It is a coherent ~+0.25 dex HARPS-specific offset
across the whole pool, on lines whose gf are pure laboratory and identical to KP's — plus
one line that is broken on both arms. The 7.966 should not anchor anything, and the reason
is not that it needs a correction: it is that a nine-line pool spanning 1.4 dex with
`red_chi2` of 68–171 on every member does not support a three-decimal value.**

## The pool

Nine lines, 4233–4584 Å. **Every one carries `PRIMARY LAB DenHartog2019` gf at tier `LAB`
with a cited σ of 0.05–0.10 dex** — the best gf we have, and the same nine lines RYA-853's
DH19 referee used.

| λ (Å) | log gf | σ | A(Fe II) | dev. from median | red_chi2 |
|---|---|---|---|---|---|
| **4303.170** | −2.52 | 0.10 | **10.011** | **+2.045** | **1212** |
| 4351.762 | −1.95 | 0.07 | 7.185 | −0.781 | 103 |
| 4385.377 | −2.64 | 0.07 | 8.560 | +0.594 | 171 |
| 4583.829 | −1.94 | 0.05 | 7.450 | −0.516 | 93 |
| 4522.628 | −2.29 | 0.06 | 8.372 | +0.406 | 87 |
| 4233.162 | −2.02 | 0.05 | 7.619 | −0.347 | 124 |
| 4549.466 | −1.92 | 0.06 | 8.081 | +0.115 | 100 |
| 4508.280 | −2.42 | 0.06 | 7.896 | −0.070 | 69 |
| 4416.819 | −2.57 | 0.09 | 7.966 | 0.000 | 68 |

The published 7.966 is this pool's **median**.

## 🔴 The ticket's framing does not survive the arithmetic

The hypothesis was that a few lines skew the value. **They do not.** Removing 4303.170 —
the worst line in the pool by a factor of four — moves the published median by **−0.035 dex**
(7.966 → 7.931). A median is robust to exactly the failure mode suspected. The pool is not
being dragged up by outliers; **the pool is high.**

## The KP cross-check settles where the offset lives

Same nine lines. Same pure-lab gf. Only the arm differs.

| comparison | median HARPS | median KP | Δ | HARPS > KP |
|---|---|---|---|---|
| vs KP molecfit | 7.966 | 7.569 | **+0.249** | **8 of 9** |
| vs KP kurucz2005 | 7.966 | 7.592 | **+0.263** | **8 of 9** |

The one exception is **4303.170**, which reads 10.011 on HARPS and 10.075 on KP — high on
*both*.

Per the ticket's own diagnostic split, that is decisive:

* **4303.170 is high on both arms** → a line-physics problem, not HARPS-specific.
* **The other eight are high on HARPS and normal on KP** → **HARPS-specific**, and it is an
  offset affecting the whole pool rather than a handful of members. A sign test on 8 of 9 in
  one direction is *p* ≈ 0.018.

**gf is excluded as the cause.** Both arms use the identical `PRIMARY LAB DenHartog2019`
values; a gf error moves both arms together and cannot open a gap between them.

## It is *not* the RYA-911 residual

RYA-911's pathology is in the **EW step**, and it ran **0.34 dex low**. This product has
**no EW leg at all** — `ew_mA`, `rew` and `observed_depth` are NaN on every row, and
`ew_method` records a synthesis flux-fit. The route is different and the sign is opposite.
So this is a **separate synthesis-route HARPS-specific offset**, and routing it to RYA-911
would attach it to the wrong mechanism.

## ⚠️ Two of this ticket's asks cannot be executed today

Declared, not silently skipped (RYA-833).

**1. There is no GRADED (non-deep) Fe II tier to compare against.** Every Fe II product in
the feed — every instrument, every band — is `tier=DEEPGRADED`. So the DEEPGRADED-minus-
GRADED deep contribution has **no comparand**. Splitting this pool on a depth proxy to
manufacture one would be precisely the conflation the ticket names CRITICAL. **What is
needed is a GRADED-tier Fe II run.**

**2. The ratified per-line COG covers Fe I only.** `rya1041_perline_cog_4200_6910.csv`
holds 626 lines over 4204–6903 Å and **zero** of this pool's nine. So the saturation question
is **UNANSWERED, not answered negatively**. Asserting saturation from EW magnitude instead
is a CRITICAL substitution — and moot regardless, since this route measures no EW.

## Disposition: one line, on measured evidence, and it is not the fix

**Fe II 4303.170 → `problem_children`, `BAD_FIT`, `required_treatment=exclude`.**

Evidence carried with it: A = 10.011 (HARPS) / 10.075 (KP), **+2.0 dex above its own pool
median on both arms**, with `red_chi2` 1212 / 1133 against a pool median of ~99, on an
identical pure-lab gf.

Three things it is explicitly **not**:

* not a saturation claim — the COG cannot rule on it, and inventing the verdict is the
  CRITICAL failure this ticket names;
* not a diagnosis — blend, damping and continuum are all untested, so the **mechanism is
  undiagnosed** and recorded as such;
* not what closes the arm gap — it moves the median by 0.035 dex. **Anyone reading this
  exclusion as the repair has misread it.**

The other eight are **not** excluded. They are a coherent systematic, and excluding lines
until an arm agrees with another arm is the tuning RYA-161/RYA-523 forbid.

## ⚠️ A correction to this audit's own first pass

The first cut reported `A` vs `red_chi2` correlating at **+0.86** and read it as "the high
values are a fit-quality artifact". That correlation is **driven entirely by 4303.170**. On
the other eight it is **+0.31** (HARPS) and +0.59 (KP). **Fit quality does not explain the
offset**, and the stronger claim would have been wrong.

## What the trustworthy HARPS Fe II VIS value is

**There isn't one yet, and that is the finding.** Not because the number needs correcting,
but because the pool cannot support one:

* n = 9, spanning **1.375 dex** even with the broken line removed;
* `red_chi2` of **68–171 on every surviving member** — no line in this pool is well fitted;
* a **+0.25 dex arm offset** whose mechanism is unidentified;
* the deep contribution **unquantified**, because no GRADED Fe II tier exists;
* saturation **unruled**, because the COG has no Fe II coverage.

The residual offset is, on this evidence, an **artifact rather than physics** — it is
arm-specific on identical gf, which physics cannot be — but its mechanism is not identified
here, and the honest next steps are a GRADED-tier Fe II run and a per-line COG pass over
Fe II, not a correction to 7.966.

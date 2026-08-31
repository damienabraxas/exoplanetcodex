# RYA-1162 — per-line artifacts for the 3 live IR Fe products

27 artifacts landed (9 per arm: `_lines.csv`, `_budgets.txt`, `_provenance.txt` for
1D-LTE / ENGINE-A / ENGINE-B). **No published value moved.**

## Reproduction — the gate the ticket set

Every published field reproduced **byte-identically**. All nine committed
`_products.csv` carry the same sha256 after the runs as before (CRIRES+ and KP rewrote
their own files with identical bytes, which IS the proof; RYA-1084 determinism holds):

| arm | A | n | excluded | scatter | vs published |
| -- | --: | --: | --: | --: | -- |
| IAG 1D-LTE | 7.546 | 25 | 3 | 0.794 | exact |
| IAG ENGINE-A | 7.599 | 6 | 19 not-served | — | exact |
| KP 1D-LTE | 7.621 | 26 | 2 | 0.890 | exact |
| KP ENGINE-A | 7.503 | 7 | 19 not-served | — | exact |
| CRIRES+ 1D-LTE | 7.551 | 5 | 0 | 0.138 | exact |
| CRIRES+ ENGINE-A | 7.492 | 2 | 3 not-served | — | exact |

## 🔴 F1 — NO TELLURIC LINE FILTER RUNS ON THE BAND-PRODUCT ROUTE

This is the finding. It is a live defect, it inflates every IR bar, and it explains
RYA-1114's F4 completely.

**The mechanism.** `telluric_policy.exclusion()` (band membership) was **superseded as a
decision input** by RYA-1079, which replaced it with `telluric_observability`'s measured
three-state per-line policy (CLEAN / RECOVERABLE / SATURATED). The successor was **never
wired into the measurement path.**

    partition_pool  production callers : scripts/rya1079_observability_census.py   (a CENSUS)
    derive_band_products.py            : does not import telluric_observability at all

So the old policy no longer decides, the new policy is not called, and **nothing filters
telluric-contaminated lines out of a synthesis band product.** `telluric_observability`'s
own docstring asserts the opposite — *"`partition_pool` runs BEFORE tier assignment, and
a SATURATED line never reaches the grader"* — which is not true on this route.

**Measured consequence.** Splitting each pool by membership of an enumerated
`TELLURIC_BANDS` interval:

| arm | inside a telluric band | outside |
| -- | -- | -- |
| IAG | **n=11, median 7.436, std 1.120, range 3.629** | n=14, median 7.549, **std 0.194** |
| KP | **n=13, median 7.598, std 1.136, range 4.301** | n=13, median 7.690, std 0.600 |
| CRIRES+ | **n=0** | n=5, median 7.551, **std 0.138** |

Every in-band line is in the H₂O **9280–9600 Å** complex, and all of them are
`in_aggregate`. Individual kept lines reach **A = 4.866** (IAG 9454.194 Å) and
**A = 10.937** (KP 9437.793 Å) — 2.7 and 3.3 dex from solar — with `red_chi2` up to
**309**. These are not measurements; they are failed fits being averaged in.

⚠️ IAG's holding is `telluric_basis = corrected` (Baker+2020, corrected at source), and
its in-band lines still scatter **1.120 dex**. So the RYA-1114 brief's Part E-a premise
— *telluric-CORRECTED is not telluric-CLEAN* — is now MEASURED, not asserted.

## 🔴 F2 — RYA-1114's F4 and F6 are REFUTED by this data

RYA-1114 concluded the wide IR bars "track line/gf quality, not n". **Wrong.** The
deriver reports for all three arms:

> gf rung 3 (gf scale (cited lab)): **every one of the 25 / 26 / 5 Fe I lines is GF-LAB
> and 100% carry a published per-line sigma; RMS 0.0537 / 0.0538 / 0.0537 dex**

All three pools are **100% laboratory gf with identical RMS**, yet their dispersions are
0.794 / 0.890 / 0.138. gf quality cannot explain a 6x difference it does not vary across.

**F6 is refuted outright.** RYA-1114 reported IR as "5.9% LAB, 22 of 34 ungraded" — that
described `Fe_perline.csv`, a stale 2026-08-18 PROJECTION which RYA-1114 itself flagged
as INDICATIVE and which reconciles with no live product. The live IR pools are the
**exhaustive-LAB** case after all, like near-UV.

**The real discriminator is band placement.** CRIRES+ spans 9801–10794 Å and contains
**no** telluric interval; IAG and KP both span across the H₂O 9280–9600 complex. That,
not gf and not n, is the 6x.

## ⚠️ F3 — the IAG published stem predates RYA-1046 and overstates its red edge

The IAG run wrote stem `9199_11081`, not the published `9199_11083`. Cause traced:
RYA-1046 narrows a band so every fitting window lies wholly inside the data
(`9199.0–11083.0 → 9199.0–11081.7`). That commit (`08c65bf`) landed
**2026-08-25T19:53Z**; the committed IAG artifact was built **18:40Z** — 73 minutes
EARLIER. So the published stem records the REQUESTED range from before the narrowing
existed and overstates the red edge by ~1.3 Å.

Content is byte-identical (same 25 lines, same 7.546, same bars). Only the name differs.

**Not resolved here, deliberately.** The artifacts are landed under `9199_11081` — the
range actually covered — because renaming them to `11083` would assert coverage the run
did not have. Re-stemming the PUBLISHED product changes the name the feed points at, so
it is a publication decision, not this ticket's (the ticket says "do NOT re-cut the band
edges"). The three `9199_11081_*_products.csv` this run produced were deliberately NOT
committed, to avoid a duplicate product under a second stem.

## ⚠️ F4 — ENGINE-A drops are `not-served`, measured against the LIVE service

Both ENGINE-A legs queried MPIA live (*"Bergemann MPIA per-line delta_nlte (live query,
solar node)"*) and reported the dropped lines as **`not-served`**: IAG 19, KP 19,
CRIRES+ 3. This answers RYA-1165's open F7 question toward **service reach**, not an
RYA-1050 labelling defect — measured, not assumed. The per-line files carry it.

## Verification

* 9/9 committed `_products.csv` sha256 unchanged; `git status` shows no tracked
  modification outside the added artifacts and `GENERATORS.yaml`.
* `check_result_generators.py` green (516 → 543 tracked artifacts, all registered).
* Runs executed with the xi-campaign's proven Mac env (`CODEX_KP_ATLAS`, `ISPEC_DIR`,
  `SIRIUS_DATA_ROOT`, `CODEX_SRV_SOLAR_ROOT`); the Baker FITS is mirrored locally so IAG
  needed no Sirius access.

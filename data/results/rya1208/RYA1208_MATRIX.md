# RYA-1208 — the Gerber matrix, cell by cell

Every (holding × band × Gerber-treatment) cell is a published product or an N/A with a
stated reason. No value is calibrated (RYA-161).

The Gerber route has two legs:

* **`synth-1D-LTE-gerber`** — Turbospectrum LTE on the Gerber deck's *own* atmosphere
  (MARCS.GES), departures withheld. **Needs no NLTE labels**, so it is runnable wherever
  the synthesis route reaches. Verified in source: `assert_linelist_supports_nlte` is
  called only under `if _nlte:`.
* **`ENGINE-B-NLTE`** — the same deck *with* departures. Needs per-line NLTE labels,
  because bsyn falls back to departure = 1 for any unlabelled line and would publish an
  LTE spectrum under an NLTE name (RYA-764/1050).

## Applicability, measured before running anything

### NLTE label coverage per band

| band | line list | Fe lines | NLTE-labelled | ENGINE-B-NLTE |
| --- | --- | --- | --- | --- |
| near-UV | `ispec_nearuv_3000_3780` | 5,038 | **0** | **N/A — no labels** |
| NIR | `ispec_ir_9200_13000` | 237 | **0** | **N/A — no labels** (confirms RYA-1203) |
| H | `ispec_h_15007_17494` | 1,935 | 502 | available — published by RYA-1206, n=9 |

So the ticket's open question ("near-UV: DETERMINE label coverage first; if 0 labels,
document as N/A") is answered: **0 of 5,038. N/A.**

### Fe II NLTE — N/A everywhere, and it is not a labelling problem

`gerber_nlte.nlte_ion_capability('Fe', 2)` returns `False` with its own reason:
`atom.fe607a` carries 58 Fe II levels but **zero** Fe II bound-bound transitions
(12,635 of 12,635 are Fe I → Fe I), so bsyn applies departure = 1 to every Fe II line
and **no line list can enable it** (RYA-1055). Labels are irrelevant to the outcome.

### Band coverage — the two cells that cannot be reached at all

| cell | measured | verdict |
| --- | --- | --- |
| HARPS near-UV | holding starts **3782.6 Å**; the near-UV product window is 3000–3780 Å | **N/A — overlap is empty.** (19 Fe lines sit in 3782.6–3800 Å, between the near-UV product window's red edge and VIS's 3800 Å start — reachable by neither band as currently defined.) |
| ~~KP-kurucz2005 NIR~~ | ~~10 Å of overlap, zero Fe lines~~ | 🔴 **WRONG — see the correction below. It is runnable and was measured.** |

## 🔴 A blocker found on the first run, fixed here

`--engine-b-deck gerber-1d-lte` died immediately with
`UnboundLocalError: cannot access local variable '_unlabelled'`. RYA-1206 bound that name
inside `if _nlte:` and reads it on every deck path. **Every non-NLTE Engine-B deck was
broken on main, including `ts-lte`, the default.** It survived because the only runs since
either used an NLTE deck or passed `--skip-engine-b`. Fixed, with a static regression test
whose control is the pre-fix source.


## The ENGINE-B-NLTE refusals, in the harness's own words

Not inferred from the label count — `gerber_nlte.assert_linelist_supports_nlte` was called
with the arguments `derive_band_products` passes it, and these are its verbatim replies.
Captured directly rather than by running a full cell, which reaches the guard only after
~45 minutes of LTE fitting.

**near-UV, Fe I**
> Fe: 4364 lines in 3000.0-3780.0 A but NONE carry NLTE level labels. bsyn sets departure = 1 for an unidentified level and iSpec drops an unlabelled element into `nlte_ignored` — either way the synthesis runs in LTE WITHOUT RAISING (RYA-534 Co/Ni, RYA-764). Refusing: an LTE spectrum under an NLTE label is worse than no product. This is a LINE-LIST coverage gap, not a deck failure — the atom carries

**NIR, Fe I**
> Fe: 235 lines in 9199.0-12976.0 A but NONE carry NLTE level labels. bsyn sets departure = 1 for an unidentified level and iSpec drops an unlabelled element into `nlte_ignored` — either way the synthesis runs in LTE WITHOUT RAISING (RYA-534 Co/Ni, RYA-764). Refusing: an LTE spectrum under an NLTE label is worse than no product. This is a LINE-LIST coverage gap, not a deck failure — the atom carries

**Fe II, any band** — refused one level earlier, on the atom rather than the list:
> Fe II: NLTE was requested, but this deck's model atom carries NO bound-bound transitions for this ionisation stage, so bsyn would apply departure = 1 to every line and the product would be LTE under an NLTE label (RYA-764/RYA-1050). NO LINE LIST CAN FIX THIS — label coverage is irrelevant when there is nothing in the atom to label against. Fe II NLTE unavailable on the Gerber deck: atom.fe607a car

## ⚠️ The Gerber deck is Sirius-only

`pipeline.gerber_nlte` reads `/mnt/codex-ext/codex-grids/nlte/gerber_ts/auxData_Fe_MARCS_*.dat`
and raises `GerberDeckError: missing Gerber aux table (Sirius-only)` anywhere else. **No
Gerber leg — LTE or NLTE — can run on the Mac.** Worth knowing before anyone tries: a Mac
attempt completes the entire 1D-LTE synthesis first and only fails when Engine-B loads the
deck, so the failure arrives ~45 minutes in and looks like a late crash rather than a
missing prerequisite.


## 🔴 A correction: KP-kurucz2005 NIR is NOT N/A, and the error was mine

I first ruled it out by testing the holding's span (2990–10010 Å) against the **band-policy**
NIR window (10000–24000 Å): 10 Å of overlap holding zero Fe lines. That is the wrong window.

Products do not run to the policy edges. The NIR products run **9199–12976 Å**, and
kurucz2005 overlaps *that* by **811 Å carrying 92 Fe lines, 24 of them lab-tier**. Measured:

    A(Fe I; NIR, synth-1D-LTE-gerber) = 7.515   (n=23, excluded 0)

This is precisely RYA-1069's lesson — *ask at the reader's overlap, not at the band
definition* — which this very document quotes two sections above. It surfaced only because
the cell was **attempted** rather than asserted N/A, which is the whole reason the two
coverage cells were probed on Sirius instead of argued on paper.

⚠️ **This holding had no NIR product of any kind**, so the same run produces the first
1D-LTE (7.520, n=23) and ENGINE-A (7.385, n=6) NIR cells for kurucz2005. Those are outside
the Gerber matrix this ticket asks for; they are reported here rather than quietly kept or
quietly dropped.

### The corrected coverage table — product windows, every holding × band

| holding | band | overlap (Å) | Fe lines | had a product? |
| --- | --- | --- | --- | --- |
| HARPS-molecfit | VIS | 3800–6910 | 2434 | yes |
| KP-kurucz2005 | near-UV | 3000–3780 | 987 | yes |
| KP-kurucz2005 | VIS | 3800–6910 | 2434 | yes |
| KP-kurucz2005 | red-optical | 6910–9199 | 474 | yes |
| **KP-kurucz2005** | **NIR** | **9199–10010** | **92** | **NO — new cell** |
| KP-molecfit | near-UV / VIS / red-optical / NIR | — | — | yes |
| IAG | VIS / red-optical / NIR | — | — | yes |
| CRIRES+ y_wide | NIR | 9800–10796 | 70 | yes |
| CRIRES+ h | H | 15009–17491 | 0* | yes |

\* `per_line.csv` does not reach past ~13000 Å (the RYA-1094 stranded-lines gap), so the
zero is an accounting-table limit, not an absence — the H products exist and carry n=19–20.

**HARPS near-UV is the only genuine coverage N/A**, and the harness says so itself:

> harps/solar_harps_molecfit_corrected spans 3782.6-6910.0 A and the request is
> 3000.0-3780.0 A — they do not overlap by even one fittable window (half-width 0.40 A).
> This arm cannot serve this band at all (RYA-832).

## 🔴 The two Fe I near-UV cells are BLOCKED, and the block is a finding

`publish_product` refuses them:

> tier=DEEPGRADED claims laboratory gf, but the budget's own verdict is gf rung 1
> (UNGRADED). A tier is a claim about the gf SCALE, and the stem is not evidence for it.

The budget is right. Fe I near-UV is a **mixed pool** — 56 of 57 lines are GF-LAB and one
is `systematic:K07`, and *"a pool is graded only if every line in it is."* `canonical_gf`
starts at 3780.04 Å, so the band has no lab-gf single-sourcing at all (RYA-759's declared
provenance gap).

**The existing near-UV Fe I products carry the same tier against the same verdict.** Their
budget also reads `gf rung 1 (UNGRADED)` — 54 of 55 lab, one `systematic:K07`. They exist
only because they were ingested **83 minutes before the guard that would refuse them**:

| | timestamp |
| --- | --- |
| the four near-UV Fe I products ingested | 2026-08-24T23:30:21Z |
| RYA-1034's tier/rung guard committed (`c870de94`) | 2026-08-25T00:53:04Z |

So four published products claim a laboratory-gf tier their own budget refuses, and were
never checked. My identical new cells ARE checked, and correctly refused.

**Not resolved here, because both routes change a published claim.** Forcing mine through
under DEEPGRADED would replicate an unsupported pedigree; publishing them as UNGRADED
would put them in a different cell from the 1D-LTE sibling they exist to be compared
against; re-tiering the existing four changes four published identities. **Ryan's call.**

The two cells are measured and the values are in hand either way:

| cell | A(Fe I) | n | blocked by |
| --- | --- | --- | --- |
| near-UV KP-kurucz2005 Fe I | **7.580** | 57 | tier/rung guard |
| near-UV KP-molecfit Fe I | **7.475** | 57 | tier/rung guard |

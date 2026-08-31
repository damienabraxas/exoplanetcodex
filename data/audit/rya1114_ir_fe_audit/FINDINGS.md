# RYA-1114 — IR Fe xi-aware uncertainty re-verification

READ-ONLY audit. **No value changed, no product moved, no tuning (RYA-161).**
Regenerate: `python3 scripts/rya1114_ir_fe_audit.py`

Scope: the **6 live IR Fe products** in `data/products/solar/Fe.json` (`band == NIR`),
read from the FEED and not from the `band_products/` directory (RYA-1097 lesson: a
directory listing is not the live set). The 3 `ENGINE-B` NIR rows are already
`superseded[]` and are excluded; they are byte-identical duplicates of the `1D-LTE`
rows, which is RYA-1100's point.

## Pre-flight

**The RYA-920 gate was dropped, with cause.** The brief blocks this audit on RYA-920
("IR runs on IAG, and auditing on a lying IAG registry inherits the lie"). RYA-920 has
not landed — zero commits anywhere, and `instrument_catalog.csv` still carries IAG at
`404.7 – 1065.0` nm. But the gate does not bite this audit:

* NIR is 9199–13000 Å, entirely above Baker's 5001.1 Å floor, so the Reiners blue arm
  is span-capped out of the IR path and cannot reach it (register v116).
* The live IAG IR product is `FeI_9199_11083` — cut to the FITS maximum, which I
  re-measured on Sirius through the repo's own `vac_to_air`: **air 5001.1076 –
  11083.4318 Å** (HDU 1 BinTable, `TTYPE1='v'`, vacuum wavenumber cm⁻¹, 728,122 rows;
  Reiners for contrast: air 4047.4255 – 10649.8663, 4,057,462 rows). The product rests
  on the FILE, not on the wrong catalogue row.

⚠️ The catalogue row is still wrong and is still read by `preflight_check.py:1264` and
`measure_reference_lineset.py:90`. It did not bite here; it is not fixed.

**Telluric verification — all three arms pass, no HOLD.** Through
`telluric_policy.applied_state` / `gate_holding`, not from a flag:

| holding | state | basis |
| -- | -- | -- |
| `solar_iag` | applied | `corrected` — Baker+2020 telluric-corrected at source |
| `solar_kpno_molecfit_corrected` | applied | **`line_selection`** — avoidance, not correction |
| `solar_crires_plus_y_wide_rya1054` | applied | `telluric_applied` — cr2res+molecfit |

⚠️ KP's basis is `line_selection`, and **RYA-928 is open**: KP declares that strategy and
then keeps 61 red-optical lines inside enumerated telluric bands. Measured here, the NIR
region does **not** inherit that failure — 0 of 34 line centres fall inside a band (but
see F5).

## Findings

### F1 🔴 `sigma_syst` omits the stellar-parameter term on all 6 IR products — CONFIRMED

Rebuilding the budget through the same `error_budget.build()` call
`derive_band_products.py:1556` makes returns:

```
line-to-line scatter   dex=<scatter>   random
gf scale (NIST-graded) dex=0.041
harness residual       dex=0.0         (SynthesisHandler)
stellar parameters     dex=None        <-- NOT MEASURED
telluric residual      dex=0.03
```

A term with `dex=None` contributes nothing to `total()` and leaves **no trace in the
published float**. `derive_band_products.py:1556` never passes
`stellar_param_sigma_dex`. This is RYA-1112's F1 (→ RYA-1120) live in IR — and IR was
**not** covered by the RYA-1120 fix (F2).

### F2 🔴 All 6 IR products are ABSENT from the RYA-1120 `sigma_reported` harness

The brief requires "Type B with per-pool dA/dxi (RYA-1120, NOT the single RYA-1089
number)". That harness holds 46 products and **none is IR**;
`solar_crires_plus_y_wide_rya1054` is not among its holdings at all.

Cause, traced: the xi campaign (`~/xi_campaign/out`, 30 runs = 15 pools × 2 xi) keys
pools as `Fe{I,II}_<holding>_<TIER>_<ROUTE>` — **with no band**. `FeI_solar_iag_GRADED_SYNTH`
resolves to `FeI_5002_6910…` — the **VIS** product. IR was never perturbed. The harness
refuses to borrow a derivative ("it will not difference two aggregates… the derivative
is a per-line paired differential or it is nothing"), and IR has no per-line data (F3),
so `sigma_params` is not merely missing for IR — it is currently **unobtainable**.

⇒ The xi-aware bar this ticket asks for **cannot be assembled for any IR product today.**
Every IR bar below is a FLOOR.

### F3 🔴 No per-line artifact exists for ANY live IR product — the audit is BLIND

None of the six has a `_lines.csv`, `_budgets.txt` or `_provenance.txt`. The only
per-line NIR data in the repo, `data/products/solar/Fe_perline.csv` (34 rows, **kpno
only**), was generated **2026-08-18** from `rya847/gated` + `rya877` — a week before the
live NIR artifacts (mtime **2026-08-25**) — and its 34 rows reconcile with no live
product (KP NIR is n=26 + 2 excluded = 28). So IAG (n=25) and CRIRES+ (n=5) have **zero**
per-line visibility, and the KP numbers below are INDICATIVE of the region, not the pool.

⚠️ This is a committed-artifact gap, **not** a structural one: the xi-campaign outputs
show `derive_band_products` does emit `_lines.csv`. Re-running would recover them.

### F4 🔴 The >0.1 bars are NOT a small-n story — the ticket's canonical framing is refuted

`sigma_stat` is a standard error (`error_budget.py:609`), so the underlying dispersion is
`sigma_stat x sqrt(n)`:

| arm | treat | n | sigma_stat (SE) | raw scatter | >0.1 |
| -- | -- | --: | --: | --: | :--: |
| CRIRES+ | 1D-LTE | 5 | 0.0616 | **0.138** | |
| CRIRES+ | ENGINE-A | 2 | 0.0700 | **0.099** | |
| IAG | 1D-LTE | 25 | 0.1588 | **0.794** | ✔ |
| IAG | ENGINE-A | 6 | 0.0721 | **0.177** | |
| KP | 1D-LTE | 26 | 0.1746 | **0.890** | ✔ |
| KP | ENGINE-A | 7 | 0.4970 | **1.315** | ✔ |

The ticket names the KP n=7 / 0.497 cell as "the canonical case — verdict = irreducible
small-n, fix = more/better NIR Fe lines". **Measured, that is the wrong diagnosis.**

* The two other flagged products carry the **largest** n in the set (25, 26) and still
  the widest bars. More lines did not help them.
* The arm with the **fewest** lines, CRIRES+ (n=5), has the **tightest** dispersion —
  0.138 dex, ~6x tighter than IAG at n=25.
* The n=7 cell's problem is a **1.315 dex line-to-line dispersion**. Adding lines shrinks
  the SE while leaving that dispersion untouched; the bar would fall without the
  measurement improving.

RCA verdict: the flagged bars track **line/gf quality, not n**. CRIRES+ is a curated
GRADED lab-gf set; the IAG/KP NIR pools are dominated by ungraded VALD3 (F6). Adding more
VALD3-grade NIR lines is predicted to *widen* the dispersion, not narrow it.

### F5 ⚠️ Residual-telluric (Part E-a) — one line, INDICATIVE

0 of 34 KP NIR line **centres** fall inside a `TELLURIC_BANDS` interval. But
**11119.795 Å sits 0.205 Å below the H₂O band edge (11120.0)** and is `flagged_kept`.
The NIR synthesis half-width is **1.40 Å**, so that line is measured over
[11118.40, 11121.20] — overlapping the H₂O band by ~1.2 Å. Band membership is tested on
the centre; the measurement integrates over the window. IAG and CRIRES+ cannot be checked
for this at all (F3).

### F6 ⚠️ gf floor (Part E-c) — IR is NOT the near-UV's exhaustive-LAB case

In the indicative KP NIR set: **2 of 34 lines (5.9%) are PRIMARY LAB** (Ruffoni2014,
per-line sigma 0.02); 32 are VALD3, and **22 of 34 are `ungraded`** while the product
publishes at tier `GRADED`. RYA-1113 found near-UV was 1% LAB in the list but 100% LAB in
the product; IR is the opposite shape. This is the reducibility story: IR **is** genuinely
gf-floored, which is also the mechanism behind F4.

### F7 🔴 Part E-b's premise is refuted — the state is REACH-UNKNOWN, a third case

The brief asks to classify each IR ENGINE-B by whether `atom.fe607a` "ceilings at 7.505 eV,
below IR upper levels". Three independent checks say that framing does not hold:

1. **7.50 is an ABUNDANCE, not an energy.** `gerber_nlte.deck_abundance` documents it:
   *"Fe's grid was computed at A(Fe) = 7.50: that is what `atom.fe607a` declares on its
   own second line (`7.50  55.85`)"* — dex, not eV.
2. **Measured upper levels do not approach it anyway.** Over the 34 KP NIR lines,
   `E_upper = EP + hc/lambda` spans **3.211 – 7.351 eV**; **0 of 34** exceed 7.505.
3. **The atom is not the IR limit.** `coverage.py:393`: *"The 9199.9 Å wall is a linelist
   limit, not physics (atom.fe607a reaches 20000 Å)."*

The recorded state in `data/catalog/engine_coverage.csv` is neither "cannot reach" nor
"can reach but unlabelled" but **`REACH-UNKNOWN`**, for both engines in NIR:

* Fe I / Engine A / NIR — *"no local Engine-A level table (no `label_Fe.txt`) — this
  element's Engine A is a web-service supplier, so its reach is not locally decidable"*.
  Fe's ENGINE-A is the **Bergemann/MPIA web service**; there is no local grid to interrogate.
* Fe I / Engine B / NIR — *"no catalogued line in band — a zero here measures the
  LINELIST's span, not the grid's, so absence cannot be asserted from it"*.

⇒ The classification the ticket asks for **cannot be made from local files by
construction**, and RYA-776 built `REACH-UNKNOWN` precisely so that a zero is not read as
an absence. Answering it needs the live MPIA service, not an atom file.

⚠️ Note the ENGINE-A pools drop most lines (IAG 6 of 25, KP 7 of 26 — 19 excluded each).
Whether that is service reach or a labelling defect (RYA-1050) is **exactly** what F3's
missing per-line artifacts would decide, and it is undecidable today.

### F8 🔴 ENGINE-A breaches the holding-spread bound

| treatment | CRIRES+ | IAG | KP | spread | vs 0.10 (RYA-1133) |
| -- | --: | --: | --: | --: | -- |
| 1D-LTE | 7.551 | 7.546 | 7.621 | 0.0750 | within |
| ENGINE-A | 7.492 | 7.599 | 7.503 | **0.1070** | **BREACH** |

Three arms measuring the same species in the same band disagree by 0.107 dex under NLTE,
driven by IAG (7.599) vs CRIRES+ (7.492). Under LTE the same three agree to 0.075. The
disagreement is therefore **introduced by the NLTE leg**, on pools of 6, 7 and 2 lines —
and F7 says we cannot currently tell whether those pools are even the same physical
subset. The published IR bars do not cover this spread.

### F9 ⚠️ Two band tables still disagree (RYA-1094 / register v128 shape)

`band_policy.POLICIES` NIR = **10000–24000 Å**; `config/synth_bands.yaml` NIR =
**9199–13000 Å**. The IAG product spans 9199–11083, so its sub-10000 Å lines fall in
`band_policy`'s *red-optical*. Both carry `telluric_required=True`, so the budget term is
unaffected here — but `engine_coverage.csv` is keyed on the `band_policy` boundaries, so
its NIR row does not describe the 9199–10000 Å part of the product at all.

## Verdict

All 6 products pass the telluric gate; **no HOLDs**. Three exceed the 0.1 dig-in
threshold, and their RCA is **gf/line quality, not small-n** (F4/F6) — a documented
irreducible under the present line pool, and one that "more NIR lines" would not fix.

But the audit this ticket specifies **cannot be completed** as written: the xi-aware
budget is unobtainable for IR (F2), the per-line classifications are blind (F3), and the
model-domain question is not locally decidable (F7). The published IR bars are FLOORS
that omit the stellar-parameter term entirely (F1) and do not cover the 0.107 dex
inter-arm NLTE spread (F8).

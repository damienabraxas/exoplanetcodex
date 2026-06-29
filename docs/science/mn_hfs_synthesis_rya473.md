# Solar Mn via HFS-resolved synthesis — RYA-473

**The real unlock.** RYA-468 settled the Mn gf disposition and found — empirically, on
the solar spectrum — that the gf grade was never Mn's final blocker; **saturation is**.
The Den Hartog e6S→z6P triplet (6013.51 / 6016.67 / 6021.82) now carries a graded gf
(Den Hartog+2011, MED tier) yet SAT-culls through the EW path (measured EW ~90–100 mÅ,
REW −4.78..−4.82, over the −4.90 saturation knee; each line is HFS-split, hfs_n=6). A
single-profile EW fit over-saturates them → KEPT=0, Mn stays no-value. The gf is correct;
the lines needed the right **method**. This ticket measures Mn by **HFS-resolved
synthesis**, the same path RYA-466 used for Cu and RYA-411 validated for Mn.

## Method (reuse, not rebuild)

- **Engine:** `pipeline.cno_synthesis._fit_element` (synthesis v2, RYA-338) — flux-space
  χ²ᵣ fit against the observed HARPS solar spectrum, varying only A(Mn) per line.
- **HFS + gf:** the GES atomic line list the production synthesis already loads
  (`GESv6_atom_hfs_iso`, RYA-353) carries the ⁵⁵Mn hyperfine pattern of the triplet
  **natively — 6 resolved components per feature**, gf-summing to exactly the Den Hartog
  values (−0.354 / −0.181 / −0.054). Those components are the Den Hartog Lawler Sobeck
  Sneden Cowan 2011 (ApJS 194,35) hfs pattern (GES tag `DLSSC` = **Den Hartog Table 4**) —
  **cited, not invented here**. The script verifies the GES gf-sum equals the RYA-468
  canonical value before trusting the fit.
- **NLTE:** the MPIA/Bergemann Mn NLTE grid (`data/nlte_grids/Mn_Bergemann_MPIA.csv`,
  on disk) is the registered Mn source. The script first **tries the live Amarsi
  HFS-resolved delta** (RYA-411 path); the Amarsi GALAH PySME departure grid
  (`nlte_Mn_scatt_pysme.grd`) is offline, so it falls back to the **MPIA grid solar δ
  (+0.108), LOUDLY FLAGGED**.

## NLTE caveat (RYA-411) — carried loud

The MPIA grid's nodes are the **high-EP** lines 4998 / 6304 / 6306 / 6867. The low-EP,
strongly-HFS triplet **6013/16/21 is NOT a grid node**, and its HFS-resolved Amarsi delta
**differs** (HFS desaturation, RYA-411). So the applied +0.108 is the grid's high-EP
reference value used as a vendored approximation — the triplet-exact NLTE δ requires the
offline Amarsi grid. Mn is therefore **not certified PASS on the vendored δ**.

## Result (validate-don't-tune)

| quantity | value |
|---|---|
| A(Mn)_LTE median | **5.446** (mean 5.357, σ 0.154, n=3) |
| per-line A(Mn)_LTE | 6013.51 → 5.447 · 6016.67 → 5.179 · 6021.82 → 5.446 |
| MPIA NLTE δ (vendored, flagged) | +0.108 |
| A(Mn)_NLTE | **5.554** |
| vs Asplund 2021 (5.42) | **+0.134** |

A(Mn)_LTE sits only +0.026 from Asplund; the vendored MPIA NLTE then lifts it to +0.134.
Nothing was pulled toward the anchor (the per-line values genuinely scatter; the reported
value is the fitted median, not 5.42). The high per-line χ²ᵣ (~120) is the expected
signature of fitting strong saturated features with a tiny photometric σ — a finding to
note, not a fit failure.

## Verdict delta

Mn moves **no-value → MEASURED** via HFS synthesis (folded into `phase_c_verdict_rya371`
through `_apply_mn_synthesis`). It **stays CURATION-OWED** — off no-value, but not PASS,
because (a) the +0.134 offset survives and (b) the NLTE is the vendored MPIA high-EP value,
not the triplet-exact Amarsi δ. **The MEASUREMENT-TOOL blocker is fixed**; what remains is
the line-exact NLTE (RYA-411 Amarsi grid) and fuller curation — a finding, not a tune.
Counts unchanged: **4 / 1 / 21 / 0**; Fe anchor 7.516.

## Owed

- Line-exact HFS-resolved Mn NLTE on 6013/16/21 (Amarsi GALAH `nlte_Mn_scatt_pysme.grd`,
  offline here) — would replace the vendored MPIA high-EP +0.108 and could move Mn to PASS.

**Related:** RYA-468 (the graded gf — input), RYA-411 (HFS-synth machinery + the NLTE
caveat), RYA-466 (the Cu precedent — same pattern), RYA-371 (the verdict).

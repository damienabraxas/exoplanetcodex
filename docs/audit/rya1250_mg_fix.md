# RYA-1250 fix: solar Mg lines as the synthesis reads them

Reopened 2026-10-23 after the post-merge review of PR #585 (main `19b1b4c8`). Script:
`scripts/rya1250_mg_reference_lines.py` (idempotent). Tests:
`tests/test_mg_reference_lines_rya1250_fix.py`.

## Why "right in the table" was wrong in the synthesis

`gf_resolver.apply_to_synth_array` clusters each synthesis list into physical lines (same lower
level, gaps < 0.10 A) and scales each cluster to the ONE canonical row nearest its centroid
(0.02 A / 0.02 eV). A component-level published set is honoured only if each cluster's row carries
the published total of the components inside it. Resolving main's table the way the synthesis
does:

| cluster | main | fixed | what was wrong |
|---|---|---|---|
| Mg I 5183.604 | **+0.180** | -0.168 (PR17 exp) | NIST -0.180 with the sign dropped (every list has -0.24); predates RYA-1250 |
| Mg I 5167.322 | -1.031 | -0.854 (PR17 exp) | |
| Mg I 4571.096 | -5.623 | -5.397 (B17) | |
| Mg I 12083.28 {.278,.346} | -1.347 | -1.116 | one component's gf on a two-component cluster |
| Mg I 10965.45 {.414,.450} | -0.260 | -0.184 | same |
| Mg I 10811.11 (7 VALD comps) | 0.500 (VALD) | 0.401 | B17 rows added beside it, never resolved |
| Mg II 10092.217 | 1.283 | 1.033 (+ .095 at 0.923) | feature total on one of two clusters 0.122 A apart (~0.25 dex excess) |
| K 19723 / 19733, near-UV 3092.98 / 3096.89 | -0.864 / -0.593 / -0.572 / -0.301 | VALD totals | NIST grades only the dominant component, whose value VALD already carries |

The `+0.180` on 5183 alone shifts that line's abundance by roughly -0.3 dex: the main reason the
HARPS VIS Deep product scattered by 0.278 dex.

**12083 identity (correction to the first review).** VALD lists 12083.278 at +0.45. PR17 measured
it (3d 1D2 - 4f 3F3 intercombination) at -1.347 +/- 0.05; the strong line is 12083.662 (4f 1F3,
+0.377 +/- 0.04). The large drop at .278 was a correction, not a defect.

## Reference Grade lines (Asplund + lab gf), by band

Ryan's definition (2026-10-21): Reference Grade = the AGSS21 set plus any lab-referenced line.

* **AGSS21** (`MG_AGSS21`): 8 Mg I lines, Table A.1 gf = PR17 theoretical. Sigma 0.05 dex is
  Bergemann+2017's own convention for PR17 theoretical gf (Table 2 note a; PR17 Fig. 4: theory
  agrees with experiment to a few percent). 15879: AGSS21 feature total, sigma 0.06 (PR17 exp
  error of the dominant component). Mg II: 7 rows (10092 = 2 components), NIST A/A+.
* **Bergemann+2017 Table 2** (`MG_BERGEMANN2017`, was `_IR`): now includes the optical lines
  4571, 4702, 5172, 5183, 5528 and 5711. One row per synthesis feature, with the components
  summed. 15748 is culled for the Sun (blends carry ~half the EW, as B17 do).
* **Lab pool** (`gf_tier=LAB`, `PR2017_exp`): PR17 Table 2 experimental gf for 5167, 5172, 5183,
  5528, 7657.6, 7659.15, 7659.90, 8806, 9255.8, 11828, 12083.66, 15025, 15040, 15048 and 15889.
* gf correlation: `PR2017_exp`, `MG_AGSS21_PR17calc` and `MG_BERGEMANN2017` share one 10% lifetime
  scale (`PR17_LIFETIMES`, scale_only 0.0414), the PR2024 precedent.

| band | holdings | Reference lines |
|---|---|---|
| near-UV | KP x2 | none. The 6 NIST B/C lines are saturated (Codex held) -- declared gap |
| VIS | KP x2, HARPS, IAG, IAG-R16 | 4571 4702 5167 5172 5183 5528 5711 6318 6319 |
| red-optical | IAG, KP x2 | 7657.6 7659.15 7659.90 8712 8717 8806 8923; Mg II 7877 7896 |
| NIR | IAG (<10000), KP, CRIRES Y-wide/Y/J | 9255.8 9983 10312 10811 10965 11828 12083.66; Mg II 9218 9244 10092x2 10914 |
| H | CRIRES H x2 | 15025 15040 15048 15765.6 15765.8 15879 15886 15889 (15748 culled) |
| K | CRIRES K | none -- declared gap, see below |

Coverage notes: 10811/10914/10965 fall between CRIRES Y-wide (to 10796 A) and J (from 11160 A),
so only the KP atlas serves them. PR17's 14877-14878 lab lines lie between J and H, where no
holding exists.

**K, attempted.** PR17 Table 2 measures one K line, 4d 1D2 - 5f 1F3 at 24572.92 A vacuum
(24566.2 A air), log gf -0.391 +/- 0.06. VALD carries three components on the same lower level and
a single, near-degenerate 5f upper level (24565.585 +0.11, 24565.821 -0.87, 24566.296 -0.225),
and PR17 measured only one of them. Its identity among them is ambiguous, so it is not adopted.
Civis+2013 K lines carry no per-line gf error.

## Code guards (all elements)

* `element_prepare`: a published line with a canonical row inside the resolver tolerance is
  reported `NEAR ... -- join explicitly`. No duplicate row is planned (PR #585 left 14).
* `published_line_add.plan`: refuses a feature whose components span more than one synthesis
  cluster; register its components instead.
* Census: the 31 other resolver-tolerance pairs in canonical_gf (Fe/Cr/Ti/Mn/Ce, mostly near-UV
  VALD) are distinct transitions, not add-path duplicates.

## Not changed, for Ryan

Commits `5566b66d` (blocking.nlte_missing forced to 0) and `f5e723c0`/`cddfd84b` (gf-sigma-less
lines admitted) changed the orchestrator gate for every element. Left as is.

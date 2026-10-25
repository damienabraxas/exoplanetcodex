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
| K | CRIRES K | 24566.30 (PR17 exp, see below) |

Coverage notes: 10811/10914/10965 fall between CRIRES Y-wide (to 10796 A) and J (from 11160 A),
so only the KP atlas serves them. PR17's 14877-14878 lab lines lie between J and H, where no
holding exists.

**K.** PR17 Table 2 measures one K line, 4d 1D2 - 5f 1F3 at 24572.92 A vacuum, log gf
-0.391 +/- 0.06 (experimental; no PR17 theory value). NIST ASD (T5539) identifies the same
transition as 24566.30 A air, grade A, -0.225, which settles the identity among VALD's three
components on that level pair. Adopted at PR17's value (Reference Grade = lab gf). The other
two VALD rows are Kurucz intercombination lines to 5f 3F; the strongest, 24565.585 at +0.110,
would dominate the feature, but the CRIRES+ K spectrum shows only a 9% shoulder there against
22% at 24566.37 -- the 12083 pattern (PR17 refuted VALD's +0.45 at 12083.278). Kept at VALD,
since no measured value exists; read the fit residual. Civis+2013 K lines carry no per-line gf
error.

**Literature (2026-10-24).** Amarsi & Grevesse 2026 (arXiv 2608.23155) recommend
A(Mg) = 7.56 +/- 0.03: A21's Mg I lines plus only the three Mg II lines with the least severe
departures from LTE. They add no IR line. Nandakumar+2023 (IGRINS I, A&A 675 A23) measure Mg
from K 21059.76 / 21060.89 / 21458.87 and drop H 15740 / 15748 as saturated.

**NIST grades for the IR (Codex / Deep pools).** The IR Mg rows were seeded from VALD and never
matched to NIST, so K offered two gf-graded lines. `scripts/rya1250_mg_ir_nist_grades.py`
writes NIST's accuracy class onto 49 Y/J/H/K rows whose gf already IS NIST's value (or whose
dominant cluster component is, the 19723 precedent); gf unchanged; NIST rows frozen in
`data/reference/nist/mg_ir_rya1250.csv`. gf-graded lines: Y/J 29, H 23, K 23 (was 2).

## Code guards (all elements)

* `element_prepare`: a published line with a canonical row inside the resolver tolerance is
  reported `NEAR ... -- join explicitly`. No duplicate row is planned (PR #585 left 14).
* `published_line_add.plan`: refuses a feature whose components span more than one synthesis
  cluster; register its components instead.
* Census: the 31 other resolver-tolerance pairs in canonical_gf (Fe/Cr/Ti/Mn/Ce, mostly near-UV
  VALD) are distinct transitions, not add-path duplicates.

## NLTE level labels for the VALD-built lists (Ryan 2026-10-23: "tag the lines")

The 2026-10-23 run FAILED 11 Mg Gerber NLTE cells in NIR/H: no Mg line in the IR/H/K lists
carried an NLTE level label (`nearuv_linelist` writes 'none'). Al had none there either,
and Fe none in K; only RYA-1094's Fe IR/H labels existed. The VIS/red-optical GESv6 list
ships its own labels (Mg 244/245, Al 84/84).

`pipeline/nlte_label_sidecar.py` + `scripts/rya1250_nlte_label_sidecars.py` (Sirius) write
`nlte_labels.tsv` beside each list for every element with a Gerber deck (Mg, Al, Fe). Two
routes per level, unioned, both endpoints or neither: (energy, J) as RYA-1094, and VALD's
own LS designation from the raw extract. The second is what Mg needs: `atom.mg86d` is
term-resolved above 3p and its energies above 4d sit ~0.02 eV below VALD's, so (energy, J)
alone labels 3 of 51 IR lines.

| list | Mg I | Al I | Fe I |
|---|---|---|---|
| IR | 46 / 51 | 31 / 31 | (RYA-1094: 191 in list) |
| H | 37 / 137 | 75 / 78 | (RYA-1094: 502 in list) |
| K | 67 / 126 | 62 / 66 | 695 / 1369 |

Checks: the routes reproduce RYA-1094's Fe labels (IR 191/191, H 498/502, zero
disagreements); 12083.346 -> 4f3F*, 12083.649 -> 4f1F*. Unlabelled Mg H/K lines go to
n >= 10 levels above the atom's ceiling (9l, 7.458 eV). Mg II: the atom holds only its
ground level, nothing to label. LTE flux with and without the labels is identical
(max |diff| 0.0 on Mg H, Fe K, Al K windows).

Why a sidecar, not the list: `run_matrix` fingerprints list ROWS, so labelling in place would
re-run every LTE cell of every element in those bands. The sidecar is merged at load
(`_load_synth_resources`) and fingerprinted only for `gerber-nlte`/`gerber-mean3d` cells of
an element that has rows there.

Not done: Fe IR/H could gain 12 / 296 lines -- that changes published Fe products (Ryan's
call). The near-UV list is regenerated per checkout and not built in the production clone,
so it has no sidecar.

## Parallel Gerber runs shared one departure file

Run 2's budget assembly died on IAG red-optical `SET-MG_BERGEMANN2017` gerber-nlte (8806 only):
its xi-0.1 leg reported "interpol_modeles_nlte produced no departure file"; re-run alone it gives
7.496. `gerber_nlte.for_node` defaulted to a host-wide `/tmp/rya798_gerber`, where every process
deletes, writes and reads `solar_<El>_coef.dat`; under `--jobs 4` a sibling deleted the file, and
could as easily have been read at its own trial abundance with no error. Exposed: the
abundance-axis 1D decks (Mg, Al). Fe 1D interpolates at one fixed A(X), so a sibling's file is
identical; the <3D> decks read directly. Now one directory per process (`default_workdir`).

## Not changed, for Ryan

Commits `5566b66d` (blocking.nlte_missing forced to 0) and `f5e723c0`/`cddfd84b` (gf-sigma-less
lines admitted) changed the orchestrator gate for every element. Left as is.

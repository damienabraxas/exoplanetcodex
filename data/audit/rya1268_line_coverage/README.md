# RYA-1268: Reference Grade line lists for every remaining element

Scope: Li, Na, K, Ca, P, S, Sc, Ti, V, Cr, Mn, Co, Ni, Cu, Zn, Sr, Y, Zr, Ba, Eu (everything
not done: Fe, C, N, O, Si, Al, Mg). Reference Grade = **the Asplund (AGSS21) line set + any
lab-referenced line** (Ryan 2026-10-21). Codex/Deep (VALD) candidates are listed only for the
band cells where Reference Grade is thin.

## How it was built (rerun in this order)

| step | script | writes |
|---|---|---|
| 1 NIST ASD, 32 species, 3000-25000 Å | `rya1268_pull_nist_species.py` | `data/linelists/primary_gf/nist_asd_<El>_<Ion>_3000_25000.tsv` |
| 2 Asplund tables (Scott+2015 x2, Grevesse+2015 from the authors' LaTeX; AGSS21 Table A.1; Gallagher/Storm/Wang by hand) | `rya1268_transcribe_agss21_cited.py`, `rya1268_transcribe_agss21_tableA1.py`, `../../reference/agss21_cited_sets/extra_solar_sets.csv` | `data/reference/agss21_cited_sets/` |
| 3 lab papers' own tables (VizieR MRTs, free; four more parsed from the papers' free PDFs) | `rya1268_parse_lab_pdfs.py <pdf dir>` -> `transcribed/`, then `rya1268_transcribe_lab_gf.py` | `data/reference/lab_gf_tables/lab_gf_lines.csv` (3,458 lines, A +/- dA or stated dex) |
| 4 VALD solar extracts' gf reference codes | `rya1268_vald_gf_refs.py` | `vald_solar_gf_refs.csv` |
| 5 lab / stated-sigma / NIST gf into canonical_gf | `rya1268_adopt_lab_gf.py` | `canonical_adoption.csv`, `held_large_gf_moves.csv` |
| 6 Asplund set + sigma + canonical join, coverage matrix, line lists | `rya1268_coverage_audit.py` | `agss21_reference_lines.csv`, `coverage_matrix.csv`, `line_lists_by_band.csv` |
| 7 line sets | `python -m pipeline.line_sets --element X` | `data/reference/line_sets/*_agss21_*.csv`, `REGISTRY.csv` |
| 8 orchestrator prepare -> apply (no culls) | `rya1268_apply_sets.py` (`--report` re-runs prepare only) | `apply_sets_summary.json`, `held_set_adoptions.csv`, `data/results/orchestrator/prepare/solar_<El>.*` |
| litscans (AGSS21 Table 2 + Lodders 2025 Table 2) | `rya1268_write_litscans.py` | `data/reference/litscan/<El>.yaml` |

Reruns 5-8 on main. `gf_error_model.csv` gets one row per new gf source (step 8).

## The Asplund set

271 lines, **190 with a per-line gf sigma** (55 before the lab tables; 178 in the first cut of this PR). The sigma comes from,
in this order: the source table itself (Gallagher Ba II, Storm Eu/Y), the lab paper's A +/- dA
for the same gf (Wood 2013/2014, Lawler 2014/2015/2017, Sobeck 2007, Den Hartog 2011/2021), a stated
sentence in the text (Grevesse Zr II 5-10%, Scott [Ca II] 0.03 dex), or NIST's class for the
same gf (|d log gf| <= 0.03). Where AGSS21 replaced the printed gf (Sc, V, Cr II, Co) the lab
gf is used.

- **Eu:** AGSS21 and Lodders 2025 both take Storm+2024 (3D NLTE), one diagnostic line,
  Eu II 6645.104, Lawler+2001 gf (0.12) at the 12% Storm states. The canonical row already
  carries RYA-1075's isotope adjudication with the same gf (0.1198) but no sigma. Finished
  work is not overwritten, so Eu still has no graded line and prepare refuses it. Lawler+2001's
  own table (VizieR J/ApJS/182/51, table 6) gives A = 17.9 +/- 2.1 e6/s for this line, i.e.
  **sigma 0.048 dex**. **Decision for Ryan:** attach that sigma to the RYA-1075 row without
  changing its gf. The other 16 Lawler Eu II lines it covers (3688-7426 Å) are now LAB rows.
  Eu 6437.64 is also an RYA-1075 row and stays unpriced the same way.
- **Y:** two sets are registered. Y_AGSS21 (Grevesse, what Asplund adopts) adopts. Y_STORM2024
  (Lodders' source) is registered for comparison only, because it shares 3 lines with a different gf.
- **Unpriced: 81 lines** -- Ti I 22, P I 8, S I 8, Sc II 8, Zr I 8, Mn I 7, Zn I 5, Cr I 3,
  and 1-2 each for Ca II, Cu I, Y II, Ca I, Sc I, Cr II, Co I, Ni I, Sr II. Their lab
  sources (Berzinsh 1997 P, Biémont 1981 Zr I, Lawler 2013 Ti I tables, the NBS compilations)
  print no free per-line error, and NIST grades them D/E or for a different gf. Den Hartog+2021
  priced 8 of the 9 Smith & Raggett Ca I lines (5%); the Sc HFS fix below priced 4 Sc I.
  They are in the line sets and get the published gf, but they are not graded.

## Reference lines per band

Read each cell as **RG / RG + NIST-graded / clean VALD candidates** (blend ratio < 0.5, not
held). RG = Asplund lines + lab-tier lines. A line is one transition (VALD's hyperfine
components are merged and their gf summed). Full table: `coverage_matrix.csv`.

| species | near-UV | VIS | red-opt | NIR | H | K |
|---|---|---|---|---|---|---|
| Li I | - | 2/3/0 | - | - | - | - |
| Na I | - | 4/31/0 | 0/3/0 | 1/4/6 | 0/2/1 | 0/11/2 |
| K I | 0/4/0 | 0/3/0 | 1/5/0 | 2/7/3 | 0/2/0 | 0/1/1 |
| Ca I | 1/1/17 | 40/47/67 | 3/3/62 | 0/0/61 | 0/0/22 | 0/0/59 |
| Ca II | - | 2/2/3 | 3/6/3 | 0/0/4 | 0/0/3 | 0/0/4 |
| P I | - | - | - | 8/8/1 | 0/0/9 | - |
| S I | - | 3/3/16 | 2/2/25 | 3/3/29 | 0/2/29 | 0/0/35 |
| Sc I | 0/5/0 | 5/59/6 | - | - | - | 0/2/1 |
| Sc II | 0/2/0 | 9/18/8 | - | - | - | - |
| Ti I | 2/24/39 | 42/175/201 | 20/33/55 | 3/7/51 | 7/7/36 | 7/10/13 |
| Ti II | 11/11/17 | 70/71/37 | 5/5/5 | 2/2/2 | 0/0/5 | - |
| V I | 108/108/7 | 353/373/48 | 12/14/0 | 0/0/1 | 0/1/2 | 0/0/3 |
| V II | 31/148/10 | 32/55/23 | - | - | - | - |
| Cr I | 6/7/66 | 123/163/321 | 0/3/39 | 3/3/39 | 0/0/9 | 0/0/14 |
| Cr II | 5/5/40 | 51/51/14 | 0/0/1 | - | - | - |
| Mn I | 0/19/13 | 21/57/76 | 0/0/23 | 2/2/6 | 0/0/20 | 0/0/12 |
| Mn II | 0/9/4 | 0/7/2 | 0/0/5 | 0/0/4 | - | - |
| Co I | 161/161/32 | 165/191/197 | 15/17/28 | 0/0/14 | 0/0/25 | 0/0/4 |
| Co II | 5/5/6 | 0/0/3 | - | - | - | - |
| Ni I | 3/4/39 | 106/109/196 | 11/11/87 | 0/0/43 | 0/0/66 | 0/0/32 |
| Cu I | 0/2/2 | 3/14/4 | 2/2/2 | - | 0/0/5 | - |
| Zn I | - | 3/3/3 | - | 1/1/0 | - | 0/0/2 |
| Sr I | - | 1/17/1 | 1/2/0 | - | - | - |
| Sr II | - | 0/0/1 | - | 3/6/0 | - | - |
| Y I | 0/4/0 | 0/22/6 | - | 0/0/2 | - | - |
| Y II | 0/4/0 | 13/43/1 | - | 0/0/1 | - | - |
| Zr I | 0/0/2 | 7/7/7 | 1/1/1 | - | - | - |
| Zr II | 4/4/19 | 6/6/16 | 0/0/2 | - | - | - |
| Ba I | 0/2/0 | 0/2/1 | 0/2/0 | - | - | - |
| Ba II | - | 4/10/0 | - | - | - | - |
| Eu II | 2/2/0 | 8/8/0 | 4/4/0 | - | - | - |

Verdicts over 216 species x band cells: 45 REFERENCE (>= 3 RG lines), 16 THIN, 155 EMPTY.
VIS is covered for every element. **The IR (NIR/H/K) is the gap**: P I (8 NIR), S I (3 NIR),
Cr I (3 NIR), Sr II (3 NIR) and now **Ti I (3 NIR, 7 H, 7 K -- Blackwell-Whitehead+2006)** reach
3 RG lines; Mn I has 2 NIR (Blackwell-Whitehead+2011). **No other element has an RG line in H
or K.** Those cells are where Codex/Deep candidates come in. They are listed in
`line_lists_by_band.csv` with `pool = vald`, VALD gf and gf-reference code, and the blend ratio.

## "Super evil" lines (Ryan 2026-10-23: super bad dex is not included)

- **200 lines held, not adopted** (`held_large_gf_moves.csv`): a lab table or NIST gf that
  would move the canonical gf by more than 0.5 dex (Co I 73, V I 47, Y I 14, Na I 11, K I 10,
  Zn I 9, ...; up to 4.2 dex, e.g. V I 7936.845 -6.07 vs NIST -1.83, Cr II 6421.699 -4.99 vs
  Lawler -1.72). These keep their old gf, are flagged `gf_conflict_held`, and count in **no**
  pool. Each one is a review item, not a fix.
- **Finished work is never overwritten.** A row another ticket adjudicated (`nist`,
  `adjudicated_kr1968_rya466`, `isotope_rya1075`, `adjudicated_rya468`,
  `flagged_manual_pull_rya354`, ...) is skipped by the lab adoption. A published-set
  adoption that lands on one is held: **9 held** in `held_set_adoptions.csv` (Li I 6707.76/.91
  NIST rows at +0.18 dex vs Wang/Yan, Na I 6154/6160 identical gf, Cu I 5105/5218/5220 KR1968
  rows up to 0.10 dex, Ba II 5853.675 0.013 dex, Eu II 6645.104 identical gf). 0 set adoptions exceed 0.5 dex.
- **32 NIST rows refused as corrupt**: log gf disagrees with NIST's own A-value by > 0.05 dex
  (e.g. Zn I 11054.25: fik = 16 gives log gf +1.68, but its A gives -0.36). **3 refused** because
  the upper level does not match VALD's (a different transition sharing the lower level).
- Blends are not cut. `blend_ratio` (summed depth of other VALD lines within 0.1 Å / this line's
  depth) is reported for the synthesis to resolve.

## Orchestrator state (post-apply prepare)

See `data/results/orchestrator/prepare/solar_<El>.md`.
- canonical_gf: **4,040 rows adopted** (lab table 2,429, NIST A..C 1,545, VALD lab code +
  NIST class 52, published stated sigma 14); 805 move by more than 0.1 dex, all within 0.5.
  Published-set gf adopted by apply: Cr 18, Ti 7, Y 6, Sc 5, Mn 4, Co 4, Sr 3, Na 1, Ca 1. Added 7 missing P I lines (NIR) to canonical_gf and the IR synthesis list.
  Apart from those 7 P rows, the IR list is byte-identical in content.
- **READY to measure: P, Zr.** No UNREVIEWED gf source is left. The only adoptions still
  listed (Li 2, Na 2, Cu 3, Ba 1, Eu 1) are the held prior-adjudication rows above.
- **3 Asplund lines the orchestrator cannot add** (Ca, S stay NOT READY on them): Ca I 6169.3 is
  AGSS21's one fitted feature for 616.90 + 616.95 nm (no single transition, so no VALD row);
  [Ca II] 7323.89 and S I 8694.626 fall in iSpec's own GES synthesis list, which is outside the repo.
- **Culls not applied**: every element lists saturated non-literature cull candidates (V 594,
  Ti 429, Co 405, Cr 261, Ni 202, ...). These are Ryan's yes/no.
- **Zn is not one of the 27 orchestrator targets** (`data/config/elements_master.json`), so
  prepare refuses it. Its set (5 Zn I lines) and litscan exist; adding the target is a decision.
- **Eu now prepares** (18 Lawler+2001 graded lines) but is NOT READY: its one Storm/AGSS21 line
  is the protected RYA-1075 row (sigma decision above), plus 9 saturated near-UV/blue cull candidates.

## Findings

1. **The RYA-1250 merge (#585, `ee13da88`) silently reverted RYA-1232 work.** It dropped
   FE_AGSS21 from SOURCES.yaml, the Fe REGISTRY rows, and 7 gf_error_model rows. All three
   were restored here (`99f517da`). It also reverted code in `pipeline/element_publish.py`
   (reference-pool budgeting, supersede `scope` guard, `legs.reindex`). That code is **not**
   restored here; it belongs on RYA-1232.
2. The lab papers that AGSS21 cites as paywalled (Lawler 2014 V I, 2015 Co I, 2017 Cr II) have
   **free machine-readable tables on CDS VizieR** that include A +/- dA. The arXiv copies carry
   only a stub. Lawler 2019 Sc, 2013 Ti I and 2001 Eu II reach us through VALD's gf codes.
3. A wavelength + lower-EP match alone is not a transition: V I 7936.845 matched an unrelated
   line sharing its lower level (-6.07 -> -1.83). Every adoption now also checks the upper level
   against VALD.
4. The Mac fails `check_result_generators.py` with ~2,400 pre-existing violations, and the
   tracked Al/Fe/Mg/Si prepare reports have no GENERATORS entry on main. The 38 new reports
   are registered.
5. NIST pull: the "FAILED" chunks in `data/output/rya1268/nist_pull.log` are "Result did not
   contain a table". Files are named `<El>_<Ion>` because the macOS file system is
   case-insensitive (`NiI` == `NII`).
6. **Literature sweep, all bands (RYA-1268 follow-up, Ryan: "any literature we can get").**
   Checked against Heiter+2021 (the Gaia-ESO atomic-data paper, local `aa36291-19.pdf`, and its
   full VizieR line list J/A+A/645/A106), which names the preferred lab gf source per element.
   Added, all free:
   - **Ti I IR -- Blackwell-Whitehead+2006** (MNRAS 373, 1603; VALD `BLNP`): 78 lines 0.47-2.4 um,
     0.04-0.09 dex. Scott+2015 ranks Lawler+2013 > Nitz+1998 > BW2006, so BW2006 yields wherever
     VALD carries Lawler/Nitz (38 lines; e.g. 5252.100 is -2.36 Lawler vs -2.63 BW2006).
   - **Mn I IR -- Blackwell-Whitehead+2011** (A&A 525, A44; `B-WPNP`): 20 lines, 0.02-0.05 dex.
   - **Ca I -- Den Hartog+2021** (ApJS 255, 27): 77 lines with an experimental reference
     (Smith & Raggett 1981, Smith 1988, new UW lifetimes), 2200-7326 Å. Its 125 other rows are
     theory (M17) and stay out. No lab Ca I gf exists beyond 7326 Å.
   - **Co II -- Lawler+2018** (ApJS 238, 7): 84 lines 2000-3621 Å, A +/- dA.
   - **Eu II -- Lawler+2001** via VizieR J/ApJS/182/51 table 6: 24 lines, A +/- dA.
   - VALD lab codes added to the audit's lab list: `BLNP`, `B-WPNP`, `BXPNL` (BW2005 Mn I),
     `NWL` (Nitz 1998 Ti I), `ZCBS` (Zerne 1997 S I). `PTP` (Pickering 2001 Ti II) stays out
     (Scott+2015 rejected it).
   Not obtainable free: Pehlivan+2015 Sc I IR (A&A 582, A98; tables online-only, no VizieR),
   Ni II 2026 (arXiv 2609.12726, table is a stub), Lawler+2013 Ti I tables (its gf reach us via
   VALD `LGWSC`). The CNO-era IR papers (Smith+2013 APOGEE, Nandakumar+2024, Elgueta/Lim/
   Chiappino 2026) use astrophysical (solar-fitted) gf, which is circular for a solar abundance.
   They are line-selection evidence only.
7. **Two counting bugs in this PR's first cut, both fixed.** (a) The Sc swap took one hyperfine
   component's gf from VALD instead of the sum, so every Asplund Sc line was set 0.4-0.7 dex too
   low (Sc II 5641.000: -1.513 instead of -0.99; the printed value is -1.131). (b) The band
   census joined a merged hyperfine group to canonical_gf at its strongest component within
   0.01 Å. A canonical row at the centroid was then missed, which undercounted every strongly
   split species (Co I VIS 61 -> 165, V I VIS 187 -> 353, Eu II 1 -> 14 RG lines). Lines now sit
   at the gf-weighted centroid and join within their own component spread.

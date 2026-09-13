# RYA-1218 — preflight checkpoint, not a completed Si campaign

Baseline: `origin/main` at `c64eccfe` (2026-09-13). RYA-1169 is merged
(PR #495); its freeze gate remains BLOCKED. This checkpoint does not establish
any abundance, grade pool, ionization balance, correction, or total uncertainty.
No Si authoritative product feed exists at this baseline (`data/products/solar`
contains Fe only). Historical Si abundances remain diagnostic history.

## Reproduce

```sh
python scripts/build_si_intake_rya1169.py
python scripts/run_si_asplund_ew_rya1169.py --out data/results/rya1218/source_set_ew_diagnostic
python -m pipeline.si_protocol_audit
python -m pytest tests/test_si_intake_rya1169.py tests/test_si_protocol_audit_rya1218.py tests/test_si_asplund_grading_rya1169.py tests/test_si_asplund_ew_run_rya1169.py
```

`summary.json` pins source-file SHA256s. The census contains 2,251 Si I and
222 Si II canonical rows. The external replication set retains nine used Si I
rows, one used Si II row, and two rejected Si I rows. All 20 DH23 canonical gaps
remain under the existing species/wavelength/EP candidate join. These are
canonical gaps, not proof of absent observations or unusable transitions.

`line_holding_eligibility.csv` enumerates 32,409 cells: all canonical Si rows
plus the 20 unmatched DH23 rows, across all 13 registered Solar holdings.
It deliberately distinguishes declared spans from untested pixel coverage.
`solar_holdings.csv` preserves raw/corrected siblings and registry conditioning
claims. Nineteen cells link newly served corrected-atlas EW windows; all other cells retain explicit untested coverage.
Grades are separate columns in `canonical_si_census.csv`; unresolved grade
construction is explicit. `historical_grade_delta.csv` preserves external
membership without promoting it to Reference Grade.

`historical_engine_reach.csv` retains the existing Si-specific reach evidence.

The expanded holding run is recorded under
`data/results/rya1218/si_multiholding_ew/`. It measures the ten RYA-1169 source
lines on every registered solar holding: corrected Kitt Peak **9/10** (the
7226.2079 A H2O window is still refused), Kurucz 2005 residual **10/10**,
HARPS raw and Molecfit siblings **8/10 each** (both stop at 6910 A), and IAG
Baker **10/10**. The IAG Reiners blue sibling has no source lines in its
4047--5001 A span. CRIRES+ Y/H holdings have no RYA-1169 source lines because
the source set is optical; their NIR Si line pool must be constructed
separately before a CRIRES abundance route is valid.

For the NIR request, `si_nir_line_pool.csv` records the four published Si I
J-band lines from Bergemann et al. (2013, ApJ 764, 115, Table 1): 11984.20,
11991.57, 12031.50, and 12103.54 A. The CRIRES audit is in
`data/results/rya1218/si_crires_nir/`. None of the corrected solar holdings
reaches J: the available corrected products are Y (9800--10796 A) and H
(15007--17494 A). The raw Vesta IDP does reach these wavelengths, but the
registry correctly refuses it because telluric correction and rest-frame
conditioning are absent. Thus all four CRIRES J lines remain HOLD_MEASUREMENT;
no EW or abundance is reported from that raw holding.

The Kitt Peak result is intentionally split: the Molecfit sibling uses the
RYA-940 product, and its 7160--7340 A H2O band has **no admissible correction**
(the Fe work did not produce a valid corrected H2O product there). Therefore
7226.2079 A remains refused on that holding. The Kurucz residual holding is
telluric-corrected at source and serves that line.

Engine applicability remains separate from observed EW reach. The current
coverage ledger reports Si I Engine A served in VIS and reachable-but-not-extracted
in red-optical; Si I Engine B served in VIS/red-optical; Si II Engine B
uncovered in VIS but served in red-optical; and no Si II Engine A route.
Unsupported engine cells remain HOLD.
`model_inventory_cells.csv` enumerates the nine roster models against both
ions, holdings, and configured synthesis bands. Every cell is currently HOLD,
not executed. It is an applicability work queue, not a completed engine audit;
grade-specific pools and EW/synthesis execution still need adjudication.

## Source evidence and corrections

Inspected local PDFs in `Documents/Exoplanet Codex/Reference documents`:

- Scott et al. 2015, *The elemental composition of the Sun. I. The intermediate
  mass elements Na to Ca*, DOI [10.1051/0004-6361/201424109](https://doi.org/10.1051/0004-6361/201424109),
  `1405.0279v2.pdf`, section 5.4 and Table 2 notes 7–8.
- Amarsi & Asplund 2017, *The solar silicon abundance based on 3D non-LTE
  calculations*, DOI [10.1093/mnras/stw2445](https://doi.org/10.1093/mnras/stw2445),
  `1609.07283v1.pdf`, Tables 1 and 2. RYA-725's old A&A citation is incorrect;
  the existing RYA-1169 MNRAS DOI is the correct lineage.
- Den Hartog et al. 2023, *Atomic Transition Probabilities for Transitions of
  Si I and Si II and the Silicon Abundances of Several Very Metal-poor Stars*,
  DOI [10.3847/1538-4365/acb642](https://doi.org/10.3847/1538-4365/acb642),
  `2301.11391v1.pdf`, section 3.2 and Tables 4–5.
- Kelleher & Podobedova 2008, *Atomic Transition Probabilities of Silicon*,
  DOI [10.1063/1.2734566](https://doi.org/10.1063/1.2734566),
  `jpcrd3720081501p.pdf`, Table 8 and the Si II wavelength/energy tables.
- Pehlivan Rhodin et al. 2024, A&A 682, A184, `aa45686-22.pdf`, modern
  experimentally validated Si I/Si II transition data.
- Bergemann et al. 2013, ApJ 764, 115, `1212.2649v1.pdf`, NLTE Si J-band
  application and atom-model context.

The Si II 6371.370 A source is Scott's mean of Schulz-Gulde (1969), Blanco et al.
(1995), and Matheron et al. (2001), with a reported 0.02 dex uncertainty on that
mean. The old builder incorrectly assigned it the Si I Garz/lifetime provenance.
The generator and generated reference table are corrected; log gf stays -0.044.
That reported uncertainty is not a complete product uncertainty and is not
silently assigned to the Si I rows.

**DH23 lower-level identity is now resolved:** Table 4 identifies 2334.407 A
as the Si II ground-term J=1/2 transition and 2350.172 A as the J=3/2
transition. NIST gives the corresponding lower energies as 0 and 287.24 cm-1;
the intake now records 0.0000 and 0.0356 eV. Table 5 rounds both values to
0.036 eV, so the rounded table value is not used for physical matching.

## Exact source requests / stopped evidence branches

Correction after a file-signature/content scan: Garz is held as the extensionless
`1973A&A....26..471G` (seven scanned pages; title page visually verified).
It is no longer a download request. See `reference_folder_review.md` and the
96-file PDF inventory for the broader applicable holdings.

The following other primary papers were not located in the checked library;
their references are verified in Scott's bibliography or the project bibliography.
Supply the papers and their transition/lifetime tables for a primary-source audit:

- O'Brian & Lawler 1991a, Physics Letters A **152**, 407,
  [DOI 10.1016/0375-9601(91)90834-U](https://doi.org/10.1016/0375-9601(91)90834-U).
- O'Brian & Lawler 1991b, Physical Review A **44**, 7134,
  [DOI 10.1103/PhysRevA.44.7134](https://doi.org/10.1103/PhysRevA.44.7134).
- Schulz-Gulde 1969, JQSRT **9**, 13; Si I/Si II oscillator strengths.
- Blanco et al. 1995, Physica Scripta **52**, 628; Si II experimental data.
- Matheron et al. 2001, JQSRT **69**, 535; Si II experimental data.

These are the remaining specific primary-lineage gaps; they do not block
use of the now-held NIST compilation, Pehlivan Rhodin 2024 data, or Bergemann
2013 model context. AGSS21, Elgueta, later theoretical Si work, historical Si
NLTE papers, NIST ASD row queries, and VALD physical level crossmatches still
require a complete dossier. The NIST compilation is held locally as an
evaluated cross-check and has not been substituted for experimental primary
evidence.

The standing bibliography audit fails on pre-existing unindexed documents,
including duplicate-named copies of the Si papers. No library files were deleted
or renamed. Global bibliography cleanup is not claimed complete.

## Framework findings and remaining gates

| Class | Finding | Disposition |
|---|---|---|
| PIPELINE DEFECT | Si II inherited Si I gf source in RYA-1169 builder | Fixed with source-based regression test |
| ELEMENT-SPECIFIC | DH23 Table 5 rounds the two Si II lower energies to the same 0.036 eV | Resolved using Table 4 identities and NIST energies; intake now stores 0.0000/0.0356 eV |
| PIPELINE DEFECT | Canonical schema has no lower/upper level or J fields | Preserve candidate matches; full physical identity remains HOLD |
| DOCUMENTATION GAP | Hardened phases in RYA-1218 are absent from current ELEMENT_PROTOCOL | Ticket governs this campaign; protocol synchronization remains owed |
| DOCUMENTATION GAP | RYA-725 cites the Amarsi Si paper as A&A and contains approximate historical identities | Use verified source tables and DOI above |
| PIPELINE DEFECT | Importing measure_band_ew performs atlas availability checks even for metadata inventory | Audit reads literal declarations without importing driver; shared metadata extraction remains owed |
| DOCUMENTATION GAP | Si modern lab intake exists but Si is absent from gf_grades.LAB_TABLES | Wiring requires source/identity adjudication first |

Phase 0 is **incomplete**: current registers, canonical data, holdings, models,
Si intake and historical reach were inspected, but runtime/grid/source readiness
has not been fully checked on Sirius. Phases 1–4 are partial evidence inventories.
Phase 5 has only a source-set EW diagnostic rerun; abundance inversion and
Phases 6–10 have not run. The unmodified RYA-1169 runner was executed into
`data/results/rya1218/source_set_ew_diagnostic` after its focused tests passed.
It reproduces 9/10 served lines on the corrected 1984 composite and 10/10 on
Kurucz 2005. The 7226.2079 A corrected-1984 window again raises a missing
correction error; no raw fallback occurs. Si II is served on both holdings.
The output's `ticket: RYA-1169` identifies the reused runner, not this campaign;
the generator registration records the RYA-1218 invocation and output location.
The 13-holding matrix is not an assertion that every holding is science-ready.

No current Si I/Si II balance, gf zero-point cap, LTE→NLTE/3D shift, product
uncertainty, abundance feed, appendix, or PDF was produced. No measurement gate
is signed off. The next executable work is source/identity adjudication and
exact-holding pixel/conditioning checks, followed by independently constructed
grade pools and validated per-ion model routes. Missing primary sources stop
their evidence branches, not unrelated inventory work.

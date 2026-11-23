# RYA-1397 Solar patch campaign

Phase 1 started from merged main at `4f1e69e91d7eee1ae99f034a3f9a7271b262669f`, on 2026-11-23. The [baseline census](solar_product_baseline.csv) joins the canonical element list to committed Solar product files and records SHA-256, version, profile count, bands, treatments and routes. There are 12 element files, with Fe I and Fe II accounted for separately. Rows without a committed file remain visible for reconciliation against active runs and unpublished products. The [input snapshot](baseline_inputs.json) identifies the literature, citation, line-set, stellar and model inputs used at intake.

## Oxygen initial comparison

[Table 2 of Pereira et al. (2009)](https://doi.org/10.1051/0004-6361/200912840), inspected in the local PDF, supplies five oxygen entries. The [comparison](oxygen_pereira_table2.csv) cross-references candidate identities from the [RYA-1353 census](../rya1353_line_census/v4/canonical_census.csv.gz) using species, wavelength, excitation and diagnostic grouping. The source table labels wavelengths in nm, and the census explicitly uses air Å. Fine-structure term and J verification remains open. The four matching oscillator strengths support retaining the existing entries for this source. Section 4.3.1 identifies three unresolved components at 615.8149, 615.8172 and 615.8187 nm. The Table 2 oscillator strength and the nearby VALD row still require component-level verification. Exclude this feature from the initial patch while evaluating other oxygen lines independently. The table's NIST-sourced triplet values are corroboration, and the paper's center-to-limb measurements provide observational validation evidence.

These are source-level decisions. Oxygen still needs the remaining literature, blend, model applicability and observational review before element-level patch adjudication. The [RYA-1352 governance dossier](../../../docs/audit/rya1352/README.md) describes Theoretical Codex as a held proposal. Preserve existing product grades and keep theoretical evidence separate pending scientific review. No patch runs have been launched.

## Batch line inventory

The [line review](solar_line_review.csv) covers all 13 product targets across the 12 committed elements. It contains 9,762 records, including 2,506 baseline indicator records, source-specific alternatives and held candidates. A row is an evidence record, so multiple rows can refer to the same transition. The [coverage counts](line_review_coverage.csv) and [element decisions](element_review_gaps.csv) describe what is available and what still prevents a run-ready list. All 13 targets are deferred at the element level. A wavelength and lower-energy match is provisional until the full level and component identity is verified. Product inclusion is based on the IDs supplied by the baseline product records, plus nine recovered Fe I VIS pools. Those nine pools contain 515 accepted line records, counted with repetition across products. Fourteen Fe I VIS summaries still need per-line reconciliation.

The newly examined sources include the full Fe II laboratory Table 6, calculated Al and Si transition tables, experimental Si I entries and the published source rows transcribed in [inspected_source_rows.csv](inspected_source_rows.csv). The [Den Hartog et al. source](https://par.nsf.gov/servlets/purl/10168066) reports 11 near-UV solar exclusions at 3170–3227 Å. Nine rows of the held laboratory table fall in that interval, so the nine are flagged by interval pending reconciliation with the paper's description. The Biémont and Wood V II mirrors were compared by both level energies and J. Eleven Biémont rows lacked a Wood match. Ten are in the current VIS band and one is outside the current band policy. Confirm the mirror identities against the primary table before adoption.

Calculated rows are held separately as theory candidates. Gauge differences and published fractional diagnostics require their own interpretation. They carry no abundance-line eligibility and the original classifications remain in the baseline files. Molecular rows retain their published vacuum wavelength and require quantum-identity and window matching.

To regenerate from this worktree and the verified local source files:

```sh
python3 scripts/build_rya1397_line_review.py \
  --reference-documents "$HOME/Documents/Exoplanet Codex/Reference documents"
```

The acquired Fe, Si and Al supplements are checked against RYA-1396 hashes before processing. The local mirror and other source artifacts are pinned in [line_review_inputs.json](line_review_inputs.json). Review redistribution terms before publishing the combined export.

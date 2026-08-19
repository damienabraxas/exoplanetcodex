# Run bug ledger

`run_bug_ledger.csv` is the canonical committed history of concrete engineering or
science-pipeline defects encountered during ticket execution and review. It is global
and filtered by `ticket_id`, not split into independent per-ticket truths. Generated
slices are views only.

It is distinct from `debug/intake/*.jsonl`: intake traces are local, git-ignored,
chronological run evidence and may contain routine INFO/WARN events. A trace event is not
automatically a bug. Promote only a real defect, failure, or scientifically material
anomaly; link `run_id` when available without committing the raw trace.

## Ownership and workflow

The ticket that discovers a defect appends the row and owns its initial evidence. Later
tickets amend that row when root cause, disposition, regression test, or fix SHA becomes
known; do not add a duplicate row for the same defect. Unknown facts stay blank or
`UNKNOWN`, never inferred. `fix_ticket_id` links split work to Linear. Git preserves row
history; Linear preserves mission and discussion.

Every implementation/science session close answers: **Did this ticket discover any new
bugs, silent fallbacks, data/provenance defects, or reproducibility hazards?** If yes,
append/amend the ledger and link follow-up Linear work when work remains. If no, say
`none discovered` in the Linear report; never create an empty placeholder row.

## Schema

- Identity: `ticket_id`, `run_id`, `discovered_at_utc`, `element`, `system`,
  `instrument`, `holding_id`, `product_identity`.
- Classification: `stage`, `severity` (`CRITICAL`, `WARNING`, `SUGGESTION`), and a
  concise `category` such as `config_leak`, `hardcode`, `wrong_holding`, `stale_cache`,
  `missing_asset`, `silent_fallback`, `units_frame`, `line_accounting`, `model_domain`,
  `provenance`, `numerical`, `performance`, or `docs_state`.
- Evidence: `symptom`, `expected_behavior`, `actual_behavior`, `root_cause`,
  `detected_by`, `affected_files`, `affected_artifacts`, `science_impact`.
- Disposition: `status` (`open`, `characterized`, `fixed`, `wont_fix`,
  `expected_behavior`, `duplicate`), `fix_ticket_id`, `regression_test`, `base_sha`,
  `fix_sha`, and `notes_provenance`.

Use semicolon-separated values inside a cell for multiple paths/IDs. Do not store mutable
configuration or duplicate full trace contents here. A row describes the defect and
points to evidence; it is not a second results ledger or a replacement for Linear.

The initial historical seed is the Kitt Peak→HARPS wrong-holding regression documented
by RYA-911/913. Fields not recoverable from committed evidence remain explicitly unknown.

# Science CI replay for PR #630

[The GitHub job](https://github.com/damienabraxas/exoplanetcodex/actions/runs/45503121726) stopped at the 45-minute workflow limit. Pytest reached 66% and reported failures and errors before cancellation. The later guards did not run.

Nine early failing or errored test positions were mapped to the Sirius collection order and replayed with the same CI interpreter. The pre-release worktree `925caac8` carries the original Solar science feed. The GitHub checkout was synthetic merge `85852d9f`, containing the updated V I feed. Both selected runs returned eight failures and one error, with identical test IDs:

- `test_al_intake_rya1132.py::test_manifest_is_reproducible_and_never_derives_abundance`
- `test_atomic_intake_rya1129.py::test_builds_cno_atomic_manifests_without_abundances`
- `test_atomic_intake_rya1129.py::test_ledger_is_live_registry_driven`
- `test_bride_h2_tables_rya1281.py::test_published_knots_units_and_strict_domain`
- `test_bride_h2_tables_rya1281.py::test_network_independent_root_and_molecular_limit`
- `test_bride_h2_tables_rya1281.py::test_source_tampering_is_refused`
- `test_bride_native_lte_rya1281.py::test_native_saha_boltzmann_and_debye_lowering` (setup error)
- `test_bride_thermal_sweep_rya1281.py::test_analysis_rejects_changed_fields_and_diagnostic_column_even_with_new_digest`
- `test_chiappino_digest_rya1059.py::test_reproducible`

The replay establishes that these nine failures predate the proposed release. The remaining full-suite failures were not individually adjudicated. The required check remains non-green.

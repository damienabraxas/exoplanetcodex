# Orchestrator run reports (RYA-1222)

Written by `pipeline/run_matrix.py`, driven by
`python run_pipeline.py --star <id> --element <X>`.

| file | what it is |
|---|---|
| `<star>_<El>_latest.json` | the most recent run's report. Every cell of the (band x holding x ion x engine) matrix appears with a terminal status and a reason. A cell missing from here is a build defect, not an empty. |
| `<star>_<El>_<UTC>.json` | the same report, timestamped. Gitignored -- scratch. |
| `inputs_hashes.json` | the idempotency ledger: what each completed cell was built from. |

## Why the hash is here and not on the product

The RYA-1222 spec asked for `inputs_hash` to be stored on the product in
`data/products/<star>/<El>.json`. It cannot be. `publish_product.write_feed`
calls `uncertainty_contract.assert_publication_feed`, which re-validates every
product that is not byte-identical to the row already on disk. The committed
products are legacy rows with no `uncertainty` block and pass only via that
retention clause, so adding one bookkeeping field to one of them gets the whole
feed refused:

```
RYA-587 refuses incomplete live products:
Fe|II|VIS|harps|solar_harps_molecfit_corrected|DEEPGRADED|...: star: provenance is required
```

Writing it here instead keeps the guarantee (same key, same comparison, tracked
in git) and keeps this layer's writes out of the science feed entirely.
`tests/test_rya1222_run_matrix.py::test_the_orchestrator_never_writes_into_the_science_feed`
pins that measurement, so whoever changes the contract finds out this depended on it.

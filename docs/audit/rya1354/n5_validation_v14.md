# N5 loading-control replay and independent terminal validator

The n5 solver is still running; no n5 state or rate is qualified. A read-only
native log prefix captured at 2026-11-12T03:58:08.753720Z contains **24 complete
ground cycles**, with maximum residual **7.58e-6** in the last one, above original
ACCY **5.960464477539063e-8**. A later live observation reached cycle26 with
maximum residual4.48e-6. Neither is a terminal result. The native solver has CPU
activity, not a stalled input prompt; the all-fixed loading control completed.

The lossless `n5_ground_snapshot_v14` archive contains22 native paths in one
424774-byte part: captured log prefix, completed same-input loading control and
the stage's immutable inputs. `scripts/rya1354_extension_progress.py` independently
replays original hashes/history, verifies loading controls and H/N, all25 loaded
orbitals, source n4 wave/initialization controls, and target n5 CSF/nucleus. It
does not compare an unfinished optimization wave against fixed-record acceptance
or admit a state. Machine report: `source_curation_v2/n5_progress_v14.json`.

`scripts/rya1354_extension_audit.py` prepares independent terminal replay. It
requires a terminal n5 frozen-model receipt and matching source/target receipt
hashes, revalidates the n4 model, then checks all6 native n5 stages for strict
residuals, physical labels, exact controls, CSF/nucleus, adopted native outputs,
25-orbital inventory, all16 earlier-layer records and actual initializer/fresh
restart identity. It verifies odd-state ordering, transition inputs and operators,
and the complete13-channel decay inventory before returning exploratory rates.
It has not yet replayed a successful n5 campaign because none is available.

The focused118-test run initially reported116passed and2failed: two new test
expectations incorrectly assumed the earlier23-cycle observation rather than the
actual24-cycle capture, and expected a residual failure instead of the earlier
missing-final-table refusal. These test expectations were corrected to match the
preserved native evidence, without changing parser or scientific criteria. The
affected3-test suite then passed in1.40s. It replays the immutable checkpoint,
rejects live logs as a closed state and rejects replaced native loading inputs.

Production/reference/catalog/star comparisons against original task base and
whitespace checks pass; new archive files are below5MiB. The standing bibliography
audit still exits1 on pre-existing civis2013/mallinson2024 verification vocabulary
before library traversal. No library-success or remote-CI claim.

## Resume after the running campaign terminates

Read the existing `evidence/frozen_collection_v13/status.json` on Sirius. Do not
relaunch the native driver or collector documented in `frozen_extension_v13.md`.
If terminal failure is archived, retain and audit that failure; do not call the
successful terminal validator or use its unqualified wave as a new seed. If the
frozen-model completion status is archived, retrieve all parts to a new local
bundle, then run:

```sh
/Users/ryanschmitt/codex/venv/bin/python scripts/rya1354_extension_audit.py --bundle <retrieved-bundle-directory> --out <new-audit-json>
```

Independently review the resulting six-stage and13-channel evidence before any
n5 comparison. A passed exploratory model still does not establish full core
relaxation, accuracy, calibrated uncertainty or a grade. Earlier full-release
failure and historical n3 component constraint limitation remain preserved.
EXPOSED_NONBLIND, scientific HOLD, adoption false; accepted lifetime, uncertainty,
covariance and calibrated sigma null. No production or scientific state change.

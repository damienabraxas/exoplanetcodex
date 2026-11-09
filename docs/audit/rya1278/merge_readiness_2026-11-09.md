# RYA-1278 / RYA-1281 integration review

Ryan explicitly authorized merging completed work on2026-11-09. This integrates
published foundation PR597/head4635c742 and experimental oxygen PR598/head319819b3
with current main078cf7de. The original active Bride checkout, unpublished2429a892
OpenMULTI3D commit and uncommitted files are excluded and untouched.

The isolated integration branch contains both published heads and current main.
A separate integration PR avoids updating or rewriting the actively developed
original source branch. Original PR598 will target main so both original heads
can be recognized as integrated. Branches are retained for ongoing work.

No physical oxygen-model,1.5D/full3D, abundance, opacity, continuum or theoretical
grade adoption is implied. RYA-1281 remains In Progress; its current physical
qualification tasks and evidence remain visible in the existing reviews.

## Conflict disposition

Four files conflicted. Append-only run-bug records, bibliography rows and artifact
generator declarations from both histories are retained. Register source rows merge
without conflicting scientific values. Both histories independently used v173-v180:
canonical main v173-v180 retain their identifiers; former Bride v173-v243 are
renumbered v181-v251, preserving every narrative. New v252 journals this integration.
Historical ticket comments retain their original branch-local version references;
this explicit mapping resolves their meaning. SEQUENCE receives the paired entry.
Canonical target mirrors check equal to merged constants. Production products,
canonical/reference line lists and stellar config/constants equal current main.

## Validation

Validation is performed on the isolated tree after integrating current main, with
clean current-main baselines for failure-set comparison. Mac held LW0.16.1 supplies
compiled fixture/constructor controls; Sirius validation uses an isolated environment
without changing production solver installations. No production RT/SE or abundance
run is authorized or executed by this merge review. Scientific acceptance remains
HOLD. Exact results and merge receipts are appended after checks complete.


The first review found one new test-collection defect: the FeII-width test imported
an optional LW backend before checking availability, aborting the standard project
suite without LW. Its module now explicitly skips when LW is absent; compiled
width assertions remain unchanged and run in the held LW environment. This is
recorded in the run-bug ledger. Initial disposable-checkout attempts also lacked
the host holdings profile, and root-wide pytest collected an archived command-line
smoke script; validation uses the documented CI `pytest tests/` scope and identical
explicit host profile. Held LW lacks pandas, so iSpec/policy post-effects run in the
project interpreter rather than changing the existing held environment.


Linux validation (isolated Sirius directory `rya1278-1281-merge-20261109`):
LW0.16.1, Python3.12.13, NumPy2.2.6, pytest9.2.0, numba0.68.1, scipy1.18.1.
Production environments are unchanged. Fixture/constructor/converter/synthetic
RT/SE checks:137passed/11skipped in101.45s. All11skips came from a pre-existing
Mac-only vendor path despite the registered Sirius6.4MiB Fe atom being present.
The test now resolves `grids.gerber_ts`; those11 checks pass in6.98s against the
existing read-only deck, with no copied alias/download. This coverage defect is
also recorded in the bug ledger. Vendor SHA256:
53446fe1db5e22e35ad972c3d36b990d27bcf6b13c92c95c00aa57cf1f6d7243.
Native fixture runtime code is unchanged from8a0cd4fa; canonical test-path fix is
0a452a18. The isolated environment initially lacked pytest, and a compact source
subset lacked a helper script; final tests use the complete transferred snapshot.
No failed setup is presented as a completed check.

Mac compiled controls:398passed in53.04s using the existing held LW0.16.1
interpreter and explicit host profile. Three iSpec/policy checks pass in5.12s in
the project interpreter (plus an explicit LW-module skip). After the path/guard
fixes, both changed compiled test files pass14checks in5.39s in the held backend.
The generic CI-style project suite is run separately against clean current main
and the integrated tree, with the same profile and interpreter.


Full integration regression exposed a critical catalog-type defect: three model/
source assets had been appended to the observational holdings registry with blank
instruments and invalid telluric states. Instrument and matrix setup then failed.
All three original records are preserved byte-for-byte in
`data/audit/rya1278/merge_review_20261109/non_observational_source_assets.csv`;
the observational registry now equals current main exactly. Existing path-register
source discovery remains intact, and no pipeline/script/test consumes these IDs as
observations. No fake instrument or telluric status is assigned. This is corrected
before merge and recorded as a critical wrong-holding defect.

The same full run identified40 backend-absence failures beyond the collection
issue. Only28 native test functions across13 mixed test modules now carry an
explicit optional-LW requirement; their pure mathematical/metadata controls still
run without LW. The compiled numerical assertions remain unchanged. Both the
standard project controls and complete held-LW controls are rerun after correction.
Initial full-run failures remain in the validation receipt and are never described
as a passing suite. The final effective failure-set comparison incorporates the
recorded corrected-case results against clean current main.


## Final comparison and merge scope

Clean current main:43failed/4697passed/51skipped/7xfailed in2088.97s.
Initial integration:122failures/4errors/4873passed (the receipt retains all results).
After the catalog/dependency corrections, complete affected-file reruns give
391passed/68skipped/3failures in724.55s; all three failures are baseline failures.
The remaining catalog regression group gives77passed in13.25s. Repeated compiled
LW controls give398passed in52.86s. Linux catalog/telluric contracts give16passed
in5.54s. All numerical test function bodies are unchanged by dependency markers.

Joining the full initial run with the latest executed corrected-case dispositions
gives4916passed/114skipped/7xfailed and the exact same43failures as clean main:
**zero added or removed failures,219 additional passes**. This is a full-suite plus
targeted-correction comparison, not a second fresh full-suite run or a global-green
claim. Details, node identities and receipt hashes:
`data/audit/rya1278/merge_review_20261109/validation.json`.

Four derived pytest fixture outputs in this disposable checkout were restored
byte-for-byte from its committed tree after testing; none is staged as a science
change. Production products, canonical/reference lines, stellar config/constants
and the observational holdings registry equal current main. Source-asset records
are preserved separately. Register/target/model checks pass. Global pre-existing
failures remain visible; no science threshold is relaxed.

While validation ran, the original PR598 advanced from reviewed319819b3 to
b81cfd5d with three explicitly RYA-1371 commits (2429a892,d52b74a9,b81cfd5d).
Those newer reference-engine/EOS changes are not part of this scoped review or
integration. PR598 now tracks active1281/1371 continuation; merging the completed
checkpoint through PR619 leaves only that newer work outstanding. Original source
branches and the actively edited checkout are retained. Only1278/1281 completed
work is authorized for this integration; other branches are untouched.

Push uses the documented draft --no-verify pre-push exception for unrelated merged
sequence entries; no file-size-hook bypass. Hosted CI is not asserted passed.

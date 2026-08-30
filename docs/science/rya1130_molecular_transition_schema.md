# RYA-1130 — dedicated canonical molecular-transition schema

## Decision

Molecular transitions live in `codex.molecular_transition/1`, independently of
atomic `canonical_gf.csv`. They do not receive Codex Grade or Deep Grade: those
grades require a primary-laboratory atomic gf and measured Solar feature depth.
Asplund reference-set membership is also a separate downstream join and is not
inferred from presence in a molecular list.

The machine contract is
[`data/schemas/canonical_molecular_transition.schema.json`](../../data/schemas/canonical_molecular_transition.schema.json).
It represents, independently:

- molecule and isotopologue;
- raw source transition identity plus typed electronic/vibrational/rotational
  identity where a release-specific decoder can do so unambiguously;
- wavelength value, medium, frame, and conversion provenance;
- transition probability and exact source identity;
- database/release/table, DOI/ADS/catalog identifiers, retrieval date, license,
  and source-file checksum;
- partition-function and dissociation-energy sources;
- component and isotopologue normalization;
- precedence/conflict decisions, band reach, and intake status.

Null provenance means “not yet established”, never “not applicable”. A source
label such as `VALD3`, `CH PGopher`, or `Sneden web` is not promoted into an exact
release, DOI, license, thermochemistry source, or normalization convention.

## Identity and migration

`pipeline/molecular_canonical.py` streams Turbospectrum `.bsyn`/`.dat` rows into
the schema. The stable `transition_id` is derived from molecule, isotopologue,
source checksum, source line number, and the raw transition label. The source
label is always retained. Only the unambiguous ExoMol CO `v…_J…` convention is
typed today; heterogeneous PGopher/Masseron labels remain lossless raw text until
their exact release-specific grammar is ratified.

Some legacy C2 records contain an empty quoted transition label. Those rows are
retained as `SOURCE_LABEL_NOT_SUPPLIED`; they are never silently dropped or given
an invented identity. With otherwise complete provenance they route to
`CROSSMATCH_REVIEW`.

Run the audit without writing a transition product:

```bash
python scripts/migrate_molecular_transitions_rya1130.py \
  --report data/results/rya1130/molecular_migration_audit.json
```

An incomplete audit exits 2. JSONL output is refused while any row is blocked
unless `--allow-blocked` is explicitly supplied; this produces a review artifact,
not a frozen-ready product.

## Current gate

The initial audit traverses every vendored CH, C2, CN, NH, OH, CO, and newly held
rovibrational source row. It proves they are representable, while also proving the
v1 manifest is not provenance-complete enough to freeze. The missing typed fields
are recorded per source family in
`data/results/rya1130/molecular_migration_audit.json`. RYA-1131 source work must
fill those fields from the actual releases; this schema prevents a delivery label
from hiding that debt.

No abundance is derived and no oscillator strength is tuned by this work.

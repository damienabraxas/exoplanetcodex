# RYA-1129 Wave 1 — C/N/O atomic-data intake

## Gate verdict

**C, N, and O remain `CROSSMATCH_REVIEW`; none is
`FROZEN_READY_FOR_MEASUREMENT`.**  This pass inventories the live atomic store
without changing a gf or deriving an abundance.  Molecular intake stopped at a
real schema limitation and is now RYA-1130.

The authoritative `data/config/elements_master.json` currently declares **27
target rows / 26 unique atomic symbols**, not 28 elements.  It contains no Zn.
That conflicts with RYA-757-era artifacts describing Zn as the 28th canonical
element.  The generated ledger follows the registry and does not invent a Zn
row; the discrepancy must be reconciled at the registry source.

## Atomic inventory

Counts below are generated from `data/linelists/canonical_gf.csv`.  “Fallback”
means the existing `gf_tier` is `KURUCZ` or `VALD3`; “evaluated” means an
explicit NIST grade/tier is stored.  `OTHER` is not promoted to laboratory
quality merely because a short source code exists.

| element | required atomic species present | rows | primary lab | evaluated/NIST | fallback | other | explicit source DOI | wavelength range (air A) | gate |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| C | C I–V (pipeline abundance indicator: C I) | 1,788 | 0 | 1 | 1,392 | 395 | 0 | 4209.717–24443.430 | review |
| N | N I–V (pipeline abundance indicator: N I) | 474 | 0 | 0 | 14 | 460 | 0 | 3466.497–15682.875 | review |
| O | O I–VI (pipeline abundance indicators: O I and [O I]) | 559 | 0 | 4 | 158 | 397 | 0 | 3692.395–21330.817 | review |

The current pipeline requirements come from
`config/physics_regime_rya400.yaml`: C I 5052/5380 plus IR C I; N I red
multiplets 7468/8216/8683; O I 777 plus [O I] 6300.  Higher ions are retained in
the manifest because this is an atomic-store inventory, not a Solar line cut.

### Carbon

The dominant stored source codes are K10 (799), NIST10 (326), K14 (315), and
VALD3 (278).  No row is tagged `LAB`; one row carries an evaluated NIST grade,
and no C atomic row has `gf_source_doi`.  The repository therefore has broad C
I–V coverage but not a provenance-complete frozen C atomic foundation.

Required molecular indicators are CH, C2, CN, and CO.  The current store has
CH 2,265, C2 2,654, and CN 2,766 rows, but the sampled/current representation is
VALD3-only with no source DOI.  CO has zero canonical rows.

### Nitrogen

The dominant stored source codes are NIST10 (333) and KP (110), followed by
VALD3 (14) and FTa (13).  No row is tagged `LAB`, no explicit NIST grade is
stored, and no N atomic row has `gf_source_doi`.  N I is present, including the
red/near-IR reach required by the current pipeline, but its source chain is not
frozen to the RYA-1129 standard.

Required molecular indicators are CN and NH.  NH has 114 canonical rows; like
CN, it lacks the typed molecular provenance needed to distinguish transition
data, isotopologue/component normalization, partition functions, and
dissociation-energy inputs.

### Oxygen

The dominant stored source codes are NIST10 (193), KP (164), K11 (129), VALD3
(29), and WSG (28).  Four rows carry evaluated NIST grades, including the
Storey & Zeippen forbidden-line provenance already represented in the store;
none of the O atomic rows has `gf_source_doi`.  The [O I] 6300 blend also
depends on the separately canonical Ni I blend gf and cannot be frozen by an O
row in isolation.

Required molecular indicators are OH and CO.  The canonical store contains one
OH row and no CO rows, so molecular O coverage is structurally incomplete.

## Source verification boundary

No citation was created from memory.  This pass relies only on already verified
repository bibliography records, particularly the NIST C/N compilation context
(`wiese_fuhr2006`), the forbidden oxygen source (`storey_zeippen2000`), VALD3
(`ryabchikova2015` + `vald3_tool`), and the molecular line-list context
(`tennyson2016`).  Those records establish source identity; they do **not**
magically attach a source to individual rows.  Row-level DOI/catalog mapping is
still missing and is why the gate stays closed.

## Schema stop and next action

`canonical_gf` mixes atomic and molecular species but supplies the latter no
typed fields for isotopologue, electronic/vibrational/rotational identity,
line-list release/table, partition function, dissociation energy, license, or
normalization convention.  RYA-1130 now owns that general defect.  Forcing
molecular data into atomic `gf_tier` semantics would erase rather than preserve
provenance.

After RYA-1130, resume Wave 1 with exact row-level source resolution for the
pipeline C I, N I, O I/[O I] indicators, then ingest the independently typed
molecular sources.  Only after zero unresolved crossmatches and complete source
identifiers should these ledger rows advance from `CROSSMATCH_REVIEW`.

## Reproduction

```bash
python3 scripts/build_atomic_intake_rya1129.py
python3 -m pytest -q tests/test_atomic_intake_rya1129.py
```

Outputs are the three atomic manifests, the live-registry-driven intake ledger,
and metadata under `data/audit/rya1129_atomic_intake/`.  The metadata explicitly
asserts `abundances_generated: false`.

# RYA-1134 provisional Al adjudication

This is a reviewable checkpoint, **not a scientifically frozen measurement pool**.
No abundance was measured, adopted, or published. RYA-1134 remains open.

`verified_v2/line_dispositions.csv` gives every frozen candidate a source disposition.
`pool_memberships.csv` preserves all 505 rows for each of Reference (`reference`),
Codex (`our-graded`), Deep (`our-deep-graded`), and external replication (`asplund-al`).
Consumers can call `pipeline.al_grade_verification.load_pool_dispositions(line_set)`;
they must retain HOLD/NONMEMBER rows and separately pass exact-holding eligibility.
The replication set continues to use its own published gf via `reference_lineset`.

## Results

| Historical source class | Laboratory | Mixed laboratory/theory | Evaluated theory | Unresolved |
|---|---:|---:|---:|---:|
| Laboratory (18) | 15 | 3 | 0 | 0 |
| Evaluated (19) | 0 | 0 | 12 | 7 |
| Theory (1) | 0 | 0 | 0 | 1 |
| Fallback (467) | 0 | 0 | 24 | 443 |

Reference has 54 provisional members; Codex 4; Deep 7; external replication 6
plus its explicitly excluded seventh line. These axes overlap. Of 55 source-qualified
rows, 50 pass the atomic handoff checks; four are source-only mid-IR controls without
raw holding components. Johnson Al II is retained from the accessible Johnson
measurement, with the unavailable Träbert successor recorded as a limitation.
451 rows lack a qualified, finite-uncertainty source match in the inspected assets.
This is a bounded source verdict, not proof that no stronger evidence exists anywhere.

The source and candidate ledgers retain alternate evidence and its limitations.
Burheim's three downgraded claims use theoretical 4f lifetimes. The 3092.839 line
does not inherit 3092.710's uncertainty; the 6696.185 neighbor receives no laboratory
gf from 6696.015; the 11254.9 unresolved feature cannot adopt one component's gf.
Uncertainties are **source bounds with their stated confidence semantics**, not
silently converted to Gaussian one-sigma errors. NIST E gets no finite placeholder.

`component_proof.csv` recovers raw VALD term/J identities and upper energies with
file/line references. Numeric component fractions and rescaled sums are proof
artifacts; no production component deck was replaced. See the independent source
review in `data/reference/al_grade_verification_rya1134/source_review.md`.

## Holding and engine gates

The matrix contains 6,565 line × holding cells across all 13 registered Solar
holdings. Coverage is currently a **declared-range audit**: in-range pixels,
order gaps and usable fitting windows have not been independently inspected.
Registry telluric and normalization declarations remain distinct from exact-window
telluric evidence and the still-unestablished observed-conditioning field.

The 5,616 engine/route cells contain 676 HOLD and 4,940 N/A, with reasons in every
cell. They cover the nine registry models, four pools, two routes and every
represented band/holding. Models 5/6 bind to `Al@mean3D`, keeping LTE and NLTE
separate; they are not full 3D. Engine A's actual Al source is the registered
GALAH/PySME-derived additive 1D-NLTE grid. Its transition identity still needs
revalidation before applying corrections to the refreshed pool.

**No scientifically eligible run cell is established.** RYA-1176's product-identity
and conditioning propagation, exact-holding pixel/context validation, and relevant
model/atom/grid validation must pass before these planning cells become RUN.
The matrix is not an assertion that those gates have been executed or passed.
RYA-1217 Gate 0 remains closed.

## Bibliography

The reconciliation inventories 135 top-level local files (93 PDF first-page scope
reviews), including alternate copies and relevant supplements. Browser asset folders
are excluded. It does not claim full-text review of every contextual paper.
Five canonical records were added, four updated, 11 alternate source copies
reconciled, and 37 citation keys checked. Nine cited sources lack a Reference
Documents copy. Träbert et al. 1999, J. Phys. B 32, 537–552,
DOI `10.1088/0953-4075/32/2/031`, remains unavailable for primary uncertainty
review; Johnson 1986 is retained with that limitation.
The two Buurman papers have verified ADS identities; their DOIs remain unestablished.
Exact counts and evidence are in `bibliography_summary.json` and
`bibliography_reconciliation.csv` under the source-review directory.

## Reproduction

Use the scientific environment with NumPy/Pandas/SciPy and a configured Kitt Peak
atlas path (`CODEX_KP_ATLAS`): importing the existing holding-reader registry
requires that staging, although this generator does not measure spectra.

```sh
python -m pipeline.al_grade_verification --out data/results/rya1134/verified_v2 --check
```

To generate again, use a **new** output directory without `--check`; existing
outputs are never overwritten. `input_hashes.json` pins code, catalogs and source
assets. The original 505-row manifest and canonical gf remain unchanged.

Remaining closeout work is explicit: independent
pixel/context audit of the holdings; validated correction/model applicability;
RYA-1176 propagation; and execution/comparison of any subsequently eligible cells.

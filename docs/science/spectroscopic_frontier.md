# The Spectroscopic Frontier — Atomic & Molecular Model Atlas

**Draft inventory for scientific review (RYA-1373).** No new scientific validation,
line grade, atmosphere/model admission or execution permission is conferred. Public
searchable HTML follows approval of the RYA-1370 Roadmap; this PR delivers the
independent data inventory, schema, reproducible matrix and contribution process.

## Canonical ownership

- Editorial scope and source selectors: `config/frontier/atlas.v1.json`.
- Generated atlas: `data/products/frontier/coverage.v1.json`.
- Contract: `data/schemas/spectroscopic_frontier.schema.json`.
- Generated human-readable matrix: `spectroscopic_frontier_matrix.md`.
- Builder/semantic validation: `pipeline/frontier_atlas.py` and
  `scripts/generate_frontier_atlas.py`.

The seed owns atlas descriptions and review dispositions. It points to existing
model/star/element registries, atomic/molecular schemas, correction tables and
bibliography; it does not replace them. Source file hashes and immutable Git URLs
are generated from the recorded science base. The builder rejects a working file
that differs from that pinned Git source. Derived transition and parameter summaries
are reproducible views, not hand-maintained alternate scientific definitions.

## What the seed establishes

O I and C I have distinct held 2019 1D NLTE and 3D NLTE correction tables. These are
abundance corrections from published source calculations, not installed full-3D
transfer solvers. Numeric stellar-axis envelopes are measured from their stored
nodes using the primary CDS fixed-width table contract. A rectangular envelope
is not the valid interpolation hull: it does not establish every cross-combination,
line, chemistry, atmosphere or observed-star application.

The CDS release explicitly labels Fe II table4 **3D LTE versus 1D LTE**. That table
is described upstream but is intentionally not vendored here. The historical RYA-1035
decision text calling its Fe II leg “3D-NLTE” conflicts with the primary release and is
recorded as an actionable atlas gap. No affected Fe II product was changed or
asserted to be wrong. Primary release: Amarsi, Nissen & Skúladóttir (2019),
DOI [10.1051/0004-6361/201936265](https://doi.org/10.1051/0004-6361/201936265),
retained `data/nlte_grids/amarsi2019_cno/ReadMe` and `provenance.json`.

Fe I/II 2022 neural-network resources are **correction surrogates**, not coupled
RT executables. Training-line excitation, energy and wavelength envelopes are derived
from the existing recovered training set. Full stellar, abundance and line-domain
checks remain with the original loader/audit. Neither VIS proximity nor a finite
network value authorizes IR extrapolation. Citation: Amarsi, Liljegren & Nissen
(2022), [10.1051/0004-6361/202244542](https://doi.org/10.1051/0004-6361/202244542).

The measured `atom.fe607a` receipt distinguishes Fe I bound-bound support from Fe II
reservoir levels with no Fe II bound-bound transitions. Its Fe II NLTE application
is **incompatible for this specific atom**, not absent from science generally.
Independent Fe II correction resources remain separate records.

Held CO/OH/CN/NH/CH/C2 list families and additional OH/NH/CH rovibrational releases
retain their manifest source/isotopologue distinctions. A legacy manifest's wavelength
medium, missing CO isotope metadata, uncertain release, partition function,
dissociation energy or normalization are not guessed. Primary molecular transition
identity and source joins remain with the existing RYA-1130/1136 products.
H2O/TiO/FeH/CaH/MgH/SiO records are discovery pointers, explicitly **not yet searched
at release level**, with null literature-audit dates. Resolved database landing pages
are not evidence of a complete search or scientific validation.

## Independent scientific dimensions

Every record distinguishes component type, treatment, local availability, coverage
status, independent verification, license, version and audit scope. `components`
links solver/atmosphere/model-atom records by type; no default mixing is implied.
Atmospheres, model atoms, line lists, departures/corrections and synthesis/SE solvers
cannot substitute for each other. Mean-3D, independent-column 1.5D and coupled
full-3D NLTE remain distinct treatment values. A 3D atmosphere file does not prove
horizontal radiative transport.

A physically useful domain needs Teff, gravity, metallicity, abundance/chemistry,
C/O and molecule formation, activity, magnetic sensitivity, pressure/broadening,
line/atom completeness and matched observational conditioning as applicable. Spectral
class is a coarse discovery label. Unassessed constraints are explicit, and all
inventory records have `execution_approved: false`. The atlas cannot approve a
scientific production run. IR abundance requires the exact holding's telluric
verification under the standing science standards; a line-list wavelength envelope
does not authorize an extension of the current project science window.

Bands are display categories in nanometers: UV below 400, VIS 400 to below 780,
IR from 780. These display boundaries are not applicability cuts. Each grid's
actual transition labels, air/vacuum values and source checks govern interpretation.
O I multiplet-average labels without a single component wavelength are retained in
`labels_without_component_wavelength`, never dropped or assigned an invented line
center. Training identities retain their exact CSV row under the source checksum.

## Coverage vocabulary

| State | Meaning |
|---|---|
| not_yet_searched | Release-level search/audit remains owed. A catalog gap is not a literature absence. |
| not_published | Requires an explicitly scoped completed search and evidence; no such claim is made in this draft. |
| unavailable | A named resource cannot be obtained under documented circumstances, not that no resource exists. |
| incompatible | A specific resource fails the stated physical/interface requirement. |
| unvalidated | A resource or claim has not received the independent scientific review required by this atlas. |
| validated | Reserved for independent sign-off with scope/evidence. This draft validator rejects this promotion. |

Availability is separate: `held`, `published_not_acquired`, `not_checked` and
`unavailable`. Published literature validation, historical numerical controls and
current project adoption are different facts. All current entries retain independent
review as owed; the atlas does not downgrade or overwrite existing adopted products.
Dates record the stated scoped source inspection, not a comprehensive literature audit.
Discovery records retain a null date rather than invent a literature review.

## Citations and software credits

Exact citation metadata is joined from `data/refs/bibliography.csv`; source records
retain file hashes and evidence URLs. The Fe surrogate's retained license is MIT,
copyright 2022 Sofie Liljegren; preserve its copyright/permission notice with source
copies or substantial portions and cite the associated paper. No upstream solver,
network or large database is redistributed by this atlas addition.

Turbospectrum is external software credited to Gerber et al. (2023),
[10.1051/0004-6361/202243673](https://doi.org/10.1051/0004-6361/202243673);
per-install version, atom, departure grid, atmosphere and reuse terms still belong to
its exact run/source records. Published 3D LTE molecular work is distinct from
molecular NLTE: Amarsi et al. (2021),
[10.1051/0004-6361/202141384](https://doi.org/10.1051/0004-6361/202141384).
The ExoMol data index and database paper
[10.1016/j.jms.2016.05.002](https://doi.org/10.1016/j.jms.2016.05.002)
are discovery resources, not proof of local use or release-specific reuse rights.
Exact database/extraction license fields remain pending until reviewed upstream.

## Reproduction and review

```sh
python scripts/generate_frontier_atlas.py
python scripts/generate_frontier_atlas.py --check
python -m unittest tests.test_frontier_atlas
```

The dependency `jsonschema>=4.18,<5` is declared in requirements. Existing CI pytest
discovery includes the new unittest cases. No raw spectra, solver models or protected
reference products are written. See the contribution guide and audit dossier for
review gates, limitations, commands, source conflicts and handoff.

Roadmap integration: the atomic laboratory's model/property scope (RYA-1356),
independent cross-engine research product (RYA-1368), molecular framework (RYA-1314)
and cool-star benchmark work (RYA-751/970) link into the atlas gap registry. The
Roadmap PR's public claims and the atlas scientific inventory need separate review.
A later website integration should consume the reviewed generated JSON, retain the
source/status distinctions above and preserve missingness and null results.

# O777 omitted-opacity source inventory

2026-11-07, RYA1281/1278, scientific HOLD. Worktree `/Users/ryanschmitt/codex/rya1279`,
branch `ryandamienschmitt/rya-1279-1280-bride-oxygen-gates`, stacked base4635c742,
starting evidence a3f476b7. Origin refreshed, main de96ae51 unchanged. No other
worker changes. This follows the tiny fresh-O feedback result, not a new O fit.

## Predeclared design

Plan `config/experiments/rya1281_missing_opacity.json`: inventory the whole held
7763-7783 consumer-AIR-A region, fixed O +/-0.5A windows and prior fixed FeII
shoulders at both fixed RVs. Counts are geometry only, not transition identity,
line strengths, missing-absorption attribution or molecular EOS admission.
Rank remaining atomic rows from both existing frozen LTE screens at the originally
declared logtau500~0depth; exclude only modeled41FeI and already tested gf_079010.
Preserve all other candidates and current-main mismatch flags. No new transfer.

Reuse existing owners, not a parallel molecular schema/parser:
`pipeline/molecular_lists.py`, `pipeline/molecular_canonical.py` and
`scripts/migrate_molecular_transitions_rya1130.py`. Read RYA1130 schema and1136
primary-intake documentation first. Molecular provenance is separate from atomic
canonical_gf; source labels such as Sneden web/PGopher are not full provenance.
Nine routed750-800nm files are byte-identical between canonical vendor and
installed iSpec; verify again at runtime, along with declared/parsed row counts.
Native TS babsma.f/hydropac.f explicitly state bsyn uses AIR wavelengths, and
the wrapper converts input nm to A. Preserve unknown primary-release wavelength
metadata: consumer-coordinate membership does not qualify the primary frame.
Stable source transition IDs come from the canonical owner; do not match transitions
on wavelength alone or invent rotational/component labels.

No gf reweighting, guessed isotopic/partition factors, molecular populations,
pressure widths or absorption installed. Existing BC2016neutral-H2 qualification
does not automatically qualify CN/CH/C2/NH/OH chemical equilibrium, spin conventions,
coupled depletion or line strengths. Existing primary328row CNO work is specific
to its matched abundance-feature census; do not extend its admission by molecule name.

Precommit plan/runner before source evaluation. Exact command: project Mac
interpreter after `source config/hosts/mac.env`,
`python scripts/audit_bride_missing_opacity.py --output <fresh-directory>`.
No Sirius/CI, production/default/product/adoption/abundance/merge/publication.

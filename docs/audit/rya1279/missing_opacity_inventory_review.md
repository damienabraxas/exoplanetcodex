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

## Measured inventory

Code/plan precommitted102d8ff6; inventory command exits0. Nine routed files parse
19504data rows, each exactly matching its declared header count. Canonical vendor
and installed files match byte-for-byte.651source rows fall within7763-7783A in
the TS consumer coordinate; raw labels, stored energies/loggf, isotopologues and
canonical source IDs retained. No decoding or normalization invented.

| Molecule | Rows in20Aregion | O7771window | O7774window | O7775window | FeIIblue shoulder | FeIIred shoulder |
|---|---:|---:|---:|---:|---:|---:|
| CH |26|1|1|3|0|0|
| CN |355|15|20|18|4|5|
| C2 |149|3|2|4|0|0|
| OH |121|5|11|6|0|0|
| NH |0|0|0|0|0|0|

Window/shoulder counts in this table use the fixed plus-gravity convention;
both RVs are in geometry.csv. O apertures remain the predeclared observed
+/-0.5A windows; FeIIshoulders are the prior fixed0.35-0.55Aoffsets about its
fixed observer centre. Rows are not observed-feature identifications or independent
physical transitions; different isotopologues/components can contribute multiple
rows. Counts are not abundance-weighted or strengths. NH0describes these held
lists/region only, not absence of all possible NH opacity. Family/isotopologue
or wavelength coverage does not establish line-list completeness.

All651retain `BLOCKED_PROVENANCE` from the existing canonical molecular owner.
The legacy manifest lacks its typed release/retrieval/license/partition/
dissociation/component/isotopologue fields; primary medium/frame and probability
kind remain unknown. Source labels retained as labels. Native TS expects AIR,
but that consumer behavior does not independently verify each primary source's
frame or authorize changing the canonical metadata to AIR. No molecular opacity,
partition/chemical population or C/N abundance has been calculated here. Existing
primary CNO matched-feature work remains valid in its own scope; these rows have
not been promoted by association. Current production products remain unchanged.

Remaining atomic ranking retains200rows per frozen column (400case/line rows),
all400current-main canonical rows exactly match their frozen source rows. The
earlier changed TiI row remains withheld in the source screen, not reinstated.
Highest canonical natural-LTE coefficient ratio at the reference depth:

| Candidate | AIR A | Lower E(eV) | Stored loggf / reference | Ratio at logtau500~0 | Strongest historical footprint |
|---|---:|---:|---|---:|---|
| V I gf_137515 |7768.902|3.632|-0.306 / K09|0.001857479|7771 IAG |
| Ti II gf_130078 |7776.895|6.504|-3.011 / K10|0.000025644|7775 synthesis |
| Ti I gf_127776 |7768.632|3.702|-1.414 / K10|0.000019289|7771 IAG |
| Fe II gf_078998 |7773.443|9.849|-2.574 / RU|0.000006642|7775 IAG |

These are previous natural-only local extinction/absorption ratios on historical
anchor footprints, not flux, isolated EW, physical line-grade or formation-depth
predictions. The V I same-footprint maximum over all depths is0.04065near the
outermost layer; no claim that this depth dominates the emergent spectrum.
gf accuracy, full component structure and pressure/Stark source qualification
are separate from natural-only eligibility. No V/Ti/Cr line installed here.

## Validation and evidence

Reused canonical molecular parser and primary-intake tests13passed4.83s, exit0:
`python -m pytest -q tests/test_molecular_canonical_rya1130.py tests/test_cno_primary_ingest_rya1136.py`
with project Mac interpreter. Test log `/private/tmp/rya1281-missing-opacity-tests.log`;
inventory log `/private/tmp/rya1281-missing-opacity.log`. Runtime source/count/
numeric checks are additional inventory evidence, not physical abundance tests.

Four result copies under `data/results/rya1281/missing_opacity_v1` byte-match the
immutable registered `experiments.bride_missing_opacity` originals; largest95780bytes.
Actual generator/invocation registered. No source, old output or molecular
canonical product overwritten, no raw spectral/binary copy, no5MiBguard exemption.
Mac one-thread source audit; no new RT/SE, Sirius, CI, installs/downloads or benchmark.
No new source/parser defect discovered; known molecular provenance debt retained.

Fresh library audit:196documents/127accounted/69unindexed/12off-Mac declarations,
exit1. **Correction:** preceding coupled-phase review/Linear text quoted195/68,
but its saved `/private/tmp/rya1281-feii-coupled-library.log` already recorded196/69.
That was a carried-forward documentation error. Corrected the prior native review,
recorded fixed `rya1281-library-count-carryforward` in bug ledger and explicitly
corrected it in this handoff. No spectrum/population/rate/abundance/source changed.
No bibliography provenance invented to silence the audit.

Statev222/SEQUENCE record the scoped inventory, not molecular-strength or
completeness admission. Self-review checked canonical ownership, source/vendor
equality, row counts, stable IDs/raw labels, AIR consumer versus unknown primary
frame, both RV geometries, isotope/component-count interpretation, atomic source
screen lineage and current-main identities. Scientific HOLD; no independent
physical molecular blend identification or source/chemistry closure claimed.

The separate held CO.dat source is audited in `co_source_addendum_review.md`:
117783declared rows parse,191enter the AIR interval after the manifest-documented
vacuum conversion. Raw numeric interval192rows shares only171IDs;41membership
differences. It is outside the inspected automatic *.bsyn glob, so it was not
part of the predeclared five-family651row census. This is not a verdict on
dedicated production CO routes or molecular chemistry. Both sources remain
separate in the evidence; do not conflate source counts or frames.

The CO supplement adds2byte-exact outputs,6total new inventory results; original
COv1preserved and extendedv2rows.csvbyte-identical. No additional primary-release
or physical-opacity claim. No new bug found in source parsing/routing by this
limited audit; the library reporting error is the fixed ledger item.

Final artifact checks: exactly2436existing generator diagnostic strings, zero
added/removed versus preceding phase, exit1retained; both actual writer CLI
invocations accepted. MIRROR/version/SEQUENCE/whitespace pass. Global ledger
tests1failed/1passed (historical RYA587science_impact blank); new documentation
row has required evidence, no unrelated repair. Postcommit freshness and normal
pre-push guards run before handoff, no bypass.

## Exact next bounded task

Prepare a matched native LTE molecules-on/off control at the already-held native
reference atmosphere, fixed composition/O8.70/xi and identical atomic inputs,
using the full nine source files selected by native routing (651centres in the
20Aoutput interval), preserving native treatment of out-of-window wings. Such an explicitly
labelled legacy-source sensitivity test can progress while primary qualification
remains open; it does not authorize installing molecular opacity in the Bride,
multiplying independent normalized spectra, or adopting CNO abundances. Retain
native5750/4.5versus canonical5772/4.438mismatch if using the native reference.

Read-only preflight finds installed TS `bin/bsyn_lu` and `bin/babsma_lu` are
MacARM64executables; no synthesis executed yet. Source/build lineage and runtime
dependencies still need verification; binary presence is not qualification.
An isolated source-receipted build is an option if existing lineage is insufficient,
without replacing the installed engine. Existing bound-opacity point producer
is not automatically a full molecular-spectrum producer. This separates native
LTE sensitivity from the later same-transfer/matched EOS/partition/isotope/
normalization bridge work. Additional O diagnostics, C/N then Fe/reference star,
wider interactions and actual adaptive3D/all-band gates remain OPEN.

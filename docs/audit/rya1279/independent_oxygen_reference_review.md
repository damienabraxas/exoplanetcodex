# Independent Oxygen 1D reference checks — RYA-1281 / RYA-1371

2026-11-09. Scientific **HOLD**.

After PR619/main d1728c92 integration, the original continuation register versions
v244/v245/v246 map to v253/v254/v255; v256 records the clean main-based branch.
Original journal and review subsection IDs remain historical branch identifiers.
RYA-1278 and RYA-1281 are now Done for their reviewed checkpoint scopes.
Further validation continues under RYA-1371, with RYA-1279/1280 still In Progress;
these issue statuses do not certify a physical Oxygen model.
 This is independent numerical corroboration,
not a completed solar abundance model or 3D validation. RYA-1374 owns the later
band and stellar-class expansion; its dependencies are recorded in Linear.

## Reference engine and input identity

OpenMULTI3D source is pinned at `32085d26c46e9919a986397b5514a9728abe2074`.
The experimental registry is `config/reference_engines.yaml`; it does not add a
production abundance token. Private Mac GFortran build: 3,614,040 bytes,
SHA256 `7cefed6ecadb563c1e519fd4ace116fa680d8aa64ed42128f818466950494bf4`.
Framework BSD3 and application GPL3 licensing/attribution remain distinct.
The supplied atom inventory also contains Fe607, Mg86 and Ba06, so the adapter
accepts explicit element/atom identities rather than hard-coding Oxygen.
Those atoms and 3D modes have not been scientifically qualified here.

The held O41 atom and MARCS 5750 K / logg4.5 / A(O)8.70 source atmosphere are
byte-identical to the held MULTI reference. The supplied OpenMULTI3D O41 file is
not substituted: it differs by UV973A transition sampling NQ9 versus101.
Native OpenMULTI3D consumes a pressure-graded, interpolated mesh. Foreign
populations are compared by physical height, never by row index. The endpoint
roundoff clamp is0.0162m, inside the independently calculated float32 accumulation
bound4.761m; this is not atmosphere extrapolation.

## Numerical and spectral results

The dense native solve took5.6794s, frequency refinement10.7498s and six-inclination
control6.6364s. All use one Mac thread and the original600s ceiling. All depths
reach the same final iteration13. Independently reconstructed **stored** physical
rate residuals are at most3.96e-8 and number conservation2.23e-16. This is not a
fresh fixed-population radiation replay of the independent native engine.

On the source logtau500 interval[-3,0], the current Bride versus held MULTI
triplet-level maximum absolute departure difference is0.006935; OpenMULTI3D
versus held MULTI is0.009253 after explicit height interpolation. This interval
is a diagnostic selection, not a contribution function. Large differences in
other/outer populations remain (Bride versus held MULTI all-level maximum8.155).

The fixed Bride local/held-MULTI-shape/LTE transfer controls completed88.663s on
the dense grid, without updating SE populations. The original seven-ray flux
endpoint replay is exactly unchanged in normalized flux and continuum.
The native and Bride spectra share exact vacuum sample values, generated from
explicit AIR windows with the canonical converter once, at spacing<=0.005A.
The canonical inverse is a one-step approximation; round-trip differences are
not corrected by new shifts. All measurements below use the actual returned
canonical AIR coordinates and reject gaps>0.0051A.

Finite +/-0.5A residual-intensity areas, mA (not full EWs or abundance corrections):

| AIR centre A | Bride local | held MULTI b on Bride LTE | Bride LTE | OpenMULTI3D NLTE | OpenMULTI3D LTE |
|---|---:|---:|---:|---:|---:|
|6158.170|5.4554|5.4572|5.4376|5.2420|5.2235|
|6300.304|3.4322|3.4322|3.4322|3.1900|3.1900|
|7771.944|84.8787|85.1145|78.4150|84.3318|77.7050|
|7774.166|73.8893|74.0730|68.5960|73.3677|67.9452|
|7775.388|58.7403|58.8620|54.9671|58.2440|54.3801|
|8446.400|105.5778|105.7858|98.7287|101.2925|94.2396|

Triplet inter-engine area differences are0.4963–0.5469mA. Native frequency
sampling doubled2->4 changes the selected areas by at most0.00330mA; Lobatto
inclinations5->6 change them by at most0.02742mA. Both pass the predeclared
diagnostic0.1mA screen; this is not a scientific error budget. No default is
selected from these controls. The initial seven-inclination request failed
because the pinned native Lobatto implementation supports at most six;
the failed run and original plan are preserved and the adapter now refuses it.

After the native per-A to per-Hz Jacobian, the absolute continua differ by about
0.5percent near the triplet,0.96percent near8446 and1.3–1.5percent near6158–6363.
Thus agreement in normalized profiles alone does not establish background parity.
EOS, collider populations/interpolation, depth representation, finite top boundary
and bound-free output handling remain distinct. Both native controls exclude
background metal line blanketing. No observational blends or instrument response
are included in this theoretical intensity comparison.

A further predeclared output-only control explicitly enables the native atom's
bound-free opacity in the output pass. It completed6.1671s; saved populations,
LTE populations and height grid are bit-identical, and all six selected NLTE
and LTE finite-window areas are exactly unchanged. This measured switch effect
does not explain these optical discrepancies; no general opacity parity follows.
The supplement `oxygen_output_opacity_v1` receipts188 further raw files and12
exact compact copies.

## Failures and corrections retained

1. Single-ray Lightweaver output squeezes the ray dimension. The old cap helper
   broadcast this into an N-by-N spectrum. Original transfer-v1 arrays are invalid
   and remain external with receipts. Actual-backend regression fails for the old
   helper and passes for the corrected helper (five cap tests). Multi-angle flux
   calculations are unaffected by this shape defect.
2. The initial new reference reader labelled native saved per-A intensity as
   per-Hz. Direct `rt_patch_mod.f90::save_intensity` inspection establishes
   erg/s/cm2/A/sr. The corrected reader preserves native values and converts via
   lambda_vac_A^2/c_A_s; an independent Planck-law regression checks the units.
   Earlier absolute-intensity interpretations are invalid; no production outputs
   were changed. Dense-analysis-v3 and comparison-v3 are the current evidence.
3. The native atom-grid output contains no6363 samples. A provisional integration
   interpolated across that absence; it is inadmissible. Dense sampling confirms
   the held atom lacks the6363 transition (6300 is present). A nearly zero signal
   in this window is a missing-transition control, not validation of[OI]6363.
4. Earlier exploratory reader failures (inactive continuum sentinel, Python import,
   numpy JSON serialization) remain preserved outside the compact archive. The
   strict binary reader requires exact byte counts/axes and accepts upper-level0
   only for inactive continuum slots with exactly zero rates.

The two software defects are in `data/audit/run_bug_ledger.csv` with tests and
invalid artifact identities. No failed/invalid arrays are promoted to science.

## Evidence and next work

### Boundary and molecular-reservoir follow-up (v245)

Pinned native `rt_patch_mod.f90::set_upstream_boundary` sets top incoming
radiation to zero. `HOP/interpol_mod.f90::optical_depth` instead defines the
reference optical-depth coordinate; its extrapolated top optical depth is not
the RT incoming boundary. A separate zero-cap Bride fixed-population transfer
control changes local finite-window areas by less than0.0003mA at6158/777/8446,
and-0.012825mA at6300. This does not measure the SE response to a boundary change.

At the actual native consumed T/ne/rho/height, Bride LTE triplet level fractions
agree with native values to0.0125–0.0379percent, while raw total atomic O differs
by up to21.43percent in cool upper layers. Native species exports show that
selected CO/OH/H2O contain up to17.624percent of the Oxygen accounted for by
OI+OII+those molecules. This is not a complete molecular/ion budget.

The native baseline species-export replay is bit-identical in its solved O
populations to the previous dense run. A separate native molecular-equilibrium
suppression control retains metadata and sets `tmolim=1K`, below all atmospheric
temperatures. Consumed T/ne/rho remain bit-identical; selected molecule densities
become zero; the matched raw atomic-total discrepancy drops below0.075percent.
The fractional LTE-level difference is unchanged. This establishes a bulk
chemistry distinction, rather than a triplet-level energy/degeneracy discrepancy.

| AIR centre A | native normal chemistry mA | native suppressed chemistry mA | Bride empty-molecule context mA |
|---|---:|---:|---:|
|6158.170|5.2420|5.2505|5.4554|
|6300.304|3.1900|3.3645|3.4322|
|7771.944|84.3318|84.9547|84.8787|
|7774.166|73.3677|73.8917|73.8893|
|7775.388|58.2440|58.6359|58.7403|
|8446.400|101.2925|101.9450|105.5778|

Native chemistry suppression changes the triplet areas by0.3919–0.6229mA and
leaves its triplet difference from Bride below0.105mA. It does not resolve8446
or establish which combined physics is accurate. Chemistry changes the native
atomic reservoir, colliders and background together; this is not a unique
single-rate attribution, a fitted parameter choice or a recommended default.
The control takes6.0881s, reaches iteration14 and stored-rate residual4.06e-8.

Two further exploratory failures are retained. `tmolim=0`, documented upstream
as disabling molecule loading, segfaults during pinned native initialization;
the adapter now refuses that configuration. No source/binary patch or relaxed
budget is used. Initial consumed-state continuum diagnostics added scattering
twice; pinned LW `Source/LwMiddleLayer.pyx` already adds it to background chi.
Only consumed-state-v2 reports are used for extinction. Corrected chi500 ratios
span0.938–1.009 with normal native chemistry and0.938–1.010 under suppression;
these are depth-dependent source-table comparisons, with the native internal
versus canonical AIR-converter distinction retained. Earlier v1 continuum
reports are invalid, preserved and excluded.

`oxygen_boundary_eos_v1` receipts538 further raw files and44 compact copies.
The next experimental step is a fully explicit molecular-LTE equilibrium
preflight on the Bride side, with inspected source fits, range validity and
element conservation, before coupling it into Oxygen SE. This is distinct from
molecular NLTE or production adoption.

That preflight has now begun (`oxygen_molecular_preflight_v1`). Eighteen installed
RH definitions were inspected and receipted. Fifteen have unsupported depths in
the full column: most end9000K (five deepest cells), H2O ends6000K (24 cells).
Full-network range admission remains HOLD. A standalone CO-only attempt is
refused by LW because its chemistry solver requires H in the molecular nuclei.
The separately declared CO+H2 LTE control runs on51 temperature-supported
independent EOS states, explicitly excluding original top-to-bottom indices51–55
only from this diagnostic. It is not a truncated RT atmosphere or full EOS.
It completes in three chemical iterations, conserves O/C to8.89e-16/5.56e-16,
and puts up to16.998percent O in CO. The native chemistry and this incomplete
network are not interchangeable. The LW warning in Newton relative-change
division is preserved in the log; the zero-valued component has not been
attributed, and saved quantities are finite.
No coupled SE, molecular line opacity, molecular NLTE or model admission follows.

A full-column HOP standalone EOS preflight now succeeds (`hop_provider_preflight_v1`).
It binds the existing pinned shared library through signatures parsed from its
committed C header, using standard Python ctypes because CFFI is absent from both
existing environments. No package/environment install or native binary change.
Exact source abundance/data identities and all native reference files are checked.
All56 input T/ne/rho values are bit-preserved. Saved OI/OII/HI/HII/CO/OH/H2O
densities agree with native outputs within5.67e-8 relative, consistent with native
float32 storage. Reference chi500 differs9.05e-5 relative; trace-continuum ownership
and exact opacity-source configuration still need measurement. This is an EOS
interface preflight, not a coupled Bride background provider or scientific model.
It provides a concrete route to complete chemical-state packets without extending
the held RH fits beyond their declared ranges.

### Trace ownership interface (continuation v246 before main reconciliation)

The9.05e-5 reference-opacity gap is resolved by exposing the native non-hydrogen
atomic `exclude_double_continua` selection in an isolated GPL derivative. The
4874-byte patch, COPYING and reproducible fresh-worktree builder are committed;
original reference source/build remain unchanged. A fresh44.760s single-thread
build produces1,164,328-byte libhop.so, SHA256
`306e3cf806ab335e29a1b685a8de05e708942e9ac588ec385f44ad90d41df2d0`.
The loader requires its exact upstream pin, patch, tracked diff, header and
binary receipt. This does not replace the independent reference engine.

Selection removes one composite background record (`O_I OP/DK65LEV`) already
represented by the model's continua. Stored native chi500 then matches exactly,
while T/ne/rho and exported EOS species remain unchanged. Before-init,
hydrogen, molecular-model, negative-excitation, invalid-stage and after-load
operations all refuse. Hydrogen/molecular special-case semantics are not
generalized from atomic Oxygen.

The prospective seven-wavelength by56-depth opacity/scattering check uses
2500/5000/6158.17/6300.304/7774.166/8446.4/16000A in AIR, converted canonically
once. All extinction/scattering values and residual absorption are finite and
nonnegative; the independent scattering call is bit-identical and opacity calls
leave species unchanged. Canonical5000A opacity agrees with stored native chi500
within9.74e-8. This verifies a non-square sampled packet, not UV/IR spectra or
all-band physical applicability. `hop_trace_provider_v1` receipts13 raw files
and7 exact compact copies, including the final fresh build; prototype build
failures remain external with receipts. Six redirected diagnostic logs are
durably copied under external `rya1371-hop-trace-diagnostic-logs-v1`.

The first binding generation refused a noncanonical character declaration; the
next compile found a non-TARGET pointer alias assumption. The final patch uses
the correct interoperable declaration and an ASSOCIATE alias. The initial
opacity probe failed on a relative plan path after cwd changed; resolving the
plan before initialization fixes it. No failed output is promoted. Next is an
explicit complete-state fixed-population transfer check; native and LW historical
physical constants/source functions must be distinguished before SE coupling.

### Clean main-based continuation

The newer reference/provider changes are applied to a fresh branch/worktree from
main `d1728c92` (PR619), preserving the original branch and main's reviewed
source-assets/observation separation and optional-backend gates. RYA-1278/1281
are Done for that checkpoint; RYA-1371 owns the active continuation. Version
mapping is explicit above and in register v256; no prior journal is rewritten.

All64 focused compiled/reference/register tests pass on this base. Every array
of a fresh seven-by56 HOP packet replay is bit-identical to the prior packet.
All131 compact copied artifacts and raw receipts verify. Standard artifact-guard
diagnostics match the clean main-tree baseline exactly:5120, zero additions or
removals. Current path check146 present/18 broken, with exactly the same18 broken
keys as main; reference library212 documents/129 indexed/83 unindexed. These are
known global debts, not a green full-suite claim. Target mirror/version checks
pass. Products, gold/reference, stellar parameters/constants, model registry and
observational holdings remain byte-identical to main.

Sources-page regeneration exposed two pre-existing bibliography schema blockers:
Civis2013's unsupported `vizier` verification tag and Mallinson2024's prose in the
verification-state field. Metadata/notes are preserved and normalized to the
supported publisher state after primary A&A/arXiv checks; the current direct
Mallinson A&A403 is recorded. No line/grid/gf data or roles change. Sources-page
generation then passes (219 rows). This reproducibility defect is ledgered.
Native logs, GPL text and the exact unified patch retain their source whitespace;
the remaining continuation whitespace check passes.

`data/results/rya1371/oxygen_reference_v1/archive_receipt.json` receipts948 raw
files and64 exact compact copies. Raw native binaries and invalid large arrays
remain external; no5MiB guard is bypassed. Generator ownership is registered.
Earlier dense-analysis-v1 intensity labels are retained solely as wavelength
input provenance for the Bride transfer; only their wavelength axes were consumed.

Next: measure native top-boundary effects separately;
compare consumed EOS/continuum/collider packets on a matched physical mesh;
verify reference native fresh-field balance; then return to the solar intensity
and flux benchmarks with continuum, blends, instrument response and uncertainty
handling qualified. Extend the atom for6363 only with inspected atomic data.
Keep molecular OH/CO/H2O and UV/IR/cool-star applicability as explicit separate
gates under RYA-1374. Oxygen 1D scientific completion remains open.

### Complete native-state fixed transfer (v257)

The isolated HOP provider now binds the complete native56-node T/ne/rho geometry,
41-level Oxygen n/nStar/totals, native H/He perturber densities, trace-continuum
ownership and native historical-constant background source. LW atomic/broadening
formulas remain explicit; no SE/population update occurs. Five native Lobatto
angles retain the source default-REAL rounding. Native and LW lower-boundary,
formal-solver, profile and atomic-constant conventions remain distinct.

Both fixed cases complete in7.1776s, with n/nStar/totals bit-identical before/after.
The common dense vacuum wavelengths are identical; AIR conversion occurs once.
Finite +/-0.5A intensity-area differences (LW/HOP minus native, mA):

| AIR centre | native NLTE populations | native LTE populations |
|---|---:|---:|
|6158.170|-.161453|-.161441|
|6300.304|+.002317|+.002319|
|7771.944|-.889573|-.872779|
|7774.166|-.778204|-.764277|
|7775.388|-.606745|-.596575|
|8446.400|-.710438|-.691925|

6363 remains an absent-transition control, not a measured Oxygen line.
The8446 difference decreases from4.285mA to0.710mA,6300 from0.242mA to0.0023mA.
Triplet differences are similar for LTE/NLTE populations, suggesting remaining
transfer/profile conventions rather than a population-specific discrepancy.
Raw continuum differences remain0.109-0.142percent. Visual inspection confirms
small profile differences, especially6158 wings; no continuum or line adjustment
is fitted. These are finite-window vertical intensities, not full EW/flux/abundance.

Native iteration uses linear formal transfer, native emergent output quadratic;
LW uses BESSER. Next compare predeclared linear output on both engines, preserving
all population/background inputs, then inspect mesh/lower ghost/profile effects.
Do not choose a production default by best agreement. Coupled SE, fresh native
field/rate balance, observed solar qualification and RYA-1374 expansion stay open.

Evidence: data/results/rya1371/hop_fixed_transfer_v1, exact compact copies plus raw
receipts. Three constructor preflights refused invalid option types/empty filename/
perturber lookup; corrected before the completed run, no failed output promoted.
Five focused packet/state/units/purpose tests pass in1.28s. Scientific HOLD.

### Matched formal solver and Stark semantics (v258)

Native linear output completes2.433s. Every native analysis array is bit-identical
to quadratic output, including full n/nStar/T/ne/rho/height. LW linear paired
transfer completes5.137s with frozen populations. Continuum differences fall to
1.5e-6--1.95e-5relative (below0.002percent); triplet area differences remain
-.825/-.722/-.564mA. Formal order explains most continuum offset, not line wings.

Source inspection identifies a real recipe difference: pinned HOP
bbtrans_mod::prepare_stark treats gstk_in=0 as a Cowley estimate, whereas the
Bride converter currently omits Stark broadening for held gStark=0. No claim is
made yet that either interpretation is the original MULTI convention. A native
linear Stark-off run (disable_stk=T) completes2.283s; thermodynamics/mesh/nStar
stay bit-identical, but maximum n changes0.21023percent. Accordingly the next LW
paired transfer consumes that run's complete populations, rather than attributing
all spectral change to transfer at the old state. It completes4.911s.

| AIR centre | matched linear NLTE difference mA | Stark-off matched-state NLTE difference mA |
|---|---:|---:|
|6158.170|-.155750|+.004719|
|6300.304|+.000112|+.000090|
|7771.944|-.824604|+.345538|
|7774.166|-.722228|+.274136|
|7775.388|-.563708|+.185684|
|8446.400|-.624791|+.078609|

Differences are LW/HOP minus native, finite +/-0.5A disk-centre intensity areas.
The LTE residual pattern is similar. Visual inspection confirms improved6158
wings and retained small triplet/8446 differences. No physically preferred
Stark/default/abundance choice follows from reducing these differences.

Further source difference: LW Unsold/Barklem broadening includes a He term;
native active-line damping uses total HI only. Next measure that term with an
explicit broadening-only He-off control, keeping EOS/continuum and native state
unchanged, then inspect remaining profile/constant effects and original MULTI
recipe semantics. Source files/fingerprints and all comparisons are in the
formal/Stark audit. 44 focused tests pass1.73s;389 raw files receipted and37
compact copies verified under hop_formal_stark_v1. Coupled SE/observed solar/3D
and RYA-1374 gates remain HOLD.

### He broadening-only control (v259)

At the Stark-off native state, disable the He contribution solely in LW elastic
broadening. All native populations/nStar, wavelength axes and continuum intensities
remain bit-identical to He-on, both LTE/NLTE cases. Actual perturber reads verify
the He-off route; EOS/He abundance/background opacity are untouched. Remaining
NLTE intensity-area differences (LW/HOP minus native, mA):6158 -.000320,
6300 +.000123,7771 +.031445,7774 +.021511,7775 +.011687,8446 -.084768.
All six finite windows agree within0.085mA under these diagnostic conventions.
This is a transfer comparison milestone, not an Oxygen physical-model admission.
Pointwise normalized differences remain up to0.00130 in triplet profiles; integrated
area agreement cannot certify exact profile/centre/Voigt-kernel parity. No empirical
shift or scale was fitted. Next audit original MULTI Stark semantics and native/LW
profile/Unsold conventions before implementing qualified coupled-state physics.

Initial He-off v1 was invalid: the prospective flag was recorded but not forwarded
to Perturbers(seed), so actual reads still used He. Caught in source/read inspection,
not used in comparisons; preserved manifest/report and external array receipts.
Corrected v2 is separately fresh and independently audited. Audit asserts actual
He-off reads and exact unchanged continuum/population/axes; no invalid promotion.
Evidence: hop_helium_v1,12 raw receipts/10 exact compact copies. Scientific HOLD;
source/native model/production defaults unchanged. RYA-1374 remains downstream.

### Original semantics and consumed table audit (v260)

Pinned pymulti4d6aa3c MULTI2.3 mul23_sub.f:1439-1446 sets GS=GQ*ne for
nonnegative GQ, with a special zero-input automatic branch only for Hydrogen.
Thus Oxygen GQ0 means no Stark in that source, agreeing with Bride's interpretation;
native HOP Cowley is a recipe difference, not evidence that Bride dropped the
original Oxygen Stark contribution. Producer identity/parity of the held MULTI
record is still a separate open gate. Original source fingerprints are retained.

Native atom_convert_mod writes and reads back its own table. The actual consumed
O.txt retains256 line identities but rounds AIR line centres to0.001A and loggf to
0.001dex. Triplet table-minus-source AIR positions are+.000131/+.000455/+.000476A;
8446 components-.000287/-.000322/-.000249A. These are source-determined differences,
not fitted wavelength corrections. The6158 native f differs+.0563percent.
They offer a concrete next test for pointwise residuals; no exact causal amount
is assigned before a controlled replay. Source/atom inputs remain unchanged.

The project canonical coordinate remains vacuum below2000A. External table
conversion uses its own refractive index and emits different numerical table
coordinates there; the largest table-coordinate difference is1.0448A near972A.
This is NOT a demonstrated RT wavelength error: the engine converts that coordinate
back internally. Roundtrip and FUV physical conventions need separate qualification
before RYA-1374/all-band claims. Initial table-audit-v1 labels called the canonical
UV coordinate AIR; retained as superseded/invalid labels, v2 explicitly distinguishes
canonical vacuum/AIR from native table coordinates. No UV result admitted.

He-switch routing now has a runtime guard checking actual compiled-context
perturber reads. A negative regression refuses the initial miswire; the full fresh
v3 pair is bit-identical to v2 for every output array. The committed replay harness
also reproduces its report byte-exact. Defect added to run_bug_ledger; failed v1
not promoted.45 focused tests pass1.70s. Source/profile/guard evidence has9 raw
receipts/8 exact compact copies in hop_profile_source_v1. Path register147 present,
18 pre-existing broken. Next: source-determined table precision/profile controls,
then qualified collider/rate/EOS integration and fresh native balance, before the
observed solar intensity/flux and broader parameter/element validation. All
physical-model, abundance, production and3D gates remain HOLD.

### Source-determined optical table controls and kernel probe (v261)

Three fresh paired fixed-population transfers separately use native consumed table
centres, strengths, and both. Native table receipt is checked before/after; canonical
AIR-to-vacuum conversion once. Support is ten transitions in the predeclared windows,
including two neighbouring6156 lines. Original levels, atom file, continua/collisions,
all other transitions and fixed n/nStar remain unchanged. Transfer-only wavelength
properties decouple photon and level energies exactly as this diagnostic intends;
this object is explicitly NOT qualified for coupled SE. No fitted shift/strength.

The three pairs finish5.650/7.794/5.977s. An independent audit proves populations,
nStar, emergent continuum and output wavelength arrays bit-identical to the original
He-off baseline in every case. Centre-only versus prior maximum NLTE normalized
profile differences:7771 .000273 vs.000528;7774 .000250 vs.001299;
7775 .000266 vs.001256. Combined centre+strength differences are.000269/.000253/
.000260. Thus known table rounding explains most triplet pointwise mismatch.
Strength-only effects are small and do not uniformly improve agreement.

| AIR centre | centre-only area difference mA | centres+strengths area difference mA |
|---|---:|---:|
|6158.170|-.000319|+.002155|
|6300.304|+.000124|+.000060|
|7771.944|+.031403|+.029012|
|7774.166|+.021298|+.023027|
|7775.388|+.012191|+.015370|
|8446.400|-.079535|-.085083|

These are finite +/-0.5A vertical intensity areas (LW/HOP minus native), not flux or
an abundance. Combined plot inspected; remaining small shape differences visible
only on a residual scale.6363 still absent. Constructor preflights expecting eight
immediate contributors refused the ten-window support; source table inspection
corrected that explicit identity before successful runs. Failed manifests preserved.

Pinned native SRC/rt_patch_mod::bound_bound path directly calls qlg_bot_profile,
which calls QLG-bot fp32; this is independent of standalone line-list Voigt options.
A fresh isolated GPL driver compiles the exact pinned qlg_mod source and compares
its QLG-bot and fp64 reference with independent SciPy wofz.61 damping values1e-5..10
and801 offsets-20..20 yield48,861 samples. QLG-bot max absolute H error.00335024;
max relative error.00641260 where H>1e-8. Native fp64 reference itself has max absolute
error1.82053e-6. These are kernel-grid errors, NOT an Oxygen emergent-intensity bound
or proof of the remaining spectral discrepancy. Probe v2 completes.966s; original
checkout/binaries untouched. v1's ambiguous threshold key is retained external with
receipts; v2 explicitly records the1e-8 floor. No profile default selected or changed.

47 focused compiled/reference/state tests pass2.45s.41 raw receipts/30 exact compact
copies verified in hop_table_profile_v1. Reference library remains212documents,
129indexed/83unindexed,17non-documents; known global debt not claimed green. Next:
measure the kernel on actual matched atmospheric damping support, then an isolated
native accurate-profile transfer control with explicit population response, before
qualified collider/EOS/SE and fresh-field validation. Physical1D/3D/solar/all-band
admission remains HOLD; no production or source atom changes.

### Measured damping-support probe and main reconciliation (v262)

Sample the actual LW damping formulas used by the matched table/He-off transfer:
ten optical transitions,56 consumed native T/ne/rho states. The isolated HOP packet
supplies the same native HI perturbers; original broadening formulas and source
constants remain declared. Damping spans3.50742e-5..0.131501. External LW atomic,
broadening, table, constants and compiled module fingerprints are bound; no RT/SE
or population update. These are measured LW damping inputs at native atmospheric
states, NOT a claim to reproduce native stored damping arrays.

The actual pinned QLG module is then probed at all560 distinct damping values and
401 offsets over[-20,20] (224,560 samples). QLG-bot max absolute H difference from
SciPy is.00326579; max relative.00651143 where H>1e-8. Finite-grid profile integrals
differ at most.000241852relative (0.024185percent). Native fp64 reference max absolute
error on this support is1.73230e-8. Native probe v3 takes1.581s. This still does not
bound an emergent Oxygen profile/area/abundance or demonstrate the full residual's
cause. Accurate-kernel RT control must measure its own SE/population response and
boundary/broadening differences, not silently replace the original reference.

Earlier sample/probe reports remain receipted external; final sample v2 binds the
external modules and final probe v3 includes the damping input in its initial
manifest and correctly distinguishes measured support from a generic grid.28raw
files/6 exact compact copies verified in hop_atmospheric_profile_v1.47 focused tests
pass (see session journal for elapsed time). Kernel/transfer/state/observation
qualification remain separate; physical1D/3D/all-band gates stay HOLD.

While pushing v261, the ordinary pre-push gate correctly refused a newly unrecorded
main landing: RYA-1373 PR618 at10f0dc95. Inspected merged diff, atlas audit/methodology
and both chronological ticket comments (8990d0ae, a1c68001). It adds the30-entry
Spectroscopic Frontier inventory/schema/source-hashed generation and coverage/gap
matrix;12 focused tests pass per merge evidence. It confers no scientific model
admission or public HTML publication; RYA-1377 owns the latter after1370. The FeII
historical LTE/NLTE source conflict is recorded, no abundance affected asserted.

Current worktree remains based on d1728c92; this records newer main reality without
silently importing its code or rewriting our experiments. No original/other worker
source modified. SEQUENCE stamp advances to10f0dc95 and records the actual landing;
register includes its scoped Done state. No hook bypass. Future integration must
retain both independent bug-ledger additions. Artifact-guard comparison continues
to name the clean d1728c92 tree baseline rather than claiming a full test of newer
main code. Next: isolated accurate-kernel native transfer, fresh-field balance and
qualified collider/EOS/SE, then observed solar validation.

### Isolated accurate-kernel native control (v263)

A fresh private pinned worktree/build changes only three active non-H profile call
sites and adds a helper invoking the existing fp64 reference. Input/offset arithmetic,
final fp32 conversion, normalization and native frequency Jacobian are retained.
Original source, executable, HOP provider and default reference registry unchanged.
Opt-in variant uses a separate experiment spec, exact patch/HEAD diff, executable/
library receipts and original-reference identity. Default dispatch is unchanged;
unknown variants and unexpected including staged source changes refuse.

Build v3 passes98.882s within one-thread600s budget. Two earlier compiled prototypes
remain unqualified: v1 canonical Git header/receipt mismatch; v2 started during the
packaging correction and correctly refused changed patch identity. Both logs/receipts
remain external and receipted, not numerical evidence. Canonical patch provenance
proved by exact reverse-application/two-file scope and final full Git diff equality.
No original mutation, cached-failure adoption or hook bypass.

Native accurate-profile solve completes6.275s, all eight slabs iteration13, reconstructed
stored-rate residual3.53223e-8, conservation2.22e-16. Thermodynamics, mesh, LTE populations,
consumed atom table and all collision matrices are bit-identical to original-profile
native control; n changes at most0.0103702percent. Paired transfer consumes its own
complete new n/nStar, not the old state. No fresh-native-field certification.

| AIR centre | own-state area difference mA | max NLTE normalized-profile difference |
|---|---:|---:|
|6158.170|+.002496|.0001644|
|6300.304|+.000024|.00001774|
|7771.944|+.003511|.00008342|
|7774.166|+.003542|.00007834|
|7775.388|+.003396|.00006991|
|8446.400|-.092746|.00038599|

Differences remain finite +/-0.5A vertical intensity windows, LW/HOP minus native.
Triplet pointwise residuals decrease substantially from table-matched QLG-bot;
8446 persists, so the kernel is not its full explanation. Source broadening conventions
and constants remain distinct. Next isolate van der Waals contribution, with its
own native population response and matched-state pair; never select no-broadening
as a physical default. Then qualify original-source widths, collider/EOS/SE and
fresh-field balance before observed solar fitting.

Actual default native replay finishes2.488s and every exported population, atmospheric
and spectrum array is bit-identical to the original pre-variant result.49 focused tests
pass; original/default preservation is an actual execution check as well as a unit
regression.392 raw receipts/34 compact copies verified in accurate_profile_v1.
Original nine-model registry and production/source atom unchanged; physical1D/3D,
observational and all-band gates HOLD. Explicit profile variant is diagnostic only.

### VDW mechanism and typed native wavelength roundtrip (v264)

An explicit native VDW-off control (accurate-kernel, zero-source-Stark, linear RT)
completes5.911s. Saved native params confirm disable_vdw=T and disable_stk=T.
Mesh/thermodynamics/nStar and all collision matrices remain bit-identical; n changes
at most0.852342percent, stored-rate residual3.95880e-8. The paired LW transfer uses
that complete new state; all elastic VDW terms and perturber reads are demonstrably
absent, natural widths retained. No omission of physical VDW is adopted.

| AIR centre | own-state VDW-off area difference mA | max NLTE normalized difference |
|---|---:|---:|
|6158.170|+.000976|.00002386|
|6300.304|-.00000105|.00001540|
|7771.944|+.004886|.00009297|
|7774.166|+.004620|.00008655|
|7775.388|+.004136|.00007696|
|8446.400|+.007140|.00005469|

Finite +/-0.5A vertical intensities, not flux/abundance. The8446 discrepancy falls
from-.092746 to+.007140mA; source VDW conventions explain most of its remaining
area difference. All six windows now within.0072mA under declared diagnostic
conventions. A residual-only six-panel plot is inspected and archived: smooth
small antisymmetric triplet/6300 residuals persist after the QLG structure is removed.
Each curve is its own native population comparison, not an attribution at fixed old n.

A typed Fortran probe links the verified original HOP library/module and directly
calls native air2vac/vac2air; no converter formula is copied.59 FUV line centres:
consumed rounded table maps back to source vacuum within.000492563A, unrounded
native roundtrip within6.97e-7A.197 other lines: consumed-table difference within
.000497446A, unrounded roundtrip within9.1e-13A. The prior~1.04A FUV table-coordinate
difference is NOT a photon-frequency error; internal inverse mapping recovers the
source vacuum frequency to table precision. This qualifies numerical BB line-centre
encoding only; canonical observation coordinates remain vacuum below2000A. UV
opacity, bound-free/rate accuracy and all-band applicability remain separate gates.

Initial probe SIGABRT from unresolved native library install-name build/libhop.so
is preserved with an authored failure note. Fresh v2 invokes pure helpers from their
required read-only library cwd, no original/binary modification or dependency install.
Runtime stderr now saved. Native saved nml has embedded NUL bytes; read without
altering its source bytes. No invalid output promoted.

A further measured precision difference: canonical project's documented one-step
AIR inverse versus native iterative inverse differs by2.23e-5A at the triplet,
2.81e-5A6300 and2.04e-5A8446. This is a small documented approximation difference,
not a wrong-medium finding or a reason to change global production conversion.
Next use an explicit source-bound native vacuum photon-centre packet in this
fixed-population diagnostic only (observation axes retain canonical conversion
once), then qualify source-aware full widths and collider/EOS/SE. No empirical shift.

50 focused tests pass2.43s.210 raw receipts/25 exact compact copies verified in
vdw_wavelength_v1, including residual plot. Generated archive was expanded before
commit to include that plot; initial compact-copy directory retained external under
/private/tmp/rya1371-vdw-wavelength-initial-archive. Physical1D/observed/3D/all-band
admission remains HOLD.

### Explicit vacuum photon-centre packet (v265)

Final prospective fixed-transfer control uses actual original-native vacuum photon
centres, qualified by the typed roundtrip packet and its source/library/compiler/
atom/table receipts. No new converter formula, fitted shift or global SSOT change.
The ten optical transfer properties consume already-vacuum centres directly; output/
observation axes retain their canonical conversion once. All other transitions and
original levels remain unchanged; this transfer object remains unqualified for SE.
Loader refuses changed producer/report/inputs/table/schema/precision; runtime audit
checks every actual optical photon centre against its source-bound packet.

Both native NLTE/LTE population cases finish5.857s, preserving n/nStar/continuum/
output axes bit-identically to the preceding VDW-off control. No native SE rerun or
population fitting occurs. One-axis control changes only the small inverse precision.

| AIR centre | fixed-state area difference mA | max NLTE normalized difference |
|---|---:|---:|
|6158.170|+.000976|.00000835|
|6300.304|-.00000112|.000000821|
|7771.944|+.004885|.00002326|
|7774.166|+.004619|.00002230|
|7775.388|+.004134|.00002173|
|8446.400|+.007663|.00001741|

All six pointwise differences below2.4e-5; all finite-window areas within.0077mA.
The8446 finite-area difference slightly increases despite better pointwise agreement:
blended truncated windows are not full equivalent widths or an agreement optimizer.
6363 remains a missing-transition control. Four-stage residual plot inspected: the
antisymmetric remainder falls with the verified vacuum packet; no visual continuum
adjustment or result-selected wavelength/abundance shift was introduced.

51 focused tests pass2.75s, including no-double-conversion and strict purpose routing.
10 raw receipts/10 compact copies verify in native_vacuum_centres_v1. Path register
153 present/18 pre-existing broken. Target mirror and artifact guards checked against
explicit clean d1728c92 baseline; source constants, stars, model registry, observations,
production products and gold/reference remain unchanged.

This completes a controlled transfer-comparison milestone, not a completed physical
Oxygen1D model. Next stop stripping physics: qualify/restore full original-source
width prescriptions, collider and molecular/atomic EOS population interfaces; then
coupled SE with fresh radiative/collisional balance and solar intensity/flux benchmarks,
continuum/bracketing, blends, instrumental response and uncertainties. Parameter/xi
sweeps, CNO then Fe/reference-star and RYA-1374 expansion remain downstream. Original
reference and default dispatch remain intact; no production or scientific admission.

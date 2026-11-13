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

### Current-main integration for review (v266)

PR620 was reported mergeable=false/dirty after RYA-1373 landed. Merge origin/main
10f0dc95 into this continuation, without merging/publishing the PR. Sole content
conflict is the append-only bug ledger. Resolve using full-row multiset union and
original record bytes:96 own rows,91 main rows,97 combined rows; the added main row
is RYA-1373. No row dropped/reworded, no historical verdict rewritten.

Incoming Frontier package generation initially reports a stale coverage snapshot,
as this branch already corrected two bibliography verification-state fields.
Freeze its source-evidence pin at immutable pre-integration ebbaf109 and regenerate
only provenance-bearing draft metadata. All30 entries, gaps, citations, release/
publication states equal main exactly; no scientific inventory fact is changed.
The pin now points to actual retained bibliography bytes rather than older metadata.

Use the existing Python3.9.6/jsonschema environment /private/tmp/rya1370-qa read-only
(PYTHONDONTWRITEBYTECODE=1), without installs or environment changes. Atlas generation,
--check and all12 unittest mutations pass2.739s. Repeat the51 compiled RT/reference/
register tests after integration; no scientific gate promoted. Source constants,
stars, catalogs, graded products, observations and gold/reference remain equal to
current main; only the draft Frontier provenance metadata differs outside experimental
Oxygen changes. No original/other worker worktree modified.

All15 Oxygen archive groups verify:2794 raw-file receipts and303 exact compact copies.
Current main integration keeps the global artifact diagnostic debt explicit and
separate from focused passes. Next resume at full source-aware widths, collider/
atomic-molecular EOS interface qualification, coupled SE/fresh balance, and observed
solar benchmarks. The controlled numerical transfer milestone is complete; physical
Oxygen1D,3D and all-band model completion stays open. PR620 remains draft.

### Original full-width oracle and prescription meaning (v267)

Compile the actual pinned pymulti4d6aa3c MULTI2.3 DAMP routine and BCONST, following
its native DP recipe (PREC.dp, REAL-to-REAL(8) include conversion). Only workspace
capacities grow MK30->64/MLINE150->512. Recording WDAMP shim observes natural/VDW/
Stark/total widths and ADAMP without replacing formulas. ATOMID hard O; NH aggregate
is used solely by SUM NH1:5, not represented as H level populations. No RT/EOS/SE.
Fresh native build.945s; all41levels/256lines/56 matched native states evaluated.

Interpretation clarification: the original DAMP header explicitly describes a fixed
He/H0.1 assumption in its historical Unsold recipe; the current executable uses
8.411 and1.283984e-12. Earlier HI-only wording describes runtime density input,
NOT proof that the physical coefficient excludes helium. Comments and executable
prescription are kept distinct, not used to infer an extra He term automatically.

Measured actual LW/original VDW ratios across all matched states:

| LW prescription |253 Unsold lines|3 ABO lines|
|---|---:|---:|
|without added LW He term|.96528903-.96528908|.9999844340|
|with actual native HeI term|1.0223318-1.0245011|1.0326208-1.0343547|

No empirical coefficient adjustment or preferred physical recipe follows. This
quantifies why dropping the added LW He term alone was not exact source-width
reproduction. Original natural widths match source GA exactly, all source O GQ0
Stark widths zero, total equals natural+VDW. Source mass15.999amu versus LW15.9994
is recorded, not silently changed. Native default-REAL literal promotion is visible
in exported constants (e.g. BK1.38066198349e-16, PI3.14159274101); this oracle does
not establish the unavailable held MULTI producer's binary/build identity.

Five fresh source-law controls retain EOS densities explicitly: halve HI, +10K,
xi.25/2/4km/s. Half-H VDW ratio is exactly.5. +10K follows T^.3 Unsold/T^((1-alpha)/2)
ABO with max relative error2.94e-11. All three xi changes leave every natural/VDW/
Stark width bit-identical; native ADAMP responds through Doppler width. Independent
SI normalization differs at most4.24e-8, below an8-default-REAL-epsilon bound9.54e-7.
These are width-component controls, not new physical stellar atmospheres, HSE/EOS
or full SE parameter sweeps. v1 report serialization refused a NumPy scalar;
corrected fresh v2 preserves all budgets and source physics; failed attempt retained.

59 raw receipts/12 compact copies verified in multi_width_v1. Next export the exact
native source-law coefficients from one source evaluation (no fitted regression),
verify them on all remaining states and controls, and restore full positive widths
through an explicitly named experimental component before coupled SE qualification.
Source-law temperature/density response must be owned and tested; no extra He factor
may be silently stacked on a source prescription. Physical1D/solar/3D/all-band HOLD.

Standing bibliography refresh found three newly held CrII PDFs. First two pages
inspected and metadata indexed only (Gurell2010 preprint, Ward2023, Engstrom2014).
Gurell's2018 template-header/date versus2010 arXiv stamp caveat is explicit; no final
journal/DOI guessed, no numerical table/gf/grade adoption. Library now215 documents,
132 indexed/83 old unindexed,17 non-documents. Generated sources page222 rows passes.
Atlas source evidence pin refreshed to metadata commit a0bbd0ec; all30 entries/gaps/
citations/states remain unchanged. No dependency install or website publication.


### Restored positive source widths and radiative-input audit (v268)

Exact source-law coefficients exported algebraically from native sample0, never
regressed or fitted to spectra. All remaining states and five independent native
controls reconstruct within6.66e-16 relative. Explicit numeric T support
4045.89892578125--9779.8701171875K is a component validation domain, not a physical
stellar-grid hull. Source owns its historical He convention; extra He is refused.
Immutable neutral-H total density comes from HOP, never artificial H level storage.
The component is opt-in and fixed-population only; coupled SE remains refused.

Actual full-positive-VDW transfer completes5.646s, preserving native full-physics
accurate-kernel n/nStar/totals, native vacuum photon centres and explicit table
rounding controls. Six finite +/-0.5A NLTE area differences (Bride-native,mA):
6158 +.00557010,6300 +.00000391,7771 +.00366057,7774 +.00366909,
7775 +.00345414,8446 +.00748132. Max normalized profile error2.40849e-5.
No VDW omission, empirical width multiplier, wavelength shift or continuum fit.
These are controlled theoretical vertical intensities, not solar validation or
complete equivalent widths. Missing6363 remains an absent-transition control.

Compiled-backend audit inspects all256x56 actual transition Qelast/aDamp arrays:
all positive VDW, original DAMP gamma agreement5.55e-16, independent SI aDamp
agreement7.77e-16. Source15.999 versus LW15.9994amu remains explicit. Width law,
constructor, transfer, comparison and backend proof are preserved byte-for-byte.

Radiative-input qualification remains open: source GA is below the same model's
own Einstein A for source transitions1->6 (1355.597669vacA:2951 versus4196.02943/s)
and1->4 (6302.046237vacA:2.11e-5 versus.006400562/s). NIST's lifetime definition
is inverse sum of outgoing spontaneous rates:
https://physics.nist.gov/Pubs/AtSpec/node18.html . The incomplete model graph is
not a complete independent lifetime and no source correction is adopted.
The actual pressure widths dominate: replacing GA diagnostically with the greater
of GA and this partial upper+lower graph would alter total damping by at most
.000452355 fractionally for1355 and3.85932e-9 for6300 over the56 states.
This is an algebraic width sensitivity, NOT an emergent-spectrum/SE bound.
Thus the6300 input inconsistency is unlikely to explain its visible discrepancy;
independent lifetime provenance and direct response remain to test. No source atom
or radiative A modified. This caveat is data qualification, not a new code defect.

Focused current tests41pass/11skip: all11 skips are historical Fe fixture availability,
not newly bypassed Oxygen tests. Three new source-law tests exercise immutable H,
density/T/micro roles, invalid domain/species and double-He/VDW-omission refusal.
Scientific1D/observed-solar/3D/all-band HOLD; next collider/EOS and fresh coupled
balance, with radiative lifetime source qualification tracked explicitly.

All17 archive groups verify:2867 raw receipts/329 exact compact copies. Staged
plain artifact diagnostics exactly match clean baseline5120; no new diagnostic.
Path check156present/18pre-existing broken. Register pointer v268 passes.
Overlay inspected; residual magnitudes come from arrays, not visual overlap.


### Collider continuation and current ORCH integration (v269)

Fresh origin/main45901e76 (PR621, RYA1378) integrated without conflict in88d9dc64;
merge retains both sequence histories. This adds the current atomic gap register
and ORCH step without changing any experimental input. Eight incoming gap-register
tests and41 focused Oxygen tests pass locally after merge. No CI/Sirius run.

Source inspection establishes the next collision contract, not a new discovery of
an untested converter fix. Existing bride_gencol_controls/native_spline already
record source direction differences. Pinned mul23_opacu.f GENCL4660--4740:
CH uses ground-H NH(1) with downward coefficient, CH0 uses raw directed endpoints;
CP uses proton NH(6). Source DAMP instead reduces SUM NH(1:5). Thus the neutral-H
width packet must NOT be substituted for the full collision H-level packet.
Installed LW collisional_rates.py150--180 uses H[0] and CH upward coefficients;
its generic CH0 is downward only. Pinned OpenMULTI3D atom_mod.f90 CH is downward,
CH0 ion-pair upward plus LTE inverse (optional H-minus energy control), and its
nh1 comes from native composition. These are distinct hypotheses, not a universal
class-name mapping. Reuse the existing named direction/interpolation controls;
explicit ground-H/proton/H-minus roles and native EOS reservoir consistency are
required before the HOP provider's active-SE refusal can be removed. No new
rate direction, density substitution, reverse channel or producer identity adopted.
Raw inspected source hashes: /private/tmp/rya1371-collider-next-audit.txt (local
resume aid); authoritative source pins/build receipts remain registered. Physical
1D/solar and1374 gates stay open.


### Occupation-aware native Hydrogen collider packet (v270)

Separate GPL-derived HOP interface exports exact native collision nh1=
comp%elements(1)%ions(1)%levels(1)%n and its occupation probability separately.
Existing species exporter supplies total HI, HII and H-minus. No synthetic H level
array, rate-direction change, coefficient adjustment or active-SE admission.
Original provider/source/executable unchanged. New build is opt-in by explicit
provider_variant='collider', default remains trace; unknown variants refuse.

Important source semantics: composition_mod::init_H_levels2694--2759 explicitly
leaves occupation probabilities OUTSIDE stored level%n; its partition includes
weights and opacity later multiplies them. Thus raw collision ground/HItotal can
exceed1 without violating the source's chosen storage convention. Current56-state
range0.999943958573--1.000003871441; occupation-weighted ground/HItotal lies
0.999939903668--0.999999999999179. Never clamp or renormalize the raw collision
field, nor substitute it for total neutral-H pressure broadening. This is source
convention qualification, not independent physical endorsement of collision use.

Fresh source/build v4 succeeds28.625s,1thread,600s unchanged cap; exact canonical
patch equals tracked diff; library12591629257b87eb5e7c637fd6cc4756f08d8f958f3d6ad312df89b27296e38a.
Raw packet audit v4 compares ALL original state/species/opacity/emissivity/scattering
arrays bit-for-bit with the separate exporter. All match exactly. Three unloaded/
invalid-size calls(-1,0,56) refuse status3 and leave both output sentinels intact.
A subsequent fresh packet reproduces every collider array bit-for-bit. This does
not test wrong-size calls while a valid state is live; no such claim is made.
Native/exported fields finite/positive; occupation weights within[0,1], weighted
ground bounded by total. Unweighted raw values preserved even above total.

Failures preserved: buildv1 canonical-diff gate refused generated-header changes;
v2 original raw-only interface built but remains superseded/unadmitted. v3 adds
explicit occupation export but canonical header whitespace gate refused. Freshv4
includes actual generated header, passes without guard bypass. Packetattemptv1
refused missing OMP_NUM_THREADS=1 before native load. v2 and diagnosticv3 refused
an incorrectly assumed raw-ground<=total invariant; inspected actual init_H_levels,
corrected the semantic guard to occupation*raw<=total, not an empirical tolerance.
No source physics altered. Diagnostic v3 log retained /private/tmp/rya1371-collider-v3.log;
v4 raw packet/report under bride_mac_checks/rya1371-hop-collider-packet-v4.
14 raw receipts/5 compact files in hop_collider_v1 preserve build/prototype evidence.

42 focused Oxygen tests pass, including unknown-variant refusal before native
access. Next compare native collision matrix using exact native nStar, ground-H,
electron fields and source interpolation/directions; then qualify active O HOP
coupling and fresh R+C balance, followed by solar observations. Native CH0 effective
LTE inverse and original directed-only/author MN hypotheses remain distinct.
Scientific Oxygen1D/3D/all-band/RYA1374 gates remain open.

Final checks:18 archive groups verify2881 raw receipts/334 compact copies. Artifact
diagnostic set remains exactly clean-baseline5120;158 paths present/18 existing
broken. Register pointer v270 passes. Exact canonical patch retains blank context
lines reported by whitespace check; no stripping of source-diff bytes. All other
staged whitespace passes. No new adopted pipeline defect; source storage convention
and failed experimental assumptions are explicit above, not invented bug-ledger rows.


### Native collision-matrix numerical milestone (v271)

Read actual native O.txt/params and all eight final depth-rate slabs. All1106
records supported explicitly:741CE,39CI,302CH,24CH0,14 common temperature nodes.
Native classical electron/H fill recipes blank, all five rate scale factors1,
CH0_HMIP false. These checked prerequisites forbid silently adding other recipes.
Native nStar all41x56 is identical to the pinned analysis; T/ne collider state
identical to the qualified occupation-aware HOP packet. H raw ground and totalHI
remain distinct. No generic H-population array or density inference.

New GPL read-only driver calls actual compiled HOP hermite_derivative and
hermite_interpolation (log coefficients) and exports its native e_k. Qualified
existing collider-provider objects/modules/sources are all receipted; native
source/object files unchanged. Fresh standalone executable builds.609s using
one-thread600s cap; no C ABI guessing or reimplementation of interpolation.
Numeric Maxwellian table warnings are source checks, not coefficient corrections.

Reusable explicit collision_matrix computes C[destination,source,depth] from
source-native interpolated q, T/ne, raw ground-H, O nStar/energies/g and e_k.
CE downward, CI upward, CH downward, CH0 upward with O LTE inverse follow the
specified OpenMULTI3D contract. This is a named experiment, not adoption over
original directed-only CH0 or independent author MN source hypotheses.

All45,920 stored upward cells (including zeros) become BIT-IDENTICAL after the
same float32 output quantization; no tolerance-selected adjustment. Full double
reconstruction versus stored upward/inferred downward matrix max relative
5.95564208e-8, zero extra nonzero cells. Native largest rate1.28354353e11/s;
absolute discrepancy4016/s is output quantization at large rates, not a physical
uncertainty. Downward reference rates are reconstructed from saved upward rates
and nStar, NOT independently exported reverse rates. No fresh radiation or
nonlinear solve, physical collision-data or stellar model validation follows.

Labelled controls use the SAME native q/state, altering one convention only:

| control | max relative difference vs stored native nonzero matrix |
|---|---:|
| totalHI replaces raw collision ground-H |5.60315483e-5|
| directed-only CH0, omit effective LTE inverse |.999995572|
| upward instead of downward CH convention |9150.59662|

These maxima describe individual rate entries, not spectra/abundance effects or
whole-engine disagreements; alternate controls do not recreate the other engine's
interpolation/EOS. They discriminate conventions without empirical scales.

v1 standalone audit, v2 shared calculation and finalv3 quantized-output proof
all retained; finalv3 is qualified evidence. Native interpolation/rate inputs,
matrices and receipts live in native_collisions_v1 (22raw/7compact). Three meaningful
synthetic tests discriminate CH/CH0 orientation and reverse channel, Htotal versus
raw ground, electron-only roles, CE/CI Boltzmann and LTE ratios, and invalid states/
unsupported controls.45 focused Oxygen tests pass. No new production bug or silent
fallback discovered; previously documented source direction hypotheses remain.

Next implement an explicitly opt-in fixed-background collision adapter from this
verified calculation, bind immutable T/ne/H/nStar and atom/table identities, prove
actual backend consumes it before permitting the HOP provider's active-O guard to
change. Then run a bounded coupled O solve and fresh physical R+C/conservation/
restart checks. Table-rounded transfer-only line controls are NOT a self-consistent
SE atom; reconcile full atom/UV/BF support before interpreting SE differences.
Independent lifetime/collision evidence and observed solar intensity/flux validation
remain required; physical Oxygen1D,3D and1374/all-band gates stay open.

All19 archive groups verify2903 raw receipts/341 exact compact copies; plain
artifact diagnostics exactly match existing baseline5120.159 paths present/18
existing broken; staged whitespace and v271 version pointer pass. No hosted CI,
Sirius work, source/default/abundance/product/gold change, or physical gate sign-off.


### Actual compiled collision binding milestone (v272)

Opt-in FixedNativeCollisionBinding loads pinned matrix/report/atom identities,
checks source provenance, recomputes rates from native interpolation coefficients
and explicit native H/T/ne/nStar fields, and requires original source matrix parity
within256float64eps for SI/CGS round-trip arithmetic. No stored-reference C array
is copied blindly into the calculation. Source recipe remains explicit, no
empirical adjustment or new rate direction adopted. Every cached input is copied
read-only. Runtime checks T/ne/nHTot/height/vturb, O model identity/nStar/total,
level E/g/stage identities and whole-matrix non-stacking. NLTE populations may
change; collision parameters remain fixed. No fake H-level population lookup.

make_context adds a default-None collision_binding keyword; original branch is
unchanged. HOP Context remains detailed-static, with ZERO active atoms; the active
provider guard is intact. Backend proof uses the actual compiled LwAtom constructor
with detailed=False solely to allocate a detached computational atom's C storage.
It is never registered in Context or run through RT/SE. This directly exercises
the same compiled compute_collisions path used for active atoms without pretending
a detailed-static object owns C or enabling a premature coupled solve.

Final auditv3 proves compiled C exactly equals the source-recomputed matrix.
Six actual compiled calls with changed T/ne/nHTot/vturb/height/O nStar all refuse;
restore/replay gives exactly the original C, native n/nStar preserved. Compiled
fillDiagonal conservation normalized column residual5.02142e-16. Source-contract
LTE flow residual6.67876e-16; CH0 inverse closure is by construction, not independent
physical reverse-rate validation. Three successful callback calls including the
diagonal check. All256x56 source widths remain positive, compiled gamma agrees
5.55e-16 and independent aDamp7.78e-16; no broadening regression.

Failures recorded: first actual-backend probe correctly encountered unallocated
C storage on a detailed-static atom (Memoryview is not initialized). Failed log
and INVALID metadata are retained; exploratory source snapshot unavailable, so
never claim this is qualified reproducible numeric evidence. Source inspection
identified active-only allocation. v2 detached-storage proof passes and is
superseded by finalv3 conservation/closure proof. Initial synthetic mutation test
added1 to1e16/1e20 densities, which rounds away; test now uses an actual1% mutation
so the negative canary discriminates. This was a test-fixture mistake, not a guard
relaxation or source physics change. No production bug/silent fallback discovered.

47 focused Oxygen tests pass; new tests cover input immutability, stale fields,
O nStar/total/level identity, rate stacking and allowed NLTE population variation.
collision_binding_v1 archives6raw receipts/2compact copies; finalv3 is evidence,
v1invalid/v2superseded remain receipt-only. Explicit plan
config/experiments/rya1371_hop_collision_binding.json is backend-audit-only.

Next reconcile a full self-consistent SE atom with full source line/BF/UV support:
transfer-only optical table overrides are not suitable as a silent SE atom. Then
permit explicitly guarded active-O HOP binding, run a bounded coupled solve, fresh
R+C/conservation/restart and intensity/flux checks. This checkpoint validates actual
collision plumbing, not fresh radiation, statistical equilibrium, observed solar,
physical1D or3D/all-band completion. RYA1374 and CNO/Fe expansion remain gated.

Review hardening: plan purpose is now fixed-native-oxygen-collision-backend-audit,
explicitly checked by the audit and refused by the generic fixed-transfer runner.
This prevents the generic runner from ignoring a collision-binding field. Freshv4
replays identical rate/conservation/refusal results under that dedicated purpose;
v3 remains a numerical precursor, finalv4 is the executable plan closure. Archive
collision_binding_v2 adds8raw/4compact, preserving both versions without raw overwrite.

Final verification:21 archive groups/2917 raw receipts/347 exact compact copies;
plain generator diagnostic set remains baseline5120 exactly, no added/removed.
160 registered paths present/18 existing broken. Staged whitespace/version272
checks pass. Existing N refinement debt697 and missing engine-B atom remain open;
O/C tracker debt cells empty mean no known path, not physical completion.


### Full-source coupled numerical and restart milestone (v273)

First genuinely evolving full-source O/HOP diagnostic at the held5750K/logg4.5,
A8.70 reference node (not canonical solar parameters).41levels/256BB/22BF,
6795SE wavelengths90.0064--6876.574nm. This wavelength extent is internal atom
support, NOT validated dense full-spectrum/all-band/observational coverage.
All source E/g/stages/f values preserved; source photon/level identity exact.
All22 BF grids retain sorted wavelength A->nm and alpha cm2->m2 units, at least
3/4/26 usable source nodes apiece, none redward of source edge. No optical table
or photon-centre overrides, fitted wavelength/continuum/abundance or omitted widths.
One active O atom only; fixed native T/ne/rho/O LTE/total and trace HOP molecular
reservoir/EOS. No charge/HSE/EOS/chemistry feedback update or molecular NLTE.

Default active guard remains; explicitly labelled coupled_diagnostic path requires
qualified immutable collision binding, full source widths and no transfer-only or
continuum-only atom. Actual collision C matches source calculation before solve.
Original5rays/qscale2/Ng settings, rtol1e-6,1500outer,600s caps inherited from the
reviewed full-cap baseline. Top incoming zero (cap length0), LW thermalized lower
boundary; native ghost/boundary and constants remain distinct.

Thermodynamic qualification: imported nStar vs source/LW Boltzmann BB ratios
max relative.000450847471 (~.045%). Native table energies and native e_k reduce
residual to1.44250e-6. This exposes historical-constant/energy conventions; the
source atom is internally photon/energy consistent, but imported LTE plus HOP
historical-constant background is NOT a complete self-consistent LTE claim.
No artificial energy/temperature shift or fitted scale is applied.

The initial new harness omitted the previously required iteration plugin, treating
it as optional while inheriting scalar controls. Default runv1 reaches tiny n/J
changes (~1e-13) yet physical residual.00393860444;1500iterations exhausted in
419.30s. No computational/scientific pass. Exact old driver source hash preserved
under raw driver_source.py and compact archive; history/final populations retained.
This is a harness routing omission, not a new upstream-method discovery. New
validate_operator_contract requires the exact inherited plugin receipt before any
solve; missing/wrong receipt or resource refuses. Canonical bug-ledger row added.
No threshold/budget relaxation, clipping or reuse of a failed seed.

Full-coupling-cap process-local plugin is the existing pinned reviewed build;
zero-cap derivative explicit. Operator action identity is independently checked
against physical R+C at EACH SE step (256float64eps bound). Complete source/helper/
Python-constant/backend/plugin fingerprints retained. Final native-seedv3 passes:
11outer iterations,7.133s, max fresh physical residual4.89214990e-10. Final operator
action discrepancy2.88928e-14, conservation6.66e-16. Subsequent fresh-process
restartv2 passes3steps/2.902s, fresh residual4.28118e-11, maximum population change
7.23447e-9. Previous successful precursor versions and exact driver snapshots are
retained, finalv3/restartv2 are qualified source-closure evidence.

Independent check rebuilds R from278 saved per-transition up/down arrays plus C,
using separate incoming/outgoing sums. Physical residuals4.89214879e-10 and
4.28119096e-11; saved-net disagreement normalized below5.63e-16. This does not use
MALI Gamma as physical rates, and proves numerical closure/restart, not source-data
accuracy. All layers remain in SE. Native population differences remain large in
outer layers: minimum own/native.00182799 at level7/logtau500=-4.9177, maximum
1.00405121. Triplet-level max relative difference in the PREDECLARED diagnostic
logtau[-3,0] interval.00117150; interval is not a contribution-function measurement.
UV/BF/background/radiative-rate and outer-layer qualification remains owed.

Forward own converged populations through the full SOURCE atom and exact native
dense optical grid, with separate intrinsic continuum and no SE updates in the
postprocess. Finite +/-0.5A vertical intensity areas (mA; own minus native):

| AIR A | native | independently coupled | difference |
|---|---:|---:|---:|
|6158.170|5.103614|5.106943|+.003328|
|6300.304|3.184512|3.184579|+.0000676|
|7771.944|83.441868|83.461379|+.019511|
|7774.166|72.604638|72.617630|+.012993|
|7775.388|57.665659|57.672769|+.007111|
|8446.400|100.876253|100.894958|+.018705|

Max pointwise normalized difference.00126295; residual plots inspected, including
antisymmetric6300/triplet structure. Earlier source-table precision controls make
rounding a plausible contributor, but this own-state result alone is not a unique
attribution.6158 has small sampling/support steps retained.6363 absent-transition
control remains unqualified (3.72e-6mA only tails). Finite areas are not complete EW,
flux or abundance corrections. No observed Kitt/IAG/HARPS/CRIRES fit performed here.

49focused tests pass; new guards distinguish operator receipt/routing and explicitly
reject unqualified active binding/optical overrides. hop_coupled_v1 captures51raw
receipts/22compact copies, including exact exploratory generator snapshots. Standing
library audit still has83existing unindexed documents/12absent local_file pointers;
no new bibliography adoption or dependency install. Physical Oxygen1D/3D/all-band/
1374 and element dispositions/gold/defaults unchanged. Numerical coupling is now
working under the declared fixed-background conventions.

Next: matched-population UV/BF radiative-rate replay against native exports to
isolate large outer-layer differences, reconcile LTE/constant conventions without
source tuning, quantify effect on optical intensity and flux; then canonical-solar
parameter/control and ORCH continuum/telluric observation benchmarks with independent
O literature/atomic lifetime evidence. CNO then Fe/reference-star remains the order;
N697/missing engine-B atom debt remains open, not cleared by this O milestone.

Final verification:22archive groups verify2968raw receipts/369exact compact copies.
Plain artifact diagnostics remain exactly baseline5120;161paths present/18existing
broken. v273pointer and staged whitespace pass. Historical default-operator plan
is retained as failure evidence; current driver refuses it without the required
pinned plugin. Only the explicitly pinned MALI plan/restart are executable current
coupled controls. No full-suite green claim or source/default/gold/product adoption.


### Native zero-ghost boundary and numerical operator milestone (v274)

Matched native n/nStar/T/ne/rho replay covers every278 radiative transition
(556positive directional pairs), no populations updated. Source topology agrees;
no missing BB/BF graph links. The direct LW physical imbalance at these native
populations is.331134968 for a zero-incoming boundary placed AT the sampled top
node. OI1302/1304/1306 resonance channels dominate level7/top-cell differences.

Important representation correction: native PRECONDITIONING=T. Pinned
atom_mod::line_rate_coefficients850--900 and rt_patch::bound_bound_rates2803--3058
show saved coefficients reconstructed from thermal/deviation integrals containing
preconditioning corrections. They are effective SE weights, not gross Einstein
rates. Their small downward values are not evidence of a bad A coefficient.
The initial generic coefficient ratio (~2019max) is inadmissible as a physical
gross-rate comparison. Prototypes retained, canonical bug-ledger row added;
rate_semantics guard now refuses that mixed comparison and final replay labels
native-effective versus replay-physical quantities separately. Earlier stored
native balance remains a stored-SE diagnostic only. Native fresh unpreconditioned
radiation/rate export and the linearization population are still owed; no claim
that its saved3.53e-8 residual is fresh physical certification.

Source-defined geometry distinction: rt_patch::set_upstream_boundary3188--3215
zeros chi,S,I in the upstream GHOST plane; long_char_linear184--275 integrates
to the first sampled node using dtau=0.5*ds*(chi_ghost+chi_node)/abs(mu).
With ghost chi,S,I=0:

    tau = chi_node*ds/(2*abs(mu))
    I_node = (eta_node/chi_node)*(1-(1-exp(-tau))/tau)

Actual consumed mesh coordinates are verified, ds_top20656.773022242953m. No
fitted height, opacity or source coefficient. Python fp64 implementation follows
this discrete native formula, including thin-limit series; independent formal
integral tests cover thin/thick limits. This is a numerical ghost-plane match,
not a newly adopted physical atmospheric extension. Lower boundary stays LW
thermalized versus native extrapolated ghost, and all other physics unchanged.

At identical native populations, that source-defined ghost interval reduces own
physical imbalance.331134968 -> .00317959534 (about104times), without changing
rates/data/thermodynamic fields or excluding layers. Focus resonance-channel flow
differences shrink by orders of magnitude; individual native/own coefficients
remain semantically noncomparable. This establishes a major boundary-location
contribution, not unique attribution of every remaining UV/BF difference.

A coupled ghost control with boundary derivative held zero preserves operator
action but oscillates and exhausts the original600s budget (600.278s including
final atomic call). Failure/history/populations preserved and never used as seeds.
A fresh private numerical plugin adds the analytic FIXED-EXTINCTION emissivity
response F(tau)/chi, limit ds/(4*abs(mu)); physical boundary/rates are unchanged.
Source-based compiled build2.464s at clean held LW d33058f, installed source/backend
untouched, one thread and5MiB size guard preserved. Actual33975-entry compiled
response matches independent Python source formula2.22e-16 relative; frozen-pop
Gamma action agrees with physical R+C2.97e-14. No empirical preconditioner scale.

With that operator and ORIGINAL native initial populations, full-source ghost
coupled run passes8iterations/6.652s. Fresh physical residual1.36442330e-9;
independent per-transition R+C/incoming-outgoing reconstruction1.36442321e-9,
conservation6.66e-16. Fresh-process restart passes (fresh residual8.36701e-10),
max population change7.07550e-8. Original rtol1e-6/rays/qscale/Ng/600s budgets
unchanged; physical action checked at every step, no clipping or tolerance tuning.

All-level/all-depth own/native population ratios now .99107826--1.01237371,
maximum relative difference1.23737percent, versus99.8172percent in the zero-at-node
case. Minimum ratio level11 at outer logtau=-4.9177. Triplet max relative in inherited
diagnostic logtau[-3,0] band.00117136 (not a measured contribution function); all
layers remain in SE. Imported source/native constant/LTE convention mismatch and
remaining UV/BF/native gross-rate/lower-boundary checks remain open.

Forward spectra use the ghost incoming boundary on the ACTUAL dense diagnostic
wavelength grid, holding parent full-angle J for scattering. Actual coefficients
recomputed before the final ray pass; no sharp-opacity interpolation, no one-ray
angular-equilibrium replacement, no population update. I is evaluated at the
sampled top node, matching the native output pivot, not propagated through an
invented external atmosphere. No optical table/centre/strength overrides or fit.
Finite±0.5A disk-centre area differences(mA), own-native:
6158+.00332851;6300+.000067586;7771+.01953462;7774+.01301063;
7775+.00712215;8446+.01871824. Pointwise precision/shape discrepancies retained;
6363 absent-transition control still not validated. No observed-solar/flux/full-EW/
abundance result. Physical1D/3D/all-band/1374 gate remains HOLD.

52focused tests pass, including ghost formal-integral/invalid-input and explicit
rate-semantics discrimination. Native_ghost_v1 receipts109raw/37compact files,
including failure/prototype source snapshots and isolated operator build. Current
replay plans deliberately separate zero-node, source-ghost fixed populations and
coupled operator controls. Default make_context boundary remains zero-at-node.

Next: obtain genuine fresh unpreconditioned native rates or matching stage-specific
net-flow exports; inspect remaining interface/path, UV/BF and LTE/constant effects;
then canonical-solar parameters and ORCH intensity/flux continuum/telluric benchmarks
with independent Oxygen evidence. Reference node5750/4.5/A8.70 remains distinct
from canonical Sun. CNO→Fe/reference star and broader parameter sweeps follow the
Oxygen physical qualification; no stellar parameters chosen to hide boundary errors.

Final v274 archive verification:23 groups/3077 raw receipts/406 compact copies;
all SHA256/byte counts and compact identities pass. Register pointer passes;
path register164 present/18 pre-existing broken. Result-generator diagnostic set
is identical to baseline5120 violations; no global green claim. Staged diff check passes.


## v275 — stage-matched physical native rate observer (2026-11-11)

Worktree codex/rya1371, branch ryandamienschmitt/rya-1371-oxygen-reference-provider,
base fb5a7f85; fetched origin/main remains45901e76. No other worker changes present.
Mac-only one-thread source-derived numerical observer, no source/default/gold/product
or scientific admission change. Scientific HOLD remains; N697/missing engine-B atom
debt unchanged. Empty C/O debt cells do not imply scientific completion.

New isolated source patch repeats ONLY BB/BF rate integration against the existing
formal radiation field, using an optional LOCAL no-preconditioning argument. Solver
rate accumulators are saved/restored exactly. Original radiation/populations/opacities,
global preconditioning, maxiter1500/rtol1e-6/600s controls and Ng are untouched.
Additional sidecars capture the exact radiation-population stage, source normalization,
thermal/deviation BB integrals and BF integrals before the SE population step.
Pinned source32085d26c46e9919a986397b5514a9728abe2074 plus existing accurate-profile
patch. GPL notices preserved, original installed/reference tree unchanged.

Private build93.761s; source/input/executable/HOP/diff and patch receipts preserved.
Native run9.276s, every one of145 ORIGINAL saved binary outputs byte-identical to
qualified accurate-profile baseline (including final n/nStar/rates/spectra).104 gross
and104 population-state sidecars retained raw. Final8 slabs cover all41levels/56depths;
192 rate-integration contributions reconstructed at the final radiation sweep.
Native RT sweep12 feeds SE population update13, not sweep13. Initial new analysis
requested the nonexistent sweep13 sidecar and failed FileNotFoundError BEFORE any
measurement/output; corrected to explicit source-stage mapping, additionally requiring
exact bit identity of n/nStar in every contribution. No nearest-stage substitution.
No failed-state/source snapshot is claimed for that initial exploratory reader.

Independent reader sums native FP32 contributions in FP64 (not an exact native
accumulator-rounding-order replay), applies actual profile normalization and BB
statistical-weight/nStar thermal/deviation conversion. BF downward factors use actual
nStar ratio. Inactive sentinels, negative/nonfinite physical rates, invalid axes/maps,
truncation, trailing state bytes and mismatched population stages refuse.

Native physical R+C residual at the EXACT PRE-update radiation population:
7.58259721e-8. Following native population update changes populations by at most
5.60163388e-7. This is now genuine unpreconditioned physical-rate evidence, distinct
from the historical effective SE coefficients. It is NOT a fresh additional native
formal solution at the POST-update population. Native physical atomic data, boundary
and LTE conventions are still not independently validated by numerical closure.

Exact-state Bride replay with source ghost geometry:2.232s, fixed native PRE-update
n/nStar, no SE step; own physical residual.003179595408. Native physical residual in
the same replay7.582597216e-8. Gross physical directional rate graph/support agrees.
Largest relative individual rate contrast.04242321193 (4.2423percent), upward level2
→19 at second outer depth: native1.36935577584e-13/s, replay1.31126330556e-13/s.
This is an extremely weak pair, not the largest population/abundance effect.
Largest normalized difference in net flow is at level13/outermost depth (~.003180).
At inherited level7/top focus, net-flow contrasts are driven by1302/1304/1306A links
and the8448.568A vacuum component. Those now have physically comparable coefficients;
we no longer attribute meaning to mixed effective/gross individual ratios.
No fitted wavelength, continuum, abundance, energy or atmosphere adjustment.

55focused tests pass (previous52 plus thermal/BF normalization, truncated export and
stage/population canaries). Archive native_gross_v1:397raw receipts/13compact copies.
Standing bibliography audit remains83 pre-existing unindexed documents; no new
literature adoption/download. Initial stage reader refusal disclosed above; no new
production bug or silent scientific fallback discovered, no new bug-ledger row.

Exact reproduction, after `source config/hosts/mac.env`, with OMP_NUM_THREADS=1 and
OPENBLAS_NUM_THREADS=1. Project Python=/Users/ryanschmitt/codex/venv/bin/python;
LW Python=/Users/ryanschmitt/bride_rya1118/.venv/bin/python. Fresh output paths required:

```text
PROJECT scripts/build_native_gross_rates.py --source /Users/ryanschmitt/codex/bride_inputs/rya1371/openmulti3d-gross-rates-v1 --output /Users/ryanschmitt/codex/bride_mac_checks/rya1371-native-gross-build-v1
PROJECT scripts/run_native_gross_rates.py --source /Users/ryanschmitt/codex/bride_inputs/rya1371/openmulti3d-gross-rates-v1 --build /Users/ryanschmitt/codex/bride_mac_checks/rya1371-native-gross-build-v1 --plan config/experiments/rya1371_openmulti3d_oxygen_accurate_profile.json --output /Users/ryanschmitt/codex/bride_mac_checks/rya1371-native-gross-observer-v1 --execute
LW scripts/analyze_native_gross_rates.py --source /Users/ryanschmitt/codex/bride_mac_checks/rya1371-native-gross-observer-v1 --baseline /Users/ryanschmitt/codex/bride_mac_checks/rya1371-oxygen-accurate-profile-v1 --output /Users/ryanschmitt/codex/bride_mac_checks/rya1371-native-gross-analysis-v1
LW scripts/replay_hop_oxygen_radiative_rates.py --plan config/experiments/rya1371_hop_native_physical_rates.json --output /Users/ryanschmitt/codex/bride_mac_checks/rya1371-hop-native-physical-rates-v1
PROJECT scripts/archive_oxygen_reference_checks.py --native-gross-supplement --raw-root /Users/ryanschmitt/codex/bride_mac_checks --output data/results/rya1371/native_gross_v1
```

Next: partition the actual physical-rate/net-flow discrepancies by transition/depth
and isolate profile/LTE constants, cell-interface transfer and UV/BF contributions.
Then canonical-Sun intensity/flux and independently sourced observed Oxygen benchmarks
with ORCH continuum/telluric handling. Prior coupled/restart milestone stands, no
new coupled solve here. CNO→Fe/reference star/3D/all-band and1374 remain gated.

Final v275 verification:24 archive groups/3474 raw receipts/419 compact byte-identical
copies pass hash/size verification. Register pointer passes;165 registered paths
present/18 pre-existing broken. Generator diagnostics exactly equal baseline5120.
Full staged whitespace check reports only verbatim native.log trailing spaces and
unified-patch blank context/EOF; source/log/patch bytes preserved, no guard changed.
Handwritten code/config/docs whitespace check passes. No global-green claim.

## v276 — physical-rate decomposition and frequency-grid controls (2026-11-11)

Base d0253a61; same codex/rya1371 worktree/branch, fetched main remains45901e76.
Scientific HOLD unchanged. Mac-only, serial one-thread fixed-population RT; original
600s/rtol/rays/Ng controls retained. No new SE/population solve, adopted parameters,
source/default/gold/product changes, or observed solar/3D/all-band qualification.

Full physical-flow audit covers41levels x56depths across SEVEN controls (16072
level-depth rows;1946 distinct pair summaries=278pairs/case). All n/nStar/native R/C
and replay C arrays remain BIT-IDENTICAL across cases. Separate incoming/outgoing
reconstruction verifies delta-R plus delta-C plus native net equals replay net within
256fp64 eps on the gross-flow scale. Collision-output quantization contributes at
most1.164e-9 on this normalized scale, far below the remaining discrepancy.

Worst default residual is level13/outermost depth. The level7->13 component at
8448.568284A vacuum (near8446A AIR) contributes +.00318551 to that cell's gross-normalized
net difference, partly cancelled by other partners. Its upward rate differs.820515%.
This is distinct from the previously largest4.24% coefficient ratio on an extremely
weak transition. BF contributions become the worst cells in some alternative grids;
no element or layer is omitted or reweighted to manufacture agreement.

| Numerical control | Points | RT seconds | Max physical residual at exact native state | Top7->13 upward rate contrast |
|---|---:|---:|---:|---:|
| Default LW/source grid |6795|2.232|.0031795954|+.820515%|
| Default plus all-grid midpoint insertion |13589|3.758|.0029892427|+.772389%|
| Then local8446-window midpoint insertion |13691|3.523|.0029390292|+.759632%|
| Native wing endpoints only, LW node shape |6795|2.155|.0032284679|+.833126%|
| Native fractional support AND spacing |6795|2.136|.0007287804|+.018720%|
| Native-style grid plus all-grid midpoints |13589|3.569|.0025799191|+.613471%|
| Then same local8446-window midpoints |13691|3.592|.0028102399|+.726628%|

The apparently closer COARSE native-style result is NOT numerical convergence or
physical validation. Refining it moves back toward the independently refined default
grid family. Do not adopt.000729 as a scientifically preferred result. These tests
measure Bride frequency-grid sensitivity with fixed native populations; they do not
establish that the native RT engine itself is under-resolved, nor uniquely attribute
remaining differences. A native-side refinement and BF-radiation audit are next.

The grid oracle links read-only qualified HOP objects, verifies linked trans_mod bytes
against pinned source32085d26, and calls ACTUAL get_lam for256source transitions
(6288line nodes). Build.428s. Actual native QNORM10km/s, NQFAC2, legacyF confirmed.
Fractional offsets are applied around ORIGINAL source photon centres; no native-table
centre/f-strength replacement, E/g/C/width change, or empirical grid fit. This is a
native-STYLE fractional grid, not a claim of bit-identical native global axes or BF
radiation handling. Optional explicit numerical control only; defaults remain unchanged.
Support-only stretches the original doubled LW node shape to exact source blue/red
endpoints; full mode replaces spacing too. Tests enforce centroid/f/endpoints/support.

Both local refinements use the SAME inherited8443--8449A AIR diagnostic window,
converted ONCE to vacuum nm [844.5320129375,845.1321751772]. Global wavelengths remain
90.0064--6876.574194nm; this is internal SE integration support, not dense observed
all-band validation. Adding nodes preserves original line endpoints and frozen state.

Failed control: prospective full two-pass midpoint grid27177 exceeded the EXISTING
HOP packet cap20000 and refused before formal RT/SE (.482s). Failure/manifest preserved.
New helper had initially allowed30000, an authoring bound inconsistency; fixed helper
now imports the original20000 bound from its canonical packet owner, with discriminating test. Local refinement remains below
cap; no size/tolerance/budget guard bypass. Canonical bug-ledger row appended for this
new harness defect. Exact first failed helper source was not captured; no complete
failed-driver reproducibility claim. Two other exploratory prep attempts failed before
outputs (namelist whitespace regex, an unused incorrect path key); both corrected,
without substituted data or scientific results.

The simple cell-interface hypothesis was checked against ACTUAL params: AMR='P', and
long_char_linear already uses upstream ds for the first downward cell when AMR is
nonblank. Thus the guessed missing-upstream-spacing mechanism does not apply. Broader
interface/precision/lower-boundary effects remain open, not uniquely ruled out.

66focused tests pass (three existing upstream LW configuration/deprecation
warnings). New tests cover midpoint support/node retention, invalid/bounded grids,
native fractional offsets, full model line identity, and photon-centre/f invariance.
Heatmap residual_maps.png inspected; CSV retains every level/depth and every pair.
Raw archive supplement35files/30compact copies; four over5MiB NPZ files and the compiled oracle executable remain raw-only
with receipts, preserving the existing compact-artifact guard. Old native-grid module
and replay driver snapshots captured before later extensions where available; no
invented historical snapshots. No literature/source acquisition or science admission.

Reproduction: source config/hosts/mac.env; OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1;
PROJECT=/Users/ryanschmitt/codex/venv/bin/python;
LW=/Users/ryanschmitt/bride_rya1118/.venv/bin/python. Fresh owned output paths required.

```text
PROJECT scripts/build_hop_line_grid_oracle.py --output /Users/ryanschmitt/codex/bride_mac_checks/rya1371-native-line-grid-build-v1
PROJECT scripts/audit_native_line_grid.py --build /Users/ryanschmitt/codex/bride_mac_checks/rya1371-native-line-grid-build-v1 --output /Users/ryanschmitt/codex/bride_mac_checks/rya1371-native-line-grid-v1
LW scripts/replay_hop_oxygen_radiative_rates.py --plan config/experiments/PLAN.json --output /Users/ryanschmitt/codex/bride_mac_checks/RAW
LW scripts/audit_oxygen_physical_rate_controls.py --raw-root /Users/ryanschmitt/codex/bride_mac_checks --output /Users/ryanschmitt/codex/bride_mac_checks/rya1371-physical-rate-control-audit-v1
PROJECT scripts/archive_oxygen_reference_checks.py --physical-rate-controls-supplement --raw-root /Users/ryanschmitt/codex/bride_mac_checks --output data/results/rya1371/physical_rate_controls_v1
```

PLAN/RAW in dependency order (prefix PLAN with rya1371_hop_physical_):
midpoint1 /rya1371-hop-physical-midpoint1-v1;
midpoint2 /rya1371-hop-physical-midpoint2-v1 (EXPECTED refusal, no rerun needed);
midpoint8446 /rya1371-hop-physical-midpoint8446-v1;
native_grid /rya1371-hop-physical-native-grid-v1;
native_support /rya1371-hop-physical-native-support-v1;
native_midpoint /rya1371-hop-physical-native-midpoint-v1;
native_midpoint8446 /rya1371-hop-physical-native-midpoint8446-v1.
Plans contain prerequisite packet hashes and registered paths. Independent audit uses
all successful cases plus original v275 baseline; no failed control is used as a seed.

Next: predeclare a native NQFAC4 accurate-profile reference solve with the same held
O41/MARCS node, unchanged physics/rtol/budgets and serial Mac execution. Measure full
population/optical sensitivity; obtain independently stage-qualified physical native
rate exports for that refinement. Trace native separate BB/BF radiation/continuum
frequency handling before attributing cross-engine mismatches to atom constants or
stellar parameters. Then revisit source constants/lifetimes and canonical-Sun observed
intensity/flux with ORCH continuum/tellurics. CNO->Fe/reference star/1374 expansion
remain gated; existing N697/missing engine-B atom debt stays open.

Final v276 verification:25archive groups/3509raw receipts/449compact copies verify
SHA256/byte counts and exact compact identity.170registered paths present/18existing
broken; register pointer passes. Generator diagnostics exactly match baseline5120,
no new violations. Handwritten code/config/docs staged whitespace check passes;
full check reports only standard CRLF endings in the two verbatim generated CSVs,
which are preserved byte-identically. No guard/whitespace setting changed and no
repository-wide green claim. Final66focused tests pass in.91s, three upstream warnings.

## v277 — native refinement and explicit all-instrument UV/VIS/IR target (2026-11-11)

Base9e4f2d7e; same worktree/branch. Ryan explicitly confirmed the target must cover
all possible Exoplanet Codex instruments, not just VIS. The current canonical catalog
has27instruments/reference atlases and11modes, nominal90--5500nm. Target every catalog
row and regenerate the coverage matrix as it grows; unknown/new physical applicability
must remain explicit. Catalog statuses/dates are inherited metadata, not newly verified
manufacturer capabilities. Rejected catalog entries remain rejected. No element or
production/instrument-disposition change. Scientific1D/3D/all-band gate remains HOLD.

### Coverage is now recorded at distinct evidence levels

Internal O41 SE support reaches90.0064--6876.574194nm,41levels/256BB/22BF. That is
sparse atomic integration support, not continuous complete spectra. Existing detailed
emergent comparisons cover five narrow AIR output windows (including absent6363
control), six actual theoretical O-only line areas around6158--8446A. Prior work is not
full UV/VIS/IR observed validation. Internal vacuum coordinates are correct for atomic
radiation; AIR-observed products such as Kitt Peak require canonical conversion ONCE
at the observation comparison boundary, not an artificial change of atomic radiation
physics to an Earth refractive medium.

New dynamic audit maps ALL27catalog entries and11modes, preserving actual catalog
frames (AIR, vacuum, raw-product vacuum-wavenumber convention). IAG catalog bounds are
already wavelength-nm metadata; they must NOT be inverted as if they were raw
wavenumbers. Mode/setting/order gaps remain product-specific: a nominal CRIRES+ envelope
is not a continuously observed spectrum. The matrix distinguishes source atomic line
centres, internal SE envelope, intersecting diagnostic windows, numerical smoke output,
and independently validated full-band spectra. ALL full-band validation cells remain
NOT_ESTABLISHED. Current O41 has0atomic O centres in APOGEE1510--1700nm: catalog OH/CO/CN
features require molecular line formation; EOS molecule reservoirs do not supply those
lines. Background line completeness, chemistry feedback, stellar-parameter coverage,
LSF/R/rotation, flux/intensity, continuum/tellurics and UV/chromospheric applicability
must qualify separately. Existing admitted Codex products are not re-dispositioned by
this experimental O-only matrix.

### Coordinate and real forward-output implementation

Reference adapter previously accepted only AIR-A windows above200nm. Added explicit
`frame: vacuum_nm`/`step_nm`, scaling10 to native vacuum-A with NO AIR conversion.
AIR behavior remains unchanged and conversion is canonical once; source default
input namelists remain unchanged when new controls are omitted. Unknown/ambiguous
frame/unit keys, nonfinite/overlapping/nonpositive grids and over-budget output refuse
before allocation. Original fine sampling limits and20000-output-point bound retained.
Spectral analysis keeps vacuum primary; AIR coordinates at/below200nm are explicitly
unsupported (NaN with air_frame_supported=false, not missing intensity). New coordinate
helper lives outside the frozen native collision-reader module.

Actual output-only runs, preserving the same NQFAC4 SE model:

- UV/NIR/IR samples:130.15--130.25,1128.8--1129.3,3098.2--3099.0nm;10.061s.
- Full-BF/BB catalog-edge samples:90--90.02,130.15--130.25,1128.8--1129.3,
  3098.2--3099.0,5499.8--5500nm;3246points,10.186s. Bound-free output explicitly ON.

Each probe preserves140 original internal native binary files EXACTLY, including n,
nStar, all saved SE rates and internal radiation. Only output-axis/intensity/moment
files differ. Actual source save_intensity writes flux-labelled output-window moments
as well as intensity; under disk-centre-only quadrature they are NOT a certified
angularly integrated stellar flux. The paired audit excludes precisely those output
labels (out_lam/ie/flux/Jnu_out) and verifies all remaining binary identities. The first
preparation of this proof incorrectly treated output flux files as internal; refusal
was corrected after source/meta inspection, no physical discrepancy hidden.

UV/NIR/IR plot inspected. O-only absorption appears near130.2nm,1.129um and3.099um;
90nm/5.5um endpoints are line-free coordinate/continuum smoke checks. These are NOT
observed benchmarks, dense full-band products, molecular spectra or solar/chromospheric
validation. The90nm endpoint is below the current internal SE grid minimum90.0064nm;
that lack of internal-grid exercise remains explicit, not silently called qualified.
Complete human table and figure: native_refinement_coverage_v1/
rya1371-multiband-coverage-view-v1/{coverage.md,uv_ir_numerical_probes.png}.

### Genuine native frequency and radiation-field controls

Changed only native NQFAC2->4 on the exact held O41/MARCS5750/logg4.5/A8.70 reference,
then separately DECOUPLE_CONTINUUM=T->F at NQFAC4. All source atomic data, EOS, molecular/
background/width/collision conventions,5rays/rtol1e-6/Ng/maxiter1500/600s/one-thread Mac
limits retained. These are explicitly labelled numerical/radiation controls, no fitted
abundance/T/energy/continuum or automatic preferred-physics adoption. Native source
allocates separate BB/BF atom-transition views when decoupling is enabled; the shared
control removes that approximation. Its physical impact is measured, not assumed.

| Native case | Solve seconds | Observer seconds | All original binary outputs preserved | Independent PRE-update physical residual |
|---|---:|---:|---:|---:|
| NQFAC2 separated (previous) |6.275|9.276|145|7.58260e-8|
| NQFAC4 separated |10.079|15.065|145|6.85759e-8|
| NQFAC4 shared BB/BF field |10.088|15.123|112|6.85618e-8|

Observers all export192integration contributions/final8slabs at radiation sweep12 /
population update13, with exact n/nStar identity required. Native stored default rates
remain PRECONDITIONED SE coefficients; physical observer is a separate certificate at
the actual PRE-update field population, not a new post-update formal solve. New replay
qualification derives the expected binary count from the exact selected native producer
instead of hardcoding145; it still requires identity of EVERY original binary file.

NQFAC4 versus2 changes some populations by1.76923% (level35/top), while six theoretical
finite optical areas change by at most.00287880mA. Visible appearance is therefore not
sufficient evidence of full-model convergence. Shared versus separated NQFAC4 changes
populations at most.186477% (level41/second outer cell); six optical area shifts at most
.0000218013mA. Heights/T/ne/rho/source-logtau/nStar/C remain bit-identical across cases.
Native refinement is not yet an asymptotic convergence certificate.

Bride source-grid fixed-population physical replays at each EXACT PRE-update state:
NQFAC4-separated residual.000755275863 /2.245s; shared4 residual.000582205815 /1.969s,
versus.003179595408 at old separated NQFAC2. No SE update in these replays. Native
physical residuals remain below6.86e-8, source/native radiative graph unchanged. The
previous independently closed own coupled population solution differs from these native
post-update populations by at most.9144%/.9163%, respectively (old native2 max1.2374%);
this is a reference diagnostic, not a new own coupled solve or physical abundance error.
Shared4 worst fixed-state imbalance returns to level13/top. Constants/LTE conventions,
source table rounding, numerical grids, atomic lifetimes and observation handling remain
open; closer agreement alone does not choose a physics/default change.

### Failures and provenance, retained explicitly

- Legacy analysis v1 inherited misleading "stored physical-rate balance" wording;
  numerical arrays remain useful but report semantics withdrawn. Final v3 reads actual
  PRECONDITIONING and explicitly labels effective coefficients. Canonical existing
  matched-rates semantics bug row amended to include this analysis exporter.
- New v2 analysis used wrong report variable (`report` instead of `doc`) at finalization:
  NameError after arrays/plots, no complete report. Exact driver/partial output and manual
  failure record preserved; fixed v3 and all real multiband analyses complete. New ledger row.
- Coordinate helper initially added to source-pinned openmulti3d_output.py caused existing
  collision-binding fingerprint guard to refuse before RT/output allocation. Qualified
  reader restored BYTE-IDENTICALLY, helper moved to spectral_frames.py; actual frequency4/
  shared4 replays then pass. Refused reader and manual exception record retained. Guard
  never suppressed/relaxed; new provenance ledger row. No failed seed used.
- Initial new coordinate test compared decimal53000.02 by exact float equality; corrected
  to FP64 precision, without changing conversion, sampling or scientific tolerances.

85focused tests pass in1.36s, three existing LW configuration/deprecation warnings. Real
native output probes validate final report emission and explicit frame masks. New archive
supplement1494raw receipts/78compact copies, including failures; all originals remain.
No new dependency/download, upstream-source/gold/observational overwrite, Sirius/hosted-CI
workload, merge/publication or element-disposition change.

### Reproduction and exact next step

Environment: source config/hosts/mac.env; OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1.
PROJECT=/Users/ryanschmitt/codex/venv/bin/python;
LW=/Users/ryanschmitt/bride_rya1118/.venv/bin/python. Fresh output paths required.
Run native solves and observers serially. Full commands/input/output identities are in
per-run manifests/receipts and the following executable map:

```text
PROJECT scripts/run_reference_engine.py --plan config/experiments/PLAN.json --output /Users/ryanschmitt/codex/bride_mac_checks/RUN --execute
PROJECT scripts/run_native_gross_rates.py --source /Users/ryanschmitt/codex/bride_inputs/rya1371/openmulti3d-gross-rates-v1 --build /Users/ryanschmitt/codex/bride_mac_checks/rya1371-native-gross-build-v1 --plan config/experiments/PLAN.json --output /Users/ryanschmitt/codex/bride_mac_checks/OBSERVER --execute
LW scripts/analyze_native_gross_rates.py --source /Users/ryanschmitt/codex/bride_mac_checks/OBSERVER --baseline /Users/ryanschmitt/codex/bride_mac_checks/RUN --output /Users/ryanschmitt/codex/bride_mac_checks/PHYSICAL_ANALYSIS
PROJECT scripts/analyze_openmulti3d_oxygen_reference.py --source /Users/ryanschmitt/codex/bride_mac_checks/RUN --output /Users/ryanschmitt/codex/bride_mac_checks/ANALYSIS
LW scripts/replay_hop_oxygen_radiative_rates.py --plan config/experiments/rya1371_hop_physical_CASE.json --output /Users/ryanschmitt/codex/bride_mac_checks/rya1371-hop-physical-CASE-v1
LW scripts/audit_native_oxygen_refinement.py --raw-root /Users/ryanschmitt/codex/bride_mac_checks --output /Users/ryanschmitt/codex/bride_mac_checks/rya1371-native-refinement-audit-v1
PROJECT scripts/audit_oxygen_instrument_coverage.py --probe-analysis /Users/ryanschmitt/codex/bride_mac_checks/rya1371-oxygen-catalog-edges-analysis-v1 --probe-run /Users/ryanschmitt/codex/bride_mac_checks/rya1371-oxygen-catalog-edges-v1 --output /Users/ryanschmitt/codex/bride_mac_checks/rya1371-oxygen-instrument-coverage-v1
PROJECT scripts/plot_oxygen_multiband_coverage.py --coverage /Users/ryanschmitt/codex/bride_mac_checks/rya1371-oxygen-instrument-coverage-v1 --analysis /Users/ryanschmitt/codex/bride_mac_checks/rya1371-oxygen-catalog-edges-analysis-v1 --output /Users/ryanschmitt/codex/bride_mac_checks/rya1371-multiband-coverage-view-v1
PROJECT scripts/archive_oxygen_reference_checks.py --native-refinement-coverage-supplement --raw-root /Users/ryanschmitt/codex/bride_mac_checks --output data/results/rya1371/native_refinement_coverage_v1
```

PLAN prefix rya1371_openmulti3d_oxygen_: frequency4_accurate ->RUN
rya1371-oxygen-frequency4-accurate-v1 /OBSERVER rya1371-native-gross-frequency4-v1 /
PHYSICAL_ANALYSIS rya1371-native-gross-frequency4-analysis-v1 /ANALYSIS
rya1371-oxygen-frequency4-accurate-analysis-v3. Shared4 uses PLAN shared4 /RUN
rya1371-oxygen-shared4-v1 /OBSERVER rya1371-native-gross-shared4-v1 /
PHYSICAL_ANALYSIS rya1371-native-gross-shared4-analysis-v1 /ANALYSIS
rya1371-oxygen-shared4-analysis-v1. CASE for own replay:frequency4 or shared4.
Output-only PLAN uv_ir_output and catalog_edges ->RUN rya1371-oxygen-uv-ir-output-v1
and rya1371-oxygen-catalog-edges-v1; ANALYSIS names append-analysis before-v1.

NEXT: predeclare a genuine native NQFAC8 convergence control (within original runtime/
thread/tolerance budgets) for shared-field case; compare populations, physical R+C and
all radiative groups before broad T/xi/g/abundance drift sweeps. Resolve remaining
thermodynamic/LTE/source-constant and lifetime evidence without tuning. Add true angularly
integrated flux and independent UV/IR diagnostics to the catalog evidence ladder, then
canonical-solar observed comparisons using ORCH continuum/tellurics and explicit AIR/
vacuum/velocity/LSF handling. O41 cannot stand in for OH/CO/CN or complete other-element
blanketing; CNO then Fe/reference star and adaptive3D expansion remain the order. RYA-1374
stays gated; N697/missing engine-B atom debt remains open.

Final v277 verification:26archive groups/5003raw receipts/527compact copies verify
SHA256/byte counts and exact compact identity.181registered paths present/18existing
broken. Register pointer passes. Generator diagnostic set EXACTLY matches baseline5120.
Handwritten code/config/docs/ledger whitespace passes; full check flags only verbatim
native-log trailing spaces and generated coverage.md final blank line, retained exactly.
No guard setting changed, no repository-wide green claim. Complete file inventory is
the v277 commit; experimental output samples are not production or science admission.


## 2026-11-12 — Shared-grid convergence and thermal convention canary

Native shared-radiation NQFAC4→8→16 controls completed on the Mac under unchanged rtol1e-6, 600s, one-thread, five-ray and source/physics guards. Worst population changes fall from 0.12580024% (8/4) to 0.00330666% (16/8); six optical finite-area changes fall from 0.00481671 to 0.000195171 mA. Independent PRE-update physical residuals are 2.93908e-7 and 7.95078e-7; both within the original tolerance but the increasing residual warrants a precision investigation before further doubling. Each observer preserves all112 original binary outputs. UV/IR output-only probes preserve107 internal binaries each. At130.15–130.25,1128.8–1129.3,3098.2–3099nm, 16/8 maximum normalized-intensity differences are3.40887e-6,4.03285e-7,3.28632e-7; edge probes90–90.02 and5499.8–5500nm are unchanged. Sparse vacuum probes are not full-band/observed/flux validation.

Own fixed-native-population replay residual improves from0.000582206 (4) to0.000224179 (8) and0.000221166 (16); no fresh own coupled solve is claimed. An analytic LTE convention canary reads the actual HOP and installed LW constants and verifies the held compiled LW header/build receipt. HOP h=6.626176e-34,k=1.380662e-23; LW h=6.6260755e-34,k=1.380658e-23 SI. Self-consistent BB/BF controls close below2.20e-16. Imported LTE populations with LW blackbody yield line-centre defect0.000225475; HOP blackbody with LW BF coefficients yields per-frequency defect0.000234811. These are analytic detailed-balance defects, not integrated RT or a unique causal attribution of the remaining NLTE residual. No constants, temperatures, rates or adopted physics changed.

All raw receipts and compact archived bytes verified. Two tooling failures retained: convergence plot v2 shadowed its fine array with a CSV handle; analytic preflight attempted abs(list). Fixed variable names/array conversion; final convergence v4 metrics exactly match v1 and thermal v2 adds provenance checks. 89 focused tests pass (three existing LW warnings). Scientific HOLD remains: held5750K/logg4.5/A8.70 is not canonical solar; no true flux, observed Kitt Peak/ORCH, molecular/full-instrument or3D admission. Next: qualify a coherent opt-in thermodynamic convention jointly across LTE, background, coefficients and collisions; investigate native residual floor; then flux/canonical solar and observed ORCH checks before parameter sweeps, CNO→Fe and reference-star validation. RYA-1374 remains gated.

Evidence: data/results/rya1371/shared_convergence_v1; final raw audits rya1371-shared-convergence-audit-v4 and rya1371-thermodynamic-canary-v2. Drivers: run_reference_engine.py, run_native_gross_rates.py, analyze_native_gross_rates.py, analyze_openmulti3d_oxygen_reference.py, replay_hop_oxygen_radiative_rates.py, audit_native_oxygen_refinement.py --shared-convergence, audit_oxygen_thermodynamic_conventions.py. All source identities and failure snapshots retained in receipts.


## 2026-11-12 — Exact same-stage LTE and collision reconstruction

RYA-1371 continuation: exact LTE/collision convention reconstruction (Mac, no new RT or SE)

Branch ryandamienschmitt/rya-1371-oxygen-reference-provider; worktree /Users/ryanschmitt/codex/rya1371; starting9d9cee99ecfdfaaee78376f9e9e61a89dd0c50be; freshly fetched origin/main45901e76 unchanged. Original source/build, physical constants, frozen gold, observed spectra and solver inputs unchanged. Existing N697/engine-B atom debt and RYA-1374 gates remain open.

New source-derived audit reconstructs all742 same-ion-stage level pairs at56 matched depths (41552 ratio checks). Native atom_mod computes e_k/T in REAL(4), then level exponents/ratios in REAL(8), from the consumed six-decimal eV table. Same-stage ratios cancel the stage normalization and ionization terms. Maximum native ratio reconstruction error9.3258734e-15. Omitting the source precision round gives1.443342e-6; using unrounded source energies with HOP h*c/k gives3.745093e-6; using original source cm-1 energies with installed LW h*c/k gives0.00045084747. This isolates precision/consumed-energy effects from the larger thermal convention difference. It does not reconstruct Saha balance, stage totals, ionization lowering or EOS.

All780 connected collision pairs /43680 pair-depths satisfy native-LTE reciprocity: bound FP64 recipe max signed pair defect2.43986e-16. Native file semantics are explicitly verified: upward triangular C serialized FP32, inverse rates reconstructed in FP64 from original cgs nStar. Reproducing this exact representation gives bit-identical C to the saved reader result; reciprocity defect2.20845e-16. The independent recipe coefficients are not globally bit-identical to stored C (maxrelative5.96e-8), because storage rounds upward values. No false exact-FP64 identity claim remains.

A diagnostic LW excitation redistribution preserving EACH ion-stage total, never supplied to a solver, breaks reciprocity with the unchanged collision matrix by0.00022547456. This confirms an LTE-only replacement is thermodynamically inconsistent; it does not prove that this explains the NLTE residual or validate microscopic collision data (inverse rates themselves use native LTE). No correction, fitted temperature, coefficient replacement, model adoption or observed/flux/full-band/3D admission.

Verification: OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1, held LW Python ran scripts/audit_oxygen_lte_reconstruction.py --native-analysis <raw>/rya1371-oxygen-shared16-analysis-v1 --native-run <raw>/rya1371-oxygen-shared16-v1 --replay <raw>/rya1371-hop-physical-shared16-v1 --output <raw>/rya1371-lte-reconstruction-v1. Source/library/state/replay hashes checked before/after, native producer linkage and exact state checked. All tests/test_*rya1371.py:74 passed in1.44s. New discriminating collision test checks reciprocal, disconnected, one-way, changed-population and invalid-diagonal cases. Disconnected entries are explicitly NaN with a presence mask, not successful zero measurements. Existing library audit still83 unindexed documents; no new documents acquired.

Failures preserved: three diagnostic preflights refused before output allocation. Initial full-matrix exact-C assumption failed at5.96e-8; whole-matrix FP32 interpretation also wrong because inverse entries are reconstructed; final inverse-bit check exposed4.42e-16 differences caused by taking ratios after SI conversion. Fixed by reproducing actual upward-only storage and cgs inverse arithmetic, without relaxing a threshold. Exact failed driver/failure metadata retained; final exact representation assertion passes. These are audit implementation errors, not upstream physical defects.

Evidence: data/results/rya1371/lte_reconstruction_v1, source driver scripts/audit_oxygen_lte_reconstruction.py, raw <raw>/rya1371-lte-reconstruction-v1. Archive receipts9 raw files/8 compact copies; full5.4MiB pair-depth CSV stays external under unchanged5MiB archive cap. State Registerv279 and SEQUENCE updated; scientific HOLD retained. Next: controlled joint thermal contract and integrated blackbody BB/BF rate checks across the actual quadrature, including collision inverse-rate semantics, before any fresh coupled solve or solar/Kitt Peak/ORCH validation. Preserve all-band goal, molecular coverage, flux and CNO→Fe plans.


## 2026-11-12 — Integrated thermal rate canaries

Continuation milestone: integrated thermal BB/BF canaries (Mac only)

Branch ryandamienschmitt/rya-1371-oxygen-reference-provider, worktree /Users/ryanschmitt/codex/rya1371, startinge082e7294b35a7747100eb2938eb6c05b5684b5a. Fresh origin/main remains45901e76. Source binaries, source atom, source energies, background, collisions, populations and physical constants unchanged; no installs/downloads, formal solve, field injection or SE iteration.

Independent integral audit constructs the qualified detailed-static O context at the existing56-depth native state with exact source width law and saved own wavelength grid. Uses actual compiled phi (five native angular rays/two directions), continuum cross sections, trapezoid wavelength measures, BB lambda0/lambda and BF1/lambda rate measures; common Einstein/profile-normalization/depth factors cancel in each normalized pair defect. This is not absolute rate validation or direct comparison with the compiled rate accumulator. All256 lines and22continua are accounted for. Context n/nStar and exact wavelength identity preserved. Held compiled backend and qualified source headers (Constants.hpp,LwTransition.hpp,LwAtom.hpp,FormalScalar.cpp,LwMiddleLayer.pyx) are checked against the existing ghost-operator build receipt; source native atom_mod/rt_patch_mod/const_mod and input hashes also retained and rechecked.

Results (maximum absolute normalized up/down defect over the actual quadrature/depths):
- BB constant line-centre blackbody/consistent LTE identity:3.6950441e-16.
- BB frequency-varying LW Planck field with consistent LW line-centre LTE:1.4611875e-6. Constant Einstein ratios across a finite CRD profile do not impose exact detailed balance for a frequency-varying Planck field.
- BB LW Planck with imported native LTE:0.000225481761.
- BB HOP Planck with imported native LTE and LW Einstein coefficient convention:9.3515624e-6 (still mixed, not adopted).
- BF consistent LW per-frequency Planck/exponential:2.8547076e-16.
- BF HOP Planck with LW per-frequency exponential/prefactor:0.000233076453.

Important source finding: native rt_patch_mod::bound_bound_rates and atom_mod::line_rate_coefficients explicitly split a thermal integral and deviation from Planck. At J=B, dev=0 and rates use the imported native LTE ratio, so the thermal contribution balances irrespective of finite profile, quadrature and rounding of the thermal integral. Source-formula canaries using the actual integrated thermal contribution close at2.1048361e-16, or2.1245922e-16 after FP32 rounding. This is a source-derived formula check on the own profile measure, not an additional native RT run. It explains why importing populations/constants into a different line-rate formulation is not automatically equivalent. No unique attribution of the actual NLTE residual is established and neither formulation is declared physically superior here.

Commands: source config/hosts/mac.env; OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 held LW Python scripts/audit_oxygen_integrated_thermal_rates.py --plan config/experiments/rya1371_hop_physical_shared16.json --replay /Users/ryanschmitt/codex/bride_mac_checks/rya1371-hop-physical-shared16-v1 --output /Users/ryanschmitt/codex/bride_mac_checks/rya1371-integrated-thermal-v3. Final result~1.3s (exact duration in report). v1 initial canaries completed; v2 adds native split identities; v3 adds transitive compiled-source receipt qualification. All v2/v3 numerical metrics exactly identical. No failed attempts or new bug-ledger defects this continuation; superseded receipts and driver snapshots preserved. Existing LW no-config warning/stdout trace preserved in stdout_stderr.log. No guard relaxation.

Verification: tests/test_*rya1371.py76 passed in2.05s;5 targeted thermal tests passed. Scientific self-review checks units, transition graph, exact context-state preservation, source/profile/kernel provenance and distinction between analytical identity, independent quadrature, compiled accumulator and actual coupled RT. Integrated canary test deliberately broadens occupation variation to discriminate finite-profile effect; native split test includes FP32-rounded thermal integrals. Data at data/results/rya1371/integrated_thermal_v1, final raw rya1371-integrated-thermal-v3. New archive11 raw receipts/5 compact copies; no5MiB cap change.

Scientific HOLD retained: actual matched-population NLTE residual remains0.000221166; no improvement claimed from these canaries. All-band capability/instrument coverage, true flux, canonical solar atmosphere, observed Kitt Peak/ORCH, molecular opacity, parameter sweeps and3D gates remain open. Next concrete task: independently reconstruct absolute own BB/BF rates from the saved J field and actual profiles, verify against saved compiled rates, then isolate the native thermal/deviation formulation at that fixed field before any opt-in coupled operator change. Joint population/background/collision contract still required; no LTE-only or single-constant patch. CNO then Fe and RYA-1374 remain sequenced behind Oxygen qualification.

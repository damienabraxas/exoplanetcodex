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

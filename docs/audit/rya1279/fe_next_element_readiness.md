# Fe handoff after the Oxygen foundation

Ryan's 2026-11-14 direction makes Fe the next Bride expansion target. C/N remain
queued. This is readiness work alongside Oxygen qualification; no Fe RT/SE run,
adopted abundance or production Fe disposition changed. Existing ownership is
RYA-1371 for the reference adapter/Oxygen foundation and RYA-1372/1368 for the
subsequent cross-element comparison. RYA-1137/1139 are historical atom/line-pool
and compression precedents, whose later comments supersede their original briefs.

## What was checked now

The supplied, pinned OpenMULTI3D `atom.fe607a` was parsed in full. It contains
607 levels (548 Fe I, 58 Fe II, one Fe III), 12,635 bound-bound transitions **all
within Fe I**, 606 continua and 29,983 collision records. Its declared abundance
is **7.50**, an input to match in a controlled experiment, not a solar measurement.
The transition centres span approximately 170.08–1999.91 nm vacuum; there are no
bound-bound transitions at or beyond 2 microns in this source atom. Fe II levels
do not establish Fe II bound-bound NLTE capability.

Full `AtomicModel` construction and a dry wavelength union completed, without
allocating a Context, atmosphere, profiles or rate matrices. Under the inherited
Oxygen diagnostic sampling prescription, the candidate union has **242,487
points**, exceeding the current **20,000-point** HOP packet scope. Its projected
BB profile array alone is **2,321,285,120 bytes** at 56 depths, five rays and two
directions; a single dense 607-square matrix at those depths is 165,065,152 bytes.
These are arithmetic projections from the dry grid, not measured peak memory or
runtime. Full Fe needs a separate resource design; reusing the O41 allocation is insufficient.

The supplied atom's SHA differs from the earlier RYA-1055 legacy atom receipt.
The registered local Gerber directory is absent, so physical equivalence has not
been established. Different bytes alone do not prove different physics. No
unregistered atom or unavailable grid was substituted. Source attribution/rights
still need explicit audit; the solver's code license does not settle atomic-data
redistribution rights.

Machine evidence: `data/results/rya1371/uv_background_fe_v1/rya1371-fe-readiness-v1/report.json`.

## Execution ladder

1. **Close the Oxygen numerical foundation.** Resolve the UV continuum rate
   integration uncertainty and recheck physical closure/restart with the qualified
   discretization. Retain the full source graph, source states and geometry. Fe
   metadata/coverage preparation can proceed, but this does not clear the gate for
   a Fe abundance result.

2. **Pin the Fe atom and actual assigned line pool.** Reconcile the supplied atom
   with the legacy lineage when its source is available; record data rights and
   all source differences. Remeasure both-level and actual-transition coverage on
   the chosen first-light set, starting from the held AGSS21 Fe I reference set.
   Do not substitute an ungraded canonical superset, rely on old coverage counts,
   or treat term-averaged and fine-structure levels as interchangeable. Record
   eligible/gap cases by band, including the absent Fe II lines and long-IR reach.

3. **Choose a resource-safe representation with a physics check.** Compare a
   measured, bounded native reference-column route with a frequency-windowed
   full-coupling bridge and, if needed, the existing qualified-compact-atom route.
   Windowing must accumulate all transitions and correct whole-line quadrature
   normalization across windows. A restricted/super-level atom changes coupling
   and needs a separate full-versus-restricted departure/line test. RYA-1139's
   latest scope is compression with inherited rates, not arbitrary collision
   filling or ceiling extension. No packet/memory/tolerance bypass, silent atom
   pruning, or new ab-initio data. Old FAL-C timings are historical evidence,
   not current matched-solar resource measurements.

4. **Qualify Fe-specific physics and one matched 1D column.** Generalize the
   currently O41-bound collision, width and background-ownership contracts only
   with independent Fe references. Check every TEMP block, rate direction and
   detailed-balance convention; exclude exactly the active Fe continua from
   background ownership. Match atmosphere, abundance, depth coordinates, hydrogen
   populations, velocities, radiation geometry and all numerical conventions.
   Compare populations/departure structure and independently reconstructed rates,
   then convergence/refinement and fresh restart. Keep the source-input 7.50
   distinct from any later observed abundance.

5. **Validate solar Fe I in VIS, with a separate Fe II diagnostic.** Use canonical
   solar parameters and a small preregistered clean line set with the actual
   Kitt Peak observable: AIR wavelengths, verified continuum/tellurics, matched
   disk-centre intensity or integrated flux, blends and instrument convolution.
   An Fe II ionization check requires an independently eligible model or an
   explicitly labelled LTE limitation. This atom's Fe II levels alone do not
   qualify general Fe II NLTE line synthesis. Do not force a 1D calculation to
   the familiar 3D 7.46 reference; validate matched physics and geometry.

6. **Stress parameters and expand supported bands.** Sweep temperature,
   microturbulence, gravity, metallicity/abundance, damping and continuum choices;
   compare line-by-line excitation/strength trends and Fe I/Fe II consistency.
   Extend to eligible UV/red/NIR/H diagnostics and then a known reference star.
   Track long-IR/atomic completeness and molecular backgrounds explicitly. Keep
   1D, mean-3D, 1.5D and full-3D claims separate.

Each gate needs reproducible input/output identities, quantitative comparisons,
uncertainties and an explicit PASS/HOLD. Parsing, a finite result or agreement
with one abundance is insufficient. All new findings remain experimental until
reviewed; existing Codex Fe products, frozen gold and publications are untouched.

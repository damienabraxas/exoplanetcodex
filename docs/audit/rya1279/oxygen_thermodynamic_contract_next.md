# Next Oxygen control: one declared thermodynamic contract

This is an implementation handoff under RYA-1371, not a physics adoption. Fe
remains the next element after Oxygen qualification. C/N remain queued.

The held 5750 K / logg 4.5 / A(O) 8.70 case is not the canonical Sun. Its
39,841 -> 59,003 -> 86,048 native grids preserve every inherited node. The
second population step is about nine times smaller than the first; optical
changes are small. This supports numerical stabilization but does not certify
full rate/population convergence or observed-solar accuracy.

## What the new evidence distinguishes

The remaining bridge peak is excited O II level 41 at the top cell, connected
by source continuum 4 -> 41 over **94.3253–106.4179 nm vacuum**. It remains near
2.0e-4 after the second native refinement. On the same saved bridge field,
native BF thermal coefficients reduce this cell's residual to about 1.0e-6.
Changing cross-section interpolation alone barely changes this particular cell.

These constant profiles were inspected directly:

| Profile | h (J s) | k (J/K) | Owner |
|---|---:|---:|---|
| Held Lightweaver | 6.6260755e-34 | 1.380658e-23 | Pinned `lightweaver/constants.py` |
| Native HOP | 6.626176e-34 | 1.380662e-23 | Pinned `HOP/const_mod.f90`, converted from cgs |

The legacy `modern_thermal` diagnostic keys denote the **held Lightweaver**
profile, not current SI values. Earlier immutable reports are retained; their
terminology is clarified separately. Neither package nor either profile has
been modified.

The six combinations of three alpha choices (LW WENO wavelength, native
Hermite frequency on raw tables, and native Hermite on consumed tables) and two
BF thermal profiles are all reported in `grid-refined-factorial-v1`. No variant
was selected by reference closeness. Native-alpha plus native-thermal shadow
matrices actually worsen the global maximum: a rate-only substitution at an
unchanged radiation field is not a self-consistent model. Measured interaction
terms are reported rather than silently adding separate effects.

## Concrete execution ladder

1. **Declare coefficient ownership.** Trace the background chi/eta/scattering
   packet in `NativeHopBackground`, the active continuum alpha implementation in
   `MultiExplicitContinuum`, compiled radiative thermal/Einstein coefficients,
   imported LTE populations and the collision binding's reverse rates. Keep
   interpolation choice, constant profile and data-table rounding as separate
   axes. Source forward rates/data must not be adjusted to reproduce an abundance.

2. **Test a consistent alpha control first.** An explicitly experimental native
   frequency-Hermite alpha route must affect both active opacity/emissivity and
   the rate integrals at the same fixed populations. Validate it against the
   compiled native oracle, including table endpoints, units, support, positivity,
   all 22 continua and source-frame identities. Preserve the existing default.
   Recompute the radiation field; do not inject a native J or interpret a shadow
   rate matrix as a new RT solution.

3. **Specify the complete thermal profile before a coupled change.** A candidate
   profile must define background emissivity/Kirchhoff behavior, active BB/BF
   thermal factors and Einstein relations, LTE/partition/energy conventions,
   collider populations and collision detailed balance together. Tests must
   distinguish thermodynamic equilibrium from an ordinary non-LTE solution.
   Changing Python constants alone does not change a compiled backend. Importing
   foreign LTE populations or a reverse-rate matrix does not establish a common
   convention. Resolve implementation authority through the existing science
   standards and input owners; do not choose a constant set by the smallest fit.

4. **Qualify the joint operator.** Require independent physical rate and particle
   conservation checks, source/own fixed-state comparisons, numerical refinement,
   a fresh coupled solve and restart. Preserve every level and the top cell,
   including the weak excited-ion population. Keep the original one-thread,
   600-second, 2 GB and HOP 20,000-point constraints; frequency tiling or changed
   quadrature needs its own normalization/coupling validation.

5. **Return to the solar observable.** Only after the operator is qualified,
   compare canonical-solar spectra using the actual Kitt Peak observable, air
   wavelengths where supported, intensity/flux geometry, continuum, tellurics,
   blends and instrumental response. Retain parameter sensitivities and the
   broader UV–VIS–IR/instrument-coverage gates. Internal FUV rate diagnostics are
   not observed FUV validation.

Step 2 now has a qualified fixed-population diagnostic: all22 continua, fresh
radiation, exact default reproduction and independent rates (joint_alpha_v1).
The residual remains2.00412e-4 versus default2.00421e-4; interpolation alone does
not resolve the dominant cell. The next work is the complete thermal-profile
contract in step 3, followed by joint operator qualification. The thermal profile work is a contract to establish and test,
not authorization to overwrite constants, defaults, source tables, gold or
production products. Scientific admission remains HOLD.

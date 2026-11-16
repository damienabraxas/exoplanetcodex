# Full oxygen LTE population and inverse-input milestone

RYA-1371 continuation from53055422. Explicitly read-only prerequisite for a
common active thermal convention; no coupled solve or physical admission.
Earlier source reconstruction was same-stage only. This audit covers all41
levels, both OI/OII stages and56depths, including the disconnected excited OII
level41. Source `atom_mod.lte_pop` consumes FP32 e_kT, lsaha and Debye lowering,
with FP64 consumed level energies/statistical weights and exponentiation.
`composition_mod.set_aux_data` owns lsaha and lowering. Source constants and
actual DEBYE_SHIELD/LEVEL_PARTF flags are pinned. No parameter is inferred by
fitting nStar. The1e-12 arithmetic screen is inherited from the prior same-stage
audit and declared before this execution; it is not a physical error budget.

Independently reconstruct all normalized source populations and ground ratios;
compare a separate held-LW log formula against actual held `lte_pops`. Retain
native FP32, native FP64, no-shield negative control and held-LW alternatives.
Keep the actual total oxygen density at every depth. This qualifies finite-atom
shape ONLY: absolute element density, external partition-function normalization
(LEVEL_PARTF=F), chemistry and EOS remain open. No substituted LTE population
is supplied to a solver.

Audit fixed collision reciprocity under each diagnostic population, and a
clearly algebraic reverse-ratio projection preserving source upward coefficients.
The latter tests ownership coupling, not independent microscopic inverse physics
and never becomes solver input. No missing level41 collisions are invented.
All prior files, model/production defaults, observations/gold, source tables,
resource constraints and scientific HOLD remain unchanged.

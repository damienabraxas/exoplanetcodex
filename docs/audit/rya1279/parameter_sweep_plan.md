# Bride solar parameter exploration: staged sweep

Ryan explicitly requested a broader sweep on2026-11-03. The0.9/1.0/1.1km/s
microturbulence slice is the first local check, not the complete promised sweep.
Every subsequent slice must state its actual inputs, changed axes and coverage;
there is no automatic claim that earlier controls used the current physics.
All science stays on the Mac under the existing bounded, serial authorization.
Scientific HOLD; fixedO8.70 baseline is not a measured solar abundance.

## Current executed slice

`config/experiments/rya1281_mac_xi_slice.json` and `rya1281_mac_xi_wide.json` declare named cases,
fixed solar5772/g4.438/full56-depth thermodynamic structure, O8.70, bound FeI
continuum,41 source FeI pressure lines and Cowley Stark, native atomic H/EOS.
H2 pressure-only result is separate; full chemistry/continuum/O collisions are
not admitted, so H2 stays omitted in this coupled slice.

Each requested xi has a fresh native producer evaluation, output-xi check and
complete certificate. Cross-case and one-ULP changes refuse. Adding more xi
amplitudes requires an explicit plan and fresh inputs, not a runner rewrite or
certificate relabel. The current implementation refuses changing another axis
inside a microturbulence-only phase. Future axes need their own validated input
contracts. A malformed or out-of-coverage case fails rather than being dropped.

## Broader sweep inventory

| Axis or control | Evidence already present | Next coverage / dependencies |
|---|---|---|
| Microturbulence | Earlier native-node0.5/0.9/1/1.1/1.5 controls use older physics | Seven fresh distinct solar points0.5/0.75/0.9/1/1.1/1.25/1.5 completed (xi_wide_control_review.md). Representative numerical refinements and interactions remain owed; no xi fitted to Kitt. |
| Effective temperature | Five fresh solar5752-5792K at10K spacing with bound continuum; sparse Fe-line widths added later | Wider source-matched slice owed. Fresh matched5572/5672/5772/5822 cases converge (thermal_control_review.md);5872/5972 remain refused at976tables because deep T exceeds10000K. Resolve source coverage, no extrapolation.10K sampling is not10K physical accuracy. |
| Surface gravity | Earlier native-node checks, not fresh source-matched solar series | Fresh matchedg4.338/4.438/4.538 cases converge with complete bindings;4.388/4.488 remain domain-only preflight. Composition/interactions and model accuracy remain owed. |
| Composition: metallicity, alpha and individual CNO/Fe | Native composition explicit; O8.70 fixed; earlier diagnostic O perturbations | Fresh producer/EOS/opacity inputs for each declared composition. Distinguish active-O response from atmosphere mixture response. Passive Fe for oxygen is not the Fe abundance milestone. |
| H/electron reservoirs and chemistry | Native charge/EOS audits, full-domain H2 source and partial collider check | Full H/Hminus/H2, electron/charge and metal/molecular consistency; separate source prescriptions, never arbitrary ne/rho scalings promoted as physical atmospheres. |
| Atomic data and collisions | Source tables/conventions and numerous numerical controls preserved | Source-backed alternatives and uncertainty inputs; no fitted rate/gf multiplier chosen to match Kitt. Independent matched1D oxygen benchmark remains indispensable. |
| Broader line/molecular opacity |41 FeI lines; combined components and other species omitted | Source-qualified components, species and later CNO molecular opacity; preserve accounting and fixed-physics counterpart. |
| Continuum/tellurics/observation | Shipped Kitt baseline plus ORCH envelope/synthesis/IAG reference legs; reference/anchor sensitivity remains | Qualify continuum geometry, masks and independent-reference spread. These are conditioning/reference choices, not free parameters to optimize. No second telluric correction. |
| Broadening/velocity/instrument response | Canonical and zero/gravitational RV variants retained | Source/uncertainty-backed vmac/vsini/FTS-LSF controls; wavelength frame fixed by source, not a fitted vacuum/air choice. |
| Numerical resolution and boundaries | Finite cap, wavelength/angular/depth and convergence controls exist at specified earlier physics | Combined nominal rays3->5/qscale1->2/flux5->7 converges; O-area differences<=0.002881mA. Matched5572/5822/g4.538 refinements pass with max area0.0048073mA; g4.338 hits original600s limit, not accepted (thermal_extreme_refinement.md). Xi extremes/individual contributions remain owed; numerical checks do not validate physical boundaries. |
| Model family/dimensionality | Experimental static1D solar O777 | Independent matched1D reference, atmosphere-family spread, then1.5D/full3D with actual validated grids/backend. No claim1D parameter sweeps establish3D NLTE. |

Begin with one axis at a time and a matching nominal control. After those
responses and coverage are understood, select combined parameter points to test
interactions (e.g. Teff-xi and gravity-composition); do not assume separability.
Predeclare those combinations before looking for a favorable residual. Preserve
failed/out-of-domain cases, not just converged or visually attractive spectra.
Record timing/memory for actual cases before scheduling any larger lattice.

The immediate species sequence remains oxygen, C/N (CNO), then Fe abundance,
then a reference star with verified spectra/parameters/abundance. Additional
oxygen lines and bands need appropriate atoms, blends and observations; O777
VIS controls are not full-spectrum/all-band validation. HARPS does not cover
O777; any HARPS oxygen test must use a separately qualified overlapping line.

Full-profile finite geometry now separates xi response more than three areas do; T/g remain similar. Both fixed RV variants checked; source-backed broadening/continuum/abundance nuisance and noise controls are required before any inverse interpretation. Next numerical work must instrument actual Ng flags and compare declared controls for the unqualified lower-gravity endpoint.

Native Ng control (ng_control_review.md): exact221-row on replay, actual43Ng updates; disabled control smoother but both original600s attempts refuse. Lower-gravity refinement remains unqualified. Next cost/acceleration controls must retain all numerical accuracy gates; no universal off policy or accepted spectrum inferred.

Cost-only replay (ng_extended_review.md) closes lower-g combined refinement:233steps/630.71s with1200s declared allowance, unchanged accuracy and exact221-row prefix. All four declared thermal/gravity endpoints plus nominal now checked, max area0.0048073mA. Original600s failures remain. Two refined continuum legs retain numerical/geometry identity and4.93-10.88mA reference spread, so independent continuum/opacity/LSF/reference work is next; no broader-grid/physical-error inference.

Frozen continuum support audit (continuum_anchor_review.md): denominator changes remain6-11mA with anchor membership fixed; one-sided synthesis7771/7774 extrapolate>4.10A. Native7774 fifth-pass coefficient/support-count mismatch is characterized, with~-0.927mA counterfactual refit; no estimator adoption. Next bounded implementation: expose coefficient-support/survivor counts and cap status as separate audit provenance while proving byte-identical native continuua, then source-backed reference/FTS-LSF and opacity checks. Broader parameter/interactions and O->C/N->Fe->reference-star sequence remain owed.

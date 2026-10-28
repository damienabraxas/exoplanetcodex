# Adaptive Bride direction: experimental roadmap, not a certified capability

Ryan reiterated the destination during the Mac controls: an adaptive full-3D NLTE
model that can ingest useful atmosphere grids, calculate spectra across bands,
and eventually generate grids of its own. Existing packages are components, not
the limit of the design. This proposal records that direction; it does not
ratify a new treatment token or claim that the current1D tests establish3D.

## The contract is physics and provenance, not a particular file format

An atmosphere adapter should provide explicit geometry/coordinates, depth and
horizontal topology, units, T/rho/ne/nH, velocity vectors, chemical composition,
EOS/opacity lineage, boundaries and time/snapshot identity. Its verified domain
and missing fields must be machine-readable. MARCS, STAGGER and another simulation
family should enter through their own validated adapters. Unsupported geometry,
rates or chemistry cannot be made equivalent by reshaping arrays or silently
inventing thermodynamics. An independent snapshot/mesh and unit canary belongs
to each adapter, including non-square multidimensional cases.

An atomic/molecular bundle separately owns levels, radiative/collisional rates,
cross-sections, broadening, temperature ranges, licenses and transition reach.
The current41-level O atom is one experimental bundle; it is not a complete
CNO/molecular/full-spectrum background model. Its10000K collision-table ceiling
already exposed a real domain boundary in the6000K MARCS atmosphere.

The backend must declare actual capabilities:1D, independent-column1.5D,2D, or
horizontally coupled full3D. Current Lightweaver provides1D/2D building blocks;
its layout API is not proof of a full3D solver. A future3D implementation needs
an explicit backend choice or new formal-transfer/SE implementation, then
analytic limits, homogeneous/horizontal-coupling tests, independent published
benchmarks and conservation/error/resource checks. RYA-1282 remains the owner
of that next dimensionality gate. No horizontally averaged structure or solar
correction table may acquire a full3D label merely through an adapter.

## Original layer to build: error-controlled adaptive refinement

The useful new component is a controller that decides when available grids are
insufficient, rather than relying on increasingly dense interpolation:

1. Record the distinction between native/solved, interpolated and corrected
   nodes, their source hashes and their model/chemistry/geometry domains.
2. Hold out native nodes and predict them. Compare thermodynamic fields as well
   as NLTE populations and spectra; a smooth spectrum alone is insufficient.
3. Use structure-preserving coordinates: positive pressure/density and positive
   geometric layer thickness are explicit candidate transforms. Reject extrapolation,
   composition changes without an owner, nonpositive fields and inverted meshes.
4. Measure parameter response and interpolation error separately.10K output
   sampling is not a10K physical-accuracy statement. Current experiments compare
   linear fields, log-positive fields and log-positive layer-thickness coordinates
   against native atmospheres, independently of Kitt Peak fits.
5. Refine where held-out error, curvature or a physical domain boundary exceeds
   the applicable predeclared scientific budget. Choose new calculations to reduce
   that error, not to drive an observed abundance toward a familiar answer.
6. Reuse verified nodes and bounded job manifests; record failed or unsupported
   cells explicitly. Benchmark wall time, memory and scaling before scheduling
   large parameter/snapshot/wavelength campaigns on the incoming hardware.

The first concrete demonstration is the Mac interpolation check: a basic gravity
interpolation can give large geometric thickness errors despite a close temperature
profile. The positive-layer prototype is original engineering under test, not
newly established atmospheric physics or an EOS solver. Its independent holdout
evidence determines what can be trusted next.

## What generating our own grids would mean

- **Spectra/departure grids on existing atmospheres:** run an independently
  validated NLTE solver at selected Teff/logg/composition/xi/abundance nodes,
  preserve convergence/provenance, and refine adaptively. This is the nearest
  feasible grid-generation milestone; it does not require pretending interpolated
  atmospheres are newly solved structures.
- **New atmosphere grids:** require an atmosphere-equilibrium or radiation-
  hydrodynamics code plus qualified EOS/opacities/chemistry and reuse rights.
  More RT samples alone do not solve a new atmosphere. A matched solar structure
  is a prerequisite before an experimental spectrum becomes a solar-abundance fit.
- **New3D time-dependent grids:** add simulation volumes, boundaries, temporal
  relaxation/snapshot representativeness and horizontal transfer. This requires
  resource and independent physical validation, beyond a larger processor count.

An expanded compute budget makes these campaigns possible; it cannot supply missing
atomic rates, opacity coverage, licensing or independent validation. The controller
should say which resource is limiting each cell instead of silently substituting
the nearest model or a different physics regime.

## Full-spectrum and observation boundary

All-band coverage is a sequence of justified capabilities, not a switch. The coupled
NLTE solution needs the relevant UV/optical/IR radiation field; the emergent spectrum
also needs qualified weak blends, molecular chemistry/opacity and continuum data.
Atomic777nm validation does not certify OH/CO/CN/NH or CRIRES H/K. Each spectral
window needs an actual held product, coverage/quality/error evidence and supported
line/atom/molecule physics. Existing RYA-1214/1220/587 retain their owners.

Reuse orchestrator holding policies: corrected/normalized products are frozen and
served unchanged, not reprocessed for each model. Preserve instrument/holding
identity, LSF, observation geometry, wavelength/velocity conventions and covariance.
The current777nm comparison is air Angstrom; Kitt Peak2005 is delivered in vacuum
nm and converted once. RT uses vacuum frequencies. FUV below2000Angstrom stays
vacuum under the canonical converter because an air wavelength is undefined there.
Earth's telluric absorption and solar gravitational/convective shifts are distinct
effects, not alternate meanings of air/vacuum.

## External stellar benchmark after the solar controls

Select a known benchmark with a verified Gaia crossmatch, independent Teff/logg/
composition and an independently sourced oxygen abundance with stated treatment
and uncertainty. Gaia identity or metallicity alone is not an oxygen reference.
Audit the actual held spectrum at the intended transitions; HARPS cannot test777nm.
Prefer existing qualified holdings before acquiring data. Keep flux vs intensity,
air/vacuum, RV, resolution, rotation, telluric and continuum identity explicit.

Freeze that benchmark's inputs before running the model. Use the solar-tested
paired parameter responses and independent line/spectral residuals; do not tune
the star's parameters or reference abundance to make the Bride pass. A second
temperature/gravity/composition regime then tests whether agreement generalizes.

The immediate next steps remain a matched solar-atmosphere/interpolation check,
resolution of the O model/MULTI background/EOS discrepancy, then a qualified
stellar benchmark. Full3D and own-grid generation remain explicit subsequent
engineering/science gates, with no merge or publication implied by this roadmap.

## Solar element sequence (Ryan, 2026-10-28)

Finish oxygen, then carbon and nitrogen: CNO is the first multi-element milestone.
Iron follows after shared reaction/rate provenance, population benchmarks and
solar observation checks are established across CNO. Another-star validation
follows the solar element work. This order does not imply C/N/Fe models are
already ingested or scientifically admitted.

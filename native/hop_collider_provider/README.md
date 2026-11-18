# Experimental native Hydrogen collision-field export

GPL-3.0-only derivative of OpenMULTI3D at
`32085d26c46e9919a986397b5514a9728abe2074`, retaining the trace-continuum
ownership patch and upstream AUTHORS/COPYING. No coefficients change.

`hop_get_h_ground` reads the exact `levels(1)%n` used by native Oxygen collision
`nh1` and separately exports `levels(1)%w`. With the pinned occupation-probability
recipe, `init_H_levels` leaves the occupation factor outside stored populations.
Thus the raw collision field is not a normalized probability partition of H I;
only its product with the occupation weight is bounded by total H I. Neither
field is rescaled, clamped or substituted for total neutral-H pressure broadening.
H II and H-minus use the existing separately named species export.

The function refuses unloaded/incorrect-size calls and missing level/weight storage
before writing either output. This derived interface deliberately requires the
explicit occupation-weighted Hydrogen recipe, not other unqualified EOS modes.
It is read-only and grants no coupled-SE, collision-physics or solar admission.
Build using `scripts/build_hop_trace_provider.py --patch` with fresh source/output
paths. Generated headers are included in the canonical patch; original source,
trace provider and binaries remain intact.

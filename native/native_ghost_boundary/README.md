# Native linear zero-ghost geometry control

Source: OpenMULTI3D at `32085d26c46e9919a986397b5514a9728abe2074`,
`SRC/rt_patch_mod.f90::set_upstream_boundary` and
`SRC/rt_solvers_mod.f90::long_char_linear/linear_weights`.

The incoming ghost plane has zero extinction, source and intensity. The sampled
top node is one native mesh spacing below it. Thus

```
t = chi_node * native_step / (2 * abs(mu))
I_node = (eta_node / chi_node) * [1 - (1 - exp(-t))/t]
```

This differs from putting zero incoming intensity at the sampled top node. The
spacing is read from the actual consumed native mesh, with all node coordinates
verified; it is never fitted. The fixed-extinction emissivity derivative is
`[1-(1-exp(-t))/t]/chi_node`, with zero-extinction limit `step/(4*abs(mu))`.
Small-t evaluation uses the Taylor series. Python boundary and compiled numerical
operator are checked against independent formal integration and actual responses.

The isolated operator retains copied Lightweaver header notices, preserves the
installed source/backend, and changes preconditioning only. Both physical-action
identity and fresh R+C/conservation/restart tests are required. The geometry
control does not itself grant physical, solar, all-band or production admission.

Native saved radiative coefficients contain SE preconditioning corrections.
They are not gross Einstein rates; individual up/down ratios against direct LW
rates are refused. Saved net-rate balance is distinct from a fresh unpreconditioned
native radiation/rate replay, which remains a separate qualification step.

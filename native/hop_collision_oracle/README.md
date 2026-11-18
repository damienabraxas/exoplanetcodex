# Native collision interpolation oracle

GPL-3.0-only driver of OpenMULTI3D HOP `interpol_mod` and `const_mod`, pinned at
`32085d26c46e9919a986397b5514a9728abe2074`. Upstream AUTHORS/COPYING are retained
in the registered source checkout. The wrapper adds no interpolation formula.

`build_hop_collision_oracle.py` links the existing qualified collider provider's
native objects and receipts every object/module/source, wrapper and command in a
fresh output directory. No native source, object or original executable changes.
Input has bounded record/node/query counts, increasing temperature nodes, positive
coefficients and in-domain queries. Every coefficient is interpolated in log space,
matching the supported CE/CI/CH/CH0 native table reader. The native `e_k` constant
is exported explicitly. This is a numerical component oracle, not independent
physical collision evidence or a coupled statistical-equilibrium solution.

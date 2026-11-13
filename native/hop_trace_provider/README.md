# Experimental HOP trace-continuum interface

This GPL-3.0-only patch derives from OpenMULTI3D's GPL application at
`32085d26c46e9919a986397b5514a9728abe2074`. Upstream authors and attribution are
in that application's `AUTHORS`; `COPYING` preserves its license text.

It exposes the non-hydrogen atomic selection in
`SRC/nlte_patch_mod.f90::exclude_double_continua` through a C ABI before loading
an EOS state. Atomic species, mass, stage and lower-level excitation matching
follow the native method. It changes active background ownership flags, not
cross sections, rates, abundances or EOS coefficients. Hydrogen and molecular
trace-model special cases are explicitly refused. No production admission.

`scripts/build_hop_trace_provider.py` applies it only in a fresh detached
worktree, builds one thread with the original compiler/options, and receipts
the upstream pin, exact patch/diff, generated header and shared library. The
original reference checkout and binaries remain unchanged.

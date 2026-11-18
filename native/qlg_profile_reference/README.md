# Isolated active-profile reference control

GPL-3.0-only derivative of pinned OpenMULTI3D application commit
32085d26c46e9919a986397b5514a9728abe2074. License text is retained at
[../hop_trace_provider/COPYING](../hop_trace_provider/COPYING); original authorship
remains in the upstream application AUTHORS and QLG module.

A new QLG helper retains qlg_bot_profile's fp32 input/offset arithmetic, final
fp32 conversion, normalization and frequency Jacobian, calling the existing fp64
reference kernel in the middle. Only the three active non-H RT call sites select
this helper. Existing QLG routines, standalone background/linelist behavior,
atomic/EOS inputs and the original source/executable remain unchanged.

This is an opt-in diagnostic, not an accurate-model/default adoption. Its native
population response, rates and emergent spectra must be measured independently;
fixed-population equality is never assumed. Fresh private source/build receipts
and a one-thread600s budget are required. No production model registry entry.

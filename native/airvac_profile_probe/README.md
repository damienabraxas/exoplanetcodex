# Native wavelength encoding probe

GPL-3.0-only typed Fortran driver of the pinned OpenMULTI3D handyfuncs_mod
functions, linked against the verified original HOP library and compiler module.
License text: [../hop_trace_provider/COPYING](../hop_trace_provider/COPYING).
No vacuum/air formula is copied or introduced into the project converter.

The native internal table encoding is distinguished from the canonical observation
coordinate, which remains vacuum below2000A. A numerical internal roundtrip does
not certify an atmospheric-air interpretation, FUV synthesis, rates, or all-band
model applicability. Original source, library, atom and binary remain unchanged.

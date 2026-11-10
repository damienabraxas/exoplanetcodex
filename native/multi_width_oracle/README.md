# Original MULTI width oracle

Standalone O-only recording driver of the exact `DAMP` routine at pymulti commit
4d6aa3c39d634b38da772e469a5692435eba1354. Source/includes are reconstructed from
verified external pinned bytes, not vendored here. No upstream code/data rights
are inferred from this driver or another solver's license.

The build follows the original DP recipe: PREC.dp and REAL→REAL(8) include
conversion, unchanged source and BCONST literals. Only capacities increase:
MK30→64, MLINE150→512. WDAMP records original computed components; no width
formula is replaced. Oxygen only; no RT/EOS/collision/SE calculation.

Native DAMP sums NH(1:5) for neutral-H broadening. The driver stores the supplied
aggregate only to service that reduction; it does not supply fictitious H level
populations. Hydrogen-specific branches are unreachable (ATOMID hard O). vturb
is explicitly m/s input mapped to the original km/s/QNORM storage. AWGT is amu
mapped through the original UU, energies cm^-1 through the original HCE.

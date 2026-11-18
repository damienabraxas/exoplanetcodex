# Source-derived line-frequency-grid oracle

Read-only caller of pinned OpenMULTI3D HOP `trans_mod::line_type%get_lam`,
source32085d26c46e9919a986397b5514a9728abe2074. Upstream GPL notices preserved;
driver SPDX GPL-3.0-only. Links existing qualified collider-provider HOP objects in
a fresh owned output directory; verifies original and linked trans_mod source bytes
are identical. No source/build overwrite or global/default installation.

Consumes ALL256 held-source transitions with actual NQFAC2/QNORM10/legacyF controls.
Uses lam0=1A to export dimensionless source offsets, then applies those offsets to
ORIGINAL source photon centres. This tests numerical support/spacing independently
of native table photon-centre/loggf rounding. No E/g/f/population/C/width override.

The support-only control stretches the original doubled LW node shape to exact
source blue/red endpoints; the full control replaces both support and spacing.
Neither result certifies physical accuracy or numerical convergence. In fact the
apparently closer coarse native-grid rate comparison deteriorates under refinement;
see v276 audit. Generic coupled runners/defaults do not request these controls.

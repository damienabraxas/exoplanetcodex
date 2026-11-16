# Oxygen emission and diffusion-boundary profile axis

RYA-1371 continuation from c603edc55ee74df23bfab8bd55b40e7afb7ae063.
Worktree `codex/rya1371`, branch `ryandamienschmitt/rya-1371-oxygen-reference-provider`.
Refreshed origin/main remains 0523151e. Draft PR620, scientific HOLD.

## Predeclared bounded contract

The next fixed-population diagnostic selects either held-LW or source-native
Planck emission for BOTH background thermal absorption/emission and bottom
thermal diffusion inputs. Both alternatives must be executed and retained;
neither is selected by smaller residual. Background chi/scattering, colliders,
active BB/BF coefficients and rates, imported n/nStar, geometry, wavelengths,
upper ghost, resource limits and tolerances remain unchanged. The existing
native-Hermite alpha control is retained in every comparison.

`EmissionProfileBoundary` uses a Python callable with the exact linear-solver
bottom condition inspected in held `FormalScalar.cpp`: B_bottom minus the
upstream-to-bottom B difference divided by trapezoidal optical depth along the
ray. Total bottom opacity is captured from actual depthData after the first
iteration and required invariant thereafter. First capture starts at B_bottom;
convergence includes the boundary refresh and may not stop at that initialization.
Wavelength/state extrapolation and negative/invalid coefficients/intensities
are refused, never clipped. Backend/default files and constants are unchanged.
The source inspection is distinct from a source-to-installed-binary attestation;
a disposable continuum context checks the callable against the actual compiled
ThermalisedRadiation boundary, including gradient and angular dependence.

Before the held source case: test invalid contracts/arrays and compiled-boundary
reproduction. Then fresh default reproduction against prior native-alpha replay,
followed by BOTH profiles at the exact inherited 19,963-node state. Preserve all
failed attempts, input/helper snapshots, source/build receipts and raw outputs.
Use existing one-thread 600s/20,000-node constraints, Mac only, no SE update.
Independently reconstruct all source radiative rates from each new saved J.
Report total-opacity/scattering/population/collision identity separately from
changed background eta, lower intensity, J and rates. No observation, abundance,
canonical-solar, UV/fullband/3D, coupled closure or adoption claim follows.

## Intake and verification

Repository task-intake, implementation-brief, scientific-code-review and
state-register procedures apply. Existing Solar Beta refinement report has
10 ticketed unresolved entries and five without resolving tickets; none is
cleared by this diagnostic. O model admission remains open. No new library
acquisition or physics source needed: held source/backend/atom receipts own
this control. Full source and frozen products remain protected.

Initial disposable test attempt failed before Context construction because a
passive-only RadiativeSet cannot generate its wavelength configuration. The
fixture now generates the configuration while H is detailed-static, then marks
H passive before construction (the existing frozen-total-provider pattern).
No native/solar run occurred on that fixture failure. No pipeline defect or
source-data change discovered at this point.

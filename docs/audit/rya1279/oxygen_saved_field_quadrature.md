# Joint oxygen saved-field quadrature budget

RYA-1371 continuation from00579232, branch
ryandamienschmitt/rya-1371-oxygen-reference-provider, worktreecodex/rya1371,
main4cdc449a, draftPR620. Scientific HOLD; latest Linear handoff unchanged.
Task-intake/implementation-brief/scientific-review/state-register procedures
continue; no new acquisition, model or observational production run.

## Predeclared bounded milestone

Audit BOTH existing joint17263/19963-node radiation fields. Preserve all41levels,
56depths,22BF/256BB, physical populations/LTE/C and source alpha/thermal factors.
Use physical angle-weighted J_nu linearly interpolated in VACUUM wavelength;
this is a declared saved-field representation, not new radiation evaluations.
Integrate its22BF rates at Gaussian orders4/8/16 on panels split at every saved
J node and raw source-alpha knot. Actual qualified compiled frequency-Hermite
FP32 alpha is queried, not replaced by a Python interpolation. Preserve actual
angular-weight sum in spontaneous terms. Native alpha queries remain bounded
by20000; stream256panels per request (maximum4096nodes). No larger HOP/RT/SE
packet, field/population update, backend change or guard bypass. Original
one-thread600s deadline retained; raw requests private/receipt-only.

Before Gaussian comparisons, source-alpha/explicit wavelength trapezoid must
reconstruct every saved compiled BF directional rate within the existing
16*maxnodes*5rays*2directions FP64 positive-summation screen. Every BF/depth
included. Independent tests verify positive Gaussian weights and exact known
polynomial moments, interpolation/depth axes, no extrapolation and size refusal.
Gaussian agreement is reported numerically, never promoted into physical error
or convergence policy. FP32 alpha rounding remains part of the queried function.

Separate the observed BF refinement into the exact telescoping identity:
Tfine-Tcoarse=(Gfine-Gcoarse)+[(Tfine-Gfine)-(Tcoarse-Gcoarse)]. G16 integrates
each saved-J representation. Report roundoff closure relative to summed positive
coefficients under128eps, plus changes separately; maxima are not additive.
Keep whole-graph shadow matrices with BB untouched for every order/grid. These
are NOT a modified physical RT/SE operator or a selected best-fitting result.
Unknown between-node radiation accuracy, finite-CRD effects, absolute EOS and
coupled/observational admission remain open regardless of Gaussian agreement.

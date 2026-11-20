# Where stock Lightweaver differs from the codes our references were made with

RYA-1282, 2026-11-20. The Bride runs on Lightweaver 0.16.1. Our oxygen atom (atom.o41f)
was written for MULTI 2.3, the reference departures come from the Gerber, Bergemann et
al. 2023 MULTI grid (A&A 669, A43), and the atmospheres and background opacities come
from MARCS / Turbospectrum (TS). Four places were found where Lightweaver, used as
shipped, does something different from those codes, plus one solver weakness (section 5) and a set of missing ultraviolet sources (section 6).
Each one is now fixed in the Bride and pinned by tests. Anyone using Lightweaver elsewhere in the Codex should read this.

Benchmark: the solar node of the Gerber+2023 grid (MARCS 5750/4.5, A(O) = 8.70),
graded by `scripts/bride_qa_benchmark.py`.

## 1. Collision directions (fixed, c364add1)

The MULTI atom format says what each collision row means. Lightweaver reads two row types
differently.

* **CH** (hydrogen collisions). MULTI: the coefficient times n(H I) is the *downward*
  rate. Lightweaver reads it as the *upward* rate.
* **CH0** (charge transfer with H). MULTI: the rate runs in the direction the row is
  written (here O I to O II). Lightweaver sorts the two levels and makes it a
  downward-only rate, with no reverse.

The CH0 misreading drained O II and broke the LTE limit at depth (|b - 1| = 1). The CH
misreading cannot be seen by the LTE limit, because the misread rate is still balanced;
it shows up against the Gerber departures instead (triplet error 0.088 dex, against
0.0047 dex once fixed). The fix is the classes `MultiCHDownward` and `MultiCH0Directed`
in `pipeline/multi_atom.py`. Tests: `tests/test_bride_collision_semantics_rya1282.py`,
and negative controls in `tests/test_bride_qa_rya1282.py`.

## 2. How collision tables are interpolated in temperature (fixed, 54aafcfb)

MULTI interpolates every collision table with a taut spline (TAUTSP, tension 5.5). It
works on log(rate) for CE, CI and CH0 and on the rate itself for CH. Outside the table it
holds the end value. Lightweaver uses a different interpolator (WENO4) on the rate
itself.

On the oxygen atom the two disagree, per coefficient, by up to 12% (CE), 42% (CI) and
29% (CH0, median 11%). No standard interpolator gets within 2.5%. We ported MULTI's
routines to Python (`pipeline/multi_tautsp.py`). The port agrees with the compiled
MULTI Fortran to 1.4e-14 on all 1106 table rows. Tests: `tests/test_multi_tautsp_rya1282.py`.

## 3. The top boundary (fixed, 30fb2144)

Lightweaver's usual top boundary assumes no light comes in from above the model. A MARCS
model stops at a finite height, and the gas above it still shines in strong ultraviolet
lines. With an empty sky, the populations in the top layers collapse.

The Bride now places a slab of the overlying gas on top of each column. Its thickness
comes from the model's own column mass (or the pressure scale height for 3D columns),
and its light comes from the top layer's own opacity and emission. The solar worst-case
departure error fell from 2.98 to 0.36 dex. The 777 triplet itself did not change.
Code: `ColumnCap` in `pipeline/bride_top_cap.py`, solved with J and the cap converged
before each population step. Tests: `tests/test_bride_column_cap_rya1282.py`.

## 4. Ultraviolet background opacity (fixed, this change)

This is the largest difference. Lightweaver's passive atoms (from RH) give each metal
only a few bound-free continua. Its Fe I ground state is also one merged level (g = 25),
which makes its own Fe I partition function about 15% too large at 4000 K. MARCS and TS
use full Opacity Project tables for C, Mg, Al, Si, Ca and Fe.

We compared the Bride's background absorption with a native TS run on the same solar
model. TS reproduces MARCS to 0.12% in this region. In the top layer, Lightweaver gave:

| Wavelength (nm) | 205 | 215 | 240 | 260 | 277 | 300 | 355 |
|---|---|---|---|---|---|---|---|
| Lightweaver / TS | 0.64 | 0.16 | 0.42 | 0.07 | 0.27 | 0.61 | 0.55 |
| Bride now / TS | 0.98 | 0.96 | 0.98 | 0.96 | 0.96 | 0.96 | 0.86 |

The missing absorption was mainly Fe I, Al I and Mg I. It matters for oxygen because
light at these wavelengths ionises O I: the 3s 5S level (the lower level of the 777
triplet) ionises below 277 nm, and the 3p 5P levels below 431 nm. Too little background
absorption lets too much of this light through, so too much O I is ionised.

What the Bride now does (`pipeline/bride_ts_background.py`, `pipeline/ts_jonabs.py`):
for each passive metal ion that TS tabulates, Lightweaver's bound-free opacity is removed
and TS's is put in its place. TS's table is read and interpolated exactly as TS does it,
the partition functions are TS's own (Irwin), and the number of atoms comes from the
Bride's own equation of state. Ca I, which Lightweaver's calcium model lacks, is added
from Ca II by the Saha equation.

Checks against native TS:

* Fe I opacity per atom matches to 1.3e-4, and the Fe I partition function to 2.6e-6.
* Total absorption over 200-360 nm at all 56 depths: typical difference 0.226 dex before,
  0.025 dex after; worst case 1.28 dex before, 0.10 dex after.

Effect on the solar benchmark (top cap and MULTI interpolation in both cases):

| Grade | Lightweaver background | TS background | Reference |
|---|---|---|---|
| Triplet levels, photosphere | 0.0045 dex | 0.0022 dex | |
| All levels, photosphere | 0.224 dex | 0.057 dex | |
| All levels, all depths | 0.356 dex | 0.122 dex | |
| O II departure, top layer | 18.5 | 10.8 | 8.17 |
| Iterations to converge | 533 | 336 | |
| 777 NLTE EWs (mA) | 70.73 / 60.16 / 47.08 | 70.65 / 60.09 / 47.01 | |

The 777 lines barely move (about 0.1%); the atom's populations, especially O II, move a
lot. Tests: `tests/test_bride_ts_background_rya1282.py`, including a check that the old
background fails the same comparison.

## 5. Background scattering is iterated slowly (fixed, 554c6b73 + c6e94706)

This is a solver difference rather than a physics difference. Lightweaver puts background
scattering (Rayleigh and electron scattering) into the formal solution using J from the
previous solution, a plain "Lambda iteration". Each pass shrinks the error by roughly the
fraction of the opacity that is scattering. In the solar model that fraction is modest and
J converges in about 90 passes (0.78 per pass). In metal-poor models the upper layers are
almost pure scattering near 122 nm (99.8% at 5750/4.5/-2), so each pass removes only
0.5% of the error. J had not converged after 3000 passes.

Before the top cap, J was updated once per population step, so this slowness was hidden
inside the population solve. That is the likely reason three of the [Fe/H] -2 benchmark
models never converged.

The Bride now uses accelerated Lambda iteration (ALI) for the scattering term, the
standard method: each update divides the change by (1 - Lambda* x albedo), with Lambda*
the local part of the transfer operator, including the cap's own scattering. Ng
acceleration is applied on top. Passes needed to reach 1e-10 on 5750/4.5/-2:

| Method | Vacuum top | Top cap |
|---|---|---|
| Lightweaver's plain iteration | > 3000 | - |
| Ng only | 491 | - |
| ALI only | 268 | 432 |
| ALI + Ng (Bride) | 66 | 132 |

All methods reach the same J (to 2.4e-8 or better), so the answer is unchanged; only the
speed differs. (These counts predate section 6. With the full ultraviolet background the
extra absorption lowers the 122 nm albedo of 5750/4.5/-2 to 0.949 and ALI + Ng needs 49
passes; the most scattering-dominated model held is now 5750/3.5/-2, albedo 0.982, 63
passes, while plain iteration is still far off after the same number.) Code: `scattering_lambda_star` and `_fixed_population_radiation` in
`pipeline/bride_columns.py`. Tests: `tests/test_bride_ng_conservation_rya1282.py`.

Related: Lightweaver's own Ng acceleration of the populations can occasionally break
conservation of the atom's total number (2.7e-10 on one 3D column, against 1.8e-15
without it). The Bride rejects such a step and redoes it without acceleration
(f3bdc4f9).

## 6. Missing ultraviolet sources in metal-poor stars (fixed, ffbcd4f1 + 789fd36c)

Section 4 matched TS at solar metallicity. A native TS run for all 14 benchmark models
(150-450 nm) showed that in metal-poor photospheres most of the far-UV absorption comes
from sources Lightweaver does not have at all. At [Fe/H] = -2, Lightweaver's background
was only 8-10% of TS's at 150-200 nm and 63-69% at 200-277 nm.

| Source | Share of TS absorption, [Fe/H] -2 photosphere | In Lightweaver? |
|---|---|---|
| Dissolved Lyman continuum (hbop, n=1 bound-free extrapolated past 91 nm) | up to 90% at 150-200 nm | no |
| H + H collision-induced absorption (Doyle 1968) | 39-66% at 175-230 nm | no |
| OH continuum | 20-55% at 150-230 nm (cool stars) | no |
| H2+ | 4-16% | no |
| CH continuum | 0-3% (12% near 355 nm at solar) | no |

The Bride now adds all five, read and interpolated exactly as TS does
(`pipeline/ts_jonabs.py`, `pipeline/ts_hbop.py`, `pipeline/bride_ts_background.py`).
CH and OH amounts come from Lightweaver's molecular equilibrium, run for the
background only; the oxygen atom keeps its own populations. The hbop port agrees with
TS's compiled Fortran to 3e-6 (median).

Agreement with TS (photosphere, oxygen excluded, typical difference in dex):

| Band | Lightweaver background | Bride now |
|---|---|---|
| 150-200 nm | 0.14-1.21 | 0.007-0.106 |
| 200-277 nm | 0.15-0.48 | 0.005-0.021 |
| worst point, 200-431 nm | 0.61-1.50 | 0.017-0.17 |

Effect on the benchmark (all 14 models converge in 6-31 steps): the LTE-limit error at
depth in metal-poor models fell from up to 2.3e-4 to 2e-6 or less, and the O II top
departure moved toward the reference in every model. The 777 triplet error at
[Fe/H] -1 and -2 did not change (see Still open).

Two related bugs were fixed along the way: Lightweaver's `compute_rays` deep-copies the
Context, which the molecular background could not survive (27b35842); and TS prints
`T=10191.4` without a space in hot layers, which broke our reader for 6500 K models
(118f5eb5).

## Still open

* **Triplet away from solar metallicity.** Metallicity ladder (2026-11-20, c6abf9ca,
  `bride_mac_checks/rya1282-qa-zladder-v1`): 120 exact Gerber nodes, Teff 4500/5000/5750/6500,
  log g 3.5/4.5, all 15 [Fe/H] from +1 to -5, A(O) = 8.70 + [Fe/H] + [O/Fe]. All 120
  converge (5-33 steps); LTE limit <= 6e-6, conservation 2e-15, detailed balance 5e-16.
  Triplet photosphere max |dlog b| (5750/4.5): +1 0.088, +0.5 0.040, 0 0.002, -0.5 0.018,
  -1 0.026, -3 0.023, -5 0.022. Worst rows: 6500/3.5 (0.03-0.12 everywhere), 4500/4.5 at
  -3..-5 (0.07-0.09). The error is smallest at solar and grows on BOTH sides: at +1 the
  Bride's 3p 5P levels are further below LTE than MULTI's, at -1..-5 closer to LTE.
  Ruled out by test: MULTI's O count ignoring the A(O) node (the record responds to the
  node); a shifted depth axis (record tau = MARCS tau500 exactly); the background metal mix
  (Bride with solar metals moves the error 0.088 -> 0.083 at +1, 0.0225 -> 0.0216 at -3);
  microturbulence (grid Vturb 1 km/s = MARCS t01). Earlier note: against the Gerber+2023 MULTI grid the triplet error is
  0.001-0.005 dex at [Fe/H] 0 but 0.009-0.046 dex at -1 and -2. The Bride's 3p 5P
  levels sit closer to LTE than MULTI's, and its three J sublevels are nearly equal where
  MULTI's spread. Scaling the hydrogen collision (CH) rates would fit it, but the best
  factor depends on metallicity (1 at 0, about 0.5 at -1, about 0.3 at -2), so it is not a
  units or direction error. MULTI 2.3 applies CH exactly as the Bride does (down = n(H I)
  x coefficient, up by detailed balance; mul23_opacu.f gencol). Not a ray-count effect
  (3 vs 6 rays: 0.026 vs 0.028). Next suspects: MULTI's own background opacity
  (Uppsala package, not TS) and the hydrogen populations it used.
* **Atom counts.** The Bride's own equation of state gives 3-4% fewer neutral Fe atoms
  than TS in the upper layers. We keep the Bride's own count rather than borrow TS's.
* **The remaining O II gap** (10.8 against 8.17 in the top layer) is not yet explained.

## For other Codex work

Anything that uses Lightweaver's built-in passive atoms for an ultraviolet background,
or reads a MULTI atom with Lightweaver's own converter, has these differences. Use the
Bride's converter (`pipeline/multi_atom.py`) and background (`background='ts_metals'` in
`pipeline/bride_columns.column_context`), or check the result against these numbers.

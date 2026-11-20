# Where stock Lightweaver differs from the codes our references were made with

RYA-1282, 2026-11-20. The Bride runs on Lightweaver 0.16.1. Our oxygen atom (atom.o41f)
was written for MULTI 2.3, the reference departures come from the Gerber, Bergemann et
al. 2023 MULTI grid (A&A 669, A43), and the atmospheres and background opacities come
from MARCS / Turbospectrum (TS). Four places were found where Lightweaver, used as
shipped, does something different from those codes. Each one is now fixed in the Bride
and pinned by tests. Anyone using Lightweaver elsewhere in the Codex should read this.

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

## Still open

* **CH continuum.** TS includes ultraviolet absorption by the CH molecule (up to 12% of
  the top-layer total near 355 nm). The Bride does not yet, because it needs the amount
  of CH, which needs a small molecular equilibrium (carbon is largely locked in CO in
  the upper layers). This is the main remaining background difference.
* **Atom counts.** The Bride's own equation of state gives 3-4% fewer neutral Fe atoms
  than TS in the upper layers. We keep the Bride's own count rather than borrow TS's.
* **The remaining O II gap** (10.8 against 8.17 in the top layer) is not yet explained.

## For other Codex work

Anything that uses Lightweaver's built-in passive atoms for an ultraviolet background,
or reads a MULTI atom with Lightweaver's own converter, has these differences. Use the
Bride's converter (`pipeline/multi_atom.py`) and background (`background='ts_metals'` in
`pipeline/bride_columns.column_context`), or check the result against these numbers.

# VoronoiRT rung-1 evaluation — RYA-1012

**Status: rung 1 PASSED. Continuum now runs end-to-end on a cube we built.
The blocker is a cube with real horizontal structure, not the toolchain.**
Date 2026-08-29. Ran on Sirius (Ubuntu 26.04, x86_64), never on the Mac; all work
at `nice -n 19`, two CI jobs completed normally inside the window.

VoronoiRT is **external MIT code we are evaluating, not adopting**. Per RYA-1012 G3
nothing from it is committed here — this repo carries only our own tooling
(`scripts/rya1012_*.py`) and this record. The evaluation tree lives outside every
codex worktree at `/home/damienabraxas/scratch/rya1012_voronoirt/` on Sirius.

---

## The one-line answer

**Yes, VoronoiRT stands up on our hardware and reproduces its own published test
output — to 9 ULP.** But the searchlight is the only arm that could ever have
passed: it is the one arm that touches no opacity, and *the continuum and line
paths have been broken against every released Transparency.jl since October 2022.*

---

## Results

| test | verdict | numbers |
|---|---|---|
| `compare_searchlight.jl` | **PASS** | `I_160_45_regular` bitwise identical; `I_20_15_regular` 1750/2401 elements bitwise, max **9 ULP**, max rel **1.16e-15** |
| flux invariant | **PASS** | 80.00000000000027 / 80.00000000000016 vs exact 80; 79.999964714–80.000000000001 across all 360 scanned angles |
| `compare_continuum.jl` (regular) | **RUNS** | **42.29 kW nm⁻¹ m⁻²** at 500 nm disk centre — solar reference ≈40.8, **within 3.6%**; uniform across map to 6.4e-14 |
| `compare_continuum.jl` (irregular vs regular) | **RUNS** | max rel diff **0.50%** — *the paper's own criterion* |
| `compare_line.jl` | not run | same cube wall; needs **no** atom choice (see below) |

The paper itself judges the searchlight **qualitatively** ("visual comparison",
peak reduced to ~70%); no numerical tolerance exists in the paper or the repo.
Our 9-ULP match is a stricter test than the authors applied.

---

## Environment (reproducible)

Julia **1.10.10**, voro++ `b0dac57`, VoronoiRT `39c0c68` (2023-06-05, MIT © 2021
Udnaes), gcc 15.2.0, PyPlot 2.11.6 against system matplotlib 3.10.7.

🔴 **Julia 1.9.4 — the era-appropriate version — cannot run on Ubuntu 26.04.** Its
bundled `libopenlibm.so` is `GNU_STACK RWE`; the kernel refuses it
(`cannot enable executable stack`). 1.10.10 ships it `RW`. Do not re-diagnose this.

🔴 **Pin `Transparency` to rev `v0.1.8` for anything touching opacity** (drags
Interpolations back to 0.13.6). See the dependency finding below.

---

## Findings about the code

**F1 — the environment is unpinned.** No `Project.toml`, no `Manifest.toml`. The
13-package set had to be reverse-engineered from `using` statements.

**F2 — the continuum and line paths are broken against modern Transparency.jl.**
VoronoiRT calls eight opacity functions — `hminus_ff`, `hminus_bf`,
`hydrogenic_ff`, `hydrogenic_bf`, `h2plus_ff`, `h2plus_bf`, `rayleigh_h`,
`thomson`. **Every one was removed at Transparency v0.2.0 (2022-10-23) — eight
months BEFORE VoronoiRT's last commit (2023-06-05).** They survive only in
v0.1.8. Current `main` retains just `calc_hminus_density`, `gaunt_bf`, `gaunt_ff`,
`n_eff`; the opacities migrated into Muspel.jl (`background.jl`). So the published
repository could not have run its own continuum or line test against any
contemporaneous release — the author was on an unrecoverable local state. This is
exactly what F1 costs.

**F3 — the shipped scripts do not run their own tests.** All three test calls at
the bottom of `compare_searchlight.jl` are commented out; the live entry is
`do_timing()`, which needs the cube. The `npzwrite` calls that produced the
bundled references are commented out too.

**F4 — `read_quadrature` parses the angle count out of the FILENAME** (digits
after the first `n`), not the file. A sensibly-named file dies with an opaque
`ArgumentError: input string is empty`.

**F5 — no shipped quadrature contains the reference angles.** `ul7n12.dat` never
yields (20,15) or (160,45). `data/searchlight_data/plots/` is labelled `160_330`
against an array named `I_160_45`.

**F6 — `I_20_15_regular.npy` was generated at φ=195°, not 15°** (195 = 15+180).
Recovered by scanning all 360 integer φ at θ=20: φ=195 matches to 2.08e-17, φ=15
is off by 8.07e-01 — the full beam amplitude. Its sibling `I_160_45` *does* match
at its labelled angle. **A reference file's name is not its provenance.**

**F7 — the Voronoi arm is not reproducible, even against itself.**
`searchlight_irregular()` draws `rand(3, n_sites)` with **no `Random.seed!`**
(seeds exist in `compare_line.jl` = 2022 and `compare_continuum.jl` = 1998). Two
runs, same machine, minutes apart: **30.6% / 15.1%** flux spread — as large as the
28.7%/20.1% gap to the published arrays. The published-vs-ours mismatch is
therefore fully explained by RNG and is **not** a portability failure.

**F8 — only the searchlight has numeric references.** `data/` holds 8 `.npy`, all
searchlight. Continuum and line ship PDFs only. **But this does not block them**:
the continuum criterion is *internal* (irregular vs regular on the same cube).

**F9 — `plot_utils.jl` calls `pyplot()` at module load**, so PyPlot + a working
matplotlib is a hard dependency of the entire module, even headless.

**F10 — VoronoiRT has x and y transposed.** Muspel's `read_atmos_rh` reads the
identical file layout as `nz, ny, nx, nhydr, nt`; `get_atmos` labels those dims
`z, x, y`. Invisible on a square box (their 256×256 cube, the 51³ searchlight) —
**wrong for any non-square cube or x/y-asymmetric diagnostic.**

**F11 — `compare_line.jl` needs no science decision.** `test_atom()` is hardcoded
hydrogen Ly-α (χu = 82258.211 cm⁻¹, gl=2, gu=8, f=0.4162, Z=1) and `nλ_bb`/`nλ_bf`
are literals. Its only blocker is the cube.

---

## The cube

`bifrost_qs006023_s525.hdf5` is **not public**. The A&A paper (arXiv 2306.01041)
gives it as 256×256×430 over 6×6 Mm² — matching the naming convention (`qs006023`
= quiet sun, 6 Mm, ~23 km) — and carries **no data-availability statement**.

🔴 **`get_atmos` reads the RH 1.5D HDF5 format** — Pereira's own, and he
co-authored VoronoiRT. Identical dataset names and SI units, and the cryptic
`[:,:,:,1,1]` indexing lands exactly on RH's `(nt, nhydr, nx, ny, nz)`. **RH 1.5D
HDF5 is the lingua franca** across VoronoiRT, Muspel.jl, RH 1.5D and helita —
build a cube once, reuse it everywhere. Converters already exist:
`helita.sim.bifrost.BifrostData(...).write_rh15d()` (one call), or
`helita.sim.rh15d.make_xarray_atmos(...)` from plain arrays.

**What we built and validated:** `scripts/rya1012_build_rh15d_cube.py`, seeded
from Muspel's bundled `FALC.hdf5` (itself made by `make_xarray_atmos`). Two edits
are required and both cost a run to find — `velocity_x`/`velocity_y` must be added
(`get_atmos` reads them unconditionally; RH omits them), and `x`/`y` must be given
real spacing. ⚠️ **A zero-extent horizontal axis makes VoronoiRT exit 0 while
writing an entirely-NaN intensity map.** Read the values, never the exit code.

⚠️ **Honest limit on the 0.50% number:** our cube is FALC tiled 3×3, i.e.
horizontally **homogeneous**. That validates plumbing and 1D consistency, not 3D
structure. On a real granulation cube the spread will be larger, and only then is
it a physics test. The 48 km horizontal spacing is a free parameter we chose.

### Public cubes

🔴 **`http://sdc.uio.no/vol/simulations/` is a plain Apache index** — six Bifrost
runs, no registration (the JS search page exposes no links; `curl` needs `-L`).

| run | grid | box | spacing |
|---|---|---|---|
| `qs006023_s525` (theirs) | 256×256×430 | 6×6×8.7 Mm | 23 km — **not public** |
| `qs006005_dyc` | 1200×1200×1736 | 6×6×10.5 Mm | 5 km — same box, ~60 GB/snap |
| `qs024048_by3363` | 504×504×496 | 24×24×17 Mm | 48 km |
| `en024048_hion` | 504×504×496 | 24×24×17 Mm | 48 km — **non-equilibrium H** |

`en024048_hion/atmos/` ships per-variable-per-snapshot FITS at 481 MB each:
`lgtg` (log₁₀ T/K), `lgne`, **`lgn1..lgn6`** (6-level non-equilibrium hydrogen —
`lgn1` is the ground state, exactly the `[:,:,:,1,1]` slice), `ux/uy/uz`. Minimum
for VoronoiRT ≈ **2.9 GB**. Headers are self-describing: `DO_HION=1`,
`HATOMFIL=H_6.atom.ccpol`, `PERIOD_X/Y=1 PERIOD_Z=0` (matching VoronoiRT's
periodicity exactly), `TEFF=5773`.

⚠️ **Blocker on that route: z is non-uniform and the z array is not in the FITS.**
No extension (8 header blocks + padded data, verified against Content-Length); it
lives in `cb24bih.mesh`, which is not on the server, and `atmos_derived/` is empty.
`CDELT3` is only nominal. A wrong z corrupts optical depth.

---

## Muspel.jl — investigated and ruled out as a vehicle

Actively maintained (v0.2.6, last commit 2026-04-17, full compat bounds), unlike
VoronoiRT. **But it does not solve NLTE**: no statistical equilibrium, no lambda
iteration, no rate matrix, no ALI/MALI anywhere in `src/`. It has formal solvers
(`lte.jl`, `saha_boltzmann`, `feautrier`, `piecewise_1D_*`) plus `read_pops_rh` /
`read_pops_multi3d` / `ExtinctionItpNLTE` — it **consumes** NLTE populations
computed by RH or Multi3D. So it does not close the RYA-1008 gap; it independently
**confirms** it, from a new direction: Pereira's own maintained package stops at
the formal solution.

---

## Recommended next steps

1. **Ask ITA/Oslo for `cb24bih.mesh`** (the z-axis for `en024048_hion`). This is a
   small, well-defined ask — a mesh file, not a simulation — and it unblocks the
   only public cube that ships non-equilibrium hydrogen populations. Cheapest path
   to a cube with real granulation structure. The author's contact is in the
   VoronoiRT README; the same message can ask about `qs006023_s525` itself.
2. **In parallel, convert our own STAGGER cube** (`t5777g44m00`, 80×80×240, already
   on disk at `/srv/codex/grids-overflow/`). Needs an EOS lookup (rho, e → T, nₑ;
   STAGGER's EOS python is broken-as-shipped — put both `post_processed/EOS/` and
   `stagger/` on `sys.path`) then `rya1012_build_rh15d_cube.py`. Independent of any
   external party, and it makes the cube axis element-agnostic for later rungs.
3. **Re-run the continuum comparison on a structured cube.** The 0.50% figure is
   currently measured on a homogeneous atmosphere and is a plumbing result. Only a
   granulation cube turns it into the physics test the paper claims. Use the
   `skip=n` decimation already in `get_atmos` to control cost.
4. **Fix the x/y transposition (F10) before trusting any non-square cube.** Both
   public options are square (504×504), which hides it — but our STAGGER cube is
   80×80 (also square) and any subcube may not be. One-line fix in a local patch;
   do not report results from a non-square cube until it is done.
5. **Do not spend effort chasing a published numeric reference for continuum/line.**
   None exists (F8), and none is needed: the criterion is internal. Treat "get a
   cube" and "get a reference" as separate problems — only the first is real.
6. **Leave rung 2 gated.** It remains behind Solar-done+presented and a read of
   Bergemann & Hoppe (arXiv 2511.04254, 143 pp, Living Reviews in Computational
   Astrophysics) — confirmed to exist. Nothing above requires opening rung 2.

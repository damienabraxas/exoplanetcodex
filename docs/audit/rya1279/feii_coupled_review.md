# Coupled oxygen response to the FeII neighbor

2026-11-07, RYA1281/1278, scientific HOLD. Historical worktree
`/Users/ryanschmitt/codex/rya1279`, branch
`ryandamienschmitt/rya-1279-1280-bride-oxygen-gates`, stacked base4635c742.
Starting evidence d7b834bf; latest ORCH pulled separately at main de96ae51 and
all48native continuum comparisons already exact. No other worktree changes.

## Predeclared scope

`config/experiments/rya1281_feii_coupled.json`: canonical column only,
three modes natural/cowley/provided, identical saved O populations as initial
guess for every leg. Fresh O SE plus radiation, unchanged atmosphere/EOS/ne/H/He,
bound background and41FeI. No charge/temperature/hydrodynamic feedback. Original
600s/1500iterations/1e-6outer convergence; original5rays/qscale2/7flux rays and
1e-8inner radiation. No compute extension or tolerance relaxation. Serial Mac,
one thread; pre-run load2.52/2.63/2.44. No Sirius/CI or priority claim.

Use unchanged source-matched collision spline and author-MN receipts. Require
joint population/radiation/cap/physical-rate convergence, then evaluate individual
collision classes afresh and compare to native C; independently reconstruct
transition R plus class C net/gross flows and conservation without accelerated
Gamma. Residual must remain<1e-6, conservation<=1e-10. Fixed O transfer is not
accepted as a new SE result. Each new final O state gets its own standalone
continuum, rather than reusing an old continuum after changing populations.

Source/canonical/main target, native widths and pinned source receipts, atmosphere,
saved populations, inherited provider/cap/solver and binding identities checked.
Natural and both ionic Stark prescriptions retained without fit/preference.
Preserve old outputs. Populations/rates/J/history stay in immutable registered
Mac store; no binary rename or5MiBguard bypass. Compact reports/flux/history can
be copied with actual generators. Downstream continuum uses verified latest1315
implementation; both references/RVs retained when comparing spectra.

No new adopted abundance/opacity/default/product, production/merge/publication,
full3D or all-band claim. Broader source/model/EOS/reference and missing neighboring
opacity gates remain OPEN. Precommit plan/runner before observational evaluation.

Exact runner: held Bride interpreter after `source config/hosts/mac.env`, with
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`,
`python scripts/run_bride_feii_coupled.py --case canonical_refined --mode natural --output <fresh-directory>/natural`;
repeat serially for cowley and provided. Do not reuse final populations between modes.

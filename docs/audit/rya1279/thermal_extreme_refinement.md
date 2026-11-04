# RYA-1281: numerical refinement at thermal-axis endpoints

Predeclared 2026-11-04, before new RT. Scientific HOLD. Continue the accepted
six-case thermal experiment without changing its source plan or physical inputs.
Run thermal_T5572, thermal_T5822, thermal_g4338 and thermal_g4538, in that order,
serially on the Mac using the already committed run_bride_thermal_sweep.py
--refine-numerics. Reuse each exact producer certificate; original results stay
immutable. Rays3->5, qscale1->2, flux rays5->7 only; retain original1500iteration/
600second limits and all residual/conservation tolerances. No retries with relaxed
budgets or favorable case selection. Retain any refusal/failure.

Compare each refined run with its exact unrefined physical counterpart, including
complete source/binding/controls checks, saved independent R+C residual and
conservation, raw spectra, three O-aperture areas and maximum profile differences.
The solar nominal refinement is an existing control, not a new solve. This checks
combined sampling sensitivity at these endpoints, not individual error sources,
depth convergence, atmosphere interpolation or physical correctness. Numerical
differences will be reported directly; no new accuracy threshold invented after
seeing results. Continuum/observation choices and adopted solar parameters stay
unchanged. Full EOS, independent O reference/atom and continuum gates stay open.

New immutable Mac store root: rya1281-thermal-extremes-refined-v1. Full R/C/J
archives remain there if over the repository size ceiling, with verified receipts.
Base main f8c723a85f905425b1f7ae189f15a8763998ffc0; stacked base
4635c7421a7dee8671ffefe05d771c68da1cd24e; initial HEAD
d61e4ab9b8034d337b2bf2564688a4da4e46094e. No merge/production/publication.

## Conditional initialization improvement, declared before the first outcome

The legacy runner starts every thermal refinement from the old solar control,
including off-solar cases. If this bounded replay exhausts its budget, retain
that complete failed attempt and stop the legacy serial batch. The next bounded
attempt may instead initialize each refinement from its own already converged,
exact physical baseline. Validate case/source/certificate/local fields and
original controls, saved population/rate/history hashes, positive nuclei totals
and exact LTE state before using that guess. Re-solve with all original joint,
fresh-rate, conservation, iteration and time limits unchanged. This improves an
initial guess, not the atmosphere, atom or fitted parameters. Record the alternate
initialization and code lineage explicitly; do not call the failed replay passed.
Use a separate immutable root rya1281-thermal-extremes-matched-v1 and keep the
legacy attempt root intact. No increase of the600-second per-attempt budget.

## Independent light analysis while the serial RT batch runs

Use only the already archived thermal/xi profile products to compare finite
temperature5672->5822K, gravity4.338->4.538 and xi0.9->1.1km/s response directions.
Check exact shared nominal arrays, observation, broadening and source receipts.
Compare unit Euclidean directions over the union of the same three AIR O777
apertures, retaining both broadening/RV variants. Report signed cosines and
singular values of the unit-column profile matrix, alongside existing area
directions. This is unweighted geometry on a fixed grid: no noise model,
covariance, inversion, uncertainty estimate or observational identifiability
claim. No chosen favorable variant, new RT, continuum correction or fit. It
asks whether line shape carries distinctions lost by integrating three areas.

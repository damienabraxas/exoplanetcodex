# Crossed continuum controls

2026-11-07, RYA-1281 / RYA-1278, draft PR598, HOLD.
Worktree `/Users/ryanschmitt/codex/rya1279`, branch
`ryandamienschmitt/rya-1279-1280-bride-oxygen-gates`; stacked base4635c742,
starting evidence3f0fb5ab, main459f5614. No unrelated worktree changes.

## Predeclared design

`config/experiments/rya1281_crossed_continuum.json`: both FeII Stark alternatives,
both columns, all three O lines, synthesis and IAG references. For each set,
cross the before/after denominator with before/after support in three designs:
original selection with clipping; original selection without clipping; actual
native coefficient support without clipping. All four corners retained. Final
coefficient support means the pixels used by the polynomial, not the postfit
survivor count (known7774actual73/reported70 distinction remains intact).

Native-clipping diagonal endpoints and final-coefficient diagonal endpoints
must reproduce archived level, slope, indices and observed area **exactly**.
Off-diagonal values are labelled counterfactual diagnostics, never adopted
continua. Finite path averages, interactions and closure are arithmetic accounting
of correlated changes, not physical causal fractions or independent uncertainties.
No new RT, O populations, widths/gf/RV fits, moved windows, reference preference,
source/observed conditioning, product/default/adoption/production/merge/publication.
Historical3percent rule retained by hash; main15percent rule not exercised.
The matched continuum audit uses the existing plus-gravity reference; both RV
variants remain in the preceding frozen transfer analysis, not newly refitted.

## Initial validation

Mac held Bride interpreter, one thread:
`python -m pytest -q tests/test_bride_crossed_continuum_rya1281.py tests/test_bride_continuum_audit_rya1281.py`
10passed in0.14s, exit0. Tests cover finite paths/interaction/closure/refusal,
analytic denominator effects, exact fixed supports, immutable arrays and a
discriminating clipped-versus-fixed control. Initial synthetic outlier fixture
triggered the native minimum-anchor refusal rather than clipping as intended
(1failed/9passed); sparse outliers corrected the test setup, no native-rule change.
No scientific implementation defect or silent fallback discovered.

Precommit plan/kernel/tests/runner before evaluation. Exact next command:

```sh
source config/hosts/mac.env
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/codex/venv/bin/python scripts/audit_bride_crossed_continuum.py --output <fresh-registered-crossed-directory>
```

Results, artifact registrations, review and final guards remain owed at this checkpoint.

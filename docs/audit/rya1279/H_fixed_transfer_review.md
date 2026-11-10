# RYA-1281 fixed-fraction hydrogen transfer interventions

Scientific HOLD. This follows the [completed coupled pair](H_nlte_restart_review.md).
It measures a conditional transfer intervention, not a unique decomposition into
continuum and collision physics or a preferred H/EOS.

## Predeclared method

The first anchor is the completed v2 model-H oxygen fraction field. Map that SAME
41level/56depth fraction field to each of the two explicit H-linked O totals,
then hold it fixed through fresh full-line and standalone-continuum radiation
and emergent-ray solves. The reciprocal anchor uses the completed v2 native-H
fractions under BOTH totals. Plans were committed before their respective runs.
Native-H never receives a newly solved population from the other frozen case.

Both retain the MARCS5750/logg4.5/xi1 column, T/Pe/ne/geometry/composition,
source O8.70/41-level atom, bare FastBackground and passive atoms, empty molecular
table, experimental geometric cap and original5rays/q2/6796SE-grid layout.
There is no added Fe/CN provider, EOS/chemistry/charge update, statistical-equilibrium
step, or imported J/rate field. Collision conventions affect the inherited coupled
populations; collisions are not separately diagnosed or varied here.

Fresh fixed-population radiation must converge in J AND cap to1e-8 (max1000
iterations); original600second per-case limits remain. Emergent flux uses7rays
and4001AIR samples,7763–7783A/0.005A, one canonical AIR-to-vacuum conversion.
Saved O populations must equal the imposed field exactly after full/continuum
ray solves; fractions across H cases agree at roundoff and totals conserve to1e-10.
Positive/finite flux, source/backend/state/common-parent guards and no-CN lookup
are mandatory. No rate-balance screen is assigned to deliberately frozen
populations, and `SE_updated`/`SE_balance_claimed` are explicitly false.

The same-H anchor case must replay its own coupled spectrum within a prospective
1e-8absolute normalized-flux and1e-8relative continuum screen. This screens
active-versus-detailed-static transfer interoperability, not physical accuracy.
Comparison uses the same ±0.5A windows at AIR7771.944/7774.166/7775.388.
No solar observation, continuum fit, broadening/RV or telluric adjustment occurs.

## Reproduction

Load `config/hosts/mac.env`; held Bride virtualenv with OMP_NUM_THREADS=1 and
OPENBLAS_NUM_THREADS=1. Runs execute their two cases serially, on the Mac:

```
python scripts/run_bride_H_fixed_transfer.py --output <fresh-model-anchor>
python scripts/run_bride_H_fixed_transfer.py --plan config/experiments/rya1281_H_fixed_transfer_reciprocal.json --output <fresh-native-anchor>
```

Original model-anchor generator is commit807f82f8. Commitfe917467 adds explicit
plan selection/reciprocal anchor while preserving the original default. All source
and output hashes/plan controls are retained in input/report manifests. Raw logs:
`/private/tmp/rya1281-H-fixed-transfer-v1.log` and
`/private/tmp/rya1281-H-fixed-transfer-reciprocal-v1.log`.

## Model-H anchor results

Both cases exit0. Model-H10.646059874910861s, native-H10.100358916912228s;
process peak669761536/758775808bytes, one thread. Peaks are whole-process
high-water marks in a serial two-case process, not isolated per-context memory.
Full radiation converges63iterations in both cases, final J changes
9.31869e-9/9.31590e-9 and cap6.52582e-10/6.46976e-10. Standalone continua
converge35iterations in both, final J7.08436e-9/7.08235e-9 and
cap2.75565e-9/2.74071e-9. No SE step occurs.

Fraction mismatch max1.1102230246251565e-15relative; O total conservation
2.220446049250313e-16/8.881784197001252e-16. Same-H normalized replay
max4.9960036108132044e-15; continuum replay exactly equal.

|AIR centre(A), ±0.5A|Coupled area contrast(mA)|Fixed model-anchor contrast(mA)|Conditional remainder(mA)|
|---|---:|---:|---:|
|7771.944|0.19432663|0.13913615|0.05519048|
|7774.166|0.17707735|0.12970829|0.04736906|
|7775.388|0.14926984|0.11304469|0.03622515|

Peak absolute normalized-flux contrasts: coupled0.0008426290723997365,
fixed0.0005782634362950789, conditional remainder0.0002643656361046576.
Maximum relative fixed continuum contrast0.0023478576895832326, essentially
the same continuum response as the coupled pair under these assumptions.
The model-anchor plot was inspected: all three contrasts have the same core
sign, and fixed-fraction transfer retains most of each finite-window area change.
These signed-area ratios are not fractions of uniquely attributable physical
variance, collision contributions or abundance corrections.

## Reciprocal anchor and intervention dependence

Both reciprocal cases exit0. Model-H10.492319208104163s and native-H
10.14155716681853s; whole-process peaks675790848/762167296bytes, one thread.
Full radiation again63iterations in both, continua35iterations; final J/cap
values remain below1e-8. O conservation8.88178e-16/2.22045e-16, common
fractions agree9.99201e-16relative; populations remain exact after transfer.
The native-H anchor replay matches its coupled normalized flux within
2.9976021664879227e-15, with exactly equal continuum. No SE solve.

|AIR centre(A), ±0.5A|Fixed native-anchor contrast(mA)|Conditional remainder(mA)|
|---|---:|---:|
|7771.944|0.13908874|0.05523790|
|7774.166|0.12966997|0.04740739|
|7775.388|0.11301478|0.03625506|

Peak absolute fixed/remainder normalized-flux contrasts are
0.0005783692230932891/0.0002642598493064474. Maximum relative continuum
contrast0.0023478576895767933. Both anchor plots were inspected; qualitative
core signs and contrast ordering agree. The largest wavelength-wise difference
between the two fixed contrasts is1.9599682254956008e-7. The three fixed-area
contrasts differ by0.00004741108/0.00003832651/0.00002990307mA between anchors.
These small differences do not establish physical interaction precision below
the inherited coupled joint1e-6accuracy. They support the qualitative attribution's
stability for these two interventions, not general path independence.

All four radiation controls took41.3803seconds in total, excluding setup/tests/
plotting and earlier coupled attempts. Each original600second case cap remained
unchanged. Twenty-four original-name files byte-match external originals;
maximum400167bytes. No compression, split/rename or5MiBguard bypass.

## Interpretation and next task

The intervention includes changing the H-dependent background AND rescaling the
O total while holding its fractions fixed. The coupled-minus-frozen remainder
measures the additional response conditional on this intervention. It must not
be relabeled a pure collision effect or a reference-calibrated error budget.
The reciprocal anchor tests dependence on the chosen fixed fraction field;
it is not an independent atmosphere/EOS/source benchmark.

Next return both completed coupled spectra and fixed controls to the previously
pinned current-ORCH solar diagnostic, holding the AIR coordinate ownership,
continuum/telluric policy and RV/broadening choices fixed. First verify that the
bare-background diagnostic is compared on equal processing terms with its own
baseline; do not silently compare it as a replacement Fe/CN-blanketed product.
Archive residual changes with no observational tuning or abundance adoption.
Primary molecular/EOS/native-state provider qualification remains required before
inserting real CN; synthetic purpose guards remain intact. C/N then Fe and the
reference-star/adaptive3D/all-band work remain open.

## Verification and provenance limits

Nine focused H-fraction/cap/binding tests pass1.58s with the Mac profile. The
end-to-end experiment additionally checks both anchor replay screens, frozen
identity, radiation convergence and matched source/state/flux guards. No numerical
kernel or adopted abundance/default/product changes. Mac only, no Sirius/CI.
Main remains7f023fd7, the inspected PR615Ca notes; no merge/rebase.
Fresh library212documents/127accounted/85unindexed/12off-Mac,exit1; no new
paper used. Primary atom lineage, independent MULTI departure discrepancy,
collision/cap hypotheses, chemistry/EOS and reference-data gates remain open.
No model/abundance/production/publication/merge admission.

MIRROR/version-pointer/sequence/staged-whitespace guards pass. Generator invocation
audit remains exactly2436historical diagnostics, zero added/removed. This is not
a global-green claim; prior RYA587blank science_impact and tracker/gold debts
remain. No new source/numerical defect was identified for the bug ledger.
Postcommit stacked-base freshness and ordinary push hooks are required; no bypass.

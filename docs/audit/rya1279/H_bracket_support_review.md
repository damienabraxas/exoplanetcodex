# RYA-1281 fixed bracketing support and reference controls

Scientific HOLD. Fixed bracketing support removes the large synthesis-reference
extrapolation in the [preceding Kitt replay](H_kitt_orch_review.md), but leaves a
6–7mA conditional reference-denominator difference on exactly the SAME pixels.
This diagnoses support versus reference sensitivity; no continuum is adopted.

## Predeclared controls

The committed plan fingerprints the previous report, support traces, summary
and profiles/reference/continuum arrays. All six spectra, both RV controls,
three AIR±0.5A triplet windows, original broadening, observation and IAG holdings
are unchanged. No new spectrum/RT/SE/EOS, observation shift, telluric correction
or abundance fit occurs. AIR ownership remains as documented in the parent:
delivered observed vacuum nm converted once; model outputs already AIR.
Current ORCH checkout/source still equals pinned de96ae51, with actual model
level bound0.15, selector90th percentile, degree1 and native minimum-pixel count.
The new fixed-support fits are explicitly diagnostic, not the native estimator.

For each of6models x2RVs x3lines, cross2reference denominators with5supports:

- Native synthesis initial selected pixels, with clipping disabled.
- Native synthesis exact final coefficient pixels, with clipping disabled.
- Native IAG exact final coefficient pixels, with clipping disabled.
- Shared intersection of baseline synthesis/IAG global top-decile pixels.
- IAG top-decile pixels selected separately in the blue and red sidebands.

The last two masks are shared across all six models and are selected without
observed residuals. Sidebands are0.5<|AIR offset|<=5A. Shared masks must have
at least8pixels on EACH side and bracket the centre; insufficient geometry is
refused, never relaxed. Native support masks are retained as deliberate
extrapolation controls and are not mislabeled as geometrically validated.
Fixed coefficient supports never change through clipping. The same pinned0.15
level bound applies to all fits. Nonpositive continua are refused.

An independent linear-prediction reconstruction checks the centre coefficient.
Exact native coefficient-support/reference replays must match the previously
recorded native level, slope and application decision under a prospective
floating-point screen. Shared-mask IAG continua must be exactly equal across
model labels. Source arrays remain exactly unchanged after all calculations.

## Execution and outcome

With `config/hosts/mac.env`, project analysis environment, serial Mac threads:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/codex/venv/bin/python scripts/analyze_bride_H_bracket_support.py --output <fresh-bracketing-support>
```

Both attempts exit0. First v1contains the complete crossed fits; v2adds36finite
2x2accounting quartets comparing exact native synthesis support with balanced
support, and synthesis with IAG denominators. summary.csv, contrasts.csv,
fits.json and the plot are byte-identical across attempts. Original logs and
output directories are preserved; no overwrite or source/kernel fix required.
Accounting helper is fingerprinted; maximum arithmetic closure1.42109e-14mA.
Finite path averages/interactions are correlated intervention summaries, not
unique physical causes, independent uncertainties or abundance corrections.

There are360fits and360reference/H-residual contrasts.276fits apply and84are
refused as predeclared. All72exact native coefficient replays have ZERO level
and slope difference. Refusals are controlled experimental outcomes, not failed
engine runs or missing observations.

|Support|Applied|Refused|Reason|
|---|---:|---:|---|
|Native synthesis initial selection, no clipping|60|12|7775synthesis level~1.2399exceeds pinned0.15bound|
|Native synthesis exact coefficient support|72|0|Original coefficient fits replay exactly|
|Native IAG exact coefficient support|72|0|Original coefficient fits replay exactly|
|Shared global top-decile intersection|0|72|One-sided or insufficient overlap|
|Balanced IAG per-side top decile|72|0|45blue/45red pixels per window|

Plus-gravity baseline shared intersections contain36blue/0red pixels at7771.944,
0blue/4red at7774.166, and no pixels at7775.388. They do not yield usable
bracketing support. The prescribed fallback experiment is the independently
predeclared balanced mask, not a threshold lowered after seeing a refusal.

## Balanced-support results

Plus-gravity coupled model-H baseline (all finite areas in mA):

|AIR centre(A)|Synthesis level|IAG level|Synthesis-conditioned Kitt area|IAG-conditioned Kitt area|Model area|
|---|---:|---:|---:|---:|---:|
|7771.944|0.99696599|1.00375779|74.74933|81.00986|68.74211|
|7774.166|0.99699096|1.00404689|64.38518|70.96012|59.20498|
|7775.388|0.99671086|1.00401443|50.67184|57.57759|46.28401|

The synthesis/IAG conditioned-observation area differences on these identical
90pixel supports are6.26053/6.57495/6.90575mA. Remaining model-minus-observed
area residuals are-6.00722/-5.18020/-4.38783mA on synthesis denominators and
-12.26775/-11.75515/-11.29358mA on IAG denominators. The plotted cores remain
deeper than the model under both diagnostic placements. No smaller residual is
used to choose a reference, new continuum or abundance.

Earlier native synthesis levels0.99677/1.03447/1.13185become approximately
0.997throughout on balanced support. At7774/7775, switching support while
holding the synthesis denominator fixed changes conditioned observed area by
-33.90064/-113.30479mA. Switching denominator AFTER bracketing still changes it
by+6.57495/+6.90575mA. The large native apparent discrepancy is therefore
sensitive to support/extrapolation, while denominator sensitivity persists once
the line is bracketed. This does not adjudicate whether the residual denominator
difference is missing opacity, atlas normalization, line formation or another
source convention; those require separate source/state controls.

Coupled-H residual contrasts remain small. On balanced IAG support they are
0.194282/0.177044/0.149270mA; on balanced synthesis support
0.193812/0.176633/0.148687mA. Neither explains the larger residuals.
Both original RVs and all six source controls remain in the machine tables;
no optimized RV/broadening, all-band validation or noise/significance claim.

The plot was inspected: bracketing removes the large extrapolated synthesis
levels, and same-pixel observed/reference ratios retain an offset between the
two denominators. The balanced mask is selected from IAG and therefore remains
a reference-conditioned intervention, not an independent solar calibration.

## Artifacts, verification and next task

Six original-name outputs byte-match external v2; largest662316bytes. All masks,
indices, fits, refusal reasons and finite accounting are retained; parent arrays
are referenced by hash rather than duplicated. Actual writer/invocation and both
attempt resource roots are registered. No5MiBsplit/rename/guard bypass.
Twelve focused fixed-support/crossed/trace tests pass0.20s; actual end-to-end
replay additionally checks72coefficient identities, geometry refusals, fixed
support, common IAG arrays, source immutability and36accounting closures.
No new numerical/source defect was identified for the bug ledger.

Main remains7f023fd7, already inspected PR615Ca notes; no merge/rebase.
Fresh library212documents/127accounted/85unindexed/12off-Mac,exit1; no new
paper used. Mac only, no Sirius/CI/install/download; original source and adopted
abundance/default/product/continuum remain unchanged. Earlier ledger/tracker/gold
debts and all scientific source/atom/MULTI/collision/cap/EOS gates remain open.

Next bounded task: audit the SAME balanced pixels' synthesis, Kitt and IAG values
and source conventions directly. Quantify how reference flux-level/line-blanketing
differences map into the6–7mA denominator contrast before using a continuum to
interpret a solar abundance. Retain both references and current geometry; do not
force atlas ratios to unity or tune support thresholds to reduce residuals.
Verify flux/intensity, frame/resolution and normalization lineage from held
sources before attributing the difference. Primary/EOS/native-state CN provider
qualification remains a separate gate before real missing-opacity insertion;
no synthetic purpose guard bypass. Additional oxygen indicators, C/N then Fe,
reference-star and adaptive3D/all-band goals remain open. Scientific HOLD;
no production/publication/merge or model/abundance admission.

Final MIRROR/version-pointer/sequence/staged-whitespace checks pass. Generator
invocation audit retains exactly2436historical diagnostics, zero added/removed.
Previously verified global ledger RYA587blank science_impact and tracker/gold
debts remain; this is not a global-green claim. Postcommit stacked-base freshness
and normal push hooks are required, with no bypass.

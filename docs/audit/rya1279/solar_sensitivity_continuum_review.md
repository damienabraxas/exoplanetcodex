# Solar oxygen sensitivity and current ORCH continuum replay

Scientific HOLD. The controls below separate numerical convergence, physical input sensitivity and observation conditioning. They do not choose an abundance, best-fitting velocity, atmosphere or continuum reference. All work is local to the Mac, serial and one thread; no Sirius/CI or product/publication changes.

## Exact continuum policy comparison

The actual `pipeline/local_continuum.py` blobs at de96ae51fd8c52f79e5457fc0205e4baaa4596fd (RYA-1315) and 078cf7dedfa961452c8bd9af43c4bf83c5463fdf (merged RYA-1350) are loaded independently and copied into the evidence package. The new rule requires neighbourhood corroboration for a model-guided level more than 3% from unity, retaining the existing 15% outer limit. The executable callback requires at least three complete neighbour windows and a centre level within 0.02 of their median. This is an actual estimator replay, not a full ORCH synthesis or abundance run.

The nominal matrix has 96 rows: four qualified columns x two existing model RV controls x three O777 lines x Kitt/native IAG observations x synthesis/IAG continuum references. Continuum windows remain centred on the rest AIR line, matching the earlier ORCH replay. A separate 96-row diagnostic shifts only the window origin by the already declared gravitational model RV; observations are never shifted. It preserves the same refusal count and conclusion. The source Kurucz2005 corrected holding is vacuum-to-AIR converted once by its canonical reader; model outputs and IAG are already AIR. Actual holding gates and pre-normalisation guards run for both `solar_kpno_kurucz2005_corrected` and `solar_iag`. The raw sources remain unchanged. No new telluric operation is performed.

The nominal new policy changes 24/96 rows, refusing precisely the beyond-3% synthesis-reference corrections: Kitt 7774/7775 and native IAG 7775, across four columns and two RVs. All remaining 72 rows are byte-identical between policy versions, including fitted records and corrected arrays. None of the 48 IAG-reference legs changes. Nominal fitted levels span 0.990462-1.131856; largest old/new corrected-pixel difference is 0.225303. This is not a validated uncertainty or an abundance change.

All 48 synthesis-guided fits have one-sided coefficient support. For the solar/plus-gravity Kitt case, the centre lies 4.102,4.103,4.101 A outside the anchor hull; centre-prediction absolute weight sums are14.24,15.18,16.69. The fitted levels are0.996768,1.034470,1.131853. The last two are now refused. The IAG-reference levels for the same holding are1.004551,1.003892,1.004148; the first IAG-reference fit is also one-sided (1.352 A extrapolation), so near-unity level is not proof of good geometry.

The small-source synthesis spans7763-7783 AIR A, insufficient for the required independent +/-50 A neighbourhood. Missing wide synthesis is explicitly recorded (zero complete neighbour windows), never padded, extrapolated or replaced by an atlas under a synthesis label. The IAG reference has all ten neighbour windows. A beyond-3% synthesis correction is therefore refused here; obtaining a physically complete wider synthesis would be a new experiment, not permission to bypass the gate.

The safeguard reduces extreme corrections but does not solve continuum placement. In native IAG/solar/plus-gravity at7774, the synthesis-reference level1.024867 is within3% and remains applied; the model-minus-observation finite +/-0.5 A area changes to -31.570 mA. An accepted bound does not certify the one-sided linear extrapolation. These finite areas differ from the earlier fixed-shoulder metrics and must not be substituted into their factorial table. The bare O/background spectrum lacks the complete neighbouring line opacity; no preferred continuum/reference is adopted to improve agreement.

The independent trace records coefficient-support count separately from the estimator's reported n_bins and last-pass clipping status. This preserves the previously characterized final-clipping provenance issue; no postclip refit is silently applied.

Synthetic controls on the actual pinned modules show a smooth10% sag is retained when corroborated, a local10% dip is refused, missing neighbours refuse, and sub3% fits never request a wide reference and remain exact. The plotted nominal comparison preserves shipped, previous/current synthesis-reference and previous/current IAG-reference curves separately.

## Prospective solar parameter matrix

Plan `config/experiments/rya1281_solar_sensitivity.json` was committed at3cb22111 before the RT matrix completed. It declares a fresh canonical-solar replay,5762/5782 K at logg4.438/xi1, and xi0.25/0.5/0.75/1.25/1.5/2.0 km/s at5772/logg4.438. The baseline context is dry-run only. Ten dry contexts precede nine serial real solves. Original budgets, thresholds, atom, background, cap and physical-rate admission remain unchanged; each starts independently from the original common fractions, not a preceding case's endpoint.

Xi changes only the depth-constant RT Doppler velocity, retaining the same xi1 MARCS temperature/pressure/density/height structure. This is a controlled transfer experiment, not a turbulent-pressure atmosphere grid. The source helper tests verify exact field isolation and reject invalid negative velocity. Temperature uses the existing audited source-cell interpolation;10 K sampling does not establish10 K physical accuracy.

## Reproduction and evidence

Use `config/hosts/mac.env`; analysis interpreter is the registered project Python (NumPy2.2.6), while actual RT uses the held LW0.16.1 environment (NumPy2.5.2). No environment dependencies were installed or changed. Run the continuum script with either declared continuum plan and a fresh output directory. For the parameter matrix, run `scripts/run_bride_solar_matched_mali.py --plan config/experiments/rya1281_solar_sensitivity.json --case <declared-name> --output <fresh-directory>`, first with `--dry-run`. The sensitivity analyzer requires every declared real case and independently replays physical rates, conservation, common starts and solar endpoint parity before producing comparisons.

Earlier continuum preflights stopped on a missing provenance field and an incorrect IAG holding identifier; both occurred before corrected-spectrum output. The correct existing receipt and canonical holding identity are used in the qualified replay. Intermediate shifted-origin and reporting-label attempts remain external with exact receipts; final nominal and shifted-control outputs are archived separately. No failed result is promoted.

[Continuum archive](../../../data/results/rya1281/orch_continuum_parity_v1/receipt.json). Full external source outputs remain under the registered Mac experiment root. The register and numerical-risk ledger record the one-sided continuum finding; scientific admission stays HOLD.

## Completed sensitivity results

All nine independently initialized real cases pass the original physical-rate, conservation and resource guards. Total measured setup/solve/synthesis time is460.631 seconds; each uses8 or9 outer iterations and42.97-57.48 seconds. Maximum independently replayed physical-rate residual is1.950e-8; conservation is at most1.554e-15. The solar endpoint spectrum is byte-identical to the previous four-case solar result. The analyzer emits108 fixed-observable metric rows,324 observation comparisons and12 temperature/velocity finite-difference rows, preserving both RVs and both shoulder schemes.

For the predeclared plus-gravity/inner-shoulder convention:

| AIR line (A) | Lower10K strengthening to solar (mA) | Upper10K strengthening from solar (mA) | Central xi response (mA per km/s) |
| --- | ---: | ---: | ---: |
|7771.944|0.711501|0.714120|4.400779|
|7774.166|0.643041|0.645863|3.672294|
|7775.388|0.539287|0.542164|2.687009|

Temperature second differences are0.00262/0.00282/0.00288 mA over the20 K bracket: smooth, slightly curved conditional response. This does not establish interpolation or atmosphere accuracy.

| xi (km/s) |7771 area change (mA)|7774 area change (mA)|7775 area change (mA)|
| ---: | ---: | ---: | ---: |
|0.25|-2.177735|-1.820403|-1.335810|
|0.50|-1.724507|-1.441071|-1.056887|
|0.75|-0.989479|-0.826412|-0.605562|
|1.00|0|0|0|
|1.25|1.210911|1.009735|0.737943|
|1.50|2.607932|2.172449|1.585019|
|2.00|5.822882|4.839240|3.517097|

The velocity response is curved, so the local derivative is not extrapolated across the full sweep. Values are changes from the fixed xi1 solar anchor, not corrections chosen to close the atlas residual. No best velocity or inferred abundance is selected.

The unweighted three-line finite response vectors for T and xi have direction cosine0.996630 for inner shoulders and0.997126 for outer shoulders (almost unchanged between RVs). This measures near-parallel observable responses only, not a statistical covariance, confidence interval or proof that full line profiles contain no distinguishing information. It explains why the three feature areas alone cannot safely identify an incorrect input factor. Continuum conditioning remains a separate, much larger possible systematic in these bare-opacity diagnostic spectra.

[Exact sensitivity archive](../../../data/results/rya1281/solar_sensitivity_v1/receipt.json), [analysis report](../../../data/results/rya1281/solar_sensitivity_v1/rya1281-solar-sensitivity-analysis-v1/report.json), [velocity response plot](../../../data/results/rya1281/solar_sensitivity_v1/rya1281-solar-sensitivity-analysis-v1/microturbulence_response.png). Full J/rate arrays remain external and SHA-receipted.

## Next milestone

Before an inversion, compare the parameter responses with an oxygen diagnostic of different excitation/strength and test complete neighbouring line opacity and bracketed continuum support on the canonical solar column. Preserve the actual current ORCH rule as the nominal policy and any geometry/placement alternatives as explicit diagnostics. Run profile-width/broadening and abundance/gravity interventions separately, followed by joint checks that measure interactions; do not infer a preferred parameter from the triplet-area residual. Extend oxygen across additional available bands, then C and N before Fe and the reference star. True3D NLTE and adaptive cross-grid coverage remain foundational work, not validated by this column sweep.

Validation:119 focused numerical/control tests pass in6.63s (666 existing Lightweaver warnings);29 register tests plus the4 new/extended controls pass in the final33-test run. Target mirror and register structure pass. All118 new compact evidence files match their SHA/size receipts; four copied estimator blobs are byte-identical to pinned Git. Normal whitespace checks report two existing blank-at-EOF notices in the exact archived old source copies; checks pass outside those two immutable copies, which were not trimmed. Path guard has139 present and15 pre-existing broken paths. Generator guard retains2436 prior diagnostic strings, zero added/removed. Library audit remains212 documents/127 accounted/85 unindexed/12 off-Mac; refinement debt reports21 rows with10 unfired tickets. Existing global ledger/tracker/gold debts remain. No global-green claim.

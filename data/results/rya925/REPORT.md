# RYA-925 — Al solar Kitt Peak all-model sweep

Base: `origin/main` at `909b6395937440fc2ac14d3d2b2b5a4cc39afc02`. Compute: Sirius, dedicated worktree `/mnt/codex-data/codex/rya925`. Holding: `solar_kpno` (2960–13000 Å); no acquisition performed. Frozen gold was not read.

## Outcome

This is the **NSO/Kitt Peak FTS disk-integrated solar flux atlas** (`solar_kpno`), not HARPS. Its held coverage is 2960–13000 Å. The instrument and holding are now explicit columns in both matrices.

For the calibration question — can our engines recover the accepted abundance within combined measurement + literature uncertainty? — **all six VIS/red-optical products replicate the literature within 1σ**. They are not failures. The old `fail-with-reason` strings in `validation.json` answer a narrower central-value rule (“is the point estimate inside 6.39–6.47?”); that rule is retained as diagnostic history but is not the RYA-925 calibration verdict.

| Kitt Peak band | independent product | Codex A(Al) | literature | combined σ | abs(Δ)/σ | calibration |
|---|---|---:|---:|---:|---:|---|
| VIS | EW · 1D-LTE | 6.531 | 6.430 ± 0.040 | 0.212 | 0.48 | replicated |
| VIS | EW · 1D-NLTE · Amarsi | 6.351 | 6.430 ± 0.040 | 0.176 | 0.45 | replicated |
| VIS | Synth · 1D-LTE | 6.271 | 6.430 ± 0.040 | 0.176 | 0.90 | replicated |
| red-optical | EW · 1D-LTE | 6.470 | 6.430 ± 0.040 | 0.185 | 0.22 | replicated |
| red-optical | EW · 1D-NLTE · Amarsi | 6.380 | 6.430 ± 0.040 | 0.182 | 0.28 | replicated |
| red-optical | Synth · 1D-LTE | 6.387 | 6.430 ± 0.040 | 0.179 | 0.24 | replicated |
| near-UV | Synth · 1D-LTE (1/3 lines) | 4.198 | 6.430 ± 0.040 | 0.201 | 11.08 | not replicated |

The accepted 6.43 ± 0.04 is **3D+NLTE**, not 1D-NLTE. Current products are therefore diagnostic cross-scale comparisons; a scale-matched 3D-NLTE product remains owed. The near-UV one-line result is a genuine non-replication and cannot be used for calibration.

### How the literature value was obtained

Asplund et al. (2021) retain the Scott et al. (2015) result: screened Al I lines, a 3D hydrodynamic solar atmosphere, NLTE line formation, vetted transition probabilities/HFS, blend rejection, and three 1D atmospheres plus mean-3D comparisons to price systematic uncertainty. The result is 6.43 ± 0.04 (±0.01 statistical, ±0.04 systematic).

Andrievsky et al. (2020) used MULTI departure coefficients inside SYNTHV blended-spectrum synthesis across optical and IR lines. Crucially, they **adopted** solar Al=6.43 and adjusted IR atomic data to reproduce the solar spectrum. Their optical/IR agreement is therefore valuable line/calibration evidence, but not an independent blind recovery of 6.43.

The baseline does not reproduce the ticket's “6.443 ± 0.068 / 6 lines” claim. Repository history shows 6.443 was the first four-IR-line integrated-absorption result (`b09cff2`), later superseded by a corrected four-line profile result 6.415 (`09f1619`). The present six-line red-optical EW baseline is 6.470 ± 0.052 (stat) ± 0.173 (syst). Its changed line count comes from the current whole-line/profile-width gates; several recovered lines pin to the width floor, and the dominant systematic remains UNGRADED gf. This is traced divergence, not tuning.

## Model/route disposition

- Successful: VIS and red-optical `EW · 1D-LTE`, `EW · 1D-NLTE · Amarsi`, and `Synth · 1D-LTE`; near-UV `Synth · 1D-LTE` with one constrained line of three.
- Near-UV: A(Al)=4.198 from one constrained 3057 Å line; 3050/3066 are non-minima. It is a loud frontier discrepancy, not a promotion candidate. The 3961 resonance line lies in the VIS band and is skipped because it requires dedicated blend-aware optical synthesis; 3944 remains CH-blend “do not attempt.”
- NIR: no product. 10872 is present but its theoretical depth 0.061 is below the predeclared 0.15 selection floor; 11254 is inside the registered H2O interval.
- Gerber 1D-NLTE: attempted/unsupported. The committed registry says reachable-not-extracted. Sirius holds an Al deck, but production supports Fe only and assumes no abundance axis; Al's deck has an abundance axis. Emitting a Gerber number would be mislabeled, so the path remains loud-unsupported.
- Nordlander–Lind: the primary paper supports 6696/6698/7835/10872 and reports the solar independent `<3D>` NLTE aggregate 6.461 ± 0.022, but the actual mean-3D EW grid is not held or wired in the current repo. It is therefore an independent literature cross-check, honestly named `MEAN3D_NLTE`, not a rerun product and never `FULL_3D_NLTE`.

## Implementation findings

The non-Fe correction helper returned an empty map instead of consulting the registered Al grid, and the frontier synthesis fitter selected Al lines while varying Fe abundance. Both are fixed with regression tests. Physics-axis emission now identifies Al as model `amarsi`, not the legacy Bergemann family default. The minimum Al-only RYA-813 litscan projection is `data/reference/litscan/Al.yaml`; `validation.json` is produced through `validate_element('Al')`, whose structural no-gold guard runs first.

## Evidence

- `line_disposition.csv`: full measured/quarantined/skipped accounting, including 3944, 3961, 10872, and telluric 11254.
- `per_line_matrix.csv`: every emitted per-line treatment result and failure reason.
- `product_matrix.csv`: stored axes, display identity, bars, same-line LTE deltas, and domain status.
- `calibration_matrix.csv`: Kitt Peak provenance, combined-uncertainty replication statistic, and calibration verdict.
- `validation.json`: legacy RYA-813 central-value classifier by product/band; not the calibration verdict.

The line-specific gf conflict remains visible: RYA-835 supports the 7835/7836 pair, while 8772/8773 did not confirm into the same provenance class. No gf, continuum, or model parameter was adjusted toward 6.43. HARPS was not touched.

## Continuation owed: all regions, data, and models

- Repeat the matrix on the independent telluric-corrected IAG solar atlas (4047–10650 Å) wherever it covers the Al lines; do not mix it with Kitt Peak.
- Run Al on the landed HARPS solar loader (`solar_harps`, RYA-911 / PR #313) across its full 3780–6910 Å VIS arm. The loader is ready; no HARPS Al value exists yet.
- Build a telluric-safe IR solar arm for 10872/11254 Å. Kitt Peak retains tellurics; the held CRIRES+ Y product ends at 10680 Å and does not reach 10872 Å. An IR abundance cannot be manufactured from either gap.
- Wire the Al abundance axis in the Gerber 1D-NLTE adapter and acquire/register the Nordlander–Lind mean-3D grid before emitting either model as a Codex product.
- Run every supported route as a separate band × holding × model product. HARPS becomes the future stellar calibration arm; this solar Kitt Peak run is the reference-atlas calibration, not a HARPS proxy.

# RYA-925 — Al solar Kitt Peak all-model sweep

Base: `origin/main` at `909b6395937440fc2ac14d3d2b2b5a4cc39afc02`. Compute: Sirius, dedicated worktree `/mnt/codex-data/codex/rya925`. Holding: `solar_kpno` (2960–13000 Å); no acquisition performed. Frozen gold was not read.

## Outcome

The independent products are in `product_matrix.csv`; no products are averaged together. VIS validation uses the RYA-716 machine-readable 1D-NLTE comparator 6.43 [6.39, 6.47]. All three VIS products are `fail-with-reason`: EW LTE 6.531 is also scale-mismatched, Amarsi 6.351 misses the lower bound by 0.039 dex, and synthesis LTE 6.271 is scale-mismatched. Red-optical and near-UV are report-and-note only.

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
- `validation.json`: RYA-813 verdict by product/band.

The line-specific gf conflict remains visible: RYA-835 supports the 7835/7836 pair, while 8772/8773 did not confirm into the same provenance class. No gf, continuum, or model parameter was adjusted toward 6.43. HARPS was not touched.

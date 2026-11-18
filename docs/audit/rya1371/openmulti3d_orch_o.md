# OpenMULTI3D in the ORCH: solar O I, red-optical (RYA-1371)

Ryan, 2026-11-18: OpenMULTI3D is a publication engine wherever it has a model atom.
This note covers the first ORCH products made with it: O I, 6910-9199 Å, on the three
holdings the RYA-1230 O runs used.

## What was wired

- `pipeline/openmulti3d_nlte.py` supplies departure coefficients b(level, depth). They are
  solved by OpenMULTI3D (pinned build) on the MPG Keeper O atom, at the same solar MARCS
  model the Gerber 1D decks use, over an abundance grid A(O) = 8.30-9.10.
  `scripts/build_openmulti3d_deck.py` builds the deck. It takes about 6 s and verifies
  itself on load (atom hash, engine build, atmosphere, deck hash).
- The departures are matched to the MARCS layers **by column mass**. The OpenMULTI3D and
  MARCS optical-depth scales differ by up to 74 K in T(τ). Column mass agrees to 10.5 K at
  most. One MARCS bottom layer lies outside the engine mesh and is held at b = 1.
- Everything else is the ORCH's own: the Turbospectrum fit, the GESv6 line list, the
  continuum rule, the windows and the line set (AGSS21 O grid).
  `derive_band_products.py --engine-b-deck openmulti3d-nlte | openmulti3d-1d-lte`
  produces the pair:
  - `synth-1D-NLTE-openmulti3d` (model 11);
  - `synth-1D-LTE-openmulti3d` (model 10). This is the same setup with the departures
    switched off, so the NLTE effect is measured on one atmosphere.
- The orchestrator dispatches the pair. `pipeline/run_matrix.DECK_EMITS` binds
  `openmulti3d-nlte` / `openmulti3d-1d-lte` to models 11/10. `DECK_PAIRS` makes them a
  mandatory pair. A cell is held, with its reason, when the element has no OpenMULTI3D
  deck, the star is not the Sun, or the host has no atom/engine to build the deck. The
  EW route refuses the deck rather than run in LTE under its name.
- `--openmulti3d-deck-key O@o41j` runs the Belyaev H-collision atom as a named variant. It
  is never blended into the nominal product.

## Results

Per-line values, A(O):

| Holding | Treatment | 7771.94 | 7774.17 | 7775.39 | 8446.36 | Product |
|---|---|---|---|---|---|---|
| IAG | 1D-LTE OM | 8.923 | 8.888 | 8.866 | 8.793 | 8.877 (n=4) |
| IAG | 1D-NLTE OM o41f | 8.738 | 8.737 | 8.752 | excluded | **8.738** (n=3) |
| IAG | 1D-NLTE OM o41j | 8.656 | 8.659 | 8.677 | excluded | 8.659 (n=3) |
| KP kurucz2005 | 1D-LTE OM | 8.928 | 8.885 | 8.871 | 8.784 | 8.878 |
| KP kurucz2005 | 1D-NLTE OM o41f | 8.738 | 8.738 | 8.746 | excluded | **8.738** |
| KP kurucz2005 | 1D-NLTE OM o41j | 8.658 | 8.656 | 8.682 | excluded | 8.658 |
| KP molecfit | 1D-LTE OM | 8.935 | 8.889 | 8.867 | 8.777 | 8.878 |
| KP molecfit | 1D-NLTE OM o41f | 8.747 | 8.744 | 8.747 | excluded | **8.747** |
| KP molecfit | 1D-NLTE OM o41j | 8.665 | 8.660 | 8.682 | excluded | 8.665 |

Three things are worth noting:

- **NLTE pulls the triplet into agreement.** In LTE the three lines spread over 0.057 dex;
  in NLTE they spread over 0.015.
- **The NLTE effect matches the direct fit.** On IAG the paired triplet effect is −0.150
  dex; the direct OpenMULTI3D fit (`openmulti3d_solar_o_direct_fit.md`) gave −0.149. The
  absolute levels differ by about 0.045 dex (ORCH 8.74 vs direct 8.79). That gap comes from
  where the continuum is placed — the ORCH synthesis reference versus the observed upper
  envelope — not from the departures.
- **The H-collision systematic is 0.08 dex** (o41f − o41j). The direct fit gave 0.07.

## 8446 is excluded, and why the guard missed it at first

GESv6 flags O I 8446 `nlte=T` but writes its level labels as `x`. iSpec passes `x` to
Turbospectrum unchanged. No atom level is called `x`, so Turbospectrum applies b = 1 and
the line runs in LTE. Every label guard counted only `none` as "no label", so the first
run counted 8446 as labelled and put an LTE line into the NLTE aggregate. Its NLTE − LTE
shift was −0.007 dex, against −0.15 for the triplet.

The fix is `gerber_nlte.is_label`: `""`, `none` and `x` all mean no label. It is used at
every label check (`gerber_nlte`, `nlte_label_sidecar`, `run_matrix`). 8446 is now
dropped from the NLTE aggregate and kept in the per-line record. The window coverage
report changed from 21/46 labelled O I lines to 3/46.

### The same defect affects published Mg products (finished element; not rerun)

GESv6 has 235 Mg I lines flagged `nlte=T` with `x` labels. Mg I lines from that set sit
inside the published Gerber `ENGINE-B-NLTE` aggregates:
- 6318.717 and 6319.237 in the 4200-6910 products, on all four holdings;
- 7657.603, 7659.152, 7659.901, 8712.689, 8717.825 and 8923.569 in the 6910-9199 products.

In KP kurucz2005 6910-9199 GRADED:
- 7659.152, 7659.901 and 8717.825 give exactly the same value in NLTE and in
  `synth-1D-LTE-gerber`. They are LTE lines under the NLTE label.
- The labelled line 8806.757 shifts by −0.077 dex, for comparison.

Whether to re-derive is Ryan's call. The guard is fixed, so any future Mg run excludes
these lines.

## Artifacts

- `data/results/rya1371/openmulti3d_orch_o_v1/<holding>_<deck>/`: lines, products,
  provenance and budgets for each run.
- Decks: `grids.openmulti3d_decks/O_5cc8a84d`, `O_o41j_9dc06ba0` (Mac). Rebuild with
  `scripts/build_openmulti3d_deck.py --key O|O@o41j`.

## Open

- **Headline continuum.** The ORCH rule gives 8.74; the envelope gives 8.79. This is Ryan's call.
- **Other run elements.** Na, Mg, Al, Si, Ca, Ti, Cr and Fe still need decks. Their Keeper
  atoms have to be copied to `engines.openmulti3d_atoms` first.

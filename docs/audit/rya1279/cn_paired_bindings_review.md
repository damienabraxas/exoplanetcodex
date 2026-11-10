# Paired native-column H controls and background component closure

2026-11-08; RYA-1281/1278; scientific HOLD. Start07e5c41e; worktree codex/rya1279,
branch ryandamienschmitt/rya-1279-1280-bride-oxygen-gates, stacked base4635c742.
Paired plan/code5e993aa4; source-component check predeclared2227dccb; input-domain
review fix5201fabc. No formal solution, RT/SE or abundance fit.

## Explicit paired state

Two named controls retain56depths from the held native5750/logg4.5column:

- model_mass_H: original MARCS/LW mass-density/H-fraction normalization;
- native_network_H: exported native full-implemented-network H nuclei normalization.

Common fields are exactly equal: native exported T and Pe, original geometry,
source xi1000m/s,zero velocity and elemental composition. Electron density uses
native Pe converted to SI divided by explicit installed LW k and native T in both
cases. This keeps the electron-pressure convention common; it does not silently
substitute native k in LW. Density owners remain metadata, not a chosen correction.
Original model92abundances are retained with O8.70as in the native experiment.
The same held41level O atom and passive atom factories are read anew per case;
source/model/state/backend/constants fingerprints are recorded and protected.

The LW molecular-population table is explicitly empty. This does **not** remove
all molecule-associated continua: the built-in H/proton free-free pair process
still uses atomic H populations. The actual baseline lookup inventory is CH,H,H-,
H2,OH, never CN. No new CN population provider is attached; the synthetic-only
provider guard remains unchanged.

LW Context initializes background and LTE starting populations on the same4001
vacuum wavelength points converted once from AIR7763–7783Aat0.005A. No radiative
transfer/formal solve or statistical-equilibrium iteration is called. Background
arrays only; not all VIS or other bands, not a new oxygen spectrum.

## Population/background response

| Diagnostic | Result |
|---|---:|
|Native-network/model-mass H and total LTE O ratio|0.99196759–1.01906076|
|O/H ratio identity maximum relative error|1.11022e-15|
|Normalized LTE O level fraction max absolute change|6.66134e-16|
|Total background extinction ratio over window/depth|0.99204504–1.01931938|
|Background thermal emission ratio|0.99196750–1.01932185|
|Background scattering ratio|0.99371445–1.00024798|

Each O particle budget agrees with abundance*nH within1.23e-15relative against
predeclared16*Nlevels*float64eps arithmetic screen. LTE fractions are unchanged
because T/ne are held, while the total O budget follows H. This does not establish
unchanged NLTE departures or justify a solar abundance correction.

Extinction/emission depart from simple H-proportional scaling by up to0.0253791%/
0.0256213%relative. Scattering responds differently because fixed-ne electron
scattering does not follow total H. LW chi includes scattering; eta is thermal
emission. These were measured, not forced to scale linearly.

## Frozen external CN diagnostic

Only three previously sampled native CN transitions are evaluated, using unchanged
absolute lower populations per m3, source constants and the declared exact-vacuum
width/a/stim rescaling policy from the coordinate check. No isotope/population
factor is applied again. Their common extinction is fixed across both controls;
LTE diagnostic emission uses explicit LW Planck convention. Postprocessed sums are
saved as samples; they never become Context backgrounds or molecular equilibrium.

Max summed sampled-CN/background extinction ranges1.85586e-11–2.24579e-4across
depths (maximum over wavelength at each depth). This is three source lines only,
not all CN/molecular opacity or an EW/flux prediction. No H-dependent CN chemistry
feedback is inferred. Plot source logtau500 is a depth coordinate, not a computed
contribution function. The normalization PNG was visually inspected.

## Source-based component closure without fitted coefficients

Inspected LW Source/Background.cpp:359–384defines the H0*proton free-free term;
its tabulated wavelength/temperature coefficient is interpolated directly from
that source. Source/Background.cpp:125–134defines Thomson scattering; the source
Cython wrapper adds scattering into total chi. No amplitude is inferred/fitted
from the paired-array difference. Tables remain within their supplied domains.
Source H populations are separately reconstructed by the actual LW LTE routine.

Let r=nH_B/nH_A, P_A=the source-table H/proton pair opacity and E=fixed-ne Thomson
opacity. Predictions use:

```
chi_B = r*chi_A + (r*r-r)*P_A + (1-r)*E
eta_B = r*eta_A + (r*r-r)*P_A*B_nu
sca_B = r*sca_A + (1-r)*E
```

Maximum relative closure differences: chi3.33067e-15,eta3.10862e-15,
sca8.88178e-16. All pass the prospective512float64eps=1.13687e-13 arithmetic
screen. The pair term contributes~0.00108%–2.07117%of baseline total extinction;
Thomson contributes~21.52%–99.33%of baseline scattering across sampled depths/
wavelengths. The quadratic and fixed-electron terms explain the array response.

This is source/runtime numerical interoperability, not independent continuum-data
accuracy, EOS validation or source-to-binary rebuild attestation. Inspected source
hashd4d0f92e28ebf0dc5f82f4ffa7e6798acee9363023e2727689ca6e875031f048;
installed backend/constants hashes match the parent control. No engine source or
physical background component was altered.

## Attempts, validation, resources and reproduction

Both pairedv1/v2and componentv1/v2commands exit0. Source review found the new state
helper wrongly rejected valid zero microturbulence, and a private component-domain
guard needed explicit finite-coordinate rejection. Fix5201fabcadds targeted tests;
actual paired xi remains1km/s. Freshv2numeric CSVs, background files, sampled
arrays and figure are byte-identical to v1. Inputs/results from both preserved.
No threshold, H choice, CN strength, rates or physical input tuned. Ledger records
the validation-domain defect as fixed; nonlinear response is not classified a bug.

Preflight2state tests pass0.10s; component parser/interpolation1passes0.10s;
expanded input-domain preflight5passes0.15s. Final31held-LW tests pass3.08s, exit0,
including new controls/components plus prior native inventory/closure, CN coefficient,
synthetic provider, SI/frame and opacity-reader tests. Runtime common-field,
composition, particle-budget and component closure checks are additional evidence.

Acceptedv2main-block2.26576s includes LW import/arrays/figure/output work but excludes
process startup/top-level NumPy import; individual Context setup0.52632/0.13594s.
Darwin whole-process peak322830336bytes, not per-Context or RT memory. OMP/OPENBLAS1,
Context Nthreads1. No Sirius/CI/benchmark/install/download; existing LW config warning
retained. No photons propagated, rate balance tested or model admitted.

After source config/hosts/mac.env:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/bride_rya1118/.venv/bin/python scripts/check_bride_cn_paired_bindings.py --output <fresh-paired-bindings>
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/bride_rya1118/.venv/bin/python scripts/check_bride_H_background_components.py --source <paired-bindings> --output <fresh-components>
```

Final pytest command lists tests/test_bride_binding_controls_rya1281.py,
test_bride_H_background_components_rya1281.py,test_bride_cn_native_binding_rya1281.py,
test_bride_cn_binding_closure_rya1281.py,test_bride_cn_coefficients_rya1281.py,
test_bride_cn_provider_rya1281.py,test_bride_cn_diagnostic_rya1281.py,
test_bride_cn_coordinates_rya1281.py,test_bride_ts_budget_rya1281.py, using the held
LW interpreter and thread environment above. Saved log
/private/tmp/rya1281-cn-paired-bindings-final-tests.log records31passed.

11compact artifacts archive paired/component reports,2summaries,2sampled NPZs,
2bindings,2full O LTE population tables and PNG byte for byte. Largest98351bytes.
Each original full background_window.npz is5410354bytes and remains external under
its original name/hash; not staged, compressed, split or renamed to bypass5MiBguard.
Predeclared eight-point samples are diagnostic products, not replacements for those
raw files. Actual writers/invocations registered; all source/old results preserved.

Fresh library203documents/127accounted/76unindexed/12off-Mac,exit1. No new paper
used to validate the continuum or source tables; no fabricated bibliography entries.

## Next bounded test and gates

Prepare paired O NLTE sensitivity manifests under these same named H controls,
first isolating H without attaching the frozen CN diagnostic sum as an admitted
provider. Preserve source atom/rate/background lineage and independently rebuilt
rate/conservation checks, joint tolerance and budgets. An explicit CN-on/off transfer
intervention requires its own external-state wrapper/identity test; the synthetic
provider's native-purpose refusal is not bypassed.

No preferred normalization, new atmosphere/EOS, measured solar abundance or
independent physical reference is established here. Native41level atom lineage,
primary molecular provenance, chemistry, C/N then Fe/reference star and actual
adaptive3D/all-band goals remain OPEN. No model/default/product/production,
publication or merge admission. Scientific HOLD remains.

Exact final test command:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/bride_rya1118/.venv/bin/python -m pytest -q tests/test_bride_binding_controls_rya1281.py tests/test_bride_H_background_components_rya1281.py tests/test_bride_cn_native_binding_rya1281.py tests/test_bride_cn_binding_closure_rya1281.py tests/test_bride_cn_coefficients_rya1281.py tests/test_bride_cn_provider_rya1281.py tests/test_bride_cn_diagnostic_rya1281.py tests/test_bride_cn_coordinates_rya1281.py tests/test_bride_ts_budget_rya1281.py
```

Integration audit: main7f023fd7merged RYA1311PR615Ca product/budget notes.
Inspected delta changes no pipeline/config/scripts/linelist source used here.
No merge/rebase of this stacked branch. Final MIRROR/version/sequence/whitespace
guards pass. Generator set exactly2436historical diagnostics, zero added/removed;
existing global ledger RYA587blank science_impact and tracker/gold debts remain,
not a global-green claim. Postcommit freshness/ordinary push checks required.

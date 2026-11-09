# Native hydrogen, density and composition binding inventory

2026-11-08, RYA-1281/1278, scientific HOLD. Start60a37a5a, branch
ryandamienschmitt/rya-1279-1280-bride-oxygen-gates, worktree codex/rya1279,
stacked base4635c742. Native inventory plan/codee444c402; separate read-only
unit/convention decomposition predeclared232f8d33. No new oxygen NLTE solve.

## Reporting parity and implemented network

Private reporting-only additions to pinned native bsyn/eqmol source
21269948750e2177fe32c5df38e1d9bb81363374 expose two existing density owners:

- bsyn retains atmospheric ro(k) and gas pressure read from its continuum input;
- eqmol computes rhotsuji and gas pressure from its implemented chemical network
  at the supplied temperature/electron pressure; bsyn does not replace ro(k).

Both are recorded, never silently reconciled. Original/reporting binaries each
run off,CN,all with the same original inputs except result paths: all six spectrum
files are byte-identical to the accepted original spectra. All three reporting
inventories have identical EOS/atmosphere/species/composition values. Native
physical expressions/constants, installed binaries and old inputs remain unchanged.

56depths,83implemented elements,627molecules,959species slots per depth
(4atomic/ion stages per element plus molecules). All53704pressure/mass/H-count
records retained in an uncompressed NPZ. This is the implemented network, not a
claim of complete physical chemistry, isotope-specific EOS or primary validation.
Composition JSON records every included atomic number, mean mass, relative
abundance and native input abundance. Runtime C/N/O/Fe dex approximately
8.390000667/7.780000588/8.700000105/7.450000217; original serialized inputs retained.

Hydrogen nuclei pressure is independently summed over all atomic/ionic H slots
and every molecule's H stoichiometry. Neutral-H fraction is0.958823–0.999823;
molecular H nuclei fraction4.58894e-5–0.00166861. Neutral H or H2 alone cannot
supply total-H normalization. Hsum agrees with the native default-REAL fictitious
H pressure within1.51922e-7relative. Independent pressure/mass aggregation
reconstructs native EOS density within9.83261e-7relative. Prospective positive
summation arithmetic screen gamma_(2N+8)=0.0002296498 for959slots, not a physical
EOS accuracy budget. Independently summed gas pressure agrees with native EOS
pressure within5.14e-8relative. No tolerance or source changed after inspection.

## Actual normalization differences

| Native quantity / original comparison | range across56depths |
|---|---:|
|EOS density / retained atmosphere density|1.00004172–1.02757878|
|EOS gas pressure / input gas pressure|1.00004174–1.02653917|
|Network H nuclei / current MARCS-LW nH|0.99196759–1.01906076|
|Native electron density / MARCS-LW ne|0.999998488–0.999998631|

Original CN coefficient-state T/Pe/Pg/rho fields match this atmosphere export
exactly at all56depths for all three sampled lines. The density/H difference is
therefore in the same native reference experiment, not a substituted atmosphere.
Reference remains5750/logg4.5, not canonical solar5772/4.438. No option is adopted.

## Separate legacy-density/unit decomposition

Native eqmol source uses rho=sum(pressure*mass)*1.2123e-8/T after adding the
free-electron mass term; exported species masses are in amu. Its source comment
identifies the factor with mH/k. The explicit existing LW atomic-mass-unit-over-k
reference gives a different scale. This is recorded as a convention difference,
not silently corrected or declared an independently proved upstream physics bug.

Current comparison constants: native k1.38066e-16erg/K (default REAL), LW
Amu1.6605402e-27kg, LW k1.380658e-23J/K. These are the actual installed LW values,
not substituted modern CODATA constants. Native mean H mass1.007859230amu differs
from LW's1.007975amu. Both remain explicit.

| Factor | range/value |
|---|---:|
|Native legacy density coefficient / explicit LW Amu/k_native|1.00796967670|
|Pressure-weighted mass-unit density / original model density|0.99236428–1.01948510|
|Network H mass fraction|0.73786795–0.73788023|
|Original model H mass fraction|0.73826|
|Network/model H mass fraction|0.99946896–0.99948558|
|LW/native H mean mass|1.00011486719|

Native composition-only1/xmy is0.7378804574; pressure-based H mass fraction and
composition-only value remain distinct, including native convergence/mass/electron
conventions. Factorization of network/model nH into density, H mass fraction and
H mean-mass terms reproduces the reported ratio within4.44e-16relative. This is
same-state dimensional bookkeeping, not an independent chemical solution.

Native CN per-gram populations already divide by retained atmosphere density.
The absolute extinction calculation multiplies by that same density, cancelling
it. Multiplying an existing absolute CN opacity by EOS/atmosphere density would
introduce an inconsistent new scaling. The normalization issue matters for binding
external populations, H-relative abundances and the background; it is not an
abundance correction inferred from a profile residual.

## Verification, resources and reproduction

Native compile7.55291s. Original/reporting off0.65432/1.05935s,
CN0.62582/0.72659s,all0.63376/0.71692s. Native timings exclude preparation/hash
preflight. Caps90scompile/180seach program retained; serial OMP/OPENBLAS1. No new
peak-memory measurement, no benchmark, Sirius, CI, installation or download.
137protected-input check slots pass before/after (not137distinct sources).

Preflight2tests pass0.10s. Final26held-LW tests pass2.59s, exit0: unique/fixed-form
source insertion, H2stoichiometry/incomplete inventory refusal, prior CN coefficient,
synthetic provider, external-state kernel, coordinate and opacity-reader checks.
Inventory and both decomposition commands exit0. A post-processing dtype gap
was fixed as described below; no native execution failure. Density/normalization difference is characterized in the
run bug ledger; no native physics fix or chosen normalization.

After source config/hosts/mac.env, OMP_NUM_THREADS=1/OPENBLAS_NUM_THREADS=1:

```
/Users/ryanschmitt/bride_rya1118/.venv/bin/python scripts/probe_bride_cn_native_binding.py --output <fresh-native-binding>
/Users/ryanschmitt/bride_rya1118/.venv/bin/python scripts/analyze_bride_cn_binding_closure.py --source <native-binding> --output <fresh-binding-closure>
/Users/ryanschmitt/bride_rya1118/.venv/bin/python -m pytest -q tests/test_bride_cn_native_binding_rya1281.py tests/test_bride_cn_coefficients_rya1281.py tests/test_bride_cn_provider_rya1281.py tests/test_bride_cn_diagnostic_rya1281.py tests/test_bride_cn_coordinates_rya1281.py tests/test_bride_ts_budget_rya1281.py
```

Registered native/closure outputs retain source/build/six replays and raw logs.
CN reporting stdout4861921bytes remains under its original name externally, with
report receipt.15compact artifacts archive source block,56row state,83element
composition, complete species NPZ,6inputs and4closure outputs across preservedv1/typedv2 byte for byte.
Largest1074886bytes. NPZ is the producer's new uncompressed species representation,
not a renamed raw log or source binary; original raw files remain untouched.
Actual generators/invocations registered. No existing results/source/FITS/gold changed.

Fresh library audit now203documents/127accounted/76unindexed/12off-Mac, exit1.
Seven additional unindexed documents appeared since the preceding196/69audit;
old audits are historical, not rewritten. No new document was used to claim
partition/EOS/physical correctness or given fabricated bibliography metadata.

## Next bounded diagnostic

Build an explicit56depth external binding with **both** original MARCS-LW H and
native-network H normalizations retained as named controls at the same source T/Pe
and geometry. Keep density owner, constants, composition and absolute CN population
units explicit; never select whichever closes an observed residual. First compare
background arrays and O population normalization under the two bindings, without
changing adopted atmosphere/abundance or treating this as EOS parity. A new O NLTE
solve follows matched-state array/identity checks, not an arbitrary density rescaling.

Current adapter remains synthetic-only/H-only. Extending it requires source receipts
and explicit native-packet scope; no native CN population was attached to a Bride
column in this phase. Primary molecular frame/isotope/partition provenance, complete
chemistry/reference validation, solar abundance and adaptive3D/all-band goals remain
OPEN. C/N then Fe/reference-star sequence preserved; no adoption/production/merge.

Integration audit: shared origin/main advanced toadefcaf0during this phase,
merging RYA1353O/Fe transition census PR614. Inspected delta adds the census
script/dossier, no pipeline/config/linelist source changes used here. This native
H/EOS diagnostic does not duplicate or adopt census grading/theory decisions.
No branch merge/rebase or source-input replacement.

Post-processing self-review: NumPy default-float CSV reading represented the unused
boolean arithmetic_passed column as NaN, although the native report flag was checked
and all computed physical columns/outputs were finite. Fixb32543aa preserves boolean
type and refuses nonfinite/false flags;2regression tests pass. Fresh closurev2
CSV is byte-identical to v1; both reports/outputs retained. This is a parser defect,
not a changed native state, physical conclusion, threshold or spectrum. Final26test
suite includes the two typed-CSV tests; original24test run also passed.

Final repository guards: MIRROR/version/sequence/whitespace pass after census
integration record. Generator diagnostic set exactly2436historical, zero added/
removed. Existing global ledger RYA587blank science_impact and tracker/gold debts
remain; no global-green claim. Postcommit freshness and ordinary pre-push guards
are required before handoff.
